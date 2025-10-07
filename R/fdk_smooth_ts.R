#' Smooth MODIS time series
#'
#' Smooth and gapfill modis time series, and insert
#' results into FluxnetLSM netcdf files for data fusion.
#'
#' @param df vector of dates
#' @param variable what MODIS variable is processed (LAI or FPAR)
#' @param nc_file FluxnetLSM netcdf ERA Met file to amend, scaling
#'  the data to the respective time step required
#' @param csv_file a FLUXNET standard CSV file (full path) into which the
#' MODIS data will be included.
#'
#' @return a data frame with smoothed gapfilled time series
#' @export

fdk_smooth_ts <- function(
    df,
    variable,
    nc_file = NA,
    csv_file = NA) {

  if (is.na(nc_file) && is.na(csv_file)) {
    warning("Arguments nc_file and csv_file are both missing. At least one must be provided.")
  } else if (!(is.na(nc_file)) && !(is.na(csv_file))) {
    warning("Arguments nc_file and csv_file are both non-missing. Specify only one.")
  } else {
    # create a new data frame based on the variable selected
    df <- df |>
      rename(
        "values" = tolower(variable)
      )

    # detect outliers
    df <- suppressWarnings(
      fdk_detect_outliers(
        df = df,
        sigma = 1,
        plot = FALSE
      )
    )

    # get mean annual cycle
    # and linearly interpolate and
    # carry forward or backward any
    # remaining gaps
    mean_df <- df |>
      mutate(
        values = ifelse(outlierflag == 1, NA, values)
      ) |>
      group_by(doy) |>
      summarize(
        mean_values = median(values, na.rm = TRUE)
      ) |>
      mutate(
        mean_values = zoo::na.approx(mean_values, na.rm = FALSE),
        mean_values = zoo::na.locf(mean_values, na.rm = FALSE),
        mean_values = zoo::na.locf(mean_values, fromLast = TRUE, na.rm = FALSE)
      )

    # join with original data, and create
    # combined results (filling NA values in the
    # original)
    df <- left_join(df, mean_df, by = "doy") |>
      mutate(
        combined_values = ifelse(is.na(values) | outlierflag == 1, mean_values, values),
        weights = ifelse(is.na(values), 0.1, 1)
      )

    # calculate the optimal span based on
    # a BIC metric (phenocam approach)
    span <- suppressWarnings(
      fdk_optimal_span(
        y = df$combined_values,
        x = as.numeric(df$date),
        weights = df$weights
      )
    )

    # fit the model with the optimal
    # span and export the smooth
    # time series (see modis_tseries below)
    fit <- stats::loess(
      combined_values ~ as.numeric(date),
      span = span,
      weights = weights,
      data = df
    )

    if (!is.na(nc_file)){
      #---- netcdf file to get time vector ----
      # open netcdf file for writing
      site_nc <- ncdf4::nc_open(nc_file, write = TRUE)

      # Get timing info for site
      site_start_time <- ncdf4::ncatt_get(site_nc, "time")$units
      site_time <- ncdf4::ncvar_get(site_nc, "time")
      site_tstep_size <- 86400 / (site_time[2] - site_time[1])

      # Extract year
      start_year <- as.numeric(substr(site_start_time, start = 15, stop = 18))
      site_time <- as.Date(as.POSIXct(site_time, origin = sprintf("%s-01-01", start_year)))

    } else {
      site_time <- read.csv(csv_file) |>
        dplyr::select(TIMESTAMP) |>
        dplyr::distinct() |>
        dplyr::arrange(TIMESTAMP) |>
        pull(TIMESTAMP) |>
        as.Date()
    }

    # predict values using the fit
    # model
    modis_smoothed <- suppressWarnings(
      stats::predict(
        fit,
        as.numeric(site_time),
        se = FALSE
      )
    )

    # this data frame can be appended to CSV or NetCDF file
    df_modis_smooth_raw <- tibble(
      date = site_time,
      modis_smoothed = modis_smoothed
    ) |>
      left_join(
        df |>
          select(date, combined_values),
        by = join_by("date")
      ) |>
      rename(modis_raw = combined_values)

    # # check visually
    # df_modis_smooth_raw |>
    #   ggplot(aes(x = date)) +
    #   geom_point(aes(y = modis_raw)) +
    #   geom_line(aes(y = modis_smoothed))

    # Check that the number of time steps match
    if (length(modis_smoothed) != length(site_time)) {
      stop("MODIS and site time steps don't match")
    }

    # Also check that no missing values
    if (any(is.na(modis_smoothed))) {
      stop("Missing values in final MODIS time series")
    }

    variable_raw <- paste0(variable, "_raw")
    variable_short <- paste0(variable, "_MODIS")
    variable_short_raw <- paste0(variable, "_MODIS_raw")

    #---- write data to netcdf file ----
    if (!is.na(nc_file)){
      message("writing MODIS data to netcdf file")

      # Define variable:
      var <- ncdf4::ncvar_def(
        variable_short,
        "-",
        list(
          site_nc$dim[[1]],
          site_nc$dim[[2]],
          site_nc$dim[[3]]
        ),
        missval = -9999,
        longname = paste0("MODIS 8-daily (LOESS smoothed)", variable)
      )

      var_raw <- ncdf4::ncvar_def(
        variable_short_raw,
        "-",
        list(
          site_nc$dim[[1]],
          site_nc$dim[[2]],
          site_nc$dim[[3]]
        ),
        missval = -9999,
        longname = paste0("MODIS 8-daily (raw)", variable)
      )

      # Add variable and then variable data:
      site_nc <- ncdf4::ncvar_add(site_nc, var)
      ncdf4::ncvar_put(
        site_nc,
        variable_short,
        df_modis_smooth_raw$modis_smoothed
      )

      # Add raw variable and then variable data:
      site_nc <- ncdf4::ncvar_add(site_nc, var_raw)
      ncdf4::ncvar_put(
        site_nc,
        variable_short_raw,
        df_modis_smooth_raw$modis_raw
      )

      # Close file handle
      ncdf4::nc_close(site_nc)

    } else {
      message(paste0("writing MODIS data to CSV file: ", variable_short))
      tmp <- readr::read_csv(csv_file)

      if (variable == "FPAR"){
        if (variable %in% names(tmp)){
          tmp <- tmp |>
            dplyr::select(-FPAR)
        }
        if (variable_raw %in% names(tmp)){
          tmp <- tmp |>
            dplyr::select(-FPAR_raw)
        }
      } else if (variable == "LAI"){
        if (variable %in% names(tmp)){
          tmp <- tmp |>
            dplyr::select(-LAI)
        }
        if (variable_raw %in% names(tmp)){
          tmp <- tmp |>
            dplyr::select(-LAI_raw)
        }
      }

      tmp |>
        dplyr::left_join(
          df_modis_smooth_raw |>
            rename(TIMESTAMP = date),
          by = join_by("TIMESTAMP")
        ) |>
        dplyr::rename(
          !!variable := modis_smoothed,
          !!variable_raw := modis_raw
          ) |>
        readr::write_csv(file = csv_file)
    }

    return(modis_tseries)
  }
}
