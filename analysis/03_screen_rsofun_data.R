library(tidyverse)

df <- readRDS("data/rsofun_driver_data.rds")

# screen for missing gpp data
screening <- df |>
  group_by(sitename) |>
  unnest(forcing) |>
  summarize(
    keep = ifelse(all(is.na(gpp)), FALSE, TRUE)
  )

df <- left_join(df, screening)
df <- df |>
  filter(
    keep
  )

# screen manually, drop funky sites

screen <- readODS::read_ods("data/rsofun_sites.ods") |>
  filter(
    drop
  )

df <- df |>
  filter(
    !(tolower(sitename) %in% tolower(screen$sitename))
  ) |>
  select(
    -keep
  )

# screen manually, drop funky years
screen <- readODS::read_ods("data/rsofun_sites.ods") |>
  filter(
    is.na(drop)
  )

data_fix_years <- df |>
  filter(tolower(sitename) %in% tolower(screen$sitename)) |>
  group_by(sitename) |>
  unnest(forcing) |>
  left_join(screen) |>
  mutate(
    year = as.numeric(format(date, "%Y"))
  ) |>
  filter(
    (year >= start & year <= end)
  ) |>
  select(
    -params_siml,
    -site_info,
    -year,
    -end,
    -start,
    -drop,
    -notes
  ) |>
  nest() |>
  rename(
    forcing = data
  ) |>
  ungroup()

df1 <- df |>
  filter(tolower(sitename) %in% tolower(screen$sitename)) |>
  select(
    -forcing
  ) |>
  left_join(data_fix_years)

df2 <- df |>
  filter(!(tolower(sitename) %in% tolower(screen$sitename)))

data <- bind_rows(df1, df2)

# save data
saveRDS(data, "data/rsofun_driver_data_clean.rds", compress = "xz")

data |>
  group_by(sitename) |>
  do({

    tmp <- .$forcing[[1]] |>
      tidyr::pivot_longer(
        col = !contains("date"),
        names_to = "measurement",
        values_to = "value"
      )

    sitename <- .$sitename[1]

    l <- seq(as.Date("1990/1/1"), as.Date("2023/1/1"), "years")

    p <- ggplot(data = tmp) +
      geom_line(
        aes(
          date,
          value
        ),
        colour = "red"
      ) +
      geom_vline(xintercept = l) +
      labs(
        title = sitename
      ) +
      theme_bw() +
      theme(panel.grid.minor = element_line()
      ) +
      facet_grid(
        measurement ~ .,
        scales = "free"
        )

    ggsave(
      paste0("manuscript/rsofun_input/",sitename,".png"),
      width = 12,
      height = 14,
      dpi = 175
    )

  })