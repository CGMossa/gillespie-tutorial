library(jsonlite)
library(tidyverse)

sir_df <- jsonlite::fromJSON("../gillespie-tutorial/rust/output/rejection_sampling_sir_mf.json" %>% here::here())

sir_df %>%
  list_transpose() %>%
  as_tibble() %>%
  rowid_to_column("repetition") %>%
  unnest() %>%
  pivot_longer(cols = susceptible:recovered) %>%
  identity() %>%
  {
    ggplot(.) +
      geom_step(aes(tick, value, color = name, group = str_c(repetition, name)))
  }
