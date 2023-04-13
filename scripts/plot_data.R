# plot cases when things look whacky
library(tidyverse)

data <- read_csv(here::here("data", "cases_hospitalizations_by_county.csv"))

county_region_key <- read_csv(here::here("data", "county_region_key.csv"))

la = data %>% filter(county == "Los Angeles")

looksee = la %>% 
          dplyr::select(date, cases, tests) %>%
          pivot_longer(-date) %>%
          filter(date >= "2023-03-01") %>%
          ggplot(aes(x = date, y = value)) + 
          geom_point() +
          geom_line() +
          facet_wrap(vars(name))



# when we're suspicious of the data ---------------------------------------

library(tidyverse)
library(ckanr)
library(lubridate)
library(splines)
library(fs)

results_dir <- "results"

case_reporting_delay_ecdf <- read_rds("data/case_reporting_delay_ecdf.rds")
death_reporting_delay_ecdf <- read_rds("data/death_delay_ecdf.rds")

prop_omicron_county_dat <-
  read_csv("data/prop_omicron_county_dat.csv") %>%
  # Repeat last prop if data isn't updated
  bind_rows(., crossing(group_by(., county) %>%
                          filter(date == max(date)) %>%
                          select(-date),
                        date = seq(max(.$date), today(), 1)))

# variants_dat <- read_tsv("https://raw.githubusercontent.com/blab/rt-from-frequency-dynamics/master/data/omicron-us/omicron-us_location-variant-sequence-counts.tsv") %>%
#   filter(location == "California") %>%
#   select(-location) %>%
#   distinct() %>%
#   pivot_wider(names_from = variant, values_from = sequences, values_fill = 0) %>%
#   pivot_longer(-date, names_to = "variant", values_to = "sequences") %>%
#   mutate(sequences = sequences + 1) %>%
#   group_by(date) %>%
#   summarize(variant = variant,
#             prop = sequences / sum(sequences)) %>%
#   filter(variant == "Omicron") %>%
#   select(date, prop_omicron = prop)

quiet <- function(x) {
  sink(tempfile())
  on.exit(sink())
  invisible(force(x))
}


ckanr_setup(url="https://data.ca.gov")
ckan <- quiet(ckanr::src_ckan("https://data.ca.gov"))

# get resources
resources <- rbind(resource_search("name:covid-19", as = "table")$results,
                   resource_search("name:hospitals by county", as = "table")$results)


cases_deaths_url <- resources %>% filter(name == "Statewide COVID-19 Cases Deaths Tests") %>% pull(url)
hosp_url <- resources %>% filter(name == "Statewide Covid-19 Hospital County Data") %>% pull(url)

cases <-
  read_csv(cases_deaths_url) %>%
  mutate(date = lubridate::ymd(date),
         deaths = as.integer(deaths),
         reported_cases = as.integer(reported_cases),
         cases = as.integer(cases),
         positive_tests = as.integer(positive_tests),
         total_tests = as.integer(total_tests)) %>%
  select(date,
         cases = cases,
         tests = total_tests,
         deaths,
         county = area) %>%
  arrange(date, county)

la_cases = cases %>% filter(county == "Los Angeles" & date >= "2023-03-01")
