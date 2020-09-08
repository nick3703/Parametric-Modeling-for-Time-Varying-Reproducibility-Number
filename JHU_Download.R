if (!require(tidyverse)) install.packages('tidyverse')
if (!require(lubridate)) install.packages('lubridate')
if (!require(janitor)) install.packages('janitor')
if (!require(openintro)) install.packages('openintro')
if (!require(usdata)) install.packages('usdata')


library(tidyverse)
library(lubridate)
library(janitor)
library(usdata)
options(scipen = 999)
options(tibble.width = Inf)

read_jhu_file_US <- function(jhu_directory, file_name){
  daily <- read_csv(str_c(jhu_directory, file_name),
                    col_types = "cccccddddddcdd") %>%
    mutate(date = mdy(file_name %>% str_replace(fixed(".csv"), ""))) %>%
    filter(Country_Region == "US") %>%
    mutate(Province_State = case_when(
      Admin2 == "American Samoa" ~ "American Samoa",
      TRUE ~ Province_State
    )) %>% 
    mutate(FIPS = case_when(
      Province_State == "American Samoa" ~ "60999", 
      Province_State == "Puerto Rico" ~ "72999",
      Province_State == "Northern Mariana Islands" ~ "69999", 
      Province_State == "Guam" ~ "66999", 
      Province_State == "Virgin Islands" ~ "78999",
      TRUE ~ FIPS
    )) %>% 
    filter(!(Province_State %in% c("Diamond Princess", "Recovered", "Grand Princess"))) %>%
    # filter(str_sub(FIPS, 1, 2) != "80") %>%
    arrange(Province_State) %>%
    mutate(Admin2 = ifelse(is.na(Admin2), Province_State, Admin2)) %>%
    dplyr::filter(Admin2 != "Unassigned") %>%
    # mutate(FIPS = ifelse(Admin2 == "Unassigned", lag(FIPS) %>% str_sub(1, 2) %>% str_c("000"), FIPS)) %>%
    mutate(FIPS = ifelse(Admin2 == "Kansas City" & Province_State == "Missouri", "29095", FIPS)) %>%
    mutate(FIPS = str_pad(FIPS, 5, "left", "0"))
  
  return(daily)
}


read_jhu_file_pacom <- function(jhu_directory, file_name){
  
  pacom_countries <- read_csv("pacom_countries.csv", col_types = cols()) 
  
  daily <-
    read_csv(str_c(jhu_directory, file_name),
             col_types = "cccccddddddcdd") %>% 
    mutate(date = mdy(file_name %>% str_replace(fixed(".csv"), ""))) %>%
    mutate(Country_Region = case_when(
      Province_State == "Hong Kong" ~ "Hong Kong",
      Province_State == "" ~ Country_Region,
      TRUE ~ Country_Region
    )) %>%
    inner_join(pacom_countries, by = c("Country_Region" = "Country")) %>% 
    mutate(FIPS = fips) %>% 
    mutate(Admin2 = Country_Region, 
           Province_State = Country_Region) %>% 
    group_by(FIPS, Admin2, Province_State, Country_Region, date) %>%
    summarize_at(.vars = c("Confirmed", "Deaths", "Recovered", "Active"), .funs = sum) %>% 
    ungroup()
  
  return(daily)
}

jhu_add_empty_counties <- function(jhu_daily_df, county_df){
  jhu_daily_df %>%
    full_join(county_df, by = c("FIPS" = "county_fips")) %>% 
    replace_na(list(Confirmed = 0, Deaths = 0, Recovered = 0, Active = 0)) %>%
    mutate(date = max(date, na.rm = T))
}

# jhu_files[1] %>%
#   purrr::map(~read_jhu_file_pacom(jhu_directory, .x)) %>%
#   purrr::map_df(~jhu_add_empty_counties(.x, county_pop_pacom)) %>%
#   filter(!complete.cases(.))


get_jhu_daily_data <- function(){
  
  jhu_url_root <- "https://raw.githubusercontent.com/CSSEGISandData/COVID-19/master/csse_covid_19_data/csse_covid_19_daily_reports/"
  filenames <- seq(mdy("03-22-2020"), today() - 1, by = 1) %>% format("%m-%d-%Y") %>% str_c(".csv")
  local_data_folder <- "01_data/jhu_daily_data/"
  
  jhu_data_downloaded <- list.files(local_data_folder)
  
  jhu_data_needed <- filenames[!(filenames %in% jhu_data_downloaded)]
  
  jhu_data_needed %>% walk(~download.file(url = str_c("https://raw.githubusercontent.com/CSSEGISandData/COVID-19/master/csse_covid_19_data/csse_covid_19_daily_reports/", .x),
                                          destfile = str_c("01_data/jhu_daily_data/", .x))                        
  )
  return(length(jhu_data_needed))
}

get_jhu_daily_data()

# load_data_files(): 
#   INPUT: nothing 
# STUFF IT DOES: loads all csvs used in sections 1, 3, and 4
# OUTPUT: named list of dataframes used in sections 1, 3, and 4
# load_data_files()


load_data_files <- function(max_date = lubridate::today()){
  
  # 2. county level census data (not by age)
  county_pop_US <- read_csv("County_Level_Census_Data_territories.csv", col_types = "cccd") %>%
    janitor::clean_names() %>%
    distinct() %>%
    mutate(county_fips = str_pad(county_fips, width = 5, side = "left", pad = "0")) %>%
    dplyr::select(county_fips, state, county_name_x, population = popest2018)
  
  # county_pop_total <- read_csv("01_data/pacom_data/pacom_total_population.csv") %>% 
  #   bind_rows(county_pop_US)
  
  county_pop_pacom <- read_csv("pacom_total_population.csv") %>%
    group_by(county_fips, state, county_name_x) %>%
    summarize(population = sum(population))
  
  county_pop_total <- county_pop_US %>% bind_rows(county_pop_pacom)
  
  # 1. jhu data files
  jhu_directory <- "01_data/jhu_daily_data/"
  jhu_files <- list.files(jhu_directory)
  jhu_latest_file <- jhu_files[length(jhu_files)]
  # territory_fips <- tibble(fips = c("78999","60999","66999","69999","72999"), 
  #                          province_state = c("Virgin Islands", "American Samoa", "Guam", "Northern Mariana Islands", "Puerto Rico"))
  
  jhu_US <- jhu_files %>% 
    purrr::map(~read_jhu_file_US(jhu_directory, .x)) %>%
    purrr::map(~jhu_add_empty_counties(.x, county_pop_US)) %>% 
    bind_rows() %>%
    janitor::clean_names() %>%
    mutate(state = abbr2state(state)) %>%
    mutate(admin2 = ifelse(is.na(admin2), county_name_x, admin2),
           province_state = ifelse(is.na(province_state), state, province_state)) %>%
    dplyr::select(fips, date, county_name = admin2, province_state, country_region, confirmed, deaths, recovered, active)
  
  
  jhu_pacom <- jhu_files %>% 
    purrr::map(~read_jhu_file_pacom(jhu_directory, .x)) %>% 
    purrr::map(~jhu_add_empty_counties(.x, county_pop_pacom)) %>%
    bind_rows() %>%
    janitor::clean_names() %>%
    dplyr::select(fips, date, county_name = admin2, province_state, country_region, confirmed, deaths, recovered, active) %>% 
    mutate(county_name = if_else(is.na(county_name),province_state,county_name),
           country_region = if_else(is.na(country_region),province_state,country_region)) %>%
    left_join(county_pop_pacom, by = c("fips" = "county_fips")) %>%
    mutate(county_name = ifelse(is.na(county_name), state, county_name),
           country_region = ifelse(is.na(country_region), state, country_region),
           province_state = ifelse(is.na(province_state), state, province_state)) %>%
    dplyr::select(-state, -county_name_x, -population)
  
  jhu_all <- jhu_US %>% bind_rows(jhu_pacom)
  
  # jhu_pacom %>% 
  #   # filter(is.na(country_region)) %>% arrange(desc(date))
  #   filter(!complete.cases(.)) 
  #   
  #   # as.data.frame() 
  
  if(!is.null(max_date)){
    
    jhu_all <- jhu_all %>%
      filter(date <= as_date(max_date))
  }
  
  # 3. covid outcomes (section 3)
  

  
  
  data_list <- list(jhu_all = jhu_all)
  
  return(data_list)
}  

jhu_all <- load_data_files()$jhu_all

write_csv(path = str_c(lubridate::today(),"_JHU.csv"),jhu_all)

load_data_files()

source("https://raw.githubusercontent.com/nick3703/Parametric-Modeling-for-Time-Varying-Reproducibility-Number/master/GitPushFunctions.R")
git2r::config(user.name = "nick3703",user.email = "nick3703@hotmail.com")

gitadd()
gitcommit()
gitpush()
