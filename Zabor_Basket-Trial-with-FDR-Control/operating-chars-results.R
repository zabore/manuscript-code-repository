library(tidyverse)
library(magrittr)
library(here)


# this function takes the name of a list (i.e. simres), the element of a list (i.e. "basket"), and the value in that list (i.e. "pp") and returns a 5000x10 dataframe
get_val_df <- function(list, element, val) {
  list %>% 
    map(element) %>% 
    map(val) %>%  
    unlist() %>% 
    matrix(nrow = 5000, ncol = 10, byrow = T) %>% 
    set_colnames(c("Lung", "Breast", "Bladder", "Colorectal",
                   "Biliary tract", "Endometrial", "Cervical",
                   "Gastroesophageal", "Ovarian", "Other")) %>% 
    as_tibble() %>% 
    mutate(simnum = 1:5000) %>% 
    pivot_longer(
      cols = -simnum,
      names_to = "Basket", 
      values_to = val) %>% 
    mutate(Basket = 
             forcats::fct_relevel(
               Basket, 
               "Lung", "Breast", "Bladder", "Colorectal",
               "Biliary tract", "Endometrial", "Cervical",
               "Gastroesophageal", "Ovarian", "Other"))
}

### Null setting
load(here(
  "simulations",
  "operating-chars",
  "results",
  "scenario1-null.rda"))

pp_scenario1_null <- get_val_df(sim_res, "basket", "pp")
pnull_scenario1_null <- get_val_df(sim_res, "frequentist_res", "p_null") %>% 
  group_by(simnum) %>% 
  mutate(p_null_adj = p.adjust(p_null, method = "fdr")) %>% 
  ungroup()

### Alt0 setting
load(here(
  "simulations",
  "operating-chars",
  "results",
  "scenario1-alt0.rda"))

pp_scenario1_alt0 <- get_val_df(sim_res, "basket", "pp")
pnull_scenario1_alt0 <- get_val_df(sim_res, "frequentist_res", "p_null")%>% 
  group_by(simnum) %>% 
  mutate(p_null_adj = p.adjust(p_null, method = "fdr")) %>% 
  ungroup()

### Alt1 setting
load(here(
  "simulations",
  "operating-chars",
  "results",
  "scenario1-alt1.rda"))

pp_scenario1_alt1 <- get_val_df(sim_res, "basket", "pp")
pnull_scenario1_alt1 <- get_val_df(sim_res, "frequentist_res", "p_null") %>% 
  group_by(simnum) %>% 
  mutate(p_null_adj = p.adjust(p_null, method = "fdr")) %>% 
  ungroup()

### Alt2 setting
load(here(
  "simulations",
  "operating-chars",
  "results",
  "scenario1-alt2.rda"))

pp_scenario1_alt2 <- get_val_df(sim_res, "basket", "pp")
pnull_scenario1_alt2 <- get_val_df(sim_res, "frequentist_res", "p_null") %>% 
  group_by(simnum) %>% 
  mutate(p_null_adj = p.adjust(p_null, method = "fdr")) %>% 
  ungroup()

### Alt3 setting
load(here(
  "simulations",
  "operating-chars",
  "results",
  "scenario1-alt3.rda"))

pp_scenario1_alt3 <- get_val_df(sim_res, "basket", "pp")
pnull_scenario1_alt3 <- get_val_df(sim_res, "frequentist_res", "p_null") %>% 
  group_by(simnum) %>% 
  mutate(p_null_adj = p.adjust(p_null, method = "fdr")) %>% 
  ungroup()


# save the results
save(pp_scenario1_null, pnull_scenario1_null, 
     pp_scenario1_alt0, pnull_scenario1_alt0, 
     pp_scenario1_alt1, pnull_scenario1_alt1, 
     pp_scenario1_alt2, pnull_scenario1_alt2, 
     pp_scenario1_alt3, pnull_scenario1_alt3, 
     file = here("simulations",
                 "operating-chars",
                 "results",
                 "scenario1-summary.rda"))
  