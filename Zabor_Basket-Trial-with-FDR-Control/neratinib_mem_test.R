# use the GitHub version from Michael Kane's account - most recent
remotes::install_github("kaneplusplus/basket")
library(basket)
library(dplyr)

# Set up the data based on the 2018 NEJM article
nerat_dat <- 
  tibble(
    cancer_type = 
      c("Lung",
        "Breast",
        "Bladder",
        "Colorectal",
        "Biliary tract",
        "Endometrial",
        "Cervical",
        "Gastroesophageal",
        "Ovarian",
        "Other"),
    enrolled = 
      c(26, 25, 18, 17, 11, 8, 5, 7, 5, 19),
    responses = 
      c(1, 8, 0, 0, 2, 0, 1, 0, 0, 1)
  )

# Perform the analysis.
# Doesn't work with mem_exact - too many baskets
# took ~30 minutes to run
for(i in c(0.5, 0.4, 0.3, 0.2, 0.1)) {
  
  nerat_basket <- 
    mem_mcmc(responses = nerat_dat$responses,
             size = nerat_dat$enrolled,
             name = nerat_dat$cancer_type,
             p0 = 0.10, 
             prior = diag(10)/(1/(1 - i)) + matrix(i, nrow = 10, ncol = 10)
    )
  
  save(nerat_basket, file = here::here("results", 
                                       paste0("nerat_basket_", i, ".rda")))
  
}

i <- 0.25
nerat_basket <- 
  mem_mcmc(responses = nerat_dat$responses,
           size = nerat_dat$enrolled,
           name = nerat_dat$cancer_type,
           p0 = 0.10, 
           prior = diag(10)/(1/(1 - i)) + matrix(i, nrow = 10, ncol = 10)
  )

save(nerat_basket, file = here::here("results", 
                                     paste0("nerat_basket_", i, ".rda")))






# Read in the actual patient level data
# accessed from http://www.cbioportal.org/study/clinicalData?id=summit_2018 on 2019-09-19
pt_nerat_dat <- read.csv(file = here::here("data", "summit_2018_clinical_data.tsv"), header = TRUE, sep = "\t")





