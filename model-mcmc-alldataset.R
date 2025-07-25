library(dplyr)
library(purrr)
library(tidyr)
library(MCMCglmm)
source("Response-functions.R")

d<-read.csv("tmp/d.csv") 
dfpici<-read.csv("tmp/dfpici.csv")
dfmci<-read.csv("tmp/dfmci.csv")

set.seed(100)
iterations=100000

# Step 1: Define models
model_specs <- 
  load_mcmc_model_designs()


#####*********************************************************************######
#####* Effect size all dataset ####

mydata<-d

# Step 2: Fit models and store
model_list <- model_specs %>%
  mutate(model = pmap(list(fixed_effects, random_effects, prior), function(fe, re, pr) {
    metamcmc(mydata,
             variable = "HedgesG",
             se = "HedgesG.err",
             fixed_effects = fe,
             random_effects = re,
             prior = pr,              nitt=iterations)
  }))


# Step 3 (replacing full summary parsing):
model_results <- model_list %>%
  dplyr::mutate(posterior = map(model, posteriormcmc)) %>%
  dplyr::mutate(
    DIC = map_dbl(posterior, ~.x$DIC),
    center = map_dbl(posterior, ~.x$center),
    se = map_dbl(posterior, ~.x$se),
    hpi_lower = map_dbl(posterior, ~.x$hpi_lower),
    hpi_upper = map_dbl(posterior, ~.x$hpi_upper),
    hpi_half_width = map_dbl(posterior, ~.x$hpi_half_width),
    significant = hpi_lower > 0 | hpi_upper < 0,
    summary_string = map_chr(posterior, ~.x$summary_string),
    # Convert fixed/random/prior to character strings
    fixed_str = map_chr(fixed_effects, ~ if (is.null(.x)) "none" else as.character(.x)),
    random_str = map_chr(random_effects, ~ if (is.null(.x)) "none" else paste(.x, collapse = " + ")),
    prior_str = map_chr(prior, function(p) {
      if (is.null(p)) return("none")
      components <- c()
      if (!is.null(p$R)) components <- c(components, paste0("R(nu=", p$R$nu, ")"))
      if (!is.null(p$G)) {
        g_terms <- paste0(names(p$G), "(nu=", sapply(p$G, `[[`, "nu"), ")")
        components <- c(components, g_terms)
      }
      paste(components, collapse = "; ")
    })
  ) %>%
  dplyr::select(
    label, fixed_str, random_str, prior_str,
    DIC, center, se, hpi_lower, hpi_upper, hpi_half_width,
    significant, summary_string
  )

# Step 4: View ranked modelsmodel_results %>%
model_results<-
  model_results %>%
  arrange(DIC)

write.csv(file = "tmp/metareg_model_results_hedge_alldata.csv",model_results)


#####*********************************************************************######
#### Exactly shaw

mshaw<-metamcmc(d,variable = "HedgesG", se = "HedgesG.err", fixed_effects= "Z.Gen.interval +Z.Year.Midpoint", random_effects=c("PaperID"))
# 
# mshaw<-metamcmc(
#   na.omit(dplyr::select(d, HedgesG, HedgesG.err, Z.Gen.interval ,Z.Year.Midpoint ,Early.Sample.N , Recent.Sample.N , N.Loci,PaperID,Binomen.rotl)),
#                 variable = "HedgesG", se = "HedgesG.err",
#                 fixed_effects= "Z.Gen.interval +Z.Year.Midpoint +Early.Sample.N + Recent.Sample.N + N.Loci", 
#                 random_effects=c("PaperID","Binomen.rotl"))
mshaw %>% summary
save(file = "tmp/mshaw.rda",mshaw)


mshawnopim<-metamcmc(
                     dplyr::filter(d,!(Diversity.Type %in% c("He","pi","A","AR"))),
                     variable = "HedgesG", se = "HedgesG.err", fixed_effects= "Z.Gen.interval +Z.Year.Midpoint", random_effects=c("PaperID"))
save(file = "tmp/mshawnopim.rda",mshawnopim)
