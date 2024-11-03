##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##
##              Combination no year
##
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
rm(list = ls())
gc()

library(ggplot2)
library(openxlsx)
library(dplyr)
library(ipmr)
library(gcookbook)
library(lme4)

## function
inv_logit <- function(int, slope, sv) {
  1/(1 + exp(-(int + slope * sv)))
}

sens_ltre <- function(species = "Bro_ere")
{

spec <- species  
## read data
dat <- read.xlsx("C://Users/mandrzej24/Documents/01_revision_demography/04_data/03_original_data/Demography_data_2018_2022.xlsx") %>% 
  subset(species == spec)

## make the seedling columns numeric
dat$sl_spring_t0 <- as.numeric(dat$sl_spring_t0)
dat$sl_fall_t0 <- as.numeric(dat$sl_fall_t0)
dat$sl_spring_t1 <- as.numeric(dat$sl_spring_t1)

## set the number of seedlings to 0 if they are NA as Na means 0 
dat$sl_spring_t0[is.na(dat$sl_spring_t0)] <- 0
dat$sl_fall_t0[is.na(dat$sl_fall_t0)] <- 0 
dat$sl_spring_t1[is.na(dat$sl_spring_t1)] <- 0 

## load seed data
dat_seeds <- read.csv("C://Users/mandrzej24/Documents/01_revision_demography/04_data/03_original_data/Seed_number_transects.csv") %>% 
  subset(., year != 2022 & species == species) 
dat_seeds$unviable_seeds <- as.numeric(dat_seeds$unviable_seeds)

## combine the unviable and viable column to one... thats how we did it with julia
dat_seeds$nr_seeds <- rowSums(dat_seeds[,c ("Seed_count", "unviable", "unviable_seeds")], na.rm = T)

## calculate the mean seeds per species per year in each treatment
mean_seed <- aggregate(nr_seeds ~ year + management + climate, FUN = mean, na.rm = T, data = dat_seeds)
mean_seed$treatment <- paste(mean_seed$climate, mean_seed$management, sep = "_")

## merge the seed data to the entire data frame
dat_comp <- merge(dat, mean_seed)

## caclulate how many seeds each individual produced at any given time (seeds * flower_all)
dat_comp$seeds <- round(dat_comp$fl_all * dat_comp$nr_seeds)

## to a column which tells us if a ind had flowered 
dat_comp$flower_pr <- ifelse(dat_comp$fl_all >= 1, 1, 0)

##Models for different vital rates that are continuous 
library(emmeans)
surv_mods <- tryCatch(glmer(survival ~ log(size_t0) * climate * management + (1|plot:climate), data = dat_comp, family = "binomial"), 
                      error = function(e){cat("ERROR :", conditionMessage(e), "\n")}) %>% summary()

grow_mods <- tryCatch(lmer(log(size_t1) ~ log(size_t0) * climate * management + (1|plot:climate), data = dat_comp),
                      error = function(e){cat("ERROR :", conditionMessage(e), "\n")}) %>% summary() 
  
flow_mods <- tryCatch(glmer(flower_pr ~ log(size_t0) * climate * management + (1|plot:climate), data = dat_comp, family = "binomial"),
                      error = function(e){cat("ERROR :", conditionMessage(e), "\n")}) %>% summary() 
  
dat_seed <- subset(dat_comp, flower_pr == 1)
seed_mods <- tryCatch(glmer(round(seeds) ~ log(size_t0) * climate * management + (1|plot:climate), data = dat_comp, family = "poisson"),
                      error = function(e){cat("ERROR :", conditionMessage(e), "\n")}) %>% summary()

dat_comp$Treatment <- paste(dat_comp$management, dat_comp$climate)

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##
##         reproduction 2018
##
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## seeds that would become seedlings
## reproduction stuff for 2018 (only subplot till 3)
dat_2018 <- subset(dat_comp, year == 2018 & subplot <= 3)

## number of seedlings for first 3 subplots
sl_spring_t0_18 <- unique(dat_2018[, c("plot", "subplot", "treatment", "sl_spring_t0")]) %>% aggregate(sl_spring_t0 ~ 
                                                                                                          treatment + plot + subplot,
                                                                                                        FUN = sum, na.rm = T, data = .)
sl_fall_t0_18 <- unique(dat_2018[, c("plot", "subplot", "treatment", "sl_fall_t0")]) %>% aggregate(sl_fall_t0 ~ 
                                                                                                      treatment + plot + subplot, 
                                                                                                    FUN = sum, na.rm = T, data = .)
sl_spring_t1_18 <- unique(dat_2018[, c("plot", "subplot", "treatment", "sl_spring_t1")]) %>% aggregate(sl_spring_t1 ~ 
                                                                                                          treatment + plot + subplot,
                                                                                                        FUN = sum, na.rm = T, data = .)
sl_18 <- merge(sl_spring_t0_18, sl_fall_t0_18) %>% merge(., sl_spring_t1_18)

## number of seeds for first 3 subplots
repr_18 <- aggregate(seeds ~ treatment + plot + subplot, FUN = sum, na.rm = T, data = dat_2018) %>% merge(., sl_18)

### temporaly solution for plantago...
  if(species == "Pla_lan")
  {
    pla_rep_gra_fut <- subset(dat_comp, treatment == "future_grazing" & year == 2018)
    pla_rep_gra_amb <- subset(dat_comp, treatment == "ambient_grazing" & year == 2018)

    pla_rep_gra_fut <- aggregate(seeds ~ plot, FUN = sum, na.rm = T, data = pla_rep_gra_fut)
    pla_rep_gra_amb <- aggregate(seeds ~ plot, FUN = sum, na.rm = T, data = pla_rep_gra_amb)

    pla_rep_mow_fut <- subset(dat_comp, treatment == "future_mowing" & year == 2018)
    pla_rep_mow_amb <- subset(dat_comp, treatment == "ambient_mowing" & year == 2018)

    pla_rep_mow_fut <- aggregate(seeds ~ plot, FUN = sum, na.rm = T, data = pla_rep_mow_fut)
    pla_rep_mow_amb <- aggregate(seeds ~ plot, FUN = sum, na.rm = T, data = pla_rep_mow_amb)

  ## again lucky that both have the same number of rows (I mean almost impossible they wouldn't have the same but you never know)
for(i in 1:nrow(repr_18))
  for(j in 1:nrow(pla_rep_gra_amb))
  {
    if(repr_18[i, "plot"] == pla_rep_gra_amb[j, "plot"] & 
       repr_18[i, "treatment"] == "ambient_grazing")
      repr_18[i, "seeds"] <- pla_rep_gra_amb[j, "seeds"]
    
    if(repr_18[i, "plot"] == pla_rep_gra_fut[j, "plot"] & 
       repr_18[i, "treatment"] == "future_grazing")
      repr_18[i, "seeds"] <- pla_rep_gra_fut[j, "seeds"]
    
    if(repr_18[i, "plot"] == pla_rep_mow_amb[j, "plot"] & 
       repr_18[i, "treatment"] == "ambient_mowing")
      repr_18[i, "seeds"] <- pla_rep_mow_amb[j, "seeds"]
    
    if(repr_18[i, "plot"] == pla_rep_mow_fut[j, "plot"] & 
       repr_18[i, "treatment"] == "future_mowing")
      repr_18[i, "seeds"] <- pla_rep_mow_fut[j, "seeds"]
    
    else
      next
  }
}

## seed to seedling
repr_18$se_sl_fall_t0 <- (repr_18$sl_fall_t0 / repr_18$seeds) %>% ifelse(. > 1, 1, .)
repr_18$se_sl_spring_t1 <- (repr_18$sl_spring_t1 /repr_18$seeds) %>% ifelse(. > 1, 1, .)

## seedling survival
repr_18 <- aggregate(new_individual ~ treatment + plot + subplot, FUN = sum, na.rm = T, data = dat_2018) %>% merge(., repr_18)

repr_18$sl_all_t0 <- repr_18$sl_fall_t0 + repr_18$sl_spring_t0

repr_18$new_ind_per_sl <- (repr_18$new_individual/repr_18$sl_all_t0) %>% ifelse(. > 1, 1, .)

repr_18$year <- 2018

repr_18 <- aggregate(repr_18[, c("se_sl_fall_t0", "se_sl_spring_t1", "new_ind_per_sl")], by = list(repr_18$treatment), 
                     FUN = mean, na.rm = T) %>% rename(treatment = Group.1)

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##
##             reproduction for the REST
##
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
dat_rest <- subset(dat_comp, year != 2018)

## think about the way i caclulate the mean maybe do it by plot too
sl_spring_t0_rest <- unique(dat_rest[, c("plot", "subplot", "treatment", "sl_spring_t0")]) %>% aggregate(sl_spring_t0 ~ 
                                                                                                            treatment + plot + 
                                                                                                            subplot,
                                                                                                          FUN = sum, na.rm = T, data = .)
sl_fall_t0_rest <- unique(dat_rest[, c("plot", "subplot",  "year", "treatment", "sl_fall_t0")]) %>% aggregate(sl_fall_t0 ~
                                                                                                                 treatment + plot + 
                                                                                                                 subplot,
                                                                                                               FUN = sum, na.rm = T, data = .)
sl_spring_t1_rest <- unique(dat_rest[, c("plot", "subplot",  "year", "treatment", "sl_spring_t1")]) %>% aggregate(sl_spring_t1 ~
                                                                                                                     treatment + 
                                                                                                                     plot + subplot,
                                                                                                                   FUN = sum, na.rm = T, data = .)
sl_rest <- merge(sl_spring_t0_rest, sl_fall_t0_rest) %>% merge(., sl_spring_t1_rest)

## aggregate for number of seeds
repr_rest <- aggregate(seeds ~ treatment + plot + subplot, FUN = sum, na.rm = T, data = dat_rest) %>% merge(., sl_rest)

##  seed to seedling
repr_rest$se_sl_fall_t0 <- (repr_rest$sl_fall_t0 / repr_rest$seeds) %>% ifelse(. > 1, 1, .)
repr_rest$se_sl_spring_t1 <- (repr_rest$sl_spring_t1 / repr_rest$seeds) %>% ifelse(. > 1, 1, .)

## seedling to new individual
repr_rest <- aggregate(new_individual ~ treatment + plot + subplot, FUN = sum, na.rm = T, data = dat_rest) %>% merge(., repr_rest)

repr_rest$sl_all_t0 <- repr_rest$sl_spring_t0 + repr_rest$sl_fall_t0

repr_rest$new_ind_per_sl <- (repr_rest$new_individual / repr_rest$sl_all_t0) %>% ifelse(. > 1, 1, .)

## aggregate for mean
repr_rest <- aggregate(repr_rest[, c("se_sl_fall_t0", "se_sl_spring_t1", "new_ind_per_sl")], by = list(repr_rest$treatment), 
                       FUN = mean, na.rm = T) %>% rename(treatment = Group.1)

repr_dat <- rbind(repr_18, repr_rest)

### calculate the mean of each reproduction parameter and merge

## get the max sizes and the sd of it for all new individuals
new_ind_size <- aggregate(log(size_t1) ~ treatment, FUN = mean, na.rm = T, 
                          data = dat_comp %>% subset(., new_individual == 1 & log(size_t1) < 2.5))

new_ind_size$sd <- aggregate(log(size_t1) ~ treatment, FUN = sd, na.rm = T, 
                             data = dat_comp %>% subset(., new_individual == 1 & log(size_t1) < 2.5))$'log(size_t1)'

new_ind_size$sd[is.na(new_ind_size$sd)] <- 0 

## get the uper and lower limits
max_t1 <- aggregate(log(size_t1) ~ treatment, FUN = max, na.rm = T, data = dat_comp) %>% rename(size = 'log(size_t1)')

U <- aggregate(log(size_t0) ~ treatment, FUN = max, na.rm = T, data = dat_comp) %>% rename(size = 'log(size_t0)') %>% 
  rbind(., max_t1) %>% aggregate(size ~ treatment, FUN = max, data = .)

min_t1 <- aggregate(log(size_t1) ~ treatment, FUN = min, na.rm = T, data = dat_comp) %>% rename(size = 'log(size_t1)')

L <- aggregate(log(size_t0) ~ treatment, FUN = min, na.rm = T, data = dat_comp) %>% rename(size = 'log(size_t0)') %>% 
  rbind(., min_t1) %>% aggregate(size ~ treatment, FUN = min, data = .)

repr_agg <-  aggregate(repr_dat[, c("se_sl_fall_t0", "se_sl_spring_t1", "new_ind_per_sl")], by = list(repr_dat$treatment), 
                       FUN = mean, na.rm = T) %>% rename(treatment = Group.1)
repr_agg[,c("se_sl_fall_t0_sd", "se_sl_spring_t1_sd", "new_ind_per_sl_sd")] <-  aggregate(repr_dat[, c("se_sl_fall_t0", "se_sl_spring_t1", 
                                                                                                       "new_ind_per_sl")], 
                                                                                          by = list(repr_dat$treatment), 
                                                                                          FUN = sd, na.rm = T)[,c("se_sl_fall_t0", 
                                                                                                                  "se_sl_spring_t1", 
                                                                                                                  "new_ind_per_sl")]
## parameter for ambient grazing
pars_amb_gra <- tryCatch(list(s_int = coef(surv_mods)[1], ## survival intercept
                              s_slope = coef(surv_mods)[2], ## survival slope
                                   
                              g_int = coef(grow_mods)[1], ##  growth intercept
                              g_slope = coef(grow_mods)[2], ## growth slope
                              g_sd = grow_mods$sigma, ## growth sd 
                                   
                              r_r_int = coef(flow_mods)[1], ## flower probability intercept
                              r_r_slope = coef(flow_mods)[2], ## flower probability slope
                                   
                              r_s_int = coef(seed_mods)[1], ## number of seeds intercept
                              r_s_slope = coef(seed_mods)[2], ## number of seeds slope
                                   
                              r_d_mu = new_ind_size[1, "log(size_t1)"], ## size distribution of new plants
                              r_d_sd = new_ind_size[1, "sd"], ## sd of size distribution of new plants
                                   
                              g_i1 = repr_agg[1, "se_sl_fall_t0"], ## number of seeds that turn to seedling in fall
                              g_i2 = repr_agg[1, "se_sl_spring_t1"], ## number of seeds that turn to seedling in spring
                                   
                              e_p = repr_agg[1, "new_ind_per_sl"], ## number of new individual that emerge from 
                                   
                              L = L[1, "size"], ## lower limit size (no worries the size is logged)
                              U = U[1, "size"], ## upper limit size (no worries the size is logged)
                                   
                              n = 200), error = function(e){cat("ERROR :", conditionMessage(e), "\n")})

## parameter for ambient mowing
pars_amb_mow <- tryCatch(data.frame(s_int = coef(surv_mods)[1] + coef(surv_mods)[4], ## survival intercept
                                   s_slope = coef(surv_mods)[2] + coef(surv_mods)[6], ## survival slope
                                   
                                   g_int = coef(grow_mods)[1] + coef(grow_mods)[4], ##  growth intercept
                                   g_slope = coef(grow_mods)[2] + coef(grow_mods)[6], ## growth slope
                                   g_sd = grow_mods$sigma, ## growth sd 
                                   
                                   r_r_int = coef(flow_mods)[1] + coef(flow_mods)[4], ## flower probability intercept
                                   r_r_slope = coef(flow_mods)[2] + coef(flow_mods)[6], ## flower probability slope
                                   
                                   r_s_int = coef(seed_mods)[1] + coef(seed_mods)[4], ## number of seeds intercept
                                   r_s_slope = coef(seed_mods)[2] + coef(seed_mods)[6], ## number of seeds slope
                                   
                                   r_d_mu = new_ind_size[2, "log(size_t1)"], ## size distribution of new plants
                                   r_d_sd = new_ind_size[2, "sd"], ## sd of size distribution of new plants
                                   
                                   g_i1 = repr_agg[2, "se_sl_fall_t0"], ## number of seeds that turn to seedling in fall
                                   g_i2 = repr_agg[2, "se_sl_spring_t1"], ## number of seeds that turn to seedling in spring
                                   
                                   e_p = repr_agg[2, "new_ind_per_sl"], ## number of new individual that emerge from 
                                   
                                   L = L[2, "size"], ## lower limit size (no worries the size is logged)
                                   U = U[2, "size"], ## upper limit size (no worries the size is logged)
                                   
                                   n = 200), error = function(e){cat("ERROR :", conditionMessage(e), "\n")})

## future grazing parameter
pars_fut_gra <- tryCatch(data.frame(s_int = coef(surv_mods)[1] + coef(surv_mods)[3], ## survival intercept
                                   s_slope = coef(surv_mods)[2] + coef(surv_mods)[5], ## survival slope
                                   
                                   g_int = coef(grow_mods)[1] + coef(grow_mods)[3], ##  growth intercept
                                   g_slope = coef(grow_mods)[2] + coef(grow_mods)[5], ## growth slope
                                   g_sd = grow_mods$sigma, ## growth sd 
                                   
                                   r_r_int = coef(flow_mods)[1] + coef(flow_mods)[3], ## flower probability intercept
                                   r_r_slope = coef(flow_mods)[2] + coef(flow_mods)[5], ## flower probability slope
                                   
                                   r_s_int = coef(seed_mods)[1] + coef(seed_mods)[3], ## number of seeds intercept
                                   r_s_slope = coef(seed_mods)[2] + coef(seed_mods)[5], ## number of seeds slope
                                   
                                   r_d_mu = new_ind_size[3, "log(size_t1)"], ## size distribution of new plants
                                   r_d_sd = new_ind_size[3, "sd"], ## sd of size distribution of new plants
                                   
                                   g_i1 = repr_agg[3, "se_sl_fall_t0"], ## number of seeds that turn to seedling in fall
                                   g_i2 = repr_agg[3, "se_sl_spring_t1"], ## number of seeds that turn to seedling in spring
                                   
                                   e_p = repr_agg[3, "new_ind_per_sl"], ## number of new individual that emerge from 
                                   
                                   L = L[3, "size"], ## lower limit size (no worries the size is logged)
                                   U = U[3, "size"], ## upper limit size (no worries the size is logged)
                                   
                                   n = 200), error = function(e){cat("ERROR :", conditionMessage(e), "\n")})

## future mowing parameter
pars_fut_mow <- tryCatch(data.frame(s_int = coef(surv_mods)[1] + coef(surv_mods)[3] + coef(surv_mods)[4] + coef(surv_mods)[7], ## survival intercept
                                   s_slope = coef(surv_mods)[2] + coef(surv_mods)[5] + coef(surv_mods)[6] + coef(surv_mods)[8], ## survival slope
                                   
                                   g_int = coef(grow_mods)[1] + coef(grow_mods)[3] + coef(grow_mods)[4] + coef(grow_mods)[7], ##  growth intercept
                                   g_slope = coef(grow_mods)[2] + coef(grow_mods)[5] + coef(grow_mods)[6] + coef(grow_mods)[8], ## growth slope
                                   g_sd = grow_mods$sigma, ## growth sd 
                                   
                                   r_r_int = coef(flow_mods)[1] + coef(flow_mods)[3] + coef(flow_mods)[4] + coef(flow_mods)[7], ## flower probability intercept
                                   r_r_slope = coef(flow_mods)[2] + coef(flow_mods)[5] + coef(flow_mods)[6] + coef(flow_mods)[8], ## flower probability slope
                                   
                                   r_s_int = coef(seed_mods)[1] + coef(seed_mods)[3] + coef(seed_mods)[4] + coef(seed_mods)[7], ## number of seeds intercept
                                   r_s_slope = coef(seed_mods)[2] + coef(seed_mods)[5] + coef(seed_mods)[6] + coef(seed_mods)[8], ## number of seeds slope
                                   
                                   r_d_mu = new_ind_size[4, "log(size_t1)"], ## size distribution of new plants
                                   r_d_sd = new_ind_size[4, "sd"], ## sd of size distribution of new plants
                                   
                                   g_i1 = repr_agg[4, "se_sl_fall_t0"], ## number of seeds that turn to seedling in fall
                                   g_i2 = repr_agg[4, "se_sl_spring_t1"], ## number of seeds that turn to seedling in spring
                                   
                                   e_p = repr_agg[4, "new_ind_per_sl"], ## number of new individual that emerge from 
                                   
                                   L = L[4, "size"], ## lower limit size (no worries the size is logged)
                                   U = U[4, "size"], ## upper limit size (no worries the size is logged)
                                   
                                   n = 200), error = function(e){cat("ERROR :", conditionMessage(e), "\n")})
  
## create "look up table" for species names
look_up_species <- data.frame(abriv = c("Bro_ere", "Dia_car", "Pla_lan", "Sca_och", "Tra_ori"), 
                              fulln = c("Bromus erectus", "Dianthus carthusianorum", "Plantago lanceolata", 
                                        "Scabiosa ochroleuca", "Tragopogon orientalis"))

###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
###             LTRE sensitivity
###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Preparation
amb_diff <- pars_amb_gra[1:14] - pars_amb_mow[1:14] 
amb_diff <- t(amb_diff) %>% as.data.frame() %>% rename(diff = V1)
amb_diff$pars_man <- rownames(amb_diff)

gra_diff <- pars_amb_gra[1:14] - pars_fut_gra[1:14]
gra_diff <- t(gra_diff) %>% as.data.frame() %>% rename(diff = V1)
gra_diff$pars_man <- rownames(amb_diff)

mow_diff <- pars_amb_mow[1:14] - pars_fut_mow[1:14]
mow_diff <- t(mow_diff) %>% as.data.frame() %>% rename(diff = V1)
mow_diff$pars_man <- rownames(mow_diff)

fut_diff <- pars_fut_gra[1:14] - pars_fut_mow[1:14]
fut_diff <- t(fut_diff) %>% as.data.frame() %>% rename(diff = V1)
fut_diff$pars_man <- rownames(fut_diff)

##~~~~~~~~~~~~~~~~~
##  Ambient
##~~~~~~~~~~~~~~~~~
## IPM from the mean of the pars for ambient Bromus
amb_mean <- bind_rows(pars_amb_gra, pars_amb_mow) %>% summarise_all(mean) %>% unlist() %>% relist(., skeleton = pars_amb_gra) %>% as.list()

my_ipm <- init_ipm(sim_gen = "general", di_dd = "di", det_stoch = "det") %>%
  
  define_kernel(name = "P",
                formula = s * g * d_ht,
                family = "CC",
                
                g = dnorm(ht_2, g_mu, g_sd),
                g_mu = g_int + g_slope * ht_1,
                s = inv_logit(s_int, s_slope, ht_1),
                
                data_list = amb_mean,
                states = list(c("ht")),
                uses_par_sets = F,
                evict_cor = T, 
                evict_fun = truncated_distributions("norm",
                                                    "g")
  ) %>%
  define_kernel(name = "go_discrete",
                formula = r_r * r_s * g_i1 * d_ht,
                
                family = "CD",
                
                r_r = inv_logit(r_r_int, r_r_slope, ht_1),
                r_s = exp(r_s_int * r_s_slope * ht_1),
                
                data_list = amb_mean,
                states = list(c("ht", "b")),
                usues_par_sets = F
  ) %>%
  define_kernel(name = "stay_discrete",
                formula = 0,
                family = "DD",
                states = list(c("b")),
                evict_cor = F)%>%
  define_kernel(name = "leave_discrete",
                formula = e_p * g_i2 * r_d,
                
                r_d = dnorm(ht_2, r_d_mu, r_d_sd),
                family = "DC",
                
                data_list = amb_mean,
                states = list(c("ht", "b")),
                uses_par_set = F, 
                evict_cor = T,
                evict_fun = truncated_distributions("norm", "r_d")
  )

my_ipm <- my_ipm %>%
  define_impl(
    list(
      P              = list(int_rule    = "midpoint",
                            state_start = "ht",
                            state_end   = "ht"),
      go_discrete    = list(int_rule    = "midpoint",
                            state_start = "ht",
                            state_end   = "b"),
      stay_discrete  = list(int_rule    = "midpoint",
                            state_start = "b",
                            state_end   = "b"),
      leave_discrete = list(int_rule    = "midpoint",
                            state_start = "b",
                            state_end   = "ht")
    )
  )


init_pop_vec   <- runif(200)
init_seed_bank <- 20

general_ipm <- my_ipm %>%
  define_domains(
    
    ht = c(amb_mean$L, amb_mean$U, amb_mean$n)
    
  ) %>%
  define_pop_state(
    
    pop_vectors = list(
      n_ht = init_pop_vec,
      n_b  = init_seed_bank
    )
  ) %>%
  make_ipm(iterations = 100,
           usr_funs = list(inv_logit   = inv_logit))

amb_ori_lambda <- lambda(general_ipm)

##~~~~~~~~~~~~~~~~~~
## Sensetivity IPMs
##~~~~~~~~~~~~~~~~~~
sens_amb <- c()

for(i in 1:14)
{
  pars_man <- amb_mean
  pars_man[i] <- pars_man[[i]] + 0.001
  
  my_ipm <- init_ipm(sim_gen = "general", di_dd = "di", det_stoch = "det") %>%
    
    define_kernel(name = "P",
                  formula = s * g * d_ht,
                  family = "CC",
                  
                  g = dnorm(ht_2, g_mu, g_sd),
                  g_mu = g_int + g_slope * ht_1,
                  s = inv_logit(s_int, s_slope, ht_1),
                  
                  data_list = pars_man,
                  states = list(c("ht")),
                  uses_par_sets = F,
                  evict_cor = T, 
                  evict_fun = truncated_distributions("norm",
                                                      "g")
    ) %>%
    define_kernel(name = "go_discrete",
                  formula = r_r * r_s * g_i1 * d_ht,
                  
                  family = "CD",
                  
                  r_r = inv_logit(r_r_int, r_r_slope, ht_1),
                  r_s = exp(r_s_int * r_s_slope * ht_1),
                  
                  data_list = pars_man,
                  states = list(c("ht", "b")),
                  usues_par_sets = F
    ) %>%
    define_kernel(name = "stay_discrete",
                  formula = 0,
                  family = "DD",
                  states = list(c("b")),
                  evict_cor = F)%>%
    define_kernel(name = "leave_discrete",
                  formula = e_p * g_i2 * r_d,
                  
                  r_d = dnorm(ht_2, r_d_mu, r_d_sd),
                  family = "DC",
                  
                  data_list = pars_man,
                  states = list(c("ht", "b")),
                  uses_par_set = F, 
                  evict_cor = T,
                  evict_fun = truncated_distributions("norm", "r_d")
    )
  
  my_ipm <- my_ipm %>%
    define_impl(
      list(
        P              = list(int_rule    = "midpoint",
                              state_start = "ht",
                              state_end   = "ht"),
        go_discrete    = list(int_rule    = "midpoint",
                              state_start = "ht",
                              state_end   = "b"),
        stay_discrete  = list(int_rule    = "midpoint",
                              state_start = "b",
                              state_end   = "b"),
        leave_discrete = list(int_rule    = "midpoint",
                              state_start = "b",
                              state_end   = "ht")
      )
    )
  
  init_pop_vec   <- runif(200)
  init_seed_bank <- 20
  
  general_ipm <- my_ipm %>%
    define_domains(
      
      ht = c(pars_man$L, pars_man$U, pars_man$n)
      
    ) %>%
    define_pop_state(
      
      pop_vectors = list(
        n_ht = init_pop_vec,
        n_b  = init_seed_bank
      )
    ) %>%
    make_ipm(iterations = 100,
             usr_funs = list(inv_logit   = inv_logit))
  
  sens_amb[[i]] <- data.frame(pars_man = names(pars_man)[i], value_man = pars_man[[i]], lambda_man = lambda(general_ipm)) 
  
}

sens_amb <- do.call("rbind", sens_amb)  
sens_amb$lambda_ori <- amb_ori_lambda
sens_amb$sensetivity <- (sens_amb$lambda_man - sens_amb$lambda_ori) / 0.001

sens_amb <- merge(sens_amb, amb_diff, by = "pars_man")

sens_amb$LTRE <- sens_amb$sensetivity * sens_amb$diff

spec_title <- look_up_species[which(look_up_species$abriv == species), 2]
sens_amb$spec <- species
sens_amb$treat <- "amb"

## Sensitivity plot ambient
amb_sens <- ggplot(sens_amb, aes(x = pars_man, y = sensetivity)) + geom_bar(stat = "identity") + theme_bw() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1.15)) + ggtitle(bquote(italic(.(spec_title))~ "grazing ambient * mowing ambient")) + 
  xlab("") + ylab("")

##~~~~~~~~~~~~~~~~~
##  Future
##~~~~~~~~~~~~~~~~~
## IPM from the mean of the pars for futient Bromus
fut_mean <- bind_rows(pars_fut_gra, pars_fut_mow) %>% summarise_all(mean) %>% unlist() %>% relist(., skeleton = pars_fut_gra) %>% as.list()

my_ipm <- init_ipm(sim_gen = "general", di_dd = "di", det_stoch = "det") %>%
  
  define_kernel(name = "P",
                formula = s * g * d_ht,
                family = "CC",
                
                g = dnorm(ht_2, g_mu, g_sd),
                g_mu = g_int + g_slope * ht_1,
                s = inv_logit(s_int, s_slope, ht_1),
                
                data_list = fut_mean,
                states = list(c("ht")),
                uses_par_sets = F,
                evict_cor = T, 
                evict_fun = truncated_distributions("norm",
                                                    "g")
  ) %>%
  define_kernel(name = "go_discrete",
                formula = r_r * r_s * g_i1 * d_ht,
                
                family = "CD",
                
                r_r = inv_logit(r_r_int, r_r_slope, ht_1),
                r_s = exp(r_s_int * r_s_slope * ht_1),
                
                data_list = fut_mean,
                states = list(c("ht", "b")),
                usues_par_sets = F
  ) %>%
  define_kernel(name = "stay_discrete",
                formula = 0,
                family = "DD",
                states = list(c("b")),
                evict_cor = F)%>%
  define_kernel(name = "leave_discrete",
                formula = e_p * g_i2 * r_d,
                
                r_d = dnorm(ht_2, r_d_mu, r_d_sd),
                family = "DC",
                
                data_list = fut_mean,
                states = list(c("ht", "b")),
                uses_par_set = F, 
                evict_cor = T,
                evict_fun = truncated_distributions("norm", "r_d")
  )

my_ipm <- my_ipm %>%
  define_impl(
    list(
      P              = list(int_rule    = "midpoint",
                            state_start = "ht",
                            state_end   = "ht"),
      go_discrete    = list(int_rule    = "midpoint",
                            state_start = "ht",
                            state_end   = "b"),
      stay_discrete  = list(int_rule    = "midpoint",
                            state_start = "b",
                            state_end   = "b"),
      leave_discrete = list(int_rule    = "midpoint",
                            state_start = "b",
                            state_end   = "ht")
    )
  )


init_pop_vec   <- runif(200)
init_seed_bank <- 20

general_ipm <- my_ipm %>%
  define_domains(
    
    ht = c(fut_mean$L, fut_mean$U, fut_mean$n)
    
  ) %>%
  define_pop_state(
    
    pop_vectors = list(
      n_ht = init_pop_vec,
      n_b  = init_seed_bank
    )
  ) %>%
  make_ipm(iterations = 100,
           usr_funs = list(inv_logit   = inv_logit))

fut_ori_lambda <- lambda(general_ipm)

##~~~~~~~~~~~~~~~~~~
## Sensetivity IPMs
##~~~~~~~~~~~~~~~~~~
sens_fut <- c()

for(i in 1:14)
{
  pars_man <- fut_mean
  pars_man[i] <- pars_man[[i]] + 0.001
  
  my_ipm <- init_ipm(sim_gen = "general", di_dd = "di", det_stoch = "det") %>%
    
    define_kernel(name = "P",
                  formula = s * g * d_ht,
                  family = "CC",
                  
                  g = dnorm(ht_2, g_mu, g_sd),
                  g_mu = g_int + g_slope * ht_1,
                  s = inv_logit(s_int, s_slope, ht_1),
                  
                  data_list = pars_man,
                  states = list(c("ht")),
                  uses_par_sets = F,
                  evict_cor = T, 
                  evict_fun = truncated_distributions("norm",
                                                      "g")
    ) %>%
    define_kernel(name = "go_discrete",
                  formula = r_r * r_s * g_i1 * d_ht,
                  
                  family = "CD",
                  
                  r_r = inv_logit(r_r_int, r_r_slope, ht_1),
                  r_s = exp(r_s_int * r_s_slope * ht_1),
                  
                  data_list = pars_man,
                  states = list(c("ht", "b")),
                  usues_par_sets = F
    ) %>%
    define_kernel(name = "stay_discrete",
                  formula = 0,
                  family = "DD",
                  states = list(c("b")),
                  evict_cor = F)%>%
    define_kernel(name = "leave_discrete",
                  formula = e_p * g_i2 * r_d,
                  
                  r_d = dnorm(ht_2, r_d_mu, r_d_sd),
                  family = "DC",
                  
                  data_list = pars_man,
                  states = list(c("ht", "b")),
                  uses_par_set = F, 
                  evict_cor = T,
                  evict_fun = truncated_distributions("norm", "r_d")
    )
  
  my_ipm <- my_ipm %>%
    define_impl(
      list(
        P              = list(int_rule    = "midpoint",
                              state_start = "ht",
                              state_end   = "ht"),
        go_discrete    = list(int_rule    = "midpoint",
                              state_start = "ht",
                              state_end   = "b"),
        stay_discrete  = list(int_rule    = "midpoint",
                              state_start = "b",
                              state_end   = "b"),
        leave_discrete = list(int_rule    = "midpoint",
                              state_start = "b",
                              state_end   = "ht")
      )
    )
  
  init_pop_vec   <- runif(200)
  init_seed_bank <- 20
  
  general_ipm <- my_ipm %>%
    define_domains(
      
      ht = c(pars_man$L, pars_man$U, pars_man$n)
      
    ) %>%
    define_pop_state(
      
      pop_vectors = list(
        n_ht = init_pop_vec,
        n_b  = init_seed_bank
      )
    ) %>%
    make_ipm(iterations = 100,
             usr_funs = list(inv_logit   = inv_logit))
  
  sens_fut[[i]] <- data.frame(pars_man = names(pars_man)[i], value_man = pars_man[[i]], lambda_man = lambda(general_ipm)) 
  
}

sens_fut <- do.call("rbind", sens_fut)  
sens_fut$lambda_ori <- fut_ori_lambda
sens_fut$sensetivity <- (sens_fut$lambda_man - sens_fut$lambda_ori) / 0.001

sens_fut <- merge(sens_fut, fut_diff, by = "pars_man")

sens_fut$LTRE <- sens_fut$sensetivity * sens_fut$diff
sens_fut$spec <- species
sens_fut$treat <- "fut"

## Sensitivity plot futient
fut_sens <- ggplot(sens_fut, aes(x = pars_man, y = sensetivity)) + geom_bar(stat = "identity") + theme_bw() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1.15)) + ggtitle(bquote(italic(.(spec_title))~ "grazing future * mowing future")) + 
  xlab("") + ylab("")

##~~~~~~~~~~~~~~~~~
##  Grazing
##~~~~~~~~~~~~~~~~~
## IPM from the mean of the pars for futient Bromus
gra_mean <- bind_rows(pars_amb_gra, pars_fut_gra) %>% summarise_all(mean) %>% unlist() %>% relist(., skeleton = pars_amb_gra) %>% as.list()

my_ipm <- init_ipm(sim_gen = "general", di_dd = "di", det_stoch = "det") %>%
  
  define_kernel(name = "P",
                formula = s * g * d_ht,
                family = "CC",
                
                g = dnorm(ht_2, g_mu, g_sd),
                g_mu = g_int + g_slope * ht_1,
                s = inv_logit(s_int, s_slope, ht_1),
                
                data_list = gra_mean,
                states = list(c("ht")),
                uses_par_sets = F,
                evict_cor = T, 
                evict_fun = truncated_distributions("norm",
                                                    "g")
  ) %>%
  define_kernel(name = "go_discrete",
                formula = r_r * r_s * g_i1 * d_ht,
                
                family = "CD",
                
                r_r = inv_logit(r_r_int, r_r_slope, ht_1),
                r_s = exp(r_s_int * r_s_slope * ht_1),
                
                data_list = gra_mean,
                states = list(c("ht", "b")),
                usues_par_sets = F
  ) %>%
  define_kernel(name = "stay_discrete",
                formula = 0,
                family = "DD",
                states = list(c("b")),
                evict_cor = F)%>%
  define_kernel(name = "leave_discrete",
                formula = e_p * g_i2 * r_d,
                
                r_d = dnorm(ht_2, r_d_mu, r_d_sd),
                family = "DC",
                
                data_list = gra_mean,
                states = list(c("ht", "b")),
                uses_par_set = F, 
                evict_cor = T,
                evict_fun = truncated_distributions("norm", "r_d")
  )

my_ipm <- my_ipm %>%
  define_impl(
    list(
      P              = list(int_rule    = "midpoint",
                            state_start = "ht",
                            state_end   = "ht"),
      go_discrete    = list(int_rule    = "midpoint",
                            state_start = "ht",
                            state_end   = "b"),
      stay_discrete  = list(int_rule    = "midpoint",
                            state_start = "b",
                            state_end   = "b"),
      leave_discrete = list(int_rule    = "midpoint",
                            state_start = "b",
                            state_end   = "ht")
    )
  )


init_pop_vec   <- runif(200)
init_seed_bank <- 20

general_ipm <- my_ipm %>%
  define_domains(
    
    ht = c(gra_mean$L, gra_mean$U, gra_mean$n)
    
  ) %>%
  define_pop_state(
    
    pop_vectors = list(
      n_ht = init_pop_vec,
      n_b  = init_seed_bank
    )
  ) %>%
  make_ipm(iterations = 100,
           usr_funs = list(inv_logit   = inv_logit))

gra_ori_lambda <- lambda(general_ipm)

##~~~~~~~~~~~~~~~~~~
## Sensetivity IPMs
##~~~~~~~~~~~~~~~~~~
sens_gra <- c()

for(i in 1:14)
{
  pars_man <- gra_mean
  pars_man[i] <- pars_man[[i]] + 0.001
  
  my_ipm <- init_ipm(sim_gen = "general", di_dd = "di", det_stoch = "det") %>%
    
    define_kernel(name = "P",
                  formula = s * g * d_ht,
                  family = "CC",
                  
                  g = dnorm(ht_2, g_mu, g_sd),
                  g_mu = g_int + g_slope * ht_1,
                  s = inv_logit(s_int, s_slope, ht_1),
                  
                  data_list = pars_man,
                  states = list(c("ht")),
                  uses_par_sets = F,
                  evict_cor = T, 
                  evict_fun = truncated_distributions("norm",
                                                      "g")
    ) %>%
    define_kernel(name = "go_discrete",
                  formula = r_r * r_s * g_i1 * d_ht,
                  
                  family = "CD",
                  
                  r_r = inv_logit(r_r_int, r_r_slope, ht_1),
                  r_s = exp(r_s_int * r_s_slope * ht_1),
                  
                  data_list = pars_man,
                  states = list(c("ht", "b")),
                  usues_par_sets = F
    ) %>%
    define_kernel(name = "stay_discrete",
                  formula = 0,
                  family = "DD",
                  states = list(c("b")),
                  evict_cor = F)%>%
    define_kernel(name = "leave_discrete",
                  formula = e_p * g_i2 * r_d,
                  
                  r_d = dnorm(ht_2, r_d_mu, r_d_sd),
                  family = "DC",
                  
                  data_list = pars_man,
                  states = list(c("ht", "b")),
                  uses_par_set = F, 
                  evict_cor = T,
                  evict_fun = truncated_distributions("norm", "r_d")
    )
  
  my_ipm <- my_ipm %>%
    define_impl(
      list(
        P              = list(int_rule    = "midpoint",
                              state_start = "ht",
                              state_end   = "ht"),
        go_discrete    = list(int_rule    = "midpoint",
                              state_start = "ht",
                              state_end   = "b"),
        stay_discrete  = list(int_rule    = "midpoint",
                              state_start = "b",
                              state_end   = "b"),
        leave_discrete = list(int_rule    = "midpoint",
                              state_start = "b",
                              state_end   = "ht")
      )
    )
  
  init_pop_vec   <- runif(200)
  init_seed_bank <- 20
  
  general_ipm <- my_ipm %>%
    define_domains(
      
      ht = c(pars_man$L, pars_man$U, pars_man$n)
      
    ) %>%
    define_pop_state(
      
      pop_vectors = list(
        n_ht = init_pop_vec,
        n_b  = init_seed_bank
      )
    ) %>%
    make_ipm(iterations = 100,
             usr_funs = list(inv_logit   = inv_logit))
  
  sens_gra[[i]] <- data.frame(pars_man = names(pars_man)[i], value_man = pars_man[[i]], lambda_man = lambda(general_ipm)) 
  
}

sens_gra <- do.call("rbind", sens_gra)  
sens_gra$lambda_ori <- gra_ori_lambda
sens_gra$sensetivity <- (sens_gra$lambda_man - sens_gra$lambda_ori) / 0.001

sens_gra <- merge(sens_gra, gra_diff, by = "pars_man")

sens_gra$LTRE <- sens_gra$sensetivity * sens_gra$diff

sens_gra$spec <- species
sens_gra$treat <- "gra"

gra_sens <- ggplot(sens_gra, aes(x = pars_man, y = sensetivity)) + geom_bar(stat = "identity") + theme_bw() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1.15)) + ggtitle(bquote(italic(.(spec_title))~ "grazing ambient * grazing future")) + 
  xlab("") + ylab("")

##~~~~~~~~~~~~~~~~~
##  Mowing
##~~~~~~~~~~~~~~~~~
## IPM from the mean of the pars for futient Bromus
mow_mean <- bind_rows(pars_amb_mow, pars_fut_mow) %>% summarise_all(mean) %>% unlist() %>% relist(., skeleton = pars_amb_mow) %>% as.list()

my_ipm <- init_ipm(sim_gen = "general", di_dd = "di", det_stoch = "det") %>%
  
  define_kernel(name = "P",
                formula = s * g * d_ht,
                family = "CC",
                
                g = dnorm(ht_2, g_mu, g_sd),
                g_mu = g_int + g_slope * ht_1,
                s = inv_logit(s_int, s_slope, ht_1),
                
                data_list = mow_mean,
                states = list(c("ht")),
                uses_par_sets = F,
                evict_cor = T, 
                evict_fun = truncated_distributions("norm",
                                                    "g")
  ) %>%
  define_kernel(name = "go_discrete",
                formula = r_r * r_s * g_i1 * d_ht,
                
                family = "CD",
                
                r_r = inv_logit(r_r_int, r_r_slope, ht_1),
                r_s = exp(r_s_int * r_s_slope * ht_1),
                
                data_list = mow_mean,
                states = list(c("ht", "b")),
                usues_par_sets = F
  ) %>%
  define_kernel(name = "stay_discrete",
                formula = 0,
                family = "DD",
                states = list(c("b")),
                evict_cor = F)%>%
  define_kernel(name = "leave_discrete",
                formula = e_p * g_i2 * r_d,
                
                r_d = dnorm(ht_2, r_d_mu, r_d_sd),
                family = "DC",
                
                data_list = mow_mean,
                states = list(c("ht", "b")),
                uses_par_set = F, 
                evict_cor = T,
                evict_fun = truncated_distributions("norm", "r_d")
  )

my_ipm <- my_ipm %>%
  define_impl(
    list(
      P              = list(int_rule    = "midpoint",
                            state_start = "ht",
                            state_end   = "ht"),
      go_discrete    = list(int_rule    = "midpoint",
                            state_start = "ht",
                            state_end   = "b"),
      stay_discrete  = list(int_rule    = "midpoint",
                            state_start = "b",
                            state_end   = "b"),
      leave_discrete = list(int_rule    = "midpoint",
                            state_start = "b",
                            state_end   = "ht")
    )
  )


init_pop_vec   <- runif(200)
init_seed_bank <- 20

general_ipm <- my_ipm %>%
  define_domains(
    
    ht = c(mow_mean$L, mow_mean$U, mow_mean$n)
    
  ) %>%
  define_pop_state(
    
    pop_vectors = list(
      n_ht = init_pop_vec,
      n_b  = init_seed_bank
    )
  ) %>%
  make_ipm(iterations = 100,
           usr_funs = list(inv_logit   = inv_logit))

mow_ori_lambda <- lambda(general_ipm)

##~~~~~~~~~~~~~~~~~~
## Sensetivity IPMs
##~~~~~~~~~~~~~~~~~~
sens_mow <- c()

for(i in 1:14)
{
  pars_man <- mow_mean
  pars_man[i] <- pars_man[[i]] + 0.001
  
  my_ipm <- init_ipm(sim_gen = "general", di_dd = "di", det_stoch = "det") %>%
    
    define_kernel(name = "P",
                  formula = s * g * d_ht,
                  family = "CC",
                  
                  g = dnorm(ht_2, g_mu, g_sd),
                  g_mu = g_int + g_slope * ht_1,
                  s = inv_logit(s_int, s_slope, ht_1),
                  
                  data_list = pars_man,
                  states = list(c("ht")),
                  uses_par_sets = F,
                  evict_cor = T, 
                  evict_fun = truncated_distributions("norm",
                                                      "g")
    ) %>%
    define_kernel(name = "go_discrete",
                  formula = r_r * r_s * g_i1 * d_ht,
                  
                  family = "CD",
                  
                  r_r = inv_logit(r_r_int, r_r_slope, ht_1),
                  r_s = exp(r_s_int * r_s_slope * ht_1),
                  
                  data_list = pars_man,
                  states = list(c("ht", "b")),
                  usues_par_sets = F
    ) %>%
    define_kernel(name = "stay_discrete",
                  formula = 0,
                  family = "DD",
                  states = list(c("b")),
                  evict_cor = F)%>%
    define_kernel(name = "leave_discrete",
                  formula = e_p * g_i2 * r_d,
                  
                  r_d = dnorm(ht_2, r_d_mu, r_d_sd),
                  family = "DC",
                  
                  data_list = pars_man,
                  states = list(c("ht", "b")),
                  uses_par_set = F, 
                  evict_cor = T,
                  evict_fun = truncated_distributions("norm", "r_d")
    )
  
  my_ipm <- my_ipm %>%
    define_impl(
      list(
        P              = list(int_rule    = "midpoint",
                              state_start = "ht",
                              state_end   = "ht"),
        go_discrete    = list(int_rule    = "midpoint",
                              state_start = "ht",
                              state_end   = "b"),
        stay_discrete  = list(int_rule    = "midpoint",
                              state_start = "b",
                              state_end   = "b"),
        leave_discrete = list(int_rule    = "midpoint",
                              state_start = "b",
                              state_end   = "ht")
      )
    )
  
  init_pop_vec   <- runif(200)
  init_seed_bank <- 20
  
  general_ipm <- my_ipm %>%
    define_domains(
      
      ht = c(pars_man$L, pars_man$U, pars_man$n)
      
    ) %>%
    define_pop_state(
      
      pop_vectors = list(
        n_ht = init_pop_vec,
        n_b  = init_seed_bank
      )
    ) %>%
    make_ipm(iterations = 100,
             usr_funs = list(inv_logit   = inv_logit))
  
  sens_mow[[i]] <- data.frame(pars_man = names(pars_man)[i], value_man = pars_man[[i]], lambda_man = lambda(general_ipm)) 
  
}

sens_mow <- do.call("rbind", sens_mow)  
sens_mow$lambda_ori <- mow_ori_lambda
sens_mow$sensetivity <- (sens_mow$lambda_man - sens_mow$lambda_ori) / 0.001

sens_mow <- merge(sens_mow, mow_diff, by = "pars_man")

sens_mow$LTRE <- sens_mow$sensetivity * sens_mow$diff
sens_mow$spec <- species
sens_mow$treat <- "mow"

mow_sens <- ggplot(sens_mow, aes(x = pars_man, y = sensetivity)) + geom_bar(stat = "identity") + theme_bw() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1.15)) + ggtitle(bquote(italic(.(spec_title))~ "mowing ambient * mowing future")) + 
  xlab("") + ylab("")

## combine parameters for each treatment
rownames(sens_amb) <- sens_amb$pars_man
rownames(sens_fut) <- sens_fut$pars_man
rownames(sens_gra) <- sens_gra$pars_man
rownames(sens_mow) <- sens_mow$pars_man

LTRE_combined_fut <- data.frame(pars = c("Survival", "Growth", "Reproduction", "Recruitment", "Establishment"), 
           LTRE = c(sum(sens_fut[c("s_int", "s_slope"), "LTRE"]),
                    sum(sens_fut[c("g_int", "g_slope", "g_sd"), "LTRE"]),
                    sum(sens_fut[c("r_r_int", "r_r_slope", "r_s_int", "r_s_slope"), "LTRE"]),
                    sum(sens_fut[c("g_i1", "g_i2"), "LTRE"]),
                    sum(sens_fut[c("r_d_mu", "r_d_sd", "e_p"), "LTRE"])))

LTRE_combined_fut$scaled <- LTRE_combined_fut$LTRE / sum(abs(LTRE_combined_fut$LTRE))

LTRE_combined_amb <- data.frame(pars = c("Survival", "Growth", "Reproduction", "Recruitment", "Establishment"), 
                                LTRE = c(sum(sens_amb[c("s_int", "s_slope"), "LTRE"]),
                                         sum(sens_amb[c("g_int", "g_slope", "g_sd"), "LTRE"]),
                                         sum(sens_amb[c("r_r_int", "r_r_slope", "r_s_int", "r_s_slope"), "LTRE"]),
                                         sum(sens_amb[c("g_i1", "g_i2"), "LTRE"]),
                                         sum(sens_amb[c("r_d_mu", "r_d_sd", "e_p"), "LTRE"])))

LTRE_combined_amb$scaled <- LTRE_combined_amb$LTRE / sum(abs(LTRE_combined_amb$LTRE))

LTRE_combined_gra <- data.frame(pars = c("Survival", "Growth", "Reproduction", "Recruitment", "Establishment"), 
                                LTRE = c(sum(sens_gra[c("s_int", "s_slope"), "LTRE"]),
                                         sum(sens_gra[c("g_int", "g_slope", "g_sd"), "LTRE"]),
                                         sum(sens_gra[c("r_r_int", "r_r_slope", "r_s_int", "r_s_slope"), "LTRE"]),
                                         sum(sens_gra[c("g_i1", "g_i2"), "LTRE"]),
                                         sum(sens_gra[c("r_d_mu", "r_d_sd", "e_p"), "LTRE"])))

LTRE_combined_gra$scaled <- LTRE_combined_gra$LTRE / sum(abs(LTRE_combined_gra$LTRE))

LTRE_combined_mow <- data.frame(pars = c("Survival", "Growth", "Reproduction", "Recruitment", "Establishment"), 
                                LTRE = c(sum(sens_mow[c("s_int", "s_slope"), "LTRE"]),
                                         sum(sens_mow[c("g_int", "g_slope", "g_sd"), "LTRE"]),
                                         sum(sens_mow[c("r_r_int", "r_r_slope", "r_s_int", "r_s_slope"), "LTRE"]),
                                         sum(sens_mow[c("g_i1", "g_i2"), "LTRE"]),
                                         sum(sens_mow[c("r_d_mu", "r_d_sd", "e_p"), "LTRE"])))

LTRE_combined_mow$scaled <- LTRE_combined_mow$LTRE / sum(abs(LTRE_combined_mow$LTRE))
LTRE_combined_mow$treat <- "mow"
LTRE_combined_gra$treat <- "gra"
LTRE_combined_fut$treat <- "fut"
LTRE_combined_amb$treat <- "amb"

LTRE_comb <- rbind(LTRE_combined_amb, LTRE_combined_fut, LTRE_combined_gra, LTRE_combined_mow)
LTRE_comb$spec <- species

ltre_amb <- ggplot(LTRE_combined_amb, aes(x = pars, y = scaled)) + geom_bar(stat = "identity") + theme_bw() + 
  labs(y = expression(paste(Delta, lambda, " (grazing - mowing)")), x = "") + geom_text(aes(x = 0.75, y = 1, label = "A.)"), size = 5) + 
  ggtitle("Ambient") + scale_y_continuous(limits = c(-1, 1)) + 
  theme(axis.text.x = element_text(angle = 45, vjust = 0.65))
ltre_fut <- ggplot(LTRE_combined_fut, aes(x = pars, y = scaled)) + geom_bar(stat = "identity") + theme_bw() + 
  labs(y = expression(paste(Delta, lambda, " (grazing - mowing)")), x = "") + geom_text(aes(x = 0.75, y = 1, label = "B.)"), size = 5) + 
  ggtitle("Future") + scale_y_continuous(limits = c(-1, 1)) + 
  theme(axis.text.x = element_text(angle = 45, vjust = 0.65))
ltre_gra <- ggplot(LTRE_combined_gra, aes(x = pars, y = scaled)) + geom_bar(stat = "identity") + theme_bw() + 
  labs(y = expression(paste(Delta, lambda, " (ambient - future)")), x = "") + geom_text(aes(x = 0.75, y = 1, label = "C.)"), size = 5) + 
  ggtitle("Grazing") + scale_y_continuous(limits = c(-1, 1)) + 
  theme(axis.text.x = element_text(angle = 45, vjust = 0.65))
ltre_mow <- ggplot(LTRE_combined_mow, aes(x = pars, y = scaled)) + geom_bar(stat = "identity") + theme_bw() + 
  labs(y = expression(paste(Delta, lambda, " (ambient - future)")), x = "") + geom_text(aes(x = 0.75, y = 1, label = "D.)"), size = 5) + 
  ggtitle("Mowing") + scale_y_continuous(limits = c(-1, 1)) + 
  theme(axis.text.x = element_text(angle = 45, vjust = 0.65))

ltre_pl <- ggpubr::ggarrange(ltre_amb, ltre_fut, ltre_gra, ltre_mow) %>% 
  ggpubr::annotate_figure(., top = ggpubr::text_grob(spec_title, face = "italic"))

sens_pl <- ggpubr::ggarrange(amb_sens, fut_sens, gra_sens, mow_sens)

sens_dat <- rbind(sens_amb, sens_fut, sens_gra, sens_mow)

return(list(ltre = ltre_pl,sens = sens_pl, data = sens_dat, data_comb = LTRE_comb, ori_dat = dat_comp))

}

bro <- sens_ltre()

dia <- sens_ltre(species = "Dia_car")

pla <- sens_ltre(species = "Pla_lan")

sca <- sens_ltre(species = "Sca_och")

tra <- sens_ltre(species = "Tra_ori")

plot_dat <- rbind(bro$data_comb, dia$data_comb, pla$data_comb, sca$data_comb, tra$data_comb)
sens_dat <- rbind(bro$data, dia$data, pla$data, sca$data, tra$data)

## use new names for the parameter
## survival
sens_dat$pars_man[sens_dat$pars_man == "s_int"] <-      "S Int"
sens_dat$pars_man[sens_dat$pars_man == "s_slope"] <-    "S Slope"
sens_dat$pars_man[sens_dat$pars_man == "g_int"] <-      "G Int"
sens_dat$pars_man[sens_dat$pars_man == "g_slope"] <-    "G Slope"
sens_dat$pars_man[sens_dat$pars_man == "g_sd"] <-       "G SD"
sens_dat$pars_man[sens_dat$pars_man == "r_r_int"] <-    "P Int"
sens_dat$pars_man[sens_dat$pars_man == "r_r_slope"] <-  "P Slope"
sens_dat$pars_man[sens_dat$pars_man == "r_s_int"] <-    "F Int  "
sens_dat$pars_man[sens_dat$pars_man == "r_s_slope"] <-  "F Slope"
sens_dat$pars_man[sens_dat$pars_man == "g_i1"] <-       "rec fall"
sens_dat$pars_man[sens_dat$pars_man == "g_i2"] <-       "rec spring"
sens_dat$pars_man[sens_dat$pars_man == "e_p"] <-        "B"
sens_dat$pars_man[sens_dat$pars_man == "r_d_mu"] <-     "new size"
sens_dat$pars_man[sens_dat$pars_man == "r_d_sd"] <-     "new size sd"

## now structure it like previous paper
sens_dat$pars_man <- factor(sens_dat$pars_man, levels = c("S Int", "S Slope", "G Int", "G Slope", "G SD", "P Int", "P Slope", 
                                     "F Int  ", "F Slope", "rec fall", "rec spring", "B", "new size", "new size sd"))

labeller_treat <- c('amb' = "(a) Ambient",
                       'fut' = "(b) Future",
                       'gra' = "(c) Grazing",
                       'mow' = "(d) Mowing")

sens <- ggplot(sens_dat, aes(x = pars_man, y = sensetivity, fill = spec)) + 
  geom_bar(stat = "identity", position = "dodge", width = 1) + 
  facet_wrap(~ treat, labeller = labeller(treat = labeller_treat)) +
  theme_bw() + 
  scale_fill_viridis_d(labels = c("Bromus erectus", "Dianthus carthusianorum", 
                                  "Plantago lanceolata", "Scabiosa ochroleuca", "Tragopogon orientalis")) + 
  labs(y = expression(paste("Sensitivity of ", lambda)), fill = "Species", x = "") + 
  theme(text = element_text(size = 20),  
        axis.text.x = element_text(angle = 90, vjust = 0.35), 
        strip.text = element_text(size = 20)) + 
  scale_x_discrete(labels = c("rec fall" = expression(paste(theta, "f")),
                              "rec spring" = expression(paste(theta, "s")),
                              "new size" = expression(paste(eta, "")), 
                              "new size sd" = expression(paste(eta, " SD")))) 

#ggsave(filename = "C://Users/marti/Documents/Martin/01_PhD/02_GCEF_demografie_ipm/01_revision/03_Figures/sens_full_model.jpeg", plot = sens, 
 #      dpi = 300, device = "jpeg", width = 4000, height = 2342, units = "px")

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##                        LTRE
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
plot_dat$pars <- factor(plot_dat$pars, levels = c("Survival", "Growth", "Reproduction", "Recruitment", "Establishment"))

ltre_amb <- ggplot(plot_dat %>% subset(treat == "amb"), aes(x = pars, y = scaled, fill = spec)) + 
  geom_bar(stat = "identity", position = "dodge") + 
    theme_bw() + scale_fill_viridis_d(labels = c("Bromus erectus", "Dianthus carthusianorum",
                                                 "Plantago lanceolata", "Scabiosa ochroleuca", 
                                                 "Tragopogon orientalis")) + 
  labs(fill = "Species", y = expression(paste(Delta, lambda, " (grazing - mowing)")),
       x = "", title = "(a) Ambient") + 
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        text = element_text(size = 20)) + 
  scale_y_continuous(limits = c(-1, 1)) 

ltre_fut <- ggplot(plot_dat %>% subset(treat == "fut"), aes(x = pars, y = scaled, fill = spec)) + 
  geom_bar(stat = "identity", position = "dodge") + 
  theme_bw() + scale_fill_viridis_d(labels = c("Bromus erectus", "Dianthus carthusianorum",
                                               "Plantago lanceolata", "Scabiosa ochroleuca", 
                                               "Tragopogon orientalis")) + 
  labs(fill = "Species", y = expression(paste(Delta, lambda, " (grazing - mowing)")), x = "", 
       title = "(b) Future") + 
  theme(axis.text.x = element_blank(), axis.text.y = element_blank(), axis.ticks = element_blank(),
        text = element_text(size = 20)) + 
  scale_y_continuous(limits = c(-1, 1)) 

ltre_gra <- ggplot(plot_dat %>% subset(treat == "gra"), aes(x = pars, y = scaled, fill = spec)) + 
  geom_bar(stat = "identity", position = "dodge") + 
  theme_bw() + scale_fill_viridis_d(labels = c("Bromus erectus", "Dianthus carthusianorum",
                                               "Plantago lanceolata", "Scabiosa ochroleuca", 
                                               "Tragopogon orientalis")) + 
  labs(fill = "Species", y = expression(paste(Delta, lambda, " (ambient - future)")), x = "", 
       title = "(c) Grazing") + 
  theme(axis.text.x = element_text(angle = 90, size = 20, vjust = 0.25),
        text = element_text(size = 20)) + 
  scale_y_continuous(limits = c(-1, 1)) 

ltre_mow <- ggplot(plot_dat %>% subset(treat == "mow"), aes(x = pars, y = scaled, fill = spec)) + 
  geom_bar(stat = "identity", position = "dodge") + 
  theme_bw() + scale_fill_viridis_d(labels = c("Bromus erectus", "Dianthus carthusianorum",
                                               "Plantago lanceolata", "Scabiosa ochroleuca",
                                               "Tragopogon orientalis")) + 
  labs(fill = "Species", y = expression(paste(Delta, lambda, " (ambient - future)")), x = "", 
       title = "(d) Mowing") + 
  theme(axis.text.x = element_text(angle = 90, size = 20, vjust = 0.25), axis.text.y = element_blank(), axis.ticks.y = element_blank(),
        text = element_text(size = 20)) + 
  scale_y_continuous(limits = c(-1, 1))

ltre_fig_comp <- ggpubr::ggarrange(ltre_amb, ltre_fut, ltre_gra, ltre_mow, common.legend = T, 
                                   legend = "top", align = "h")
ltre_fig_comp

ggsave(filename = "C://Users/mandrzej24/Documents/01_revision_demography/03_figures/tre_full_model.jpeg", plot = ltre_fig_comp, dpi = 300,
       device = "jpeg", width = 5000, height = 2928, units = "px")


