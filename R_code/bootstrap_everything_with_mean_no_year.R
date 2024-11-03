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
library(parallel)
library(tidyverse)
library(lme4)

## color platte
rbPalette <- c("#0072B2", "#D55E00")

standard_error <- function(x) sd(x, na.rm = T) / sqrt(length(x))

## read data
dat <- read.xlsx("C://Users/mandrzej24/Documents/01_revision_demography/04_data/03_original_data/Demography_data_2018_2022.xlsx") %>% 
  subset(species %in% c("Bro_ere", "Dia_car", "Pla_lan", "Sca_och", "Tra_ori"))

## get sample size per treatment and species
aggregate(dat$species ~ dat$climate + dat$management, FUN = table)

unique(dat$year)

dat$treatment <- paste(dat$climate, dat$management, sep = "_")
dat$spec_treat <- paste(dat$species, dat$treatment, sep = "_")

## make the seedling columns numeric
dat$sl_spring_t0 <- as.numeric(dat$sl_spring_t0)
dat$sl_fall_t0 <- as.numeric(dat$sl_fall_t0)
dat$sl_spring_t1 <- as.numeric(dat$sl_spring_t1)

## set the number of seedlings to 0 if they are NA as Na means 0 
dat$sl_spring_t0[is.na(dat$sl_spring_t0)] <- 0
dat$sl_fall_t0[is.na(dat$sl_fall_t0)] <- 0 
dat$sl_spring_t1[is.na(dat$sl_spring_t1)] <- 0 

## figure out 
dat$ind_t0 <- ifelse(dat$size_t0 > 0, 1, 0)

## create new fl_pr
dat$flower_pr <- ifelse(dat$fl_all > 0, 1, 0)

## load seed data
dat_seeds <- read.csv("C://Users/marti/Documents/Martin/01_PhD/02_GCEF_demografie_ipm/Seed_count_data/Seed_number_transects.csv") %>% 
  subset(year != 2022) %>% subset(species %in% c("Bro_ere", "Dia_car", "Pla_lan", "Sca_och", "Tra_ori"))
dat_seeds$treatment <- paste(dat_seeds$climate, dat_seeds$management, sep = "_")

## replace the nas
dat_seeds$unviable[is.na(dat_seeds$unviable)] <- 0
dat_seeds$Seed_count[is.na(dat_seeds$Seed_count)] <- 0
dat_seeds$unviable_seeds[is.na(dat_seeds$unviable_seeds)] <- 0

dat_seeds$unviable_seeds <- as.numeric(dat_seeds$unviable_seeds)

## combine the unviable and viable column to one... thats how we did it with julia
dat_seeds$nr_seeds <- rowSums(dat_seeds[,c ("Seed_count", "unviable", "unviable_seeds")], na.rm = T)

suber <- unique(dat$species)

boot_lambda <- function(ii)
{
  library(dplyr)
  library(ipmr)
  library(lme4)
  
  inv_logit <- function(int, slope, sv) {
    1/(1 + exp(-(int + slope * sv)))
  }
  
  ## create empty data frame 
  lamb <- data.frame(lambda = NA, species = NA, treatment = NA)
  all_lamb <- c()

  for(i in 1:length(suber))
  {
    
    set.seed(ii)
    
    ## subset data into the different species, climate and years
    dat_mod <- subset(dat, species == suber[i])
    dat_mod_s <- subset(dat_seeds, species == suber[i])

    dat_mod <- sample_n(dat_mod, nrow(dat_mod), replace = T)
    dat_mod_s <- sample_n(dat_mod_s, nrow(dat_mod_s), replace = T)
    
    dat_mod_s$unviable_seeds <- as.numeric(dat_mod_s$unviable_seeds)
    ## combine the unviable and viable column to one... that's how we did it with Julia
    dat_mod_s$nr_seeds <- rowSums(dat_mod_s[,c ("Seed_count", "unviable", "unviable_seeds")], na.rm = T)
    
    ## calculate the mean seeds per species per year in each treatment
    mean_seed <- aggregate(nr_seeds ~ species + management + climate, FUN = mean,
                           na.rm = T, data = dat_mod_s)
    mean_seed$treatment <- paste(mean_seed$climate, mean_seed$management, sep = "_")
    
    ## merge the seed data to the entire data frame by species and climate
    dat_mod <- merge(dat_mod, mean_seed[,c("species", "climate", "nr_seeds")], 
                     by = c("species", "climate"))
    
    ## calculate number of seeds
    dat_mod$seeds <- dat_mod$fl_all * dat_mod$nr_seeds 
    
    surv_mods <- tryCatch(glmer(survival ~ log(size_t0) * climate * management + (1|plot:climate), data = dat_mod, family = "binomial"), 
                          error = function(e){cat("ERROR :", conditionMessage(e), "\n")}) %>% summary() 
    
    grow_mods <- tryCatch(lmer(log(size_t1) ~ log(size_t0) * climate * management + (1|plot:climate), data = dat_mod),
                          error = function(e){cat("ERROR :", conditionMessage(e), "\n")}) %>% summary() 
    
    flow_mods <- tryCatch(glmer(flower_pr ~ log(size_t0) * climate * management + (1|plot:climate), data = dat_mod, family = "binomial"),
                          error = function(e){cat("ERROR :", conditionMessage(e), "\n")}) %>% summary() 
    
    dat_seed <- subset(dat_mod, flower_pr == 1)
    seed_mods <- tryCatch(glmer(round(seeds) ~ log(size_t0) * climate * management + (1|plot:climate), data = dat_seed, family = "poisson"),
                          error = function(e){cat("ERROR :", conditionMessage(e), "\n")}) %>% summary()
    
        dat_2018 <- subset(dat_mod, year == 2018 & subplot <= 3)
    
    ## number of seedlings for first 3 subplots
    sl_spring_t0_18 <- unique(dat_2018[, c("plot", "subplot", "treatment", "sl_spring_t0")]) %>% 
      aggregate(sl_spring_t0 ~ treatment + plot + subplot, FUN = sum, na.rm = T, data = .)
    
    sl_fall_t0_18 <- unique(dat_2018[, c("plot", "subplot", "treatment", "sl_fall_t0")]) %>% 
      aggregate(sl_fall_t0 ~ treatment + plot + subplot, FUN = sum, na.rm = T, data = .)
    
    sl_spring_t1_18 <- unique(dat_2018[, c("plot", "subplot", "treatment", "sl_spring_t1")]) %>% 
      aggregate(sl_spring_t1 ~ treatment + plot + subplot, FUN = sum, na.rm = T, data = .)
    
    sl_18 <- merge(sl_spring_t0_18, sl_fall_t0_18) %>% merge(., sl_spring_t1_18)
    
    ## number of seeds for first 3 subplots
    repr_18 <- aggregate(seeds ~ treatment + plot + subplot, FUN = sum, na.rm = T, data = dat_2018) %>% 
      merge(., sl_18)
    
    ## seed to seedling
    repr_18$se_sl_fall_t0 <- (repr_18$sl_fall_t0 / repr_18$seeds) %>% ifelse(. > 1, 1, .)
    repr_18$se_sl_spring_t1 <- (repr_18$sl_spring_t1 /repr_18$seeds) %>% ifelse(. > 1, 1, .)
    
    ## seedling survival
    repr_18 <- aggregate(new_individual ~ treatment + plot + subplot, FUN = sum, na.rm = T, data = dat_2018) %>% 
      merge(., repr_18)
    
    repr_18$sl_all_t0 <- repr_18$sl_fall_t0 + repr_18$sl_spring_t0
    
    repr_18$new_ind_per_sl <- (repr_18$new_individual/repr_18$sl_all_t0) %>% ifelse(. > 1, 1, .)
    
    dat_mod_rest <- subset(dat_mod, year != 2018)
    
    ## think about the way i caclulate the mean maybe do it by plot too
    sl_spring_t0_rest <- unique(dat_mod_rest[, c("plot", "subplot", "treatment", "sl_spring_t0")]) %>% 
      aggregate(sl_spring_t0 ~ treatment + plot + subplot, FUN = sum, na.rm = T, data = .)
    
    sl_fall_t0_rest <- unique(dat_mod_rest[, c("plot", "subplot",  "year", "treatment", "sl_fall_t0")]) %>%
      aggregate(sl_fall_t0 ~ treatment + plot + subplot, FUN = sum, na.rm = T, data = .)
    
    sl_spring_t1_rest <- unique(dat_mod_rest[, c("plot", "subplot",  "year", "treatment", "sl_spring_t1")]) %>%
      aggregate(sl_spring_t1 ~ treatment + plot + subplot, FUN = sum, na.rm = T, data = .)
    
    sl_rest <- merge(sl_spring_t0_rest, sl_fall_t0_rest) %>% merge(., sl_spring_t1_rest)
    
    ## aggregate for number of seeds
    repr_rest <- aggregate(seeds ~ treatment + plot + subplot, FUN = sum, na.rm = T, data = dat_mod_rest) %>% 
      merge(., sl_rest)
    
    ## seed to seedling
    repr_rest$se_sl_fall_t0 <- (repr_rest$sl_fall_t0 / repr_rest$seeds) %>% ifelse(. > 1, 1, .)
    repr_rest$se_sl_spring_t1 <- (repr_rest$sl_spring_t1 / repr_rest$seeds) %>% ifelse(. > 1, 1, .)
    
    ## seedling to new individual
    repr_rest <- aggregate(new_individual ~ treatment + plot + subplot, FUN = sum, na.rm = T, data = dat_mod_rest) %>% 
      merge(., repr_rest)
    
    repr_rest$sl_all_t0 <- repr_rest$sl_spring_t0 + repr_rest$sl_fall_t0
    
    repr_rest$new_ind_per_sl <- (repr_rest$new_individual / repr_rest$sl_all_t0) %>% ifelse(. > 1, 1, .)
    
    repr_dat <- rbind(repr_rest, repr_18)
    ## aggregate for mean
    
    repr_dat <- aggregate(repr_dat[, c("se_sl_fall_t0", "se_sl_spring_t1", "new_ind_per_sl")],
                          by = list(repr_dat$treatment), FUN = mean, na.rm = T) %>% 
      rename(treatment = Group.1)
    
    ## calculate the mean of each reproduction parameter and merge
    
    ## get the max sizes and the sd of it for all new individuals
    new_ind_size <- tryCatch(aggregate(log(size_t1) ~ treatment, FUN = mean, na.rm = T, 
                                       data = dat_mod %>% subset(., new_individual == 1 & log(size_t1) < 2.5)), 
                             error = function(e){cat("ERROR :", conditionMessage(e), "\n")})
    
    new_ind_size$sd <- tryCatch(aggregate(log(size_t1) ~ treatment, FUN = sd, na.rm = T, data = dat_mod %>% 
                                            subset(., new_individual == 1 & log(size_t1) < 2.5))$'log(size_t1)', 
                                error = function(e){cat("ERROR :", conditionMessage(e), "\n")}) 
    
    new_ind_size$sd[is.na(new_ind_size$sd)] <- 0 
    
    ## get the upper and lower limits
    max_t1 <- tryCatch(aggregate(log(size_t1) ~ treatment, FUN = max, na.rm = T, data = dat_mod) %>% 
                         rename(size = 'log(size_t1)'), error = function(e){cat("ERROR :", conditionMessage(e), "\n")})
    
    U <- aggregate(log(size_t0) ~ treatment, FUN = max, na.rm = T, data = dat_mod) %>% rename(size = 'log(size_t0)') %>% 
      rbind(., max_t1) %>% aggregate(size ~ treatment, FUN = max, data = .)
    
    min_t1 <- tryCatch(aggregate(log(size_t1) ~ treatment, FUN = min, na.rm = T, data = dat_mod) %>% 
                         rename(size = 'log(size_t1)'), error = function(e){cat("ERROR :", conditionMessage(e), "\n")}) 
    
    L <- aggregate(log(size_t0) ~ treatment, FUN = min, na.rm = T, data = dat_mod) %>% rename(size = 'log(size_t0)') %>% 
      rbind(., min_t1) %>% aggregate(size ~ treatment, FUN = min, data = .)
    
 ## parameter for ambient grazing
 data_list_amb_gra <- tryCatch(list(s_int = coef(surv_mods)[1], ## survival intercept
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
                                    
                                    g_i1 = repr_dat[1, "se_sl_fall_t0"], ## number of seeds that turn to seedling in fall
                                    g_i2 = repr_dat[1, "se_sl_spring_t1"], ## number of seeds that turn to seedling in spring
                                    
                                    e_p = repr_dat[1, "new_ind_per_sl"], ## number of new individual that emerge from 
                                    
                                    L = L[1, "size"], ## lower limit size (no worries the size is logged)
                                    U = U[1, "size"], ## upper limit size (no worries the size is logged)
                                    
                                    n = 200), error = function(e){cat("ERROR :", conditionMessage(e), "\n")})
 
    ## parameter for ambient mowing
    data_list_amb_mow <- tryCatch(list(s_int = coef(surv_mods)[1] + coef(surv_mods)[4], ## survival intercept
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
                               
                                       g_i1 = repr_dat[2, "se_sl_fall_t0"], ## number of seeds that turn to seedling in fall
                                       g_i2 = repr_dat[2, "se_sl_spring_t1"], ## number of seeds that turn to seedling in spring
                               
                                       e_p = repr_dat[2, "new_ind_per_sl"], ## number of new individual that emerge from 
                               
                                       L = L[2, "size"], ## lower limit size (no worries the size is logged)
                                       U = U[2, "size"], ## upper limit size (no worries the size is logged)
                               
                                       n = 200), error = function(e){cat("ERROR :", conditionMessage(e), "\n")})
    
    ## future grazing parameter
    data_list_fut_gra <- tryCatch(list(s_int = coef(surv_mods)[1] + coef(surv_mods)[3], ## survival intercept
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
                                       
                                       g_i1 = repr_dat[3, "se_sl_fall_t0"], ## number of seeds that turn to seedling in fall
                                       g_i2 = repr_dat[3, "se_sl_spring_t1"], ## number of seeds that turn to seedling in spring
                                       
                                       e_p = repr_dat[3, "new_ind_per_sl"], ## number of new individual that emerge from 
                                       
                                       L = L[3, "size"], ## lower limit size (no worries the size is logged)
                                       U = U[3, "size"], ## upper limit size (no worries the size is logged)
                                       
                                       n = 200), error = function(e){cat("ERROR :", conditionMessage(e), "\n")})
    
    ## future mowing parameter
    data_list_fut_mow <- tryCatch(list(s_int = coef(surv_mods)[1] + coef(surv_mods)[3] + coef(surv_mods)[4] + coef(surv_mods)[7], ## survival intercept
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
                                       
                                       g_i1 = repr_dat[4, "se_sl_fall_t0"], ## number of seeds that turn to seedling in fall
                                       g_i2 = repr_dat[4, "se_sl_spring_t1"], ## number of seeds that turn to seedling in spring
                                       
                                       e_p = repr_dat[4, "new_ind_per_sl"], ## number of new individual that emerge from 
                                       
                                       L = L[4, "size"], ## lower limit size (no worries the size is logged)
                                       U = U[4, "size"], ## upper limit size (no worries the size is logged)
                                       
                                       n = 200), error = function(e){cat("ERROR :", conditionMessage(e), "\n")})
    
    if(is.null(data_list_amb_gra) | is.null(data_list_amb_mow) | is.null(data_list_fut_gra) | is.null(data_list_fut_mow))
     next
    
    ## Just ambient treatment
    amb_gra <- do.call("rbind", data_list_amb_gra)
    amb_mow <- do.call("rbind", data_list_amb_mow)
    ambient <- cbind(amb_gra, amb_mow) %>% rowMeans() %>% as.list()
    
    ## Just future treatment
    fut_gra <- do.call("rbind", data_list_fut_gra)
    fut_mow <- do.call("rbind", data_list_fut_mow)
    future <- cbind(fut_gra, fut_mow) %>% rowMeans() %>% as.list()
    
    ## Just grazing treatment
    amb_gra <- do.call("rbind", data_list_amb_gra)
    fut_gra <- do.call("rbind", data_list_fut_gra)
    grazing <- cbind(amb_gra, fut_gra) %>% rowMeans() %>% as.list()
  
    ## Just mowing treatment
    amb_mow <- do.call("rbind", data_list_amb_mow)
    fut_mow <- do.call("rbind", data_list_fut_mow)
    mowing <- cbind(amb_mow, fut_mow) %>% rowMeans() %>% as.list()
    
    ## combine the lists into one list for easier and shorter coding... as if that mattered by now XD
    data_list_full <- list(data_list_amb_gra, data_list_amb_mow, data_list_fut_mow, data_list_fut_gra, ambient, future, grazing, mowing)  
    names(data_list_full) <- c("amb_gra", "amb_mow", "fut_mow", "fut_gra", "amb", "fut", "gra", "mow")
    
    for(j in 1:length(data_list_full))
    {
     data_list <- data_list_full[[j]] 
      
      if(is.null(data_list))
      next
    
      if(is.null(data_list$g_i1))
        next
    
      else
      {
        my_ipm <- init_ipm(sim_gen = "general", di_dd = "di", det_stoch = "det") %>%
        
          define_kernel(name = "P",
                        formula = s * g * d_ht,
                        family = "CC",
                      
                        g = dnorm(ht_2, g_mu, g_sd),
                        g_mu = g_int + g_slope * ht_1,
                        s = inv_logit(s_int, s_slope, ht_1),
                      
                        data_list = data_list,
                        states = list(c("ht")),
                        uses_par_sets = F,
                        evict_cor = T, 
                        evict_fun = truncated_distributions("norm", "g")) %>%
        
          define_kernel(name = "go_discrete",
                        formula = r_r * r_s * d_ht,
                      
                        family = "CD",
                      
                        r_r = inv_logit(r_r_int, r_r_slope, ht_1),
                        r_s = exp(r_s_int * r_s_slope * ht_1),
                        
                        data_list = data_list,
                        states = list(c("ht", "b")),
                        usues_par_sets = F) %>%
        
          define_kernel(name = "stay_discrete",
                        formula = 0,
                        family = "DD",
                        states = list(c("b")),
                        evict_cor = F) %>%
        
          define_kernel(name = "leave_discrete",
                        formula = e_p * g_i1 * g_i2 * r_d,
                      
                        r_d = dnorm(ht_2, r_d_mu, r_d_sd),
                        family = "DC",
                      
                        data_list = data_list,
                        states = list(c("ht", "b")),
                        uses_par_set = F, 
                        evict_cor = T,
                        evict_fun = truncated_distributions("norm", "r_d"))
      
        my_ipm <- my_ipm %>%
          define_impl(
            list(P              = list(int_rule    = "midpoint",
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
                                      state_end   = "ht")))
      
        init_pop_vec   <- runif(200)
        init_seed_bank <- 20
      
        general_ipm <- tryCatch(my_ipm %>%
                                define_domains(
                                  
                                  ht = c(data_list$L, data_list$U, data_list$n)
                                  
                                ) %>%
                                define_pop_state(
                                  
                                  pop_vectors = list(
                                    n_ht = init_pop_vec,
                                    n_b  = init_seed_bank)) %>%
                                
                                make_ipm(iterate = T, iterations = 50,
                                         usr_funs = list(inv_logit   = inv_logit)), error = function(e){cat("ERROR :", conditionMessage(e), "\n")})
      
      if(is.null(general_ipm))
        next
      
      else
        lamb[j,] <- data.frame(lambda = lambda(general_ipm), species = unique(dat_mod$species), 
                               treatment = names(data_list_full)[j])
      }
   
    }
    
    all_lamb[i] <- list(lamb)
  }
  
  return(all_lamb)
}

cl <- makeCluster(detectCores() - 1)

clusterExport(cl, c("dat", "dat_seeds", "suber"))

start_time <- Sys.time()
lambdas <- parSapply(cl, 1:2, FUN = boot_lambda)
Sys.time() - start_time

stopCluster(cl)
 
booted_lambda <- do.call("rbind", lambdas)
booted_lambda <- do.call("rbind", booted_lambda)
## write the data so I dont have to rerun the whole code
#write.csv(x = booted_lambda, file = "C://Users/marti/Documents/Martin/01_PhD/02_GCEF_demografie_ipm/01_revision/04_lambdas_from_full_model/interaction_boot_full_mod_all_treats.csv")

## read booted data
booted_lambda <- read.csv("c://Users/mandrzej24/Documents/01_revision_demography/04_data/01_new_data/interaction_boot_full_mod_all_treats.csv")
mean_lambdas <- aggregate(log(lambda) ~ species + treatment, FUN = mean, data = booted_lambda, na.rm = T) 
colnames(mean_lambdas)[3] <- "lambda"

mean_lambdas$sd <- aggregate(log(lambda) ~ species + treatment, FUN = sd, 
                             data = booted_lambda, na.rm = T)$'log(lambda)'

mean_lambdas$species_prep <- c("'(a) ' * italic('Bromus erectus')",
                              "'(b) ' * italic('Dianthus carthusianorum')", 
                              "'(c) ' * italic('Plantago lanceolata')",
                              "'(d) ' * italic('Scabiosa ochroleuca')",
                              "'(e) ' * italic('Tragopogon orientalis')") 
mean_lambdas$species_prep_poster <-  c("italic('Bromus erectus')",
                                       "italic('Dianthus carthusianorum')", 
                                       "italic('Plantago lanceolata')",
                                       "italic('Scabiosa ochroleuca')",
                                       "italic('Tragopogon orientalis')") 

mean_lambdas$climate <- ifelse(mean_lambdas$treatment %in% c("amb_gra", "amb_mow"), "ambient", 
                               ifelse(mean_lambdas$treatment == "amb", "ambient", "future"))
mean_lambdas$manage <- ifelse(mean_lambdas$treatment %in% c("amb_gra", "fut_gra"), "grazing", 
                              ifelse(mean_lambdas$treatment == "gra", "grazing", "mowing"))

inter_lamb_sd_bro <- ggplot(mean_lambdas %>% subset(!(treatment %in% c("amb", "fut", "gra", "mow")) & species == "Bro_ere"), aes(x = climate, y = lambda, shape = treatment)) +
  geom_line(aes(x = climate, y = lambda, group = manage, linetype = manage), position = position_dodge(.08)) +
  geom_point(position = position_dodge(.08), size = 4, mapping = aes(color = treatment)) + 
  facet_wrap(~ species_prep, labeller = label_parsed) +
  theme_bw() +
  geom_errorbar(aes(ymin = lambda - sd, ymax = lambda + sd, color = treatment), 
                width = .15, position = position_dodge(.08)) + 
  scale_linetype_manual(values = c(2,1), guide = "none") +
  scale_shape_manual(name = "Treatment", values = c(17, 18, 17, 18), 
                     label = c("Ambient grazing", "Ambient mowing", "Future grazing", "Future mowing")) + 
  xlab("") + ylab(expression(paste("log (", lambda, ")"))) +
  theme(text = element_text(size = 20), strip.text = element_text(size = 15), axis.text.x = element_blank(), 
        axis.ticks.x = element_blank()) + 
  scale_color_manual(name = "Treatment", values = c(rep(rbPalette[1], 2), rep(rbPalette[2], 2)), 
                     label = c("Ambient grazing", "Ambient mowing", "Future grazing", "Future mowing")) +
  geom_hline(yintercept = 0, linetype = 2) + scale_y_continuous(limits = c(-2, 2))

inter_lamb_sd_dia <- ggplot(mean_lambdas %>% subset(!(treatment %in% c("amb", "fut", "gra", "mow")) & species == "Dia_car"), aes(x = climate, y = lambda, shape = treatment)) +
  geom_line(aes(x = climate, y = lambda, group = manage, linetype = manage), position = position_dodge(.08)) +
  geom_point(position = position_dodge(.08), size = 4, mapping = aes(color = treatment)) + 
  facet_wrap(~ species_prep, labeller = label_parsed) +
  theme_bw() +
  geom_errorbar(aes(ymin = lambda - sd, ymax = lambda + sd, color = treatment), 
                width = .15, position = position_dodge(.08)) + 
  scale_linetype_manual(values = c(2,1), guide = "none") +
  scale_shape_manual(name = "Treatment", values = c(17, 18, 17, 18), 
                     label = c("Ambient grazing", "Ambient mowing", "Future grazing", "Future mowing")) + 
  xlab("") + ylab("") +
  theme(text = element_text(size = 20), strip.text = element_text(size = 15), axis.text.x = element_blank(), 
        axis.ticks.x = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank()) + 
  scale_color_manual(name = "Treatment", values = c(rep(rbPalette[1], 2), rep(rbPalette[2], 2)), 
                     label = c("Ambient grazing", "Ambient mowing", "Future grazing", "Future mowing")) +
  geom_hline(yintercept = 0, linetype = 2) + scale_y_continuous(limits = c(-2, 2))

inter_lamb_sd_pla <- ggplot(mean_lambdas %>% subset(!(treatment %in% c("amb", "fut", "gra", "mow")) & species == "Pla_lan"), aes(x = climate, y = lambda, shape = treatment)) +
  geom_line(aes(x = climate, y = lambda, group = manage, linetype = manage), position = position_dodge(.08)) +
  geom_point(position = position_dodge(.08), size = 4, mapping = aes(color = treatment)) + 
  facet_wrap(~ species_prep, labeller = label_parsed) +
  theme_bw() +
  geom_errorbar(aes(ymin = lambda - sd, ymax = lambda + sd, color = treatment), 
                width = .15, position = position_dodge(.08)) + 
  scale_linetype_manual(values = c(2,1), guide = "none") +
  scale_shape_manual(name = "Treatment", values = c(17, 18, 17, 18), 
                     label = c("Ambient grazing", "Ambient mowing", "Future grazing", "Future mowing")) + 
  xlab("") + ylab("") +
  theme(text = element_text(size = 20), strip.text = element_text(size = 15), axis.text.y = element_blank(), 
        axis.ticks.y = element_blank()) + 
  scale_color_manual(name = "Treatment", values = c(rep(rbPalette[1], 2), rep(rbPalette[2], 2)), 
                     label = c("Ambient grazing", "Ambient mowing", "Future grazing", "Future mowing")) +
  geom_hline(yintercept = 0, linetype = 2) + scale_y_continuous(limits = c(-2, 2))

inter_lamb_sd_sca <- ggplot(mean_lambdas %>% subset(!(treatment %in% c("amb", "fut", "gra", "mow")) & species == "Sca_och"), aes(x = climate, y = lambda, shape = treatment)) +
  geom_line(aes(x = climate, y = lambda, group = manage, linetype = manage), position = position_dodge(.08)) +
  geom_point(position = position_dodge(.08), size = 4, mapping = aes(color = treatment)) + 
  facet_wrap(~ species_prep, labeller = label_parsed) +
  theme_bw() +
  geom_errorbar(aes(ymin = lambda - sd, ymax = lambda + sd, color = treatment), 
                width = .15, position = position_dodge(.08)) + 
  scale_linetype_manual(values = c(2,1), guide = "none") +
  scale_shape_manual(name = "Treatment", values = c(17, 18, 17, 18), 
                     label = c("Ambient grazing", "Ambient mowing", "Future grazing", "Future mowing")) + 
  xlab("") + ylab(expression(paste("log (", lambda, ")"))) +
  theme(text = element_text(size = 20), strip.text = element_text(size = 15)) + 
  scale_color_manual(name = "Treatment", values = c(rep(rbPalette[1], 2), rep(rbPalette[2], 2)), 
                     label = c("Ambient grazing", "Ambient mowing", "Future grazing", "Future mowing")) +
  geom_hline(yintercept = 0, linetype = 2) + scale_y_continuous(limits = c(-2, 2))

inter_lamb_sd_tra <- ggplot(mean_lambdas %>% subset(!(treatment %in% c("amb", "fut", "gra", "mow")) & species == "Tra_ori"), aes(x = climate, y = lambda, shape = treatment)) +
  geom_line(aes(x = climate, y = lambda, group = manage, linetype = manage), position = position_dodge(.08)) +
  geom_point(position = position_dodge(.08), size = 4, mapping = aes(color = treatment)) + 
  facet_wrap(~ species_prep, labeller = label_parsed, scale = "free") +
  theme_bw() +
  geom_errorbar(aes(ymin = lambda - sd, ymax = lambda + sd, color = treatment), 
                width = .15, position = position_dodge(.08)) + 
  scale_linetype_manual(values = c(2,1), guide = "none") +
  scale_shape_manual(name = "Treatment", values = c(17, 18, 17, 18), 
                     label = c("Ambient grazing", "Ambient mowing", "Future grazing", "Future mowing")) + 
  xlab("") + ylab("") +
  theme(text = element_text(size = 20), strip.text = element_text(size = 15)) + 
  scale_color_manual(name = "Treatment", values = c(rep(rbPalette[1], 2), rep(rbPalette[2], 2)), 
                     label = c("Ambient grazing", "Ambient mowing", "Future grazing", "Future mowing")) +
  geom_hline(yintercept = 0, linetype = 2)

int_lotte_wunsch <- ggpubr::ggarrange(inter_lamb_sd_bro, inter_lamb_sd_dia, inter_lamb_sd_pla, inter_lamb_sd_sca, inter_lamb_sd_tra, 
                  common.legend = T, ncol = 3, nrow = 2, legend = "right", align = "hv")
int_lotte_wunsch

ggsave("C://Users/mandrzej24/Documents/01_revision_demography/03_figures/interaction_boot_sd_full_model_mean_scales.jpeg", 
       dpi = 300, device = "jpeg", width = 4000, height = 2677, units = "px", plot = int_lotte_wunsch)


clim_sd <- ggplot(mean_lambdas %>% subset(treatment %in% c("amb", "fut")), aes(x = climate, y = lambda)) +
  geom_point(position = position_dodge(.08), size = 4, mapping = aes(color = climate)) + 
  facet_wrap(~ species_prep, labeller = label_parsed) +
  theme_bw() +
  geom_errorbar(aes(ymin = lambda - sd, ymax = lambda + sd, color = climate), 
                width = .15, position = position_dodge(.08)) + 
  scale_linetype_manual(values = c(2,1), guide = "none") + 
  xlab("") + ylab(expression(paste("log (", lambda, ")"))) +
  theme(text = element_text(size = 20), strip.text = element_text(size = 15)) + 
  scale_color_manual(name = "Climate", values = rbPalette,
                     label = c("Ambient", "Future")) +
  geom_hline(yintercept = 0, linetype = 2)
clim_sd

manage_sd <- ggplot(mean_lambdas %>% subset(treatment %in% c("gra", "mow")), aes(x = manage, y = lambda)) +
  geom_point(position = position_dodge(.08), size = 4, mapping = aes(shape = manage)) + 
  facet_wrap(~ species_prep, labeller = label_parsed) +
  theme_bw() +
  geom_errorbar(aes(ymin = lambda - sd, ymax = lambda + sd), 
                width = .15, position = position_dodge(.08)) + 
  scale_linetype_manual(values = c(2,1), guide = "none") + 
  xlab("") + ylab(expression(paste("log (", lambda, ")"))) +
  theme(text = element_text(size = 20), strip.text = element_text(size = 15)) + 
  scale_shape_manual(name = "Management", values = c(17, 18),
                     label = c("Grazing", "Mowing")) +
  geom_hline(yintercept = 0, linetype = 2)
manage_sd

ggsave("C://Users/mandrzej24/Documents/01_revision_demography/03_figures/interaction_boot_sd_full_model_mean.jpeg", 
       dpi = 300, device = "jpeg", width = 4000, height = 2677, units = "px", plot = inter_lamb_sd)
ggsave("C://Users/mandrzej24/Documents/01_revision_demography/03_figures/climate_boot_sd_full_model_mean.jpeg", 
       dpi = 300, device = "jpeg", width = 4000, height = 2677, units = "px", plot = clim_sd)
ggsave("C://Users/mandrzej24/Documents/01_revision_demography/03_figures/manage_boot_sd_full_model_mean.jpeg", 
       dpi = 300, device = "jpeg", width = 4000, height = 2677, units = "px", plot = manage_sd)

## statistics, permutation test
uniq_spec <- unique(booted_lambda$species)
uniq_trea <- unique(booted_lambda$treatment)

permutation <- function(species = "Bro_ere", treat1 = "amb_mow", treat2 = "amb_gra", mutations = 1000)
{
  set.seed(5)
  
  spec <- species
  
  first_lambda <- subset(booted_lambda, species == spec & treatment == treat1)
  second_lambda <- subset(booted_lambda, species == spec & treatment == treat2)
  
  n1 <- nrow(first_lambda)
  n2 <- nrow(second_lambda)
  
  lamb_mean1 <- mean(first_lambda$lambda)
  lamb_mean2 <- mean(second_lambda$lambda)
  
  savers <- c()
  for(i in 1:mutations)
  {
    if(lamb_mean2 > lamb_mean1)
      savers[i] <- second_lambda[sample(1:n2, 1), "lambda"] - first_lambda[sample(1:n1, 1), "lambda"]
    
    if(lamb_mean2 < lamb_mean1)
      savers[i] <- first_lambda[sample(1:n1, 1), "lambda"] - second_lambda[sample(1:n2, 1), "lambda"]
  }
  
  perm_out <- length(which(savers <= 0.05))/mutations
  
  return(perm_out)
}


## with permutation nothing is significant
permutation("Bro_ere", "amb_gra", "amb_mow")
permutation("Bro_ere", "amb_gra", "fut_gra")
permutation("Bro_ere", "amb_mow", "fut_mow")
permutation("Bro_ere", "fut_gra", "fut_mow")

permutation("Dia_car", "amb_gra", "amb_mow")
permutation("Dia_car", "amb_gra", "fut_gra")
permutation("Dia_car", "amb_mow", "fut_mow")
permutation("Dia_car", "fut_gra", "fut_mow")

permutation("Pla_lan", "amb_gra", "amb_mow")
permutation("Pla_lan", "amb_gra", "fut_gra")
permutation("Pla_lan", "amb_mow", "fut_mow")
permutation("Pla_lan", "fut_gra", "fut_mow")

permutation("Sca_och", "amb_gra", "amb_mow")
permutation("Sca_och", "amb_gra", "fut_gra")
permutation("Sca_och", "amb_mow", "fut_mow")
permutation("Sca_och", "fut_gra", "fut_mow")

permutation("Tra_ori", "amb_gra", "amb_mow")
permutation("Tra_ori", "amb_gra", "fut_gra")
permutation("Tra_ori", "amb_mow", "fut_mow")
permutation("Tra_ori", "fut_gra", "fut_mow")

permutation("Bro_ere", "amb", "fut")
permutation("Bro_ere", "gra", "mow")

permutation("Dia_car", "amb", "fut")
permutation("Dia_car", "gra", "mow")

permutation("Pla_lan", "amb", "fut")
permutation("Pla_lan", "gra", "mow")

permutation("Sca_och", "amb", "fut")
permutation("Sca_och", "gra", "mow")

permutation("Tra_ori", "amb", "fut")
permutation("Tra_ori", "gra", "mow")
