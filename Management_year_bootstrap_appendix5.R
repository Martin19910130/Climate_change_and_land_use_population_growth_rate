##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##  Use only the landuse data over the year
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
rm(list = ls())
gc()

## if needed use this to install packages
install.packages(c("ggplot2", "openxlsx", "dplyr", "ipmr", "parallel", "tidyverse"))

library(ggplot2)
library(openxlsx)
library(dplyr)
library(ipmr)
library(parallel)
library(tidyverse)

## function
inv_logit <- function(int, slope, sv) {
  1/(1 + exp(-(int + slope * sv)))
}

## read data, we don't use Lychnes and Trifolium anymore (not enough data)
dat <- read.xlsx("C://Users/marti/Documents/Martin/01_PhD/02_GCEF_demografie_ipm/Size_seedling_flower_data/Demography_data_2018_2022.xlsx") %>% 
  subset(species %in% c("Bro_ere", "Dia_car", "Pla_lan", "Sca_och", "Tra_ori"))
dat$spec_year_mana <- paste(dat$species, dat$year, dat$management, sep = "_")

## create column with if flowered 1 and if not 0 
dat$flower_pr <- ifelse(dat$fl_all > 0 , 1, 0)

## load seed data
dat_seeds <- read.csv("C://Users/marti/Documents/Martin/01_PhD/02_GCEF_demografie_ipm/Seed_count_data/Seed_number_transects.csv") %>% 
  subset(year != 2022 & species %in% c("Bro_ere", "Dia_car", "Pla_lan", "Sca_och", "Tra_ori")) 
dat_seeds$spec_year_mana <- paste(dat_seeds$species, dat_seeds$year, dat_seeds$management, sep = "_")
suber <- unique(dat$spec_year_mana)

aggregate(Seed_count ~ year + species + management, FUN = mean, data = dat_seeds)

for(i in 1:length(suber))
  subset(dat_seeds, spec_year_mana == suber[i]) %>% print()

boot_lambda <- function(ii)
{
  library(dplyr)
  library(ipmr)
  
  inv_logit <- function(int, slope, sv) {
    1/(1 + exp(-(int + slope * sv)))
  }
  
  ## create empty data frame 
  lamb <- data.frame(lambda = NA, species = NA,
                     year = NA, management = NA)
  
  for(i in 1:length(suber))
  {
    
    set.seed(ii)
    
    ## subset data into the different species, manaate and years
    dat_mod <- subset(dat, spec_year_mana == suber[i])
    dat_mod_s <- subset(dat_seeds, spec_year_mana == suber[i])
    
    ## sample seed and normal data and combine both data frames
    dat_mod <- sample_n(dat_mod, nrow(dat_mod), replace = T)
    dat_mod_s <- sample_n(dat_mod_s, nrow(dat_mod_s), replace = T)
    
    dat_mod_s$unviable_seeds <- as.numeric(dat_mod_s$unviable_seeds)
    
    ## combine the unviable and viable column to one... that's how we did it with Julia
    dat_mod_s$nr_seeds <- rowSums(dat_mod_s[,c ("Seed_count", "unviable", "unviable_seeds")], na.rm = T)
    
    ## calculate the mean seeds per species per year in each treatment
    mean_seed <- aggregate(nr_seeds ~ year + species + management, FUN = mean, na.rm = T, data = dat_mod_s)
    
    ## merge the seed data to the entire data frame by species and manaate
    dat_mod <- merge(dat_mod, mean_seed[,c("species",  "year", "management", "nr_seeds")], 
                     by = c("species", "year", "management"))
    
    ## calculate number of seeds
    dat_mod$seeds <- dat_mod$fl_all * dat_mod$nr_seeds 
    
    surv_mods <- tryCatch(glm(survival ~ log(size_t0), data = dat_mod, family = "binomial"), 
                          error = function(e){cat("ERROR :", conditionMessage(e), "\n")}) %>% summary() 
    grow_mods <- tryCatch(lm(log(size_t1) ~ log(size_t0), data = dat_mod),
                          error = function(e){cat("ERROR :", conditionMessage(e), "\n")}) %>% summary() 
    
    flow_mods <- tryCatch(glm(flower_pr ~ log(size_t0), data = dat_mod, family = "binomial"),
                          error = function(e){cat("ERROR :", conditionMessage(e), "\n")}) %>% summary() 
    
    dat_seed <- subset(dat_mod, flower_pr == 1)
    seed_mods <- tryCatch(glm(seeds ~ log(size_t0), data = dat_seed, family = "poisson"),
                          error = function(e){cat("ERROR :", conditionMessage(e), "\n")}) %>% summary()
    
    if(dat_mod[1, "year"] == 2018)
    {
      
      if(dat_mod[1, "species"] == "Pla_lan")
      {
        dat_2018_seeds <- subset(dat_mod, year == 2018)
        repr_18 <- aggregate(seeds ~ spec_year_mana + plot + subplot, FUN = sum, na.rm = T, data = dat_2018_seeds) %>% 
          merge(., sl_18)
        
        dat_2018 <- subset(dat_mod, year == 2018 & subplot <= 3)
        sl_spring_t0_18 <- unique(dat_2018[, c("plot", "subplot", "spec_year_mana", "sl_spring_t0")]) %>% 
          aggregate(sl_spring_t0 ~ spec_year_mana + plot + subplot, FUN = sum, na.rm = T, data = .)
        
        sl_fall_t0_18 <- unique(dat_2018[, c("plot", "subplot", "spec_year_mana", "sl_fall_t0")]) %>% 
          aggregate(sl_fall_t0 ~ spec_year_mana + plot + subplot, FUN = sum, na.rm = T, data = .)
        
        sl_spring_t1_18 <- unique(dat_2018[, c("plot", "subplot", "spec_year_mana", "sl_spring_t1")]) %>% 
          aggregate(sl_spring_t1 ~ spec_year_mana + plot + subplot, FUN = sum, na.rm = T, data = .)
        
        sl_18 <- merge(sl_spring_t0_18, sl_fall_t0_18) %>% merge(., sl_spring_t1_18)
        
        repr_18 <- merge(repr_18, sl_18)
        
      }
      
      else
      {
        dat_2018 <- subset(dat_mod, year == 2018 & subplot <= 3)
        
        ## number of seedlings for first 3 subplots
        sl_spring_t0_18 <- unique(dat_2018[, c("plot", "subplot", "spec_year_mana", "sl_spring_t0")]) %>% 
          aggregate(sl_spring_t0 ~ spec_year_mana + plot + subplot, FUN = sum, na.rm = T, data = .)
        
        sl_fall_t0_18 <- unique(dat_2018[, c("plot", "subplot", "spec_year_mana", "sl_fall_t0")]) %>% 
          aggregate(sl_fall_t0 ~ spec_year_mana + plot + subplot, FUN = sum, na.rm = T, data = .)
        
        sl_spring_t1_18 <- unique(dat_2018[, c("plot", "subplot", "spec_year_mana", "sl_spring_t1")]) %>% 
          aggregate(sl_spring_t1 ~ spec_year_mana + plot + subplot, FUN = sum, na.rm = T, data = .)
        
        sl_18 <- merge(sl_spring_t0_18, sl_fall_t0_18) %>% merge(., sl_spring_t1_18)
        
        ## number of seeds for first 3 subplots
        repr_18 <- aggregate(seeds ~ spec_year_mana + plot + subplot, FUN = sum, na.rm = T, data = dat_2018) %>% 
          merge(., sl_18)
        
        ## seed to seedling
        repr_18$se_sl_fall_t0 <- (repr_18$sl_fall_t0 / repr_18$seeds) %>% ifelse(. > 1, 1, .)
        repr_18$se_sl_spring_t1 <- (repr_18$sl_spring_t1 /repr_18$seeds) %>% ifelse(. > 1, 1, .)
        
        ## seedling survival
        repr_18 <- aggregate(new_individual ~ spec_year_mana + plot + subplot, FUN = sum, na.rm = T, data = dat_2018) %>% 
          merge(., repr_18)
        
        repr_18$sl_all_t0 <- repr_18$sl_fall_t0 + repr_18$sl_spring_t0
        
        repr_18$new_ind_per_sl <- (repr_18$new_individual/repr_18$sl_all_t0) %>% ifelse(. > 1, 1, .)
        
        repr_dat <- aggregate(repr_18[, c("se_sl_fall_t0", "se_sl_spring_t1", "new_ind_per_sl")], 
                              by = list(repr_18$spec_year_mana), FUN = mean, na.rm = T) %>%
          rename(spec_year_mana = Group.1)
      }
    }
    
    else
    {
      ## think about the way i caclulate the mean maybe do it by plot too
      sl_spring_t0_rest <- unique(dat_mod[, c("plot", "subplot", "spec_year_mana", "sl_spring_t0")]) %>% 
        aggregate(sl_spring_t0 ~ spec_year_mana + plot + subplot, FUN = sum, na.rm = T, data = .)
      
      sl_fall_t0_rest <- unique(dat_mod[, c("plot", "subplot",  "year", "spec_year_mana", "sl_fall_t0")]) %>%
        aggregate(sl_fall_t0 ~ spec_year_mana + plot + subplot, FUN = sum, na.rm = T, data = .)
      
      sl_spring_t1_rest <- unique(dat_mod[, c("plot", "subplot",  "year", "spec_year_mana", "sl_spring_t1")]) %>%
        aggregate(sl_spring_t1 ~ spec_year_mana + plot + subplot, FUN = sum, na.rm = T, data = .)
      
      sl_rest <- merge(sl_spring_t0_rest, sl_fall_t0_rest) %>% merge(., sl_spring_t1_rest)
      
      ## aggregate for number of seeds
      repr_rest <- aggregate(seeds ~ spec_year_mana + plot + subplot, FUN = sum, na.rm = T, data = dat_mod) %>% 
        merge(., sl_rest)
      
      ## seed to seedling
      repr_rest$se_sl_fall_t0 <- (repr_rest$sl_fall_t0 / repr_rest$seeds) %>% ifelse(. > 1, 1, .)
      repr_rest$se_sl_spring_t1 <- (repr_rest$sl_spring_t1 / repr_rest$seeds) %>% ifelse(. > 1, 1, .)
      
      ## seedling to new individual
      repr_rest <- aggregate(new_individual ~ spec_year_mana + plot + subplot, FUN = sum, na.rm = T, data = dat_mod) %>% 
        merge(., repr_rest)
      
      repr_rest$sl_all_t0 <- repr_rest$sl_spring_t0 + repr_rest$sl_fall_t0
      
      repr_rest$new_ind_per_sl <- (repr_rest$new_individual / repr_rest$sl_all_t0) %>% ifelse(. > 1, 1, .)
      
      ## aggregate for mean
      repr_dat <- aggregate(repr_rest[, c("se_sl_fall_t0", "se_sl_spring_t1", "new_ind_per_sl")],
                            by = list(repr_rest$spec_year_mana), FUN = mean, na.rm = T) %>% 
        rename(spec_year_mana = Group.1)
    }
    
    ## calculate the mean of each reproduction parameter and merge
    
    ## get the max sizes and the sd of it for all new individuals
    new_ind_size <- tryCatch(aggregate(log(size_t1) ~ spec_year_mana, FUN = mean, na.rm = T, 
                                       data = dat_mod %>% subset(., new_individual == 1 & log(size_t1) < 2.5)), 
                             error = function(e){cat("ERROR :", conditionMessage(e), "\n")})
    
    new_ind_size$sd <- tryCatch(aggregate(log(size_t1) ~ spec_year_mana, FUN = sd, na.rm = T, data = dat_mod %>% 
                                            subset(., new_individual == 1 & log(size_t1) < 2.5))$'log(size_t1)', 
                                error = function(e){cat("ERROR :", conditionMessage(e), "\n")}) 
    
    new_ind_size$sd[is.na(new_ind_size$sd)] <- 0 
    
    ## get the upper and lower limits
    max_t1 <- tryCatch(aggregate(log(size_t1) ~ spec_year_mana, FUN = max, na.rm = T, data = dat_mod) %>% 
                         rename(size = 'log(size_t1)'), error = function(e){cat("ERROR :", conditionMessage(e), "\n")})
    
    U <- aggregate(log(size_t0) ~ spec_year_mana, FUN = max, na.rm = T, data = dat_mod) %>% rename(size = 'log(size_t0)') %>% 
      rbind(., max_t1) %>% aggregate(size ~ spec_year_mana, FUN = max, data = .)
    
    min_t1 <- tryCatch(aggregate(log(size_t1) ~ spec_year_mana, FUN = min, na.rm = T, data = dat_mod) %>% 
                         rename(size = 'log(size_t1)'), error = function(e){cat("ERROR :", conditionMessage(e), "\n")}) 
    
    L <- aggregate(log(size_t0) ~ spec_year_mana, FUN = min, na.rm = T, data = dat_mod) %>% rename(size = 'log(size_t0)') %>% 
      rbind(., min_t1) %>% aggregate(size ~ spec_year_mana, FUN = min, data = .)
    
    data_list <- tryCatch(list(s_int = coef(surv_mods)[1], ## survival intercept
                               s_slope = coef(surv_mods)[2], ## survival slope
                               
                               g_int = coef(grow_mods)[1], ##  growth intercept
                               g_slope = coef(grow_mods)[2], ## growth slope
                               g_sd = grow_mods$sigma, ## growth sd 
                               
                               r_r_int = coef(flow_mods)[1], ## flower probability intercept
                               r_r_slope = coef(flow_mods)[2], ## flower probability slope
                               
                               r_s_int = coef(seed_mods)[1], ## number of seeds intercept
                               r_s_slope = coef(seed_mods)[2], ## number of seeds slope
                               
                               r_d_mu = new_ind_size[, "log(size_t1)"], ## size distribution of new plants
                               r_d_sd = new_ind_size[, "sd"], ## sd of size distribution of new plants
                               
                               g_i1 = repr_dat[, "se_sl_fall_t0"], ## number of seeds that turn to seedling in fall
                               g_i2 = repr_dat[, "se_sl_spring_t1"], ## number of seeds that turn to seedling in spring
                               
                               e_p = repr_dat[, "new_ind_per_sl"], ## number of new individual that emerge from 
                               
                               L = L[, "size"], ## lower limit size (no worries the size is logged)
                               U = U[, "size"], ## upper limit size (no worries the size is logged)
                               
                               n = 200), error = function(e){cat("ERROR :", conditionMessage(e), "\n")})
    
    
    if(is.null(data_list))
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
        lamb[i,] <- data.frame(lambda = lambda(general_ipm), species = unique(dat_mod$species),
                               year = unique(dat_mod$year), management = unique(dat_mod$management))
    }
  }
  
  return(lamb)
}

#start_time <- Sys.time()
#lambdas <- sapply(1:1000, FUN = boot_lambda)
#Sys.time() - start_time

cl <- makeCluster(detectCores() - 1)

clusterExport(cl, c("dat", "dat_seeds", "suber"))

start_time <- Sys.time()
lambdas <- parSapply(cl, 1:1000, FUN = boot_lambda)
Sys.time() - start_time

stopCluster(cl)

lambdas <- t(lambdas)

lamb_list <- c()
for(i in 1:nrow(lambdas))
  lamb_list[[i]] <- data.frame(lambda = lambdas[i, "lambda"], species = lambdas[i, "species"], 
                               management = lambdas[i, "management"], year = lambdas[i, "year"], iteration = i)

booted_lambda <- do.call("rbind", lamb_list) %>% subset(!is.na(species))

## save the result so I dont have to load it again
#write.xlsx(booted_lambda, 
#         "C:/Users/marti/Documents/Martin/01_PhD/02_GCEF_demografie_ipm/lambda_data/booted_original_manag_year.xlsx")

booted_lambda <- read.xlsx("C:/Users/marti/Documents/Martin/01_PhD/02_GCEF_demografie_ipm/lambda_data/booted_original_manag_year.xlsx")

mean_lambdas <- aggregate(log(lambda) ~ species + year + management, FUN = mean, data = booted_lambda, na.rm = T) %>% rename(lambda = 'log(lambda)')
mean_lambdas$ci025 <- aggregate(log(lambda) ~ species + year + management, FUN = quantile, 
                                data = booted_lambda, na.rm = T, probs = c(0.025, 0.975))$'log(lambda)' %>% .[,1]
mean_lambdas$ci975 <- aggregate(log(lambda) ~ species + year + management, FUN = quantile, 
                                data = booted_lambda, na.rm = T, probs = c(0.025, 0.975))$'log(lambda)' %>% .[,2] 

## names for the panel
mean_lambdas$species_prep <- ifelse(mean_lambdas$species == "Bro_ere", "'(A) ' * italic('Bromus erectus')",
                                    ifelse(mean_lambdas$species == "Dia_car", 
                                           "'(B) ' * italic('Dianthus carthusianorum')",
                                           ifelse(mean_lambdas$species == "Pla_lan", 
                                                  "'(C) ' * italic('Plantago lanceolata')",
                                                  ifelse(mean_lambdas$species == "Sca_och", 
                                                         "'(D) ' * italic('Scabiosa ochroleuca')",
                                                         "'(E) ' * italic('Tragopogon orientalis')")))) 

lambda <- ggplot(mean_lambdas, aes(x = year, y = lambda, shape = management)) +
  geom_point(position = position_dodge(.5), size =4) + 
  facet_wrap(~ species_prep, labeller = label_parsed) + 
  scale_shape_manual(name = "Management", values = c(18, 17)) +
  geom_errorbar(aes(ymin = 0 + ci025, ymax = 0 + ci975), width = .15, 
                position = position_dodge(.5)) + theme_bw() +
  labs(y = expression(paste("log (", lambda, ")")), 
       x = "") + 
  theme(text = element_text(size = 20))

lambda  

ggsave(plot = lambda, 
       file = "C:/Users/marti/Documents/Martin/01_PhD/02_GCEF_demografie_ipm/Paper_figures/Appendix/lambda_year_manag_boot.jpeg",
       device = "jpeg", dpi = 300, width = 3500, height = 2342, units = "px")
