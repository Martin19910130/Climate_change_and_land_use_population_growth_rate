##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##  Lambda for interaction over the years
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
rm(list = ls())
gc()

library(ggplot2)
library(openxlsx)
library(dplyr)
library(ipmr)
library(parallel)

## color platte
rbPalette <- c("#0072B2", "#D55E00")

## function
inv_logit <- function(int, slope, sv) {
  1/(1 + exp(-(int + slope * sv)))
}

## read data, we don't use Lychnes and Trifolium anymore (not enough data)
dat <- read.xlsx("C://Users/marti/Documents/Martin/01_PhD/02_GCEF_demografie_ipm/Size_seedling_flower_data/Demography_data_2018_2022.xlsx") %>% 
  subset(species %in% c("Bro_ere", "Pla_lan", "Sca_och", "Tra_ori", "Dia_car")) 
dat$spec_year_clim <- paste(dat$species, dat$year, dat$climate, dat$management, sep = "_")

## create column with if flowered 1 and if not 0 
dat$flower_pr <- ifelse(dat$fl_all > 0 , 1, 0)

## load seed data
dat_seeds <- read.csv("C://Users/marti/Documents/Martin/01_PhD/02_GCEF_demografie_ipm/Seed_count_data/Seed_number_transects.csv") %>% 
  subset(year != 2022 & species %in% c("Bro_ere", "Pla_lan", "Sca_och", "Tra_ori", "Dia_car"))
dat_seeds$spec_year_clim <- paste(dat_seeds$species, dat_seeds$year, dat_seeds$climate, dat_seeds$management, sep = "_")

dat$spec_year <- paste(dat$species, dat$year, sep = "_")
dat_seeds$spec_year <- paste(dat_seeds$species, dat_seeds$year, sep = "_")
suber <- unique(dat$spec_year)

boot_lambda <- function(ii)
{
  library(dplyr)
  library(ipmr)
  library(lme4)
  
  inv_logit <- function(int, slope, sv) {
    1/(1 + exp(-(int + slope * sv)))
  }
  
  ## create empty data frame 
  lamb <- data.frame(lambda = NA, species = NA,
                     year = NA, treatment = NA)
  all_lamb <- c()
  
  for(i in 1:length(suber))
  {
    #set.seed(ii)
    
    ## subset data into the different species, climate and years
    dat_mod <- subset(dat, spec_year == suber[i])
    dat_mod_s <- subset(dat_seeds, spec_year == suber[i])

    ## sample seed and normal data and combine both data frames
    dat_mod <- sample_n(dat_mod, nrow(dat_mod), replace = T)
    dat_mod_s <- sample_n(dat_mod_s, nrow(dat_mod_s), replace = T)
    
    if(nrow(dat_mod_s) == 0)
      next
    
    dat_mod_s$unviable_seeds <- as.numeric(dat_mod_s$unviable_seeds)
    
    ## combine the unviable and viable column to one... that's how we did it with Julia
    dat_mod_s$nr_seeds <- rowSums(dat_mod_s[,c ("Seed_count", "unviable", "unviable_seeds")], na.rm = T)
    
    ## calculate the mean seeds per species per year in each treatment
    dat_mod$nr_seeds <- mean(dat_mod_s$nr_seeds, na.rm = T) %>% round()
    
    ## calculate number of seeds
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
    
    if(dat_mod[1, "year"] == 2018)
    {
      
      if(dat_mod[1, "species"] == "Pla_lan")
      {
        dat_2018_seeds <- subset(dat_mod, year == 2018)
        repr_18 <- aggregate(seeds ~ spec_year_clim + plot + subplot, FUN = sum, na.rm = T, data = dat_2018_seeds)
        
        dat_2018 <- subset(dat_mod, year == 2018 & subplot <= 3)
        sl_spring_t0_18 <- unique(dat_2018[, c("plot", "subplot", "spec_year_clim", "sl_spring_t0")]) %>% 
          aggregate(sl_spring_t0 ~ spec_year_clim + plot + subplot, FUN = sum, na.rm = T, data = .)
        
        sl_fall_t0_18 <- unique(dat_2018[, c("plot", "subplot", "spec_year_clim", "sl_fall_t0")]) %>% 
          aggregate(sl_fall_t0 ~ spec_year_clim + plot + subplot, FUN = sum, na.rm = T, data = .)
        
        sl_spring_t1_18 <- unique(dat_2018[, c("plot", "subplot", "spec_year_clim", "sl_spring_t1")]) %>% 
          aggregate(sl_spring_t1 ~ spec_year_clim + plot + subplot, FUN = sum, na.rm = T, data = .)
        
        sl_18 <- merge(sl_spring_t0_18, sl_fall_t0_18) %>% merge(., sl_spring_t1_18)
        
        repr_18 <- merge(repr_18, sl_18)
        
      }
      
      else
      {
        dat_2018 <- subset(dat_mod, year == 2018 & subplot <= 3)
        
        ## number of seedlings for first 3 subplots
        sl_spring_t0_18 <- unique(dat_2018[, c("plot", "subplot", "spec_year_clim", "sl_spring_t0")]) %>% 
          aggregate(sl_spring_t0 ~ spec_year_clim + plot + subplot, FUN = sum, na.rm = T, data = .)
        
        sl_fall_t0_18 <- unique(dat_2018[, c("plot", "subplot", "spec_year_clim", "sl_fall_t0")]) %>% 
          aggregate(sl_fall_t0 ~ spec_year_clim + plot + subplot, FUN = sum, na.rm = T, data = .)
        
        sl_spring_t1_18 <- unique(dat_2018[, c("plot", "subplot", "spec_year_clim", "sl_spring_t1")]) %>% 
          aggregate(sl_spring_t1 ~ spec_year_clim + plot + subplot, FUN = sum, na.rm = T, data = .)
        
        sl_18 <- merge(sl_spring_t0_18, sl_fall_t0_18) %>% merge(., sl_spring_t1_18)
        
        ## number of seeds for first 3 subplots
        repr_18 <- aggregate(seeds ~ spec_year_clim + plot + subplot, FUN = sum, na.rm = T, data = dat_2018) %>% 
          merge(., sl_18)
        
        ## seed to seedling
        repr_18$se_sl_fall_t0 <- (repr_18$sl_fall_t0 / repr_18$seeds) %>% ifelse(. > 1, 1, .)
        repr_18$se_sl_spring_t1 <- (repr_18$sl_spring_t1 /repr_18$seeds) %>% ifelse(. > 1, 1, .)
        
        ## seedling survival
        repr_18 <- aggregate(new_individual ~ spec_year_clim + plot + subplot, FUN = sum, na.rm = T, data = dat_2018) %>% 
          merge(., repr_18)
        
        repr_18$sl_all_t0 <- repr_18$sl_fall_t0 + repr_18$sl_spring_t0
        
        repr_18$new_ind_per_sl <- (repr_18$new_individual/repr_18$sl_all_t0) %>% ifelse(. > 1, 1, .)
        
        repr_dat <- aggregate(repr_18[, c("se_sl_fall_t0", "se_sl_spring_t1", "new_ind_per_sl")], 
                              by = list(repr_18$spec_year_clim), FUN = mean, na.rm = T) %>%
          rename(spec_year_clim = Group.1)
      }
    }
    
    else
    {
      ## think about the way i caclulate the mean maybe do it by plot too
      sl_spring_t0_rest <- unique(dat_mod[, c("plot", "subplot", "spec_year_clim", "sl_spring_t0")]) %>% 
        aggregate(sl_spring_t0 ~ spec_year_clim + plot + subplot, FUN = sum, na.rm = T, data = .)
      
      sl_fall_t0_rest <- unique(dat_mod[, c("plot", "subplot",  "year", "spec_year_clim", "sl_fall_t0")]) %>%
        aggregate(sl_fall_t0 ~ spec_year_clim + plot + subplot, FUN = sum, na.rm = T, data = .)
      
      sl_spring_t1_rest <- unique(dat_mod[, c("plot", "subplot",  "year", "spec_year_clim", "sl_spring_t1")]) %>%
        aggregate(sl_spring_t1 ~ spec_year_clim + plot + subplot, FUN = sum, na.rm = T, data = .)
      
      sl_rest <- merge(sl_spring_t0_rest, sl_fall_t0_rest) %>% merge(., sl_spring_t1_rest)
      
      ## aggregate for number of seeds
      repr_rest <- tryCatch(aggregate(seeds ~ spec_year_clim + plot + subplot, FUN = sum, na.rm = T, data = dat_mod) %>% 
        merge(., sl_rest), 
        error = function(e){cat("ERROR :", conditionMessage(e), "\n")})
      
      ## seed to seedling
      repr_rest$se_sl_fall_t0 <- (repr_rest$sl_fall_t0 / repr_rest$seeds) %>% ifelse(. > 1, 1, .)
      repr_rest$se_sl_spring_t1 <- (repr_rest$sl_spring_t1 / repr_rest$seeds) %>% ifelse(. > 1, 1, .)
      
      ## seedling to new individual
      repr_rest <- tryCatch(aggregate(new_individual ~ spec_year_clim + plot + subplot, FUN = sum, na.rm = T, data = dat_mod) %>% 
        merge(., repr_rest), 
        error = function(e){cat("ERROR :", conditionMessage(e), "\n")})
      
      repr_rest$sl_all_t0 <- repr_rest$sl_spring_t0 + repr_rest$sl_fall_t0
      
      repr_rest$new_ind_per_sl <- (repr_rest$new_individual / repr_rest$sl_all_t0) %>% ifelse(. > 1, 1, .)
      
      ## aggregate for mean
      repr_dat <- tryCatch(aggregate(repr_rest[, c("se_sl_fall_t0", "se_sl_spring_t1", "new_ind_per_sl")],
                            by = list(repr_rest$spec_year_clim), FUN = mean, na.rm = T) %>% 
        rename(spec_year_clim = Group.1), 
        error = function(e){cat("ERROR :", conditionMessage(e), "\n")})
    }
    
    ## calculate the mean of each reproduction parameter and merge
    
    ## get the max sizes and the sd of it for all new individuals
    new_ind_size <- tryCatch(aggregate(log(size_t1) ~ spec_year_clim, FUN = mean, na.rm = T, 
                                       data = dat_mod %>% subset(., new_individual == 1 & log(size_t1) < 2.5)), 
                             error = function(e){cat("ERROR :", conditionMessage(e), "\n")})
    
    new_ind_size$sd <- tryCatch(aggregate(log(size_t1) ~ spec_year_clim, FUN = sd, na.rm = T, data = dat_mod %>% 
                                            subset(., new_individual == 1 & log(size_t1) < 2.5))$'log(size_t1)', 
                                error = function(e){cat("ERROR :", conditionMessage(e), "\n")}) 
    
    new_ind_size$sd[is.na(new_ind_size$sd)] <- 0 
    
    ## get the upper and lower limits
    max_t1 <- tryCatch(aggregate(log(size_t1) ~ spec_year_clim, FUN = max, na.rm = T, data = dat_mod) %>% 
                         rename(size = 'log(size_t1)'), error = function(e){cat("ERROR :", conditionMessage(e), "\n")})
    
    U <- tryCatch(aggregate(log(size_t0) ~ spec_year_clim, FUN = max, na.rm = T, data = dat_mod) %>% rename(size = 'log(size_t0)') %>% 
      rbind(., max_t1) %>% aggregate(size ~ spec_year_clim, FUN = max, data = .), error = function(e){cat("ERROR :", conditionMessage(e), "\n")}) 
    
    min_t1 <- tryCatch(aggregate(log(size_t1) ~ spec_year_clim, FUN = min, na.rm = T, data = dat_mod) %>% 
                         rename(size = 'log(size_t1)'), error = function(e){cat("ERROR :", conditionMessage(e), "\n")}) 
    
    L <- tryCatch(aggregate(log(size_t0) ~ spec_year_clim, FUN = min, na.rm = T, data = dat_mod) %>% rename(size = 'log(size_t0)') %>% 
      rbind(., min_t1) %>% aggregate(size ~ spec_year_clim, FUN = min, data = .), error = function(e){cat("ERROR :", conditionMessage(e), "\n")})
    
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
                                 treatment = names(data_list_full)[j], year = unique(dat_mod$year))
      }
      
    }
    
    all_lamb[i] <- list(lamb)
  }
  
  return(all_lamb)
}

cl <- makeCluster(detectCores() - 1)

clusterExport(cl, c("dat", "dat_seeds", "suber"))

start_time <- Sys.time()
lambdas <- parSapply(cl, 1:1000, FUN = boot_lambda)
Sys.time() - start_time

stopCluster(cl)

booted_lambda <- do.call("rbind", lambdas)

## write table so I dont have to run the code all the time
#write.xlsx(booted_lambda, "C:/Users/mandrzej24/Documents/01_revision_demography/04_data/01_new_data/boot_all_together_mean.xlsx")
booted_lambda <- read.xlsx("C:/Users/mandrzej24/Documents/01_revision_demography/04_data/01_new_data/boot_inter_year_full_mean.xlsx")

colnames(booted_lambda)[3:4] <- c("treatment", "year")

mean_lambdas <- aggregate(log(lambda) ~ species + year + treatment, FUN = mean, data = booted_lambda, na.rm = T) %>% 
  rename(lambda = 'log(lambda)')

mean_lambdas$sd <- aggregate(log(lambda) ~ species + year + treatment, FUN = sd, data = booted_lambda, na.rm = T)$'log(lambda)'

## names for the panel
mean_lambdas$species_prep <- ifelse(mean_lambdas$species == "Bro_ere", "'(a) ' * italic('Bromus erectus')",
                                    ifelse(mean_lambdas$species == "Dia_car", 
                                           "'(b) ' * italic('Dianthus carthusianorum')",
                                           ifelse(mean_lambdas$species == "Pla_lan", 
                                                  "'(c) ' * italic('Plantago lanceolata')",
                                                  ifelse(mean_lambdas$species == "Sca_och", 
                                                         "'(d) ' * italic('Scabiosa ochroleuca')",
                                                         "'(e) ' * italic('Tragopogon orientalis')")))) 

lambda <- ggplot(mean_lambdas %>% subset(!(treatment %in% c("amb", "fut", "gra", "mow"))), aes(x = year, y = lambda, color = treatment, shape = treatment)) + 
  geom_point(position = position_dodge(.5), size = 4) + 
  facet_wrap(~ species_prep, labeller = label_parsed) + 
  geom_errorbar(aes(ymin = lambda - sd, ymax = lambda + sd), width = .15, position = position_dodge(.5)) + 
  scale_color_manual(values = c(rep(rbPalette[1], 2), rep(rbPalette[2], 2)), name = "Treatment", 
                     labels = c("Ambient grazing", "Ambient mowing", "Future grazing", "Future mowing")) + 
  theme_bw() + 
  scale_shape_manual(values = c(17, 18, 17, 18), name = "Treatment", 
                     labels = c("Ambient grazing", "Ambient mowing", "Future grazing", "Future mowing")) + 
  labs(y = expression(paste("log (", lambda, ")")), x = "") + 
  theme(text = element_text(size = 20), axis.text.x = element_text(angle = 90, vjust = 0.5))

lambda

clim_sd <- ggplot(mean_lambdas %>% subset(treatment %in% c("amb", "fut")), aes(x = year, y = lambda, color = treatment)) +
  geom_point(position = position_dodge(.08), size = 4, mapping = aes(color = treatment)) + 
  facet_wrap(~ species_prep, labeller = label_parsed) +
  theme_bw() +
  geom_errorbar(aes(ymin = lambda - sd, ymax = lambda + sd, color = treatment), 
                width = .15, position = position_dodge(.08)) + 
  scale_linetype_manual(values = c(2,1), guide = "none") + 
  xlab("") + ylab(expression(paste("log (", lambda, ")"))) +
  theme(text = element_text(size = 20), strip.text = element_text(size = 15)) + 
  scale_color_manual(name = "Climate", values = rbPalette,
                     label = c("Ambient", "Future")) +
  geom_hline(yintercept = 0, linetype = 2)
clim_sd

manage_sd <- ggplot(mean_lambdas %>% subset(treatment %in% c("gra", "mow")), aes(x = year, y = lambda, shape = treatment)) +
  geom_point(position = position_dodge(.08), size = 4, mapping = aes(shape = treatment)) + 
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

ggsave(plot = lambda, 
      filename = "C:/Users/mandrzej24/Documents/01_revision_demography/03_figures/lambda_year_inter_booted_full_modell.jpeg",
     dpi = 300, width = 3500, height = 2342, units = "px")
ggsave(plot = clim_sd, 
       filename = "C:/Users/mandrzej24/Documents/01_revision_demography/03_figures/lambda_year_climate_booted_full_modell.jpeg",
       dpi = 300, width = 3500, height = 2342, units = "px")
ggsave(plot = manage_sd, 
       filename = "C:/Users/mandrzej24/Documents/01_revision_demography/03_figures/lambda_year_manage_booted_full_modell.jpeg",
       dpi = 300, width = 3500, height = 2342, units = "px")

