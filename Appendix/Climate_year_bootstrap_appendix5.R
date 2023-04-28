##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##        IPM only climate including years
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
rm(list = ls())
gc()

library(ggplot2)
library(dplyr)
library(openxlsx)
library(ipmr)
library(parallel)

## color platte
rbPalette <- c("#0072B2", "#D55E00")

## function
inv_logit <- function(int, slope, sv) {
  1/(1 + exp(-(int + slope * sv)))
}

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##        Read and prepare data
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## read data
dat <- read.xlsx("C://Users/ma22buky/Documents/01_PhD/02_GCEF_demografie_ipm/Size_seedling_flower_data/transect_data_rearranged.xlsx") %>% 
  subset(species %in% c("Bro_ere", "Dia_car", "Pla_lan", "Sca_och", "Tra_ori"))
dat$spec_clim_year <- paste(dat$species, dat$climate, dat$year, sep = "_")

dat$flower_pr <- ifelse(dat$fl_all > 0, 1, 0)
## load seed data
dat_seeds <- read.csv("C://Users/ma22buky/Documents/01_PhD/02_GCEF_demografie_ipm/Seed_count_data/Seed_table_18_19_20_transect.csv") %>% 
  subset(., year != 2022) 
dat_seeds$unviable_seeds <- as.numeric(dat_seeds$unviable_seeds)
dat_seeds$spec_clim_year <- paste(dat_seeds$species, dat_seeds$climate, dat_seeds$year, sep = "_")

## combine the unviable and viable column to one... thats how we did it with julia
dat_seeds$nr_seeds <- rowSums(dat_seeds[,c ("Seed_count", "unviable", "unviable_seeds")], na.rm = T)

suber <- unique(dat$spec_clim_year)

boot_lambda <- function(ii)
{
  library(dplyr)
  library(ipmr)
  
  inv_logit <- function(int, slope, sv) {
    1/(1 + exp(-(int + slope * sv)))
  }
  
  ## create empty data frame 
  lamb <- data.frame(lambda = NA, species = NA, climate = NA, year = NA)
  
  for(i in 1:length(suber))
  {
    
    set.seed(ii)
    
    ## subset data into the different species, climate and years
    dat_mod <- subset(dat, spec_clim_year == suber[i])
    dat_mod_s <- subset(dat_seeds, spec_clim_year == suber[i])
    
    ## sample seed and normal data and combine both data frames
    dat_mod <- sample_n(dat_mod, nrow(dat_mod), replace = T)
    dat_mod_s <- sample_n(dat_mod_s, nrow(dat_mod_s), replace = T)
    
    dat_mod_s$unviable_seeds <- as.numeric(dat_mod_s$unviable_seeds)
    
    ## combine the unviable and viable column to one... that's how we did it with Julia
    dat_mod_s$nr_seeds <- rowSums(dat_mod_s[,c ("Seed_count", "unviable", "unviable_seeds")], na.rm = T)
    
    ## calculate the mean seeds per species per year in each treatment
    mean_seed <- aggregate(nr_seeds ~ species + climate + year, FUN = mean, na.rm = T, data = dat_mod_s)
    
    ## merge the seed data to the entire data frame by species and climate
    dat_mod <- merge(dat_mod, mean_seed[,c("species", "climate", "year", "nr_seeds")], 
                     by = c("species", "climate", "year"))
    
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
    
    ## the reproduction for 2018 is still different due to the sampeling, problems with Plantago
    ## temporal solution for plantago...
    if(dat_mod[i, "year"] == 2018)
    {
      if(suber[i] == "Pla_lan_ambient" | suber[i] == "Pla_lan_future")
        dat_2018 <- subset(dat_mod, year == 2018)
      
      else
        dat_2018 <- subset(dat_mod, year == 2018 & subplot <= 3)
      
      ## number of seedlings for first 3 subplots
      sl_spring_t0_18 <- unique(dat_2018[, c("plot", "subplot", "spec_clim_year", "sl_spring_t0")]) %>% 
        aggregate(sl_spring_t0 ~ spec_clim_year + plot + subplot, FUN = sum, na.rm = T, data = .)
      
      sl_fall_t0_18 <- unique(dat_2018[, c("plot", "subplot", "spec_clim_year", "sl_fall_t0")]) %>% 
        aggregate(sl_fall_t0 ~ spec_clim_year + plot + subplot, FUN = sum, na.rm = T, data = .)
      
      sl_spring_t1_18 <- unique(dat_2018[, c("plot", "subplot", "spec_clim_year", "sl_spring_t1")]) %>% 
        aggregate(sl_spring_t1 ~ spec_clim_year + plot + subplot, FUN = sum, na.rm = T, data = .)
      
      sl_18 <- merge(sl_spring_t0_18, sl_fall_t0_18) %>% merge(., sl_spring_t1_18)
      
      ## number of seeds for first 3 subplots
      repr_18 <- aggregate(seeds ~ spec_clim_year + plot + subplot, FUN = sum, na.rm = T, data = dat_2018) %>% 
        merge(., sl_18)
      
      ## seed to seedling
      repr_18$se_sl_fall_t0 <- (repr_18$sl_fall_t0 / repr_18$seeds) %>% ifelse(. > 1, 1, .)
      repr_18$se_sl_spring_t1 <- (repr_18$sl_spring_t1 /repr_18$seeds) %>% ifelse(. > 1, 1, .)
      
      ## seedling survival
      repr_18 <- aggregate(new_individual ~ spec_clim_year + plot + subplot, FUN = sum, na.rm = T, data = dat_2018) %>% 
        merge(., repr_18)
      
      repr_18$sl_all_t0 <- repr_18$sl_fall_t0 + repr_18$sl_spring_t0
      
      repr_18$new_ind_per_sl <- (repr_18$new_individual/repr_18$sl_all_t0) %>% ifelse(. > 1, 1, .)
      
      repr_dat <- aggregate(repr_18[, c("se_sl_fall_t0", "se_sl_spring_t1", "new_ind_per_sl")], 
                            by = list(repr_18$spec_clim_year), FUN = mean, na.rm = T) %>%
        rename(spec_clim_year = Group.1)
    }
    else
    {
      ## think about the way i caclulate the mean maybe do it by plot too
      sl_spring_t0_rest <- unique(dat_mod[, c("plot", "subplot", "spec_clim_year", "sl_spring_t0")]) %>% 
        aggregate(sl_spring_t0 ~ spec_clim_year + plot + subplot, FUN = sum, na.rm = T, data = .)
      
      sl_fall_t0_rest <- unique(dat_mod[, c("plot", "subplot",  "year", "spec_clim_year", "sl_fall_t0")]) %>%
        aggregate(sl_fall_t0 ~ spec_clim_year + plot + subplot, FUN = sum, na.rm = T, data = .)
      
      sl_spring_t1_rest <- unique(dat_mod[, c("plot", "subplot",  "year", "spec_clim_year", "sl_spring_t1")]) %>%
        aggregate(sl_spring_t1 ~ spec_clim_year + plot + subplot, FUN = sum, na.rm = T, data = .)
      
      sl_rest <- merge(sl_spring_t0_rest, sl_fall_t0_rest) %>% merge(., sl_spring_t1_rest)
      
      ## aggregate for number of seeds
      repr_rest <- aggregate(seeds ~ spec_clim_year + plot + subplot, FUN = sum, na.rm = T, data = dat_mod) %>% 
        merge(., sl_rest)
      
      ## seed to seedling
      repr_rest$se_sl_fall_t0 <- (repr_rest$sl_fall_t0 / repr_rest$seeds) %>% ifelse(. > 1, 1, .)
      repr_rest$se_sl_spring_t1 <- (repr_rest$sl_spring_t1 / repr_rest$seeds) %>% ifelse(. > 1, 1, .)
      
      ## seedling to new individual
      repr_rest <- aggregate(new_individual ~ spec_clim_year + plot + subplot, FUN = sum, na.rm = T, data = dat_mod) %>% 
        merge(., repr_rest)
      
      repr_rest$sl_all_t0 <- repr_rest$sl_spring_t0 + repr_rest$sl_fall_t0
      
      repr_rest$new_ind_per_sl <- (repr_rest$new_individual / repr_rest$sl_all_t0) %>% ifelse(. > 1, 1, .)
      
      ## aggregate for mean
      repr_dat <- aggregate(repr_rest[, c("se_sl_fall_t0", "se_sl_spring_t1", "new_ind_per_sl")],
                            by = list(repr_rest$spec_clim_year), FUN = mean, na.rm = T) %>% 
        rename(spec_clim_year = Group.1)
    }
    ## calculate the mean of each reproduction parameter and merge
    
    ## get the max sizes and the sd of it for all new individuals
    new_ind_size <- tryCatch(aggregate(log(size_t1) ~ spec_clim_year, FUN = mean, na.rm = T, 
                                       data = dat_mod %>% subset(., new_individual == 1 & log(size_t1) < 2.5)), 
                             error = function(e){cat("ERROR :", conditionMessage(e), "\n")})
    
    new_ind_size$sd <- tryCatch(aggregate(log(size_t1) ~ spec_clim_year, FUN = sd, na.rm = T, data = dat_mod %>% 
                                            subset(., new_individual == 1 & log(size_t1) < 2.5))$'log(size_t1)', 
                                error = function(e){cat("ERROR :", conditionMessage(e), "\n")}) 
    
    new_ind_size$sd[is.na(new_ind_size$sd)] <- 0 
    
    ## get the upper and lower limits
    max_t1 <- tryCatch(aggregate(log(size_t1) ~ spec_clim_year, FUN = max, na.rm = T, data = dat_mod) %>% 
                         rename(size = 'log(size_t1)'), error = function(e){cat("ERROR :", conditionMessage(e), "\n")})
    
    U <- aggregate(log(size_t0) ~ spec_clim_year, FUN = max, na.rm = T, data = dat_mod) %>% rename(size = 'log(size_t0)') %>% 
      rbind(., max_t1) %>% aggregate(size ~ spec_clim_year, FUN = max, data = .)
    
    min_t1 <- tryCatch(aggregate(log(size_t1) ~ spec_clim_year, FUN = min, na.rm = T, data = dat_mod) %>% 
                         rename(size = 'log(size_t1)'), error = function(e){cat("ERROR :", conditionMessage(e), "\n")}) 
    
    L <- aggregate(log(size_t0) ~ spec_clim_year, FUN = min, na.rm = T, data = dat_mod) %>% rename(size = 'log(size_t0)') %>% 
      rbind(., min_t1) %>% aggregate(size ~ spec_clim_year, FUN = min, data = .)
    
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
        lamb[i,] <- data.frame(lambda = lambda(general_ipm), species = unique(dat_mod$species), climate = unique(dat_mod$climate), 
                               year = unique(dat_mod$year))
    }
  }
  
  return(lamb)
}

cl <- makeCluster(detectCores() - 1)

clusterExport(cl, c("dat", "dat_seeds", "suber"))

start_time <- Sys.time()
lambdas <- parSapply(cl ,1:1000, FUN = boot_lambda)
Sys.time() - start_time

stopCluster(cl)

lambdas <- t(lambdas)

unlist(lambdas)

lamb_list <- c()
for(i in 1:nrow(lambdas))
  lamb_list[[i]] <- data.frame(lambda = lambdas[i, "lambda"], species = lambdas[i, "species"], 
                               climate = lambdas[i, "climate"], year = lambdas[i, "year"], iteration = i)

booted_lambda <- do.call("rbind", lamb_list) %>% subset(!is.na(species))
booted_lambda[which(booted_lambda$species == "Bro_ere" & booted_lambda$year == 2018 & booted_lambda$climate == "ambient"),] 

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## run code from here if I didn't change anything in the data or code
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## write the original bootstrapped lambdas as xlsx
#write.xlsx(booted_lambda, "C://Users/marti/Documents/Martin/01_PhD/02_GCEF_demografie_ipm/data/booted_original_clim_year.xlsx")
booted_lambda <- read.xlsx("C://Users/marti/Documents/Martin/01_PhD/02_GCEF_demografie_ipm/data/booted_original_clim_year.xlsx")

mean_lambdas <- aggregate(log(lambda) ~ species + climate + year, FUN = mean, data = booted_lambda, na.rm = T) %>% rename(lambda = 'log(lambda)')
mean_lambdas$ci025 <- aggregate(log(lambda) ~ species + climate + year, FUN = quantile, 
                                data = booted_lambda, na.rm = T, probs = c(0.025, 0.975))$'log(lambda)' %>% .[,1]
mean_lambdas$ci975 <- aggregate(log(lambda) ~ species + climate + year, FUN = quantile, 
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

lambdas <- ggplot(mean_lambdas, aes(x = year, y = lambda, color = climate)) +
  geom_point(position = position_dodge(.5), size = 4) + 
  facet_wrap(~ species_prep, labeller = label_parsed) + 
  geom_errorbar(aes(ymin = 0 + ci025, ymax = 0 + ci975), width = .15, position = position_dodge(.5)) + 
  scale_color_manual(values = rbPalette, name = "Climate") +
  theme_bw() + ylab(expression(paste("log (", lambda, ")"))) + xlab("") +
  theme(text = element_text(size = 20))

lambdas

ggsave(plot = lambdas, 
       filename = "C:/Users/marti/Documents/Martin/01_PhD/02_GCEF_demografie_ipm/Paper_figures/Appendix/lambda_year_clim_boot.jpeg",
       device = "jpeg",
       dpi = 300, height = 2342, width = 3500, units = "px")

permutation <- function(species = "Bro_ere", treat1 = "ambient", treat2 = "future", y = 2018, mutations = 1000)
{
  set.seed(1)
  spec <- species
  treat1 <- treat1
  treat2 <- treat2
  y <- y
  
  first_lambda <- subset(booted_lambda, species == spec & climate == treat1 & year == y)
  second_lambda <- subset(booted_lambda, species == spec & climate == treat2 & year == y)
  
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
permutation(species = "Bro_ere", y = 2019)
