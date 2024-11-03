##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##        IPM only managment without years
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
rm(list = ls())
gc()

library(ggplot2)
library(dplyr)
library(openxlsx)
library(ipmr)
library(lme4)
library(rcartocolor)

## color platte
rbPalette <- c("#0072B2", "#D55E00")

colorbl <- carto_pal(2,"Earth")[c(2,1)]
## function
inv_logit <- function(int, slope, sv) {
  1/(1 + exp(-(int + slope * sv)))
}

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##        Read and prepare data
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## read data
dat <- read.xlsx("C://Users/mandrzej24/Documents/01_revision_demography/04_data/03_original_data/Demography_data_2018_2022.xlsx") %>% 
  subset(species %in% c("Bro_ere", "Dia_car", "Pla_lan", "Sca_och", "Tra_ori"))

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
  subset(., year != 2022) 
dat_seeds$unviable_seeds <- as.numeric(dat_seeds$unviable_seeds)

## combine the unviable and viable column to one... thats how we did it with julia
dat_seeds$nr_seeds <- rowSums(dat_seeds[,c ("Seed_count", "unviable", "unviable_seeds")], na.rm = T)

## calculate the mean seeds per species per year in each treatment
mean_seed <- aggregate(nr_seeds ~ year + species + management + climate, FUN = mean, na.rm = T, data = dat_seeds)
mean_seed$treatment <- paste(mean_seed$climate, mean_seed$management, sep = "_")

## merge the seed data to the entire data frame
dat_comp <- merge(dat, mean_seed)

## caclulate how many seeds each individual produced at any given time (seeds * flower_all)
dat_comp$seeds <- round(dat_comp$fl_all * dat_comp$nr_seeds)

## to a column which tells us if a ind had flowered 
dat_comp$flower_pr <- ifelse(dat_comp$fl_all >= 1, 1, 0)

## run the 4 easy models
suber <- unique(dat_comp$species)

surv_mods <- c()
grow_mods <- c()
flow_mods <- c()
seed_mods <- c()
lamb <- data.frame(lambda = NA, species = NA, treatment = NA)
lamb_all <- c()

for(i in 1:length(suber))
{
  #set.seed(ii)
  
  ## subset data into the different species, climate and years
  dat_mod <- subset(dat_comp, species == suber[i])
  
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
  lamb_all[[lamb[1, "species"]]] <- lamb
}

lamb_mana <- do.call("rbind", lamb_all) %>% subset(treatment %in% c("gra", "mow"))
lamb_clim <- do.call("rbind", lamb_all) %>% subset(treatment %in% c("amb", "fut"))

colnames(lamb_mana) <- c("lamb", "spec", "manag")
colnames(lamb_clim) <- c("lamb", "spec", "clima")

lamb_mana$manag <- ifelse(lamb_mana$manag == "gra", "grazing", "mowing")
lamb_clim$clima <- ifelse(lamb_clim$clima == "amb", "ambient", "future")

ggplot(lamb_mana, aes(x = spec, y = log(lamb), shape = manag)) + geom_point(position = position_dodge(.25), size = 2) + theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1.1)) + geom_hline(yintercept = 0, linetype = 3) + 
  scale_shape_manual(values = c(17, 18), name = "Management", labels = c("Grazing", "Mowing")) + ylab(expression(paste("log (", lambda, ")"))) + xlab("")

ggplot(lamb_clim, aes(x = spec, y = log(lamb), color = clima)) + geom_point(position = position_dodge(.25), size = 2) + theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1.1)) + geom_hline(yintercept = 0, linetype = 3) + 
  scale_shape_manual(values = c(17, 18), name = "Management") + ylab(expression(paste("log (", lambda, ")"))) + xlab("") + 
  scale_color_manual(values = rbPalette, labels = c("Ambient", "Future"), name = "Climate")

## load Caros flowering stuff
caro_fl <- read.csv("C://Users/mandrzej24/Documents/01_revision_demography/04_data/03_original_data/flowering times 2020_CPlos_updated.csv", sep = ";") %>% 
  rename(spec = species, manag = landuse) 
## change species names to the names I used
caro_fl$spec <- recode(caro_fl$spec, "Anth odo" = "Ant_odo", "Bro ere" = "Bro_ere", "Crep bie" = "Cre_bie", 
                       "Dia carth" = "Dia_car", "Lot cor" = "Lot_cor", "Med fal" = "Med_fal", "Pla lan" = "Pla_lan", 
                       "Sca och" = "Sca_och", "Tra orie" = "Tra_ori", "Trif prat" = "Tri pra")

caro_fl <- caro_fl %>% subset(spec %in% c("Bro_ere", "Dia_car", "Pla_lan", "Sca_och", "Tra_ori"))

caro_fl$manag <- ifelse(caro_fl$manag == "EM", "mowing", "grazing")
caro_fl$clima <- ifelse(caro_fl$climate == "amb", "ambient", "future")

caro_fl$starting_month <- as.Date(caro_fl$FFD_date, "%d.%m.%Y") %>% format("%m") %>% as.numeric()
caro_fl$ending_month <- as.Date(caro_fl$LFD_date, "%d.%m.%Y") %>% format("%m") %>% as.numeric()

caro_fl$start_month_name <- month.abb[caro_fl$starting_month]
caro_fl$ending_month_name <- month.abb[caro_fl$ending_month]

caro_fl$duration_month <- caro_fl$ending_month - caro_fl$starting_month + 1

## aggregate for the management
caro_duration <- aggregate(duration_month ~ spec + manag, FUN = mean, data = caro_fl)
caro_start <- aggregate(starting_month ~ manag + spec, FUN = mean, na.rm = T, data = caro_fl)
caro_end <- aggregate(ending_month ~ spec + manag, FUN = mean, na.rm = T, data = caro_fl ) 

## aggregate for climate
caro_duration_clim <- aggregate(duration_month ~ spec + clima, FUN = mean, data = caro_fl)
caro_start_clim <- aggregate(starting_month ~ clima + spec, FUN = mean, na.rm = T, data = caro_fl)
caro_end_clim <- aggregate(ending_month ~ clima + spec, FUN = mean, na.rm = T, data = caro_fl ) 

## Which species prefers what treatment
effect_mana <- subset(lamb_mana, manag == "grazing")[, c("lamb", "spec")] %>% rename(lamb_grazing = lamb)
effect_mana <- subset(lamb_mana, manag == "mowing")[, c("lamb", "spec")] %>% rename(lamb_mowing = lamb) %>% merge(effect_mana, ., by = "spec")

effect_clim <- subset(lamb_clim, clima == "ambient")[, c("lamb", "spec")] %>% rename(lamb_ambient = lamb)
effect_clim <- subset(lamb_clim, clima == "future")[, c("lamb", "spec")] %>% rename(lamb_future = lamb) %>% merge(effect_clim, ., by = "spec")

effect_mana$effect_size <- log(effect_mana$lamb_grazing) - log(effect_mana$lamb_mowing)
effect_mana$prefer <- ifelse(effect_mana$effect_size > 0 , "Prefers grazing", "Prefers mowing")

effect_clim$effect_size <- log(effect_clim$lamb_ambient) - log(effect_clim$lamb_future)
effect_clim$prefer <- ifelse(effect_clim$effect_size > 0 , "Prefers ambient", "Prefers future")


gra_eff <- ggplot(effect_mana, aes(x = spec, y = effect_size, fill = prefer)) + geom_bar(stat = "identity") +  
  scale_fill_manual(values = colorbl, name = "") + theme_bw() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1.05)) + 
  xlab("") + ylab(expression(paste(lambda, " grazing - ", lambda, " mowing"))) + theme(legend.position = "none") + 
  geom_text(aes(x = 0.75, y = 1, label = "B.)"), size = 5)

amb_eff <- ggplot(effect_clim, aes(x = spec, y = effect_size, fill = prefer)) + geom_bar(stat = "identity") +  
  scale_fill_manual(values = colorbl, name = "") + theme_bw() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1.05)) + 
  xlab("") + ylab(expression(paste(lambda, " grazing - ", lambda, " mowing"))) + theme(legend.position = "none") + 
  geom_text(aes(x = 0.75, y = 1, label = "B.)"), size = 5)

amb_eff
gra_eff

##  merge data with the phenology data
timing_lambda_mana <- data.frame(spec = rep(effect_mana$spec, 2),  manag = c(rep("grazing", 5), rep("mowing", 5)), 
                            lambda = c(log(effect_mana$lamb_grazing), log(effect_mana$lamb_mowing)),
                            effect_size = rep(effect_mana$effect_size, 2))

timing_lambda_clim <- data.frame(spec = rep(effect_clim$spec, 2),  clima = c(rep("ambient", 5), rep("future", 5)), 
                                 lambda = c(log(effect_clim$lamb_ambient), log(effect_clim$lamb_future)),
                                 effect_size = rep(effect_clim$effect_size, 2))


## Compare data with caros flowering stuff
timing_lambda_mana <- merge(caro_duration, caro_start) %>% merge(., caro_end) %>% merge(timing_lambda_mana)
timing_lambda_clim <- merge(caro_duration_clim, caro_start_clim) %>% merge(., caro_end_clim) %>% merge(timing_lambda_clim)


ggplot(timing_lambda_mana, aes(x = duration_month, y = lambda, label = spec)) + geom_point() + geom_smooth(method = "lm", se = F) + 
  facet_wrap(~ manag) + geom_text() + theme_bw()
ggplot(timing_lambda_clim, aes(x = duration_month, y = lambda, label = spec)) + geom_point() + geom_smooth(method = "lm", se = F) + 
  facet_wrap(~ clima) + geom_text() + theme_bw()

ggplot(timing_lambda_mana, aes(x = starting_month, y = lambda, label = spec)) + geom_point() + geom_smooth(method = "lm", se = F) + geom_text() + 
  facet_wrap(~ manag)

lm(lambda ~ duration_month, data = timing_lambda_mana) %>% summary()
lm(lambda ~ duration_month, data = timing_lambda_clim) %>% summary()

## plot for the linear model grazing - mowing using the effect size
fl_duration <- ggplot(timing_lambda_mana, aes(x = duration_month, y = effect_size)) +
  geom_point(aes(shape = manag, color = manag), size = 4, color = "black") + 
  theme_bw() + 
  ylab(expression(paste("Effect size ( ", lambda, " grazing - ", lambda, " mowing)"))) + 
  xlab("Duration of flowering (month)") + 
  scale_x_continuous(limits = c(0, 8), breaks = seq(0, 8, 1)) + 
  scale_y_continuous(limits = c(-2.5, 2)) + 
  labs(shape = "Management") + 
  scale_shape_manual(labels = c("grazing" = "Grazing", "mowing" = "Mowing"), values = c(17, 16), 
                     name = "Management") + 
  ggtitle("(a)") + 
  theme(text = element_text(size = 20))
fl_duration

fl_start <- ggplot(timing_lambda, aes(x = starting_month, y = effect_size)) + 
  geom_point(aes(shape = manag), size = 4) +
  scale_x_continuous(limits = c(4, 8), labels = c("April", "May", "June", "July", "August")) + 
  scale_y_continuous(limits = c(-2.5, 2)) + 
  theme_bw() + 
  labs(x = "Starting month of flowering", y = "",
       shape = "Management") + 
  scale_shape_manual(values = c(17, 16), labels = c("grazing" = "Grazing", "mowing" = "Mowing")) + 
  ggtitle("(b)") + 
  theme(text = element_text(size = 20))
fl_start

## Plots for climate
## Flowering duration
fl_duration_clim <- ggplot(timing_lambda_clim, aes(x = duration_month, y = effect_size)) +
  geom_point(aes(color = clima), size = 4) + 
  theme_bw() + 
  ylab(expression(paste("Effect size ( ", lambda, " ambient - ", lambda, " future)"))) + 
  xlab("Duration of flowering (month)") + 
  scale_x_continuous(limits = c(0, 8), breaks = seq(0, 8, 1)) + 
  scale_y_continuous(limits = c(-1.5, 1)) + 
  labs(shape = "Management") + 
  scale_color_manual(labels = c("Ambient", "Future"), values = rbPalette, 
                     name = "Climate") + 
  ggtitle("(a)") + 
  theme(text = element_text(size = 20))
fl_duration_clim

#Flowering start
fl_start_clim <- ggplot(timing_lambda_clim, aes(x = starting_month, y = effect_size)) + 
  geom_point(aes(color = clima), size = 4) +
  scale_x_continuous(limits = c(4, 8), labels = c("April", "May", "June", "July", "August")) + 
  scale_y_continuous(limits = c(-1.5, 1)) + 
  theme_bw() + 
  labs(x = "Starting month of flowering", y = "",
       shape = "Management") + 
  scale_color_manual(values = rbPalette, name = "Climate", labels = c("Ambient", "Future"))+ 
  ggtitle("(b)") + 
  theme(text = element_text(size = 20))
fl_start_clim


lm(effect_size ~ duration_month, data = timing_lambda_mana) %>% summary()
lm(effect_size ~ starting_month, data = timing_lambda_mana) %>% summary()


lm(effect_size ~ duration_month, data = timing_lambda_clim) %>% summary()
lm(effect_size ~ starting_month, data = timing_lambda_clim) %>% summary()

fl_plots <- ggpubr::ggarrange(fl_duration, fl_start, common.legend = T, legend = "right")
fl_plots

fl_plots_clim <- ggpubr::ggarrange(fl_duration_clim, fl_start_clim, common.legend = T, legend = "right")
fl_plots_clim

ggsave("C://Users/mandrzej24/Documents/01_revision_demography/03_figures/fl_effect_sizes_full_modell_mean_par.jpeg", 
        plot = fl_plots, device = "jpeg", dpi = 300, width = 3500, height = 2342, units = "px")
ggsave("C://Users/mandrzej24/Documents/01_revision_demography/03_figures/fl_effect_sizes_climate.jpeg", 
       plot = fl_plots_clim, device = "jpeg", dpi = 300, width = 3500, height = 2342, units = "px")
lm(effect_size ~ duration_month, data = timing_lambda) %>% summary()
lm(effect_size ~ duration_month, data = timing_lambda %>% subset(manag == "mowing")) %>% summary()
lm(effect_size ~ starting_month, data = timing_lambda) %>% summary()

## starting month 
ggplot(timing_lambda, aes(x = starting_month, y = effect_size, label = spec)) + geom_point() + 
  theme_bw() + geom_text() + geom_smooth(method = "lm", se = F)
 
lm(effect_size ~ starting_month, data = timing_lambda) %>% summary()
