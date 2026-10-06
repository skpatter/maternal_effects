
library(rethinking)
library(ggplot2)
library(dplyr)

# load data for each model 

dd <- read.csv("dd_infant_behavior_mediation.csv", header = TRUE, na.strings = "") 

# play model

dd$Focal <- as.factor(dd$Focal)
dd$Focal <- as.integer(as.factor(dd$Focal))

mediationPlay_GCeffort <- ulam(
  alist(
    
    ### 1. Mediator 1 Sub-Model: GC (Gaussian)
    gc ~ dnorm(mu_gc, sigma),
    mu_gc <- ag + ag_Focal[Focal] + bg_age * age +  bg_grpsize * grpsize + bg_rank * relrank +
      bg_female * female + bg_challenges * challenges + bg_grass * grass + bg_momage*momage + bg_momagesq*momageSq + #bg_multip*multip +
      bg_ageopuntia * ageop + bg_adv * adv + bg_adv_rank*adv*relrank,
    
    ### 3. Mediator 3 Sub-Model: Nursing (Beta)
    bn ~ beta(shape1_bn, shape2_bn),
    shape1_bn <- mu_bn * phi_bn,
    shape2_bn <- (1 - mu_bn) * phi_bn,
    logit(mu_bn) <- amn + amn_Focal[Focal] + bnm_age * age + bnm_adv * adv + 
      bnm_rank * relrank + bnm_challenge * challenges + 
      bnm_grass * grass + bnm_grpsize * grpsize + 
      bnm_opuntia * ageop + bnm_momage*momage + #bnm_multip*multip +
      bnm_female * female,
    
    ### 3. Mediator 3 Sub-Model: Carrying (Beta)
    bh ~ beta(shape1_bh, shape2_bh),
    shape1_bh <- mu_bh * phi_bh,
    shape2_bh <- (1 - mu_bh) * phi_bh,
    logit(mu_bh) <- amh + amh_Focal[Focal] + bhm_age * age + bhm_adv * adv + 
      bhm_rank * relrank + bhm_challenge * challenges + 
      bhm_grass * grass + bhm_grpsize * grpsize + 
      bhm_opuntia * ageop + bhm_momage*momage + #bhm_multip*multip + 
      bhm_female * female,
    
    ### 4. Main outcome Sub-Model: Play (ZIP)
    fc_total ~ dzipois(p, lambda),
    logit(p) <- ap + ap_Focal[Focal] + bp_age * age + bp_ageSq * ageSq + 
      bp_female * female + bp_multip * multip + bp_rank * relrank + 
      bp_grpsize * grpsize + bp_challenges * challenges + 
      bp_grass * grass + bp_opuntia * ageop + bp_adv * adv + 
      bp_gc * gc + bp_bn*bn + bp_bh*bh ,
    
    log(lambda) <- al + log(time) + al_Focal[Focal] + bl_age * age + 
      bl_ageSq * ageSq + bl_female * female + bl_multip * multip + 
      bl_rank * relrank + bl_grpsize * grpsize + 
      bl_challenges * challenges + bl_grass * grass + 
      bl_opuntia * ageop + bl_adv * adv +  
      bl_gc * gc + bl_bn*bn + bl_bh*bh,
    
    ag_Focal[Focal] ~ dnorm(0, sigma_focal),
    amn_Focal[Focal] ~ dnorm(0, sigma_focal),
    amh_Focal[Focal] ~ dnorm(0, sigma_focal),
    ap_Focal[Focal] ~ dnorm(0, sigma_focal),
    al_Focal[Focal] ~ dnorm(0, sigma_focal),
    
    c( amn, bnm_age, bnm_grpsize, bnm_opuntia, bnm_female, bnm_challenge, bnm_grass, bnm_adv, bnm_rank,bnm_momage, #bnm_multip, 
       amh, bhm_age, bhm_grpsize, bhm_opuntia, bhm_female, bhm_challenge, bhm_grass, bhm_adv, bhm_rank,bhm_momage, #bhm_multip, 
       
       ag, bg_age,bg_female,bg_grpsize, bg_rank, bg_challenges, bg_grass, bg_ageopuntia, bg_adv,bg_adv_rank,bg_momage, bg_momagesq,#bg_multip,
       ap, al, bp_age, bp_ageSq, bp_female, bp_multip, bp_rank, bp_grpsize, bp_challenges, bp_grass, bp_opuntia, 
       bp_adv, bp_gc,bp_bn,bp_bh,
       bl_age, bl_ageSq, bl_female, bl_multip, bl_rank, bl_grpsize, bl_challenges, bl_grass, bl_opuntia, 
       bl_adv, bl_gc,bl_bn,bl_bh
    ) ~ dnorm(0, 1),
    
    sigma_focal ~ dexp(1),
    c(phi_bn, phi_bh) ~ dexp(1),
    sigma ~ dexp(1)
  ),
  data = list(
    fc_total   = dd$fc_total_month,
    time       = dd$DurMins_total_month,
    Focal      = dd$Focal,
    age        = dd$s_age,   
    ageSq      = dd$s_agesquared,
    momage     = dd$s_MomAge,
    momageSq     = dd$s_MomAge^2,
    female     = dd$female,
    multip     = dd$multiparous,
    relrank    = dd$s_rank,
    ageop      = dd$s_ageopuntia,
    grpsize    = dd$s_grpsize_cur,
    challenges = dd$s_challenges,
    grass      = dd$s_grass, 
    gc         = dd$s_gc,
    bn         = dd$Bnminsrate_month0,
    bh         = dd$Bhminsrate_month0,
    adv        = dd$s_cumadvprimip
  ),
  cores = 3, chains = 1, warmup = 500, iter = 2000, log_lik = TRUE, 
  control = list(adapt_delta = 0.99)
)



# independence model

#same data as play model

mediationLeaves_GCeffort <- ulam(
  alist(
    
    ### 1. Mediator 1 Sub-Model: GC (Gaussian)
    gc ~ dnorm(mu_gc, sigma),
    mu_gc <- ag + ag_Focal[Focal] + bg_age * age + bg_grpsize * grpsize + bg_rank * relrank +
      bg_female * female + bg_challenges * challenges + bg_grass * grass + bg_momage*momage + bg_momagesq*momageSq + #bg_multip*multip +
      bg_ageopuntia * ageop + bg_adv * adv + bg_adv_rank*adv*relrank,
    
    ### 2. Mediator 2 Sub-Model: Nursing (Beta)
    bn ~ beta(shape1_bn, shape2_bn),
    shape1_bn <- mu_bn * phi_bn,
    shape2_bn <- (1 - mu_bn) * phi_bn,
    logit(mu_bn) <- amn + amn_Focal[Focal] + bnm_age * age + bnm_adv * adv + 
      bnm_rank * relrank + bnm_challenge * challenges + 
      bnm_grass * grass + bnm_grpsize * grpsize + 
      bnm_opuntia * ageop + bnm_momage*momage +
      bnm_female * female,
    
    ### 3. Mediator 3 Sub-Model: Carrying (Beta)
    bh ~ beta(shape1_bh, shape2_bh),
    shape1_bh <- mu_bh * phi_bh,
    shape2_bh <- (1 - mu_bh) * phi_bh,
    logit(mu_bh) <- amh + amh_Focal[Focal] + bhm_age * age + bhm_adv * adv + 
      bhm_rank * relrank + bhm_challenge * challenges + 
      bhm_grass * grass + bhm_grpsize * grpsize + 
      bhm_opuntia * ageop + bnm_momage*momage +
      bhm_female * female,
    
    ### 4. Main outcome Sub-Model Leaves from mother (ZIP)
    l1 ~ dzipois(p, lambda),
    logit(p) <- ap + ap_Focal[Focal] + bp_prox*prox + bp_age * age + bp_ageSq * ageSq + 
      bp_female * female + bp_multip * multip + bp_rank * relrank + 
      bp_grpsize * grpsize + bp_challenges * challenges + 
      bp_grass * grass + bp_opuntia * ageop + bp_adv * adv + 
      bp_gc * gc + bp_bn*bn + bp_bh*bh ,
    
    log(lambda) <- al + log(time) + al_Focal[Focal] + bl_prox*prox + bl_age * age + 
      bl_ageSq * ageSq + bl_female * female + bl_multip * multip + 
      bl_rank * relrank + bl_grpsize * grpsize + 
      bl_challenges * challenges + bl_grass * grass + 
      bl_opuntia * ageop + bl_adv * adv +  
      bl_gc * gc + bl_bn*bn + bl_bh*bh,
    
    ag_Focal[Focal] ~ dnorm(0, sigma_focal),
    amn_Focal[Focal] ~ dnorm(0, sigma_focal),
    amh_Focal[Focal] ~ dnorm(0, sigma_focal),
    ap_Focal[Focal] ~ dnorm(0, sigma_focal),
    al_Focal[Focal] ~ dnorm(0, sigma_focal),
    
    c( amn, bnm_age, bnm_grpsize, bnm_momage,bnm_opuntia, bnm_female, bnm_challenge, bnm_grass, bnm_adv, bnm_rank,
       amh, bhm_age, bhm_grpsize, bhm_momage,bhm_opuntia, bhm_female, bhm_challenge, bhm_grass, bhm_adv, bhm_rank,
       
       ag, bp_prox,bg_age,bg_female,bg_grpsize, bg_rank, bg_challenges, bg_grass, bg_ageopuntia, bg_adv,bg_momage,bg_momagesq,bg_adv_rank, #bg_multip,
       ap, al,bl_prox, bp_age, bp_ageSq, bp_female, bp_multip, bp_rank, bp_grpsize, bp_challenges, bp_grass, bp_opuntia, 
       bp_adv, bp_gc, bp_bn, bp_bh,
       bl_age, bl_ageSq, bl_female, bl_multip, bl_rank, bl_grpsize, bl_challenges, bl_grass, bl_opuntia, 
       bl_adv, bl_gc, bl_bn, bl_bh
    ) ~ dnorm(0, 1),
    
    sigma_focal ~ dexp(1),
    c(phi_bn, phi_bh) ~ dexp(1),
    sigma ~ dexp(1)
  ),
  data = list(
    l1         = dd$l1_total_month,
    time       = dd$DurMins_total_month,
    prox      = dd$s_prox,
    Focal      = dd$Focal,
    age        = dd$s_age,   
    ageSq      = dd$s_agesquared,
    momage     = dd$s_MomAge,
    momageSq   = dd$s_MomAge^2,
    female     = dd$female,
    multip     = dd$multiparous,
    relrank    = dd$s_rank,
    ageop      = dd$s_ageopuntia,
    grpsize    = dd$s_grpsize_cur,
    challenges = dd$s_challenges,
    grass      = dd$s_grass, 
    gc         = dd$s_gc,
    bn         = dd$Bnminsrate_month0,
    bh         = dd$Bhminsrate_month0,
    adv        = dd$s_cumadvprimip
  ),
  cores = 3, chains = 1, warmup = 500, iter = 2000, log_lik = TRUE, 
  control = list(adapt_delta = 0.99)
)


# growth model

dd <- read.csv("dd_growth_mediation.csv", header = TRUE, na.strings = "") 

mediationgrowth_effortGC <- ulam(
  alist(
    
    ### 1. Mediator 1 Sub-Model: GCs
    gc ~ dnorm(mu_gc, sigma_gc),
    mu_gc <- ag + bg_adv*adv + bg_rank*relrank + bg_grass*grass + 
      bg_ageopuntia*ageop + bg_momage*momage + bg_momageSq*momagesq + bg_female*female + 
      bg_adv_rank*adv*relrank,
    
    ### 2. Mediator 2 Sub-Model: Nursing
    bn ~ dnorm(mu_bn, sigma_bn),
    mu_bn <- amn + bnm_adv*adv + bnm_rank*relrank + bnm_grass*grass + 
      bnm_opuntia*ageop + bnm_momage*momage + bnm_female*female,
    
    ### 3. Mediator 3 Sub-Model: Carrying
    bh ~ dnorm(mu_bh, sigma_bh),
    mu_bh <- amh + bhm_adv*adv + bhm_rank*relrank + bhm_grass*grass + 
      bhm_opuntia*ageop + bhm_momage*momage + bhm_female*female,
    
    ### 4. Final Outcome Sub-Model: Growth
    growth_rate ~ dnorm(mu_gr, sigma_gr),
    mu_gr <- ap + bp_adv*adv + bp_bn*bn + bp_bh*bh + bp_gc*gc + 
      bp_opuntia*ageop + bp_multip*multip + bp_female*female + 
      bp_rank*relrank + bp_biomass*grass,
    
    c(ag, bg_adv, bg_rank, bg_adv_rank,bg_grass, bg_ageopuntia, bg_momage, bg_momageSq, bg_female, 
      amn, bnm_adv, bnm_rank, bnm_grass, bnm_opuntia, bnm_female, bnm_momage,
      amh, bhm_adv, bhm_rank, bhm_grass, bhm_opuntia, bhm_female, bhm_momage,
      ap, bp_adv, bp_bn, bp_bh, bp_gc,bp_opuntia, bp_multip, bp_female, bp_rank, bp_biomass
    ) ~ dnorm(0, 1),
    
    c(sigma_gc, sigma_bn, sigma_bh, sigma_gr) ~ dexp(1)
    
  ),
  data = list(
    growth_rate = dd$s_growth,
    bn          = dd$s_bn,
    bh          = dd$s_bh,
    gc          = dd$s_gc,
    adv         = dd$s_cumadvprimip,
    relrank     = dd$s_rank,
    multip      = dd$multiparous,
    momage      = dd$s_MomAge,
    momagesq    = dd$s_MomAge^2,
    female      = dd$female,
    grass       = dd$s_grass,
    challenges  = dd$s_challenges,
    ageop       = dd$s_ageopuntia
  ),
  cores = 3, chains = 3, warmup = 1000, iter = 3000, log_lik = TRUE, 
  control = list(adapt_delta = 0.99)
)



# trade-off model

dd <- read.csv("dd_daily_infant_behavior_tradeoff.csv", header = TRUE, na.strings = "") 

dd$Focal <- as.factor(dd$Focal)
dd$Focal <- as.integer(as.factor(dd$Focal))

fc_growth_interactions_model_Xsex <- ulam(
  alist(
    
    fc_total ~ dzipois( p , lambda ),
    
    #separated out to improve diagnostics
    logit(p) <- p_base + p_inter,
    p_base <- ap + ap_Focal[Focal] + 
      bp_age*age + bp_ageSq*ageSq + bp_female*female +bp_multip*multip + bp_rank*relrank + 
      bp_grpsize*grpsize + bp_challenges*challenges +  bp_grass*grass + bp_opuntia*opuntia + bp_adv*adv,
    p_inter <- bp_bn*bn + bp_bh*bh + bp_gc*gc + bp_growth*growth + 
      bp_bn_growth*bn*growth +  bp_bh_growth*bh*growth +  bp_gc_growth*gc*growth +
      bp_bn_female*bn*female + bp_bh_female*bh*female +   bp_gc_female*gc*female + bp_female_growth*female*growth + 
      bp_bn_growth_female*bn*female*growth +  bp_bh_growth_female*bh*female*growth + bp_gc_growth_female*gc*female*growth,
    
    log(lambda) <- l_base + l_inter,
    l_base <- al + time + al_Focal[Focal] + 
      bl_age*age + bl_ageSq*ageSq + bl_female*female + bl_multip*multip + bl_rank*relrank + 
      bl_grpsize*grpsize + bl_challenges*challenges + bl_grass*grass + bl_opuntia*opuntia + bl_adv*adv,
    l_inter <- bl_bn*bn + bl_bh*bh + bl_gc*gc + bl_growth*growth + 
      bl_bn_growth*bn*growth + bl_bh_growth*bh*growth + bl_gc_growth*gc*growth + 
      bl_bn_female*bn*female + bl_bh_female*bh*female + bl_gc_female*gc*female + bl_female_growth*female*growth + 
      bl_bn_growth_female*bn*female*growth + bl_bh_growth_female*bh*female*growth +  bl_gc_growth_female*gc*female*growth,
    
    ap_Focal[Focal] ~ dnorm(0, 1),
    al_Focal[Focal] ~ dnorm(0, 1),
    
    c(ap, al, bp_age, bp_ageSq, bp_female, bp_multip, bp_rank,
      bp_grpsize, bp_challenges, bp_grass, bp_opuntia, bp_adv, bp_bn, bp_bh, bp_gc, bp_growth,
      bp_bn_growth, bp_bh_growth, bp_gc_growth, bp_bn_female, bp_bh_female, bp_gc_female,bp_female_growth, 
      bp_bn_growth_female, bp_bh_growth_female, bp_gc_growth_female,
      bl_age, bl_ageSq, bl_female, bl_multip, bl_rank,
      bl_grpsize, bl_challenges, bl_grass, bl_opuntia, bl_adv, bl_bn, bl_bh, bl_gc, bl_growth,
      bl_bn_growth, bl_bh_growth, bl_gc_growth, bl_bn_female, bl_bh_female, bl_gc_female,bl_female_growth, 
      bl_bn_growth_female, bl_bh_growth_female, bl_gc_growth_female
    ) ~ dnorm(0, 1)
    
  ),
  data = list(
    fc_total   = dd$fc_total,
    time       = dd$s_mins,
    Focal      = as.integer(dd$Focal),
    age        = dd$s_age,  
    ageSq      = dd$s_age^2,
    female     = dd$female,
    multip     = dd$multiparous,
    relrank    = dd$s_rank,
    opuntia    = dd$s_ageopuntia,
    grpsize    = dd$s_grpsize_cur,
    challenges = dd$s_challenges,
    grass      = dd$s_grass, 
    bn         = dd$s_bnRate,
    bh         = dd$s_bhRate,
    gc         = dd$s_gc,
    adv        = dd$s_cumadvprimip,
    growth     = dd$s_growth
  ),
  cores = 3, chains = 2, warmup = 1000, iter = 3000, log_lik = TRUE
)



# Mediation Outcomes + Forest Plot

#extract posterior
post_play <- extract.samples(mediationPlay_GCeffort)

delta <- 1         # one standard deviation increase 
time_offset <- 60   # set observation time 
inv_logit <- function(x) 1 / (1 + exp(-x))


#### expected play count
calc_ey_joint_play <- function(gc_val, bn_val, bh_val, 
                               adv_val=0, rank_val=0, grass_val=0, chall_val=0) {
  
  #zip model
  # p component
  logit_p <- post_play$ap + 
    post_play$bp_rank * rank_val + 
    post_play$bp_challenges * chall_val + 
    post_play$bp_grass * grass_val + 
    post_play$bp_adv * adv_val + 
    post_play$bp_gc * gc_val + 
    post_play$bp_bn * bn_val + 
    post_play$bp_bh * bh_val
  p <- inv_logit(logit_p)
  
  # l component
  log_lambda <- post_play$al + log(time_offset) + 
    post_play$bl_rank * rank_val + 
    post_play$bl_challenges * chall_val + 
    post_play$bl_grass * grass_val + 
    post_play$bl_adv * adv_val + 
    post_play$bl_gc * gc_val + 
    post_play$bl_bn * bn_val + 
    post_play$bl_bh * bh_val
  lambda <- exp(log_lambda)
  
  return((1 - p) * lambda)
}

## decomposition across predictors
decompose_variable_play <- function(var_name) {
  
  adv_0   <- 0; rank_0  <- 0; grass_0 <- 0; chall_0 <- 0
  adv_1   <- 0; rank_1  <- 0; grass_1 <- 0; chall_1 <- 0
  
  if (var_name == "adv")        adv_1   <- delta
  if (var_name == "relrank")    rank_1  <- delta
  if (var_name == "grass")      grass_1 <- delta
  if (var_name == "challenges") chall_1 <- delta
  
  # mediator1: GC (Gaussian) *note has interaction btwn rank and adv
  gc_0 <- post_play$ag + post_play$bg_adv * adv_0 + post_play$bg_rank * rank_0 + 
    post_play$bg_adv_rank * (adv_0 * rank_0) + post_play$bg_grass * grass_0 + post_play$bg_challenges * chall_0
  
  gc_1 <- post_play$ag + post_play$bg_adv * adv_1 + post_play$bg_rank * rank_1 + 
    post_play$bg_adv_rank * (adv_1 * rank_1) + post_play$bg_grass * grass_1 + post_play$bg_challenges * chall_1
  
  # mediator2: BN = Nursing (Beta) 
  bn_0 <- inv_logit(post_play$amn + post_play$bnm_adv * adv_0 + post_play$bnm_rank * rank_0 + 
                      post_play$bnm_grass * grass_0 + post_play$bnm_challenge * chall_0)
  
  bn_1 <- inv_logit(post_play$amn + post_play$bnm_adv * adv_1 + post_play$bnm_rank * rank_1 + 
                      post_play$bnm_grass * grass_1 + post_play$bnm_challenge * chall_1)
  
  # mediator3: BH = Carrying (Beta)
  bh_0 <- inv_logit(post_play$amh + post_play$bhm_adv * adv_0 + post_play$bhm_rank * rank_0 + 
                      post_play$bhm_grass * grass_0 + post_play$bhm_challenge * chall_0)
  
  bh_1 <- inv_logit(post_play$amh + post_play$bhm_adv * adv_1 + post_play$bhm_rank * rank_1 + 
                      post_play$bhm_grass * grass_1 + post_play$bhm_challenge * chall_1)
  
  # the counterfactuals 
  Y_0000 <- calc_ey_joint_play(gc_0, bn_0, bh_0, adv_0, rank_0, grass_0, chall_0) #Baseline
  Y_1000 <- calc_ey_joint_play(gc_1, bn_0, bh_0, adv_0, rank_0, grass_0, chall_0) #shift GC
  Y_0100 <- calc_ey_joint_play(gc_0, bn_1, bh_0, adv_0, rank_0, grass_0, chall_0) #shift BN
  Y_0010 <- calc_ey_joint_play(gc_0, bn_0, bh_1, adv_0, rank_0, grass_0, chall_0) #shift BH
  Y_1111 <- calc_ey_joint_play(gc_1, bn_1, bh_1, adv_1, rank_1, grass_1, chall_1) #total 
  Y_0111 <- calc_ey_joint_play(gc_1, bn_1, bh_1, adv_0, rank_0, grass_0, chall_0) #direct 
  
  #paths
  IE_gc <- Y_1000 - Y_0000
  IE_bn <- Y_0100 - Y_0000
  IE_bh <- Y_0010 - Y_0000
  IE_total <- IE_gc + IE_bn + IE_bh
  DE<- Y_1111 - Y_0111
  TE <- Y_1111 - Y_0000
  
  data.frame(
    Variable = var_name,
    IE_gc    = IE_gc,
    IE_bn    = IE_bn,
    IE_bh    = IE_bh,
    IE_total = IE_total,
    DE       = DE,
    TE       = TE
  )
}

 
### Summary Table (89% Credible Intervals)

var_list <- c("adv", "relrank", "grass", "challenges")
all_draws_play <- do.call(rbind, lapply(var_list, decompose_variable_play))

summarize_play_89 <- function(df, var_label) {
  sub <- df[df$Variable == var_label, ]
  get_row <- function(draws, name) {
    m   <- mean(draws, na.rm = TRUE)
    te  <- mean(sub$TE, na.rm = TRUE)
    q89 <- HPDI(draws, prob =0.89) #q89 <- quantile(draws, probs = c(0.055, 0.945), na.rm = TRUE)
    
    data.frame(
      Variable                   = var_label,
      Effect                     = name,
      `Mean Count Change (fc)`   = m,
      StdDev                     = sd(draws, na.rm = TRUE),
      `5.5%`                     = q89[1],
      `94.5%`                    = q89[2],
      `Prop. Mediated (m/TE)`    = ifelse(name == "Direct Effect", NA, m / te),
      check.names                = FALSE,
      row.names                  = NULL
    )
  }
  rbind(
    get_row(sub$IE_gc,    "Indirect via GC"),
    get_row(sub$IE_bn,    "Indirect via Nursing"),
    get_row(sub$IE_bh,    "Indirect via Carrying"),
    get_row(sub$IE_total, "Combined Indirect (GC+Nursing+Carrying)"),
    get_row(sub$DE,       "Direct Effect"),
    get_row(sub$TE,       "Total Effect")
  )
}

final_play_summary_89 <- do.call(rbind, lapply(var_list, function(v) summarize_play_89(all_draws_play, v)))

### forest plot
plot_data_play <- final_play_summary_89 %>%
  filter(Effect != "Combined Indirect (GC+Nursing+Carrying)") %>%
  mutate(
    Variable = factor(Variable, 
                      levels = c("adv", "relrank", "grass", "challenges"),
                      labels = c("Maternal early adversity", "Maternal rank", "Current biomass", "Current challenges")),
    Effect = factor(Effect, 
                    levels = rev(c("Total Effect", "Direct Effect", 
                                   "Indirect via GC", "Indirect via Nursing", "Indirect via Carrying")))
  )
ggplot(plot_data_play, aes(x = `Mean Count Change (fc)`, y = Effect, color = Effect)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray50", size = 0.6) +
  geom_errorbarh(aes(xmin = `5.5%`, xmax = `94.5%`), height = 0.25, size = 0.8) +
  geom_point(size = 3) +
  facet_wrap(~ Variable, scales = "free_x", ncol = 2) +
  scale_color_manual(values = c(
    "Total Effect"     = "#2B5C8F",  # Deep Blue
    "Direct Effect"    = "#D95F02",  # Burnt Orange
    "Indirect via GC"  = "#7570B3",  # Muted Purple
    "Indirect via Nursing"  = "#1B9E77",  # Teal Green
    "Indirect via Carrying"  = "#E7298A"   # Magenta
  )) +
  labs(
    x = "Expected change in play counts",
    y = NULL,
    title = "Play: Direct and indirect pathways via GCs, Nursing, and Carrying"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    panel.border       = element_rect(color = "gray80", fill = NA, size = 0.8),
    strip.background   = element_rect(fill = "gray95", color = "gray80"),
    strip.text         = element_text(face = "bold", size = 11),
    legend.position    = "none",
    axis.text.y        = element_text(face = "bold", size = 10),
    panel.grid.minor   = element_blank(),
    plot.title         = element_text(face = "bold", size = 14),
    plot.subtitle      = element_text(color = "gray30", size = 11)
  )

