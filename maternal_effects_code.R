
library(rethinking)
library(dplyr)

#play model

fc_model <- map2stan(
  alist(
    
    fc_total  ~ dzipois( p , lambda ),
    
    logit(p) <- ap + time + ap_Focal[Focal] + bp_age*age  + bp_ageSq*ageSq + bp_female*female +
      bp_multip*multip + bp_rank*relrank + 
      bp_grpsize*grpsize + bp_challenges*challenges + bp_grass*grass + 
      bp_opuntia*opuntia + bp_adv*adv +
      bp_bn*bn + bp_bh*bh + bp_gc*gc  ,
    
    log(lambda) <- al + time + al_Focal[Focal] + bl_age*age  + bl_ageSq*ageSq + bl_female*female +
      bl_multip*multip + bl_rank*relrank + 
      bl_grpsize*grpsize + bl_challenges*challenges + bl_grass*grass +
      bl_opuntia*opuntia + bl_adv*adv +
      bl_bn*bn + bl_bh*bh + bl_gc*gc  ,
    
    c(ap_Focal,al_Focal)[Focal] ~ dmvnormNC(sigma_focal,Rho_focal) ,
    c(ap,al,bp_age,bp_ageSq,bp_female,bp_multip,bp_rank,
      bp_grpsize,bp_challenges,bp_grass,bp_opuntia,bp_adv,bp_bn,bp_bh,bp_gc,
      bl_age,bl_ageSq,bl_female,bl_multip,bl_rank,
      bl_grpsize,bl_challenges,bl_grass,bl_opuntia,bl_adv,bl_bn,bl_bh,bl_gc) ~ dnorm(0,2) ,
    
    sigma_focal  ~ dcauchy(0,2)  ,
    Rho_focal ~ dlkjcorr(3) 
    
  ),
  data=list(
    fc_total=dd$fc_total,
    time = dd$s_mins,
    Focal=dd$Focal,
    age = dd$s_age,  
    ageSq = dd$s_age^2,
    female = dd$female,
    multip = dd$multiparous,
    relrank = dd$s_rank,
    opuntia = dd$s_ageopuntia,
    grpsize = dd$s_grpsize_cur,
    challenges = dd$s_challenges,
    grass = dd$s_grass, 
    bn = dd$s_bnRate,
    bh = dd$s_bhRate,
    gc = dd$s_gc,
    adv = dd$s_cumadvprimip
    
  ),
  cores=3 , chains=2 , warmup=1000, iter=3000, WAIC=TRUE, types=list(adapt.delta=0.99)
)


# independence model

leaves_model_sexInter <- map2stan(
  alist(
    
    leaves  ~ dzipois( p , lambda ),
    
    logit(p) <- ap + time + ap_Focal[Focal] + bp_age*age  + bp_ageSq*ageSq +  bp_female*female +
      bp_multip*multip + bp_rank*relrank + bp_opuntia*opuntia +
      bp_grpsize*grpsize + bp_challenges*challenges + bp_grass*grass + bp_adv*adv +
      bp_bn*bn + bp_bh*bh + bp_gc*gc + (bp_bn_female*bn + bp_bh_female*bh + bp_gc_female*gc)*female ,
    
    log(lambda) <- al + time + al_Focal[Focal] + bl_age*age  + bl_ageSq*ageSq + bl_female*female +
      bl_multip*multip + bl_rank*relrank + bl_opuntia*opuntia +
      bl_grpsize*grpsize + bl_challenges*challenges + bl_grass*grass + bl_adv*adv +
      bl_bn*bn + bl_bh*bh + bl_gc*gc + (bl_bn_female*bn + bl_bh_female*bh + bl_gc_female*gc)*female ,
    
    c(ap_Focal,al_Focal)[Focal] ~ dmvnormNC(sigma_focal,Rho_focal) ,
    c(ap,al,bp_age,bp_ageSq,bp_female,bp_multip,bp_rank,
      bp_grpsize,bp_challenges,bp_grass,bp_opuntia,bp_adv,bp_bn,bp_bh,bp_gc,
      bp_bn_female,bp_bh_female,bp_gc_female,
      bl_age,bl_ageSq,bl_female,bl_multip,bl_rank,
      bl_grpsize,bl_challenges,bl_grass,bl_opuntia,bl_adv,bl_bn,bl_bh,bl_gc,
      bl_bn_female,bl_bh_female,bl_gc_female) ~ dnorm(0,2) ,
    
    sigma_focal  ~ dcauchy(0,2)  ,
    Rho_focal ~ dlkjcorr(3) 
    
  ),
  data=list(
    leaves=dd$l1,
    time = dd$s_mins,
    Focal=dd$Focal,
    age = dd$s_age,  
    ageSq = dd$s_age^2,
    female = dd$female,
    multip = dd$multiparous,
    relrank = dd$s_rank,
    opuntia = dd$s_ageopuntia,
    grpsize = dd$s_grpsize_cur,
    challenges = dd$s_challenges,
    grass = dd$s_grass, 
    bn = dd$s_bnRate,
    bh = dd$s_bhRate,
    gc = dd$s_gc,
    adv = dd$s_cumadvprimip
    
  ),
  cores=3 , chains=2 , warmup=1000, iter=3000, WAIC=TRUE, types=list(adapt.delta=0.99)
)


# growth model
growth_normal_biomass2 <- map2stan(
  alist(
    
    growth_rate ~ dnorm(mu , sigma ),
    
    mu <- ap + bp_bn*bn + bp_bh*bh + bp_gc*gc + bp_adv*adv + bp_opuntia*opuntia +
      bp_multip*multip + bp_female*female + bp_rank*relrank +  bp_biomass*biomass,
    
    c(ap, bp_bn, bp_bh, bp_gc, bp_adv, bp_opuntia, bp_multip, bp_female, bp_rank,  bp_biomass) ~ dnorm(0,2) ,
    sigma ~ dcauchy(0,2)
    
  ),
  data=list(
    growth_rate=d$s_growth,
    bn = d$s_bn,
    bh = d$s_bh,
    gc = d$s_gc,
    adv = d$s_cumadvprimip,
    relrank = d$s_rank,
    multip = d$multiparous,
    female = d$female,
    biomass = d$s_grass,
    opuntia = d$s_ageopuntia
  ),
  cores=3 , chains=2 , warmup=3000, iter=6000, WAIC=TRUE, types=list(adapt.delta=0.99)
)


# trade-off model

fc_growth_interactions_model <- map2stan(
  alist(
    
    fc_total  ~ dzipois( p , lambda ),
    
    logit(p) <- ap + time + ap_Focal[Focal] + bp_age*age  + bp_ageSq*ageSq + bp_female*female +
      bp_multip*multip + bp_rank*relrank + 
      bp_grpsize*grpsize + bp_challenges*challenges + bp_grass*grass + 
      bp_opuntia*opuntia + bp_adv*adv +
      bp_bn*bn + bp_bh*bh + bp_gc*gc + bp_growth*growth + 
      (bp_bn_growth*bn + bp_bh_growth*bh + bp_gc_growth*gc)*growth,
    
    log(lambda) <- al + time + al_Focal[Focal] + bl_age*age  + bl_ageSq*ageSq + bl_female*female +
      bl_multip*multip + bl_rank*relrank + 
      bl_grpsize*grpsize + bl_challenges*challenges + bl_grass*grass +
      bl_opuntia*opuntia + bl_adv*adv +
      bl_bn*bn + bl_bh*bh + bl_gc*gc + bl_growth*growth +
      (bl_bn_growth*bn + bl_bh_growth*bh + bl_gc_growth*gc)*growth,
    
    c(ap_Focal,al_Focal)[Focal] ~ dmvnormNC(sigma_focal,Rho_focal) ,
    c(ap,al,bp_age,bp_ageSq,bp_female,bp_multip,bp_rank,
      bp_grpsize,bp_challenges,bp_grass,bp_opuntia,bp_adv,bp_bn,bp_bh,bp_gc,bp_growth,
      bp_bn_growth,bp_bh_growth,bp_gc_growth,
      bl_age,bl_ageSq,bl_female,bl_multip,bl_rank,
      bl_grpsize,bl_challenges,bl_grass,bl_opuntia,bl_adv,bl_bn,bl_bh,bl_gc,bl_growth,
      bl_bn_growth,bl_bh_growth,bl_gc_growth) ~ dnorm(0,2) ,
    
    sigma_focal  ~ dcauchy(0,2)  ,
    Rho_focal ~ dlkjcorr(3) 
    
  ),
  data=list(
    fc_total=dd$fc_total,
    time = dd$s_mins,
    Focal=dd$Focal,
    age = dd$s_age,  
    ageSq = dd$s_age^2,
    female = dd$female,
    multip = dd$multiparous,
    relrank = dd$s_rank,
    opuntia = dd$s_ageopuntia,
    grpsize = dd$s_grpsize_cur,
    challenges = dd$s_challenges,
    grass = dd$s_grass, 
    bn = dd$s_bnRate,
    bh = dd$s_bhRate,
    gc = dd$s_gc,
    adv = dd$s_cumadvprimip,
    growth = dd$s_growth
    
  ),
  cores=3 , chains=2 , warmup=1000, iter=3000, WAIC=TRUE, types=list(adapt.delta=0.99)
)


# plot for play

#Offspring play

#nursing time
a_foc_z <- matrix(0,1000,length(unique(dd$Focal)))
bn.seq=seq(min(dd$s_bnRate),max(dd$s_bnRate),length=1000)

d.pred_bn <- list(
  Focal=rep(1,length(bn.seq)),
  time=rep(mean(dd$s_mins),length(bn.seq)),
  age=rep(mean(dd$s_age),length(bn.seq)),
  ageSq=rep(mean(dd$s_age^2),length(bn.seq)),
  female=rep(mean(dd$female),length(bn.seq)),
  multip=rep(mean(dd$multiparous),length(bn.seq)),  
  relrank=rep(mean(dd$s_rank),length(bn.seq)),
  opuntia=rep(mean(dd$s_ageopuntia),length(bn.seq)),
  grpsize=rep(mean(dd$s_grpsize_cur),length(bn.seq)),
  challenges=rep(mean(dd$s_challenges),length(bn.seq)),
  grass=rep(mean(dd$s_grass),length(bn.seq)),
  bh=rep(mean(dd$s_bhRate),length(bn.seq)),
  bn=bn.seq,
  gc=rep(mean(dd$s_gc),length(bn.seq)),
  adv=rep(mean(dd$s_cumadvprimip),length(bn.seq))
)

LM1 <- link(fc_model, n=1000 , data=d.pred_bn, replace=
              list(am_Focal=a_foc_z), WAIC=TRUE)
pred1 <- (1-LM1$p)*LM1$lambda
pred1.median <- apply(pred1, 2, median )
pred1.HPDI <-apply(pred1, 2, HPDI )

par(mfrow = c(2, 2), mar=c(2,0,2,0),oma=c(4,4,3,3))

xxa=bn.seq
yya=as.vector(t(pred1[1:1000,]))
smoothScatter(rep(xxa,1000),yya,xlim=c(min(dd$s_bnRate),max(dd$s_bnRate)), colramp=colorRampPalette(c("white","#33CCFF")),nbin=200,transformation = function(x) x^.5,ylab='n',xlab='n',cex=1.2, yaxt='n', xaxt='n', ylim=c(0,1), nrpoints=0)
#points(fc_total ~ s_bnRate, data=dd , col=alpha("#33CCFF",0.6),pch=16, cex=0.9)

lines( bn.seq , pred1.median , lwd=1,col="black")
lines( bn.seq , pred1.HPDI[1,],lty=2,lwd=1)
lines( bn.seq , pred1.HPDI[2,],lty=2,lwd=1)

lab.a=seq(0,1,.2)
lab.sa=(lab.a-mean(dd$Bnminsrate))/sd(dd$Bnminsrate)

axis(1,labels=NA, at=lab.sa, tck=-0.01)
mtext(lab.a,at=lab.sa,side=1, line=.5)
axis(2, at = seq(from=0 , to=3, by = .5) ,tck=-0.01 )

mtext("Nursing rate", side=1, line=2.5, cex=2,tck=-0.01)
mtext(side=2,line=2,text="Play bouts",cex=2)


#carrying time
a_foc_z <- matrix(0,1000,length(unique(dd$Focal)))
bn.seq=seq(min(dd$s_bhRate),max(dd$s_bhRate),length=1000)

d.pred_bn <- list(
  Focal=rep(1,length(bn.seq)),
  time=rep(mean(dd$s_mins),length(bn.seq)),
  age=rep(mean(dd$s_age),length(bn.seq)),
  ageSq=rep(mean(dd$s_age^2),length(bn.seq)),
  female=rep(mean(dd$female),length(bn.seq)),
  multip=rep(mean(dd$multiparous),length(bn.seq)),  
  relrank=rep(mean(dd$s_rank),length(bn.seq)),
  opuntia=rep(mean(dd$s_ageopuntia),length(bn.seq)),
  grpsize=rep(mean(dd$s_grpsize_cur),length(bn.seq)),
  challenges=rep(mean(dd$s_challenges),length(bn.seq)),
  grass=rep(mean(dd$s_grass),length(bn.seq)),
  bn=rep(mean(dd$s_bnRate),length(bn.seq)),
  bh=bn.seq,
  gc=rep(mean(dd$s_gc),length(bn.seq)),
  adv=rep(mean(dd$s_cumadvprimip),length(bn.seq))
)

LM1 <- link(fc_model, n=1000 , data=d.pred_bn, replace=
              list(am_Focal=a_foc_z), WAIC=TRUE)
pred1 <- (1-LM1$p)*LM1$lambda
pred1.median <- apply(pred1, 2, median )
pred1.HPDI <-apply(pred1, 2, HPDI )

xxa=bn.seq
yya=as.vector(t(pred1[1:1000,]))
#par(mfrow = c(1, 1), cex=1.1, mar=c(0,0,0,0), oma=c(4,3,1.1,4))
smoothScatter(rep(xxa,1000),yya,xlim=c(min(dd$s_bhRate),max(dd$s_bhRate)), colramp=colorRampPalette(c("white","#33CCFF")),nbin=200,transformation = function(x) x^.5,ylab='n',xlab='n',cex=1.2, yaxt='n', xaxt='n', ylim=c(0,1), nrpoints=0)
#points(fc_total ~ s_bhRate, data=dd , col=alpha("#33CCFF",0.6),pch=16, cex=0.9)

lines( bn.seq , pred1.median , lwd=1,col="black")
lines( bn.seq , pred1.HPDI[1,],lty=2,lwd=1)
lines( bn.seq , pred1.HPDI[2,],lty=2,lwd=1)

lab.a=seq(0,.9,.2)
lab.sa=(lab.a-mean(dd$Bhminsrate))/sd(dd$Bhminsrate)

axis(1,labels=NA, at=lab.sa, tck=-0.01)
mtext(lab.a,at=lab.sa,side=1, line=.5)
#axis(2, at = seq(from=0 , to=3, by = 1) ,tck=-0.01 )

mtext("Carrying rate", side=1, line=2.5, cex=2,tck=-0.01)
#mtext(side=2,line=2,text="Bouts of social contact play",cex=1.5)


#GCMs
a_foc_z <- matrix(0,1000,length(unique(dd$Focal)))
bn.seq=seq(min(dd$s_gc),max(dd$s_gc),length=1000)

d.pred_bn <- list(
  Focal=rep(1,length(bn.seq)),
  time=rep(mean(dd$s_mins),length(bn.seq)),
  age=rep(mean(dd$s_age),length(bn.seq)),
  ageSq=rep(mean(dd$s_age^2),length(bn.seq)),
  female=rep(mean(dd$female),length(bn.seq)),
  multip=rep(mean(dd$multiparous),length(bn.seq)),  
  relrank=rep(mean(dd$s_rank),length(bn.seq)),
  opuntia=rep(mean(dd$s_ageopuntia),length(bn.seq)),
  grpsize=rep(mean(dd$s_grpsize_cur),length(bn.seq)),
  challenges=rep(mean(dd$s_challenges),length(bn.seq)),
  grass=rep(mean(dd$s_grass),length(bn.seq)),
  bn=rep(mean(dd$s_bnRate),length(bn.seq)),
  gc=bn.seq,
  bh=rep(mean(dd$s_bhRate),length(bn.seq)),
  adv=rep(mean(dd$s_cumadvprimip),length(bn.seq))
)

LM1 <- link(fc_model, n=1000 , data=d.pred_bn, replace=
              list(am_Focal=a_foc_z), WAIC=TRUE)
pred1 <- (1-LM1$p)*LM1$lambda
pred1.median <- apply(pred1, 2, median )
pred1.HPDI <-apply(pred1, 2, HPDI )

xxa=bn.seq
yya=as.vector(t(pred1[1:1000,]))
#par(mfrow = c(1, 1), cex=1.1, mar=c(0,0,0,0), oma=c(4,3,1.1,4))
smoothScatter(rep(xxa,1000),yya,xlim=c(min(dd$s_gc),max(dd$s_gc)), colramp=colorRampPalette(c("white","#33CCFF")),nbin=200,transformation = function(x) x^.5,ylab='n',xlab='n',cex=1.2, yaxt='n', xaxt='n', ylim=c(0,1), nrpoints=0)
#points(fc_total ~ s_gc, data=dd , col=alpha("#33CCFF",0.6),pch=16, cex=0.9)

lines( bn.seq , pred1.median , lwd=1,col="black")
lines( bn.seq , pred1.HPDI[1,],lty=2,lwd=1)
lines( bn.seq , pred1.HPDI[2,],lty=2,lwd=1)


lab.b=seq(440,5400,1000)
lab.sb=(lab.b-mean(dd$monthly_mean_GC))/sd(dd$monthly_mean_GC)
axis(1,labels=NA, at=lab.sb, tck=-0.01)
mtext(lab.b,at=lab.sb,side=1, line=.5)

axis(2, at = seq(from=0 , to=3, by = .5) ,tck=-0.01 )

mtext("Maternal GCMs", side=1, line=2.5, cex=2,tck=-0.01)
mtext(side=2,line=2,text="Play bouts",cex=2)




