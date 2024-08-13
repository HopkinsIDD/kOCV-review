
# Functions for panel plots
# colored by vaccine product type 

# function assigning colors 
assign_colors <- function(fpdata){
  
  fpdata <- fpdata %>% 
    mutate(vaccine_product = case_when(
      StudyID == "Sialubanje et al, 2022" ~ "Euvichol-plus",
      StudyID == "Ali et al, 2021" ~ "Shanchol",
      StudyID == "Ferreras et al, 2018" ~ "Shanchol",
      StudyID == "Grandesso et al, 2019" ~ "Shanchol",
      StudyID == "Franke et al, 2018" ~ "Shanchol",
      StudyID == "Qadri et al, 2018" ~ "Shanchol",
      StudyID == "Malembaka et al. 2023" ~ "Euvichol-plus",
      StudyID == "Matias et al, 2023" ~ "Euvichol",
      StudyID == "Bhattacharya et al, 2013" ~ "Shanchol",
      StudyID == "Sur et al, 2009" ~ "Shanchol",
      StudyID == "Sur et al, 2011" ~ "Shanchol",
      StudyID == "Trach et al, 1997" ~ "WC",
      StudyID == "Wierzba et al, 2015" ~ "Shanchol",
      StudyID == "Ivers et al, 2015" ~ "Shanchol",
      StudyID == "Luquero et al, 2014" ~ "Shanchol",
      StudyID == "van Loon et al, 1996" ~ "WC",
      StudyID == "Clemens et al, 1990" ~ "WC",
      StudyID == "Clemens et al, 1988" ~ "WC",
      StudyID == "Qadri et al, 2015" ~ "Shanchol",
      StudyID == "Qadri et al, 2016" ~ "Shanchol",
      StudyID == "Azman et al, 2016" ~ "Shanchol"
    )) %>%
    mutate(color = case_when(vaccine_product == "Shanchol" ~ "#377EB8",
                             vaccine_product == "WC" ~ "#984EA3",
                             vaccine_product == "Euvichol-plus" ~ "#4DAF4A",
                             vaccine_product == "Euvichol" ~ "#FF7F00"))
  
  return(fpdata)
}


draw_2dose_efficacy_estimates_colored <- function(fpdata){
  
  fpdata <- fpdata %>% filter(!is.na(yearid))
  fpdata$fPaper <- factor(fpdata$Paper, levels = unique(fpdata$Paper))
  fpdata$uid <- 1:nrow(fpdata)
  fpdata$fuid <- factor(fpdata$uid)
  tau_estim_method <- "EB" # seems like a reasonable choice from DOI: 10.1002/sim.2688
  
  ##############
  # 2-dose efficacy
  ##############
  
  ## set up efficacy data for 2 dose
  fpdata_rct <- fpdata %>% mutate(se = ifelse(is.na(VE.u),
                                              (log(1-Mean.VE)-log(1-VE.l))/qnorm(0.95),
                                              (log(1-VE.l) - log(1-VE.u))/(2*qnorm(.975))),
                                  yi = log(1-Mean.VE)) %>%
    arrange((yearid)) %>%
    filter(StudyType=="Efficacy")  %>%
    filter(Dose==2) %>%
    filter(yearid!=0)
  
  # overall efficacy
  log_ve_fup_rct <- rma(yi=log(1-fpdata_rct$Mean.VE),
                        sei=fpdata_rct$se,
                        method=tau_estim_method)
  
  # estimate yearly pe
  fpdata_rct_1 <- subset(fpdata_rct,yearid==1)
  fpdata_rct_rma_1 <- rma(yi=log(1-fpdata_rct_1$Mean.VE),
                          sei=fpdata_rct_1$se,
                          method=tau_estim_method)
  fpdata_rct_2 <- subset(fpdata_rct,yearid==2)
  fpdata_rct_rma_2 <- rma(yi=log(1-fpdata_rct_2$Mean.VE),
                          sei=fpdata_rct_2$se,
                          method=tau_estim_method)
  fpdata_rct_3 <- subset(fpdata_rct,yearid==3)
  fpdata_rct_rma_3 <- rma(yi=log(1-fpdata_rct_3$Mean.VE),
                          sei=fpdata_rct_3$se,
                          method=tau_estim_method)
  fpdata_rct_4 <- subset(fpdata_rct,yearid==4)
  fpdata_rct_rma_4 <- rma(yi=log(1-fpdata_rct_4$Mean.VE),
                          sei=fpdata_rct_4$se,
                          method=tau_estim_method)
  
  # add pooled estimates
  ci.lb.rct=c(fpdata_rct_rma_1$ci.lb,fpdata_rct_rma_2$ci.lb,fpdata_rct_rma_3$ci.lb,fpdata_rct_rma_4$ci.lb)
  ci.ub.rct=c(fpdata_rct_rma_1$ci.ub,fpdata_rct_rma_2$ci.ub,fpdata_rct_rma_3$ci.ub,fpdata_rct_rma_4$ci.ub)
  yi.rct=c(fpdata_rct_rma_1$b,fpdata_rct_rma_2$b,fpdata_rct_rma_3$b,fpdata_rct_rma_4$b)
  
  ci.lb.rct=100*(1-exp(ci.lb.rct))
  ci.ub.rct=100*(1-exp(ci.ub.rct))
  yi.rct=100*(1-exp(yi.rct))
  
  # CI
  b.cr.lb.rct=100*(1-exp(log_ve_fup_rct$yi-qnorm(0.975)*sqrt(log_ve_fup_rct$vi)))
  b.cr.ub.rct=100*(1-exp(log_ve_fup_rct$yi+qnorm(0.975)*sqrt(log_ve_fup_rct$vi)))
  b.cr.rct=100*(1-exp(log_ve_fup_rct$yi))
  
  ###########################
  # Updated Horizontal Figure with just Efficacy Study for 2 dose
  ########################
  index_count=table(fpdata_rct$yearid)  # Number of studies per time period
  index_center=c(4.5,13.5,22.5,31.5,40.5) # mid point of each panel
  index_left=index_center-((index_count+1)/2-0.5)  # left most index of each panel
  index_right=index_center+((index_count+1)/2-0.5)  # right most index of each panel
  index=NULL
  for (i in 1:5){
    index=c(index,seq(index_left[i],length.out=index_count[i]))
  }
  
  # point size
  wi <- 1/sqrt(log_ve_fup_rct$vi)
  psize <- wi/sum(wi, na.rm = TRUE)
  psize <- (psize - min(psize, na.rm = TRUE))/(max(psize,na.rm = TRUE) - min(psize, na.rm = TRUE))
  psize <- (psize * 1) + 0.5
  
  # base plot
  #par(mfrow=c(2,1),mar=c(7,7,7,7))
  par(mar=c(4,4,6,6))
  plot(1:43,rep(NA,length.out=43),pch=15,
       ylim=c(-50,100),col="grey",
       xaxt = "n",ylab="kOCV Efficacy (%)",bty="n",
       xlab="")
  
  points(index,100*(1-exp(log_ve_fup_rct$yi)),pch=15,
         cex=psize,col=fpdata_rct$color,
         xaxt = "n",xlab="",ylab="kOCV Efficacy (%)",bty="n")
  
  # add CI
  nonarrowindex=which(b.cr.ub.rct>=(-50))
  arrowindex=which(b.cr.ub.rct<(-50))
  segments(index-0.3,b.cr.lb.rct,index+0.3,b.cr.lb.rct,col=fpdata_rct$color)
  segments(index[nonarrowindex],b.cr.lb.rct[nonarrowindex],
           index[nonarrowindex],b.cr.ub.rct[nonarrowindex],col=fpdata_rct$color[nonarrowindex])
  segments(index[nonarrowindex]-0.3,b.cr.ub.rct[nonarrowindex],
           index[nonarrowindex]+0.3,b.cr.ub.rct[nonarrowindex],col=fpdata_rct$color[nonarrowindex])
  # add CI for studies with large CI
  segments(index[arrowindex],b.cr.lb.rct[arrowindex],
           index[arrowindex],rep(-50,length.out=length(arrowindex)),col= fpdata_rct$color[arrowindex])
  for (i in 1:length(arrowindex)){
    polygon(y=c(-50,-46,-46,-50),
            x=c(index[arrowindex[i]],index[arrowindex[i]]-0.2 ,
                index[arrowindex[i]]+0.2,index[arrowindex[i]]),
            col = fpdata_rct$color[arrowindex[i]],border=fpdata_rct$color[arrowindex[i]])
  }
  
  # add pooled estimates
  for (i in 1:length(unique(fpdata_rct$yearid))){
    polygon(y = c(ci.lb.rct[i], yi.rct[i], ci.ub.rct[i], yi.rct[i]),
            x=c(index_right[i],index_right[i]-0.3,
                index_right[i],index_right[i]+0.3),col="black", border = "black")
  }
  
  # add axis and text
  axis(side=4, at=c(-50,0,50,100))
  axis(side=1,pos=-70,xlab="",tck=0,col="white",
       at=index_center[1:4],
       labels=c(paste(format(yi.rct,digits=2),"% [",
                      as.numeric(format(ci.ub.rct,digits=2)),
                      ",",
                      as.numeric(format(ci.lb.rct,digits=2)),"]",sep="")))
  axis(side=1,pos=-70,xlab="",tck=0,col.axis="black",col="white",
       at=index_center[5],labels=c("81% [42,94]")) # display original figures in the paper instead of reconstructed 
  # axis(side=1,pos=-75,xlab="",tck=0,col="white",
  #      at=index_center[4],labels=c("26% [-46,63]"))
  axis(side=1,pos=-50,xlab="",tck=0,lty=3,col="dimgrey",font=2,
       at=index_center,
       labels=c("0-12m","12-24m","24-36m","36-48m","48-60m"))
  axis(side=3,pos=100,xlab="",tck=0,lty=3,at=index_center[3],
       labels=c("Two-dose efficacy estimates"))
  abline(h=c(100,0,-50),lty=3,col="grey")
  title(main = "A", adj = 0, line = 1.5)
  axis(side=1,pos=-90,xlab="",tck=0,lty=3,col="white",
       at=22.5,
       # labels=c("Follow-up duration in month"),
       labels = c(""))
  
  # plot <- recordPlot()
  # return(plot)
  # 
  # # return the plot object 
  # grid.echo()
  # Fig <- grid.grab()
  
}


# For sensitive analysis (removing pre-shanchol vaccines)
draw_2dose_efficacy_estimates_colored_SSA <- function(fpdata){
  
  fpdata <- fpdata %>% filter(!is.na(yearid))
  fpdata$fPaper <- factor(fpdata$Paper, levels = unique(fpdata$Paper))
  fpdata$uid <- 1:nrow(fpdata)
  fpdata$fuid <- factor(fpdata$uid)
  tau_estim_method <- "EB" # seems like a reasonable choice from DOI: 10.1002/sim.2688
  
  ##############
  # 2-dose efficacy
  ##############
  
  ## set up efficacy data for 2 dose
  fpdata_rct <- fpdata %>% mutate(se = ifelse(is.na(VE.u),
                                              (log(1-Mean.VE)-log(1-VE.l))/qnorm(0.95),
                                              (log(1-VE.l) - log(1-VE.u))/(2*qnorm(.975))),
                                  yi = log(1-Mean.VE)) %>%
    arrange((yearid)) %>%
    filter(StudyType=="Efficacy")  %>%
    filter(Dose==2) %>%
    filter(yearid!=0)
  
  ### remove pre-shanchol vaccies 
  fpdata_rct <- fpdata_rct %>% 
    filter(Paper %in% c("Trach et al, 1997", "van Loon et al, 1996a") == FALSE)
  
  # overall efficacy
  log_ve_fup_rct <- rma(yi=log(1-fpdata_rct$Mean.VE),
                        sei=fpdata_rct$se,
                        method=tau_estim_method)
  
  # estimate yearly pe
  fpdata_rct_1 <- subset(fpdata_rct,yearid==1)
  fpdata_rct_rma_1 <- rma(yi=log(1-fpdata_rct_1$Mean.VE),
                          sei=fpdata_rct_1$se,
                          method=tau_estim_method)
  fpdata_rct_2 <- subset(fpdata_rct,yearid==2)
  fpdata_rct_rma_2 <- rma(yi=log(1-fpdata_rct_2$Mean.VE),
                          sei=fpdata_rct_2$se,
                          method=tau_estim_method)
  fpdata_rct_3 <- subset(fpdata_rct,yearid==3)
  fpdata_rct_rma_3 <- rma(yi=log(1-fpdata_rct_3$Mean.VE),
                          sei=fpdata_rct_3$se,
                          method=tau_estim_method)
  fpdata_rct_4 <- subset(fpdata_rct,yearid==4)
  fpdata_rct_rma_4 <- rma(yi=log(1-fpdata_rct_4$Mean.VE),
                          sei=fpdata_rct_4$se,
                          method=tau_estim_method)
  
  # add pooled estimates
  ci.lb.rct=c(fpdata_rct_rma_1$ci.lb,fpdata_rct_rma_2$ci.lb,fpdata_rct_rma_3$ci.lb,fpdata_rct_rma_4$ci.lb)
  ci.ub.rct=c(fpdata_rct_rma_1$ci.ub,fpdata_rct_rma_2$ci.ub,fpdata_rct_rma_3$ci.ub,fpdata_rct_rma_4$ci.ub)
  yi.rct=c(fpdata_rct_rma_1$b,fpdata_rct_rma_2$b,fpdata_rct_rma_3$b,fpdata_rct_rma_4$b)
  
  ci.lb.rct=100*(1-exp(ci.lb.rct))
  ci.ub.rct=100*(1-exp(ci.ub.rct))
  yi.rct=100*(1-exp(yi.rct))
  
  # CI
  b.cr.lb.rct=100*(1-exp(log_ve_fup_rct$yi-qnorm(0.975)*sqrt(log_ve_fup_rct$vi)))
  b.cr.ub.rct=100*(1-exp(log_ve_fup_rct$yi+qnorm(0.975)*sqrt(log_ve_fup_rct$vi)))
  b.cr.rct=100*(1-exp(log_ve_fup_rct$yi))
  
  ###########################
  # Updated Horizontal Figure with just Efficacy Study for 2 dose
  ########################
  index_count=table(fpdata_rct$yearid)  # Number of studies per time period
  index_center=c(4.5,13.5,22.5,31.5,40.5) # mid point of each panel
  index_left=index_center-((index_count+1)/2-0.5)  # left most index of each panel
  index_right=index_center+((index_count+1)/2-0.5)  # right most index of each panel
  index=NULL
  for (i in 1:5){
    index=c(index,seq(index_left[i],length.out=index_count[i]))
  }
  
  # point size
  wi <- 1/sqrt(log_ve_fup_rct$vi)
  psize <- wi/sum(wi, na.rm = TRUE)
  psize <- (psize - min(psize, na.rm = TRUE))/(max(psize,na.rm = TRUE) - min(psize, na.rm = TRUE))
  psize <- (psize * 1) + 0.5
  
  # base plot
  #par(mfrow=c(2,1),mar=c(7,7,7,7))
  par(mar=c(4,4,6,6))
  plot(1:43,rep(NA,length.out=43),pch=15,
       ylim=c(-50,100),col="grey",
       xaxt = "n",ylab="kOCV Efficacy (%)",bty="n",
       xlab="")
  
  points(index,100*(1-exp(log_ve_fup_rct$yi)),pch=15,
         cex=psize,col=fpdata_rct$color,
         xaxt = "n",xlab="",ylab="kOCV Efficacy (%)",bty="n")
  
  # add CI
  nonarrowindex=which(b.cr.ub.rct>=(-50))
  arrowindex=which(b.cr.ub.rct<(-50))
  segments(index-0.3,b.cr.lb.rct,index+0.3,b.cr.lb.rct,col=fpdata_rct$color)
  segments(index[nonarrowindex],b.cr.lb.rct[nonarrowindex],
           index[nonarrowindex],b.cr.ub.rct[nonarrowindex],col=fpdata_rct$color[nonarrowindex])
  segments(index[nonarrowindex]-0.3,b.cr.ub.rct[nonarrowindex],
           index[nonarrowindex]+0.3,b.cr.ub.rct[nonarrowindex],col=fpdata_rct$color[nonarrowindex])
  # add CI for studies with large CI
  segments(index[arrowindex],b.cr.lb.rct[arrowindex],
           index[arrowindex],rep(-50,length.out=length(arrowindex)),col= fpdata_rct$color[arrowindex])
  # for (i in 1:length(arrowindex)){
  #   polygon(y=c(-50,-46,-46,-50),
  #           x=c(index[arrowindex[i]],index[arrowindex[i]]-0.2 ,
  #               index[arrowindex[i]]+0.2,index[arrowindex[i]]),
  #           col = fpdata_rct$color[arrowindex[i]],border=fpdata_rct$color[arrowindex[i]])
  # }
  
  # add pooled estimates
  for (i in 1:length(unique(fpdata_rct$yearid))){
    polygon(y = c(ci.lb.rct[i], yi.rct[i], ci.ub.rct[i], yi.rct[i]),
            x=c(index_right[i],index_right[i]-0.3,
                index_right[i],index_right[i]+0.3),col="black", border = "black")
  }
  
  # add axis and text
  axis(side=4, at=c(-50,0,50,100))
  axis(side=1,pos=-70,xlab="",tck=0,col="white",
       at=index_center[1:4],
       labels=c(paste(format(yi.rct,digits=2),"% [",
                      as.numeric(format(ci.ub.rct,digits=2)),
                      ",",
                      as.numeric(format(ci.lb.rct,digits=2)),"]",sep="")))
  axis(side=1,pos=-70,xlab="",tck=0,col.axis="black",col="white",
       at=index_center[5],labels=c("81% [42,94]")) # display original figures in the paper instead of reconstructed 
  # axis(side=1,pos=-75,xlab="",tck=0,col="white",
  #      at=index_center[4],labels=c("26% [-46,63]"))
  axis(side=1,pos=-50,xlab="",tck=0,lty=3,col="dimgrey",font=2,
       at=index_center,
       labels=c("0-12m","12-24m","24-36m","36-48m","48-60m"))
  axis(side=3,pos=100,xlab="",tck=0,lty=3,at=index_center[3],
       labels=c("Two-dose efficacy estimates"))
  abline(h=c(100,0,-50),lty=3,col="grey")
  title(main = "A", adj = 0, line = 1.5)
  axis(side=1,pos=-90,xlab="",tck=0,lty=3,col="white",
       at=22.5,
       # labels=c("Follow-up duration in month"),
       labels = c(""))
  
  # plot <- recordPlot()
  # return(plot)
  # 
  # # return the plot object 
  # grid.echo()
  # Fig <- grid.grab()
  
}





draw_2dose_efficacy_meta <- function(df.obs = df.mr.2d.rct,
                                     df.pred = waning.ve.2dose,
                                     coeff = log.2d.rct$beta){
  
  par(mar=c(4,4,2,4))
  
  df.pred <- df.pred %>% filter(type == "Efficacy")
  length <- nrow(df.pred)
  month_index <- c(which(df.pred$month == 12), which(df.pred$month == 24), which(df.pred$month == 36),
                   which(df.pred$month == 48), which(df.pred$month == 60))
  text_labels <- paste0(round(df.pred$median[month_index]), " [", 
                        round(df.pred$ci.lb[month_index]), ",", 
                        round(df.pred$ci.ub[month_index]), "]")
  
  # waning line
  plot(x = df.pred$month, y = (1-exp(coeff[1] + coeff[2]*df.pred$month))*100,
       xlim = c(0,62),
       ylim=c(-20,100), type="l", lty=1, xaxt = "n",
       ylab="kOCV Efficacy (%)", bty="L",
       xlab="Months after vaccination", col = "black",
       lwd = 1.5
       ) 
  axis(1, at = c(12, 24, 36, 48, 60))
  
  # ribbon for CI
  polygon(c(df.pred$month, rev(df.pred$month)), c(df.pred$ci.lb, rev(df.pred$ci.ub)),
          col = adjustcolor("grey33", 0.2),  border = NA )
  # ribbon for PI
  polygon(c(df.pred$month, rev(df.pred$month)), c(df.pred$pi.lb, rev(df.pred$pi.ub)),
          col = adjustcolor("grey33", 0.2),  border = NA )
  points(x = c(12,24,36,48,60),
         y = df.pred$median[month_index], pch = 19, cex = 0.8, col = "black")
  
  # labels
  text(x = c(12,24,36,48,60)-2,
       y = df.pred$median[month_index]-5,
       labels = text_labels, col = "black")
  
  # estimates segment 
  segments(x0 = df.obs$TL, y0 = df.obs$Mean.VE * 100,
           x1 = df.obs$TR, y1 = df.obs$Mean.VE * 100, col = "grey33")
  abline(h = 0, lty = 2)
  
  title(main = "C", adj = 0, line = 1)
  
  
}


draw_2dose_effectiveness_estimates_colored <- function(fpdata){
  
  fpdata <- fpdata %>% filter(!is.na(yearid))
  fpdata$fPaper <- factor(fpdata$Paper, levels = unique(fpdata$Paper))
  fpdata$uid <- 1:nrow(fpdata)
  fpdata$fuid <- factor(fpdata$uid)
  tau_estim_method <- "EB" # seems like a reasonable choice from DOI: 10.1002/sim.2688
  
  # subset data to observational study only, and order by yearid and continent
  fpdata_obs <- fpdata %>% mutate(se = ifelse(is.na(VE.u),
                                              (log(1-Mean.VE)-log(1-VE.l))/qnorm(0.95),
                                              (log(1-VE.l) - log(1-VE.u))/(2*qnorm(.975))),
                                  yi = log(1-Mean.VE)) %>%
    arrange((yearid)) %>%
    filter(StudyType=="Effectiveness")  %>%
    filter(Dose==2) %>%
    filter(yearid!=0)
  
  # overall effectiveness
  log_ve_fup_obs <- rma(yi=log(1-fpdata_obs$Mean.VE),
                        sei=fpdata_obs$se,
                        method=tau_estim_method)
  
  # estimate yearly pe
  fpdata_obs_1 <- subset(fpdata_obs,yearid==1)
  fpdata_obs_rma_1 <- rma(yi=log(1-fpdata_obs_1$Mean.VE),
                          sei=fpdata_obs_1$se,
                          method=tau_estim_method)
  fpdata_obs_2 <- subset(fpdata_obs,yearid==2)
  fpdata_obs_rma_2 <- rma(yi=log(1-fpdata_obs_2$Mean.VE),
                          sei=fpdata_obs_2$se,
                          method=tau_estim_method)
  fpdata_obs_3 <- subset(fpdata_obs,yearid==3)
  fpdata_obs_rma_3 <- rma(yi=log(1-fpdata_obs_3$Mean.VE),
                          sei=fpdata_obs_3$se,
                          method=tau_estim_method)
  
  
  # add pooled estimates
  ci.lb.obs=c(fpdata_obs_rma_1$ci.lb, fpdata_obs_rma_2$ci.lb, fpdata_obs_rma_3$ci.lb)
  ci.ub.obs=c(fpdata_obs_rma_1$ci.ub, fpdata_obs_rma_2$ci.ub, fpdata_obs_rma_3$ci.ub)
  yi.obs=c(fpdata_obs_rma_1$b, fpdata_obs_rma_2$b, fpdata_obs_rma_3$b)
  
  ci.lb.obs=100*(1-exp(ci.lb.obs))
  ci.ub.obs=100*(1-exp(ci.ub.obs))
  yi.obs=100*(1-exp(yi.obs))
  
  # CI
  b.cr.lb.obs=100*(1-exp(log_ve_fup_obs$yi-qnorm(0.975)*sqrt(log_ve_fup_obs$vi)))
  b.cr.ub.obs=100*(1-exp(log_ve_fup_obs$yi+qnorm(0.975)*sqrt(log_ve_fup_obs$vi)))
  b.cr.obs = 100*(1-exp(log_ve_fup_obs$yi))
  
  ###########################
  # Updated Horizontal Figure with just Effectiveness Study for 2 dose
  ########################
  index_count=table(fpdata_obs$yearid)
  index_center=c(4.5,13.5,22.5,31.5) # mid point of each panel
  index_left=index_center-((index_count+1)/2-0.5)  # left most index of each panel
  index_right=index_center+((index_count+1)/2-0.5)  # right most index of each panel
  # index_left=index_center[1:length(index_count)]-((index_count+1)/2-0.5)  # left most index of each panel
  # index_right=index_center[1:length(index_count)]+((index_count+1)/2-0.5)  # right most index of each panel
  index=NULL
  for (i in 1:4){
    index=c(index,seq(index_left[i],length.out=index_count[i]))
  }
  
  # point size
  wi <- 1/sqrt(log_ve_fup_obs$vi)
  psize <- wi/sum(wi, na.rm = TRUE)
  psize <- (psize - min(psize, na.rm = TRUE))/(max(psize,
                                                   na.rm = TRUE) - min(psize, na.rm = TRUE))
  psize <- (psize * 1) + 0.5
  
  # base plot
  par(mar=c(4,6,6,4))
  plot(1:43,rep(NA,length.out=43),pch=15,
       ylim=c(-50,100),col="grey",
       xaxt = "n",ylab="kOCV Effectiveness (%)",bty="n",
       xlab="")
  
  points(index,100*(1-exp(log_ve_fup_obs$yi)),pch=15,
         cex=psize,col=fpdata_obs$color,
         xaxt = "n",xlab="",ylab="kOCV Efficacy (%)",bty="n")
  
  # add CI
  nonarrowindex=which(b.cr.ub.obs>=(-50))
  
  segments(index-0.3,b.cr.lb.obs,index+0.3,b.cr.lb.obs,col=fpdata_obs$color)
  segments(index[nonarrowindex],b.cr.lb.obs[nonarrowindex],
           index[nonarrowindex],b.cr.ub.obs[nonarrowindex],col=fpdata_obs$color[nonarrowindex])
  segments(index[nonarrowindex]-0.3,b.cr.ub.obs[nonarrowindex],
           index[nonarrowindex]+0.3,b.cr.ub.obs[nonarrowindex],col=fpdata_obs$color[nonarrowindex])
  
  # add CI for studies with large CI
  arrowindex=which(b.cr.ub.obs<(-50))
  segments(index[arrowindex],b.cr.lb.obs[arrowindex],
           index[arrowindex],rep(-50,length.out=length(arrowindex)),col=fpdata_obs$color[arrowindex])
  for (i in 1:length(arrowindex)){
    polygon(y=c(-50,-46,-46,-50),
            x=c(index[arrowindex[i]],index[arrowindex[i]]-0.2 ,
                index[arrowindex[i]]+0.2,index[arrowindex[i]]),
            col=fpdata_obs$color[arrowindex[i]],border=fpdata_obs$color[arrowindex[i]])
  }
  
  # add pooled estimates
  for(i in 1:3){
    polygon(y = c(ci.lb.obs[i], yi.obs[i], ci.ub.obs[i], yi.obs[i]),
            x=c(index_right[i],index_right[i]-0.3,
                index_right[i],index_right[i]+0.3),col="black", border = "black")
  }
  
  # add axis
  axis(side=4, at=c(-50,0,50,100))
  axis(side=1,pos=-70,xlab="",tck=0,col="white",
       at=index_center[1:3],
       labels=c(paste(format(yi.obs,digits=2),"% [",as.numeric(format(ci.ub.obs,digits=2)),
                      ",",as.numeric(format(ci.lb.obs,digits=2)),"]",sep="")))
  axis(side=1,pos=-70,xlab="",tck=0,col.axis="black",col="white",
       at=index_center[4],labels=c("94% [56,99]"))
  axis(side=1,pos=-50,xlab="",tck=0,lty=3,col="white",font=2,
       at=c(4.5,13.5,22.5,31.5, 40.5),
       labels=c("0-12m","12-24m","24-36m","36-48m","48-60m"))
  axis(side=1,pos=-90,xlab="",tck=0,lty=3,col="white",
       at=22.5,
       #labels=c("Follow-up duration in month"),
       labels = c(""))
  axis(side=3,pos=100,xlab="",tck=0,lty=3,at=index_center[3],
       labels=c("Two-dose effectiveness estimates"))
  title(main = "B", adj = 0, line = 1.5)
  abline(h=c(100,0,-50),lty=3,col="grey")
  
  
}


draw_2dose_effectiveness_meta <- function(df.obs = df.mr.2d.obs,
                                          df.pred = waning.ve.2dose,
                                          coeff = log.2d.obs$beta){
  
  par(mar=c(4,6,2,1.5))
  
  df.pred <- df.pred %>% filter(type == "Effectiveness") %>% filter(month <= 48)
  length <- nrow(df.pred)
  month_index <- c(which(df.pred$month == 12), which(df.pred$month == 24), which(df.pred$month == 36),
                   which(df.pred$month == 48))
  text_labels <- paste0(round(df.pred$median[month_index]), " [", 
                        round(df.pred$ci.lb[month_index]), ",", 
                        round(df.pred$ci.ub[month_index]), "]")
  
  # waning line
  plot(x = df.pred$month, y = (1-exp(coeff[1] + coeff[2]*log(df.pred$month)))*100,
       xlim = c(0,62),
       ylim=c(-20,100), type="l", lty=1, xaxt = "n",
       ylab="kOCV Effectiveness (%)", bty="L",
       xlab="Months after vaccination", col = "black",
       lwd = 1.5
       ) 
  axis(1, at = c(12, 24, 36, 48, 60))
  
  # ribbon for CI
  polygon(c(df.pred$month, rev(df.pred$month)), c(df.pred$ci.lb, rev(df.pred$ci.ub)),
          col = adjustcolor("grey33", 0.2),  border = NA )
  # ribbon for PI
  polygon(c(df.pred$month, rev(df.pred$month)), c(df.pred$pi.lb, rev(df.pred$pi.ub)),
          col = adjustcolor("grey33", 0.2),  border = NA )
  points(x = c(12,24,36,48),
         y = df.pred$median[month_index], pch = 19, cex = 0.8, col = "black")
  
  # labels
  text(x = c(12,24,36,48)-2,
       y = df.pred$median[month_index]-5,
       labels = text_labels, col = "black")
  
  # estimates segment 
  segments(x0 = df.obs$TL, y0 = df.obs$Mean.VE * 100,
           x1 = df.obs$TR, y1 = df.obs$Mean.VE * 100, col = "grey33")
  abline(h = 0, lty = 2)
  
  title(main = "D", adj = 0, line = 1)
  
}

# for main analysis (keep all vaccine types, including pre-shanchol vaccines)
draw_2dose_panelPlot_colored <- function(fpdata = df_FollowupPeriod){
  
  fpdata_colored <- assign_colors(fpdata)
  layout(matrix(c(1,1,2,2,2,
                  1,1,2,2,2,
                  3,3,4,4,4,
                  3,3,4,4,4), nrow = 5, ncol = 4, byrow = FALSE))
  draw_2dose_efficacy_estimates_colored(fpdata = fpdata_colored)
  draw_2dose_efficacy_meta()
  draw_2dose_effectiveness_estimates_colored(fpdata = fpdata_colored)
  draw_2dose_effectiveness_meta()
  
}

# for sensitive analysis (removing pre-shanchol vaccines)
draw_2dose_panelPlot_colored_SSA <- function(fpdata = df_FollowupPeriod){
  
  fpdata_colored <- assign_colors(fpdata)
  layout(matrix(c(1,1,2,2,2,
                  1,1,2,2,2,
                  3,3,4,4,4,
                  3,3,4,4,4), nrow = 5, ncol = 4, byrow = FALSE))
  draw_2dose_efficacy_estimates_colored_SSA(fpdata = fpdata_colored)
  draw_2dose_efficacy_meta()
  draw_2dose_effectiveness_estimates_colored(fpdata = fpdata_colored)
  draw_2dose_effectiveness_meta()
  
}


#### one dose 

# combine efficacy and effectiveness estimates
draw_1dose_estimates_colored <- function(fpdata){
  
  fpdata <- fpdata %>% filter(!is.na(yearid))
  fpdata$fPaper <- factor(fpdata$Paper, levels = unique(fpdata$Paper))
  fpdata$uid <- 1:nrow(fpdata)
  fpdata$fuid <- factor(fpdata$uid)
  tau_estim_method <- "EB" # seems like a reasonable choice from DOI: 10.1002/sim.2688
  
  ##########################
  ####### set up data for 1 dose 
  ##########################
  # (include only effectiveness estimates into the pooled estimated / meta analysis)
  temp <- fpdata %>% 
    arrange((yearid)) %>%
    filter(Dose==1) %>%
    # group by every 6 months
    mutate(TM = (TL + TR)/2) %>%
    mutate(monthid = case_when(TM > 0 & TM <= 6 ~ 1,
                               TM > 6 & TM <= 12 ~ 2,
                               TM > 12 & TM <= 18 ~ 3,
                               TM > 18 & TM <= 24 ~ 4,
                               TM > 24 & TM <= 30 ~ 5)) %>%
    filter(yearid!=0) %>% 
    group_by(monthid) %>%
    arrange(StudyType, .by_group = TRUE) %>%
    filter(Paper != "Qadri et al, 2016")
  
  # select only effectiveness studies
  fpdata_1dose <- temp %>% filter(StudyType == "Effectiveness")
  fpdata_1dose_rct <- temp %>% filter(StudyType == "Efficacy")
  
  # overall effectiveness
  log_ve_fup_1dose <- rma(yi=log(1-fpdata_1dose$Mean.VE),
                          sei=fpdata_1dose$se,
                          method=tau_estim_method)
  
  #estimate every 6 months pe
  fpdata_1dose_1 <- subset(fpdata_1dose,monthid==1)
  fpdata_rma_1dose_1 <- rma(yi=log(1-fpdata_1dose_1$Mean.VE),
                            sei=fpdata_1dose_1$se,
                            method=tau_estim_method)
  # only one monthid == 2
  fpdata_1dose_3 <- subset(fpdata_1dose,monthid==3)
  fpdata_rma_1dose_3 <- rma(yi=log(1-fpdata_1dose_3$Mean.VE),
                            sei=fpdata_1dose_3$se,
                            method=tau_estim_method)
  # there is no monthid == 4
  fpdata_1dose_5 <- subset(fpdata_1dose,monthid==5)
  fpdata_rma_1dose_5 <- rma(yi=log(1-fpdata_1dose_5$Mean.VE),
                            sei=fpdata_1dose_5$se,
                            method=tau_estim_method)
  
  
  # CI of effectiveness studies
  b.cr.lb.1dose=100*(1-exp(log_ve_fup_1dose$yi-qnorm(0.975)*sqrt(log_ve_fup_1dose$vi)))
  b.cr.ub.1dose=100*(1-exp(log_ve_fup_1dose$yi+qnorm(0.975)*sqrt(log_ve_fup_1dose$vi)))
  b.cr.1dose = 100*(1-exp(log_ve_fup_1dose$yi))
  
  # insert in the efficacy estimates into the effectiveness estimates
  b.cr.lb.1dose <- c(b.cr.lb.1dose[1:4], (fpdata_1dose_rct$VE.u*100)[1], b.cr.lb.1dose[5], (fpdata_1dose_rct$VE.u*100)[2], b.cr.lb.1dose[6:8], (fpdata_1dose_rct$VE.u*100)[3:4], b.cr.lb.1dose[9:10] )
  b.cr.ub.1dose <- c(b.cr.ub.1dose[1:4], (fpdata_1dose_rct$VE.l*100)[1], b.cr.ub.1dose[5], (fpdata_1dose_rct$VE.l*100)[2], b.cr.ub.1dose[6:8], (fpdata_1dose_rct$VE.l*100)[3:4], b.cr.ub.1dose[9:10] )
  b.cr.1dose <- c(b.cr.1dose[1:4], (fpdata_1dose_rct$Mean.VE*100)[1], b.cr.1dose[5], (fpdata_1dose_rct$Mean.VE*100)[2], b.cr.1dose[6:8], (fpdata_1dose_rct$Mean.VE*100)[3:4], b.cr.1dose[9:10])
  
  #pooled estimates
  ci.lb.1dose=c(fpdata_rma_1dose_1$ci.lb, fpdata_rma_1dose_3$ci.lb, fpdata_rma_1dose_5$ci.lb)
  ci.ub.1dose=c(fpdata_rma_1dose_1$ci.ub, fpdata_rma_1dose_3$ci.ub, fpdata_rma_1dose_5$ci.ub)
  yi.1dose=c(fpdata_rma_1dose_1$b, fpdata_rma_1dose_3$b, fpdata_rma_1dose_5$b)
  
  ci.lb.1dose=100*(1-exp(ci.lb.1dose))
  ci.ub.1dose=100*(1-exp(ci.ub.1dose))
  yi.1dose=100*(1-exp(yi.1dose))
  
  ###########################
  # Updated Horizontal Figure with all estimates for 1 dose
  ########################
  index_count=table(temp$monthid)
  index_center=c(4.5,13.5,22.5,31.5, 40.5) # mid point of each panel
  index_left=index_center-((index_count+1)/2-0.5)  # left most index of each panel
  index_right=index_center+((index_count+1)/2-0.5)  # right most index of each panel
  index=NULL
  for (i in 1:length(index_count)){
    index=c(index,seq(index_left[i],length.out=index_count[i]))
  }
  
  # point size
  wi <- 1/sqrt(log_ve_fup_1dose$vi)
  psize <- wi/sum(wi, na.rm = TRUE)
  psize <- (psize - min(psize, na.rm = TRUE))/(max(psize,na.rm = TRUE) - min(psize, na.rm = TRUE))
  psize <- (psize * 1) + 0.5
  
  # index of effectiveness studies
  obs_position <- c(index[1:4], index[6], index[8:10], index[13:14])
  rct_position <- c(index[5], index[7], index[11:12] )
  obs_index <- c(1:4, 6, 8:10, 13:14)
  rct_index <- c(5,7,11,12)
  
  # base plot
  par(mar=c(6,6,3,5))
  plot(1:43,rep(NA,length.out=43),pch=15,
       ylim=c(-50,100),col="grey",
       xaxt = "n",ylab="kOCV Efficacy/Effectiveness (%)",bty="n",
       xlab="")
  
  # points for effectiveness studies
  points(obs_position,100*(1-exp(log_ve_fup_1dose$yi)),pch=15,
         cex=psize,col=fpdata_1dose$color,
         xaxt = "n",xlab="",ylab="kOCV Efficacy/Effectiveness (%)",bty="n")
  
  # points for efficacy studies
  points(rct_position, fpdata_1dose_rct$Mean.VE*100,pch=20,
         cex=1,col=fpdata_1dose_rct$color,
         xaxt = "n",xlab="",ylab="kOCV Efficacy/Effectiveness (%)",bty="n")
  
  
  # add CI 
  nonarrowindex= c(1,2,4,5,6,7,8,9,11,12,13)
  
  # draw CI for efficacy studies in dashed lines
  segments(index[rct_index]-0.3,b.cr.lb.1dose[rct_index],
           index[rct_index]+0.3,b.cr.lb.1dose[rct_index],col=fpdata_1dose_rct$color)
  segments(index[rct_index],b.cr.lb.1dose[rct_index],
           index[rct_index],b.cr.ub.1dose[rct_index],col=fpdata_1dose_rct$color, lty = 2)
  segments(index[rct_index]-0.3,b.cr.ub.1dose[rct_index],
           index[rct_index]+0.3,b.cr.ub.1dose[rct_index],col=fpdata_1dose_rct$color)
  
  # draw CI for effectiveness studies in solid lines
  obs_nonarrowindex <- c(1,2,4,6,8,9,13)
  segments(index[obs_nonarrowindex]-0.3,b.cr.lb.1dose[obs_nonarrowindex], 
           index[obs_nonarrowindex]+0.3,b.cr.lb.1dose[obs_nonarrowindex],
           col=temp$color[obs_nonarrowindex])
  segments(index[obs_nonarrowindex],b.cr.lb.1dose[obs_nonarrowindex],
           index[obs_nonarrowindex],b.cr.ub.1dose[obs_nonarrowindex],col=temp$color[obs_nonarrowindex])
  segments(index[obs_nonarrowindex]-0.3,b.cr.ub.1dose[obs_nonarrowindex],
           index[obs_nonarrowindex]+0.3,b.cr.ub.1dose[obs_nonarrowindex],col=temp$color[obs_nonarrowindex])
  
  # add CI for studies with large CI
  arrowindex=c(3,10,14)
  segments(index[arrowindex]-0.3,b.cr.lb.1dose[arrowindex], index[arrowindex]+0.3,b.cr.lb.1dose[arrowindex],col=temp$color[arrowindex])
  segments(index[arrowindex],b.cr.lb.1dose[arrowindex],
           index[arrowindex],rep(-50,length.out=length(arrowindex)),col= temp$color[arrowindex])
  for (i in 1:length(arrowindex)){
    polygon(y=c(-50,-46,-46,-50),
            x=c(index[arrowindex[i]],index[arrowindex[i]]-0.2 ,
                index[arrowindex[i]]+0.2,index[arrowindex[i]]),
            col=temp$color[arrowindex[i]],border=temp$color[arrowindex[i]])
  }
  
  # add pooled estimates
  for(i in c(1,3,5)){
    polygon(y = c(ci.lb.1dose[ceiling(i/2)], yi.1dose[ceiling(i/2)], ci.ub.1dose[ceiling(i/2)], yi.1dose[ceiling(i/2)]),
            x=c(index_right[i],index_right[i]-0.3,
                index_right[i],index_right[i]+0.3),col="black", border = "black")
  }
  
  # add axis
  axis(side=4, at=c(-50,0,50,100))
  
  # label effectiveness studies
  axis(side=1,pos=-70,xlab="",tck=0,col="white", col.axis="black",
       at=index_center[1],
       labels=c(paste(format(yi.1dose[[1]],digits=2),"% [",
                      as.numeric(format(ci.ub.1dose[[1]],digits=2)),
                      ",",
                      as.numeric(format(ci.lb.1dose[[1]],digits=2)),"]",sep="")))
  axis(side=1,pos=-70,xlab="",tck=0,col="white", col.axis="black",
       at=index_center[3],
       labels=c(paste(format(yi.1dose[[2]],digits=2),"% [",
                      as.numeric(format(ci.ub.1dose[[2]],digits=2)),
                      ",",
                      as.numeric(format(ci.lb.1dose[[2]],digits=2)),"]",sep="")))
  axis(side=1,pos=-70,xlab="",tck=0,col="white", col.axis="black",
       at=index_center[5],
       labels=c(paste(format(yi.1dose[[3]],digits=2),"% [",
                      as.numeric(format(ci.ub.1dose[[3]],digits=2)),
                      ",",
                      as.numeric(format(ci.lb.1dose[[3]],digits=2)),"]",sep="")))
  axis(side=1,pos=-70,xlab="",tck=0,col="white", col.axis="black",
       at=index_center[2],
       labels="92% [66,98]")
  
  # label efficacy studies
  axis(side=1,pos=-80,xlab="",tck=0,col="white", col.axis="grey61",
       at=index_center[1],
       labels="58% [24,76]")
  axis(side=1,pos=-80,xlab="",tck=0,col="white", col.axis="grey61",
       at=index_center[2],
       labels="37% [-20,67]")
  axis(side=1,pos=-80,xlab="",tck=0,col="white", col.axis="grey61",
       at=index_center[3],
       labels="62% [34,78]")
  axis(side=1,pos=-80,xlab="",tck=0,col="white", col.axis="grey61",
       at=index_center[4],
       labels="67% [30,84]")
  
  
  axis(side=1,pos=-50,xlab="",tck=0,lty=3,col="white",font=2,
       at=c(4.5,13.5,22.5,31.5,40.5), # mid point of each panel
       labels=c("0-6m","6-12m","12-18m","18-24m","24-30m"))
  abline(h=c(100,0,-50),lty=3,col="grey")
  #abline(h=c(0,50),lty=3,col=AddAlpha("grey",0.4))
  axis(side=3,pos=110,xlab="",tck=0,lty=3,at=22.5,
       labels=c("One-dose estimates"))
  axis(side=1,pos=-100,xlab="",tck=0,lty=3,col="white",
       at=22.5,
       #labels=c("Follow-up duration in month"),
       labels = "")
  title(main = "A", adj = 0, line = 1.5)
  
  
}


draw_1dose_effectiveness_meta <- function(df.obs = df.mr.1d.obs,
                                          df.pred = pred.est.1d %>% filter(month_6m <= 6.2),
                                          coeff = log.1d$beta){
  
  par(mar=c(4,6,2,1.5))
  
  length <- nrow(df.pred)
  # month_index <- c(which(df.pred$month_6m == 1), which(df.pred$month_6m == 2), which(df.pred$month_6m == 3),
  #                  which(df.pred$month_6m == 4), which(df.pred$month_6m == 5))
  month_index <- c(1,2,3,4,5)*6000
  text_labels <- paste0(round(df.pred$median[month_index]), " [", 
                        round(df.pred$ci.lb[month_index]), ",", 
                        round(df.pred$ci.ub[month_index]), "]")
  
  # waning line
  plot(x = df.pred$month_6m*6, y = (1-exp(coeff[1] + coeff[2]*log(df.pred$month_6m * 6)))*100,
       xlim = c(0,32),
       ylim=c(-20,100), type="l", lty=1, xaxt = "n",
       ylab="kOCV Effectiveness (%)", bty="L",
       xlab="Months after vaccination", col = "black",
       lwd = 1.5
       ) 
  axis(1, at = c(6, 12, 18, 24, 30))
  
  # estimates segment (effectiveness studies)
  segments(x0 = df.obs$TL, y0 = df.obs$Mean.VE * 100,
           x1 = df.obs$TR, y1 = df.obs$Mean.VE * 100, col = "grey33")
  # estimates segment (efficacy studies)
  segments(x0 = c(0,6,12,18), y0 = c(58,37,62,67),
           x1 = c(6,12,18,24), y1 = c(58,37,62,67), col = "grey33", lty = 2) # from qadri 2018
  
  
  abline(h = 0, lty = 2)
  
  # ribbon for CI
  polygon(c(df.pred$month_6m * 6, rev(df.pred$month_6m * 6)), c(df.pred$ci.lb, rev(df.pred$ci.ub)),
          col = adjustcolor("grey33", 0.1),  border = NA )
  # ribbon for PI
  polygon(c(df.pred$month_6m * 6, rev(df.pred$month_6m * 6)), c(df.pred$pi.lb, rev(df.pred$pi.ub)),
          col = adjustcolor("grey33", 0.1),  border = NA )
  points(x = c(6,12,18,24,30),
         y = df.pred$median[month_index], pch = 19, cex = 0.8, col = "black")
  
  # labels
  text(x = c(6,12,18,24,30)-2,
       y = df.pred$median[month_index]-5,
       labels = text_labels, col = "black")
  
  
  
  title(main = "B", adj = 0, line = 1)
  
}


draw_1dose_panelPlot_colored <- function(fpdata = df_FollowupPeriod){
 
  fpdata_colored <- assign_colors(fpdata)
  
   layout(matrix(c(1,1,2,2,2), ncol = 1, byrow = FALSE))
  
  draw_1dose_estimates_colored(fpdata = fpdata_colored)
  draw_1dose_effectiveness_meta()
}






