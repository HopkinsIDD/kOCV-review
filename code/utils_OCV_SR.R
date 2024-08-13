# Functions to be used in updated OCV review

#################################################
# Function to make simple forest plot of VE
#################################################
ve_forest <- function(dat,
                      dose = 2,
                      vaccine_type = c("WC", "WC-rBS"),
                      tau_estim_method = "EB"
                      ){
  
  par(mar=c(4,4,1,2))
  
  # subset data rows
  dat <- dat %>% filter(Doses == dose, Vaccine %in% vaccine_type) %>%
    arrange(desc(StudyDesign),desc(Duration.months))
  
  
  log_ve <- rma(yi=yi,
                sei=se, method=tau_estim_method,
                slab=paste0(dat$Location,", ",dat$Paper," (",dat$Duration.months," months)"),
                dat=dat)
  
  # set positions
  rct.l <- 2
  rct.u <- rct.l + dat %>% filter(StudyDesign == "Randomized") %>% nrow() - 1
  obs.l <- rct.u + 4
  obs.u <- obs.l + dat %>% filter(StudyDesign == "Observational") %>% nrow() - 1
  ylim.u <- obs.u + 4
  
  forest(log_ve,
         slab=dat$Paper,
         transf=function(x) 1-exp(x),
         refline=NA,
         xlab=paste0(dose, "-dose Efficacy/Effectivness", " (", paste0(vaccine_type, collapse = " or "), ")"),
         ilab=data.frame(dat$Location,dat$Duration.months),
         ilab.xpos=c(-1.7,-.75), cex=.75,
         addcred=F,
         ylim=c(0.2,ylim.u),
         xlim=c(-4.5,3),
         rows=c(rct.l:rct.u, obs.l:obs.u),
         alim=c(-.5,1),
         order=1:nrow(dat))
  
  # add text for RCTs and Observational
  text(-4.5, c(rct.u + 1, obs.u + 1), pos=4,  c("Efficacy","Effectivness"), cex=.8, font=2)
  text(c(-1.7,-0.75), obs.u + 1.6, c("Location","Duration"),cex=.8)
  text(2.8, obs.u + 1.6, " VE [95% CI]", pos=2, cex=0.8)
  
  # add summary stats for rct/obs
  effectivness_only <- subset(dat,StudyDesign=='Observational')
  effectivness_rma <- rma(yi=yi,
                          sei=se,
                          method=tau_estim_method,
                          dat=dat,
                          subset=StudyDesign=='Observational')
  
  efficacy_only <- subset(dat,StudyDesign=='Randomized')
  effiacy_rma <- rma(yi=yi,
                     sei=se,
                     method=tau_estim_method,
                     dat=dat,
                     subset=StudyDesign=='Randomized')
  
  addpoly(effectivness_rma, row= obs.l - 1, cex=.75, transf=function(x) 1-exp(x), mlab="Mean Effectivness",  font = 4)
  addpoly(effiacy_rma, row= rct.l - 1, cex=.75,transf=function(x) 1-exp(x), mlab="Mean Efficacy", font =4)
  text(-0.75, obs.l-0.3,round(weighted.mean(effectivness_only$Duration.months,
                                      w=1/effectivness_only$se^2),0),pos=1,cex=0.75)
  text(-0.75, rct.l-0.4,round(weighted.mean(efficacy_only$Duration.months,
                                      w=1/efficacy_only$se^2),0),pos=1,cex=0.75)
  
}

ve_forest_byage <- function(dat,
                      dose = 2,
                      vaccine_type = c("WC", "WC-rBS"),
                      tau_estim_method = "EB"
){
  
  par(mar=c(4,4,1,2))
  
  # subset data rows
  dat <- dat %>% filter(Doses == dose, Vaccine %in% vaccine_type) %>%
    arrange((age.group),desc(Duration.months))
  
  
  log_ve <- rma(yi=yi,
                sei=se, method=tau_estim_method,
                slab=paste0(dat$Location,", ",dat$Paper," (",dat$Duration.months," months)"),
                dat=dat)
  
  # set positions
  o5a.l <- 2
  o5a.u <- o5a.l + dat %>% filter(age.group == "Above 5 yo") %>% nrow() - 1
  u5a.l <- o5a.u + 4
  u5a.u <- u5a.l + dat %>% filter(age.group == "Under 5 yo") %>% nrow() - 1
  ylim.u <- u5a.u + 4.5
  
  if(dose == 2){
    dat <- dat %>%
      mutate(StudyID = case_when(StudyID == "Sur et al, 2011" & Mean.VE == 0.88 ~ "Sur et al, 2011a",
                                 StudyID == "Sur et al, 2011" & Mean.VE == 0.61 ~ "Sur et al, 2011b",
                                 StudyID == "Qadri et al, 2015" & Mean.VE == 0.33 ~ "Qadri et al, 2015a",
                                 StudyID == "Qadri et al, 2015" & Mean.VE == 0.56 ~ "Qadri et al, 2015b",
                                 TRUE ~ StudyID))
  }
  if(dose == 1){
    dat <- dat %>%
      mutate(StudyID = case_when(StudyID == "Qadri et al, 2016" & Mean.VE == 0.63 ~ "Qadri et al, 2016a",
                                 StudyID == "Qadri et al, 2016" & Mean.VE == 0.56 ~ "Qadri et al, 2016b",
                                 StudyID == "Qadri et al, 2018" & Mean.VE == 0.52 ~ "Qadri et al, 2018a",
                                 StudyID == "Qadri et al, 2018" & Mean.VE == 0.59 ~ "Qadri et al, 2018b",
                                 TRUE ~ StudyID))
  }
  
  forest(log_ve,
         slab=dat$StudyID,
         transf=function(x) 1-exp(x),
         refline=NA,
         xlab=paste0(dose, "-dose Efficacy", " (", paste0(vaccine_type, collapse = " or "), ")"),
         ilab=data.frame(dat$Location,dat$Duration.months),
         ilab.xpos=c(-1.7,-.75), cex=.75,
         addcred=F,
         ylim=c(0.2,ylim.u),
         xlim=c(-4.5,3),
         rows=c(o5a.l:o5a.u, u5a.l:u5a.u),
         alim=c(-.5,1),
         order=1:nrow(dat))
  
  # add text for RCTs and Observational
  text(-4.5, c(o5a.u + 1, u5a.u + 1), pos=4,  c("5 and older","Under 5"), cex=.8, font=2)
  text(c(-1.7,-0.75), u5a.u + 1.6, c("Location","Duration"),cex=.8)
  text(2.8, u5a.u + 1.6, " VE [95% CI]", pos=2, cex=0.8)
  
  # add summary stats for rct/obs
  u5a_only <- subset(dat,age.group == "Under 5 yo")
  u5a_rma <- rma(yi=yi,
                          sei=se,
                          method=tau_estim_method,
                          dat=dat,
                          subset=age.group == "Under 5 yo")
  
  o5a_only <- subset(dat,age.group == "Above 5 yo")
  o5a_rma <- rma(yi=yi,
                     sei=se,
                     method=tau_estim_method,
                     dat=dat,
                     subset=age.group == "Above 5 yo")
  
  addpoly(u5a_rma, row= u5a.l - 1, cex=.75, transf=function(x) 1-exp(x), mlab="Mean Efficacy",  col = "steelblue", border = "steelblue")
  addpoly(o5a_rma, row= o5a.l - 1, cex=.75,transf=function(x) 1-exp(x), mlab="Mean Efficacy", col = "steelblue",  border = "steelblue")
  addpoly(u5a_rma, row= u5a.l - 1, cex=.75, transf=function(x) 1-exp(x), mlab="Mean Efficacy",  col = "steelblue",  border = "steelblue")
  addpoly(o5a_rma, row= o5a.l - 1, cex=.75,transf=function(x) 1-exp(x), mlab="Mean Efficacy", col = "steelblue", border = "steelblue")
  
  text(-0.75, u5a.l-0.3,round(weighted.mean(u5a_only$Duration.months,
                                            w=1/u5a_only$se^2),0),pos=1,cex=0.75)
  text(-0.75, o5a.l-0.4,round(weighted.mean(o5a_only$Duration.months,
                                            w=1/o5a_only$se^2),0),pos=1,cex=0.75)
  
}


#################################################
# Function to make 4-panel plot of VE over time (1 dose only figure, combining efficacy and effectiveness)
#################################################
make_fp_ve_1dose <- function(fpdata = df_FollowupPeriod){
  
  fpdata <- fpdata %>% filter(!is.na(yearid))
  fpdata$fPaper <- factor(fpdata$Paper, levels = unique(fpdata$Paper))
  fpdata$uid <- 1:nrow(fpdata)
  fpdata$fuid <- factor(fpdata$uid)
  tau_estim_method <- "EB" # seems like a reasonable choice from DOI: 10.1002/sim.2688
  
  par(mfrow=c(1,1),mar=c(4,7,7,7))
  
  ##########################
  ####### set up data for 1 dose 
  ##########################
  # (include only effectiveness estimates into the pooled estimated / meta analysis)
  temp <- fpdata %>% 
    arrange((yearid)) %>%
    filter(Dose==1) %>%
    filter(yearid!=0) %>% 
    group_by(yearid) %>%
    arrange(StudyType, .by_group = TRUE)
  
  # select only effectiveness studies
  fpdata_1dose <- temp %>% filter(StudyType == "Effectiveness")
  fpdata_1dose_rct <- temp %>% filter(StudyType == "Efficacy")
  
  # overall effectiveness
  log_ve_fup_1dose <- rma(yi=log(1-fpdata_1dose$Mean.VE),
                              sei=fpdata_1dose$se,
                              method=tau_estim_method)
  
  #estimate yearly pe
  fpdata_1dose_1 <- subset(fpdata_1dose,yearid==1)
  fpdata_rma_1dose_1 <- rma(yi=log(1-fpdata_1dose_1$Mean.VE),
                                sei=fpdata_1dose_1$se,
                                method=tau_estim_method)
  fpdata_1dose_2 <- subset(fpdata_1dose,yearid==2)
  fpdata_rma_1dose_2 <- rma(yi=log(1-fpdata_1dose_2$Mean.VE),
                                sei=fpdata_1dose_2$se,
                                method=tau_estim_method)
  fpdata_1dose_3 <- subset(fpdata_1dose,yearid==3)
  fpdata_rma_1dose_3 <- rma(yi=log(1-fpdata_1dose_3$Mean.VE),
                            sei=fpdata_1dose_3$se,
                            method=tau_estim_method)
  
  
  # CI of effectiveness studies
  b.cr.lb.1dose=100*(1-exp(log_ve_fup_1dose$yi-qnorm(0.975)*sqrt(log_ve_fup_1dose$vi)))
  b.cr.ub.1dose=100*(1-exp(log_ve_fup_1dose$yi+qnorm(0.975)*sqrt(log_ve_fup_1dose$vi)))
  b.cr.1dose = 100*(1-exp(log_ve_fup_1dose$yi))
  
  # insert in the efficacy estimates into the effectiveness estimates
  b.cr.lb.1dose <- c(b.cr.lb.1dose[1:5], (fpdata_1dose_rct$VE.u*100)[1], b.cr.lb.1dose[6:8], (fpdata_1dose_rct$VE.u*100)[2], b.cr.lb.1dose[9:10])
  b.cr.ub.1dose <- c(b.cr.ub.1dose[1:5], (fpdata_1dose_rct$VE.l*100)[1], b.cr.ub.1dose[6:8], (fpdata_1dose_rct$VE.l*100)[2], b.cr.ub.1dose[9:10])
  b.cr.1dose <- c(b.cr.1dose[1:5], (fpdata_1dose_rct$Mean.VE*100)[1], b.cr.1dose[6:8], (fpdata_1dose_rct$Mean.VE*100)[2], b.cr.1dose[9:10])
  
  #pooled estimates
  ci.lb.1dose=c(fpdata_rma_1dose_1$ci.lb, fpdata_rma_1dose_2$ci.lb, fpdata_rma_1dose_3$ci.lb)
  ci.ub.1dose=c(fpdata_rma_1dose_1$ci.ub, fpdata_rma_1dose_2$ci.ub, fpdata_rma_1dose_3$ci.ub)
  yi.1dose=c(fpdata_rma_1dose_1$b, fpdata_rma_1dose_2$b, fpdata_rma_1dose_3$b)

  ci.lb.1dose=100*(1-exp(ci.lb.1dose))
  ci.ub.1dose=100*(1-exp(ci.ub.1dose))
  yi.1dose=100*(1-exp(yi.1dose))
  
  ###########################
  # Updated Horizontal Figure with all estimates for 1 dose
  ########################
  index_count=table(temp$yearid)
  # index_center=c(4.5,13.5,22.5,31.5,40.5) # mid point of each panel
  # index_left=index_center[length(index_count)]-((index_count+1)/2-0.5)  # left most index of each panel
  # index_right=index_center[length(index_count)]+((index_count+1)/2-0.5)  # right most index of each panel
  index_center=c(4.5,13.5,22.5) # mid point of each panel
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
  obsindex <- c(index[1:5], index[7:9], index[11:12])
  rctindex <- c(index[6], index[10])
  
  # base plot
  par(mar=c(6,6,3,6))
  plot(1:43,rep(NA,length.out=43),pch=15,
       ylim=c(-50,100),col="grey",
       xaxt = "n",ylab="kOCV Efficacy/Effectiveness (%)",bty="n",
       xlab="")
  
  # points for effectiveness studies
  points(obsindex,100*(1-exp(log_ve_fup_1dose$yi)),pch=15,
         cex=psize,col="grey",
         xaxt = "n",xlab="",ylab="kOCV Efficacy/Effectiveness (%)",bty="n")
  
  # points for efficacy studies
  points(rctindex, fpdata_1dose_rct$Mean.VE*100,pch=20,
         cex=1,col="grey",
         xaxt = "n",xlab="",ylab="kOCV Efficacy/Effectiveness (%)",bty="n")
  
  
  # add CI 
  nonarrowindex= c(1,2,3,5,6,8,10,11)
  
  segments(index-0.3,b.cr.lb.1dose,index+0.3,b.cr.lb.1dose,col="grey")
  segments(index[nonarrowindex],b.cr.lb.1dose[nonarrowindex],
           index[nonarrowindex],b.cr.ub.1dose[nonarrowindex],col="grey")
  segments(index[nonarrowindex]-0.3,b.cr.ub.1dose[nonarrowindex],
           index[nonarrowindex]+0.3,b.cr.ub.1dose[nonarrowindex],col="grey")
  
  # change line color of effectiveness to black, keep the efficacy estimate as grey
  obs_nonarrowindex <- c(1,2,3,5,8,11)
  points(index[obs_nonarrowindex], b.cr.1dose[obs_nonarrowindex], pch=15,
         cex=psize[c(1,2,3,5,7,9)],col="black",
         xaxt = "n",xlab="",ylab="kOCV Efficacy/Effectiveness (%)",bty="n")
  segments(index[obs_nonarrowindex]-0.3,b.cr.lb.1dose[obs_nonarrowindex], index[obs_nonarrowindex]+0.3,b.cr.lb.1dose[obs_nonarrowindex],col="black")
  segments(index[obs_nonarrowindex],b.cr.lb.1dose[obs_nonarrowindex],
           index[obs_nonarrowindex],b.cr.ub.1dose[obs_nonarrowindex],col="black")
  segments(index[obs_nonarrowindex]-0.3,b.cr.ub.1dose[obs_nonarrowindex],
           index[obs_nonarrowindex]+0.3,b.cr.ub.1dose[obs_nonarrowindex],col="black")
  
  # add CI for studies with large CI
  arrowindex=c(4,7,9,12)
  points(index[arrowindex], b.cr.1dose[arrowindex], pch=15,
         cex=psize[c(4,6,8,10)],col="black",
         xaxt = "n",xlab="",ylab="kOCV Efficacy/Effectiveness (%)",bty="n")
  segments(index[arrowindex]-0.3,b.cr.lb.1dose[arrowindex], index[arrowindex]+0.3,b.cr.lb.1dose[arrowindex],col="black")
  segments(index[arrowindex],b.cr.lb.1dose[arrowindex],
           index[arrowindex],rep(-50,length.out=length(arrowindex)),col="black")
  for (i in 1:length(arrowindex)){
    polygon(y=c(-50,-46,-46,-50),
            x=c(index[arrowindex[i]],index[arrowindex[i]]-0.2 ,
                index[arrowindex[i]]+0.2,index[arrowindex[i]]),
            col="black",border="black")
  }

  # add pooled estimates
  for(i in 1:3){
    polygon(y = c(ci.lb.1dose[i], yi.1dose[i], ci.ub.1dose[i], yi.1dose[i]),
            x=c(index_right[i],index_right[i]-0.3,
                index_right[i],index_right[i]+0.3),col="steelblue", border = "steelblue")
  }
  
  # add axis
  axis(side=4, at=c(-50,0,50,100))
  # axis(side=3,pos=100,xlab="",tck=0,col="white",
  #      at=index_center[1],
  #      labels=c(paste(format(yi.obs,digits=2),"% [",as.numeric(format(ci.ub.obs,digits=2)),
  #                     ",",as.numeric(format(ci.lb.obs,digits=2)),"]",sep="")))
  axis(side=1,pos=-70,xlab="",tck=0,col="white", col.axis="black",
       at=index_center[1],
       labels=c(paste(format(yi.1dose[[1]],digits=2),"% [",
                      as.numeric(format(ci.ub.1dose[[1]],digits=2)),
                      ",",
                      as.numeric(format(ci.lb.1dose[[1]],digits=2)),"]",sep="")))
  axis(side=1,pos=-70,xlab="",tck=0,col="white", col.axis="black",
       at=index_center[2],
       labels=c(paste(format(yi.1dose[[2]],digits=2),"% [",
                      as.numeric(format(ci.ub.1dose[[2]],digits=2)),
                      ",",
                      as.numeric(format(ci.lb.1dose[[2]],digits=2)),"]",sep="")))
  axis(side=1,pos=-70,xlab="",tck=0,col="white", col.axis="black",
       at=index_center[3],
       labels=c(paste(format(yi.1dose[[3]],digits=2),"% [",
                      as.numeric(format(ci.ub.1dose[[3]],digits=2)),
                      ",",
                      as.numeric(format(ci.lb.1dose[[3]],digits=2)),"]",sep="")))
  # axis(side=1,pos=-70,xlab="",tck=0,col="white", col.axis="dimgrey",
  #      at=index_center[3],
  #      labels="33% [-318,89]")
  
  # axis(side=3,pos=100,xlab="",tck=0,col="white",
  #      at=index_center[1],labels=c("45% [26,59]"))
  # axis(side=3,pos=100,xlab="", col.axis = "dimgray",tck=0,col="white",
  #      at=index_center[2],labels=c("64% [43,77]"))
  axis(side=1,pos=-50,xlab="",tck=0,lty=3,col="white",font=2,
       at=c(4.5,13.5,22.5,31.5,40.5), # mid point of each panel
       labels=c("0-12m","12-24m","24-36m","36-48m","48-60m"))
  abline(h=c(100,0,-50),lty=3,col="grey")
  #abline(h=c(0,50),lty=3,col=AddAlpha("grey",0.4))
  axis(side=3,pos=110,xlab="",tck=0,lty=3,at=22.5,
       labels=c("One-dose estimates"))
  axis(side=1,pos=-100,xlab="",tck=0,lty=3,col="white",
       at=22.5,
       labels=c("Follow-up duration in month"))
  title(main = "A", adj = 0, line = 2)
  # add estimate for each study
  # axis(side=3, pos=90,tck=0,col="white",las=2, cex.axis=0.6,
  #      at=index,labels= paste0(
  #        round(b.cr.rct.1dose), " [", round(b.cr.ub.rct.1dose), ", ", round(b.cr.lb.rct.1dose), "]"
  #      ))
  
  # add label for each study
  # axis(side=3, pos=90,tck=0,col="white",las=2, cex.axis=0.5,
  #      at=index,labels= fpdata_rct_1dose$Label)
  
  # change the label color of new studies to steelblue
  # axis(side=3, pos=90,tck=0,col="white",las=2, cex.axis=0.5,col.axis="black",
  #      at=index[which(fpdata_rct_1dose$new.old == "old")],labels=fpdata_rct_1dose$Paper[which(fpdata_rct_1dose$new.old == "old")]
  # )
  # axis(side=3, pos=90,tck=0,col="white",las=2, cex.axis=0.5,col.axis="blue",
  #      at=index[which(fpdata_rct_1dose$new.old == "new")],labels=fpdata_rct_1dose$Paper[which(fpdata_rct_1dose$new.old == "new")]
  # )
}

#################################################
# Function to make 4-panel plot of VE over time (1 dose only figure, combining efficacy and effectiveness), but grouped every 6 months
#################################################
make_fp_ve_1dose_6m <- function(fpdata = df_FollowupPeriod){
  
  fpdata <- fpdata %>% filter(!is.na(yearid))
  fpdata$fPaper <- factor(fpdata$Paper, levels = unique(fpdata$Paper))
  fpdata$uid <- 1:nrow(fpdata)
  fpdata$fuid <- factor(fpdata$uid)
  tau_estim_method <- "EB" # seems like a reasonable choice from DOI: 10.1002/sim.2688
  
  par(mfrow=c(1,1),mar=c(4,7,7,7))
  
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
    arrange(StudyType, .by_group = TRUE)
  
  # remove the duplicate estimate: Qadri 2016 (0-6 months)
  temp <- temp %>%
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
  b.cr.lb.1dose <- c(b.cr.lb.1dose[1:4], (fpdata_1dose_rct$VE.u*100)[1], b.cr.lb.1dose[5], (fpdata_1dose_rct$VE.u*100)[2], b.cr.lb.1dose[6:8], (fpdata_1dose_rct$VE.u*100)[3], (fpdata_1dose_rct$VE.u*100)[4], b.cr.lb.1dose[9:10])
  b.cr.ub.1dose <- c(b.cr.ub.1dose[1:4], (fpdata_1dose_rct$VE.l*100)[1], b.cr.ub.1dose[5], (fpdata_1dose_rct$VE.l*100)[2], b.cr.ub.1dose[6:8], (fpdata_1dose_rct$VE.l*100)[3], (fpdata_1dose_rct$VE.l*100)[4], b.cr.ub.1dose[9:10])
  b.cr.1dose <- c(b.cr.1dose[1:4], (fpdata_1dose_rct$Mean.VE*100)[1], b.cr.1dose[5], (fpdata_1dose_rct$Mean.VE*100)[2], b.cr.1dose[6:8], (fpdata_1dose_rct$Mean.VE*100)[3], (fpdata_1dose_rct$Mean.VE*100)[4], b.cr.1dose[9:10])
  
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
  # index_center=c(4.5,13.5,22.5,31.5,40.5) # mid point of each panel
  # index_left=index_center[length(index_count)]-((index_count+1)/2-0.5)  # left most index of each panel
  # index_right=index_center[length(index_count)]+((index_count+1)/2-0.5)  # right most index of each panel
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
  obsindex <- c(index[1:4], index[6], index[8:10], index[13:14])
  rctindex <- c(index[5], index[7], index[11], index[12])
  
  # base plot
  par(mar=c(6,6,3,6))
  plot(1:43,rep(NA,length.out=43),pch=15,
       ylim=c(-50,100),col="grey",
       xaxt = "n",ylab="kOCV Efficacy/Effectiveness (%)",bty="n",
       xlab="")
  
  # points for effectiveness studies
  points(obsindex,100*(1-exp(log_ve_fup_1dose$yi)),pch=15,
         cex=psize,col="grey",
         xaxt = "n",xlab="",ylab="kOCV Efficacy/Effectiveness (%)",bty="n")
  
  # points for efficacy studies
  points(rctindex, fpdata_1dose_rct$Mean.VE*100,pch=20,
         cex=1,col="grey",
         xaxt = "n",xlab="",ylab="kOCV Efficacy/Effectiveness (%)",bty="n")
  
  
  # add CI 
  nonarrowindex= c(1,2,4,5,6,7,8,9,11,12,13)
  
  segments(index-0.3,b.cr.lb.1dose,index+0.3,b.cr.lb.1dose,col="grey")
  segments(index[nonarrowindex],b.cr.lb.1dose[nonarrowindex],
           index[nonarrowindex],b.cr.ub.1dose[nonarrowindex],col="grey")
  segments(index[nonarrowindex]-0.3,b.cr.ub.1dose[nonarrowindex],
           index[nonarrowindex]+0.3,b.cr.ub.1dose[nonarrowindex],col="grey")
  
  # change line color of effectiveness to black, keep the efficacy estimate as grey
  obs_nonarrowindex <- c(1,2,4,6,8,9,13)
  points(index[obs_nonarrowindex], b.cr.1dose[obs_nonarrowindex], pch=15,
         cex=psize[c(1,2,4,5,6,7,9)],col="black",
         xaxt = "n",xlab="",ylab="kOCV Efficacy/Effectiveness (%)",bty="n")
  segments(index[obs_nonarrowindex]-0.3,b.cr.lb.1dose[obs_nonarrowindex], index[obs_nonarrowindex]+0.3,b.cr.lb.1dose[obs_nonarrowindex],col="black")
  segments(index[obs_nonarrowindex],b.cr.lb.1dose[obs_nonarrowindex],
           index[obs_nonarrowindex],b.cr.ub.1dose[obs_nonarrowindex],col="black")
  segments(index[obs_nonarrowindex]-0.3,b.cr.ub.1dose[obs_nonarrowindex],
           index[obs_nonarrowindex]+0.3,b.cr.ub.1dose[obs_nonarrowindex],col="black")
  
  # add CI for studies with large CI
  arrowindex=c(3,10,14)
  points(index[arrowindex], b.cr.1dose[arrowindex], pch=15,
         cex=psize[c(3,8,10)],col="black",
         xaxt = "n",xlab="",ylab="kOCV Efficacy/Effectiveness (%)",bty="n")
  segments(index[arrowindex]-0.3,b.cr.lb.1dose[arrowindex], index[arrowindex]+0.3,b.cr.lb.1dose[arrowindex],col="black")
  segments(index[arrowindex],b.cr.lb.1dose[arrowindex],
           index[arrowindex],rep(-50,length.out=length(arrowindex)),col="black")
  for (i in 1:length(arrowindex)){
    polygon(y=c(-50,-46,-46,-50),
            x=c(index[arrowindex[i]],index[arrowindex[i]]-0.2 ,
                index[arrowindex[i]]+0.2,index[arrowindex[i]]),
            col="black",border="black")
  }
  
  # add pooled estimates
  for(i in c(1,2,3)){
    if(i == 1){
      polygon(y = c(ci.lb.1dose[i], yi.1dose[i], ci.ub.1dose[i], yi.1dose[i]),
              x=c(index_right[i],index_right[i]-0.3,
                  index_right[i],index_right[i]+0.3),col="steelblue", border = "steelblue")
    }
    if(i == 2){
      polygon(y = c(ci.lb.1dose[i], yi.1dose[i], ci.ub.1dose[i], yi.1dose[i]),
              x=c(index_right[3],index_right[3]-0.3,
                  index_right[3],index_right[3]+0.3),col="steelblue", border = "steelblue")
    }
    if(i == 3){
      polygon(y = c(ci.lb.1dose[i], yi.1dose[i], ci.ub.1dose[i], yi.1dose[i]),
              x=c(index_right[5],index_right[5]-0.3,
                  index_right[5],index_right[5]+0.3),col="steelblue", border = "steelblue")
    }
    
  }
  
  # add axis
  axis(side=4, at=c(-50,0,50,100))
  # axis(side=3,pos=100,xlab="",tck=0,col="white",
  #      at=index_center[1],
  #      labels=c(paste(format(yi.obs,digits=2),"% [",as.numeric(format(ci.ub.obs,digits=2)),
  #                     ",",as.numeric(format(ci.lb.obs,digits=2)),"]",sep="")))
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
  axis(side=1,pos=-70,xlab="",tck=0,col="white", col.axis="dimgrey",
       at=index_center[2],
       labels="93% [69,98]")
  axis(side=1,pos=-70,xlab="",tck=0,col="white", col.axis="dimgrey",
       at=index_center[4],
       labels="67% [30,84]")
  
  # axis(side=3,pos=100,xlab="",tck=0,col="white",
  #      at=index_center[1],labels=c("45% [26,59]"))
  # axis(side=3,pos=100,xlab="", col.axis = "dimgray",tck=0,col="white",
  #      at=index_center[2],labels=c("64% [43,77]"))
  axis(side=1,pos=-50,xlab="",tck=0,lty=3,col="white",font=2,
       at=c(4.5,13.5,22.5,31.5,40.5), # mid point of each panel
       labels=c("0-6m","6-12m","12-18m","18-24m","24-30m"))
  abline(h=c(100,0,-50),lty=3,col="grey")
  #abline(h=c(0,50),lty=3,col=AddAlpha("grey",0.4))
  axis(side=3,pos=110,xlab="",tck=0,lty=3,at=22.5,
       labels=c("One-dose estimates"))
  axis(side=1,pos=-100,xlab="",tck=0,lty=3,col="white",
       at=22.5,
       labels=c("Follow-up duration in month"))
  title(main = "A", adj = 0, line = 2)
  # add estimate for each study
  # axis(side=3, pos=90,tck=0,col="white",las=2, cex.axis=0.6,
  #      at=index,labels= paste0(
  #        round(b.cr.rct.1dose), " [", round(b.cr.ub.rct.1dose), ", ", round(b.cr.lb.rct.1dose), "]"
  #      ))
  
  # add label for each study
  # axis(side=3, pos=90,tck=0,col="white",las=2, cex.axis=0.5,
  #      at=index,labels= fpdata_rct_1dose$Label)
  
  # change the label color of new studies to steelblue
  # axis(side=3, pos=90,tck=0,col="white",las=2, cex.axis=0.5,col.axis="black",
  #      at=index[which(fpdata_rct_1dose$new.old == "old")],labels=fpdata_rct_1dose$Paper[which(fpdata_rct_1dose$new.old == "old")]
  # )
  # axis(side=3, pos=90,tck=0,col="white",las=2, cex.axis=0.5,col.axis="blue",
  #      at=index[which(fpdata_rct_1dose$new.old == "new")],labels=fpdata_rct_1dose$Paper[which(fpdata_rct_1dose$new.old == "new")]
  # )
}



#################################################
# Function to make 4-panel plot of VE over time (2-dose only)
#################################################
make_fp_ve_2dose <- function(fpdata = df_FollowupPeriod){
  
  fpdata <- fpdata %>% filter(!is.na(yearid))
  fpdata$fPaper <- factor(fpdata$Paper, levels = unique(fpdata$Paper))
  fpdata$uid <- 1:nrow(fpdata)
  fpdata$fuid <- factor(fpdata$uid)
  tau_estim_method <- "EB" # seems like a reasonable choice from DOI: 10.1002/sim.2688
  

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
  par(mfrow=c(2,1),mar=c(7,7,7,7))
  par(mar=c(2,6,6,6))
  plot(1:43,rep(NA,length.out=43),pch=15,
       ylim=c(-50,100),col="grey",
       xaxt = "n",ylab="kOCV Efficacy (%)",bty="n",
       xlab="")
  
  points(index,100*(1-exp(log_ve_fup_rct$yi)),pch=15,
         cex=psize,col="black",
         xaxt = "n",xlab="",ylab="kOCV Efficacy (%)",bty="n")
  
  # add CI
  nonarrowindex=which(b.cr.ub.rct>=(-50))
  arrowindex=which(b.cr.ub.rct<(-50))
  segments(index-0.3,b.cr.lb.rct,index+0.3,b.cr.lb.rct,col="black")
  segments(index[nonarrowindex],b.cr.lb.rct[nonarrowindex],
           index[nonarrowindex],b.cr.ub.rct[nonarrowindex],col="black")
  segments(index[nonarrowindex]-0.3,b.cr.ub.rct[nonarrowindex],
           index[nonarrowindex]+0.3,b.cr.ub.rct[nonarrowindex],col="black")
  # add CI for studies with large CI
  segments(index[arrowindex],b.cr.lb.rct[arrowindex],
           index[arrowindex],rep(-50,length.out=length(arrowindex)),col="black")
  for (i in 1:length(arrowindex)){
    polygon(y=c(-50,-46,-46,-50),
            x=c(index[arrowindex[i]],index[arrowindex[i]]-0.2 ,
                index[arrowindex[i]]+0.2,index[arrowindex[i]]),
            col="black",border="black")
  }
  
  # change line color of new studies to blue
  # newstudyindex = which(fpdata_rct$new.old == "new")
  # points(index[newstudyindex],100*(1-exp(log_ve_fup_rct$yi[newstudyindex])),pch=15,
  #        cex=psize[newstudyindex],col="blue",
  #        xaxt = "n",xlab="",ylab="kOCV Efficacy (%)",bty="n")
  # segments(index[newstudyindex]-0.3,b.cr.lb.rct[newstudyindex],index[newstudyindex]+0.3,b.cr.lb.rct[newstudyindex],col="blue")
  # segments(index[newstudyindex],b.cr.lb.rct[newstudyindex],
  #          index[newstudyindex],b.cr.ub.rct[newstudyindex],col="blue")
  # segments(index[newstudyindex]-0.3,b.cr.ub.rct[newstudyindex],
  #          index[newstudyindex]+0.3,b.cr.ub.rct[newstudyindex],col="blue")
  
  # add pooled estimates
  for (i in 1:length(unique(fpdata_rct$yearid))){
    polygon(y = c(ci.lb.rct[i], yi.rct[i], ci.ub.rct[i], yi.rct[i]),
            x=c(index_right[i],index_right[i]-0.3,
                index_right[i],index_right[i]+0.3),col="steelblue", border = "steelblue")
  }
  
  # add axis and text
  axis(side=4, at=c(-50,0,50,100))
  axis(side=1,pos=-70,xlab="",tck=0,col="white",
       at=index_center[1:4],
       labels=c(paste(format(yi.rct,digits=2),"% [",
                      as.numeric(format(ci.ub.rct,digits=2)),
                      ",",
                      as.numeric(format(ci.lb.rct,digits=2)),"]",sep="")))
  axis(side=1,pos=-70,xlab="",tck=0,col.axis="dimgrey",col="white",
       at=index_center[5],labels=c("81% [41,94]"))
  # axis(side=1,pos=-75,xlab="",tck=0,col="white",
  #      at=index_center[4],labels=c("26% [-46,63]"))
  axis(side=1,pos=-50,xlab="",tck=0,lty=3,col="dimgrey",font=2,
       at=index_center,
       labels=c("0-12m","12-24m","24-36m","36-48m","48-60m"))
  axis(side=3,pos=100,xlab="",tck=0,lty=3,at=index_center[3],
       labels=c("Two-dose estimates"))
  abline(h=c(100,0,-50),lty=3,col="grey")
  title(main = "B", adj = 0, line = 2.5)
  #abline(h=c(0,50),lty=3,col=AddAlpha("grey",0.4))
  
  # add estimate for each study
  # axis(side=3, pos=90,tck=0,col="white",las=2, cex.axis=0.6,
  #      at=index,labels= paste0(
  #        round(b.cr.rct), " [", round(b.cr.ub.rct), ", ", round(b.cr.lb.rct), "]"
  #      ))
  
  # # add label for each study
  # axis(side=3, pos=90,tck=0,col="white",las=2, cex.axis=0.5,
  #      at=index,labels=fpdata_rct$Label
  #      )
  
  # change the label color of new studies to blue
  # axis(side=3, pos=90,tck=0,col="white",las=2, cex.axis=0.5,col.axis="black",
  #      at=index[which(fpdata_rct$new.old == "old")],labels=fpdata_rct$Paper[which(fpdata_rct$new.old == "old")]
  # )
  # axis(side=3, pos=90,tck=0,col="white",las=2, cex.axis=0.5,col.axis="blue",
  #      at=index[which(fpdata_rct$new.old == "new")],labels=fpdata_rct$Paper[which(fpdata_rct$new.old == "new")]
  # )
  
  
  ##########################
  ####### set up effectiveness data for 2 dose
  ##########################
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
  par(mar=c(6,6,3,6))
  plot(1:43,rep(NA,length.out=43),pch=15,
       ylim=c(-50,100),col="grey",
       xaxt = "n",ylab="kOCV Effectiveness (%)",bty="n",
       xlab="")
  
  points(index,100*(1-exp(log_ve_fup_obs$yi)),pch=15,
         cex=psize,col="black",
         xaxt = "n",xlab="",ylab="kOCV Efficacy (%)",bty="n")
  
  # add CI
  nonarrowindex=which(b.cr.ub.obs>=(-50))
  
  segments(index-0.3,b.cr.lb.obs,index+0.3,b.cr.lb.obs,col="black")
  segments(index[nonarrowindex],b.cr.lb.obs[nonarrowindex],
           index[nonarrowindex],b.cr.ub.obs[nonarrowindex],col="black")
  segments(index[nonarrowindex]-0.3,b.cr.ub.obs[nonarrowindex],
           index[nonarrowindex]+0.3,b.cr.ub.obs[nonarrowindex],col="black")
  
  # change line color of new studies to blue
  # newstudyindex = which(fpdata_obs$new.old == "new")
  # points(index[newstudyindex],100*(1-exp(log_ve_fup_obs$yi[newstudyindex])),pch=15,
  #        cex=psize[newstudyindex],col="blue",
  #        xaxt = "n",xlab="",ylab="kOCV Efficacy (%)",bty="n")
  # segments(index[newstudyindex]-0.3,b.cr.lb.obs[newstudyindex],index[newstudyindex]+0.3,b.cr.lb.obs[newstudyindex],col="blue")
  # segments(index[newstudyindex],b.cr.lb.obs[newstudyindex],
  #          index[newstudyindex],b.cr.ub.obs[newstudyindex],col="blue")
  # segments(index[newstudyindex]-0.3,b.cr.ub.obs[newstudyindex],
  #          index[newstudyindex]+0.3,b.cr.ub.obs[newstudyindex],col="blue")
  
  
  # add CI for studies with large CI
  arrowindex=which(b.cr.ub.obs<(-50))
  segments(index[arrowindex],b.cr.lb.obs[arrowindex],
           index[arrowindex],rep(-50,length.out=length(arrowindex)),col="black")
  # for (i in 1:length(arrowindex)){
  #   polygon(y=c(-50,-46,-46,-50),
  #           x=c(index[arrowindex[i]],index[arrowindex[i]]-0.2 ,
  #               index[arrowindex[i]]+0.2,index[arrowindex[i]]),
  #           col="grey",border="grey")
  # }
  
  # add pooled estimates
  for(i in 1:3){
    polygon(y = c(ci.lb.obs[i], yi.obs[i], ci.ub.obs[i], yi.obs[i]),
            x=c(index_right[i],index_right[i]-0.3,
                index_right[i],index_right[i]+0.3),col="steelblue", border = "steelblue")
  }
  
  
  # add axis
  axis(side=4, at=c(-50,0,50,100))
  axis(side=1,pos=-75,xlab="",tck=0,col="white",
       at=index_center[1:3],
       labels=c(paste(format(yi.obs,digits=2),"% [",as.numeric(format(ci.ub.obs,digits=2)),
                      ",",as.numeric(format(ci.lb.obs,digits=2)),"]",sep="")))
  axis(side=1,pos=-75,xlab="",tck=0,col.axis="dimgrey",col="white",
       at=index_center[4],labels=c("94% [60,99]"))
  axis(side=1,pos=-50,xlab="",tck=0,lty=3,col="white",font=2,
       at=c(4.5,13.5,22.5,31.5, 40.5),
       labels=c("0-12m","12-24m","24-36m","36-48m","48-60m"))
  axis(side=1,pos=-100,xlab="",tck=0,lty=3,col="white",
       at=22.5,
       labels=c("Follow-up duration in month"))
  abline(h=c(100,0,-50),lty=3,col="grey")
  #abline(h=c(0,50),lty=3,col=AddAlpha("grey",0.4))
  
  # add estimate for each study
  # axis(side=3, pos=90,tck=0,col="white",las=2, cex.axis=0.6,
  #      at=index,labels= paste0(
  #        round(b.cr.obs), " [", round(b.cr.ub.obs), ", ", round(b.cr.lb.obs), "]"
  #      ))
  
  # # add label for each study
  # axis(side=3, pos=90,tck=0,col="white",las=2, cex.axis=0.5,
  #      at=index,labels= fpdata_obs$Label)
  
  # change the label color of new studies to blue
  # axis(side=3, pos=90,tck=0,col="white",las=2, cex.axis=0.5,col.axis="black",
  #      at=index[which(fpdata_obs$new.old == "old")],labels=fpdata_obs$Paper[which(fpdata_obs$new.old == "old")]
  # )
  # axis(side=3, pos=90,tck=0,col="white",las=2, cex.axis=0.5,col.axis="blue",
  #      at=index[which(fpdata_obs$new.old == "new")],labels=fpdata_obs$Paper[which(fpdata_obs$new.old == "new")]
  # )
  
}

# make plot for age stratified VE
plot_age_stratified_ve <- function(data = df.age.plot,
                                   VE.ratio.pooled,
                                   VE.ratio.pooled.w.ci,
                                   add_pooled_estim = FALSE){
  
  # order the studies by dose number
  data <- data %>% arrange(desc(Doses)) %>%
    rename(under5 = `Under 5 yo`, above5 = `Above 5 yo`) %>%
    mutate(StudyDesign = ifelse(StudyDesign == "Randomized", "RCT", "OBS")) %>%
    rename(Design = StudyDesign)

  par(mar=c(2,25,2,0))
  plot(c(data$VE.ratio, NA), 1:9 ,pch=19,
       xlim=c(-0.4,1.2),
       yaxt = "n",ylab="",bty="n",
       xlab="", cex = 0.5, cex.axis=0.4, tck=0,
       mgp = c(3,-0.3,0))


  # add estimates
  axis(side=2,pos=-5.8,xlab="",tck=0,lty=3,col="white",font=1,
       at=c(1:8),
       labels = data$StudyID, las = 1, cex.axis=0.7)
  axis(side=2,pos=-5.2,xlab="",tck=0,lty=3,col="white",font=1,
       at=c(1:8),
       labels = data$Design, las = 1, cex.axis=0.7)
  axis(side=2,pos=-4.7,xlab="",tck=0,lty=3,col="white",font=1,
       at=c(1:8), 
       labels = data$Doses, las = 1, cex.axis=0.7)
  axis(side=2,pos=-3.1,xlab="",tck=0,lty=3,col="white",font=1,
       at=c(1:8),
       labels = data$under5, las = 1, cex.axis=0.7)
  axis(side=2,pos= -1.6,xlab="",tck=0,lty=3,col="white",font=1,
       at=c(1:8), 
       labels = data$above5, las = 1, cex.axis=0.7)
  axis(side=2,pos= -0.1,xlab="",tck=0,lty=3,col="white",font=1,
       at=c(1:8), 
       labels = data$VE.ratio.CI, las = 1, cex.axis=0.7)
  
  # add estimates' titles
  axis(side=2,pos= -6.0,xlab="",tck=0,lty=3,col="white",font=2,
       at=9, # mid point of each panel
       labels = "Study", las = 1, cex.axis=0.7)
  axis(side=2,pos= -5.1,xlab="",tck=0,lty=3,col="white",font=2,
       at=9, # mid point of each panel
       labels = "Design", las = 1, cex.axis=0.7)
  axis(side=2,pos= -4.4,xlab="",tck=0,lty=3,col="white",font=2,
       at=9, # mid point of each panel
       labels = "Doses", las = 1, cex.axis=0.7)
  axis(side=2,pos= -3.3,xlab="",tck=0,lty=3,col="white",font=2,
       at=9, # mid point of each panel
       labels = "VE(<5)", las = 1, cex.axis=0.7)
  axis(side=2,pos= -1.8,xlab="",tck=0,lty=3,col="white",font=2,
       at=9, # mid point of each panel
       labels = "VE(>5)", las = 1, cex.axis=0.7)
  axis(side=2,pos= -0.1,xlab="",tck=0,lty=3,col="white",font=2,
       at=9, # mid point of each panel
       labels = "Relative VE", las = 1, cex.axis=0.7)
  
  lines(x=c(1, 1), y=c(0,8.5), lty = 3)
  
  if(add_pooled_estim == TRUE){
    # add line for pooled estimate 
    lines(x=c(VE.ratio.pooled, VE.ratio.pooled), y=c(0,9.5))
    
    axis(side=1,pos= 1, xlab="",tck=0,lty=3,col="white",font=2,
    at=VE.ratio.pooled, # mid point of each panel
    labels = VE.ratio.pooled.w.ci, las = 1, cex.axis=0.8)
  }
}


#################################
# Function to make 4-panel plot
#################################
make_fp_ve_4panels <- function(fpdata = df_FollowupPeriod){
  
  par(mfrow=c(2,2),mar=c(7,7,7,7))
  
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
         cex=psize,col="black",
         xaxt = "n",xlab="",ylab="kOCV Efficacy (%)",bty="n")
  
  # add CI
  nonarrowindex=which(b.cr.ub.rct>=(-50))
  arrowindex=which(b.cr.ub.rct<(-50))
  segments(index-0.3,b.cr.lb.rct,index+0.3,b.cr.lb.rct,col="black")
  segments(index[nonarrowindex],b.cr.lb.rct[nonarrowindex],
           index[nonarrowindex],b.cr.ub.rct[nonarrowindex],col="black")
  segments(index[nonarrowindex]-0.3,b.cr.ub.rct[nonarrowindex],
           index[nonarrowindex]+0.3,b.cr.ub.rct[nonarrowindex],col="black")
  # add CI for studies with large CI
  segments(index[arrowindex],b.cr.lb.rct[arrowindex],
           index[arrowindex],rep(-50,length.out=length(arrowindex)),col="black")
  for (i in 1:length(arrowindex)){
    polygon(y=c(-50,-46,-46,-50),
            x=c(index[arrowindex[i]],index[arrowindex[i]]-0.2 ,
                index[arrowindex[i]]+0.2,index[arrowindex[i]]),
            col="black",border="black")
  }
  
  # add pooled estimates
  for (i in 1:length(unique(fpdata_rct$yearid))){
    polygon(y = c(ci.lb.rct[i], yi.rct[i], ci.ub.rct[i], yi.rct[i]),
            x=c(index_right[i],index_right[i]-0.3,
                index_right[i],index_right[i]+0.3),col="steelblue", border = "steelblue")
  }
  
  # add axis and text
  axis(side=4, at=c(-50,0,50,100))
  axis(side=1,pos=-70,xlab="",tck=0,col="white",
       at=index_center[1:4],
       labels=c(paste(format(yi.rct,digits=2),"% [",
                      as.numeric(format(ci.ub.rct,digits=2)),
                      ",",
                      as.numeric(format(ci.lb.rct,digits=2)),"]",sep="")))
  axis(side=1,pos=-70,xlab="",tck=0,col.axis="dimgrey",col="white",
       at=index_center[5],labels=c("81% [42,94]")) # display original figures in the paper instead of reconstructed 
  # axis(side=1,pos=-75,xlab="",tck=0,col="white",
  #      at=index_center[4],labels=c("26% [-46,63]"))
  axis(side=1,pos=-50,xlab="",tck=0,lty=3,col="dimgrey",font=2,
       at=index_center,
       labels=c("0-12m","12-24m","24-36m","36-48m","48-60m"))
  axis(side=3,pos=100,xlab="",tck=0,lty=3,at=index_center[3],
       labels=c("Two-dose efficacy estimates"))
  abline(h=c(100,0,-50),lty=3,col="grey")
  title(main = "A", adj = 0, line = 2.5)
  axis(side=1,pos=-90,xlab="",tck=0,lty=3,col="white",
       at=22.5,
       labels=c("Follow-up duration in month"))
  
  
  ##########################
  # 2 dose effectiveness 
  ##########################
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
         cex=psize,col="black",
         xaxt = "n",xlab="",ylab="kOCV Efficacy (%)",bty="n")
  
  # add CI
  nonarrowindex=which(b.cr.ub.obs>=(-50))
  
  segments(index-0.3,b.cr.lb.obs,index+0.3,b.cr.lb.obs,col="black")
  segments(index[nonarrowindex],b.cr.lb.obs[nonarrowindex],
           index[nonarrowindex],b.cr.ub.obs[nonarrowindex],col="black")
  segments(index[nonarrowindex]-0.3,b.cr.ub.obs[nonarrowindex],
           index[nonarrowindex]+0.3,b.cr.ub.obs[nonarrowindex],col="black")
  
  # add CI for studies with large CI
  arrowindex=which(b.cr.ub.obs<(-50))
  segments(index[arrowindex],b.cr.lb.obs[arrowindex],
           index[arrowindex],rep(-50,length.out=length(arrowindex)),col="black")
  for (i in 1:length(arrowindex)){
    polygon(y=c(-50,-46,-46,-50),
            x=c(index[arrowindex[i]],index[arrowindex[i]]-0.2 ,
                index[arrowindex[i]]+0.2,index[arrowindex[i]]),
            col="black",border="black")
  }
  
  # add pooled estimates
  for(i in 1:3){
    polygon(y = c(ci.lb.obs[i], yi.obs[i], ci.ub.obs[i], yi.obs[i]),
            x=c(index_right[i],index_right[i]-0.3,
                index_right[i],index_right[i]+0.3),col="steelblue", border = "steelblue")
  }
  
  # add axis
  axis(side=4, at=c(-50,0,50,100))
  axis(side=1,pos=-70,xlab="",tck=0,col="white",
       at=index_center[1:3],
       labels=c(paste(format(yi.obs,digits=2),"% [",as.numeric(format(ci.ub.obs,digits=2)),
                      ",",as.numeric(format(ci.lb.obs,digits=2)),"]",sep="")))
  axis(side=1,pos=-70,xlab="",tck=0,col.axis="dimgrey",col="white",
       at=index_center[4],labels=c("94% [56,99]"))
  axis(side=1,pos=-50,xlab="",tck=0,lty=3,col="white",font=2,
       at=c(4.5,13.5,22.5,31.5, 40.5),
       labels=c("0-12m","12-24m","24-36m","36-48m","48-60m"))
  axis(side=1,pos=-90,xlab="",tck=0,lty=3,col="white",
       at=22.5,
       labels=c("Follow-up duration in month"))
  axis(side=3,pos=100,xlab="",tck=0,lty=3,at=index_center[3],
       labels=c("Two-dose effectiveness estimates"))
  title(main = "B", adj = 0, line = 2.5)
  abline(h=c(100,0,-50),lty=3,col="grey")
  
  
  ################## 
  # 1-dose estimates (efficacy)
  ##################
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
    arrange(StudyType, .by_group = TRUE)
  
  # remove the duplicate estimate: Qadri 2016 (0-6 months)
  temp <- temp %>%
    filter(Paper != "Qadri et al, 2016")
  
  fpdata_1dose_rct <- temp %>% filter(StudyType == "Efficacy")
  
  # there are only 4 studies in each follow up period 
  b.cr.lb.1dose <- fpdata_1dose_rct$VE.u * 100
  b.cr.ub.1dose <- fpdata_1dose_rct$VE.l * 100
  b.cr.1dose <- fpdata_1dose_rct$Mean.VE* 100
  
  ###########################
  # Updated Horizontal Figure with efficacy for 1 dose
  ########################
  index_count=table(fpdata_1dose_rct$monthid)
  index_center=c(4.5,13.5,22.5,31.5) # mid point of each panel
  index= index_center
  
  # base plot
  par(mar=c(6,4,4,6))
  plot(1:43,rep(NA,length.out=43),pch=15,
       ylim=c(-50,100),col="grey",
       xaxt = "n",ylab="kOCV Efficacy (%)",bty="n",
       xlab="")
  
  # points for efficacy studies
  points(index_center, fpdata_1dose_rct$Mean.VE*100,pch=15,
         cex=1,col="black",
         xaxt = "n",xlab="",ylab="kOCV Efficacy (%)",bty="n")
  
  # add CI 
  nonarrowindex=which(b.cr.ub.1dose>=(-50))
  
  segments(index-0.3,b.cr.lb.1dose,index+0.3,b.cr.lb.1dose,col="black")
  segments(index[nonarrowindex],b.cr.lb.1dose[nonarrowindex],
           index[nonarrowindex],b.cr.ub.1dose[nonarrowindex],col="black")
  segments(index[nonarrowindex]-0.3,b.cr.ub.1dose[nonarrowindex],
           index[nonarrowindex]+0.3,b.cr.ub.1dose[nonarrowindex],col="black")
  
  # add axis
  axis(side=4, at=c(-50,0,50,100))
  # axis(side=3,pos=100,xlab="",tck=0,col="white",
  #      at=index_center[1],
  #      labels=c(paste(format(yi.obs,digits=2),"% [",as.numeric(format(ci.ub.obs,digits=2)),
  #                     ",",as.numeric(format(ci.lb.obs,digits=2)),"]",sep="")))
  axis(side=1,pos=-70,xlab="",tck=0,col="white", col.axis="dimgrey",
       at=index_center[1],
       labels=c(paste(format(b.cr.1dose[[1]],digits=2),"% [",
                      as.numeric(format(b.cr.ub.1dose[[1]],digits=2)),
                      ",",
                      as.numeric(format(b.cr.lb.1dose[[1]],digits=2)),"]",sep="")))
  axis(side=1,pos=-70,xlab="",tck=0,col="white", col.axis="dimgrey",
       at=index_center[2],
       labels=c(paste(format(b.cr.1dose[[2]],digits=2),"% [",
                      as.numeric(format(b.cr.ub.1dose[[2]],digits=2)),
                      ",",
                      as.numeric(format(b.cr.lb.1dose[[2]],digits=2)),"]",sep="")))
  axis(side=1,pos=-70,xlab="",tck=0,col="white", col.axis="dimgrey",
       at=index_center[3],
       labels=c(paste(format(b.cr.1dose[[3]],digits=2),"% [",
                      as.numeric(format(b.cr.ub.1dose[[3]],digits=2)),
                      ",",
                      as.numeric(format(b.cr.lb.1dose[[3]],digits=2)),"]",sep="")))
  axis(side=1,pos=-70,xlab="",tck=0,col="white", col.axis="dimgrey",
       at=index_center[4],
       labels=c(paste(format(b.cr.1dose[[4]],digits=2),"% [",
                      as.numeric(format(b.cr.ub.1dose[[4]],digits=2)),
                      ",",
                      as.numeric(format(b.cr.lb.1dose[[4]],digits=2)),"]",sep="")))
  
  axis(side=1,pos=-50,xlab="",tck=0,lty=3,col="white",font=2,
       at=c(4.5,13.5,22.5,31.5,40.5), # mid point of each panel
       labels=c("0-6m","6-12m","12-18m","18-24m","24-30m"))
  abline(h=c(100,0,-50),lty=3,col="grey")
  #abline(h=c(0,50),lty=3,col=AddAlpha("grey",0.4))
  axis(side=3,pos=100,xlab="",tck=0,lty=3,at=22.5,
       labels=c("One-dose efficacy estimates"))
  axis(side=1,pos=-90,xlab="",tck=0,lty=3,col="white",
       at=22.5,
       labels=c("Follow-up duration in month"))
  title(main = "C", adj = 0, line = 2)
  
  
  ################## 
  # 1-dose estimates (effectiveness)
  ##################
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
    arrange(StudyType, .by_group = TRUE)
  
  # remove the duplicate estimate: Qadri 2016 (0-6 months)
  temp <- temp %>%
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
  
  #pooled estimates
  ci.lb.1dose=c(fpdata_rma_1dose_1$ci.lb, fpdata_rma_1dose_3$ci.lb, fpdata_rma_1dose_5$ci.lb)
  ci.ub.1dose=c(fpdata_rma_1dose_1$ci.ub, fpdata_rma_1dose_3$ci.ub, fpdata_rma_1dose_5$ci.ub)
  yi.1dose=c(fpdata_rma_1dose_1$b, fpdata_rma_1dose_3$b, fpdata_rma_1dose_5$b)
  
  ci.lb.1dose=100*(1-exp(ci.lb.1dose))
  ci.ub.1dose=100*(1-exp(ci.ub.1dose))
  yi.1dose=100*(1-exp(yi.1dose))
  
  ###########################
  # Updated Horizontal Figure with effectiveness for 1 dose
  ########################
  index_count=table(fpdata_1dose$monthid)
  index_center=c(4.5,13.5,22.5, 40.5) # mid point of each panel
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
  
  # base plot
  par(mar=c(6,6,4,4))
  plot(1:43,rep(NA,length.out=43),pch=15,
       ylim=c(-50,100),col="grey",
       xaxt = "n",ylab="kOCV Effectiveness (%)",bty="n",
       xlab="")
  
  # points for effectiveness studies
  temp_mid <- 100*(1-exp(log_ve_fup_1dose$yi))
  #temp_mid[5] <- 90
  
  points(index,temp_mid,pch=15,
         cex=psize,col="black",
         xaxt = "n",xlab="",ylab="kOCV Effectiveness (%)",bty="n")
  
  # add CI 
  nonarrowindex=which(b.cr.ub.1dose>=(-50))
  # b.cr.ub.1dose[5] <- 54
  # b.cr.lb.1dose[5] <- 98
  
  segments(index-0.3,b.cr.lb.1dose,index+0.3,b.cr.lb.1dose,col="black")
  segments(index[nonarrowindex],b.cr.lb.1dose[nonarrowindex],
           index[nonarrowindex],b.cr.ub.1dose[nonarrowindex],col="black")
  segments(index[nonarrowindex]-0.3,b.cr.ub.1dose[nonarrowindex],
           index[nonarrowindex]+0.3,b.cr.ub.1dose[nonarrowindex],col="black")
  
  # add CI for studies with large CI
  arrowindex=c(3,8,10)
  points(index[arrowindex], b.cr.1dose[arrowindex], pch=15,
         cex=psize[c(3,8,10)],col="black",
         xaxt = "n",xlab="",ylab="kOCV Effectiveness (%)",bty="n")
  segments(index[arrowindex]-0.3,b.cr.lb.1dose[arrowindex], index[arrowindex]+0.3,b.cr.lb.1dose[arrowindex],col="black")
  segments(index[arrowindex],b.cr.lb.1dose[arrowindex],
           index[arrowindex],rep(-50,length.out=length(arrowindex)),col="black")
  for (i in 1:length(arrowindex)){
    polygon(y=c(-50,-46,-46,-50),
            x=c(index[arrowindex[i]],index[arrowindex[i]]-0.2 ,
                index[arrowindex[i]]+0.2,index[arrowindex[i]]),
            col="black",border="black")
  }
  
  # add pooled estimates
  for(i in c(1,2,3)){
    if(i == 1){
      polygon(y = c(ci.lb.1dose[i], yi.1dose[i], ci.ub.1dose[i], yi.1dose[i]),
              x=c(index_right[i],index_right[i]-0.3,
                  index_right[i],index_right[i]+0.3),col="steelblue", border = "steelblue")
    }
    if(i == 2){
      polygon(y = c(ci.lb.1dose[i], yi.1dose[i], ci.ub.1dose[i], yi.1dose[i]),
              x=c(index_right[3],index_right[3]-0.3,
                  index_right[3],index_right[3]+0.3),col="steelblue", border = "steelblue")
    }
    if(i == 3){
      polygon(y = c(ci.lb.1dose[i], yi.1dose[i], ci.ub.1dose[i], yi.1dose[i]),
              x=c(index_right[4],index_right[4]-0.3,
                  index_right[4],index_right[4]+0.3),col="steelblue", border = "steelblue")
    }
    
  }
  
  # add axis
  axis(side=4, at=c(-50,0,50,100))
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
       at=index_center[4],
       labels=c(paste(format(yi.1dose[[3]],digits=2),"% [",
                      as.numeric(format(ci.ub.1dose[[3]],digits=2)),
                      ",",
                      as.numeric(format(ci.lb.1dose[[3]],digits=2)),"]",sep="")))
  axis(side=1,pos=-70,xlab="",tck=0,col="white", col.axis="dimgrey",
       at=index_center[2],
       labels="92% [66,98]")
  
  axis(side=1,pos=-50,xlab="",tck=0,lty=3,col="white",font=2,
       at=c(4.5,13.5,22.5,31.5,40.5), # mid point of each panel
       labels=c("0-6m","6-12m","12-18m","18-24m","24-30m"))
  abline(h=c(100,0,-50),lty=3,col="grey")
  axis(side=3,pos=100,xlab="",tck=0,lty=3,at=22.5,
       labels=c("One-dose effectiveness estimates"))
  axis(side=1,pos=-90,xlab="",tck=0,lty=3,col="white",
       at=22.5,
       labels=c("Follow-up duration in month"))
  title(main = "D", adj = 0, line = 2)
  
  
  # # return the plot object 
  # grid.echo()
  # Fig <- grid.grab()
  
  # return(Fig)
  
}




#################################
# Function to make 4-panel plot (colored by vaccine product)
#################################
make_fp_ve_4panels_colored <- function(fpdata = df_FollowupPeriod){
  
  colors <- brewer.pal(n = 8, name = "Dark2")[2:5]
  
  # first, assign vaccine product used for each study
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
  
  

  par(mfrow=c(2,2),mar=c(7,7,7,7))
  
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
         cex=psize,
         col = fpdata_rct$color,
         xaxt = "n",xlab="",ylab="kOCV Efficacy (%)",bty="n")
  
  # add CI
  nonarrowindex=which(b.cr.ub.rct>=(-50))
  arrowindex=which(b.cr.ub.rct<(-50))
  segments(index-0.3,b.cr.lb.rct,index+0.3,b.cr.lb.rct,col= fpdata_rct$color)
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
            col= fpdata_rct$color[arrowindex[i]],border=fpdata_rct$color[arrowindex[i]])
  }
  
  # add pooled estimates
  for (i in 1:length(unique(fpdata_rct$yearid))){
    polygon(y = c(ci.lb.rct[i], yi.rct[i], ci.ub.rct[i], yi.rct[i]),
            x=c(index_right[i],index_right[i]-0.3,
                index_right[i],index_right[i]+0.3),col="grey61", border = "grey61")
  }
  
  # add axis and text
  axis(side=4, at=c(-50,0,50,100))
  axis(side=1,pos=-70,xlab="",tck=0,col="white",
       at=index_center[1:4],
       labels=c(paste(format(yi.rct,digits=2),"% [",
                      as.numeric(format(ci.ub.rct,digits=2)),
                      ",",
                      as.numeric(format(ci.lb.rct,digits=2)),"]",sep="")))
  axis(side=1,pos=-70,xlab="",tck=0,col.axis="dimgrey",col="white",
       at=index_center[5],labels=c("81% [42,94]")) # display original figures in the paper instead of reconstructed 
  # axis(side=1,pos=-75,xlab="",tck=0,col="white",
  #      at=index_center[4],labels=c("26% [-46,63]"))
  axis(side=1,pos=-50,xlab="",tck=0,lty=3,col="dimgrey",font=2,
       at=index_center,
       labels=c("0-12m","12-24m","24-36m","36-48m","48-60m"))
  axis(side=3,pos=100,xlab="",tck=0,lty=3,at=index_center[3],
       labels=c("Two-dose efficacy estimates"))
  abline(h=c(100,0,-50),lty=3,col="grey")
  title(main = "A", adj = 0, line = 2.5)
  axis(side=1,pos=-90,xlab="",tck=0,lty=3,col="white",
       at=22.5,
       labels=c("Follow-up duration in month"))
  
  
  ##########################
  # 2 dose effectiveness 
  ##########################
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
         cex=psize,
         col= fpdata_obs$color,
         xaxt = "n",xlab="",ylab="kOCV Efficacy (%)",bty="n")
  
  # add CI
  nonarrowindex=which(b.cr.ub.obs>=(-50))
  
  segments(index-0.3,b.cr.lb.obs,index+0.3,b.cr.lb.obs,col= fpdata_obs$color)
  segments(index[nonarrowindex],b.cr.lb.obs[nonarrowindex],
           index[nonarrowindex],b.cr.ub.obs[nonarrowindex],col=fpdata_obs$color[nonarrowindex])
  segments(index[nonarrowindex]-0.3,b.cr.ub.obs[nonarrowindex],
           index[nonarrowindex]+0.3,b.cr.ub.obs[nonarrowindex],col=fpdata_obs$color[nonarrowindex])
  
  # add CI for studies with large CI
  arrowindex=which(b.cr.ub.obs<(-50))
  segments(index[arrowindex],b.cr.lb.obs[arrowindex],
           index[arrowindex],rep(-50,length.out=length(arrowindex)),col= fpdata_obs$color[arrowindex])
  for (i in 1:length(arrowindex)){
    polygon(y=c(-50,-46,-46,-50),
            x=c(index[arrowindex[i]],index[arrowindex[i]]-0.2 ,
                index[arrowindex[i]]+0.2,index[arrowindex[i]]),
            col= fpdata_obs$color[arrowindex[i]],border=fpdata_obs$color[arrowindex[i]])
  }
  
  # add pooled estimates
  for(i in 1:3){
    polygon(y = c(ci.lb.obs[i], yi.obs[i], ci.ub.obs[i], yi.obs[i]),
            x=c(index_right[i],index_right[i]-0.3,
                index_right[i],index_right[i]+0.3),col="grey61", border = "grey61")
  }
  
  # add axis
  axis(side=4, at=c(-50,0,50,100))
  axis(side=1,pos=-70,xlab="",tck=0,col="white",
       at=index_center[1:3],
       labels=c(paste(format(yi.obs,digits=2),"% [",as.numeric(format(ci.ub.obs,digits=2)),
                      ",",as.numeric(format(ci.lb.obs,digits=2)),"]",sep="")))
  axis(side=1,pos=-70,xlab="",tck=0,col.axis="dimgrey",col="white",
       at=index_center[4],labels=c("94% [56,99]"))
  axis(side=1,pos=-50,xlab="",tck=0,lty=3,col="white",font=2,
       at=c(4.5,13.5,22.5,31.5, 40.5),
       labels=c("0-12m","12-24m","24-36m","36-48m","48-60m"))
  axis(side=1,pos=-90,xlab="",tck=0,lty=3,col="white",
       at=22.5,
       labels=c("Follow-up duration in month"))
  axis(side=3,pos=100,xlab="",tck=0,lty=3,at=index_center[3],
       labels=c("Two-dose effectiveness estimates"))
  title(main = "B", adj = 0, line = 2.5)
  abline(h=c(100,0,-50),lty=3,col="grey")
  
  
  ################## 
  # 1-dose estimates (efficacy)
  ##################
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
    arrange(StudyType, .by_group = TRUE)
  
  # remove the duplicate estimate: Qadri 2016 (0-6 months)
  temp <- temp %>%
    filter(Paper != "Qadri et al, 2016")
  
  fpdata_1dose_rct <- temp %>% filter(StudyType == "Efficacy")
  
  # there are only 4 studies in each follow up period 
  b.cr.lb.1dose <- fpdata_1dose_rct$VE.u * 100
  b.cr.ub.1dose <- fpdata_1dose_rct$VE.l * 100
  b.cr.1dose <- fpdata_1dose_rct$Mean.VE* 100
  
  ###########################
  # Updated Horizontal Figure with efficacy for 1 dose
  ########################
  index_count=table(fpdata_1dose_rct$monthid)
  index_center=c(4.5,13.5,22.5,31.5) # mid point of each panel
  index= index_center
  
  # base plot
  par(mar=c(6,4,4,6))
  plot(1:43,rep(NA,length.out=43),pch=15,
       ylim=c(-50,100),col="grey",
       xaxt = "n",ylab="kOCV Efficacy (%)",bty="n",
       xlab="")
  
  # points for efficacy studies
  points(index_center, fpdata_1dose_rct$Mean.VE*100,pch=15,
         cex=1,
         col= fpdata_1dose_rct$color,
         xaxt = "n",xlab="",ylab="kOCV Efficacy (%)",bty="n")
  
  # add CI 
  nonarrowindex=which(b.cr.ub.1dose>=(-50))
  
  segments(index-0.3,b.cr.lb.1dose,index+0.3,b.cr.lb.1dose,col=fpdata_1dose_rct$color)
  segments(index[nonarrowindex],b.cr.lb.1dose[nonarrowindex],
           index[nonarrowindex],b.cr.ub.1dose[nonarrowindex],col=fpdata_1dose_rct$color[nonarrowindex])
  segments(index[nonarrowindex]-0.3,b.cr.ub.1dose[nonarrowindex],
           index[nonarrowindex]+0.3,b.cr.ub.1dose[nonarrowindex],col=fpdata_1dose_rct$color[nonarrowindex])
  
  # add axis
  axis(side=4, at=c(-50,0,50,100))
  # axis(side=3,pos=100,xlab="",tck=0,col="white",
  #      at=index_center[1],
  #      labels=c(paste(format(yi.obs,digits=2),"% [",as.numeric(format(ci.ub.obs,digits=2)),
  #                     ",",as.numeric(format(ci.lb.obs,digits=2)),"]",sep="")))
  axis(side=1,pos=-70,xlab="",tck=0,col="white", col.axis="dimgrey",
       at=index_center[1],
       labels=c(paste(format(b.cr.1dose[[1]],digits=2),"% [",
                      as.numeric(format(b.cr.ub.1dose[[1]],digits=2)),
                      ",",
                      as.numeric(format(b.cr.lb.1dose[[1]],digits=2)),"]",sep="")))
  axis(side=1,pos=-70,xlab="",tck=0,col="white", col.axis="dimgrey",
       at=index_center[2],
       labels=c(paste(format(b.cr.1dose[[2]],digits=2),"% [",
                      as.numeric(format(b.cr.ub.1dose[[2]],digits=2)),
                      ",",
                      as.numeric(format(b.cr.lb.1dose[[2]],digits=2)),"]",sep="")))
  axis(side=1,pos=-70,xlab="",tck=0,col="white", col.axis="dimgrey",
       at=index_center[3],
       labels=c(paste(format(b.cr.1dose[[3]],digits=2),"% [",
                      as.numeric(format(b.cr.ub.1dose[[3]],digits=2)),
                      ",",
                      as.numeric(format(b.cr.lb.1dose[[3]],digits=2)),"]",sep="")))
  axis(side=1,pos=-70,xlab="",tck=0,col="white", col.axis="dimgrey",
       at=index_center[4],
       labels=c(paste(format(b.cr.1dose[[4]],digits=2),"% [",
                      as.numeric(format(b.cr.ub.1dose[[4]],digits=2)),
                      ",",
                      as.numeric(format(b.cr.lb.1dose[[4]],digits=2)),"]",sep="")))
  
  axis(side=1,pos=-50,xlab="",tck=0,lty=3,col="white",font=2,
       at=c(4.5,13.5,22.5,31.5,40.5), # mid point of each panel
       labels=c("0-6m","6-12m","12-18m","18-24m","24-30m"))
  abline(h=c(100,0,-50),lty=3,col="grey")
  #abline(h=c(0,50),lty=3,col=AddAlpha("grey",0.4))
  axis(side=3,pos=100,xlab="",tck=0,lty=3,at=22.5,
       labels=c("One-dose efficacy estimates"))
  axis(side=1,pos=-90,xlab="",tck=0,lty=3,col="white",
       at=22.5,
       labels=c("Follow-up duration in month"))
  title(main = "C", adj = 0, line = 2)
  
  
  ################## 
  # 1-dose estimates (effectiveness)
  ##################
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
    arrange(StudyType, .by_group = TRUE)
  
  # remove the duplicate estimate: Qadri 2016 (0-6 months)
  temp <- temp %>%
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
  
  #pooled estimates
  ci.lb.1dose=c(fpdata_rma_1dose_1$ci.lb, fpdata_rma_1dose_3$ci.lb, fpdata_rma_1dose_5$ci.lb)
  ci.ub.1dose=c(fpdata_rma_1dose_1$ci.ub, fpdata_rma_1dose_3$ci.ub, fpdata_rma_1dose_5$ci.ub)
  yi.1dose=c(fpdata_rma_1dose_1$b, fpdata_rma_1dose_3$b, fpdata_rma_1dose_5$b)
  
  ci.lb.1dose=100*(1-exp(ci.lb.1dose))
  ci.ub.1dose=100*(1-exp(ci.ub.1dose))
  yi.1dose=100*(1-exp(yi.1dose))
  
  ###########################
  # Updated Horizontal Figure with effectiveness for 1 dose
  ########################
  index_count=table(fpdata_1dose$monthid)
  index_center=c(4.5,13.5,22.5, 40.5) # mid point of each panel
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
  
  # base plot
  par(mar=c(6,6,4,4))
  plot(1:43,rep(NA,length.out=43),pch=15,
       ylim=c(-50,100),col="grey",
       xaxt = "n",ylab="kOCV Effectiveness (%)",bty="n",
       xlab="")
  
  # points for effectiveness studies
  temp_mid <- 100*(1-exp(log_ve_fup_1dose$yi))
  # temp_mid[5] <- 90
  
  points(index,temp_mid,pch=15,
         cex=psize,
         col= fpdata_1dose$color,
         xaxt = "n",xlab="",ylab="kOCV Effectiveness (%)",bty="n")
  
  # add CI 
  nonarrowindex=which(b.cr.ub.1dose>=(-50))
  # b.cr.ub.1dose[5] <- 54
  # b.cr.lb.1dose[5] <- 98
  
  segments(index-0.3,b.cr.lb.1dose,index+0.3,b.cr.lb.1dose,col= fpdata_1dose$color)
  segments(index[nonarrowindex],b.cr.lb.1dose[nonarrowindex],
           index[nonarrowindex],b.cr.ub.1dose[nonarrowindex],col=fpdata_1dose$color[nonarrowindex])
  segments(index[nonarrowindex]-0.3,b.cr.ub.1dose[nonarrowindex],
           index[nonarrowindex]+0.3,b.cr.ub.1dose[nonarrowindex],col=fpdata_1dose$color[nonarrowindex])
  
  # add CI for studies with large CI
  arrowindex=c(3,8,10)
  segments(index[arrowindex]-0.3,b.cr.lb.1dose[arrowindex], index[arrowindex]+0.3,b.cr.lb.1dose[arrowindex],col=fpdata_1dose$color[arrowindex])
  segments(index[arrowindex],b.cr.lb.1dose[arrowindex],
           index[arrowindex],rep(-50,length.out=length(arrowindex)),col= fpdata_1dose$color[arrowindex])
  for (i in 1:length(arrowindex)){
    polygon(y=c(-50,-46,-46,-50),
            x=c(index[arrowindex[i]],index[arrowindex[i]]-0.2 ,
                index[arrowindex[i]]+0.2,index[arrowindex[i]]),
            col= fpdata_1dose$color[arrowindex[i]],border=fpdata_1dose$color[arrowindex[i]])
  }
  
  # add pooled estimates
  for(i in c(1,2,3)){
    if(i == 1){
      polygon(y = c(ci.lb.1dose[i], yi.1dose[i], ci.ub.1dose[i], yi.1dose[i]),
              x=c(index_right[i],index_right[i]-0.3,
                  index_right[i],index_right[i]+0.3),col="grey61", border = "grey61")
    }
    if(i == 2){
      polygon(y = c(ci.lb.1dose[i], yi.1dose[i], ci.ub.1dose[i], yi.1dose[i]),
              x=c(index_right[3],index_right[3]-0.3,
                  index_right[3],index_right[3]+0.3),col="grey61", border = "grey61")
    }
    if(i == 3){
      polygon(y = c(ci.lb.1dose[i], yi.1dose[i], ci.ub.1dose[i], yi.1dose[i]),
              x=c(index_right[4],index_right[4]-0.3,
                  index_right[4],index_right[4]+0.3),col="grey61", border = "grey61")
    }
    
  }
  
  # add axis
  axis(side=4, at=c(-50,0,50,100))
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
       at=index_center[4],
       labels=c(paste(format(yi.1dose[[3]],digits=2),"% [",
                      as.numeric(format(ci.ub.1dose[[3]],digits=2)),
                      ",",
                      as.numeric(format(ci.lb.1dose[[3]],digits=2)),"]",sep="")))
  axis(side=1,pos=-70,xlab="",tck=0,col="white", col.axis="dimgrey",
       at=index_center[2],
       labels="92% [66,98]")
  
  axis(side=1,pos=-50,xlab="",tck=0,lty=3,col="white",font=2,
       at=c(4.5,13.5,22.5,31.5,40.5), # mid point of each panel
       labels=c("0-6m","6-12m","12-18m","18-24m","24-30m"))
  abline(h=c(100,0,-50),lty=3,col="grey")
  axis(side=3,pos=100,xlab="",tck=0,lty=3,at=22.5,
       labels=c("One-dose effectiveness estimates"))
  axis(side=1,pos=-90,xlab="",tck=0,lty=3,col="white",
       at=22.5,
       labels=c("Follow-up duration in month"))
  title(main = "D", adj = 0, line = 2)
  
  
}


#################################################
# Function to make 4-panel plot of VE over time (both 1 dose and 2 dose VE) 
#################################################
make_fp_ve <- function(fpdata = df_FollowupPeriod){
  
  fpdata <- fpdata %>% filter(!is.na(yearid))
  fpdata$fPaper <- factor(fpdata$Paper, levels = unique(fpdata$Paper))
  fpdata$uid <- 1:nrow(fpdata)
  fpdata$fuid <- factor(fpdata$uid)
  tau_estim_method <- "EB" # seems like a reasonable choice from DOI: 10.1002/sim.2688
  
  par(mfrow=c(2,2),mar=c(7,7,7,7))
  
  
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
  par(mar=c(6,6,6,6))
  plot(1:43,rep(NA,length.out=43),pch=15,
       ylim=c(-50,100),col="grey",
       xaxt = "n",ylab="kOCV Efficacy (%)",bty="n",
       xlab="")
  
  points(index,100*(1-exp(log_ve_fup_rct$yi)),pch=15,
         cex=psize,col="grey",
         xaxt = "n",xlab="",ylab="kOCV Efficacy (%)",bty="n")
  
  # add CI
  nonarrowindex=which(b.cr.ub.rct>=(-50))
  arrowindex=which(b.cr.ub.rct<(-50))
  segments(index-0.3,b.cr.lb.rct,index+0.3,b.cr.lb.rct,col="grey")
  segments(index[nonarrowindex],b.cr.lb.rct[nonarrowindex],
           index[nonarrowindex],b.cr.ub.rct[nonarrowindex],col="grey")
  segments(index[nonarrowindex]-0.3,b.cr.ub.rct[nonarrowindex],
           index[nonarrowindex]+0.3,b.cr.ub.rct[nonarrowindex],col="grey")
  # add CI for studies with large CI
  segments(index[arrowindex],b.cr.lb.rct[arrowindex],
           index[arrowindex],rep(-50,length.out=length(arrowindex)),col="grey")
  for (i in 1:length(arrowindex)){
    polygon(y=c(-50,-46,-46,-50),
            x=c(index[arrowindex[i]],index[arrowindex[i]]-0.2 ,
                index[arrowindex[i]]+0.2,index[arrowindex[i]]),
            col="grey",border="grey")
  }
  
  # change line color of new studies to blue
  newstudyindex = which(fpdata_rct$new.old == "new")
  points(index[newstudyindex],100*(1-exp(log_ve_fup_rct$yi[newstudyindex])),pch=15,
         cex=psize[newstudyindex],col="blue",
         xaxt = "n",xlab="",ylab="kOCV Efficacy (%)",bty="n")
  segments(index[newstudyindex]-0.3,b.cr.lb.rct[newstudyindex],index[newstudyindex]+0.3,b.cr.lb.rct[newstudyindex],col="blue")
  segments(index[newstudyindex],b.cr.lb.rct[newstudyindex],
           index[newstudyindex],b.cr.ub.rct[newstudyindex],col="blue")
  segments(index[newstudyindex]-0.3,b.cr.ub.rct[newstudyindex],
           index[newstudyindex]+0.3,b.cr.ub.rct[newstudyindex],col="blue")
  
  # add pooled estimates
  for (i in 1:length(unique(fpdata_rct$yearid))){
    polygon(y = c(ci.lb.rct[i], yi.rct[i], ci.ub.rct[i], yi.rct[i]),
            x=c(index_right[i],index_right[i]-0.3,
                index_right[i],index_right[i]+0.3),col="black")
  }
  
  # add axis and text
  axis(side=4, at=c(-50,0,50,100))
  axis(side=1,pos=-75,xlab="",tck=0,col="white",
       at=index_center[1:4],
       labels=c(paste(format(yi.rct,digits=2),"% [",
                      as.numeric(format(ci.ub.rct,digits=2)),
                      ",",
                      as.numeric(format(ci.lb.rct,digits=2)),"]",sep="")))
  axis(side=1,pos=-75,xlab="",tck=0,col.axis="dimgrey",col="white",
       at=index_center[5],labels=c("81% [41,94]"))
  # axis(side=1,pos=-75,xlab="",tck=0,col="white",
  #      at=index_center[4],labels=c("26% [-46,63]"))
  axis(side=1,pos=-50,xlab="",tck=0,lty=3,col="dimgrey",font=2,
       at=index_center,
       labels=c("0-12m","12-24m","24-36m","36-48m","48-60m"))
  axis(side=3,pos=150,xlab="",tck=0,lty=3,at=index_center[3],
       labels=c("Two-dose estimates"))
  abline(h=c(100,0,-50),lty=3,col="grey")
  #abline(h=c(0,50),lty=3,col=AddAlpha("grey",0.4))
  
  # add estimate for each study
  # axis(side=3, pos=90,tck=0,col="white",las=2, cex.axis=0.6,
  #      at=index,labels= paste0(
  #        round(b.cr.rct), " [", round(b.cr.ub.rct), ", ", round(b.cr.lb.rct), "]"
  #      ))
  
  # # add label for each study
  # axis(side=3, pos=90,tck=0,col="white",las=2, cex.axis=0.5,
  #      at=index,labels=fpdata_rct$Label
  #      )
  
  # change the label color of new studies to blue
  axis(side=3, pos=90,tck=0,col="white",las=2, cex.axis=0.5,col.axis="black",
       at=index[which(fpdata_rct$new.old == "old")],labels=fpdata_rct$Paper[which(fpdata_rct$new.old == "old")]
  )
  axis(side=3, pos=90,tck=0,col="white",las=2, cex.axis=0.5,col.axis="blue",
       at=index[which(fpdata_rct$new.old == "new")],labels=fpdata_rct$Paper[which(fpdata_rct$new.old == "new")]
  )
  
  ##########################
  ####### set up efficacy data for 1 dose
  ##########################
  fpdata_rct_1dose <- fpdata %>% mutate(se = ifelse(is.na(VE.u),
                                                    (log(1-Mean.VE)-log(1-VE.l))/qnorm(0.95),
                                                    (log(1-VE.l) - log(1-VE.u))/(2*qnorm(.975))),
                                        yi = log(1-Mean.VE)) %>%
    arrange((yearid)) %>%
    filter(StudyType=="Efficacy")  %>%
    filter(Dose==1) %>%
    filter(yearid!=0)
  
  # overall efficacy
  log_ve_fup_rct_1dose <- rma(yi=log(1-fpdata_rct_1dose$Mean.VE),
                              sei=fpdata_rct_1dose$se,
                              method=tau_estim_method)
  
  # estimate yearly pe
  # fpdata_rct_1dose_1 <- subset(fpdata_rct_1dose,yearid==1)
  # fpdata_rct_rma_1dose_1 <- rma(yi=log(1-fpdata_rct_1dose_1$Mean.VE),
  #                               sei=fpdata_rct_1dose_1$se,
  #                               method=tau_estim_method)
  # fpdata_rct_1dose_2 <- subset(fpdata_rct_1dose,yearid==2)
  # fpdata_rct_rma_1dose_2 <- rma(yi=log(1-fpdata_rct_1dose_2$Mean.VE),
  #                               sei=fpdata_rct_1dose_2$se,
  #                               method=tau_estim_method)
  
  
  # CI
  b.cr.lb.rct.1dose=100*(1-exp(log_ve_fup_rct_1dose$yi-qnorm(0.975)*sqrt(log_ve_fup_rct_1dose$vi)))
  b.cr.ub.rct.1dose=100*(1-exp(log_ve_fup_rct_1dose$yi+qnorm(0.975)*sqrt(log_ve_fup_rct_1dose$vi)))
  b.cr.rct.1dose = 100*(1-exp(log_ve_fup_rct_1dose$yi))
  
  # pooled estimates
  # ci.lb.rct.1dose=c(fpdata_rct_rma_1dose_1$ci.lb)
  # ci.ub.rct.1dose=c(fpdata_rct_rma_1dose_1$ci.ub)
  # yi.rct.1dose=c(fpdata_rct_rma_1dose_1$b)
  # 
  # ci.lb.rct.1dose=100*(1-exp(ci.lb.rct.1dose))
  # ci.ub.rct.1dose=100*(1-exp(ci.ub.rct.1dose))
  # yi.rct.1dose=100*(1-exp(yi.rct.1dose))
  
  ###########################
  # Updated Horizontal Figure with just Efficacy Study for 1 dose
  ########################
  index_count=table(fpdata_rct_1dose$yearid)
  # index_center=c(4.5,13.5,22.5,31.5,40.5) # mid point of each panel
  # index_left=index_center[length(index_count)]-((index_count+1)/2-0.5)  # left most index of each panel
  # index_right=index_center[length(index_count)]+((index_count+1)/2-0.5)  # right most index of each panel
  index_center=c(4.5,13.5) # mid point of each panel
  index_left=index_center-((index_count+1)/2-0.5)  # left most index of each panel
  index_right=index_center+((index_count+1)/2-0.5)  # right most index of each panel
  index=NULL
  for (i in 1:length(index_count)){
    index=c(index,seq(index_left[i],length.out=index_count[i]))
  }
  
  # point size
  wi <- 1/sqrt(log_ve_fup_rct_1dose$vi)
  psize <- wi/sum(wi, na.rm = TRUE)
  psize <- (psize - min(psize, na.rm = TRUE))/(max(psize,na.rm = TRUE) - min(psize, na.rm = TRUE))
  psize <- (psize * 1) + 0.5
  
  # base plot
  par(mar=c(6,6,6,6))
  plot(1:43,rep(NA,length.out=43),pch=15,
       ylim=c(-50,100),col="grey",
       xaxt = "n",ylab="kOCV Efficacy (%)",bty="n",
       xlab="")
  
  points(index,100*(1-exp(log_ve_fup_rct_1dose$yi)),pch=15,
         cex=psize,col="grey",
         xaxt = "n",xlab="",ylab="kOCV Efficacy (%)",bty="n")
  
  # add CI
  nonarrowindex=which(b.cr.ub.rct.1dose>=(-50))
  
  segments(index-0.3,b.cr.lb.rct.1dose,index+0.3,b.cr.lb.rct.1dose,col="grey")
  segments(index[nonarrowindex],b.cr.lb.rct.1dose[nonarrowindex],
           index[nonarrowindex],b.cr.ub.rct.1dose[nonarrowindex],col="grey")
  segments(index[nonarrowindex]-0.3,b.cr.ub.rct.1dose[nonarrowindex],
           index[nonarrowindex]+0.3,b.cr.ub.rct.1dose[nonarrowindex],col="grey")
  
  # change line color of new studies to blue
  newstudyindex = which(fpdata_rct_1dose$new.old == "new")
  points(index[newstudyindex],100*(1-exp(log_ve_fup_rct_1dose$yi[newstudyindex])),pch=15,
         cex=psize[newstudyindex],col="blue",
         xaxt = "n",xlab="",ylab="kOCV Efficacy (%)",bty="n")
  segments(index[newstudyindex]-0.3,b.cr.lb.rct.1dose[newstudyindex],index[newstudyindex]+0.3,b.cr.lb.rct.1dose[newstudyindex],col="blue")
  segments(index[newstudyindex],b.cr.lb.rct.1dose[newstudyindex],
           index[newstudyindex],b.cr.ub.rct.1dose[newstudyindex],col="blue")
  segments(index[newstudyindex]-0.3,b.cr.ub.rct.1dose[newstudyindex],
           index[newstudyindex]+0.3,b.cr.ub.rct.1dose[newstudyindex],col="blue")
  
  # # add CI for studies with large CI
  # arrowindex=which(b.cr.ub.rct.1dose<(-50))
  # segments(index[arrowindex],b.cr.lb.rct.1dose[arrowindex],
  #          index[arrowindex],rep(-50,length.out=length(arrowindex)),col="grey")
  # for (i in 1:length(arrowindex)){
  #   polygon(y=c(-50,-46,-46,-50),
  #           x=c(index[arrowindex[i]],index[arrowindex[i]]-0.2 ,
  #               index[arrowindex[i]]+0.2,index[arrowindex[i]]),
  #           col="grey",border="grey")
  # }
  #
  # add pooled estimates
  # for(i in 1:2){
  #   polygon(y = c(ci.lb.rct.1dose[i], yi.rct.1dose[i], ci.ub.rct.1dose[i], yi.rct.1dose[i]),
  #           x=c(index_right[i],index_right[i]-0.3,
  #               index_right[i],index_right[i]+0.3),col="black")
  # }
  
  # add axis
  axis(side=4, at=c(-50,0,50,100))
  # axis(side=3,pos=100,xlab="",tck=0,col="white",
  #      at=index_center[1],
  #      labels=c(paste(format(yi.obs,digits=2),"% [",as.numeric(format(ci.ub.obs,digits=2)),
  #                     ",",as.numeric(format(ci.lb.obs,digits=2)),"]",sep="")))
  axis(side=1,pos=-75,xlab="",tck=0,col="white", col.axis="dimgrey",
       at=index_center[1],
       labels=c(paste(format(b.cr.rct.1dose[[1]],digits=2),"% [",
                      as.numeric(format(b.cr.ub.rct.1dose[[1]],digits=2)),
                      ",",
                      as.numeric(format(b.cr.lb.rct.1dose[[1]],digits=2)),"]",sep="")))
  axis(side=1,pos=-75,xlab="",tck=0,col="white", col.axis="dimgrey",
       at=index_center[2],
       labels=c(paste(format(b.cr.rct.1dose[[2]],digits=2),"% [",
                      as.numeric(format(b.cr.ub.rct.1dose[[2]],digits=2)),
                      ",",
                      as.numeric(format(b.cr.lb.rct.1dose[[2]],digits=2)),"]",sep="")))
  
  # axis(side=3,pos=100,xlab="",tck=0,col="white",
  #      at=index_center[1],labels=c("45% [26,59]"))
  # axis(side=3,pos=100,xlab="", col.axis = "dimgray",tck=0,col="white",
  #      at=index_center[2],labels=c("64% [43,77]"))
  axis(side=1,pos=-50,xlab="",tck=0,lty=3,col="white",font=2,
       at=c(4.5,13.5,22.5,31.5,40.5), # mid point of each panel
       labels=c("0-12m","12-24m","24-36m","36-48m","48-60m"))
  abline(h=c(100,0,-50),lty=3,col="grey")
  #abline(h=c(0,50),lty=3,col=AddAlpha("grey",0.4))
  axis(side=3,pos=120,xlab="",tck=0,lty=3,at=22.5,
       labels=c("One-dose estimates"))
  
  # add estimate for each study
  # axis(side=3, pos=90,tck=0,col="white",las=2, cex.axis=0.6,
  #      at=index,labels= paste0(
  #        round(b.cr.rct.1dose), " [", round(b.cr.ub.rct.1dose), ", ", round(b.cr.lb.rct.1dose), "]"
  #      ))
  
  # add label for each study
  # axis(side=3, pos=90,tck=0,col="white",las=2, cex.axis=0.5,
  #      at=index,labels= fpdata_rct_1dose$Label)
  
  # change the label color of new studies to steelblue
  axis(side=3, pos=90,tck=0,col="white",las=2, cex.axis=0.5,col.axis="black",
       at=index[which(fpdata_rct_1dose$new.old == "old")],labels=fpdata_rct_1dose$Paper[which(fpdata_rct_1dose$new.old == "old")]
  )
  axis(side=3, pos=90,tck=0,col="white",las=2, cex.axis=0.5,col.axis="blue",
       at=index[which(fpdata_rct_1dose$new.old == "new")],labels=fpdata_rct_1dose$Paper[which(fpdata_rct_1dose$new.old == "new")]
  )
  
  ##########################
  ####### set up effectiveness data for 2 dose
  ##########################
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
  par(mar=c(6,6,3,6))
  plot(1:43,rep(NA,length.out=43),pch=15,
       ylim=c(-50,100),col="grey",
       xaxt = "n",ylab="kOCV Effectiveness (%)",bty="n",
       xlab="")
  
  points(index,100*(1-exp(log_ve_fup_obs$yi)),pch=15,
         cex=psize,col="grey",
         xaxt = "n",xlab="",ylab="kOCV Efficacy (%)",bty="n")
  
  # add CI
  nonarrowindex=which(b.cr.ub.obs>=(-50))
  
  segments(index-0.3,b.cr.lb.obs,index+0.3,b.cr.lb.obs,col="grey")
  segments(index[nonarrowindex],b.cr.lb.obs[nonarrowindex],
           index[nonarrowindex],b.cr.ub.obs[nonarrowindex],col="grey")
  segments(index[nonarrowindex]-0.3,b.cr.ub.obs[nonarrowindex],
           index[nonarrowindex]+0.3,b.cr.ub.obs[nonarrowindex],col="grey")
  
  # change line color of new studies to blue
  newstudyindex = which(fpdata_obs$new.old == "new")
  points(index[newstudyindex],100*(1-exp(log_ve_fup_obs$yi[newstudyindex])),pch=15,
         cex=psize[newstudyindex],col="blue",
         xaxt = "n",xlab="",ylab="kOCV Efficacy (%)",bty="n")
  segments(index[newstudyindex]-0.3,b.cr.lb.obs[newstudyindex],index[newstudyindex]+0.3,b.cr.lb.obs[newstudyindex],col="blue")
  segments(index[newstudyindex],b.cr.lb.obs[newstudyindex],
           index[newstudyindex],b.cr.ub.obs[newstudyindex],col="blue")
  segments(index[newstudyindex]-0.3,b.cr.ub.obs[newstudyindex],
           index[newstudyindex]+0.3,b.cr.ub.obs[newstudyindex],col="blue")
  
  
  # add CI for studies with large CI
  arrowindex=which(b.cr.ub.obs<(-50))
  segments(index[arrowindex],b.cr.lb.obs[arrowindex],
           index[arrowindex],rep(-50,length.out=length(arrowindex)),col="grey")
  # for (i in 1:length(arrowindex)){
  #   polygon(y=c(-50,-46,-46,-50),
  #           x=c(index[arrowindex[i]],index[arrowindex[i]]-0.2 ,
  #               index[arrowindex[i]]+0.2,index[arrowindex[i]]),
  #           col="grey",border="grey")
  # }
  
  # add pooled estimates
  for(i in 1:3){
    polygon(y = c(ci.lb.obs[i], yi.obs[i], ci.ub.obs[i], yi.obs[i]),
            x=c(index_right[i],index_right[i]-0.3,
                index_right[i],index_right[i]+0.3),col="black")
  }
  
  
  # add axis
  axis(side=4, at=c(-50,0,50,100))
  axis(side=1,pos=-75,xlab="",tck=0,col="white",
       at=index_center[1:3],
       labels=c(paste(format(yi.obs,digits=2),"% [",as.numeric(format(ci.ub.obs,digits=2)),
                      ",",as.numeric(format(ci.lb.obs,digits=2)),"]",sep="")))
  axis(side=1,pos=-75,xlab="",tck=0,col.axis="dimgrey",col="white",
       at=index_center[4],labels=c("94% [60,99]"))
  axis(side=1,pos=-50,xlab="",tck=0,lty=3,col="white",font=2,
       at=c(4.5,13.5,22.5,31.5, 40.5),
       labels=c("0-12m","12-24m","24-36m","36-48m","48-60m"))
  axis(side=1,pos=-110,xlab="",tck=0,lty=3,col="white",
       at=22.5,
       labels=c("Follow-up duration in month"))
  abline(h=c(100,0,-50),lty=3,col="grey")
  #abline(h=c(0,50),lty=3,col=AddAlpha("grey",0.4))
  
  # add estimate for each study
  # axis(side=3, pos=90,tck=0,col="white",las=2, cex.axis=0.6,
  #      at=index,labels= paste0(
  #        round(b.cr.obs), " [", round(b.cr.ub.obs), ", ", round(b.cr.lb.obs), "]"
  #      ))
  
  # # add label for each study
  # axis(side=3, pos=90,tck=0,col="white",las=2, cex.axis=0.5,
  #      at=index,labels= fpdata_obs$Label)
  
  # change the label color of new studies to blue
  axis(side=3, pos=90,tck=0,col="white",las=2, cex.axis=0.5,col.axis="black",
       at=index[which(fpdata_obs$new.old == "old")],labels=fpdata_obs$Paper[which(fpdata_obs$new.old == "old")]
  )
  axis(side=3, pos=90,tck=0,col="white",las=2, cex.axis=0.5,col.axis="blue",
       at=index[which(fpdata_obs$new.old == "new")],labels=fpdata_obs$Paper[which(fpdata_obs$new.old == "new")]
  )
  
  
  #####################################
  ## set up effectiveness data for 1 ##
  #####################################
  
  fpdata_obs_1dose <- fpdata %>% mutate(se = ifelse(is.na(VE.u),
                                                    (log(1-Mean.VE)-log(1-VE.l))/qnorm(0.95),
                                                    (log(1-VE.l) - log(1-VE.u))/(2*qnorm(.975))),
                                        yi = log(1-Mean.VE)) %>%
    arrange((yearid)) %>%
    filter(StudyType=="Effectiveness")  %>%
    filter(Dose==1) %>%
    filter(yearid!=0)
  
  azman_id <- which(fpdata_obs_1dose$Paper=="Azman et al, 2016")
  fpdata_obs_1dose[azman_id,'se'] <- (log(1-fpdata_obs_1dose[azman_id,'Mean.VE'])-
                                        log(1-fpdata_obs_1dose[azman_id,'VE.l']))/(2*qnorm(.975))*2 #.9966
  
  # overall efficacy
  log_ve_fup_obs_1dose <- rma(yi=log(1-fpdata_obs_1dose$Mean.VE),
                              sei=fpdata_obs_1dose$se,
                              method=tau_estim_method)
  
  # estimate yearly pe
  fpdata_obs_1dose_1 <- subset(fpdata_obs_1dose,yearid==1)
  fpdata_obs_1dose_rma_1 <- rma(yi=yi,
                                sei=se,
                                method=tau_estim_method,
                                data=fpdata_obs_1dose_1)
  fpdata_obs_1dose_2 <- subset(fpdata_obs_1dose,yearid==2)
  fpdata_obs_1dose_rma_2 <- rma(yi=yi,
                                sei=se,
                                method=tau_estim_method,
                                data=fpdata_obs_1dose_2)
  
  # ## test difference between year-1 effiacy and effectivness
  # combined_design_onedose <- rbind(fpdata_obs_1dose_1 ,fpdata_rct_1)
  
  # rma(yi,
  #     sei=se,
  #     method=tau_estim_method,
  #     dat=combined_design_onedose,
  #     mods= ~ StudyType)
  
  ## one-dose combined estimate since no sig difference in short term study designs.
  # one_dose_combined <- rma(yi,
  #     sei=se,
  #     method=tau_estim_method,
  #     dat=combined_design_onedose
  #     )
  # round(100*c(1 - exp(one_dose_combined$b),1 - exp(one_dose_combined$ci.ub),1 - exp(one_dose_combined$ci.lb)),1)
  
  
  # onetwo_dose_effectivness <- rbind(fpdata_obs_1dose_1,fpdata_obs_1)
  # rma(yi,
  #     sei=se,
  #     method=tau_estim_method,
  #     dat=onetwo_dose_effectivness <- rbind(fpdata_obs_1dose_1,fpdata_obs_1),
  #     mods= ~ Dose)
  
  
  # add pooled estimates
  ci.lb.obs.1dose=c(fpdata_obs_1dose_rma_1$ci.lb, fpdata_obs_1dose_rma_2$ci.lb)
  ci.ub.obs.1dose=c(fpdata_obs_1dose_rma_1$ci.ub, fpdata_obs_1dose_rma_2$ci.ub)
  yi.obs.1dose=c(fpdata_obs_1dose_rma_1$b, fpdata_obs_1dose_rma_2$b)
  
  
  ci.lb.obs.1dose=100*(1-exp(ci.lb.obs.1dose))
  ci.ub.obs.1dose=100*(1-exp(ci.ub.obs.1dose))
  yi.obs.1dose=100*(1-exp(yi.obs.1dose))
  
  # CI
  b.cr.lb.obs.1dose=100*(1-exp(log_ve_fup_obs_1dose$yi-qnorm(0.975)*sqrt(log_ve_fup_obs_1dose$vi)))
  b.cr.ub.obs.1dose=100*(1-exp(log_ve_fup_obs_1dose$yi+qnorm(0.975)*sqrt(log_ve_fup_obs_1dose$vi)))
  b.cr.obs.1dose = 100* (1-exp(log_ve_fup_obs_1dose$yi))
  
  ###########################
  # Updated Horizontal Figure with just Effectiveness Study for 1 dose
  ########################
  index_count=table(fpdata_obs_1dose$yearid)
  #index_center=c(4.5,13.5,22.5,31.5,40.5) # mid point of each panel
  index_center=c(4.5,13.5,22.5)
  index_left=index_center[1:length(index_count)]-((index_count+1)/2-0.5)  # left most index of each panel
  index_right=index_center[1:length(index_count)]+((index_count+1)/2-0.5)  # right most index of each panel
  index=NULL
  for (i in 1:length(index_count)){
    index=c(index,seq(index_left[i],length.out=index_count[i]))
  }
  
  # point size
  wi <- 1/sqrt(log_ve_fup_obs_1dose$vi)
  psize <- wi/sum(wi, na.rm = TRUE)
  psize <- (psize - min(psize, na.rm = TRUE))/(max(psize,
                                                   na.rm = TRUE) - min(psize, na.rm = TRUE))
  psize <- (psize * 1) + 0.5
  
  # base plot
  par(mar=c(6,6,3,6))
  plot(1:43,rep(NA,length.out=43),pch=15,
       ylim=c(-50,100),col="grey",
       xaxt = "n",ylab="kOCV Effectiveness (%)",bty="n",
       xlab="")
  
  points(index,100*(1-exp(log_ve_fup_obs_1dose$yi)),pch=15,
         cex=psize,col="grey",
         xaxt = "n",xlab="",ylab="kOCV Effectiveness (%)",bty="n")
  
  # add CI
  nonarrowindex=which(b.cr.ub.obs.1dose>=(-50))
  
  segments(index-0.3,b.cr.lb.obs.1dose,index+0.3,b.cr.lb.obs.1dose,col="grey")
  segments(index[nonarrowindex],b.cr.lb.obs.1dose[nonarrowindex],
           index[nonarrowindex],b.cr.ub.obs.1dose[nonarrowindex],col="grey")
  segments(index[nonarrowindex]-0.3,b.cr.ub.obs.1dose[nonarrowindex],
           index[nonarrowindex]+0.3,b.cr.ub.obs.1dose[nonarrowindex],col="grey")
  
  # change line color of new studies to blue
  newstudyindex = which(fpdata_obs_1dose$new.old == "new")
  points(index[newstudyindex],100*(1-exp(log_ve_fup_obs_1dose$yi[newstudyindex])),pch=15,
         cex=psize[newstudyindex],col="blue",
         xaxt = "n",xlab="",ylab="kOCV Efficacy (%)",bty="n")
  segments(index[newstudyindex]-0.3,b.cr.lb.obs.1dose[newstudyindex],index[newstudyindex]+0.3,b.cr.lb.obs.1dose[newstudyindex],col="blue")
  segments(index[newstudyindex],b.cr.lb.obs.1dose[newstudyindex],
           index[newstudyindex],b.cr.ub.obs.1dose[newstudyindex],col="blue")
  segments(index[newstudyindex]-0.3,b.cr.ub.obs.1dose[newstudyindex],
           index[newstudyindex]+0.3,b.cr.ub.obs.1dose[newstudyindex],col="blue")
  
  
  # add CI for studies with large CI
  arrowindex=which(b.cr.ub.obs.1dose<(-50))
  segments(index[arrowindex],b.cr.lb.obs.1dose[arrowindex],
           index[arrowindex],rep(-50,length.out=length(arrowindex)),col="grey")
  for (i in 1:length(arrowindex)){
    polygon(y=c(-50,-46,-46,-50),
            x=c(index[arrowindex[i]],index[arrowindex[i]]-0.2 ,
                index[arrowindex[i]]+0.2,index[arrowindex[i]]),
            col="grey",border="grey")
  }
  # add pooled estimates
  i=1
  polygon(y = c(ci.lb.obs.1dose[i], yi.obs.1dose[i], ci.ub.obs.1dose[i], yi.obs.1dose[i]),
          x=c(index_right[i],index_right[i]-0.3,
              index_right[i],index_right[i]+0.3),col="black")
  i=2
  polygon(y = c(ci.lb.obs.1dose[i], yi.obs.1dose[i], ci.ub.obs.1dose[i], yi.obs.1dose[i]),
          x=c(index_right[i],index_right[i]-0.3,
              index_right[i],index_right[i]+0.3),col="black")
  
  
  # add axis
  axis(side=4, at=c(-50,0,50,100))
  axis(side=1,pos=-75,xlab="",tck=0,col="white",
       at=index_center[1:2],
       labels=c(paste(format(yi.obs.1dose,digits=2),"% [",as.numeric(format(ci.ub.obs.1dose,digits=2)),
                      ",",as.numeric(format(ci.lb.obs.1dose,digits=2)),"]",sep="")))
  # axis(side=1,pos=-75,xlab="",tck=0,col.axis="dimgrey",col="white",
  #      at=13.5,
  #      labels=c("40% [-31,73]"))
  axis(side=1,pos=-75,xlab="",tck=0,col.axis="dimgrey",col="white",
       at=22.5,labels=c("33% [-318,89]"))
  axis(side=1,pos=-50,xlab="",tck=0,lty=3,col="white",font=2,
       at=c(4.5,13.5,22.5,31.5,40.5),
       labels=c("0-12m","12-24m","24-36m","36-48m","48-60m"))
  axis(side=1,pos=-110,xlab="",tck=0,lty=3,col="white",
       at=22.5,
       labels=c("Follow-up duration in month"))
  abline(h=c(100,0,-50),lty=3,col="grey")
  #abline(h=c(0,50),lty=3,col=AddAlpha("grey",0.4))
  
  # add estimate for each study
  # axis(side=3, pos=90,tck=0,col="white",las=2, cex.axis=0.6,
  #      at=index,labels= paste0(
  #        round(b.cr.obs.1dose), " [", round(b.cr.ub.obs.1dose), ", ", round(b.cr.lb.obs.1dose), "]"
  #      ))
  # axis(side=3, pos=90,tck=0,col="white",las=2, cex.axis=0.5,
  #      at=index,labels= fpdata_obs_1dose$Label)
  
  # change the label color of new studies to blue
  axis(side=3, pos=90,tck=0,col="white",las=2, cex.axis=0.5,col.axis="black",
       at=index[which(fpdata_obs_1dose$new.old == "old")],labels=fpdata_obs_1dose$Paper[which(fpdata_obs_1dose$new.old == "old")]
  )
  axis(side=3, pos=90,tck=0,col="white",las=2, cex.axis=0.5,col.axis="blue",
       at=index[which(fpdata_obs_1dose$new.old == "new")],labels=fpdata_obs_1dose$Paper[which(fpdata_obs_1dose$new.old == "new")]
  )
  
  # # mark immunologically naive setting
  # axis(side=3, pos=90,tck=0,col="white",las=2, cex.axis=0.5,col.axis="red",
  #      at=index[4],labels= fpdata_obs_1dose$Label[4])
  
}


#### function to format the predictions of the rma output
get_rma_pred <- function(obs.pred, rct.pred){
  
  rma_pred <- 
    obs.pred %>%
    as.data.frame() %>% 
    rename(obs.median = pred, obs.lb = ci.lb, obs.ub = ci.ub) %>%
    mutate(year = 1:6) %>% dplyr::select(year, obs.median, obs.lb, obs.ub) %>%
    cbind(
      rct.pred %>% as.data.frame() %>% 
        rename(rct.median = pred, rct.lb = ci.lb, rct.ub = ci.ub) %>%
        dplyr::select(rct.median, rct.lb, rct.ub)
    ) %>%
    mutate(obs.median = 100*obs.median, obs.lb = 100*obs.lb, obs.ub = 100*obs.ub,
           rct.median = 100*rct.median, rct.lb = 100*rct.lb, rct.ub = 100*rct.ub)
  
  return(rma_pred)
}



# Make base plot of waning plot
make_base_waning_plot <- function(df.plt,
                                  est.type){
  
  y_label <- paste0(est.type, " (%)")
  
  plt.base <-
    df.plt %>%
    ggplot() +
    geom_segment(aes(x = TL, xend = TR, y = Mean.VE*100, yend = Mean.VE*100), 
                 data = df.plt, color = "grey") +
    xlab("Months after vaccination") + ylab(y_label) +
    scale_x_continuous(breaks = c(12, 24, 36, 48, 60, 72), 
                       labels = as.character(c(12, 24, 36, 48, 60, 72))) +
    scale_y_continuous(breaks = c(-100, -50, 0,25, 50, 75, 100), 
                       labels = as.character(c(-100, -50, 0, 25, 50, 75, 100))) +
    #ylim(c(-40, 100)) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"))
  
  return(plt.base)
}


















