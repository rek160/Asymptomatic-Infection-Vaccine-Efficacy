# Written by Matt Hitchings
## Edited by Rebecca Kahn, 24th February, 2018
# Function to take results of an epidemic and analyse the data to get vaccine efficacy estimates

  require(survival)
  require(coxme)
  require(plyr)
  require(icenReg)
  
  args=(commandArgs(TRUE))
  
  nsim<-500
  ave_community_size<-20000
  community_size_range<-40
  num_communities<-1
  rate_within<-0.00257 #0.00171 #0.00214 #0.00257	
  rate_between<-0
  ## betaS is force of infection in symptomatics and betaAS in asymptomatics -- simplifying assumption is that they are the same
  betaS<-0.005
  betaAS<-0.005
  ## pS is proportion of the infected in the control group that becomes symptomatic and vpS is the proportion in the vaccine group
  pS<-.2
  vpS<-.2
  ## introductions from main population
  num_introductions<-5
  direct_VE<-0.6
  incperiod_shape<-3300   #mean 6
  incperiod_rate<-550
  infperiod_shape<-1.13
  infperiod_rate<-0.188
  trial_startday<-100
  trial_length<-150
  num_sample<-3
  enrollment_period<-1
  num_ind_enrolled_per_day<-num_communities/enrollment_period ## number of clusters enrolled per day since using stratified indivdual randomization
  sample_percent<-.1
  cluster_coverage<-.075
  ave_inc_period<-ceiling(incperiod_shape/incperiod_rate)
  
  for (i in 1:length(args)) {
    eval (parse (text = args[[i]] ))
  }
  
  list <- structure(NA,class="result")
  "[<-.result" <- function(x,...,value) {
    args <- as.list(match.call())
    args <- args[-c(1:2,length(args))]
    length(value) <- length(args)
    for(i in seq(along=args)) {
      a <- args[[i]]
      if(!missing(a)) eval.parent(substitute(a <- v,list(a=a,v=value[[i]])))
    }
    x
  }
  
  # This code restricts to those who were susceptible or exposed when enrolled in the trial
  # The only people we don't recruit into the trial are those who were infectious at time of enrollment
  # This analysis assumes we restrict to those who were susceptible or exposed on enrollment.
  j<-j
  # directory where networks are stored
  directory<-'~/RCode_RK/1comm1.25/'
  filename<-paste0(directory,"comm_constant_Results_Analysis",cluster_coverage,"_",trial_length,"_",j,"_",betaAS,"_",vpS,"_",rate_within,"_",sample_percent,".csv")
  results_analysis<-read.csv(filename,header=TRUE)
  results_analysis<-results_analysis[,2:8]

  R0 <- (1 - (infperiod_rate/(infperiod_rate+betaS))^infperiod_shape) *
    (((ave_community_size-1)*(1-rate_within)*rate_within + 
        (num_communities-1)*ave_community_size*(1-rate_between)*rate_between + 
        ((ave_community_size-1)*rate_within + (num_communities-1)*ave_community_size*rate_between)^2)/
       ((ave_community_size-1)*rate_within + (num_communities-1)*ave_community_size*rate_between)-1)
  
  if (nrow(results_analysis)>1){
    
    trial_nodes <- nrow(!is.na(results_analysis[results_analysis$TrialStatus]))

    ## Calculate interval for sampling in the interval censored model later
    roundnum<-trial_length/num_sample
    
    ## get syptomatic, asymptomatic and total events in vaccinated and control groups
    numevents_vacc <- nrow(results_analysis[(results_analysis$eventstatus==1) & (results_analysis$TrialStatus==1),])
    numevents_vacc_S <- nrow(results_analysis[(results_analysis$eventstatus==1) & (results_analysis$TrialStatus==1) & 
                                                (results_analysis$Symptomatic==1),])
    numevents_vacc_AS <- nrow(results_analysis[(results_analysis$eventstatus==1) & (results_analysis$TrialStatus==1) & 
                                                 (results_analysis$Symptomatic==0),])
    num_vacc <- nrow(results_analysis[(results_analysis$TrialStatus==1),])
    
    numevents_cont <- nrow(results_analysis[(results_analysis$eventstatus==1) & (results_analysis$TrialStatus==0),])
    numevents_cont_S <- nrow(results_analysis[(results_analysis$eventstatus==1) & (results_analysis$TrialStatus==0) & 
                                                (results_analysis$Symptomatic==1),])
    numevents_cont_AS <- nrow(results_analysis[(results_analysis$eventstatus==1) & (results_analysis$TrialStatus==0) & 
                                                 (results_analysis$Symptomatic==0),])
    num_cont <- nrow(results_analysis[(results_analysis$TrialStatus==0),])
    
    
    ## VES: Vaccine efficacy looking at just symptomatic infections (note VES is not a good name and 
    ## should really be VES_sympt or something else but leaving for now)
    ## Create a datset where event status = 1 only if symptomatic infection (change asymptomatic event status from 1 to 0)
    ## As asymptomatic people are now counted as uninfected, their Day infected should be length of trial when censored
    results_analysis_VES <- results_analysis
    results_analysis_VES$eventstatus[(results_analysis_VES$eventstatus==1)&(results_analysis_VES$Symptomatic==0)]<-0
    results_analysis_VES$DayInfected[(results_analysis_VES$eventstatus==0)]<-trial_length
    
    ## Estimate vaccine efficacy against susceptibility to all infection with cox proportional hazards model
    ## If multiple communities, stratify by community to alleviate bias from heterogeneity in hazards
    if ((numevents_vacc>0)&&(numevents_cont>0)) {
      survmodel<-try(coxph(Surv(DayInfected,eventstatus)~TrialStatus,results_analysis),silent=T)
      usesurvmod <- !inherits(survmodel, 'try-error')
      if (!is.na(survmodel$coefficients[1])){
        vaccEffEst_total <- 1-exp(survmodel$coefficients + c(0, 1.96, -1.96)*sqrt(survmodel$var))
        VEC_var <-survmodel$var
        zval <- survmodel$coefficient/sqrt(survmodel$var)
        VEC_pval <- pnorm(zval, lower.tail = vaccEffEst_total[1]>0)*2
      } else {
        vaccEffEst_total<-c(NA,NA,NA)
        VEC_var <-NA
        VEC_pval <- NA
      }
    } else if ((numevents_vacc>0)&&(numevents_cont==0)) {
      # If there are no events in the control arm but events in the vaccine arm, VE estimate is -1
      vaccEffEst_total<-c(-1,-1,-1)
      VEC_var <-NA
      VEC_pval <-1
    } else if ((numevents_vacc==0)&&(numevents_cont>0)) {
      # If there are no events in the vaccine arm and events in the control arm, VE is 1, but
      # for the p-value we add one event to both arms and do a Cox regression on that data
      # I have chosen to give both events the median time among control events.
      results_analysis1<-results_analysis
      newevent_v_rownum <- min(which((results_analysis1$eventstatus==0)&(results_analysis1$TrialStatus==1)))
      newevent_c_rownum <- min(which((results_analysis1$eventstatus==0)&(results_analysis1$TrialStatus==0)))
      
      eventtime <- median(results_analysis1$DayInfected[results_analysis1$eventstatus==1])
      
      results_analysis1$DayInfected[newevent_v_rownum] <- eventtime
      results_analysis1$eventstatus[newevent_v_rownum] <- 1
      
      results_analysis1$DayInfected[newevent_c_rownum] <- eventtime
      results_analysis1$eventstatus[newevent_c_rownum] <- 1
      
      survmodel<-try(coxph(Surv(DayInfected,eventstatus)~TrialStatus,results_analysis1),silent=T)
      usesurvmod <- !inherits(survmodel, 'try-error')
      
      if (!is.na(survmodel$coefficients[1])){
        vaccEffEst_total<-c(1,1,1)
        VE_fake <- 1-exp(survmodel$coefficients)
        zval <- survmodel$coefficient/sqrt(survmodel$var)
        VEC_pval <- pnorm(zval, lower.tail = VE_fake>0)*2
        VEC_var <-survmodel$var
      } else {
        vaccEffEst_total<-c(1,1,1)
        VEC_pval <- NA
        VEC_var <-NA
      }
    } else {
      # If no events are observed in either arm, the trial has failed and no result can be obtained
      vaccEffEst_total<-c(NA,NA,NA)
      VEC_var <-NA
      VEC_pval <-NA
    }
    
    ## Now calculate vaccine efficacy when only looking at symptomatic infections
    if ((numevents_vacc_S>0)&&(numevents_cont_S>0)) {
      survmodel<-try(coxph(Surv(DayInfected,eventstatus)~TrialStatus,results_analysis_VES),silent=T)
      usesurvmod <- !inherits(survmodel, 'try-error')
      if (!is.na(survmodel$coefficients[1])){
        vaccEffEst_total_S <- 1-exp(survmodel$coefficients + c(0, 1.96, -1.96)*sqrt(survmodel$var))
        zval <- survmodel$coefficient/sqrt(survmodel$var)
        VES_pval <- pnorm(zval, lower.tail = vaccEffEst_total_S[1]>0)*2
        VES_var <-survmodel$var
      } else {
        vaccEffEst_total_S<-c(NA,NA,NA)
        VES_pval<-NA
        VES_var <-NA
      }
    } else if ((numevents_vacc_S>0)&&(numevents_cont_S==0)) {
      # If there are no events in the control arm but events in the vaccine arm, VE estimate is -1
      vaccEffEst_total_S<-c(-1,-1,-1)
      VES_pval<-1
      VES_var <-NA
    } else if ((numevents_vacc_S==0)&&(numevents_cont_S>0)) {
      # If there are no events in the vaccine arm and events in the control arm, VE is 1
      # If there are no events in the vaccine arm and events in the control arm, VE is 1, but
      # for the p-value we add one event to both arms and do a Cox regression on that data
      # What time to assign? I have chosen to give both events the median time
      # among control events.
      results_analysis_VES1<-results_analysis_VES
      newevent_v_rownum <- min(which((results_analysis_VES1$eventstatus==0)&(results_analysis_VES1$TrialStatus==1)))
      newevent_c_rownum <- min(which((results_analysis_VES1$eventstatus==0)&(results_analysis_VES1$TrialStatus==0)))
      
      eventtime <- median(results_analysis$DayInfected[results_analysis_VES1$eventstatus==1])
      
      results_analysis_VES1$DayInfected[newevent_v_rownum] <- eventtime
      results_analysis_VES1$eventstatus[newevent_v_rownum] <- 1
      
      results_analysis_VES1$DayInfected[newevent_c_rownum] <- eventtime
      results_analysis_VES1$eventstatus[newevent_c_rownum] <- 1
      
      survmodel<-try(coxph(Surv(DayInfected,eventstatus)~TrialStatus,results_analysis_VES1),silent=T)
      usesurvmod <- !inherits(survmodel, 'try-error')
      
      if (!is.na(survmodel$coefficients[1])){
        vaccEffEst_total_S<-c(1,1,1)
        VE_fake <- 1-exp(survmodel$coefficients)
        zval <- survmodel$coefficient/sqrt(survmodel$var)
        VES_pval <- pnorm(zval, lower.tail = VE_fake>0)*2
        VES_var <-survmodel$var
      } else {
        vaccEffEst_total_S<-c(1,1,1)
        VES_pval <- NA
        VES_var <-NA
      }
    } else {
      # If no events are observed in either arm, the trial has failed and no result can be obtained
      vaccEffEst_total_S<-c(NA,NA,NA)
      VES_pval <- NA
      VES_var <-NA
    }
    
    ## Caluclate sample sizes used for estimating the above vaccine efficacies
    ## Looking at everyone
    sample_size <- nrow(results_analysis)
    sample_sizeS <- nrow(results_analysis)
    
    ## Sample size for estimating VE against progression to symptoms -- only look within events
    sample_sizeP <- nrow(results_analysis[(results_analysis$eventstatus==1),])
    
    ## Calculate RR estimate 
    ## Counting all infections as events
    VEC_estimate<-1-(numevents_vacc/num_vacc)/(numevents_cont/num_cont)
    
    ## Only counting symptomatic infections as events 
    VES_estimate<-1-(numevents_vacc_S/num_vacc)/(numevents_cont_S/num_cont)
    
    ## VE against progression to symptoms -- only looking within infected
    VEP_estimate<-1-(numevents_vacc_S/numevents_vacc)/(numevents_cont_S/numevents_cont)
    
    ## Need to correct RR estimates because leaky vaccine biases these estimates down
    AR_vacc<-numevents_vacc/num_vacc
    AR_cont<-numevents_cont/num_cont
    VEC_corrected<-1-log(1-AR_vacc)/log(1-AR_cont)
    
    ## Interval censored cox model with time to event for symptomatic and at end for asymptomatic
    ## Note day enrolled is now representing beginning of interval during which they could be infected
    results_analysis_int <- results_analysis
    results_analysis_int$DayInfected <- results_analysis_int$DayInfected + results_analysis_int$DayEnrolled
    
    ## Serologic testing only done once at end of trial so interval for infection is from Day Enrolled to end of trial
    results_analysis_int$DayInfected[(results_analysis_int$Symptomatic==0)]<-trial_length + trial_startday
    
    ## Bginning of interval for symptomatically infected is the day they are infected
    results_analysis_int$DayEnrolled[(results_analysis_int$eventstatus==1) & (results_analysis_int$Symptomatic==1)]<- results_analysis_int$DayInfected[(results_analysis_int$eventstatus==1) & (results_analysis_int$Symptomatic==1)]
    
    ## Non events interval is from the end of the trial through infiinity 
    results_analysis_int$DayEnrolled[(results_analysis_int$eventstatus==0)]<- trial_length + trial_startday
    results_analysis_int$DayInfected[(results_analysis_int$eventstatus==0)]<-Inf

    ## Run interval censored model
    if (numevents_vacc>0 && numevents_cont>0){
      survmodel<- try(ic_sp(cbind(DayEnrolled, DayInfected) ~ TrialStatus, model = 'ph',
                            bs_samples = 50, data = results_analysis_int))
      
      if (!is.na(survmodel$coefficients[1])){
        # If no error was thrown, use the results of the model
        vaccEffEst_total_int <- 1-exp(survmodel$coefficients)
        vaccEffEst_total_int_low <- 1-exp(survmodel$coefficients + 1.96*sqrt(survmodel$var))
        vaccEffEst_total_int_high <- 1-exp(survmodel$coefficients - 1.96*sqrt(survmodel$var))
        VEC_int_var <- survmodel$var
        zval <- survmodel$coefficients[1]/sqrt(survmodel$var[1])
        VEC_pval_int <- pnorm(zval, lower.tail = vaccEffEst_total_int[1]>0)*2
      } else {
        vaccEffEst_total_int <-c(NA,NA,NA)
        vaccEffEst_total_int_low <- c(NA,NA,NA)
        vaccEffEst_total_int_high <- c(NA,NA,NA)
        VEC_int_var <- NA
        VEC_pval_int <- NA
      }
    } else if ((numevents_vacc>0)&&(numevents_cont==0)) {
      # If there are no events in the control arm but events in the vaccine arm, VE estimate is -1
      vaccEffEst_total_int<-c(-1,-1,-1)
      vaccEffEst_total_int_low <- c(-1,-1,-1)
      vaccEffEst_total_int_high <- c(-1,-1,-1)
      VEC_int_var <-NA
      VEC_pval_int <-1
    } else if ((numevents_vacc==0)&&(numevents_cont>0)) {
      # If there are no events in the vaccine arm and events in the control arm, VE is 1, but
      # for the p-value we add one event to both arms and do a Cox regression on that data
      # What time to assign? I have chosen to give both events the median time
      # among control events.
      results_analysis_int2<-results_analysis_int
      newevent_v_rownum <- min(which((results_analysis_int2$eventstatus==0)&(results_analysis_int2$TrialStatus==1)))
      newevent_c_rownum <- min(which((results_analysis_int2$eventstatus==0)&(results_analysis_int2$TrialStatus==0)))
      
      eventtime <- median(results_analysis_int2$DayInfected[results_analysis_int2$eventstatus==1])
      
      results_analysis_int2$DayInfected[newevent_v_rownum] <- eventtime
      results_analysis_int2$eventstatus[newevent_v_rownum] <- 1
      
      results_analysis_int2$DayInfected[newevent_c_rownum] <- eventtime
      results_analysis_int2$eventstatus[newevent_c_rownum] <- 1
      
      survmodel<- try(ic_sp(cbind(DayEnrolled, DayInfected) ~ TrialStatus, model = 'ph',
                            bs_samples = 50, data = results_analysis_int2))      
      usesurvmod <- !inherits(survmodel, 'try-error')
      
      if (!is.na(survmodel$coefficients[1])){       
        vaccEffEst_total_int<-c(1,1,1)
        vaccEffEst_total_int_low <- 1-exp(survmodel$coefficients + 1.96*sqrt(survmodel$var))
        vaccEffEst_total_int_high <- 1-exp(survmodel$coefficients - 1.96*sqrt(survmodel$var))
        VE_fake <- 1-exp(survmodel$coefficients)
        zval <- survmodel$coefficients[1]/sqrt(survmodel$var[1])
        VEC_pval_int <- pnorm(zval, lower.tail = VE_fake>0)*2
        VEC_int_var <-survmodel$var
      } else {
        vaccEffEst_total_int<-c(1,1,1)
        vaccEffEst_total_int_low <- c(1,1,1)
        vaccEffEst_total_int_high <- c(1,1,1)
        VEC_pval_int <- NA
        VEC_int_var <-NA
      }
    } else {
      # If no events are observed in either arm, the trial has failed and no result can be obtained
      vaccEffEst_total_int<-c(NA,NA,NA)
      vaccEffEst_total_int_low <- c(NA,NA,NA)
      vaccEffEst_total_int_high <- c(NA,NA,NA)
      VEC_int_var <- NA
      VEC_pval_int <- NA
    }
    
    vaccEffEst_total_int <- c(vaccEffEst_total_int[1],vaccEffEst_total_int_low[1],vaccEffEst_total_int_high[1])
  
  
    # interval censored cox model with time to event for symptomatic and sample 3 times (or however many "num_sample" is set to) for asymptomatic
    # Note day enrolled is now representing beginning of interval during which they could be infected
    results_analysis_int_1 <- results_analysis
    results_analysis_int_1$DayInfected <- results_analysis_int_1$DayInfected + results_analysis_int_1$DayEnrolled
    
    ## Those asymptomatically infected in first third of trial -- update day infected to day of testing (Day 50 here)
    results_analysis_int_1$DayInfected[(results_analysis_int_1$Symptomatic==0) & 
                                         (results_analysis_int_1$DayInfected<=trial_length/num_sample + trial_startday)]<-trial_length/num_sample + trial_startday
    
    ## Those asymptomatically infected in second third of trial -- update day infected to day of testing (Day 100 here)
    results_analysis_int_1$DayInfected[(results_analysis_int_1$Symptomatic==0) & 
                                         (results_analysis_int_1$DayInfected<=trial_length/num_sample*2 + trial_startday) & 
                                         (results_analysis_int_1$DayInfected>trial_length/num_sample + trial_startday)]<-trial_length/num_sample*2 + trial_startday 
    ## Change their day of enrollment to day of last test (Day 50 here)
    results_analysis_int_1$DayEnrolled[(results_analysis_int_1$Symptomatic==0) & 
                                         (results_analysis_int_1$DayInfected<=trial_length/num_sample*2 + trial_startday) & 
                                         (results_analysis_int_1$DayInfected>trial_length/num_sample + trial_startday)]<-trial_length/num_sample + trial_startday 
    
    ## Those asymptomatically infected in final third of trial -- update day infected to day of testing (Day 150 here)
    results_analysis_int_1$DayInfected[(results_analysis_int_1$Symptomatic==0) & 
                                         (results_analysis_int_1$DayInfected<=trial_length/num_sample*3 + trial_startday) & 
                                         (results_analysis_int_1$DayInfected>trial_length/num_sample*2 + trial_startday)]<-trial_length/num_sample*3 + trial_startday 
    
    ## Change their day of enrollment to day of last test (Day 100 here)
    results_analysis_int_1$DayEnrolled[(results_analysis_int_1$Symptomatic==0) & 
                                         (results_analysis_int_1$DayInfected<=trial_length/num_sample*3 + trial_startday) & 
                                         (results_analysis_int_1$DayInfected>trial_length/num_sample*2 + trial_startday)]<-trial_length/num_sample*2 + trial_startday 
    
    ## Bginning of interval for symptomatically infected is the day they are infected
    results_analysis_int_1$DayEnrolled[(results_analysis_int_1$eventstatus==1) & (results_analysis_int_1$Symptomatic==1)]<- results_analysis_int_1$DayInfected[(results_analysis_int_1$eventstatus==1) & (results_analysis_int_1$Symptomatic==1)]
    
    ## Non events interval is from the end of the trial through infiinity 
    results_analysis_int_1$DayEnrolled[(results_analysis_int_1$eventstatus==0)]<- trial_length + trial_startday
    results_analysis_int_1$DayInfected[(results_analysis_int_1$eventstatus==0)]<-Inf
    
    ## Run interval censored model on this dataset
    if ((numevents_vacc>0)&&(numevents_cont>0)) {
      survmodel<- try(ic_sp(cbind(DayEnrolled, DayInfected) ~ TrialStatus, model = 'ph',
                            bs_samples = 50, data = results_analysis_int_1))
      if (!is.na(survmodel$coefficients[1])){
        # If no error was thrown, use the results of the model
        vaccEffEst_total_int_1 <- 1-exp(survmodel$coefficients)
        vaccEffEst_total_int_1_low <- 1-exp(survmodel$coefficients + 1.96*sqrt(survmodel$var))
        vaccEffEst_total_int_1_high <- 1-exp(survmodel$coefficients - 1.96*sqrt(survmodel$var))
        VEC_int_1_var <- survmodel$var
        zval <- survmodel$coefficients[1]/sqrt(survmodel$var[1])
        VEC_pval_int_1 <- pnorm(zval, lower.tail = vaccEffEst_total_int_1[1]>0)*2
      } else {
        vaccEffEst_total_int_1 <-c(NA,NA,NA)
        vaccEffEst_total_int_1_low <- c(NA,NA,NA)
        vaccEffEst_total_int_1_high <- c(NA,NA,NA)
        VEC_int_1_var <- NA
        VEC_pval_int_1 <- NA
      }
    } else if ((numevents_vacc>0)&&(numevents_cont==0)) {
      # If there are no events in the control arm but events in the vaccine arm, VE estimate is -1
      vaccEffEst_total_int_1<-c(-1,-1,-1)
      vaccEffEst_total_int_1_low <- c(-1,-1,-1)
      vaccEffEst_total_int_1_high <- c(-1,-1,-1)
      VEC_int_1_var <-NA
      VEC_pval_int_1 <-1
    } else if ((numevents_vacc==0)&&(numevents_cont>0)) {
      # If there are no events in the vaccine arm and events in the control arm, VE is 1, but
      # for the p-value we add one event to both arms and do a Cox regression on that data
      # What time to assign? I have chosen to give both events the median time
      # among control events.
      results_analysis_int_12<-results_analysis_int_1
      newevent_v_rownum <- min(which((results_analysis_int_12$eventstatus==0)&(results_analysis_int_12$TrialStatus==1)))
      newevent_c_rownum <- min(which((results_analysis_int_12$eventstatus==0)&(results_analysis_int_12$TrialStatus==0)))
      
      eventtime <- median(results_analysis_int_12$DayInfected[results_analysis_int_12$eventstatus==1])
      
      results_analysis_int_12$DayInfected[newevent_v_rownum] <- eventtime
      results_analysis_int_12$eventstatus[newevent_v_rownum] <- 1
      
      results_analysis_int_12$DayInfected[newevent_c_rownum] <- eventtime
      results_analysis_int_12$eventstatus[newevent_c_rownum] <- 1
      
      survmodel<- try(ic_sp(cbind(DayEnrolled, DayInfected) ~ TrialStatus, model = 'ph',
                            bs_samples = 50, data = results_analysis_int12))    
      usesurvmod <- !inherits(survmodel, 'try-error')
      
      if (!is.na(survmodel$coefficients[1])){
        vaccEffEst_total_int_1<-c(1,1,1)
        vaccEffEst_total_int_1_low <- 1-exp(survmodel$coefficients + 1.96*sqrt(survmodel$var))
        vaccEffEst_total_int_1_high <- 1-exp(survmodel$coefficients - 1.96*sqrt(survmodel$var))
        VE_fake <- 1-exp(survmodel$coefficients)
        zval <- survmodel$coefficients[1]/sqrt(survmodel$var[1])
        VEC_pval_int_1 <- pnorm(zval, lower.tail = VE_fake>0)*2
        VEC_int_1_var <-survmodel$var
      } else {
        vaccEffEst_total_int_1<-c(1,1,1)
        vaccEffEst_total_int_1_low <- c(1,1,1)
        vaccEffEst_total_int_1_high <- c(1,1,1)
        VEC_pval_int_1 <- NA
        VEC_int_1_var <-NA
      }
    } else {
      # If no events are observed in either arm, the trial has failed and no result can be obtained
      vaccEffEst_total_int_1<-c(NA,NA,NA)
      vaccEffEst_total_int_1_low <- c(NA,NA,NA)
      vaccEffEst_total_int_1_high <- c(NA,NA,NA)
      VEC_int_1_var <- NA
      VEC_pval_int_1 <- NA
    }
    
    vaccEffEst_total_int_1 <- c(vaccEffEst_total_int_1[1],vaccEffEst_total_int_1_low[1],vaccEffEst_total_int_1_high[1])
    
    
    ## NOT STRATIFIED
    ## Estimate vacccine efficacy using Cox PH for whole trial overall but not stratified by community
    if ((numevents_vacc>0)&&(numevents_cont>0)) {
      survmodel<-try(coxph(Surv(DayInfected,eventstatus)~TrialStatus,results_analysis),silent=T)
      usesurvmod <- !inherits(survmodel, 'try-error')
      if (!is.na(survmodel$coefficients[1])){
        vaccEffEst_total_unstrat <- 1-exp(survmodel$coefficients + c(0, 1.96, -1.96)*sqrt(survmodel$var))
        VEC_var_unstrat <-survmodel$var
        zval <- survmodel$coefficient/sqrt(survmodel$var)
        VEC_pval_unstrat <- pnorm(zval, lower.tail = vaccEffEst_total[1]>0)*2
      } else {
        vaccEffEst_total_unstrat<-c(NA,NA,NA)
        VEC_var_unstrat <-NA
        VEC_pval_unstrat <- NA
      }
    } else if ((numevents_vacc>0)&&(numevents_cont==0)) {
      # If there are no events in the control arm but events in the vaccine arm, VE estimate is -1
      vaccEffEst_total_unstrat<-c(-1,-1,-1)
      VEC_var_unstrat <-NA
      VEC_pval_unstrat <-1
    } else if ((numevents_vacc==0)&&(numevents_cont>0)) {
      # If there are no events in the vaccine arm and events in the control arm, VE is 1, but
      # for the p-value we add one event to both arms and do a Cox regression on that data
      # What time to assign? I have chosen to give both events the median time
      # among control events.
      results_analysis3<-results_analysis
      newevent_v_rownum <- min(which((results_analysis3$eventstatus==0)&(results_analysis3$TrialStatus==1)))
      newevent_c_rownum <- min(which((results_analysis3$eventstatus==0)&(results_analysis3$TrialStatus==0)))

      eventtime <- median(results_analysis3$DayInfected[results_analysis3$eventstatus==1])

      results_analysis3$DayInfected[newevent_v_rownum] <- eventtime
      results_analysis3$eventstatus[newevent_v_rownum] <- 1

      results_analysis3$DayInfected[newevent_c_rownum] <- eventtime
      results_analysis3$eventstatus[newevent_c_rownum] <- 1

      survmodel<- try(ic_sp(cbind(DayEnrolled, DayInfected) ~ TrialStatus, model = 'ph',
                            bs_samples = 50, data = results_analysis3))
      usesurvmod <- !inherits(survmodel, 'try-error')

      if (!is.na(survmodel$coefficients[1])){
        vaccEffEst_total_unstrat<-c(1,1,1)
        VE_fake <- 1-exp(survmodel$coefficients)
        zval <- survmodel$coefficient/sqrt(survmodel$var)
        VEC_pval_unstrat <- pnorm(zval, lower.tail = VE_fake>0)*2
        VEC_var_unstrat <-survmodel$var
      } else {
        vaccEffEst_total_unstrat<-c(1,1,1)
        VEC_pval_unstrat <- NA
        VEC_var_unstrat <-NA
      }
    } else {
      # If no events are observed in either arm, the trial has failed and no result can be obtained
      vaccEffEst_total_unstrat<-c(NA,NA,NA)
      VEC_var_unstrat<-NA
      VEC_pval_unstrat <-NA
    }
    
    #Imputation
    ## Take a random sample of 10% of non-symptomatic (asymptomatically infected or uninfected) in vaccine group and control group
    results_analysis_impute_int <- results_analysis_int
    
    ## Add column to record imputed infection status
    results_analysis_impute_int <- cbind(results_analysis_impute_int,rep(NA,nrow(results_analysis_impute_int)))
    names(results_analysis_impute_int) <- c("InfectedNode","DayInfected","Community","TrialStatus","Symptomatic","DayEnrolled","eventstatus","imputation")
    
    ## calculate size of 10% sample and then take half of that amount from vaccine group and half from control group
    sample <- round(nrow(results_analysis_impute_int)*sample_percent,0)
    sample_v <- round(sample/2,0)
    sample_c <- sample-sample_v
    mysample <- results_analysis_impute_int
    mysample_vacc <- mysample[mysample$TrialStatus==1,]
    mysample_vacc <- mysample_vacc[sample(1:nrow(mysample_vacc), sample_v,replace=FALSE),]   
    mysample_cont <- mysample[mysample$TrialStatus==0,]
    mysample_cont <- mysample_cont[sample(1:nrow(mysample_cont), sample_c,replace=FALSE),] 
    mysample <- rbind(mysample_vacc,mysample_cont)
    
    ## Record their imputed infection status as their actual infection status as this will be known from testing the sample
    results_analysis_impute_int$imputation[(results_analysis_impute_int$InfectedNode %in% mysample$InfectedNode)] <- results_analysis_impute_int$eventstatus[(results_analysis_impute_int$InfectedNode %in% mysample$InfectedNode)]
    
    ## Calculate numbers infected in the sample
    num_vacc_sample <- nrow(mysample[mysample$TrialStatus==1,])
    num_events_vacc_sample <- nrow(mysample[(mysample$TrialStatus==1) & (mysample$eventstatus==1),])
    num_events_vacc_sample_symp <- nrow(mysample[(mysample$TrialStatus==1) & (mysample$eventstatus==1) & (mysample$Symptomatic==1),])
    num_events_vacc_sample_asymp <- nrow(mysample[(mysample$TrialStatus==1) & (mysample$eventstatus==1) & (mysample$Symptomatic==0),])

    num_cont_sample <- nrow(mysample[mysample$TrialStatus==0,])
    num_events_cont_sample <- nrow(mysample[(mysample$TrialStatus==0) & (mysample$eventstatus==1),])
    num_events_cont_sample_symp <- nrow(mysample[(mysample$TrialStatus==0) & (mysample$eventstatus==1) & (mysample$Symptomatic==1),])
    num_events_cont_sample_asymp <- nrow(mysample[(mysample$TrialStatus==0) & (mysample$eventstatus==1) & (mysample$Symptomatic==0),])

    ## Estimate VE against progression to symptoms in the sample
    VEP_impute <- 1- (num_events_vacc_sample_symp/num_events_vacc_sample)/(num_events_cont_sample_symp/num_events_cont_sample)
    
    #Imputation
    ## calculate proportion infected in vaccine and control in sample
    rv<-num_events_vacc_sample/num_vacc_sample
    rc<-num_events_cont_sample/num_cont_sample
    
    ## Imputed infection status of the symptomatics will be their true infection status
    results_analysis_impute_int$imputation[results_analysis_impute_int$Symptomatic==1] <- 1
    
    ## Save data set before imputation so will be able to reset for multiple imputation
    results_analysis_impute_int2 <- results_analysis_impute_int
    
    ## Will impute 10 times (multiple imputation)
    Qs <- c()
    Us <- c()
    for (k in 1:10){
      results_analysis_impute_int <- results_analysis_impute_int2
      
      ## Impute infection status using rv and rc
      results_analysis_impute_int$imputation[(results_analysis_impute_int$TrialStatus==1) & is.na(results_analysis_impute_int$imputation)]<-rbinom(nrow(results_analysis_impute_int[(results_analysis_impute_int$TrialStatus==1) & is.na(results_analysis_impute_int$imputation),]),1,rv)
      results_analysis_impute_int$imputation[(results_analysis_impute_int$TrialStatus==0) & is.na(results_analysis_impute_int$imputation)]<-rbinom(nrow(results_analysis_impute_int[(results_analysis_impute_int$TrialStatus==0) & is.na(results_analysis_impute_int$imputation),]),1,rc)
      
      ## Set day infected to be end of trial for those with imputed infection status of 1 (among asymptomatics)
      results_analysis_impute_int$DayInfected[(results_analysis_impute_int$imputation==1) & ((results_analysis_impute_int$Symptomatic==0) | is.na(results_analysis_impute_int$Symptomatic))]<- trial_length + trial_startday

      ## Set interval for those with non-event imputed status to be end of trial through infinity
      results_analysis_impute_int$DayEnrolled[(results_analysis_impute_int$imputation==0)]<- trial_length + trial_startday
      results_analysis_impute_int$DayInfected[(results_analysis_impute_int$imputation==0)]<-Inf
      
      ## Run interval censored model on imputed data
      if (numevents_vacc>0 && numevents_cont>0){
        survmodel<- try(ic_sp(cbind(DayEnrolled, DayInfected) ~ TrialStatus, model = 'ph',
                              bs_samples = 50, data = results_analysis_impute_int))
        usesurvmod <- !inherits(survmodel, 'try-error')
        if (usesurvmod==TRUE){
          # If no error was thrown, use the results of the model
          vaccEffEst_total_impute_int <- 1-exp(survmodel$coefficients)
          vaccEffEst_total_impute_int_low <- 1-exp(survmodel$coefficients + 1.96*sqrt(survmodel$var))
          vaccEffEst_total_impute_int_high <- 1-exp(survmodel$coefficients - 1.96*sqrt(survmodel$var))
          VEC_impute_int_var <- survmodel$var
          zval <- survmodel$coefficients[1]/sqrt(survmodel$var[1])
          VEC_pval_impute_int <- pnorm(zval, lower.tail = vaccEffEst_total_impute_int[1]>0)*2
        } else {
          vaccEffEst_total_impute_int <-c(NA,NA,NA)
          vaccEffEst_total_impute_int_low <- c(NA,NA,NA)
          vaccEffEst_total_impute_int_high <- c(NA,NA,NA)
          VEC_impute_int_var <- NA
          VEC_pval_impute_int <- NA
        }
      } else if ((numevents_vacc>0)&&(numevents_cont==0)) {
        # If there are no events in the control arm but events in the vaccine arm, VE estimate is -1
        vaccEffEst_total_impute_int<-c(-1,-1,-1)
        vaccEffEst_total_impute_int_low <- c(-1,-1,-1)
        vaccEffEst_total_impute_int_high <- c(-1,-1,-1)
        VEC_impute_int_var <-NA
        VEC_pval_impute_int <-1
      } else if ((numevents_vacc==0)&&(numevents_cont>0)) {
        # If there are no events in the vaccine arm and events in the control arm, VE is 1, but
        # for the p-value we add one event to both arms and do a Cox regression on that data
        # What time to assign? I have chosen to give both events the median time
        # among control events.
        results_analysis_impute_int2<-results_analysis_impute_int
        newevent_v_rownum <- min(which((results_analysis_impute_int2$eventstatus==0)&(results_analysis_impute_int2$TrialStatus==1)))
        newevent_c_rownum <- min(which((results_analysis_impute_int2$eventstatus==0)&(results_analysis_impute_int2$TrialStatus==0)))
        
        eventtime <- median(results_analysis_impute_int2$DayInfected[results_analysis_impute_int2$eventstatus==1])
        
        results_analysis_impute_int2$DayInfected[newevent_v_rownum] <- eventtime
        results_analysis_impute_int2$eventstatus[newevent_v_rownum] <- 1
        
        results_analysis_impute_int2$DayInfected[newevent_c_rownum] <- eventtime
        results_analysis_impute_int2$eventstatus[newevent_c_rownum] <- 1
        
        survmodel<- try(ic_sp(cbind(DayEnrolled, DayInfected) ~ TrialStatus, model = 'ph',
                              bs_samples = 50, data = results_analysis_impute_int2))
        usesurvmod2 <- !inherits(survmodel, 'try-error')
        
        if (usesurvmod2==TRUE){
          vaccEffEst_total_impute_int<-c(1,1,1)
          vaccEffEst_total_impute_int_low <- 1-exp(survmodel$coefficients + 1.96*sqrt(survmodel$var))
          vaccEffEst_total_impute_int_high <- 1-exp(survmodel$coefficients - 1.96*sqrt(survmodel$var))
          VE_fake <- 1-exp(survmodel$coefficients)
          zval <- survmodel$coefficients[1]/sqrt(survmodel$var[1])
          VEC_pval_impute_int <- pnorm(zval, lower.tail = VE_fake>0)*2
          VEC_impute_int_var <-survmodel$var
        } else {
          vaccEffEst_total_impute_int<-c(1,1,1)
          vaccEffEst_total_impute_int_low <- c(1,1,1)
          vaccEffEst_total_impute_int_high <- c(1,1,1)
          VEC_pval_impute_int <- NA
          VEC_impute_int_var <-NA
        }
      } else {
        # If no events are observed in either arm, the trial has failed and no result can be obtained
        vaccEffEst_total_impute_int<-c(NA,NA,NA)
        vaccEffEst_total_impute_int_low <- c(NA,NA,NA)
        vaccEffEst_total_impute_int_high <- c(NA,NA,NA)
        VEC_impute_int_var <- NA
        VEC_pval_impute_int <- NA
      }
      
      vaccEffEst_total_impute_int <- c(vaccEffEst_total_impute_int[1],vaccEffEst_total_impute_int_low[1],vaccEffEst_total_impute_int_high[1])
      
      # U is average variance of imputed estimates
      # m is number of times impute
      # B is between imputation variance
      # https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2727536/
      Qs <- c(Qs,vaccEffEst_total_impute_int[1])
      Us <- c(Us,VEC_impute_int_var[1])
    }
    
    Q <- mean(Qs,na.rm=TRUE)
    U <- mean(Us,na.rm=TRUE)
    Between <- c()
    for (r in 1:10){
      Between<- c(Between,(Qs[r]-Q)^2)
    }
    B <- 1/(10-1)*sum(Between)
    T = U + ((10+1)/10)*B
    vaccEffEst_total_impute_int<-Q
    vaccEffEst_total_impute_int_lower<-Q-1.96*sqrt(T)
    vaccEffEst_total_impute_int_upper<-Q+1.96*sqrt(T)
    VEC_impute_int_var<-T
  
    
    
    
  } else{
    num_vacc<- NA          
    numevents_vacc<- NA          
    numevents_vacc_S<- NA          
    numevents_vacc_AS<- NA          
    num_cont<- NA          
    numevents_cont<- NA          
    numevents_cont_S<- NA          
    numevents_cont_AS<- NA          
    sample_size<- NA          
    sample_sizeP<- NA          
    sample_sizeS<- NA          
    VEP_estimate<- NA          
    VES_estimate<- NA          
    VEC_estimate<- NA          
    VEC_corrected<- NA          
    vaccEffEst_total[1]<- NA          
    vaccEffEst_total[2]<- NA          
    vaccEffEst_total[3]<- NA          
    vaccEffEst_total_S[1]<- NA          
    vaccEffEst_total_S[2]<- NA          
    vaccEffEst_total_S[3]<- NA          
    VEC_pval<- NA           
    VES_pval<- NA           
    VEC_pval_int<- NA           
    VEC_pval_int_1<- NA           
    vaccEffEst_total_int[1]<- NA          
    vaccEffEst_total_int[2]<- NA           
    vaccEffEst_total_int[3]<- NA          
    vaccEffEst_total_int_1[1]<- NA          
    vaccEffEst_total_int_1[2]<- NA         
    vaccEffEst_total_int_1[3]<- NA          
    vaccEffEst_total_impute_int<- NA           
    vaccEffEst_total_impute_int_low<- NA         
    vaccEffEst_total_impute_int_high<- NA           
    VEC_var<- NA           
    VES_var<- NA           
    VEC_int_var<- NA          
    VEC_int_1_var<- NA           
    VEC_impute_int_var<- NA           
    VEC_pval_impute_int <- NA
  }
  
  names <- c("num_vacc","numevents_vacc","numevents_vacc_S","numevents_vacc_AS",
               "num_cont","numevents_cont","numevents_cont_S","numevents_cont_AS",
               "sample_size","sample_sizeP","sample_sizeS",
               "VEP_estimate","VES_estimate",
               "VEC_estimate","VEC_corrected",
               "vaccEffEst_total","vaccEffEst_total_low","vaccEffEst_total_high",
               "vaccEffEst_total_S","vaccEffEst_total_S_low","vaccEffEst_total_S_high",
               "VEC_pval", "VES_pval", "VEC_pval_int", "VEC_pval_int_1", 
               "vaccEffEst_total_int", "vaccEffEst_total_int_low", "vaccEffEst_total_int_high",
               "vaccEffEst_total_int_1","vaccEffEst_total_int_1_low","vaccEffEst_total_int_1_high",
               "vaccEffEst_total_impute_int",  "vaccEffEst_total_impute_int_low", "vaccEffEst_total_impute_int_high", 
               "VEC_var","VES_var","VEC_int_var","VEC_int_1_var", "VEC_impute_int_var", "VEC_pval_impute_int")
  
  results <- c(num_vacc,numevents_vacc,numevents_vacc_S,numevents_vacc_AS,
       num_cont,numevents_cont,numevents_cont_S,numevents_cont_AS,
       sample_size,sample_sizeP,sample_sizeS,
       VEP_estimate,VES_estimate,
       VEC_estimate,VEC_corrected,
       vaccEffEst_total[1],vaccEffEst_total[2],vaccEffEst_total[3],
       vaccEffEst_total_S[1],vaccEffEst_total_S[2],vaccEffEst_total_S[3],
       VEC_pval, VES_pval, VEC_pval_int, VEC_pval_int_1, 
       vaccEffEst_total_int[1], vaccEffEst_total_int[2], vaccEffEst_total_int[3],
       vaccEffEst_total_int_1[1],vaccEffEst_total_int_1[2],vaccEffEst_total_int_1[3],
       vaccEffEst_total_impute_int, vaccEffEst_total_impute_int_lower,vaccEffEst_total_impute_int_upper, 
       VEC_var, VES_var, VEC_int_var,VEC_int_1_var, VEC_impute_int_var, VEC_pval_impute_int)

  trialresults<- rbind(names,results)
  write.csv(trialresults, paste0("1comm_constant_trial_results_",cluster_coverage,"_",trial_length,"_",j,"_",betaAS,"_",vpS,"_",rate_within,"_",sample_percent,".csv"))

  
