
`%notin%` <- Negate(`%in%`)

bothModelFits <- function() {
  
  for (maxrot in c(45,60,90)) {  
    
    df <- loadSTLdata(maxrots=c(maxrot))
    
    
    df <- aggregate(response ~ participant + target + rotation, data=df, FUN=median)
    
    participants <- unique(df$participant)
    
    participant <- c()
    target <- c()
    n <- c()
    r_par <- c()
    c_par <- c()
    cappedMSE <- c()
    s_par <- c()
    w_par <- c()
    attributionMSE <- c()
    
    for (targ in c('arc', 'point')) {
      
      subdf <- df[which(df$target == targ),]
      
      for (ppid in c('all', participants)) {
        
        if (ppid == 'all') {
          ppid_idx <- 1
        } else {
          ppid_idx <- which(participants==ppid) + 1
        }
        
        cat(sprintf('group %d; target %s; participant %s (%d/%d)\n', maxrot, targ, ppid, ppid_idx, length(participants)+1))
        
        if (ppid == 'all') {
          pplist <- participants
        } else {
          pplist <- c(ppid)
        }
        
        capped_par      <- cappedFit( participants=pplist, data=subdf)
        attribution_par <- STLFit(    participants=pplist, data=subdf)
        
        participant <- c(participant, ppid)
        target <- c(target, targ)
        n <- c(n, length(which(subdf$participant %in% pplist)))
        r_par <- c(r_par, capped_par['r'])
        c_par <- c(c_par, capped_par['c'])
        s_par <- c(s_par, attribution_par['s'])
        w_par <- c(w_par, attribution_par['w'])
        cappedMSE <- c(cappedMSE, cappedMSE(par        = capped_par,
                                            rotations  = subdf$rotation,
                                            deviations = subdf$response) )
        
        attributionMSE <- c(attributionMSE, STLErrors(par        = attribution_par,
                                                      rotations  = subdf$rotation,
                                                      deviations = subdf$response) )
        
        
      }
      
    }
    
    mf <- data.frame(participant,
                     target,
                     n,
                     'r'=r_par,
                     'c'=c_par,
                     cappedMSE,
                     's'=s_par,
                     'w'=w_par,
                     attributionMSE)
    
    write.csv(mf, file=sprintf('data/modelFits_%d.csv',maxrot), quote=TRUE, row.names=FALSE)
    
  }
  
}


# getModelLikelihoods <- function() {
#   
#   
#   mf45 <- read.csv('data/modelFits_45.csv', stringsAsFactors = FALSE)
#   mf45$maxrotation <- 45
#   mf60 <- read.csv('data/modelFits_60.csv', stringsAsFactors = FALSE)
#   mf60$maxrotation <- 60
#   
#   mf <- rbind(mf45, mf60)
#   
#   mf$cappedLL <- NA
#   mf$attributionLL <- NA
#   
#   for (maxrot in c(45,60)) {
#     
#     participants <- unique(mf$participant[which(mf$maxrotation == maxrot)])
#     
#     for (target in c('point', 'arc')) {
#       
#       for (ppid in participants) {
#         
#         idx <- which(mf$participant == ppid & mf$target == target & mf$maxrotation == maxrot)
#         
#         MSEs <- c('capped'=mf$cappedMSE[idx], 'attribution'=mf$attributionMSE[idx])
#         N <- mf$n[idx]
#         
#         AICs <- AICc(MSE = MSEs,
#                      k   = c(2,2),
#                      N   = N)
#         
#         rLL <- relativeLikelihood(AICs)
#         
#         mf$cappedLL[idx] <- rLL['capped']
#         mf$attributionLL[idx] <- rLL['attribution']
#       }
#     }
#   }
#   
#   return(mf)
#   
# }

# not sure what's going on in this function: leaving it out for now

# compareModelFits <- function() {
#   
#   mf <- getModelLikelihoods()
#   
#   for (maxrot in c(45,60)) {
#     
#     participants <- unique(mf$participant[which(mf$participant != 'all' & mf$maxrotation == maxrot)])
#     
#     for (target in c('point', 'arc')) {
#       
#       idx <- which(mf$participant == 'all' & mf$target == target & mf$maxrotation == maxrot)
#       cat(sprintf('\nModel fit on all participants for %s targets with a max rotation of %d:\n',toupper(target),maxrot))
#       
#       cat('-- relative likelihoods:\n\n')
#       rLL <- c('capped'      = mf$cappedLL[idx],
#                'attribution' = mf$attributionLL[idx])
#       print(rLL)
#       
#       pval <- max(rLL[which(rLL < 1.0)])
#       
#       cat(sprintf('\n  this means the %s model is best (p = %0.3f)\n', names(rLL)[which.max(rLL)], pval))
#       
#       countwins <- c('capped'=0, 'attribution'=0)
#       for (ppid in participants) {
#         
#         idx <- which(mf$participant == 'all' & mf$target == target & mf$maxrotation == maxrot)
#         
#         rLL <- c('capped'      = mf$cappedLL[idx],
#                  'attribution' = mf$attributionLL[idx])
#         
#         countwins[names(rLL)[which.max(rLL)]] <- 1 + countwins[names(rLL)[which.max(rLL)]]
#         
#       }
#       
#       cat('\n-- count of participants for which each model is best:\n\n')
#       print(countwins)
#       
#     }
#   }
# }


# model evaluation functions ---

# These functions should only be used for models that are  fitted by
# minimizing errors, rather than maximizing likelihood.

# stats::AIC() uses standard R model objects with likelihoods

# AIC <- function(MSE, k, N) {
#   return( (N * log(MSE)) + (2 * k) )
# }
# 
# AICc <- function(MSE, k, N) {
#   return( AIC(MSE, k, N) * (((2*k^2) + 2*k) / (N - k - 1)) )
# }
# 
# relativeLikelihood <- function(crit) {
#   return( exp( ( min( crit  ) - crit  ) / 2 ) )
# }

# statistics on raw data -----

taskErrorANOVAs <- function() {
  
  for (maxrot in c(45,60,90)) {
    
    cat(sprintf('\n=== %d° GROUP\n\n',maxrot))
    
    df <- loadSTLdata(average=median,maxrots=c(maxrot))
    
    df$participant <- as.factor(df$participant)
    df$rotation <- as.factor(df$rotation)
    df$target <- as.factor(df$target)
    df$maxrotation <- as.factor(df$maxrotation)
    
    
    print(afex::aov_ez(
      id = 'participant',
      dv = 'response',
      data = df,
      within = c('target', 'rotation')
      # within = c('target', 'rotation', 'superblock')
    ))
    
  }
  
}

taskErrorPostHocs <- function(verbose=FALSE, dfout=TRUE) {
  
  group <- c()
  rotation <- c()
  pval <- c()
  respdiff <- c()
  
  for (maxrot in c(45,60,90)) {
    
    if (verbose) { cat(sprintf('\n=== %d° GROUP\n\n',maxrot)) }
    
    df <- loadSTLdata(average=median,maxrots=c(maxrot))
    
    df$participant <- as.factor(df$participant)
    df$rotation <- as.factor(df$rotation)
    df$target <- as.factor(df$target)
    df$maxrotation <- as.factor(df$maxrotation)
    
    for (rot in sort(unique(df$rotation))) {
      
      rdf  <- df[which(df$rotation == rot),]
      
      ardf <- rdf[which(rdf$target=='arc'),]
      names(ardf)[which(names(ardf)=='response')] <- 'arc_response'
      ardf <- ardf[,which(names(ardf) %notin% c('target','rotation','maxrotation'))]
      # str(ardf)
      
      prdf <- rdf[which(rdf$target=='point'),]
      names(prdf)[which(names(prdf)=='response')] <- 'point_response'
      prdf <- prdf[,which(names(prdf) %notin% c('target','rotation','maxrotation'))]
      # str(prdf)
      
      pardf <- merge(x = ardf, y=prdf, by='participant')
      # print(str(pardf))
      rtt <- t.test( x      = pardf$arc_response,
                     y      = pardf$point_response,
                     paired = TRUE)
      if (verbose) {
        cat(sprintf('\n=== %d GROUP: %d rotation\n', maxrot, as.integer(rot)))
        print(rtt)
      }
      
      group <- c(group, maxrot)
      rotation <- c(rotation, rot)
      pval <- c(pval, rtt$p.value)
      respdiff <- c(respdiff, rtt$estimate)
      
    }
    
  }
  
  if (dfout) {
    
    respdiff <- -1 * respdiff
    df <- data.frame(group, rotation, pval, respdiff)
    
    for (group in unique(df$group)) {
      idx <- which(df$group == group)
      df$pval.adj[idx] <- p.adjust(df$pval[idx], method='fdr')
    }
    
    return(df)
    
  } else {
    
    NULL
    
  }
  
}

# delayed feedback -----

terminalFeedbackTtest <- function() {
  
  df <- read.csv('data/FDaftereffects.csv', stringsAsFactors = FALSE)
  
  df <- df[which(df$delay == 0),]
  df <- aggregate(aftereffect ~ participant, data=df, FUN=mean)
  
  print(t.test(df$aftereffect, mu=0))
  
}

delayedFeedbackANOVA <- function() {
  
  df <- read.csv('data/FDaftereffects.csv', stringsAsFactors = FALSE)
  
  df$participant <- as.factor(df$participant)
  df$rotation <- as.factor(df$rotation)
  df$delay <- as.factor(df$delay)
  
  print(afex::aov_ez(
    id = 'participant',
    dv = 'aftereffect',
    data = df,
    within = c('delay', 'rotation'),
    type=3
  ))
  
}

# model comparison -----

modelLikelihood <- function(maxrots=c(45,60,90)) {
  
  MSEs <- c('attenuation'=0, 'capped'=0)
  groupN <- c()
  N <- 0
  
  all_mfs <- NA
  for (maxrot in maxrots) {
    mf <- read.csv(sprintf('data/modelFits_%d.csv', maxrot))
    mf <- mf[which(mf$participant %notin% c('all') & mf$target == 'point'),]
    if (is.data.frame(all_mfs)) {
      all_mfs <- rbind(all_mfs, mf)
    } else {
      all_mfs <- mf
    }
  }
  MSEs['attenuation'] <- MSEs['attenuation'] + sum(all_mfs$attributionMSE)
  MSEs['capped']      <- MSEs['capped']      + sum(all_mfs$cappedMSE)
  groupN <- c(groupN, length(all_mfs$participant))
  N <- N + sum(all_mfs$n)

  # print(N)
  # cat('sum of MSEs:\n')
  # print(MSEs)
  # 
  # AICs <- Reach::AICc(MSEs, N=N, k=c(2,2))
  # 
  # cat('\ncorrected AICs (AICc):\n')
  # print(AICs)
  # 
  # cat('\nrelative log-likelihoods:\n')
  # print(Reach::relativeLikelihood(AICs))
  # 
  # print(groupN)
  
  cat('mean of MSEs:\n')
  MSEs <- MSEs/sum(groupN)
  print(MSEs)
  
  AICs <- Reach::AIC(MSEs, N=sum(groupN), k=c(2,2))
  
  cat('\nMSE-based AICs:\n')
  print(AICs)
  
  cat('\nAIC-based relative log-likelihoods:\n')
  print(Reach::relativeLikelihood(AICs))
  
  
}

modelTtests <- function(maxrots=c(45,60,90)) {
  
  MSEs <- list('attribution'=c(), 'capped'=c())
  
  df <- NA
  for (maxrot in maxrots) {
    mf <- read.csv(sprintf('data/modelFits_%d.csv', maxrot))
    mf <- mf[which(mf$participant %notin% c('all')),]
    if (is.data.frame(df)) {
      df <- rbind(df, mf)
    } else {
      df <- mf
    }
  }
  
  print( BayesFactor::ttestBF(  x=df$cappedMSE,
                                y=df$attributionMSE,
                                paired=TRUE)              )
  
  print( t.test(  x=df$cappedMSE,
                  y=df$attributionMSE,
                  paired=TRUE)              )
  
}

comparePredictions <- function() {
  
  predictions <- getPredictionData()
  
  
  for (maxrot in c(45,60,90)) {
    
    preds <- predictions[[sprintf('%d',maxrot)]]
    ppreds <- preds[which(preds$target == 'point'),]
    
    linfits <- list()
    
    X = ppreds[,'capped']
    Y = ppreds$exponential
    
    capmod <- lm(Y~X+0)
    
    cat(sprintf('=== %d° group: CAPPED model\n',maxrot))
    
    # get the p-values and (adjusted) R-squared values, from here:
    print(summary(capmod))

    X = ppreds[,'attribution']
    Y = ppreds$exponential
    
    attrmod <- lm(Y~X+0)
    
    cat(sprintf('=== %d° group: ATTENUATION model\n',maxrot))
    
    # get the p-values and (adjusted) R-squared values, from here:
    print(summary(attrmod))
    
        
    # print(summary(linmod)$p.value)
    # print(summary(linmod)$adj.r.squared)
    
    
    print(stats::AIC(capmod, attrmod))
    
  }
  
}

# compare predicted and actual first long exposure trials -----

firstTrialTest <- function () {

  df <- getFirstTrialData()

  print(t.test(df$prediction, df$reachdeviation_deg, paired=TRUE))

}

firstTrialRegression <- function() {
  
  pd <- getPredictionData()
  predictions <- rbind(pd[['45']], pd[['60']], pd[['90']])
  
  ppreds <- predictions[which(predictions$target == 'point'),]
  # 
  X = ppreds$capped
  Y = ppreds$exponential
  
  # plot linear fit with intercept=0
  linmod <- lm(Y~X+0)
  
  print(summary(linmod))
  
}

# test-restest reliability:

doTestRetestRegressions <- function() {
  
  
  allPar <- read.csv('data/testRetest_fits.csv', stringsAsFactors = F)
  
  # for (parameter in c('r','c','s','w')) {
  for (parameter in c('r','c')) {
    
    
    cat(paste0(list('r'='CAPPED FIXED-RATE: RATE',
                    'c'='CAPPED FIXED-RATE: CAP', 
                    's'='ATTENUATION: SLOPE',
                    'w'='ATTENUATION: WIDTH')[[parameter]],'\n'))
    
    
    run1val <- allPar[which(allPar$epoch == 'early'),parameter]
    run2val <- allPar[which(allPar$epoch == 'late' ),parameter]
    
    
    lmfit <- lm(run2val ~ run1val)
    print(summary(lmfit))
    
  }
  
}