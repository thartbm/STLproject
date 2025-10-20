

getLongExpFits <- function() {
  
  
  
  participant <- c()
  maxrot      <- c()
  phase       <- c()
  lambda      <- c()
  N0          <- c()
  
  for (group in c(45, 60, 90)) {
    
    IDs <- getidentifiers(group=as.character(group))
    
    for (IDno in c(1:length(IDs))) {
      
      ID <- IDs[IDno]
      
      cat(sprintf('working on participant %s (%d/%d)\n', ID, IDno, length(IDs)))
      
      df <- loadPerformanceFile(group=as.character(group), id=ID)
      # df  <- read.csv(sprintf('data/%s_performance.csv',ID), stringsAsFactors = FALSE)
      # mrot <- getMaxRot(df)
      ldf <- extractLongExpData(df)
      fits <- fitLongExposure(ldf)
      # print(fits)
      
      for (mode in c('learning','washout')) {
        participant <- c(participant, ID)
        maxrot      <- c(maxrot, group)
        phase       <- c(phase, mode)
        lambda      <- c(lambda, fits[[mode]]['lambda'])
        N0          <- c(N0, fits[[mode]]['N0'])
      }
      
    }
    
  }
  
  df <- data.frame( participant,
                    maxrot,
                    phase,
                    lambda,
                    N0             )
  
  write.csv(df, 'data/exponential_fits.csv', row.names = FALSE, quote = FALSE)
  
}

fitLongExposure <- function(data) {
  
  lfit <- Reach::exponentialFit(  signal=data[['exposure']],
                                  mode='learning',
                                  asymptoteRange = c(10,30))
  # print(lfit)
  
  wfit <- Reach::exponentialFit( signal=data[['aftereffects']],
                                 mode='washout',
                                 asymptoteRange = c(0,30))
  # print(wfit)
  
  return( list( 'learning' = lfit,
                'washout'  = wfit  )   )
  
}



getPredictionData <- function() {
  
  maxrots=c(45,60,90)
  
  expfits <- read.csv('data/exponential_fits.csv', stringsAsFactors = FALSE)
  expfits <- expfits[which(expfits$phase == 'learning'),]
  expfits$pred <- expfits$lambda * expfits$N0
  
  outdfs <- list()
  
  for (maxrot in maxrots) {
    
    # groupparticipants <- getidentifiers(group=as.character(maxrot))
    
    # print(length(groupparticipants))
    
    mfits <- read.csv(sprintf('data/modelFits_%d.csv',maxrot), stringsAsFactors = FALSE)
    # mfits <- mfits[which(mfits$participant %in% groupparticipants),]
    
    
    participant <- c()
    target      <- c()
    attribution <- c()
    capped      <- c()
    exponential <- c()
    
    participants <- setdiff(unique(mfits$participant), c('all'))
    
    for (ID in participants) {
      
      for (tt in c('point','arc')) { # why do we include the arc targets here? we're not using that anywhere....
        
        mpfits <- mfits[which(mfits$participant == ID & mfits$target == tt),]
        
        # capped model fit:
        mpar <- c('r'=mpfits$r,
                  'c'=mpfits$c)
        capped_pred <- cappedModel(par=mpar, rotations=c(20))
        
        # # attribution model fit:
        # mpar <- c('s'=mpfits$s,
        #           'w'=mpfits$w)
        # attr_pred <- STLPredict(par=mpar, rotations=c(20))
        
        # exponential fits:
        exppred <- expfits$pred[which(expfits$participant == ID)]
        if (length(exppred) == 0) {
          print(ID)
        }

        participant <- c(participant, ID)
        target      <- c(target,      tt)
        # attribution <- c(attribution, attr_pred)
        capped      <- c(capped,      capped_pred)
        exponential <- c(exponential, exppred)
        
      }

    }

    outdfs[[sprintf('%d', maxrot)]] <- data.frame( participant,
                                                   target,
                                                   # attribution,
                                                   capped,
                                                   exponential)
    
  }
  
  return(outdfs)
  
}



# long exposure function -----

extractLongExpData <- function(df) {
  
  baseline_idx <- which(df$label == 'fixed-rotation-baseline')
  exposure_idx <- which(df$label == 'fixed-rotation')
  aftereffects_idx <- which(df$label == 'fixed-rotation-aftereffects')
  
  # print(c(length(baseline_idx),length(exposure_idx),length(aftereffects_idx)))
  
  baseline <- df$reachdeviation_deg[baseline_idx[31:60]]
  baseline <- baseline[which(abs(baseline) < 60)]
  baseline <- mean(baseline)
  # print(baseline)
  
  exposure <- df$reachdeviation_deg[exposure_idx] - baseline
  exposure[which(abs(exposure) > 80)] <- NA
  
  aftereffects <- df$reachdeviation_deg[aftereffects_idx] - baseline
  aftereffects[which(abs(aftereffects) > 80)] <- NA
  
  return( list( 'baseline'      = baseline,
                'exposure'      = exposure,
                'aftereffects'  = aftereffects) )
  
}

getLongTimeCourses <- function() {
  
  if (file.exists('data/longTimecourses.csv')) {
    
    df <- read.csv('data/longTimecourses.csv', stringsAsFactors = FALSE)
    dfs <- list('45' = df[which(df$maxrot == 45),],
                '60' = df[which(df$maxrot == 60),],
                '90' = df[which(df$maxrot == 90),])
    return(dfs)
    
  } else {
    
    maxrots = c(45,60,90)
    
    dfs <- list('45' = NA,
                '60' = NA,
                '90' = NA)
    
    for (mrot in maxrots) {
      
      IDs <- getidentifiers(group=as.character(mrot))
      
      for (IDno in c(1:length(IDs))) {
        
        ID <- IDs[IDno]
        
        pdf <- loadPerformanceFile(group=as.character(mrot), id=ID)
        ldf  <- extractLongExpData(pdf)
        
        pdf <- data.frame( participant = rep(ID,160),
                           trial       = c(1:160),
                           reachdeviation_deg = ldf[['exposure']],
                           maxrot      = rep(mrot,160) )
        
        idx <- sprintf('%d', mrot)
        
        if (is.data.frame(dfs[[idx]])) {
          dfs[[idx]] <- rbind(dfs[[idx]], pdf)
        } else {
          dfs[[idx]] <- pdf
        }
        
      }
      
    }
    
    write.csv( rbind(dfs[['45']], dfs[['60']], dfs[['90']]),
               'data/longTimecourses.csv',
               row.names = FALSE,
               quote = TRUE)
    
    return(dfs)
    
  }
  
}

getFirstTrialData <- function() {
  
  # we need the exponential fit data:
  exponentials <- read.csv('data/exponential_fits.csv', stringsAsFactors = FALSE)
  exponentials <- exponentials[which(exponentials$phase == 'learning'),]
  exponentials$prediction <- exponentials$lambda * exponentials$N0
  
  # and the long time course data itself:
  exposure <- getLongTimeCourses()
  exposure <- rbind(exposure[['45']], exposure[['60']], exposure[['90']])
  # exposure <- exposure[,-which(names(exposure) %in% c('maxrot'))]
  
  exposure <- exposure[which(exposure$trial == 2),]
  
  firstTrials <- merge(exponentials, exposure, by=c('participant', 'maxrot'))
  
  return(firstTrials)
  
}