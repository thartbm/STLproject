

getLongExpFits <- function() {
  
  IDs <- getidentifiers()
  
  participant <- c()
  maxrot      <- c()
  phase       <- c()
  lambda      <- c()
  N0          <- c()
  
  
  for (IDno in c(1:length(IDs))) {
    
    ID <- IDs[IDno]
    
    cat(sprintf('working on participant %s (%d/%d)\n', ID, IDno, length(IDs)))
    
    df  <- read.csv(sprintf('data/%s_performance.csv',ID), stringsAsFactors = FALSE)
    mrot <- getMaxRot(df)
    ldf <- extractLongExpData(df)
    fits <- fitLongExposure(ldf)
    # print(fits)
    
    for (mode in c('learning','washout')) {
      participant <- c(participant, ID)
      maxrot      <- c(maxrot, mrot)
      phase       <- c(phase, mode)
      lambda      <- c(lambda, fits[[mode]]['lambda'])
      N0          <- c(N0, fits[[mode]]['N0'])
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



getPredictionData <- function(maxrots=c(45,60,90)) {
  
  allparticipants <- getidentifiers()
  
  expfits <- read.csv('data/exponential_fits.csv', stringsAsFactors = FALSE)
  expfits <- expfits[which(expfits$phase == 'learning'),]
  expfits$pred <- expfits$lambda * expfits$N0
  
  outdfs <- list()
  
  for (maxrot in maxrots) {
    
    mfits <- read.csv(sprintf('data/modelFits_%d.csv',maxrot), stringsAsFactors = FALSE)
    mfits <- mfits[which(mfits$participant %in% allparticipants),]
    
    participant <- c()
    target      <- c()
    attribution <- c()
    capped      <- c()
    exponential <- c()
    
    participants <- unique(mfits$participant)
    
    for (ID in participants) {
      
      for (tt in c('point','arc')) {
        
        mpfits <- mfits[which(mfits$participant == ID & mfits$target == tt),]
        
        # capped model fit:
        mpar <- c('r'=mpfits$r,
                  'c'=mpfits$c)
        capped_pred <- cappedModel(par=mpar, rotations=c(20))
        
        # attribution model fit:
        mpar <- c('s'=mpfits$s,
                  'w'=mpfits$w)
        attr_pred <- STLPredict(par=mpar, rotations=c(20))
        
        exppred <- expfits$pred[which(expfits$participant == ID)]
        if (length(exppred) == 0) {
          print(ID)
        }

        participant <- c(participant, ID)
        target      <- c(target,      tt)
        attribution <- c(attribution, attr_pred)
        capped      <- c(capped,      capped_pred)
        exponential <- c(exponential, exppred)
        
      }

    }

    outdfs[[sprintf('%d', maxrot)]] <- data.frame( participant,
                                                   target,
                                                   attribution,
                                                   capped,
                                                   exponential)
    
  }
  
  return(outdfs)
  
}