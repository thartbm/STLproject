


getErrorResponses <- function() {
  
  allErrors <- NA
  
  for (group in c('45','60','90')) {
    
    IDs <- getidentifiers(group)
    
    for (ID in IDs) {
      
      df <- loadPerformanceFile(group, ID)
      
      df$error <- df$reachdeviation_deg + df$rotation
      
      tidx <- which(df$label == 'stl-target-rotation')
      # print(length(tidx))
      
      ep2  <- df$error[tidx+2]
      ep1  <- df$error[tidx+1]
      e_0  <- df$error[tidx]
      em1  <- df$error[tidx-1]
      em2  <- df$error[tidx-2]
      rotation <- df$rotation[tidx]
      
      errorResponses <- data.frame(
        ID       = ID,
        group    = group,
        rotation = rotation,
        em2      = em2,
        em1      = em1,
        e_0      = e_0,
        ep1      = ep1,
        ep2      = ep2
      )
      
      if (is.data.frame(allErrors)) {
        allErrors <- rbind(allErrors, errorResponses)
      } else {
        allErrors <- errorResponses
      }
    }
    
  }
  
  write.csv(allErrors, file = 'data/errorResponses.csv', row.names = FALSE, quote=TRUE)
  
}


errorPlot <- function() {
  
  df <- read.csv('data/errorResponses.csv', stringsAsFactors = FALSE)
  
  
  layout(mat=matrix(data=c(1:3), nrow=3, ncol=1, byrow=TRUE))
  
  par(mar=c(4.1,4.1,2,0.1))
  
  
  for (group in c(45,60,90)) {
    
    cat(sprintf('\n=== GROUP: %d\n', group))
    
    gdf <- df[df$group == group, ]
    
    plot(NULL,NULL,
         xlim=c(-100,100),
         ylim=c(-30,30),
         xlab='Error on trial n (deg)', ylab='Error on trial n+1 (deg)',
         asp=1,ax=F,bty='n'
         )
    
    points(x=gdf$e_0, y=gdf$ep1, pch=16, col=rgb(0,0,0,0.1))
    
    lines(x=c(-100,100), y=c(0,0), col='red', lty=1, lw=1)
    
    X <- seq(-group,group,0.1)
    Y <- ksmooth(x=gdf$e_0, y=gdf$ep1, kernel = "normal", bandwidth = 10, x.points=X)
    
    lines(X,Y$y, col='blue', lty=1, lw=2)
    
    axis(side=1, at=seq(-group,+group,15))
    axis(side=2, at=seq(-30,30,15))
    
  }
  
  
}


errorHistory <- function() {
  
  df <- read.csv('data/errorResponses.csv', stringsAsFactors = FALSE)
  
  for (group in c(45,60,90)) {
    
    cat(sprintf('\n=== GROUP: %d\n', group))
    
    gdf <- df[df$group == group, ]
    
    HE_lm <- lm(ep1 ~ rotation + em2 + em1 + e_0:rotation, data=gdf)
    
    print(summary(HE_lm))
    
    # print(car::vif(HE_lm, type='predictor'))
    
  }
  
}

plotWashoutErrors <- function() {
  
  df <- read.csv('data/errorResponses.csv', stringsAsFactors = FALSE)
  
  layout(mat=matrix(data=c(1:3), nrow=3, ncol=1, byrow=TRUE))
  
  par(mar=c(4.1,4.1,2,0.1))
  
  errorbased <- FALSE
  
  for (group in c(45,60,90)) {
    
    cat(sprintf('\n=== GROUP: %d\n', group))
    
    gdf <- df[df$group == group, ]
    
    plot(NULL,NULL,
         xlim=c(-100,100),
         ylim=c(-10,10),
         xlab='Rotation on trial n (deg)', ylab='Error on trial n+2 (deg)',
         ax=F,bty='n'
    )
    
    if (errorbased) {
      
      points(x=gdf$e_0, y=gdf$ep2, pch=16, col=rgb(0,0,0,0.1))
      
      lines(x=c(-100,100), y=c(0,0), col='red', lty=1, lw=1)
      
      X <- seq(-group,group,0.1)
      Y <- ksmooth(x=gdf$e_0, y=gdf$ep2, kernel = "normal", bandwidth = 10, x.points=X)
      
      lines(X,Y$y, col='blue', lty=1, lw=2)
      
    } else {
      
      lines(x=c(-100,100), y=c(0,0), col='gray', lty=1, lw=1)
      
      gdf$WO <- gdf$ep2 - gdf$em1 # baseline the first washout trials
      
      CI <- aggregate(WO ~ rotation, data=gdf, FUN=Reach::getConfidenceInterval, method='b', resamples=2000)
      avg <- aggregate(WO ~ rotation, data=gdf, FUN=mean)
      
      colors <- c('#0066FFFF', '#0066FF33')
      
      X <- c(CI$rotation, rev(CI$rotation))
      Y <- c(CI$WO[,1], rev(CI$WO[,2]))
      polygon(X,Y,border=NA,col=colors[2])
      
      lines(x=avg$rotation, y=avg$WO, col=colors[1], lw=1)
      
      ppdata <- aggregate(WO ~ ID + rotation, data=gdf, FUN=mean)
      print(afex::aov_ez(
        id = 'ID',
        dv = 'WO',
        data = ppdata,
        within = c('rotation')
      ))
      
    }
    
    axis(side=1, at=seq(-group,+group,15))
    axis(side=2, at=seq(-10,10,5))
    
  }
  
  
}

plotNextBaseline <- function() {
  
  df <- read.csv('data/errorResponses.csv', stringsAsFactors = FALSE)
  
  layout(mat=matrix(data=c(1:3), nrow=3, ncol=1, byrow=TRUE))
  par(mar=c(4.1,4.1,2,0.1))
  
  cat('can not baseline this data !!!\n')
  
  for (group in c(45,60,90)) {
    
    cat(sprintf('\n=== GROUP: %d\n', group))
    
    gdf <- df[df$group == group, ]
    
    plot(NULL,NULL,
         xlim=c(-100,100),
         ylim=c(-10,10),
         xlab='Rotation on previous bout (deg)', ylab='Error on current baseline (deg)',
         ax=F,bty='n'
    )
    
    
    lines(x=c(-100,100), y=c(0,0), col='gray', lty=1, lw=1)
    
    NBdf <- NA
    
    for (ID in unique(gdf$ID)) {
      
      idf <- gdf[gdf$ID == ID, ]
      
      rotations    <- idf$rotation[c(1:99)]
      nextBaseline <- idf$em1[c(2:100)]
      
      idf <- data.frame(
        ID       = ID,
        rotation = rotations,
        baseline = nextBaseline
      )

      if (is.data.frame(NBdf)) {
        NBdf <- rbind(NBdf, idf)
      } else {
        NBdf <- idf
      }
    }
    
    
  
    CI <- aggregate(baseline ~ rotation, data=NBdf, FUN=Reach::getConfidenceInterval, method='b', resamples=2000)
    avg <- aggregate(baseline ~ rotation, data=NBdf, FUN=mean)

    colors <- c('#0066FFFF', '#0066FF33')

    X <- c(CI$rotation, rev(CI$rotation))
    Y <- c(CI$baseline[,1], rev(CI$baseline[,2]))
    polygon(X,Y,border=NA,col=colors[2])

    lines(x=avg$rotation, y=avg$baseline, col=colors[1], lw=1)



    axis(side=1, at=seq(-group,+group,15))
    axis(side=2, at=seq(-10,10,5))
    
    
    
    ppdata <- aggregate(baseline ~ ID + rotation, data=NBdf, FUN=mean)
    print(afex::aov_ez(
      id = 'ID',
      dv = 'baseline',
      data = ppdata,
      within = c('rotation')
    ))
    
  }
  
}