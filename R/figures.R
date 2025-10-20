
# only for delayed feedback inset plot:

getColors <- function() {
  
  cols.op <- c(rgb(255, 147, 41,  255, max = 255), # orange:  21, 255, 148
               rgb(229, 22,  54,  255, max = 255), # red:    248, 210, 126
               rgb(207, 0,   216, 255, max = 255), # pink:   211, 255, 108
               rgb(127, 0,   216, 255, max = 255), # violet: 195, 255, 108
               rgb(0,   19,  136, 255, max = 255)) # blue:   164, 255, 68
  
  cols.tr <- c(rgb(255, 147, 41,  32,  max = 255), # orange:  21, 255, 148
               rgb(229, 22,  54,  32,  max = 255), # red:    248, 210, 126
               rgb(207, 0,   216, 32,  max = 255), # pink:   211, 255, 108
               rgb(127, 0,   216, 32,  max = 255), # violet: 195, 255, 108
               rgb(0,   19,  136, 32,  max = 255)) # blue:   164, 255, 68
  
  cols <- list()
  cols$op <- cols.op
  cols$tr <- cols.tr
  
  return(cols)
  
}


fig23_data <- function(models=FALSE, target='inline') {
  
  delays <- FALSE
  if (models == FALSE) {
    filename <- 'doc/fig2_alldata'
    width    <- 5
    height   <- 5.5
  }
  if (models == TRUE) {
    filename <- 'doc/fig3_model_fits'
    width    <- 5
    height   <- 5.5
  }

  Reach::setupFigureFile(
    target = target,
    filename = filename,
    width = 5,
    height = 5.5,
    dpi=300
  )
  
  layout(mat=matrix(c(1:3), ncol=1))
  par(mar=c(3.25,3.25,2,0.1))
  
  for (maxrot in c(45,60,90)) {
    
    plot(x=-1000,y=-1000,
         main=sprintf('%d° max rotation',maxrot),ylab='',xlab='',
         # xlim=c(0,maxrot+1),
         xlim=c(0,91),
         ylim=c(-0.5,8.5),
         ax=F,bty='n')

    title(ylab='reach aftereffect [°]',line=2.25)
    title(xlab='rotation [°]',line=2.25)
    
    df <- loadSTLdata(average=median, maxrots=c(maxrot))
    
    if (models) {
      fits <- read.csv(sprintf('data/modelFits_%d.csv', maxrot), stringsAsFactors = FALSE)
      fits <- fits[which(fits$participant=='all'),]
    }

    CIdf <- aggregate(response ~ target + rotation, data=df, FUN=Reach::getConfidenceInterval, method='b', resamples=2000)
    avgdf <- aggregate(response ~ target + rotation, data=df, FUN=mean)
    
    targettypes <- c('point','arc')
    if (models == TRUE) {
      targettypes <- c('point')
    }
    
    for (tt in targettypes) {
      
      if (tt == 'point') {
        colors <- c('#0066FFFF', '#0066FF33')
      }
      if (tt == 'arc') {
        colors <- c('#FF6600FF', '#FF660033')
      }
      
      tCI <- CIdf[which(CIdf$target == tt),]
      tavg <- avgdf[which(avgdf$target == tt),]
      
      X <- c(tCI$rotation, rev(tCI$rotation))
      Y <- c(tCI$response[,1], rev(tCI$response[,2]))
      polygon(X,Y,border=NA,col=colors[2])
      
      lines(x=tavg$rotation, y=tavg$response, col=colors[1], lw=1)
      
      if (models & (tt == 'point')) {
        fit <- fits[which(fits$target == tt),]
        # cross <- fit$c/fit$r
        
        lines(x=c(0,fit$c/fit$r,maxrot),
              y=c(0,fit$c,fit$c),
              lty=2,
              col=colors[1],
              lw=2)
        
        attpar = c('s'=fit$s, 'w'=fit$w)
        rots = seq(0,maxrot,0.5)
        attribution_predictions <- STLPredict(par=attpar,
                                              rotations=rots)

        lines(x=rots,
              y=attribution_predictions,
              lty=3,
              col=colors[1],
              lw=2)
      } else {
        # lines(x=tavg$rotation, y=tavg$response, col=colors[1], lw=2)
      }
      
    }
    
    if (maxrot == 45) {
      axis(side=1,at=c(1,5,10,15,20,25,30,35,40,45),cex.axis=0.75)
    }
    if (maxrot == 60) {
      axis(side=1,at=c(1,5,10,15,20,25,30,40,50,60),cex.axis=0.75)
    }
    if (maxrot == 90) {
      axis(side=1,at=c(1,5,10,15,20,30,40,50,70,90),cex.axis=0.75)
    }
    axis(side=2,at=c(0,2,4,6,8),cex.axis=0.75)
    
    
    if (maxrot == 45) {
      if (models) {
        legend(x=50,y=6,
               legend=c('data', 'capped fixed-rate model', 'attenuation model'),
               lwd=c(1,1,1),
               lty=c(1,2,3),
               col=c('#0066FFFF','#0066FFFF','#0066FFFF'),
               bty='n')
      } else {
        legend(x=50,y=6,
               legend=c('point target', 'arc target'),
               lwd=c(2,2),
               col=c('#0066FFFF','#FF6600FF'),
               bty='n')
      }
    }
    if (maxrot == 60) {
      
    }
    
    if (maxrot == 45) {
      if (delays) {
        
        colors <- getColors()
        
        fd <- read.csv('data/FDaftereffects.csv', stringsAsFactors = FALSE)
        
        delays = sort(unique(fd$delay))
        
        for (dn in c(1:length(delays))) {

          delay = delays[dn]
          
          fdd <- fd[which(fd$delay == delay),]
          
          avg <- aggregate(aftereffect ~ rotation, data=fdd, FUN=mean)
          CI  <- aggregate(aftereffect ~ rotation, data=fdd, FUN=Reach::getConfidenceInterval)
          
          polygon( x      = c(CI$rotation, rev(CI$rotation)) + 45,
                   y      = c(CI$aftereffect[,1], rev(CI$aftereffect[,2])),
                   col    = colors$tr[dn], 
                   border = NA,
                   xpd=TRUE)
          
          lines(x=avg$rotation + 45, y=avg$aftereffect, col=colors$op[dn])
          
          # print(CI$aftereffect[,1])
        }
        
        lines( x = c(48, 92),
              y = c(0, 0),
              col = '#999999',
              lty = 1,
              xpd=TRUE)
        
        axis( side = 1,
              at = unique(fd$rotation)+45,
              labels = sprintf('%d', unique(fd$rotation) ),
              cex.axis=0.75)
        # axis( side = 2,
        #       at = seq(-2,6,2))
        
        legend( x = 75,
                y = 10,
                # title='delay:',
                legend = sprintf('%0.1f s', delays),
                lty=1,
                # seg.len = 1.5,
                col=colors$op,
                bty='n',
                xpd=TRUE)
        
        
      }
    }
  }
  
  if (target %in% c('pdf','svg','png','tiff')) {
    dev.off()
  }
  
}


# plotDataParmin <- function(models=FALSE, target='inline') {
#   
#   width  = 8
#   height = 3
#   dpi    = 300
#   if (models) {
#     filename=sprintf('doc/fig4_model_fits_parmin.%s',target)
#   } else {
#     filename=sprintf('doc/fig2_alldata_parmin.%s',target)
#   }
#   
#   if (target == 'pdf') {
#     pdf(file   = filename, 
#         width  = width, 
#         height = height)
#   }
#   if (target == 'svg') {
#     svglite::svglite( filename = filename,
#                       width = width,
#                       height = height,
#                       fix_text_size = FALSE)
#   }
#   if (target == 'png') {
#     png( filename = filename,
#          width = width*dpi,
#          height = height*dpi,
#          res = dpi
#     )
#   }
#   if (target == 'tiff') {
#     tiff( filename = filename,
#           compression = 'lzw',
#           width = width*dpi,
#           height = height*dpi,
#           res = dpi
#     )
#   }
#   
#   left  <- 0.75
#   right <- 0.05
#   
#   deg15 <- (width - ((left + right) * 2)) / 7
#   
#   layout(mat=matrix(c(1:2), ncol=2), widths=c(deg15*3,deg15*4))
#   par(mai=c(0.75,left,0.4,right))
#   
#   
#   
#   
#   
#   for (maxrot in c(45,60)) {
#     
#     plot(x=-1000,y=-1000,
#          main=sprintf('%d° max rotation',maxrot),ylab='',xlab='',
#          xlim=c(0,maxrot+1),
#          # xlim=c(0,91),
#          ylim=c(-0.5,8.5),
#          ax=F,bty='n')
#     
#     title(ylab='initial change [°]',line=2.25)
#     title(xlab='rotation [°]',line=2.25)
#     
#     df <- loadSTLdata(average=median, maxrots=c(maxrot))
#     
#     if (models) {
#       fits <- read.csv(sprintf('data/modelFits_%d.csv', maxrot), stringsAsFactors = FALSE)
#       fits <- fits[which(fits$participant=='all'),]
#     }
#     
#     CIdf <- aggregate(response ~ target + rotation, data=df, FUN=Reach::getConfidenceInterval, method='b', resamples=2000)
#     avgdf <- aggregate(response ~ target + rotation, data=df, FUN=mean)
#     for (tt in c('point','arc')) {
#       
#       if (tt == 'point') {
#         colors <- c('#0066FFFF', '#0066FF33')
#       }
#       if (tt == 'arc') {
#         colors <- c('#FF6600FF', '#FF660033')
#       }
#       
#       tCI <- CIdf[which(CIdf$target == tt),]
#       tavg <- avgdf[which(avgdf$target == tt),]
#       
#       X <- c(tCI$rotation, rev(tCI$rotation))
#       Y <- c(tCI$response[,1], rev(tCI$response[,2]))
#       polygon(X,Y,border=NA,col=colors[2])
#       
#       if (models) {
#         fit <- fits[which(fits$target == tt),]
#         # cross <- fit$c/fit$r
#         
#         lines(x=c(0,fit$c/fit$r,maxrot),
#               y=c(0,fit$c,fit$c),
#               lty=3,
#               col=colors[1],
#               lw=2)
#         
#         attpar = c('s'=fit$s, 'w'=fit$w)
#         rots = seq(0,maxrot,0.5)
#         attribution_predictions <- STLPredict(par=attpar,
#                                               rotations=rots)
#         
#         lines(x=rots,
#               y=attribution_predictions,
#               lty=2,
#               col=colors[1],
#               lw=2)
#       } else {
#         lines(x=tavg$rotation, y=tavg$response, col=colors[1], lw=2)
#       }
#       
#     }
#     
#     if (maxrot == 45) {
#       axis(side=1,at=c(1,5,10,15,20,25,30,35,40,45),cex.axis=0.75)
#     }
#     if (maxrot == 60) {
#       axis(side=1,at=c(1,5,10,15,20,25,30,40,50,60),cex.axis=0.75)
#     }
#     if (maxrot == 90) {
#       axis(side=1,at=c(1,5,10,15,20,30,40,50,70,90),cex.axis=0.75)
#     }
#     axis(side=2,at=c(0,2,4,6,8),cex.axis=0.75)
#     
#     
#     if (maxrot == 60) {
#       legend(x=35,y=2,
#              legend=c('dot target', 'arc target'),
#              lwd=c(2,2),
#              col=c('#0066FFFF','#FF6600FF'),
#              bty='n')
#     }
#     if (maxrot == 60) {
#       if (models) {
#         legend(x=20,y=2,
#                legend=c('attenuation model', 'capped model'),
#                lwd=c(2,2),
#                lty=c(2,3),
#                col=c('#999999','#999999'),
#                bty='n')
#       }
#     }
#     
#   }
#   
#   if (target %in% c('pdf','svg','png','tiff')) {
#     dev.off()
#   }
#   
# }
# 
# 
# fig3_TEdiffs <- function(target='inline',posthocs=TRUE) {
#   
#   Reach::setupFigureFile(
#     target = target,
#     filename = 'doc/fig3_taskerror',
#     width = 4,
#     height = 8,
#     dpi=300
#   )
#   
#   layout(mat=matrix(c(1:3), ncol=1))
#   par(mar=c(3.25,3.25,2,0.1))
#   
#   pval.df <- taskErrorPostHocs()
#   
#   for (maxrot in c(45,60,90)) {
#     
#     plot(x=-1000,y=-1000,
#          main=sprintf('%d° max rotation',maxrot),ylab='',xlab='',
#          # xlim=c(0,maxrot+1),
#          xlim=c(0,91),
#          ylim=c(-2.5,4.5),
#          ax=F,bty='n')
#     
#     title(ylab='reach aftereffect [°]',line=2.25)
#     title(xlab='rotation [°]',line=2.25)
#     
#     df <- loadSTLdata(average=median, maxrots=c(maxrot))
#     
#     diffdf <- aggregate(response ~ participant + rotation, data=df, FUN=diff)
#     diffdf$response <- diffdf$response
#     
#     CIdf <- aggregate(response ~ rotation, data=diffdf, FUN=Reach::getConfidenceInterval, method='b', resamples=2000)
#     avgdf <- aggregate(response ~ rotation, data=diffdf, FUN=mean)
#     
#     colors <- c('#9900FFFF', '#9900FF33')
#     
#     lines(x=c(1,maxrot),y=c(0,0),col='#999999',lty=2)
#     
#     # tCI <- CIdf[which(CIdf$target == target),]
#     # tavg <- avgdf[which(avgdf$target == target),]
#     tCI <- CIdf
#     tavg <- avgdf
#     
#     X <- c(tCI$rotation, rev(tCI$rotation))
#     Y <- c(tCI$response[,1], rev(tCI$response[,2]))
#     polygon(X,Y,border=NA,col=colors[2])
#     
#     lines(x=tavg$rotation, y=tavg$response, col=colors[1])
#     
#     # add FDR corrected p-values:
#     mr.pval.df <- pval.df[which(pval.df$group == maxrot),]
#     
#     if (posthocs) {
#       for (rotation in mr.pval.df$rotation) {
#         
#         pval <- mr.pval.df$pval[which(mr.pval.df$rotation == rotation)]
#         pval.adj <- mr.pval.df$pval.adj[which(mr.pval.df$rotation == rotation)]
#         
#         pch <- 1
#         col <- '#9900FF33'
#         if (pval < 0.05) {
#           col <- '#9900FFFF'
#         }
#         if (pval.adj < 0.05) {
#           pch <- 16
#         }
#         
#         points( x = rotation,
#                 y = -2,
#                 pch=pch,
#                 col=col)
#         
#       }
#     }
#     
#     if (maxrot == 45) {
#       axis(side=1,at=c(1,5,10,15,20,25,30,35,40,45),cex.axis=0.75)
#     }
#     if (maxrot == 60) {
#       axis(side=1,at=c(1,5,10,15,20,25,30,40,50,60),cex.axis=0.75)
#     }
#     if (maxrot == 90) {
#       axis(side=1,at=c(1,5,10,15,20,30,40,50,70,90),cex.axis=0.75)
#     }
#     axis(side=2,at=c(-2,0,2,4),cex.axis=0.75)
#     
#     
#   }
#   
#   if (target %in% c('pdf','svg','png','tiff')) {
#     dev.off()
#   }
#   
#   
# }

fig3_learningPrediction <- function(target='inline') {
  
  Reach::setupFigureFile(
    target = target,
    filename = 'doc/fig3_prediction',
    width = 6,
    height = 6,
    dpi=300
  )
  
  layout(mat=matrix(c(1,1,2,3),ncol=2,byrow=TRUE),widths=c(1,1.25))
  
  par(mar=c(3.25,3.25,3.25,0.1))
  
  ltc <- getLongTimeCourses()
  timecourses <- rbind(ltc[['45']], ltc[['60']], ltc[['90']])
  
  pd <- getPredictionData()
  predictions <- rbind(pd[['45']], pd[['60']], pd[['90']])
  
    
  plot(NULL,NULL,
       main='20° time course', xlab='', ylab='',
       xlim=c(0,161),ylim=c(-5,25),
       bty='n',ax=F)
  title(xlab='trial',line=2.25)
  title(ylab='reach deviation [°]',line=2.25)
  
  # colt <- c('45'='#FF000033', '60'='#9900FF33', '90'='#0099FF33')[sprintf('%d',maxrot)]
  # colo <- c('45'='#FF0000ff', '60'='#9900FFff', '90'='#0099FFFF')[sprintf('%d',maxrot)]
  
  # colt <- c('45'='#FF000033', '60'='#9900FF33', '90'='#0099FF33')[sprintf('%d',maxrot)]
  # colo <- c('45'='#FF0000ff', '60'='#9900FFff', '90'='#0099FFFF')[sprintf('%d',maxrot)]
  
  colo <- '#0066FFFF'
  colt <- '#0066FF33'
  
  df <- timecourses
  participants <- unique(df$participant)
  
  CIdf <- aggregate(reachdeviation_deg ~ trial, data=df, FUN=Reach::getConfidenceInterval, method='t', resamples=1000)
  
  X <- c(c(1:160), rev(c(1:160)))
  Y <- c(CIdf$reachdeviation_deg[,1], rev(CIdf$reachdeviation_deg[,2]))
  polygon(X,Y,border=NA,col=colt)
  
  avg <- aggregate(reachdeviation_deg ~ trial, data=df, FUN=mean, na.rm=TRUE)
  lines(avg, col=colt)
  
  # exponential fit to average timecourse?
  expfit <- Reach::exponentialFit(signal=avg$reachdeviation_deg,
                                  mode='learning',
                                  asymptoteRange = c(10,30))
  
  expvals <- Reach::exponentialModel(par=expfit, timepoints=160)
  # print(expvals)
  lines(x=c(1:160), y=expvals$output, col=colo, lty=2)
  
  
  axis(side=1,at=c(1,40,80,120,160),cex.axis=0.8)
  axis(side=2,at=c(0,10,20),cex.axis=0.8)
  
  ppreds <- predictions[which(predictions$target == 'point'),]
  
  
  
  plot(NULL,NULL,
       main='predictions', xlab='', ylab='',
       xlim=c(0,15),ylim=c(0,15),
       bty='n',ax=F,asp=1)
  title(xlab='STL prediction [°]',line=2.25)
  title(ylab='exponential fit [°]',line=2.25)
  
  lines(x=c(0,15),y=c(0,15),col='#99999999')
  
  X = ppreds$capped
  Y = ppreds$exponential
  
  # plot data points:
  points( x = X,
          y = Y,
          pch=16, cex=1.5,
          col=colt)
  
  # plot linear fit with intercept=0
  linmod <- lm(Y~X+0)
  coeff <- linmod$coeff['X']
  lines( x=c(0,15),
         y=c(0,15*coeff),
         col=colo)
  
  axis(side=1,at=seq(0,15,5),cex.axis=0.8)
  axis(side=2,at=seq(0,15,5),cex.axis=0.8)
  
  
  df <- getFirstTrialData()
  
  
  plot(NULL,NULL,
       main='',xlab='',ylab='',
       xlim=c(0.5,3.5),ylim=c(-35,25),
       bty='n',ax=F)
  
  lines( x=c(0.5,3.5), y=c(0,0), lty=2, col='#999999' )
  lines( x=c(0.5,3.5), y=c(20,20), lty=2, col='#999999' )
  
  title(main='second trials')
  title(ylab='(reach) deviation [°]', line=2.25)

  scol <- '#0066FFFF'
  col <- '#0066FF33'
  
  for (depvar in c(1,2,3)) {
    
    if (depvar == 1) { vals <- df$reachdeviation_deg }
    if (depvar == 2) { vals <- df$prediction }
    if (depvar == 3) { vals <- df$prediction-df$reachdeviation_deg  }
    
    points(x = seq(depvar-0.3, depvar-0.1, length.out=dim(df)[1]),
           # x=rep(depvar-0.3,dim(df)[1]),
           y=vals,
           col=col)
    
    CI <- Reach::getConfidenceInterval(vals,method='b')
    
    polygon( x = c(depvar+0.1,depvar+0.3,depvar+0.3,depvar+0.1),
             y = c(CI[1],CI[1],CI[2],CI[2]),
             col=col,
             border=NA)
    
    lines(x=c(depvar+0.1,depvar+0.3),
          y=rep(mean(vals),2),
          col=scol)
    
  }
    
  axis(side=1,at=c(1,2,3),labels=c('reach\ndata', 'exponential\nprediction', 'difference\n '), padj=0.5, cex.axis=0.8)
  axis(side=2,at=seq(-20,20,10), cex.axis=0.8)
  
  if (target %in% c('pdf','svg','png','tiff')) {
    dev.off()
  }
  
}

# fig7_1stTrialFits <- function(target='inline') {
#   
#   Reach::setupFigureFile(
#     target = target,
#     filename = 'doc/fig7_exp1st',
#     width = 7,
#     height = 4,
#     dpi=300
#   )
#   
#   layout(mat=matrix(c(1:3),ncol=3,byrow=TRUE))
#   par(mar=c(3.25,3.25,3.25,0.1))
# 
#   df <- getFirstTrialData()
#   
#   plot(NULL,NULL,
#        main='',xlab='',ylab='',
#        xlim=c(0.5,3.5),ylim=c(-35,25),
#        bty='n',ax=F)
#   
#   lines( x=c(0.5,3.5), y=c(0,0), lty=2, col='#999999' )
#   lines( x=c(0.5,3.5), y=c(20,20), lty=2, col='#999999' )
#   
#   title(main='second trials')
#   title(ylab='(reach) deviation [°]', line=2.25)
#   # df <- firstTrials[which(firstTrials$maxrot == maxrot),]
#   
#   # col <- c('45'='#FF000033', '60'='#9900FF33', '90'='#0099FF33')[sprintf('%d',maxrot)]
#   # scol <- c('45'='#FF0000', '60'='#9900FF', '90'='#0099FF')[sprintf('%d',maxrot)]
#   
#   scol <- '#0066FFFF'
#   col <- '#0066FF33'
#   
#   for (depvar in c(1,2,3)) {
#     
#     if (depvar == 1) { vals <- df$reachdeviation_deg }
#     if (depvar == 2) { vals <- df$prediction }
#     if (depvar == 3) { vals <- df$prediction-df$reachdeviation_deg  }
#     
#     points(x = seq(depvar-0.3, depvar-0.1, length.out=dim(df)[1]),
#            # x=rep(depvar-0.3,dim(df)[1]),
#            y=vals,
#            col=col)
#     
#     CI <- Reach::getConfidenceInterval(vals,method='b')
#     
#     polygon( x = c(depvar+0.1,depvar+0.3,depvar+0.3,depvar+0.1),
#              y = c(CI[1],CI[1],CI[2],CI[2]),
#              col=col,
#              border=NA)
#     
#     lines(x=c(depvar+0.1,depvar+0.3),
#           y=rep(mean(vals),2),
#           col=scol)
#     
#     
#     
#     axis(side=1,at=c(1,2,3),labels=c('reach\ndata', 'exponential\nprediction', 'difference\n '), padj=0.5, cex.axis=0.75)
#     axis(side=2,at=seq(-20,20,10), cex.axis=0.75)
#     
# 
#     # print(range(df$reachdeviation_deg))
#     
#     # print(t.test(df$prediction, df$reachdeviation_deg, paired=TRUE))
#     
#   }
#   
#   if (target %in% c('pdf','svg','png','tiff')) {
#     dev.off()
#   }
#   
# }




fig4_retest <- function(target='inline') {
  
  
  Reach::setupFigureFile(
    target = target,
    filename = 'doc/fig4_retest',
    width = 8,
    height = 4.5,
    dpi=300
  )
  

  allPar <- read.csv('data/testRetest_fits.csv', stringsAsFactors = F)
  
  layout(mat=matrix(data=c(1,2),ncol=2, byrow=TRUE))
  
  par(mar=c(4.1,4.1,3.25,0.2))
  
  # colt <- c('45'='#FF000033', '60'='#9900FF33', '90'='#0099FF33')[sprintf('%d',maxrot)]
  # colo <- c('45'='#FF0000ff', '60'='#9900FFff', '90'='#0099FFFF')[sprintf('%d',maxrot)]
  
  
  for (parameter in c('r','c')) {
    
    lim <- list('r'=c(0,1), 'c'=c(0,16), 's'=c(0,1), 'w'=c(0,40))[[parameter]]
    
    plot(NULL, NULL, 
         xlab=sprintf('%s (block 1 & 2)', list('r'='rate', 'c'='cap', 's'='slope', 'w'='width')[[parameter]]),
         ylab=sprintf('%s (block 4 & 5)', list('r'='rate', 'c'='cap', 's'='slope', 'w'='width')[[parameter]]),
         xlim=lim, ylim=lim,
         asp=1, bty='n', ax=F,
         las=2
    )
    
    title(main=list('r'='A: capped fixed-rate: rate', 'c'='B: capped fixed-rate: cap', 's'='C: attentuation: slope', 'w'='D: attenuation: width')[[parameter]], line=2.25)
    
    lines(x=lim, y=lim, col='#999999')
    
    
    for (rot in c(45,60,90)) {
      
      run1val <- allPar[which(allPar$epoch == 'early' & allPar$group == rot),parameter]
      run2val <- allPar[which(allPar$epoch == 'late'  & allPar$group == rot),parameter]
      
      colt <- c('45'='#FF000033', '60'='#9900FF33', '90'='#0099FF33')[sprintf('%d',rot)]
      colo <- c('45'='#FF0000ff', '60'='#9900FFff', '90'='#0099FFFF')[sprintf('%d',rot)]
      
      points( x = run1val,
              y = run2val,
              pch=16, cex=1.6,
              col=colt)
    
    }
    
    run1val <- allPar[which(allPar$epoch == 'early'),parameter]
    run2val <- allPar[which(allPar$epoch == 'late' ),parameter]
    
    
    lmfit <- lm(run2val ~ run1val)
    # print(summary(lmfit))
    
    if (parameter == 'c') {
      
      # colors <- c('#FF6600FF', '#FF660033')   # orange
      # "#00CED1FF", "#00CED133"                # turquoise
      
      at <- range(run1val)
      

      coef <- lmfit$coefficients
      lines(at, coef[1]+(at*coef[2]), col="#000000", lty=2)
      
      ci <- predict( lmfit,
                     newdata=data.frame(run1val=seq(at[1],at[2],length.out=40)),
                     interval = "confidence")
      
      X <- c(seq(at[1],at[2],length.out=40),rev(seq(at[1],at[2],length.out=40)))
      Y <- c(ci[,'lwr'],rev(ci[,'upr']))
      polygon(x=X,y=Y,col="#00000033",border=NA)
      
    }
    
    axis(side=1,at=seq(lim[1],lim[2],length.out=list('r'=6, 'c'=5, 's'=6, 'w'=5)[[parameter]]),cex.axis=1.0)
    axis(side=2,at=seq(lim[1],lim[2],length.out=list('r'=6, 'c'=5, 's'=6, 'w'=5)[[parameter]]),cex.axis=1.0)
    
  }
  
  if (target %in% c('pdf','svg','png','tiff')) {
    dev.off()
  }
  
}