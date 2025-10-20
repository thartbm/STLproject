
# max w at: 20
# NOT: 25, 30, 32, 35, 40, 45

testRetestFits <- function() {
  
  group        <- c()
  participant  <- c()
  epoch        <- c()
  r            <- c()
  c            <- c()
  s            <- c()
  w            <- c()
  
  for (rot in c(45,60,90)) {
    
    X <- seq(0, rot, 0.1)
    
    target <- 'point'
    df <- read.csv(sprintf('data/%sTargetSTL_%d.csv', target, rot))
    

    participants <- unique(df$participant)
    
    for (ppid in participants) {
      
      cat(sprintf('group %d; participant %s (%d/%d)\n', rot, ppid, which(participants==ppid), length(participants)))
      
      dfp <- df[df$participant==ppid,]
    
      for (sbpair in c('early', 'late')) {
        
        sbs <- list('early'=c(1,2), 'late'=c(4,5))[[sbpair]]
        
        dfe <- dfp[dfp$superblock %in% sbs,]
        # dfe <- aggregate(response ~ rotation, data=dfe, FUN=median)
        # plot(x=dfe$rotation, y=dfe$response, main=sprintf('P%s %s', ppid, sbpair))
        
        Cpars <- cappedFit(participants=ppid, data=dfe)
        Apars <- STLFit(   participants=ppid, data=dfe)
        
        # lines(x=X,
        #       y=cappedModel(par = Cpars,
        #                     rotations = X),
        #       col='red')
        # lines(x=X,
        #       y=STLPredict(par = Apars,
        #                    rotations = X),
        #       col='blue')
        
        # all_parameters$group[row_idx] <- group
        # all_parameters$participant[row_idx] <- ppid
        # all_parameters$epoch[row_idx] <- epoch
        # all_parameters$r[row_idx] <- fit['r']
        # all_parameters$c[row_idx] <- fit['c']
        
        group <- c(group, rot)
        participant <- c(participant, ppid)
        epoch <- c(epoch, sbpair)
        r <- c(r, Cpars['r'])
        c <- c(c, Cpars['c'])
        s <- c(s, Apars['s'] )
        w <- c(w, Apars['w'] )
        
      }
      
    }
    
  }
  
  results <- data.frame(group=group, participant=participant, epoch=epoch, r=r, c=c, s=s, w=w)
  write.csv(results, 'data/testRetest_fits.csv', row.names = FALSE)
  
}

