

prepDemographics <- function() {
  
  # STL_90:
  # demo_90 <- read.csv('data_prep/demographics_STL_90.csv')
  # names(demo_90) <- c('datetime',
  #                     'consent',
  #                     'age',
  #                     'sex',
  #                     'eyesight',
  #                     'handedness',
  #                     'ID')
  
  
  dd <- read.csv('data_prep/SRC_demographics_STL_delayed.csv')
  names(dd)   <-    c('datetime',
                      'consent',
                      'age',
                      'sex',
                      'vision',
                      'handedness',
                      'participant')
  
  dd <- dd[,c('participant','consent','age','sex','vision','handedness')]
  
  dd$vision[which(dd$vision == "\"corrected to normal\" (I need contacts or glasses and I'm wearing them now)")] <- 'corrected'
  dd$vision[which(dd$vision == "\"normal\" (I don't need contacts or glasses)")] <- 'normal'
  
  
  dd$handedness[which(dd$handedness == "right handed")] <- 'right'
  dd$handedness[which(dd$handedness == "left handed")] <- 'left'
  
  write.csv(dd, 'data_prep/demographics_STL_delayed.csv', row.names=F)
  
}



prepGroupData <- function(group) {
  
  # locate sources (hard coded... because it's not super well-organized):
  folder <- list('45'      = 'data_prep/raw/data_45_60',
                 '60'      = 'data_prep/raw/data_45_60',
                 '90'      = 'data_prep/raw/data_90',
                 'delayed' = 'data_prep/raw/STL_delayed')[[group]]
  
  demographics_file <- list('45'      = 'data_prep/demographics_STL_45_60.csv',
                            '60'      = 'data_prep/demographics_STL_45_60.csv',
                            '90'      = 'data_prep/demographics_STL_45_60.csv',
                            'delayed' = 'data_prep/demographics_STL_delayed.csv')[[group]] 
  
  # get participants who are in the demographics file and who have data files:
  # (for the 45 and 60 degree groups, we will also need to verify the max rotation)
  
  dd <- read.csv(demographics_file)
  
  ppdirs <- list.dirs(folder, recursive=F, full.names=F)
  
  pp <- intersect(dd$participant, ppdirs)
  
  if (group %in% c('45','60')) {
    
    ppids <- c()
    
    for (ppid in pp) {
      
      sdf <- read.csv(sprintf('%s/%s/%s_performance.csv', folder, ppid, ppid), stringsAsFactors = F)

      if (max(sdf$rotation) == as.numeric(group)) {
        ppids <- c(ppids, ppid)
      }
      
    }
    
  } else {
    ppids <- pp
  }
  
  # now we have a list of ppids in the group:
  # move files over!
  
  dir.create(sprintf('data_prep/nice/STL_%s', group))
  dir.create(sprintf('data_prep/nice/STL_%s/full', group))
  dir.create(sprintf('data_prep/nice/STL_%s/short', group))
  
  # for (ppid in ppids) {
  #   # copy the 'full folder' to the "full" folder
  #   from <- sprintf('%s/%s', folder, ppid)
  #   to   <- sprintf('data_prep/nice/STL_%s/full', group)
  #   
  #   file.copy(from      = from,
  #             to        = to,
  #             recursive = TRUE)
  #   
  #   # then move the summary file to the "short" folder
  #   
  #   from <- sprintf('%s/%s/%s_performance.csv', to, ppid, ppid)
  #   to   <- sprintf('data_prep/nice/STL_%s/short/%s_performance.csv', group, ppid)
  #   # print(from)
  #   # print(to)
  #   file.copy(from      = from,
  #             to        = to)
  # }
  
  # if this all went allright, we should now create the demographics file:
  
  dd <- dd[which(dd$participant %in% ppids),]
  
  cols <- c('participant', 'consent', 'age', 'sex', 'handedness', 'vision')
  if (group %in% c('45', '60')) {cols <- c(cols, 'height_cm')}
  
  dd <- dd[,cols]
  
  for (colname in c('sex', 'handedness', 'vision')) {
    dd[,colname] <- tolower(dd[,colname])
  }
  
  write.csv(dd, sprintf('data_prep/nice/demographics_%s.csv', group), quote=T, row.names=F)
  
}

prepAllData <- function() {
  
  for (group in c('45', '60', '90', 'delayed')) {
    prepGroupData(group=group)
  }
  
}