# Alex on the 18/06/2018
# Script to calcultae the lmortality as a rate, record various summary and database 

SavingInfo = "MyDFs(dir='bureau' or 'home',
CODE = 'BETPEN',
seuil = 0.8
seuilC = 0.6)
This function extract plot data in three tables (one with management, one scaled and one non scaled. Mortality binary data are added, as well as SPEI indexes that were calculated previously"
bannerBreak = "\n*********************************************************************\n"
cat(paste0(bannerBreak,SavingInfo,bannerBreak,"\n"))

MyDFs <- function(dir="bureau",
                  CODE,
                  seuil = 0.8,
                  seuilC = 0.6) {
  if (dir == "home") {
    Dir = c("/Users/alexandrechangenet/Dropbox/FUNDIV/")
  } else if (dir == "bureau") {
    Dir <- c("/home/achangenet/Documents/FUNDIV - NFI - Europe/")
  }
  setwd(Dir)
  Dir = c(paste0(Dir, "our-data/species/", CODE, "/CLIMAP"))
  setwd(Dir)
  list.files(Dir, pattern = paste0(".rds"))
  df <- readRDS(paste0("Mydf2_", CODE, "_", seuil, "_", seuilC, ".rds")) #Base de données plot full
  ###################################################################
  ###                                                             ###
  ###   Here is the df (scaled plots but some are managed)        ###
  ###                                                             ###
  ###################################################################
  
  dfplot <-df[!duplicated(df$plotcode), -c(5:7, 12:27)] # Keep unique plot codes
  saveRDS(get("dfplot"), paste0(Dir, "/dfplot", CODE, seuil, ".rds")) # Work at the plot scale
  
  # Delete particular cases (100% recruitment by sp that result in NA for BA.ha.plot)
  # And put BAj to 0 instead of NA for plots whose mortality is 100%
  # Then recalculate the BA.O.plot.1
  
                              ############################
  ############################# Added on the 01/07/2018  ########################################
  dfplot <- dfplot[which(!is.na(dfplot$BA.ha.plot.1) | dfplot$sp.recruitment.plot.rate!=1),]   ##
  dfplot[which(is.na(dfplot$BAj.plot.1) & dfplot$sp.mortality.plot.rate==1),"BAj.plot.1"] <- 0 ##
  dfplot$BA.O.plot.1 <- dfplot$BA.ha.plot.1-dfplot$BAj.plot.1                                  ##
  ###############################################################################################
  
  # Mortality rate as a binomial (not really useful)
  dfplot$mort.bin <- dfplot[, "mortality.plot.rate"]
  dfplot[!is.na(dfplot$mortality.plot.rate) &
           dfplot$mortality.plot.rate != 0, 318] = 1
  dfplot[!is.na(dfplot$mortality.plot.rate) &
           dfplot$mortality.plot.rate == 0, 318] = 0
  dfplot[, 318] = as.numeric(dfplot[, 318])
  #  Sp Mortality rate as a binomial # Added on the 17/05/2018
  dfplot$sp.mort.bin <- dfplot[, "sp.mortality.plot.rate"]
  dfplot[!is.na(dfplot$sp.mortality.plot.rate) &
           dfplot$sp.mortality.plot.rate != 0, 319] = 1
  dfplot[!is.na(dfplot$sp.mortality.plot.rate) &
           dfplot$sp.mortality.plot.rate == 0, 319] = 0
  dfplot[, 319] = as.numeric(dfplot[, 319])
  
  
  # Divide by the number of years (perthousand unit)
  dfplot$country
  dfplot[dfplot$country == "FR", "yearsbetweensurveys"] <-
    5 # Add this line really matters
  dfplot$sp.mortality.plot.rate.yr <-
    round((
      dfplot$sp.mortality.plot.rate * 1000 / dfplot$yearsbetweensurveys
    ))
  dfplot$sp.mortality.plot.count.yr <-
    round((
      dfplot$sp.mortality.plot.count / dfplot$yearsbetweensurveys
    ))
  
  # Summary
  capture.output(print(summary(
    as.factor(dfplot$sp.mortality.plot.count)
  )),
  file = paste0(Dir, "dfplot.sp.mortality.plot.count.txt")) # Output as a latex wrapped in a txt file
  capture.output(print(summary(as.factor(
    dfplot$sp.mort.bin
  ))), file = paste0(Dir, "dfplot.sp.mort.bin.txt")) # Output as a latex wrapped in a txt file
  capture.output(print(summary(
    as.factor(dfplot$sp.mortality.plot.rate.yr)
  )),
  file = paste0(Dir, "dfplot.sp.mortality.plot.rate.yr.txt")) # Output as a latex wrapped in a txt file
  
  ## Added on the 22th june Add the minimum and mean SPEI
  Years <- c("spei01","spei03","spei06","spei12","spei18","spei24","spei36","spei48")
  for (i in Years){
    load(file = paste0("/home/achangenet/Documents/FUNDIV - NFI - Europe/our-data/climate/SPEI/plots_",i,"_1981_2015.RData"))
    dfplot[,paste0("min_",i)] <- spei_plots_all$min_spei_survey_years[match(dfplot$plotcode, spei_plots_all$plotcode,incomparables = NA)]
    dfplot[,paste0("mean_",i)] <- spei_plots_all$mean_spei_survey_years[match(dfplot$plotcode, spei_plots_all$plotcode,incomparables = NA)]
  }
         
  saveRDS(get("dfplot"), paste0(Dir, "/dfplot", CODE, seuil, ".rds")) # Work at the plot scale NO SCALE
  dfplot <-
    readRDS(paste0("dfplot", CODE, seuil, ".rds")) #Base de données plot full
  # Scale (based on all available data)
  
  ####### Added 05/07 # No scaled data but without NA and trasnfo ! ! ############
                                                                                ##
  dfplot <-dfplot[!is.na(dfplot$Plotcat) & dfplot$Plotcat != 10, ]              ## Unique plots that falls within the species distribution area remove NA and transition zone
  dfplot$Plotcat <- as.factor(dfplot$Plotcat)                                   ##
  dfplot <-dfplot[!is.na(dfplot$sp.mortality.plot.count.yr), ]                  ## Remove 2 NA or 14 (sp.count) 
  Transfo <- c("BA.ha.plot.1","BA.O.plot.1","BAj.plot.1",                       ##
               "dbh.plot.mean","BAIj.plot.1.mean","BAIj.plot.1")                ## Variable to transform 
  for (i in 1:length(Transfo)){                                                 ##
    dfplot[,paste0("log",Transfo[i])] <- log(dfplot[,Transfo[i]]+1)             ## Log transfo
    dfplot[,paste0("sqrt",Transfo[i])] <- sqrt(dfplot[,Transfo[i]])}            ## Sqrt transfo 
  saveRDS(get("dfplot"), paste0(Dir, "/dfplotbis", CODE, seuil, ".rds"))        ## Save this database as Dfplot normal 
                                                                                ##
  ################################################################################
  
  dfplot[, c(10:17, 38:43, 47:50, 60:317,338:(337+2*length(Transfo)))] <-scale(dfplot[, c(10:17, 38:43, 47:50, 60:317,338:(337+2*length(Transfo)))]
                                                                               , center = TRUE, scale =TRUE) #Center and reduced ALL
  dfplot2 <- dfplot
  #summary(as.factor(dfplot2[dfplot2$sp.mortality.plot.rate.yr==0,"Plotcat"])) #(763,1284,70)
  #summary(as.factor(dfplot2[dfplot2$sp.mortality.plot.rate.yr>0,"Plotcat"])) #(763,1284,70)
  capture.output(print(summary(as.factor(dfplot2[dfplot2$sp.mortality.plot.rate.yr ==
                                                   0, "Plotcat"]))),
                 file = paste0(Dir, "dfplot2.sp.mortality.rate.yr.txt")) # Output as a latex wrapped in a txt file
  capture.output(
    print(summary(as.factor(dfplot2[dfplot2$sp.mortality.plot.rate.yr > 0, "Plotcat"]))),
    file = paste0(Dir, "dfplot2.sp.mortality.rate.yr.txt"),
    append = T
  ) # Output as a latex wrapped in a txt file
  
  # At that point; there is not much plots left
  
  ###################################################################
  ###                                                             ###
  ###   Here is the dfplot2 (scaled plots but some are managed)   ###
  ###                                                             ###
  ###################################################################
  saveRDS(get("dfplot2"), paste0(Dir, "/dfplot2", CODE, seuil, ".rds"))
  dfplot2 <- readRDS(paste0("dfplot2", CODE, seuil, ".rds"))
  
  # Remove gestion(management 2 >= 1) & gest ???
  dfplot2 <-
    dfplot2[is.na(dfplot2$management2) |
              dfplot2$management2 == 0, ] # 255 NA (591,627 & 51)
  dfplot2 <-
    dfplot2[!is.na(dfplot2$management2) &
              dfplot2$management2 == 0, ] #  (581,539,51)
  summary(as.factor(dfplot2[dfplot2$sp.mortality.plot.count.yr > 0, "Plotcat"])) # Check my categories
  # If any, remove it
  dfplot2$gest <- as.factor(dfplot2$gest)
  dfplot3 <-
    dfplot2[-which(dfplot2$gest == 2 |
                     dfplot2$gest == 1), ] # Remove 459 + 537
  
  ##############################################################
  ###                                                        ###
  ###   Here is the dfplot3 (scaled plots and not managed)   ###
  ###                                                        ###
  ##############################################################
  capture.output(print(summary(as.factor(dfplot3[dfplot3$sp.mortality.plot.rate.yr ==
                                                   0, "Plotcat"]))),
                 file = paste0(Dir, "dfplot3.sp.mortality.rate.yr.txt")) # Output as a latex wrapped in a txt file
  capture.output(
    print(summary(as.factor(dfplot3[dfplot3$sp.mortality.plot.rate.yr > 0, "Plotcat"]))),
    file = paste0(Dir, "dfplot3.sp.mortality.rate.yr.txt"),
    append = T
  ) # Output as a latex wrapped in a txt file
  saveRDS(get("dfplot3"), paste0(Dir, "/dfplot3", CODE, seuil, ".rds"))
  dfplot3 <- readRDS(paste0("dfplot3", CODE, seuil, ".rds"))
}
