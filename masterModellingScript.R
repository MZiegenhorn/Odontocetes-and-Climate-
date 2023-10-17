
#This script used to produce all models for Ziegenhorn et al. 2023 "Odontocete detections are linked to oceanographic conditions in the Hawaiian Archipelago"


rm(list=ls())
graphics.off()


library(ggplot2)
library(plyr)
library(dplyr)
library(nlme)
library(parsedate)
library(geepack)
library(splines)
library(pracma)  
library(mctest)
library(lubridate)
library(gridBase)
library(grid)
library(tidyverse)  
library(ROCR)           
library(PresenceAbsence) 
library(mvtnorm)       
library(gridExtra)       
library(SimDesign)
library(regclass)
library(devtools)
library(car)
library(splines2)
library(stats)
library(mgcv) 

################ SETTINGS ####################################
site = "Manawai" 
kval = 4 #knots value for all models
familyuse = 'nb' #what family to use for modelling
infolder = file.path('/Users/maziegenhorn/Documents/SIO/thesis_chapters/ch3/manuscript/final/code_and_data/')
outfolder = file.path('/Users/maziegenhorn/Documents/SIO/thesis_chapters/ch3/modTestforPub')
acval = 0.05 #autocorrelation value
class = c("1","2","3","4","5_6","7_8","9","10")
concurval = 0.6 #concurvity threshold

###################################################################
infilename = paste(infolder,site,"_counts_finalData.csv",sep='')


####make directory if needed
if (!dir.exists(outfolder)){
  dir.create(outfolder)
}

#extract name
namesh = basename(infilename)
outnamesh = str_replace(namesh,'counts_finalData.csv','')

#run modelling for each class
for (ifi in seq_along(class)){
  
  #setup our file 
  Daily <- read.table(file = infilename,sep = ',',header = TRUE)
  
  #get focal class
  useclass = class[ifi]
  findclass = which(names(Daily)==paste0("class",useclass))
  Daily$pres = round(Daily[,findclass]) #rounding it because poisson wants counts, which should always be integers
  #then get rid of the other class variable for that class
  Daily = Daily[,-findclass]
  Daily = Daily %>% relocate(pres,.after=day)
  
  #setup outname
  outnamesh2 = paste0(outnamesh,'class',useclass,'_',familyuse,'.RData')
  outname = file.path(outfolder,outnamesh2)
  
  Dayuse = Daily 
  
  #save a variable counting the number of days with presence
  nPresDays = size(which(Dayuse$pres>0),2)
  print(nPresDays)
  
  #only run everything if there's more than 100 presence points
  if (nPresDays>=100){
    
    #sort things by time
    Dayuse = Dayuse[order(Dayuse$day),]
    
    varnames = c("enso","pdo","npgo","ssht_dep0","salt_dep0","tempt_dep0")
    
    if (isempty(varnames)){
      print(paste0('no variables provided for ',outnamesh2,'. Skipping this model!'))
    }else{
      
      ###### basic model ###################
      smvars = varnames
      smooths = paste(paste0("s(",smvars,",k=4,bs='cr')"),collapse = "+")
      myformula = as.formula(paste("pres~",smooths))
      basicMod<-gam(myformula,family=familyuse,data=Dayuse)
      summary(basicMod)
      
      ######################## CORRELATIONS  ##########################################################
      
      #look for concurvity if there's more than one variable
      if (length(varnames)>1){
        
        #also look at VIF from full model
        varnamesvif = varnames
        smvars = varnamesvif
        smooths = paste(paste0("s(",smvars,",k=4,bs='cr')"),collapse = "+")
        myformula = as.formula(paste("pres~",smooths))
        
        fullMod<-gam(myformula,family=familyuse, data=Dayuse)
        AICfull = AIC(fullMod)
        
        concurall = concurvity(fullMod)
        #get estimate values- remove any greater than 0.5
        badcc = concurall[2,concurall[2,]>concurval]
        
        #if any gvif^2 values are higher than 3, remove largest and try again
        while (!isempty(badcc)){
          ccorder = concurall[2,order(concurall[2,],decreasing=TRUE)]
          rmvvartemp = str_replace(names(ccorder)[1],'s\\(','')
          rmvvar = str_replace(rmvvartemp,'\\)','')
          
          varnamesvif = varnamesvif[-c(grep(rmvvar,varnamesvif))]
          smvars = varnamesvif [grep("enso|pdo|ssh|temp|sal",varnamesvif)]
          smooths = paste(paste("s(",smvars,",k=4,bs='cr')"),collapse = "+")
          # if (grep("climatestate",varnamesvif)){
          #   myformula = as.formula(paste("pres~",paste(c(smooths,"as.factor(climatestate)"),collapse = "+")))
          # }else{
          myformula = as.formula(paste("pres~",smooths))
          # }
          fullMod<-gam(myformula,family=familyuse, data=Dayuse)
          concurall = concurvity(fullMod)
          badcc = concurall[2,concurall[2,]>concurval]
        }
        
      }else{
        varnamesvif = varnamestype
        vartypevif = vartype
        smvars = varnamesvif [grep("enso|pdo|ssh|temp|sal",varnamesvif)]
        smooths = paste(paste0("s(",smvars,",k=4,bs='cr')"),collapse = "+")
        myformula = as.formula(paste("pres~",smooths))
        
        fullMod<-gam(myformula,family=familyuse, data=Dayuse)
        #full model
        AICfull = AIC(fullMod)
        
        print("Only one variable! No VIF analysis required.")
      }
      
      #get variable names again, excluding those that have been removed by concurvity 
      varsort = smvars
      
      ########## RETEST AUTOCORRELATION  ##################
      print("formula after concurvity analysis: ")
      print(myformula)
      
      corrMod<-gam(myformula,family=familyuse,data=Dayuse)
      #anova to check for variable significance, summary to look at other stuff
      summary(corrMod)
      acfRes = acf(residuals(corrMod), lag.max = 1000, ylim=c(0,0.5))
      acf(residuals(corrMod),lag.max = 100, ylim=c(0,0.2))
      
      #find acf value below 0.1 to use
      acuse = Position(function(x) x < acval,acfRes$acf)
      if (is.na(acuse)){#if it never drops below this value
        acvaln = acval*2
        print(paste0("Warning: autocorrelation does not drop below ",acval,"New limit = ",acvaln))
        acuse = Position(function(x) x < acvaln,acfRes$acf)
      }
      #one more in case that doesn't work either
      if (is.na(acuse)){#if it never drops below this value
        acvaln2 = acvaln*2
        print(paste0("Warning: autocorrelation does not drop below ",acvaln,"New limit = ",acvaln2))
        acuse = Position(function(x) x < acvaln2,acfRes$acf)
      }
      print(acuse)
      
      #ex: autocorrelation drops below 0.1 after 2 hours, so clusters for id should be 2 observations long
      #this is used to set id in geeglm(); each row in id is of observations that might be correlated,
      #observations in different rows are assumed to be not correlated
      Dayuse$id = 0
      deps = data.frame(unique(Dayuse$dep))
      #set up blocks at deployment level so that no block spans multiple deployments 
      startval = 1
      for (idep in 1:nrow(deps)){
        nrowdep = c(grep(deps$unique.Dayuse.dep[idep],Dayuse$dep))
        endval = (startval-1)+size(nrowdep,2)/acuse
        Dayuse$id[nrowdep] = round(linspace(startval,endval,size(nrowdep,2)))
        startval = endval+1
      }
      
      #use ID to create new table
      Daymeans = list()
      varnames = c("pres","day",varnames)
      for (iv in seq_along(varnames)){
        valtemp = aggregate(Dayuse[,colnames(Dayuse)==varnames[iv]],list(Dayuse$id),FUN=mean)
        Daymeans[[varnames[iv]]] = valtemp$x
      }
      
      #convert daymeans to a data frame
      Daymeans = as.data.frame(Daymeans)
      
      #redo full model with autocorrelation
      fullMod<-gam(myformula,family=familyuse, data=Daymeans)
      AICfull = AIC(fullMod)
      
      ########################################  FINAL MODEL  ############################################################
      library(MuMIn)
      smooths = paste(paste0("s(",smvars,",k=4,bs='cr')"),collapse = "+")
      myformulad = as.formula(paste("pres~",smooths))
      
      m1d<-gam(myformulad,family=familyuse, data=Daymeans,na.action="na.fail")
      testlist = dredge(m1d,rank="AICc")
      finalmod = testlist[1,]
      summary(finalmod)
      
      #get final mod
      n = names(finalmod)
      dfloc = which(n=="df")
      n=n[3:dfloc-1]
      listtr = finalmod[,3:dfloc-1]
      finalvars = n[!is.na(listtr)]
      
      if (!isempty(finalvars)){
        formfinal = as.formula(paste("pres~",paste(finalvars,collapse = "+")))
        m1final<-gam(formfinal,family=familyuse, data=Daymeans,na.action="na.fail")
        m1finsum = summary(m1final)
      }else{
        print("best model is null model!")
        formfinal = as.formula("pres~1")
        m1final<-gam(formfinal,family=familyuse, data=Daymeans,na.action="na.fail")
        m1finsum = list()
        fullMod = list()
      }
      
      
      # ###### additional stuff since using full model
      m1finsum = summary(m1final)
      AICfinal = AIC(m1final)
      print(m1finsum)
      
      #save model
      outname = str_replace(outname,'_finalData.RData','_ch3mod.RData')
      save(list = c("fullMod","m1final","m1finsum","Daily",
                    "Dayuse","Daymeans","nPresDays",
                    "AICfull","AICfinal","varnames"),file=outname)
      
      
      #############  PLOTTING   ############################################################################################
      ### plot using ggeffects package!
      library(ggeffects)
      #get names in usable form
      allnamesfinal = vector()
      
      for (is in seq_along(varsort)){
        ind = grep(varsort[is],finalvars)
        if (!isempty(ind)){
          allnamesfinal = c(allnamesfinal,varsort[is])
        }
      }
      
      for (iv2 in seq_along(allnamesfinal)){
        varname = allnamesfinal[iv2]
        vals = Dayuse[, which(names(Dayuse) == varname)]
        dfvar = data.frame(vals)
        
        #setup color for plots
        if(grepl("sal",varname)){
          col1 = 'sienna3'
          col2 = 'plum2'
        }else if (grepl("temp",varname)){
          col1 = 'darkgoldenrod3'
          col2 = 'khaki'
        }else if (grepl("npgo|pdo|enso",varname)){
          col1 = 'royalblue4'
          col2 = 'lightskyblue'
        }else if (grepl("ssh",varname)){
          col1 = 'darkred'
          col2 = 'firebrick1'
        }
        
        ggplot(ggpredict(m1final,terms=varname,interval="confidence"),aes(x=x,y=predicted))+
          geom_line(color="black")+
          geom_ribbon(aes(ymin=conf.low,ymax=conf.high),fill=col1,alpha=0.15)+
          theme_bw()+
          theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())+
          geom_rug(aes(x=x),sides="b",inherit.aes=FALSE)+
          xlab(varname)+
          theme(legend.position = "none")+
          ggtitle("")
        
        plotname2 = str_replace(outname,'.RData',varname)
        plotname2 = paste(plotname2,'.tiff',sep='',collapse='')
        ggsave(plotname2,device='tiff',dpi=700,compression="lzw")
        
      }
    }
  }else{
    print(paste("nPresDays only ",nPresDays,". Skipping ",outnamesh2,sep=""))
  }}


