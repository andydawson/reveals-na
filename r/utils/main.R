
#########################################
##### REVEALSinR - the wrapper function
#########################################
# a <- REVEALSinR(pollenFile = "revealsR/reveals_test.csv",
#                 pf         = "revealsR/reveals_inputs.csv",
#                 dwm        = "lsm unstable",
#                 tBasin     = "lake",
#                 dBasin     = 2*round(lake_rad_site), # diameter!
#                 regionCutoff = 100000,
#                 repeats      = 1000)
REVEALSinR = function(pollenFile, pf, filetype='csv', 
                      fileoptions=if (filetype=='csv') {list(sep=",",dec=".")} else {list(sheetIndex=1)},
                      dwm, tBasin, dBasin, repeats=1e3, distMethod="WLS", regionCutoff=1e5,ppefun=rnorm_reveals,
                      pollenfun=rmultinom_reveals, writeresults=TRUE,verbose=TRUE){
  
  # read pollen data and parameters from the respective file
  if (filetype=='csv') 
  {
    pollen<-do.call(read.csv,args=c(file=pollenFile,fileoptions, row.names=NULL)) 
    parameters<-do.call(read.csv,args=c(file=pf,fileoptions, row.names=NULL))  
  }
  else if (filetype=='xls')   ##implicite import???
  {    
    pollen<-do.call(xlsx::read.xlsx,args=c(file=pollenFile,fileoptions))
    parameters<-do.call(xlsx::read.xlsx,args=c(file=pf,fileoptions))
  }
  else 
  {
    pollen <- pollenFile
    parameters <- pf
    
    pollenFile <- 'Pollen file'
    pf         <- 'Parameter file'
  }
  
  # check main parameters
  dwm<-tolower(dwm); tBasin<-tolower(tBasin) # all in lowercase to make the test easier
  if (!CheckInput(dwm=dwm, tBasin=tBasin, dBasin=dBasin, repeats=repeats, distMethod=distMethod, parameters=parameters)) 
    stop ("Some parameter is not suitable - check and correct")
  
  # extract PPEs and vg (fall speed)
  PPEs <- parameters[,c('species','PPEs','PPE.errors')]
  vg <- parameters[,c('fallspeed')]
  
  # calculate deposition factor depending on basin size/type and dispersal model
  if(verbose) { message("REVEALS:\tCalculating deposition factor for each species ...")}
  
  deposition <- do.call(rbind,lapply(vg, DispersalFactorK, tBasin=tBasin, dBasin=dBasin, dwm=dwm, regionCutoff=regionCutoff))
  
  # call REVEALS for each time slice in the pollen record
  all<-  apply(pollen, 1, REVEALS,  repeats=repeats, distMethod= distMethod, PPEs=PPEs, ppefun=ppefun,
               deposition=deposition,pollenfun=pollenfun,writeresults=writeresults,verbose=verbose)

  #making a data.frame for the results statistics
  results_df<-data.frame(pollenFile, pf, dwm, tBasin,pollen[1],
                       t(sapply(all,FUN=function(x)x$meansim)),
                       t(sapply(all,FUN=function(x)x$mediansim)),
                       t(sapply(all,FUN=function(x)x$q90sim)),
                       t(sapply(all,FUN=function(x)x$q10sim)),
                       t(sapply(all,FUN=function(x)x$sdsim)))

#with the right names
colnames(results_df)<-c("Pollen.file","Parameter.file", "Distance.weighting", "Basin.type",names(pollen)[1],paste(names(pollen[-1]),rep(c('meansim','mediansim','q90sim','q10sim','sdsim'),each=length(pollen[1,])-1),sep="."))

#writing them if necessary
if(writeresults)   {
#   # add statistics to overview file - File l?schen zu beginn - nur wie? --not necessary anymore. 
#   if(!(file.exists("ResultStatistics_confidence_limits.csv"))) {
#     write.table(data.frame("Parameter file", "Distance weighting", "Basin type","median", t(names(pollen[-1])),"", "10% Quantile", t(names(pollen[-1])),"","90% Quantile", t(names(pollen[-1]))), file="ResultStatistics_confidence_limits.csv", append = TRUE, row.names=FALSE, col.names = FALSE, sep=";", dec=",")
#   }
#   
#   if(!(file.exists("ResultStatistics_SD.csv"))) {
#     write.table(data.frame("Parameter file", "Distance weighting", "Basin type","mean", t(names(pollen[-1])),"", "SD", t(names(pollen[-1]))), file="ResultStatistics_SD.csv", append = TRUE, row.names=FALSE, col.names = FALSE, sep=";", dec=",")
#   }
  nn<-names(results_df)
  #quantiles
  k<-c(1:4,grep(x=nn,pattern='mediansim'),grep(x=nn,pattern='q10sim'),grep(x=nn,pattern='q90sim'))
  write.table(results_df[,k],file="ResultStatistics_confidence_limits.csv", append = FALSE, row.names=FALSE, col.names = TRUE, sep=";", dec=",")
 
  #mean+sd
  k<-c(1:4,grep(x=nn,pattern='meansim'),grep(x=nn,pattern='sdsim'))
  write.table(results_df[,k],file="ResultStatistics_SD.csv",  append = FALSE, row.names=FALSE, col.names = TRUE,  sep=";", dec=",")
}


invisible(results_df)
}

#############################################
##### REVEALS - applying the REVEALS model to each time slice
#############################################
# 
# # call REVEALS for each time slice in the pollen record
# all<-  apply(pollen, 1, REVEALS,  repeats=repeats, distMethod= distMethod, PPEs=PPEs, ppefun=ppefun,
#              deposition=deposition,pollenfun=pollenfun,writeresults=writeresults,verbose=verbose)

REVEALS = function(counts, repeats=1e3, PPEs, ppefun=rnorm_reveals,deposition,pollenfun=rmultinom_reveals, 
                   distMethod="WLS",writeresults=TRUE,verbose=TRUE){
  
  # print time slice to indicate progress
if (verbose)  {message(paste('REVEALS:\tcomputing for age slice\t',counts[1]))}
  
  # start calculation round 
  coverM<-replicate(n=repeats,{
    
    # add random error to PPEs (using mean and stanadard deviation), 'while loop' is to avoid negative values
    PPEs[,4]<-ppefun(length(PPEs[,2]), PPEs[,2], PPEs[,3])
    
    
    # add random error to pollen data (by default: random drawing with rnorm from the given composition)
    pollenPer <- pollenfun(n=1,pollen=counts[-1])
    
    # prepare vector dispersal x PPE
    disp_ppe<-deposition*PPEs[,4]
    
    # calculate cover (cf. REVEALS model formula B1 in Sugita 2007 REVEALS paper)
    s <- sum(pollenPer/disp_ppe)
    round(100*pollenPer/(disp_ppe * s),digits=3)
    
  },simplify=TRUE)
  
   # convert results into a matrix
  rownames(coverM) <- names(counts[-1])
  
  # write complete data file for each time slice
  if(writeresults) {
    fileName <- paste("results_cover", as.character(pollen[1]),".csv",sep="")
    write.csv2(t(coverM), fileName)
  }
  # calculate mean
  coverMean <- apply(coverM, 1, mean)
  
  # calculate median
  coverMedian <- apply(coverM, 1, median)
  
  # calculate standard deviation
  coverSD <- apply(coverM, 1, sd)
  
  # quantiles
  coverq10 <- apply(coverM,1,function(u) quantile(probs=.1,x=u))
  coverq90 <- apply(coverM,1,function(u) quantile(probs=.9,x=u))
  
  return(invisible(list(meansim=coverMean,sdsim=coverSD,mediansim=coverMedian,q90sim=coverq90,q10sim=coverq10)))
}


##################################################################################
##### DispersalFaktorK
##################################################################################

# calculates the dispersal factor K for each species (=influx if the cover of that taxon were 100%)
DispersalFactorK=function(vg, tBasin, dBasin, dwm, regionCutoff){ 
  
  # print fall speed to indicate progress
  #print(vg)
  rBasin <- dBasin/2
  
  if (dBasin>regionCutoff) stop("Size of the region too small!")
  
  # branch with dispersal model (for parameters of GPM model see Jackson & Lyford 1999)
  if (dwm=="gpm neutral") {
    disModel<-GPM(fallSpeed=vg, cz=0.12, n=0.25, u=3)
  } 
  else if (dwm=="gpm unstable") {
    disModel<-GPM(fallSpeed=vg, cz=0.21, n=0.20, u=3)
  } 
  else if (dwm=="1overd") {
    disModel<-OneOverD.f()
  } 
  else if (dwm=="lsm unstable")   {
    vgr<-100*round(vg,2)+1 #find right column
    disModel<-lsm_unstable[[vgr]]
    # with the following option, alternative models may be used... 
    #} else if (dwm=="alternative")   {
    #  load("alternative.rda") # load model, which is prepared using a loess smoother
    #  vgr<-100*round(vg,2)+1 #find right column
    #  disModel<-alternative[[vgr]]
  } 
  else stop("no valid distance weighting method selected")
  
  
  # for lakes, add additinal deposition arriving due to lake mxing 
  if (tBasin=="peatland") {
    if (regionCutoff<= 1000000) influx <- predict(disModel,rBasin)-predict(disModel,regionCutoff) # deposition in a peatland
    else influx <- predict(disModel,rBasin) # if region > 1000 km - ignore limits of the region
  } 
  else if (tBasin=="lake") {
    
    # set ring width for lake model, 2 m in smaller basins and 10 m in larger basins
    stepSize=2 
    if (dBasin>=2000) stepSize=10 
    
    # create vector with steps along lake radius from center to margin
    rSteps <- seq(from = stepSize, to = rBasin-1, by = stepSize) 
    
    # call the lake model 
    lakeInflux <- do.call(rbind,lapply(rSteps, LakeModel, dBasin=dBasin, disModel=disModel, stepSize=stepSize, regionCutoff=regionCutoff))
    
    #devide total influx by (relevant) lake area 
    influx <- sum(lakeInflux) / (pi*(max(rSteps))^2) 
  }
  
  return(influx)    
}


##################################################
##### LakeModel
##################################################

# calculates additonal pollen deposition on lake points between the center and margin of the lake
# x is the distance between lake center and point for which deposition is calculated

LakeModel = function(x, dBasin, disModel, stepSize, regionCutoff){
  rBasin<-dBasin/2
  
  # create a sequence of radii from the smallest to the largest needed
  r2 <- seq(from = rBasin - x, to = 2*rBasin, by = 2) 
  
  # cut off high values (>rBasin+x), otherwise alpha calculation fails
  r2Cut <- r2; r2Cut[r2Cut>(rBasin+x)]=rBasin+x 
  
  # calculate the angel alpha (see Appendix)
  alpha <- acos((r2Cut^2+x^2-rBasin^2)/(2*r2Cut*x)) 
  
  # the proportion of the circle outside the lake area
  propOut <- 1-alpha/pi                             
  a <- diff(propOut)/2 ; propOut <- propOut[-1]-a   # produce mean values between two circles
  
  # pollen airborne at each circle 
  airborne <- predict(disModel, r2)                 
  
  # calculate difference -> influx from full rings
  influxRing <- abs(diff(airborne))                 
  
  # multiply with proportion of each ring that is outside the lake --> true influx from each ring
  influxRing <- influxRing*propOut    
  
  # influx is sum of all rings plus deposition from outside
  influx <- sum(influxRing) + predict(disModel,2*rBasin) - predict(disModel,regionCutoff) 
  
  ringArea <- pi*x^2-pi*(x-stepSize)^2
  
  # return total influx multiplied by ringArea to account for the increasing area with increasing x
  return(influx*ringArea)  
}


###################################################
##### Check parameters 
###################################################

CheckInput <- function(dwm, tBasin, dBasin, repeats, distMethod, parameters){
  
  # check distance weighting method
  if (!isTRUE((dwm=="lsm unstable") | (dwm=="gpm neutral") | (dwm=="gpm unstable") | (dwm=="1overd"))) stop("distance weighting method not defined; should be 'LSM unstable', 'GPM neutral', 'GPM unstable' or '1oved' ")
  
  # check basin type
  if (!isTRUE((tBasin=="peatland") | (tBasin=="lake"))) stop("basin type ('tBasin') not defined; should be 'peatland' or 'lake'")
  
  # check basin diameter
  if (!isTRUE(dBasin == floor(dBasin))) stop("basin size ('dBasin') should be an integer value")
  if (dBasin < 10) stop("basin diameter ('dBasin') should be at least 10 m")
  
  # check number if iterations (should be larger than 1000)
  if (!isTRUE(repeats == floor(repeats))) stop("'repeats' should be an integer value")
  if (repeats < 1000) stop("'repeats' should be at least 1000")
  
  # check the distance measure for calcualting the difference between modelled and empiric pollen data
  distMethod <- toupper(distMethod)
  if (!isTRUE((distMethod=="WLS") | (distMethod=="OLS"))) stop("distance measure not defined; should be WLS or OLS")
  
  # species column present?
  if(!'species' %in% tolower(names(parameters))) {
    print ("no 'species' column in 'parameters.csv'")
    return(FALSE)
  }
  
  #fallspeed column present?
  if(!'fallspeed' %in% tolower(names(parameters))) {
    print ("no 'fallspeed' column in 'parameters.csv'")
    return(FALSE)
  }
  
  # values between 0.01 and 0.15?
  if(max(parameters$fallspeed)>0.15) stop ("fallspeed(s) too high (>0.15 m/s), please check")
  if(max(parameters$fallspeed)<0.01) stop ("fallspeed(s) too low (<0.01 m/s), please check")
  
  # PPEs column present?
  if(!'ppes' %in% tolower(names(parameters))) {
    print ("no 'PPEs' column in 'parameters.csv'")
    return(FALSE)
  }
  
  # PPE.errors column present?
  if(!'ppe.errors' %in% tolower(names(parameters))) {
    print ("no 'PPE.errors' column in 'parameters.csv'")
    return(FALSE)
  }
  
  
  return(TRUE)
}

#########################################
##### add error in pollen data 
#########################################

# default: add errors to pollen data through random drawing
rmultinom_reveals<-function(n=1,pollen,...)
{
  100*rmultinom(n, sum(pollen), pollen/sum(pollen))/sum(pollen)
}

# alternative: add error to pollen data following Maher (1972)
rMaher <- function(n=1,pollen,...){
  
  countSum<-sum(pollen[-1])
  pDat <- pollen[-1]/countSum 
  
  ranVec<-rnorm(length(pollen)-1, mean=0, sd=1) #random Vector with mean of zero and 1 standard deviation
  i1 <- pDat+((ranVec^2)/(2*countSum))
  i2 <- (pDat*(1-pDat)/countSum)+((ranVec^2)/(4*countSum^2))
  denom <- 1+((ranVec^2)/countSum)
  ranPol <- 100*(i1+(ranVec*sqrt(i2)))/denom
  return(ranPol)
}

#changing the PPE's randomly

rnorm_reveals<- function(n,mean=rep(1,n),sds=rep(1,n)) 
{x<-0*mean-1
 #iterate till positivity is established
 while(any(x<0)){x<-rnorm(n=n,mean=mean,sd=sds)}    
 return(x)
}

##########################################################################
##### Gaussian plume model (GPM) - Sutton model follwing Prentice (1985) 
##########################################################################

# prepares a model with a Gaussian plume model
# Arguments:
# n: turbulence parameter
# cz: vertical diffusion coefficient
# u: mean wind speed [meter/s]

GPM = function (fallSpeed, cz, n, u){
  x <- seq(from = 0, to = 100, by = 0.1)
  x <- x^3
  y <- exp((-4*fallSpeed*x^(n/2))/(n*u*sqrt(pi)*cz))
  return(loess(y ~ x, span=0.015, model=FALSE, degree=2))
}


###################################################
##### 1overd dispersal model 
###################################################

# prepares a model with the 1/d dispersal function
# bitte ignorieren - noch nicht fertig
OneOverD.f = function (){
  x <- seq(from = 1, to = 100, by = 0.15) # klappt nicht sauber
  x <- x^3
  y <- 1/x
  # es muss ein remaining airborne vector bei raus kommen
  y<-y/sum(y)
  # plot(y~x, xlim=c(0,100))
  loess(y ~ x, span=0.015)
  return(loess(y ~ x, span=0.015, model=FALSE, degree=2))
}


#################################################################
#### for the future: Plot functions
#################################################################


# plot_compts<-function(compdata,cols=rainbow(ncol(compdata)-1))
# {if(!("Age"%in% names(compdata))) stop("no historical percentage data")
# 
# mdata <- melt(compdata, id=c("Age"))  #mydate ist der output
# 
# plot1 <- ggplot(mdata, aes(x = mdata$Age, y = mdata$value, fill = mdata$variable)) + 
#   geom_area(position = 'stack') + 
#   scale_fill_manual(values=cols)
# }
# 
# plot_compslice<-function(pollen)
#   plot(acomp(pollen))

