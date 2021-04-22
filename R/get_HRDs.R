#' New HRD scar score

get.ovaHRDscars <- function(seg, chrominfo = chrominfo_grch38, LOH_windos=c(10,50), LST_segSize=12e6, LST_mindistance=1e6){
  seg <- preparing.input(seg)

  #Calculating HRD-LOH
  HRD_LOHs <- features.LOH(seg, MbSizes=LOH_windos)
  HRD_LOHs <- HRD_LOHs[,1]

  #Calculating LSTs
  LSTs <- LSTs(seg, chrominfo=chrominfo, segsizes=LST_segSize, mindistance=LST_mindistance)

  #Calculating nTAIs
  res_ai<- calc.ai_new(seg, chrominfo)
  nTAIs <- res_ai[,1]

  ovaHRDscar <- HRD_LOHs + LSTs + nTAIs

  #Concatenating results

  HRDresulst <- cbind(HRD_LOHs, LSTs, nTAIs, ovaHRDscar)
  colnames(HRDresulst) <- cbind("nLOH","LSTs","nTAIs", "ovaHRDscar")
  assign("HRDresulst",as.data.frame(HRDresulst),envir = .GlobalEnv)
  return(HRDresulst)
}