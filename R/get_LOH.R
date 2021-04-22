#' Extract features per samples
#'
#' @param seg segmentation data

######### Function to extract the number of LOH events of a particular size ########
#Start allways MbSizes with 0
#As output will give default output will give you the LOH events of sizes 0 to 1MB, 1Mb to 4Mb, 7Mb to 10Mb, and so on... according to the windows sizes selected
features.LOH <- function(segs,  MbSizes = c(0,1,4,7,10,13,16,19,22,25,28,31,34,37,40), A_cn=7, B_cn=8){
    #For each of  the samples
    totalCountLOHs <- NULL
    for(sample in unique(segs[,1])){
        sample.seg <- segs[segs[,1] %in% sample,]
        #Selection of LOH events
        segLOH <- sample.seg[sample.seg[,B_cn] == 0 & sample.seg[,A_cn] != 0,,drop=F]
        sizesLOH <- (segLOH[,4] - segLOH[,3])/1e6
        countLoHSample <- NULL
        for (i in 2:length(MbSizes)){
              countLOHs <- length(which(sizesLOH >= MbSizes[i-1] & sizesLOH <= MbSizes[i]))
              countLoHSample <- c(countLoHSample, countLOHs)
        }
        countLOHs <- length(which(sizesLOH > MbSizes[i]))
        countLoHSample <- c(countLoHSample, countLOHs)
        totalCountLOHs <- rbind(totalCountLOHs, countLoHSample)
    }
    row.names(totalCountLOHs) <- unique(segs[,1])
    colnames(totalCountLOHs) <- c(paste("LoH", paste(MbSizes, "Mb", sep="" ), sep="_"))
    #colnames(totalCountLOHs) <- c(paste(MbSizes, "Mb", sep="" ))
    totalCountLOHs <- as.data.frame(totalCountLOHs)
    return(totalCountLOHs)
}