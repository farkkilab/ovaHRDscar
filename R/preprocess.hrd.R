#' Preprocessing for further analysis
#'
#' @param seg segmentation data
#' @return preprocessed data
#This function will select the segments bigger in A tmp[,8] > tmp[,7] that do not spam all the chromosome !all(seg[,8] <= seg[,7]
#After reducing it get the number of events using shrink.seg.ai.wrapper
#Important to note, is that columns seg[,7] and seg[,8] are interchanged


preprocess.hrd<-function(seg){
  #Will ignore chromosomes X,Y
  seg <- seg[!seg[,2] %in% c(paste('chr',c('X','Y','x','y',23,24),sep=''),c('X','Y','x','y',23,24)),]
  seg[,1] <- as.character(seg[,1])
  seg[,2] <- as.numeric(gsub("chr","",seg[,2]))
  if (any(is.na(seg[,2]))){ stop("Chromosomes have unusual notation")}
  #Make sure that copy numbers values in A, are bigger than in B, if no, then invert columns
  if(! all(seg[,8] <= seg[,7]) ){
    tmp <- seg
    seg[tmp[,8] > tmp[,7],7]  <- tmp[tmp[,8] > tmp[,7],8]
    seg[tmp[,8] > tmp[,7],8]  <- tmp[tmp[,8] > tmp[,7],7]
  }
  seg <- shrink.seg.ai.wrapper(seg)
  return(seg)

}


#Preparing input
preparing.input <- function(seg){
  segAUX <- seg

  if (ncol(segAUX) != 7){
    stop("Input must contain only the next seven columns: SampleID, Chromosome, Start_position, End_position, total_cn, A_cn, B_cn")
  }

  segAUX <- cbind(segAUX, seg[,7])
  colnames(segAUX)[8] <- colnames(seg)[7]
  segAUX[,7] <- seg[,6]
  colnames(segAUX)[7] <- colnames(seg)[6]

  if (!(all(is.numeric(segAUX[,3])) & all(is.numeric(segAUX[,4])) & all(is.numeric(segAUX[,5])) & all(is.numeric(segAUX[,6])) & all(is.numeric(segAUX[,7])))){
	      stop("Some of the input rows do not contain numeric values in the columns: Start_position, End_position, total_cn, A_cn, B_cn")
  }

  preprocessed_seg <- preprocess.hrd(segAUX)
  #Next line unnecessary, was used to save time of analysis, but this can lead to error in samples with 0 allelic imbalances
  #preprocessed_seg <- rm.chr.normals(preprocessed_seg)
  length_breaks <- preprocessed_seg[,4] - preprocessed_seg[,3]
  preprocessed_seg <- preprocessed_seg[length_breaks > 50, ] #Removing variation less than 50bp
  return(preprocessed_seg)
}
