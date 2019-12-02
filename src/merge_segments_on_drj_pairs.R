# #############################################################################
# This file contains the definition of the functions mergeAllDrjs
#
# Author: Claire Lemaitre
# #############################################################################

##################################################
#   Find all bounds of segmentation segment      #
#        for each pair of DRJs                   #
#   And format the results for Gbrowse (Fabrice) #
##################################################

############Script###########
mergeAllDrjs=function(inputPairDir,outputPairDir){
	totalPaires=list.files(inputPairDir,full.names =TRUE)
        a=lapply(totalPaires, function(x) mergeSegmentBounds(x,outputPairDir))
}

mergeSegmentBounds=function(pairFile,outputPairDir){

  fileName=tail(strsplit(pairFile,"/")[[1]],1)
  tab=read.table(pairFile,header=TRUE)

  ## par paire de drj, plutot que par drj ??
  res=by(tab,paste(tab$inf1,tab$sup1,tab$inf2,tab$sup2,sep="-"), function(x) {
    n=nrow(x)
    #sc1=sum(x$FvsR==2)/n
    sc2=sum(x$segType==1 || x$segType==3)/n
    sc3=sum(x$segType==2 || x$segType==3)/n

    return(data.frame(inf1=x$inf1[1],sup1=x$sup1[1],inf2=x$inf2[1],sup2=x$sup2[1],nb=n,nbfrag=length(unique(x$frag)),segAlt=sc2,segDiff=sc3))
  })

  mat=matrix(unlist(res),ncol=8,byrow=T)
  final=as.data.frame(mat)
  names(final)=c("inf1","sup1","inf2","sup2","nb","nbfrag","segAlt","segDiff")

  write.table(final,paste(outputPairDir,"/",fileName,sep=""),quote=F,row.names=F)

  ## if(nrow(unique(drj1))!=nrow(unique(drj2)) || nrow(unique(drj1))!=nrow(unique(cbind(drj1,drj2)))){
  ##   return(tab)
  ## }
  ## print(paste(nrow(unique(drj1)),nrow(unique(drj2)),nrow(unique(cbind(drj1,drj2))),sep="--"))

	## tab$bounds1=paste(tab$inf1,",",tab$sup1,sep="")
	## tab$bounds2=paste(tab$inf2,",",tab$sup2,sep="")

	## # create tables with bounds
	## segmentationByDRJ(tab,tab$bounds1,1,Seg,pairFile)
	## segmentationByDRJ(tab,tab$bounds2,2,Seg,pairFile)
}
