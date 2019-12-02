# #######################################################################################
# This file contains the definition of the function drawAllSummary for scaffolds analyses
#
# Authors: Claire Lemaitre and Stéphanie Robin
# ########################################################################################

###########################################################
#   For all pairs of DRJs draw segmentation on the DRJs   #
#         name = TRUE : display the id of the read        #
###########################################################

##########Main function#####################
drawAllSummary=function(pairDir,outputDir,readNames=FALSE,addSimi=NULL,merge=FALSE){
  ## addSimi : si directory avec mismatch vector pour chaque pair de DRJs (seulement sur drj1)

  total_files=list.files(pairDir,full.names =TRUE)
  draw=lapply(total_files, function(x) drawByFile(x,outputDir,readNames,addSimi,merge))
}

## TODO : addSimi
## Stéphanie : modif des coordonnes du tableau pour que ca corresponde aux scaffolds
drawByFile=function(file,outputDir,readNames=FALSE,addSimi=NULL,merge=FALSE){

  fileName=strsplit(tail(strsplit(file,"/")[[1]],1),".tab")[[1]][1]
  coord.DRJs=strsplit(fileName,"_")[[1]]
  bac=coord.DRJs[1]
  infDRJ1=as.numeric(coord.DRJs[4])
  supDRJ1=as.numeric(coord.DRJs[5])
  infDRJ2=as.numeric(coord.DRJs[6])
  supDRJ2=as.numeric(coord.DRJs[7])

  allSeg=read.table(file,header=TRUE)
  allSeg=allSeg[order(allSeg$inf1,allSeg$sup1,allSeg$read),]
  nbReads=nrow(allSeg)

  nbMer=0
  yMer=NULL
  if(merge){
    #merged=read.table(paste(mergeDir,"/",fileName,".tab",sep=""),header=TRUE)
    merged=unique(allSeg[,c("inf1","sup1","inf2","sup2")])
    merged=merged[order(merged$inf1,merged$sup1),]
    nbMer=nrow(merged)
  }


  index.doublon=which(allSeg$doublon==TRUE)

  xtitle=paste("Position on bac",bac)
  ytitle="Reads"

  pdf(paste(outputDir,"/",fileName,".pdf",sep=""),width = 11, height = 6, onefile = FALSE)

  par(mfrow=c(1,2), mar=c(5,4,1.5,1), oma=c(0,0,2,0))
  y=1:nbReads
  couleur=ifelse(allSeg$multi,"orange","red")
  couleur=ifelse(allSeg$segType>0,"grey40",couleur)

  ## DRJ1
  plot(c(infDRJ1,supDRJ1),c(-nbMer,nbReads),type="n",xlab=xtitle,ylab=ytitle,main="DRJ1")

  if(!is.null(addSimi)){
    mismatch=read.table(paste(addSimi,"/",fileName,".tab",sep=""))[,1]
    pcId=100-sum(mismatch)/(sum(mismatch+1))*100
    ## plotter les différences
    vmismatch=infDRJ1+which(mismatch!=0)
    abline(v=vmismatch,col="grey80")
    text(supDRJ1,0,paste(signif(pcId,3)," %id (drj)",sep=""),adj=c(1,0),col="grey30")
    ## fenetre glissante
    xpc=NULL
    ypc=NULL
    l=19
    for (i in 1:(length(mismatch)-l)){
      xpc=c(xpc,i+floor(l/2))
      ypc=c(ypc,100-sum(mismatch[i:(i+l)]>0)/(l+1)*100)
      }
    ypcNorm=ypc*(nbReads)/100
    lab=seq(0,100,by=10)
    at=lab*(nbReads)/100
    lines(xpc+infDRJ1,ypcNorm,col="lightsteelblue3")
    axis(4,at=at,labels=lab,col="lightsteelblue3",col.axis="lightsteelblue3",cex.axis=0.8)
    mtext("%id", col="lightsteelblue3", adj=1)
  }

  if(merge){
    yMer=-nbMer:(-1)
    segments(merged$inf1,yMer,merged$sup1,yMer, col="blue")
    abline(h=0,lty=2,col="grey")
  }
  segments(allSeg$inf1,y,allSeg$sup1,y, col=couleur)

  points(allSeg$sup1[index.doublon],index.doublon,pch="*")
  if(readNames==TRUE){
    text(allSeg$inf1, y,labels=as.character(allSeg$id),cex=0.5)
  }


  ## DRJ2
  plot(c(infDRJ2,supDRJ2),c(-nbMer,nbReads),type="n",xlab=xtitle,ylab=ytitle,main="DRJ2")
  if(merge){
    yMer=-nbMer:(-1)
    segments(merged$inf2,yMer,merged$sup2,yMer, col="blue")
    abline(h=0,lty=2,col="grey")
  }
  segments(allSeg$inf2,y,allSeg$sup2,y, col=couleur)
  points(allSeg$sup2[index.doublon],index.doublon,pch="*")
  if(readNames==TRUE){
    text(allSeg$inf2, y,labels=as.character(allSeg$id),cex=0.5)
  }

  if(addSimi==TRUE){
    ## if(length(seq(infDRJ2,supDRJ2))==length(simDRJ2[seq(32,length(simDRJ2)-30)])){
    ##   lines(seq(infDRJ2,supDRJ2),simDRJ2[seq(31,length(simDRJ2)-30)]*nbTot,
    ##         type="l",col="green")
    ## }
  }


  mtext(paste(fileName," : ",length(unique(allSeg$frag))," fragments différents",sep=""), side=3, line=0, adj=0.5, cex=1.4, outer=TRUE)

  dev.off()

}



compareSummary=function(refDir,outputDir=NULL,outputFile=NULL,additionalDir=NULL,dirNames=NULL,addSimi=NULL){
  ## pour additionalDir : permet de comparer sur un meme plot plusieurs ensembles de segmentations, on fait l'hypothese que ce sont exactement les memes coordonnees de drj
  ## output : soit outputDir : 1 fichier pdf par paire de drj
  ##          soit outputFile : 1 fichier pdf pour toutes les drjs

  if(!is.null(outputFile)){
    pdf(outputFile,width = 6, height = 8)
  }

  total_files=list.files(refDir,full.names =TRUE)
  draw=lapply(total_files, function(x) drawAndCompareByFile(x,outputDir,additionalDir,dirNames,addSimi))

  if(!is.null(outputFile)){
    dev.off()
  }

  ## formate les stats de nb de reads par paire de DRJ, sous la forme d'un tableau
  nbAdd=length(additionalDir)
  tmp=matrix(unlist(draw),ncol=3+nbAdd*2,byrow=T)
  tmp2=as.data.frame(matrix(as.numeric(tmp[,-1]),ncol=2+nbAdd*2,byrow=F))
  if(nbAdd>0){
    name=c("ref",paste("add",1:nbAdd,sep=""))
    names(tmp2)=paste(rep(name,rep(2,1+nbAdd)),c(".tot",".good"),sep="")
  }
  else{
    names(tmp2)=c("tot","good")
  }

  tmp2$name=tmp[,1]
  n=ncol(tmp2)
  stats=tmp2[,c(n,1:(n-1))]

  return(stats)

}


## Plot d'une figure de summary sur laquelle on peut superposer plusieurs jeux de données
## il faut que les fichiers de jonction aient le même nom (bac_inf1... .tab) mais dans des répertoires différents
## renvoie également un vecteur avec les nombre de read
drawAndCompareByFile=function(file,outputDir=NULL,addDir=NULL,dirNames=NULL,addSimi=NULL){

  fileName=strsplit(tail(strsplit(file,"/")[[1]],1),".tab")[[1]][1]
  coord.DRJs=strsplit(fileName,"_")[[1]]
  bac=coord.DRJs[1]
  infDRJ1=as.numeric(coord.DRJs[2])
  supDRJ1=as.numeric(coord.DRJs[3])
  infDRJ2=as.numeric(coord.DRJs[4])
  supDRJ2=as.numeric(coord.DRJs[5])

  allSeg=read.table(file,header=TRUE)
  allSeg=allSeg[order(allSeg$inf1,allSeg$sup1,allSeg$read),]
  nbReads=nrow(allSeg)

  stats=c(nbReads,nbReads-sum(allSeg$segType>0))


  nbAdd=0
  y.add=NULL
  max.y=nbReads
  x1.add=NULL
  x2.add=NULL
  couleur.add=NULL
  text.add=NULL
  text.y.add=NULL
  if(!is.null(addDir)){
    for(i in 1:length(addDir)){
      addFile=paste(addDir[i],"/",fileName,".tab",sep="")
      if(file.exists(addFile)){
        add=read.table(addFile,h=T)
        add=add[order(add$inf1,add$sup1,add$read),]
        nbAdd=nbAdd+nrow(add)+1
        x1.add=c(x1.add,add$inf1)
        x2.add=c(x2.add,add$sup1)
        new.max=max.y+1+nrow(add)
        y.add=c(y.add,(max.y+2):new.max)

        cols=ifelse(add$multi,"orange","red")
        cols=ifelse(add$segType>0,"grey40",cols)
        couleur.add=c(couleur.add,cols)

        if(!is.null(dirNames)){
          dataName=dirNames[i]
        }
        else{
          dataName=paste("additional dataset",i)
        }
        text.add=c(text.add,paste("(above) additional segmentations from",dataName))
        text.y.add=c(text.y.add,max.y+1)
        max.y=new.max

        stats=c(stats,nrow(add),nrow(add)-sum(add$segType>0))
      }
      else{
        stats=c(stats,0,0)
      }
    }
  }

  index.doublon=which(allSeg$doublon==TRUE)

  xtitle=paste("Position on bac",bac)
  ytitle="Reads"

  if(!is.null(outputDir)){
    if(!file.exists(outputDir)){
      system(paste("mkdir ",outputDir))
    }
    pdf(paste(outputDir,"/",fileName,".pdf",sep=""),width = 6, height = 8, onefile = FALSE)
  }

  ##par(mar=c(5,4,1.5,1), oma=c(0,0,2,0))
  y=1:nbReads
  couleur=ifelse(allSeg$multi,"orange","red")
  couleur=ifelse(allSeg$segType>0,"grey40",couleur)

  ## Plot only DRJ1
  plot(c(infDRJ1,supDRJ1),c(1,nbReads+nbAdd),type="n",xlab=xtitle,ylab=ytitle,main=fileName)
  mtext(paste("drj1 (",supDRJ1-infDRJ1+1, "pb)",sep=""), adj=0)

  if(!is.null(addSimi)){
    mismatch=read.table(paste(addSimi,"/",fileName,".tab",sep=""))[,1]
    pcId=100-sum(mismatch)/(sum(mismatch+1))*100
    ## plotter les différences
    vmismatch=infDRJ1+which(mismatch!=0)
    abline(v=vmismatch,col="grey80")
    text(supDRJ1,1,paste(signif(pcId,3)," %id (drj)",sep=""),adj=c(1,0),col="grey30")
    ## fenetre glissante
    xpc=NULL
    ypc=NULL
    l=19
    for (i in 1:(length(mismatch)-l)){
      xpc=c(xpc,i+floor(l/2))
      ypc=c(ypc,100-sum(mismatch[i:(i+l)]>0)/(l+1)*100)
      }
    ypcNorm=ypc*(nbReads+nbAdd)/100
    lab=seq(0,100,by=10)
    at=lab*(nbReads+nbAdd)/100
    lines(xpc+infDRJ1,ypcNorm,col="lightsteelblue3")
    axis(4,at=at,labels=lab,col="lightsteelblue3",col.axis="lightsteelblue3",cex.axis=0.8)
    mtext("%id", col="lightsteelblue3", adj=1)
  }

  segments(allSeg$inf1,y,allSeg$sup1,y, col=couleur)

  points(allSeg$sup1[index.doublon],index.doublon,pch="*")

  if(nbAdd>0){
    abline(h=text.y.add,lty=2,col="grey")
    text(infDRJ1, text.y.add,labels=text.add,cex=0.5,pos=4)
    segments(x1.add,y.add,x2.add,y.add, col=couleur.add)
  }

  mtext(paste(fileName," : ",length(unique(allSeg$frag))," fragments différents",sep=""), side=3, line=0, adj=0.5, outer=TRUE)

  if(!is.null(outputDir)){
    dev.off()
  }

  res=c(fileName,as.character(stats))
  return(res)

}


compareSummary2=function(drjFile,inputDir,outputDir=NULL,outputFile=NULL,dirNames=NULL,addSimi=NULL){
  ## ici pas de refDir, mais un fichier avec des coord de DRJ et une liste de répertoire dans inputDir
  ## permet de comparer sur un meme plot plusieurs ensembles de segmentations, on fait l'hypothese que ce sont exactement les memes coordonnees de drj
  ## output : soit outputDir : 1 fichier pdf par paire de drj
  ##          soit outputFile : 1 fichier pdf pour toutes les drjs

  drjs=read.table(drjFile,h=T)
  drjNames=paste(drjs$bac,drjs$inf1,drjs$sup1,drjs$inf2,drjs$sup2,sep="_")

  if(!is.null(outputFile)){
    pdf(outputFile,width = 6, height = 8)
  }

  draw=lapply(drjNames, function(x) drawAndCompareByFile2(x,inputDir,outputDir,dirNames,addSimi))

  if(!is.null(outputFile)){
    dev.off()
  }

  ## formate les stats de nb de reads par paire de DRJ, sous la forme d'un tableau
  nbData=length(inputDir)
  tmp=matrix(unlist(draw),ncol=1+nbData*2,byrow=T)
  tmp2=as.data.frame(matrix(as.numeric(tmp[,-1]),ncol=nbData*2,byrow=F))
  if(!is.null(dirNames)){
    names(tmp2)=paste(rep(dirNames,rep(2,nbData)),c(".tot",".good"),sep="")
  }
  else{
    name=paste("data",1:nbData,sep="")
    names(tmp2)=paste(rep(name,rep(2,nbData)),c(".tot",".good"),sep="")
  }

  tmp2$name=tmp[,1]
  n=ncol(tmp2)
  stats=tmp2[,c(n,1:(n-1))]

  return(stats)

}


## Plot d'une figure de summary sur laquelle on peut superposer plusieurs jeux de données
## il faut que les fichiers de jonction aient le même nom (bac_inf1... .tab) mais dans des répertoires différents
## renvoie également un vecteur avec les nombre de read
drawAndCompareByFile2=function(drjName,inputDir,outputDir=NULL,dirNames=NULL,addSimi=NULL){

  coord.DRJs=strsplit(drjName,"_")[[1]]
  bac=coord.DRJs[1]
  infDRJ1=as.numeric(coord.DRJs[2])
  supDRJ1=as.numeric(coord.DRJs[3])
  infDRJ2=as.numeric(coord.DRJs[4])
  supDRJ2=as.numeric(coord.DRJs[5])


  nbSeg=0
  stats=NULL
  y.seg=NULL
  max.y=-1
  x1.seg=NULL
  x2.seg=NULL
  couleur=NULL
  text.seg=NULL
  y.text=NULL
  x.doublon=NULL
  y.doublon=NULL
  for(i in 1:length(inputDir)){
    segFile=paste(inputDir[i],"/",drjName,".tab",sep="")
    if(file.exists(segFile)){
      seg=read.table(segFile,h=T)
      seg=seg[order(seg$inf1,seg$sup1,seg$read),]
      nbSeg=nbSeg+nrow(seg)+1
      x1.seg=c(x1.seg,seg$inf1)
      x2.seg=c(x2.seg,seg$sup1)
      new.max=max.y+1+nrow(seg)
      new.y=(max.y+2):new.max
      y.seg=c(y.seg,new.y)
      x.doublon=c(x.doublon,seg$sup1[seg$doublon==TRUE])
      y.doublon=c(y.doublon,new.y[seg$doublon==TRUE])

      cols=ifelse(seg$multi,"orange","red")
      cols=ifelse(seg$segType>0,"grey40",cols)
      couleur=c(couleur,cols)

      if(!is.null(dirNames)){
        dataName=dirNames[i]
      }
      else{
        dataName=paste("dataset",i)
      }
      text.seg=c(text.seg,paste("(above) segmentations from",dataName))
      y.text=c(y.text,max.y+1)
      max.y=new.max

      stats=c(stats,nrow(seg),nrow(seg)-sum(seg$segType>0))
    }
    else{
      stats=c(stats,0,0)
    }
  }

  xtitle=paste("Position on bac",bac)
  ytitle="Reads"

  if(nbSeg>0){
    if(!is.null(outputDir)){

      if(!file.exists(outputDir)){
        system(paste("mkdir ",outputDir))
      }
      pdf(paste(outputDir,"/",drjName,".pdf",sep=""),width = 6, height = 8, onefile = FALSE)

    }
    ##par(mar=c(5,4,1.5,1), oma=c(0,0,2,0))

    ## Plot only DRJ1
    plot(c(infDRJ1,supDRJ1),c(0,nbSeg),type="n",xlab=xtitle,ylab=ytitle,main=drjName)
    mtext(paste("drj1 (",supDRJ1-infDRJ1+1, "pb)",sep=""), adj=0)

    if(!is.null(addSimi)){
      mismatch=read.table(paste(addSimi,"/",drjName,".tab",sep=""))[,1]
      pcId=100-sum(mismatch)/(sum(mismatch+1))*100
      ## plotter les différences
      vmismatch=infDRJ1+which(mismatch!=0)
      abline(v=vmismatch,col="grey80")
      text(supDRJ1,0,paste(signif(pcId,3)," %id (drj)",sep=""),adj=c(1,0),col="grey30")
      ## fenetre glissante
      xpc=NULL
      ypc=NULL
      l=19
      for (i in 1:(length(mismatch)-l)){
        xpc=c(xpc,i+floor(l/2))
        ypc=c(ypc,100-sum(mismatch[i:(i+l)]>0)/(l+1)*100)
      }
      ypcNorm=ypc*(nbSeg)/100
      lab=seq(0,100,by=10)
      at=lab*(nbSeg)/100
      lines(xpc+infDRJ1,ypcNorm,col="lightsteelblue3")
      axis(4,at=at,labels=lab,col="lightsteelblue3",col.axis="lightsteelblue3",cex.axis=0.8)
      mtext("%id", col="lightsteelblue3", adj=1)
    }

    abline(h=y.text,lty=2,col="grey")
    text(infDRJ1, y.text,labels=text.seg,cex=0.5,pos=4)
    segments(x1.seg,y.seg,x2.seg,y.seg, col=couleur)
    points(x.doublon,y.doublon,pch="*")

    if(!is.null(outputDir)){
      dev.off()
    }
  }

  res=c(drjName,as.character(stats))
  return(res)

}
