source("Cassis-segmentation-light.R")

segmentASetOfSequences=function(coordFile,vectorDir,figureDir,segResultFile,zoom=T,margin=20,clean=F){

  dir.create(figureDir,showWarnings = FALSE)

  coordTab=read.table(coordFile,h=T)
  res=apply(coordTab,1, function(x) segmentASequence(x,vectorDir,figureDir,zoom,margin,clean))
  res=do.call("rbind",res)

  write.table(as.data.frame(res),segResultFile,quote=F,row.names=F)
}


segmentASequence=function(coordLine,vectorDir,figureDir,zoom=T,margin=20,clean=F){

  id=as.numeric(coordLine[1])
  bac=coordLine[3]
  read=coordLine[2]

  vectorTab=read.table(paste(vectorDir,"/",id,".tab",sep=""))

  if(zoom){
    ## on utilise les coordonnées des chevauchements des hits blast sur le read (attention : pas les coord exactes de la drj) 
    drjInf=as.numeric(coordLine[9])
    drjSup=as.numeric(coordLine[10])
    
    index.inf=max(1,drjInf-margin)
    index.sup=min(nrow(vectorTab),drjSup+margin)
    vectorTab=vectorTab[index.inf:index.sup,]
  }

  if(clean){
    ## for the "partial" method
    ## removes the parts specific to the read at the beginning or at the end
    index.inf=max(c(which(vectorTab$V3==0),which(vectorTab$V4==0),1))
    maxi1=max(vectorTab$V3)
    maxi2=max(vectorTab$V4)
    index.sup=min(c(which(vectorTab$V3==maxi1),which(vectorTab$V4==maxi2),nrow(vectorTab)))
    if(index.inf<index.sup){
      vectorTab=vectorTab[index.inf:index.sup,]
    }
    else{
      print("probleme dans clean(): sequence nulle !!!")
      return(data.frame(id,read,bac,br.beg1=0,br.end1=0,br.beg2=0,br.end2=0,stat=0,nbAlt=0,nbDiff=0))
    }
  }
   
  pdf(paste(figureDir,"/",id,".pdf",sep=""),width = 6, height = 6, onefile = FALSE)
  resSeg=segmentAndPlotABreak(vectorTab)
  dev.off()
  
  N=nrow(vectorTab)
  x1=resSeg[2]
  x2=resSeg[3]
  stat=resSeg[1]
  ## if(x1<2){
  ##   stat=0
  ##   x1=1
  ## }
  ## if(x2>(N-1)){
  ##   stat=0
  ##   x2=N
  ## }

  ## renvoie les coordonnées des points de cassure directement sur le bac
  inf.bac1=as.numeric(coordLine[4])
  ## br.beg1=inf.bac1+vectorTab[x1,3]
  ## br.end1=inf.bac1+vectorTab[x2,3]
  ## changement : a cause des gaps et des valeurs consecutives, ca evite d'avoir des decalages petits entre plusieurs reads sur une meme paire de drj (cf. merge_segments_on_drj_pairs.R
  br.beg1=inf.bac1+vectorTab[x1+1,3]-1
  br.end1=inf.bac1+vectorTab[x2-1,3]+1

  inf.bac2=as.numeric(coordLine[6])
  ## br.beg2=inf.bac2+vectorTab[x1,4]
  ## br.end2=inf.bac2+vectorTab[x2,4]
  br.beg2=inf.bac2+vectorTab[x1+1,4]-1
  br.end2=inf.bac2+vectorTab[x2-1,4]+1

  # attention, cas particuliers (de toute façon dans ces cas, stat=0)
  if(x2<2){
    br.end1=inf.bac1+vectorTab[x2,3]
    br.end2=inf.bac2+vectorTab[x2,4]
  }
  if(x1>(N-1)){
    br.beg1=inf.bac1+vectorTab[x1,3]
    br.beg2=inf.bac2+vectorTab[x1,4]
  }
    
  nbAlt=getNbOfAlternation(vectorTab[x1:x2,])
  nbDiff=getNbOfDifferences(vectorTab[x1:x2,])
  res=data.frame(id,read,bac,br.beg1,br.end1,br.beg2,br.end2,stat,nbAlt,nbDiff)
  return(res)
}


getNbOfAlternation=function(vectorTab){
	r=FALSE
        nb=0
        vec=vectorTab$V1-vectorTab$V2
        if(length(vec)>0){
          pentes=rle(vec[vec!=0])$values
          nb=sum(pentes==1)-sum(pentes[1]==1)
          #nb=sum(pentes==-1)-sum(pentes[length(pentes)]==-1)
        }
	return(nb)
}

getNbOfDifferences=function(vectorTab){
	nb=sum(vectorTab$V1==1 & vectorTab$V2==1)
	return(nb)
}
