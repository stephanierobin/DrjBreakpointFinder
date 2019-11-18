

getAllValidSegmentationsByDRJPair=function(segmentationFile,drjPairsFile,pairDir,confirmedDrjPairFile=NULL,multipleLog=NULL){

  segmentTab=read.table(segmentationFile,h=T)
  #coordTab=read.table(coordFile)
  #listeReads=unique(as.character(coordTab$V3))
  listeReads=unique(as.character(segmentTab$read))

  ## multipleLog : pour garder une trace des reads éliminés car multiples ex-aequo = bonnes segmentations mais on ne sait pas laquelle choisir
  if(!is.null(multipleLog)){
    con=file(multipleLog,"w")
    close(con)
  }
  
  listeChoosenId=NULL
  tabRead=data.frame(id=NULL,read=NULL,multi=NULL)
  
  ## Selection d'une seule segmentation par read :
  for (r in listeReads){
    ## listeId=coordTab$V1[as.character(coordTab$V3)==r]
    ## seg=segmentTab[is.element(segmentTab$id,listeId),]
    seg=segmentTab[segmentTab$read==r & segmentTab$stat==1,]
    if(nrow(seg)==1){
      listeChoosenId=c(listeChoosenId,seg$id)
      tabRead=rbind(tabRead,data.frame(id=seg$id,read=r,multi=FALSE))
    }
    if(nrow(seg)>1){
      choosenId=chooseMulti(seg,multipleLog)
      if(!is.null(choosenId)){
        listeChoosenId=c(listeChoosenId,choosenId)
        tabRead=rbind(tabRead,data.frame(id=choosenId,read=r,multi=TRUE))
      }
    }
  }

  ## Récupérer les bonnes infos à mettre dans les tableaux
  ## id frag read doublon multi inf1 sup1 inf2 sup2 stat segType FvsR
  segmentTab$segType=rep(0,nrow(segmentTab))
  segmentTab$segType[segmentTab$nbAlt>0]=1
  segmentTab$segType[segmentTab$nbDiff>0]=2
  segmentTab$segType[segmentTab$nbAlt>0 & segmentTab$nbDiff>0]=3
  
  segment2=merge(segmentTab[,c("id","br.beg1","br.end1","br.beg2","br.end2","stat","segType")],tabRead,by="id")
  names(segment2)=c("id","inf1","sup1","inf2","sup2","stat","segType","read","multi")

  ## TODO : gérer les doublons et ajouter colonne FvsR
  segment2$frag=unlist(lapply(as.character(segment2$read),function(x) {removeStrandSuffix(x)}),use.names=F)
  segment2$doublon=rep(FALSE,nrow(segment2))
  segment2$FvsR=rep(0,nrow(segment2))
  segment2=addDoublonInfo(segment2)
  
  ## ajout des debuts et fins des reads sur les drjs :
  ## names(coordTab)[1:7]=c("id","bac","read","infRead1","supRead1","infRead2","supRead2")
  ## segment3=merge(segment2,coordTab[,c("id","infRead1","supRead1","infRead2","supRead2")],by="id")
  ## segmentFinal=segment3[,c("id","frag","read","doublon","multi","infRead1","supRead1","inf1bk","sup1bk","infRead2","supRead2","inf2bk","sup2bk","stat","segType","FvsR")]

  segmentFinal=segment2[,c("id","frag","read","doublon","multi","inf1","sup1","inf2","sup2","stat","segType","FvsR")]
  
  ## découpage par paire de DRJs et écriture des résultats
  drjPairs=read.table(drjPairsFile,h=T)
  drjPairs$nbSupport=rep(0,nrow(drjPairs))
  drjPairs$nbSupportRead=rep(0,nrow(drjPairs))
  drjPairs$nbInit=rep(0,nrow(drjPairs))
  for (i in 1:nrow(drjPairs)){
    listeId=as.numeric(strsplit(as.character(drjPairs$idList[i]),",")[[1]])
    drjPairs$nbInit[i]=length(unique(listeId))
    tab=segmentFinal[is.element(segmentFinal$id,listeId),]
    if(nrow(tab)>0){
      drjPairs$nbSupport[i]=length(unique(tab$frag))
      drjPairs$nbSupportRead[i]=nrow(tab)
      fileName=paste(pairDir,"/",drjPairs$bac[i],"_",drjPairs$inf1[i],"_",drjPairs$sup1[i],"_",drjPairs$inf2[i],"_",drjPairs$sup2[i],".tab",sep="")
      write.table(tab,fileName,quote=F,row.names=F)
    }
  }

  if(!is.null(confirmedDrjPairFile)){
    conf=drjPairs[drjPairs$nbSupport>0,c("bac","inf1","sup1","inf2","sup2","nbSupport","nbSupportRead","nbInit")]
    ## nbSupport : nb de frag differents qui supportent la drj
    ## nbSupportRead : nb de reads differents qui supportent la drj >= nbSupport
    ## nbInit : nb de reads initiaux mappé sur cette paire de DRJ (attention : nb d'id, plusieurs fois le meme read sur une paire de drj ??)
    ## la différence nbInit-nbSupportRead : les reads avec stat segmentation =0 ou les reads mutliples ex-aequo
    conf=conf[order(conf$bac,conf$inf1,conf$inf2),]
    write.table(conf,confirmedDrjPairFile,quote=F,row.names=F)
  }
}


chooseMulti=function(segmentForARead,multipleLog=NULL){

  id=NULL
  sel=segmentForARead[segmentForARead$stat==1,]
  ## il faut stat =1 
  if(nrow(sel)>1){
    sc=sel$nbAlt+sel$nbDiff
    ## on veut un score mini
    min.sc=min(sc)
    index=which(sc==min.sc)

    sel2=sel[index,c("bac","br.beg1","br.end1","br.beg2","br.end2")]
    
    if(nrow(unique(sel2))==1){ ## si un seul id ou bien si une seule paire de région
      id=sel$id[index[1]]
    }
    else{## si plusieurs avec score mini, n'en prend aucun
      if(!is.null(multipleLog)){
        con=file(multipleLog,"a")
        writeLines(paste(sel$read[1],paste(sel$id[index],collapse=","),sep="\t"),con)
        close(con)
      }
    }
  }
  return(id)

}

removeStrandSuffix=function(readName){

  n=nchar(readName)
  suffix=substr(readName,n-2,n)
  if(is.element(suffix,c("RM1","FM1"))){
    frag=substr(readName,1,n-3)
  }
  else{
    frag=readName
  }
  return(frag)

}

addDoublonInfo=function(segmentTab){
  
  dupliques=segmentTab$frag[duplicated(segmentTab$frag)]
    
  for (f in dupliques){
    segmentTab$doublon[segmentTab$frag==f]=TRUE
    tmp=segmentTab[segmentTab$frag==f,c("inf1","sup1","inf2","sup2")]
    if(sum(tmp[1,]==tmp[2,])<4){
      segmentTab$FvsR[segmentTab$frag==f]=2
    }
    else{
      segmentTab$FvsR[segmentTab$frag==f]=1
    }
  }
  
  return(segmentTab)
}
