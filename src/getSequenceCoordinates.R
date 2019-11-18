## Fonction pour obtenir les coordonnées des séquences à aligner pour la segmentation à partir des paires de hits
getSequenceCoordinates=function(hitPairsFile,id2drjFile,coordFile){

  hits=read.table(hitPairsFile,h=T)
  id2drj=read.table(id2drjFile,h=T)
  
  sel=hits[(hits$orient==1 & hits$sup2.1-hits$inf2.1>0 & hits$sup2.2-hits$inf2.2>0)|(hits$orient==-1 & hits$sup2.1-hits$inf2.1>0 & hits$sup2.2-hits$inf2.2>0),]
  
  sel1=sel[sel$orient==1,] 
  coord1=data.frame(id=sel1$id,bac=sel1$bac,read=sel1$read,inf1=sel1$inf2.2-sel1$inf1.2,sup1=sel1$sup2.2+sel1$l-sel1$sup1.2,inf2=sel1$inf2.1-sel1$inf1.1,sup2=sel1$sup2.1+sel1$l-sel1$sup1.1,orient=sel1$orient,over1=sel1$inf1.2,over2=sel1$sup1.1)
  
  sel2=sel[sel$orient==-1,]
  coord2=data.frame(id=sel2$id,bac=sel2$bac,read=sel2$read,inf1=sel2$inf2.1-sel2$l+sel2$sup1.1,sup1=sel2$sup2.1+sel2$inf1.1,inf2=sel2$inf2.2-sel2$l+sel2$sup1.2,sup2=sel2$sup2.2+sel2$inf1.2,orient=sel2$orient,over1=sel2$l-sel2$sup1.1,over2=sel2$l-sel2$inf1.2)
  
  coordP=rbind(coord1,coord2)

  ## Verifier que les drjs sont bien inclues entierement dans les sequences
  names(id2drj)[3:6]=paste("drj",names(id2drj)[3:6],sep="")
  coordP=merge(coordP,id2drj,sort=F)

  coordP$inf1=pmin(coordP$inf1,coordP$drjinf1)
  coordP$sup1=pmax(coordP$sup1,coordP$drjsup1)
  coordP$inf2=pmin(coordP$inf2,coordP$drjinf2)
  coordP$sup2=pmax(coordP$sup2,coordP$drjsup2)

  coordP=coordP[,c("id","bac","read","inf1","sup1","inf2","sup2","orient","over1","over2")]
  ## return(coordP)
  write.table(coordP,coordFile,row.names=F,col.names=F,quote=F)
    
}

