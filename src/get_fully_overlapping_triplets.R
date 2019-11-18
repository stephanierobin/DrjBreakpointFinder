source("hits_functions.R")
source("intervals.R")


getFullyOverlappingTriplets=function(blastFile,readLengthFile,coordFile,drjPairsFile,id2drjFile=NULL,byScaff=FALSE){
  
  ## optim 2018 : if byScaff=TRUE, ie. the blast file contains the results for only one scaffold/bac sequence
  # then : great time speed-up !! (x1800 !)
  # function "couverture" is no longer called
  # instead : computes once for all the intervals of covered-by-hits for the whole scaffold and then calls the function couvertureOptim instead of couverture. 
  # TODO : to implement in the general case where several scaffolds.
  
  ## formatte + rajoute colonne 'l':readLength # fonction definie dans hits_functions.R
  tab=formatBlastTable(blastFile,readLengthFile)
  
  tabF=removeEmbeddedHits(tab,self=F,rename=T) # definie dans hits_functions.R
  ## enleve les hits qui sont completement inclus dans un autre
  ## attention change l'ordre des colonnes : inf2 toujours < sup2

  # print(nrow(tabF))
 
  ## Paramètres des filtres :
  ## ------------------------
  readOverlap=50
  flank=20
  readNotCovered=20
  bacCov=90
  insertSize=300

  ## Filtre 0 : enleve les hits dont le read peut s'aligner entierement
  ## ------------------------------------------------------------------
  entier=unique(tab$read[abs(tab$length-tab$l)<=20 & tab$pcid>=95]) # unique supprime les doubles
  tabF=tabF[!is.element(tabF$read,entier),]
    
  # print(paste("F0 : ",nrow(tabF)))
  if(nrow(tabF)==0){
    return(1)
  }

  ## Filtre 1 : les couples de read-bac avec au moins 2 hits de meme orientation :
  ## --------------------------------------------------------------------------
  tabF1=tabF[tabF$orient==1,]
  tabF2=tabF[tabF$orient==-1,]
  interest1=tabF1[is.element(tabF1$couple,unique(tabF1$couple[duplicated(tabF1$couple)])),] #duplicated retourne un vecteur dont lignes contiennent elements dupliques 
  interest2=tabF2[is.element(tabF2$couple,unique(tabF2$couple[duplicated(tabF2$couple)])),]

  # print(paste("F1-1 : ",nrow(interest1)))
  # print(paste("F1-2 : ",nrow(interest2)))

  if(nrow(interest1)+nrow(interest2)==0){
    return(1)
  }

  ## Filtre 1+2+3 : ordre des hits, chevauchement > 50 nt sur le read, pas de chavauchement sur le BAC
  ## --------------------------------------------------------------------------------------------------
  pair1=getWellOrderedPairs(interest1,seuil=readOverlap,orient=1)
  pair2=getWellOrderedPairs(interest2,seuil=readOverlap,orient=-1)
  totalPairs=rbind(pair1,pair2) # rbind combine vectors par ligne

  # print(paste("F3 : ",nrow(totalPairs)))

  if(nrow(totalPairs)==0){
    return(1)
  }

  ## Filtre 4+5+6+7 : region entre les 2 hits sur le BAC couverte à plus de bacCov % par des reads
  ##                  région non couverte sur le read <  readNotCovered nt
  ##                  régions alignées flanquant la drj chacune > flank nt
  ##                  région entre les 2 drj sur le BAC > insertSize nt
  ## ------------------------------------------------------------------------------------
  cercles=data.frame(bac=totalPairs$bac,inf=pmin(totalPairs$sup2.1,totalPairs$sup2.2),sup=pmax(totalPairs$inf2.1,totalPairs$inf2.2))
  
  if(byScaff){
    couv=couvertureOptim(cercles,tab[tab$pcid>90 & tab$length>50,])
  }
  else{
    couv=couverture(cercles,tab[tab$pcid>90 & tab$length>50,]) # fonction definie plus bas
  }

  totalPairs$cov=couv
  
  selectPairs=totalPairs[abs((totalPairs$sup1.2-totalPairs$inf1.1+1)-totalPairs$l)<readNotCovered & totalPairs$inf1.2-totalPairs$inf1.1>flank & totalPairs$sup1.2-totalPairs$sup1.1>flank & totalPairs$cov>bacCov & totalPairs$totalSize-2*totalPairs$drjSize>insertSize,]
  
  if(nrow(selectPairs)==0){
    return(1)
  }
  
  selectPairs$id=1:nrow(selectPairs)
  selectPairs=selectPairs[,c(ncol(selectPairs),1:(ncol(selectPairs)-1))]

  # print(paste("F7 : ",nrow(selectPairs)))


  ## merge les paires de drj

  res=getFinalDRJpairs(selectPairs,rlist=T) #fonction definie plus bas

  drjs=res$drjs
  selectPairs=res$pairHits
  id2drj=getId2drj(drjs) #fonction definie plus bas 

  coord=getSequenceCoordinates(selectPairs,id2drj) #fonction definie plus bas
  
  write.table(drjs,drjPairsFile,quote=F,row.names=F)
  #write.table(selectPairs,hitPairsFile,quote=F,row.names=F,col.names=T)
  write.table(coord,coordFile,quote=F,row.names=F,col.names=T)
  if(!is.null(id2drjFile)){
    write.table(id2drj,id2drjFile,quote=F,row.names=F,col.names=T)
  }
  
  return(0)

}

## CL 29/08/2012 MODIF :
##   - suppose que pairHits a une colonne id
##   - renvoie la liste des id plutôt que la liste des reads dans rlist
##   - teste si la paire de DRJs finale ne s'overlappe pas : si overlap, suppression de cette paire ET de toutes les paires de hits qui la définissait => UPDATE pairHits (le renvoie aussi)
###   ----------------------------------------

## ici algo mieux que getBacRegionPairs, merge seulement si chevauchement réciproque
## aussi : les régions sont les drj pas les hits totaux
## renvoie des paires de régions sur bac, avec le nombre de hits qui les confirme, le nombre de read, la médiane des drjSize
getFinalDRJpairs=function(pairHits,rlist=F){
# rlist= T : renvoie aussi la liste des reads dans le tableau
  
  tabPlus=pairHits

  ## Définition des coords des DRJs
  inf1=pmin(tabPlus$inf2.1,tabPlus$inf2.2)
  sup2=pmax(tabPlus$sup2.1,tabPlus$sup2.2)
  tabPlus$inf1=inf1
  tabPlus$sup1=inf1+tabPlus$drjSize
  tabPlus$inf2=sup2-tabPlus$drjSize
  tabPlus$sup2=sup2

  ## liste pour updater le tableau pairHits
  idToremove=NULL
    
  regBac=by(tabPlus,as.character(pairHits$bac),function(x) {
    #listeId=1:nrow(x)
    #x$id=listeId
    groupes=as.list(NULL)
    if(nrow(x)>1){
      for(i in 1:nrow(x)){
        id=x$id[i]
        rest=x[-i,]
        t=rest[pmin(rest$sup1,x$sup1[i])-pmax(rest$inf1,x$inf1[i])+1>0 & pmin(rest$sup2,x$sup2[i])-pmax(rest$inf2,x$inf2[i])+1>0,]
        if(nrow(t)==0 | length(groupes)==0){
          groupes[[length(groupes)+1]]=c(id,t$id)
        }
        else{
          g=which(unlist(lapply(groupes,function(y) {sum(is.element(y,c(id,t$id)))>0})))
          if(length(g)==1){
            groupes[[g]]=unique(c(groupes[[g]],id,t$id))
          }
          else{
            if(length(g)==0){
              groupes[[length(groupes)+1]]=c(id,t$id)
            }
            if(length(g)>1){
              groupes[[length(groupes)+1]]=unique(c(id,t$id,unlist(groupes[g])))
              groupes=groupes[-g]
            }
          }
        }
      }
    }
    if(nrow(x)==1){
      groupes[[1]]=x$id
    }

    ## verification : tous les id doivent appartenir a un unique groupe
    if(length(unique(unlist(groupes)))!=nrow(x)){
      print("pb manque id")
      print(x$bac[1])
      print(length(unique(unlist(groupes))))
      print(nrow(x))
    }
    if(length(unique(unlist(groupes)))!=length(unlist(groupes))){
      print("pb non disjoint")
    }

    ## ecriture des paires de DRJs
    inf1=NULL;sup1=NULL;inf2=NULL;sup2=NULL
    nbHits=NULL;nbReads=NULL;drjMed=NULL
    idList=NULL;#drjList=NULL
    for(i in 1:length(groupes)){
      ## ici récupere les coordonnees des deux regions, le nb de paires de hits, le nb de reads (la liste ?), la taille médiane des drj (la liste ?)
      t=x[is.element(x$id,groupes[[i]]),]
      ## coordonnees de la paire de DRJs (merge) 
      inf1=c(inf1,min(t$inf1))
      sup1=c(sup1,max(t$sup1))
      inf2=c(inf2,min(t$inf2))
      sup2=c(sup2,max(t$sup2))

      nbHits=c(nbHits,nrow(t))
      nbReads=c(nbReads,length(unique(t$read)))
      drjMed=c(drjMed,median(t$drjSize))
      ##readList=c(readList,paste(unique(t$read),collapse=","))
      idList=c(idList,paste(unique(t$id),collapse=","))
      ##drjList=c(drjList,paste(t$drjSize,collapse=",")) #inutile
    }
    res=data.frame(bac=rep(x$bac[1],length(inf1)),inf1,sup1,inf2,sup2,nbHits,nbReads,drjMed,idList) #fonction rep génère une suite d un meme nb
        
    return(res)
  })

  df=do.call("rbind", as.list(regBac))

  df$l1=df$sup1-df$inf1+1
  df$l2=df$sup2-df$inf2+1

  ## verification pas d'overlap entre les 2 DRJs et update de pairHits
  df2=df[df$inf2>df$sup1,]
  idToremove=as.numeric(unlist(strsplit(as.character(df$idList[df$inf2<=df$sup1]),","),r=T))
  pairHits2=pairHits[!is.element(pairHits$id,idToremove),]
  
  return(list(drjs=df2,pairHits=pairHits2))

}

## a partir de drjpairs et de sa colonne rList, renvoie le tableau de correspondance inverse : pour chaque id (triplet) sa paire de drj (identifiée par bac et paires de coordonnees)
getId2drj=function(drjPairs){

  res=apply(drjPairs,1,function(x) {
    ids=as.numeric(strsplit(as.character(x[9]),",")[[1]])
    n=length(ids)
    d=data.frame(id=ids,bac=rep(x[1],n),inf1=rep(as.numeric(x[2]),n),sup1=rep(as.numeric(x[3]),n),inf2=rep(as.numeric(x[4]),n),sup2=rep(as.numeric(x[5]),n))
    return(d)
  })

  df=do.call("rbind", as.list(res))
  rownames(df)=1:nrow(df)
  return(df)
  
}

## FUNCTION NOT TESTED YET !!!!

## function that returns a vector of percentages, such that the value i is the % of the regions[i] covered by  hits of tabHits
couverture=function(regions,tabHits){

  res=unlist(apply(regions,1,function(x) {
    bac=as.character(x[1])
    # print(bac)
    reg=data.frame(inf=as.numeric(x[2]),sup=as.numeric(x[3]))
    hits=tabHits[tabHits$bac==bac,c("inf2","sup2")]
    names(hits)=c("inf","sup")
    notCov=notCoveredBy(reg, hits, includeBoundaries=T) #fonction definie dans intervals.R qui renvoie reg qui n est pas inclue dans hits
    if(nrow(notCov)>0){
      frac=sum(notCov$sup-notCov$inf+1)/(reg$sup-reg$inf+1)*100
    }
    else{
      frac=0
    }
    return(100-frac)
    print(frac)
  }))

  return(res)
  print(res)
}


## Fonction pour obtenir les coordonnées des séquences à aligner pour la segmentation à partir des paires de hits
getSequenceCoordinates=function(hitPairs,id2drj){

  ## hits=read.table(hitPairsFile,h=T)
  ## id2drj=read.table(id2drjFile,h=T)
  
  sel=hitPairs[(hitPairs$orient==1 & hitPairs$sup2.1-hitPairs$inf2.1>0 & hitPairs$sup2.2-hitPairs$inf2.2>0)|(hitPairs$orient==-1 & hitPairs$sup2.1-hitPairs$inf2.1>0 & hitPairs$sup2.2-hitPairs$inf2.2>0),]
  
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

  ## Vérifier que les inf>0 (meme si ca ne bloque dans la recuperation de la sequence, on se sert de cette apres pour retrouver les coordonnees sur le bac, para rapport aux coordonnees sur la sequence
  coordP$inf1=pmax(coordP$inf1,1)
  coordP$inf2=pmax(coordP$inf2,1)

  coordP=coordP[,c("id","read","bac","inf1","sup1","inf2","sup2","orient","over1","over2")]
  return(coordP)
  ## write.table(coordP,coordFile,row.names=F,col.names=F,quote=F)
    
}


# 07/2018 New functions to speed-up in the case of one scaffold at a time
# Pre-processing function for couvertureOptim
getHitsIntervals=function(tab){
  
  maxi=max(tab[,2])
  binCov=rep(0,maxi)
  for (i in 1:nrow(tab)){
    b1=tab[i,1]
    e1=tab[i,2]
    binCov[b1:e1]=1
  }
  #binCov : 0/1 vector 1 if the position on scaffold is covered by at least one hit
  
  b=NULL
  e=NULL
  di=diff(binCov)
  if(binCov[1]==1){
    b=1
  }
  b=c(b,which(di==1)+1)
  e=which(di==-1)
  if(binCov[length(binCov)]==1){
    e=c(e,length(binCov))
  }
  return(data.frame(b,e))
}

# Replaces function "couverture" (in the case of only one bac/scaffold in regions and tabHits tables)
## function that returns a vector of percentages, such that the value i is the % of the regions[i] covered by  hits of tabHits
couvertureOptim=function(regions,tabHits){
  
  hitsIntervals=getHitsIntervals(tabHits[,c("inf2","sup2")])
  
  res=unlist(apply(regions,1,function(x) {
    inf1=as.numeric(x[2])
    sup1=as.numeric(x[3])
    if(sum(hitsIntervals$b<=inf1 & hitsIntervals$e>=sup1)>0){
      # the whole region is covered by a single interval
      return(100)
    }
    tmp = hitsIntervals[(hitsIntervals$b>=inf1 & hitsIntervals$b<=sup1) | (hitsIntervals$e>=inf1 & hitsIntervals$e<=sup1),]
    if(nrow(tmp)==0){
      # the whole region does not overlap any interval
      return(0)
    }
    tmp$b=pmax(tmp$b,inf1)
    tmp$e=pmin(tmp$e,sup1)
    covered=sum(tmp$e-tmp$b+1)
    frac=covered/(sup1-inf1+1)*100
    return(frac)
  }))
  
  return(res)
}
