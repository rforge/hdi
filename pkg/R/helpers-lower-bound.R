printClusterBound <- function(x){
  cat("\n lower l1-norm group bounds in a hierarchical clustering ")
  cat("\n lower bound on l1-norm of all regression coefficients:",
      signif(max(x$lowerBound),4))
  if(sum(x$isLeaf)==1){
    tmp <- sum(x$noMembers[which(x$isLeaf)])
    cat("\n only 1 significant non-overlapping cluster with", tmp,
      if(tmp==1) "member" else "members")
    }else{
      cat("\n number of non-overlapping significant clusters       :",
          sum(x$isLeaf))
      tmp <- range(x$noMembers[which(x$isLeaf)])
      if(diff(tmp)==0){
        cat("\n with", tmp[1],
            if(tmp[1]==1) "member in each non-overlapping cluster"
            else "members in each non-overlapping cluster")
      }else{
        cat("\n with", paste(tmp,collapse=" up to "), "members each")
    }
  }
  cat("\n ")
}


plotClusterBound <- function(x,cexfactor=1,yaxis="members",col=NULL){

    hh <- x$noMembers
    hh2 <- sqrt(x$lowerBound)

    if(yaxis!="members"){
        hh2 <- sqrt(x$noMembers)
        hh <- x$lowerBound
    }
    hh2 <- hh2/max(hh2)
    xvec <- x$position

    col <- if(yaxis=="members")  rgb(0.8,0.2,0.2,0.6) else rgb(0.2,0.2,0.8,0.6)

    plot(xvec, hh,cex=1,axes=FALSE,xlab="",ylab=if(yaxis=="members") "cluster size" else "lower l1-norm bound",pch=20,col="white",log= if(yaxis=="members")"y"else "")
    
    axis(2)
    coll <- rgb(0.1,0.1,0.1,0.7)
    for (k in length(x$position):1){
        if((nm <- x$leftChild[k])>0){
            lines( c(xvec[k],xvec[nm]),c(hh[k],hh[k]),col=coll)
            lines( c(xvec[nm],xvec[nm]),c(hh[k],hh[nm]),col=coll)
        }
        if((nm <- x$rightChild[k])>0){
            lines( c(xvec[k],xvec[nm]),c(hh[k],hh[k]),col=coll)
            lines( c(xvec[nm],xvec[nm]),c(hh[k],hh[nm]),col=coll)
        }
    }
    
    points(xvec, hh,cex=sqrt(hh2)*4*cexfactor,xlab="",col="white",pch=20)
    points(xvec, hh,cex=sqrt(hh2)*4*cexfactor,xlab="",col=col,pch=20)
    
}
