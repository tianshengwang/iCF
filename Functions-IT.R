
# ############################################################################
# Pruning and Size Selection Based on LeBlanc and Crowley (JASA, 1992)
# ############################################################################


# =================================
# METHOD I: THE TEST SAMPLE METHOD
# =================================

# -----------------------------------------------------------------------------
# THE PRUNING ALGORITHM GOOD FOR THE TREE SIZE SELECTION VIA TEST SAMPLE METHOD
# -----------------------------------------------------------------------------

prune.size.testsample <- function(tre)
{
  call <- match.call(); out <- match.call(expand = F)
  out$result <- out$size  <- out$... <- NULL
  ntest <- as.numeric(tre[1, ncol(tre)])
  if(is.null(dim(tre))) stop("No Need to Prune Further.")
  result <- NULL; n.tmnl <- sum(is.na(tre[,4])); subtree <- 1            
  a <- cbind(Ga.2=2, Ga.3=3, Ga.4=4, Ga.BIC=log(ntest))
  max.Ga <- rep(-1e20, 4); size <- rep(0, 4); btree <-as.list(1:4) 
  while (n.tmnl > 1 ) {
    # print(tre)
    internal <- tre$node[!is.na(tre$score)]; l <- length(internal); 
    r.value <- 1:l
    for(i in 1:l) {
      branch <- tre[is.element(tre$node,c(internal[i], de(internal[i], tree=tre))),]
      score <- as.numeric(as.vector(branch$score))
      r.value[i] <- sum(score, na.rm=T) / sum(!is.na(score))
    }
    alpha <- min(r.value)
    nod.rm <- internal[r.value == alpha]; 
    # if (length(nod.rm)>1) print("Multiple Nodes will be pruned. Check!")
    G <- sum(as.numeric(as.vector(tre$score.test)), na.rm=T); 
    Ga <- G - a*l 
    for (k in 1:4){if (Ga[k] > max.Ga[k]) {max.Ga[k] <- Ga[k]; size[k] <- n.tmnl; btree[[k]] <- tre}}                        
    result <- rbind(result, cbind(subtree=subtree, node.rm=nod.rm, size.tree=nrow(tre), 
                                  size.tmnl=nrow(tre)-l, alpha=alpha, G=G, Ga))
    tre <- tre[!is.element(tre$node, de(nod.rm,tre)),]
    tre[match(nod.rm, tre$node), c("bvar.col", "vname", "best.cut", "score", "score.test")] <- NA
    n.tmnl <- sum(is.na(tre$cut))
    if (n.tmnl ==1) {for (k in 1:4){if (0 > max.Ga[k]) {max.Ga[k] <- 0; size[k] <- 1; btree[[k]] <- tre}}}
    subtree <- subtree + 1          
  }
  # HANDLE THE NULL TREE WITH THE ROOT NODE ONLY
  result <- rbind(result, cbind(subtree=subtree, node.rm='NA', size.tree=nrow(tre), 
                                size.tmnl=1, alpha=9999, G=0, Ga=cbind(Ga.2=0, Ga.3=0, Ga.4=0, Ga.BIC=0)))     
  result <- as.data.frame(result)
  out$result <- result; out$size <- size; out$btree <- btree
  out 
}


# =================================
# METHOD II: THE BOOTSTRAP METHOD
# =================================

# -----------------------------------------------------------------------
# THE PRUNING ALGORITHM GOOD FOR THE BOOTSTRAP TREE SIZE SELECTION METHOD
# -----------------------------------------------------------------------

prune.size <- function(tre)
{
  if(is.null(dim(tre))) stop("No Need to Prune Further.")
  result <- NULL; n.tmnl <- sum(is.na(tre$score)); subtree <- 1            
  while (n.tmnl > 1 ) {
    # if (n.tmnl==5) {btre <- tre; print(btre)}
    internal <- tre$node[!is.na(tre$score)]; l <- length(internal); 
    r.value <- 1:l
    for(i in 1:l) {
      branch <- tre[is.element(tre$node,c(internal[i], de(internal[i], tree=tre))),]
      score <- as.numeric(as.vector(branch$score))
      r.value[i] <- sum(score, na.rm=T) / sum(!is.na(score))
    }
    alpha <- min(r.value)
    nod.rm <- internal[r.value == alpha]; 
    # if (length(nod.rm)>1) print("Multiple Nodes will be pruned. Check!")
    G <- sum(as.numeric(as.vector(tre$score)), na.rm=T);
    G.test <- sum(as.numeric(as.vector(tre$score.test)), na.rm=T)
    result <- rbind(result, cbind(subtree=subtree, node.rm=nod.rm, size.tree=nrow(tre), 
                                  size.tmnl=nrow(tre)-l, alpha=alpha, G=G, G.test=G.test))
    tre <- tre[!is.element(tre$node, de(nod.rm,tre)),]
    tre[match(nod.rm, tre$node), c("bvar.col", "vname", "best.cut", "score", "score.test")] <- NA
    n.tmnl <- sum(is.na(tre$score))
    subtree <- subtree + 1          
  }
  # HANDLE THE NULL TREE WITH THE ROOT NODE ONLY
  result <- rbind(result, cbind(subtree=subtree, node.rm='NA', size.tree=nrow(tre), 
                                size.tmnl=1, alpha=9999, G=0, G.test=0))    
  result <- as.data.frame(result)
  result
}


# ==========================================================================================  #
#  FUNCTIONS RELATED TO THE PRUNING AND THEN BOOTSTRAP FOR TREE SIZE SELECTION
# ==========================================================================================  #

# OPTION LeBlanc IS TO APPLY THE WHOLE SAMPLE (TRUE) OR THE OUT-OF-BAD SAMPLE (FALSE) IN THE BOOTSTRAP PROCEDURE
# OPTION min.boot.tree.size IS TO MAKE SURE A NON-NULL TREE IS GROWN AT EACH BOOTSTRAP 11/1/2007

bootstrap.grow.prune <- function(B=30, data, min.node.size=20, n0=5,  max.depth=10,
                                 col.y=NULL, col.trt=NULL, cols.split.var, cols.nominal=NULL,
                                 SSS=TRUE, a=10, q=0.01, n.intervals=1,
                                 mtry=length(cols.split.var), LeBlanc=TRUE, min.boot.tree.size=1)  
{
  call <- match.call(); out <- match.call(expand = F)
  out$boot.tree <- out$boot.prune <- out$... <- NULL
  time.start <- date()
  tree0 <- grow.IT(data, data, min.node.size=min.node.size, n0=n0, max.depth=max.depth, 
                   SSS=SSS, a=a, q=q, n.intervals=n.intervals,
                   col.y=col.y, col.trt=col.trt, 
                   cols.split.var=cols.split.var, cols.nominal=cols.nominal, mtry=mtry);  
  print(tree0);  
  prune0 <- prune.size(tree0); 
  boot.tree <- list(tree0); boot.prune <- list(prune0) 
  b <- 1
  while (b <= B) {
    print(paste("###################### b = ", b, " ###########################", sep=""))
    # SAMPLING OBSERVATION
    samp <- sample(1:nrow(data), size=nrow(data), replace=T) 
    dat <- data[samp, ];     
    dat.oob <- data[-unique(samp),]
    n.oob <- nrow(dat.oob); # print(n.oob)        
    if (LeBlanc) { tre <- grow.IT(dat, data, min.node.size=min.node.size, n0=n0, max.depth=max.depth, 
                                  SSS=SSS, a=a, q=q, n.intervals=n.intervals,
                                  col.y=col.y, col.trt=col.trt, 
                                  cols.split.var=cols.split.var, cols.nominal=cols.nominal, mtry=mtry)}  
    else {tre <- grow.IT(dat, dat.oob, min.node.size=min.node.size, n0=n0, max.depth=max.depth, 
                         SSS=SSS, a=a, q=q, n.intervals=n.intervals,
                         col.y=col.y, col.trt=col.trt, 
                         cols.split.var=cols.split.var, cols.nominal=cols.nominal, mtry=mtry)}
    print(tre)        
    if (nrow(tre)> min.boot.tree.size) {
      boot.tree <- c(boot.tree, list(tre)); 
      prune <- prune.size(tre); # print(prune)
      boot.prune <- c(boot.prune, list(prune));
      b <- b+1
    }
  }
  time.end <- date(); 
  print(paste("The Start and End time for ", B, "bootstrap runs is:"))
  print(rbind(time.start, time.end))
  out$boot.tree <- boot.tree
  out$boot.prune <- boot.prune
  # THE INITIAL LARGE TREE
  out$initial.tree <- tree0;    
  out
}   




# ========================================================================
# FUNCTION obtain.btree() OBTAINS THE BEST SUBTREE WITH KNOW SIZE bsize=
# ========================================================================

obtain.btree <- function(tre, bsize=3)
{	   
  if (bsize==1) { btre <- tre[1,]; btre[, c("bvar.col", "vname", "best.cut", "score", "score.test")] <- NA}  
  else if (bsize <1) stop("THE BEST TREE SIZE bsize= MUST BE >=1!")
  else {
    n.tmnl <- sum(is.na(tre$score)); indicator <- TRUE         
    while (n.tmnl > 1 && indicator ==TRUE) {
      # print(n.tmnl)
      if (n.tmnl==bsize) {btre <- tre; print(btre); indicator==FALSE}
      internal <- tre$node[!is.na(tre$score)]; L <- length(internal); 
      r.value <- 1:L
      for(i in 1:L) {
        branch <- tre[is.element(tre$node,c(internal[i], de(internal[i], tree=tre))),]
        score <- as.numeric(as.vector(branch$score))
        r.value[i] <- sum(score, na.rm=T) / sum(!is.na(score))
      }
      alpha <- min(r.value)
      nod.rm <- internal[r.value == alpha]; 
      tre <- tre[!is.element(tre$node, de(nod.rm,tre)),]
      tre[match(nod.rm, tre$node), c("bvar.col", "vname", "best.cut", "score", "score.test")] <- NA
      n.tmnl <- sum(is.na(tre$score))        
    }
  }
  btre
}


# ===============================================================
# SELECT THE BEST SUBTREE SIZE AND OBTAIN THE BEST SUBTREE MODEL
# ===============================================================

bootstrap.size <- function(bootstrap.grow.prune, n, plot.it=TRUE, filename=NULL, horizontal=T)
{   
  OUT <- as.list(NULL)
  tree0 <- bootstrap.grow.prune$boot.tree[[1]]
  boot.prune <- bootstrap.grow.prune$boot.prune
  #  COMPUTE THE ALPHA PRIME'S
  prune0 <- boot.prune[[1]] 
  n.subtree <- nrow(prune0)
  alpha <- as.numeric(as.character(prune0$alpha));
  # temp <- c(alpha[1], alpha[-length(alpha)])  	##############
  temp <- c(0, alpha[-length(alpha)])  			############## CHANGE FIRST VALUE OF ALPHA TO 0
  alpha.prime <- sqrt(alpha*temp)  
  # cbind(alpha,  alpha.prime=prune0$alpha.prime)
  b <- length(boot.prune)
  G <- as.numeric(as.character(prune0$G)); 
  ################
  size.tmnl <- as.numeric(as.character(prune0$size.tmnl)); 
  subtree <- as.numeric(as.character(prune0$subtree)); 
  # tree.penalty <- log(nrow(teeth))
  G.a <- matrix(0, n.subtree, 5)
  penalty <- c(0, 2:4, log(n))
  for (i in 1:n.subtree) {
    a.i <- alpha.prime[i]
    bias <- 0
    for (j in 2:b){
      prune.bs <- boot.prune[[j]]
      alpha.bs <- as.numeric(as.character(prune.bs$alpha)); 
      g <- as.numeric(as.character(prune.bs$G)); 
      g.test <- as.numeric(as.character(prune.bs$G.test)); 
      indx <- 1
      if (sum(alpha.bs <= a.i)>0) {          
        temp1 <- which.max(which(alpha.bs<=a.i))
        indx <- ifelse(is.null(temp1), 1, temp1)
      }
      temp2 <- (g-g.test)[indx]
      bias <- bias + temp2 
      # print(cbind(i, a.i, j, bias, indx, temp2))
    }
    G.honest <- G[i] - bias/(b-1) 
    G.a[i,] <- G.honest - penalty*(size.tmnl[i]-1)
  }
  out <- data.frame(cbind(size.tmnl, G.a))
  colnames(out) <- c("tmnl.size", "G", "G.2", "G.3", "G.4", "G.log(n)")
  G.a <- out
  
  # PLOT THE G.a WITH DIFFERENT CHOICES OF a
  subtree.size <- G.a[,1]
  if (plot.it) {
    if (!is.null(filename)) postscript(file=filename, horizontal=horizontal)
    par(mfrow=c(1, 1), mar=rep(4, 4))   ##################### SET THE PLOTTING PARAMETERS
    n.subtrees <- nrow(G.a)
    min.x <- min(subtree.size); max.x <- max(subtree.size)
    min.y <- min(G.a$"G.log(n)"); max.y <- max(G.a$G.2)
    plot(x=c(min.x, max.x), y=c(min.y, max.y), type="n", xlab="tree size", ylab="G(a)")
    for (j in 3:6) lines(subtree.size, G.a[,j], lty=j-1, col=j-1, lwd=2)
    legend(x=min.x, y=(max.y+min.y)/2, lty=2:5, col=2:5, legend=c("G(2)", "G(3)", "G(4)", "G(ln(n))")) 
    if (!is.null(filename)) dev.off()
  }   
  # OBTAIN THE BEST TREE SIZE 
  bsize <- btree <- as.list(NULL)
  Ga.cols <- c("G.2", "G.3", "G.4", "G.log(n)")
  for (j in Ga.cols) {
    best.size <- subtree.size[which.max(G.a[,j])]
    bsize[[j]] <- best.size 
    # print(j); print(best.size); print(tree0); print(G.a)    
    btree[[j]] <- obtain.btree(tree0, bsize=best.size)
  }
  OUT$G.a <- G.a; OUT$bsize <- bsize; OUT$btree <- btree
  return(OUT)
}	




#