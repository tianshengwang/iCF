
# ####################################################################
# FUNCTIONS FOR RANDOM FORESTS OF INTERACTION TREES
# ####################################################################

require(randomForest)
require(MASS)
require(dplyr)

expit <- function(x) (tanh(x/2)+1)/2  	# CORRECT & VERIFIED


# ======================================
# FUNCTION THAT GENERATES DATA FOR ITE
# ======================================

rdat <- function(n=300, p=5, rho=0.5, Model="IV", sigma=1, digit=2){
  # GENERATE COVARIATES
  if (rho==0){
    for (j in 1:p) assign(paste("x", j, sep=""), round(runif(n), digits=digit))
  } else {
    Sigma0 <- matrix(rho, p, p)
    # REFERENCE http://comisef.wikidot.com/tutorial:correlateduniformvariates
    Sigma0 <- 2*sin(pi*Sigma0/6) 
    diag(Sigma0) <- 1
    mu0 <- rep(0, p)
    X <- pnorm(mvrnorm(n=n, mu=mu0, Sigma=Sigma0))	
    for (j in 1:p) assign(paste("x", j, sep=""), X[,j])
  }
  # GENERATE TREATMENT 
  trt <- sample(c(0,1), size=n, replace = T)
  
  # THE POTENTIAL OUTCOMES
  # MODEL I: THE PURE NOISE MODEL
  mu0 <- -2 - 2 *x1 - 2*x2^2 + 2*x3^3  # STRONG / WEAK PROGNOSTIC EFFECTS
  alpha <- rnorm(n, sd=1)   # COMMON RANDOM EFFECT 
  y0 <- mu0 + alpha + rnorm(n, sd=sigma) 
  if (Model=="I") {delta <- 5}  # NULL MODEL I 
  else if (Model=="II") {delta <- -5 + 5*x1 + 5*x2} # MODEL II: LINEAR - MADE STRONGER  
  else if (Model=="III") {delta <- -5 + 5*x4 + 5*x5} # MODEL III: LINEAR - SEPARATE PROGNOSTIC WITH PREDICTIVE 
  else if (Model=="IV") {delta <- -2 + 2* sign(x1 <= 0.5) + 2*sign(x2 <= 0.5) * sign(x3 <= 0.5)}  # MODEL IV TREE  
  else if (Model=="V") {delta <-  -6 + 0.1*exp(4*x1) + 4*expit(20*(x2-0.5)) + 3*x3 + 2*x4 + x5}	  # MARS MODEL V - MARS 
  else if (Model=="VI") {delta <- -10 + 10*sin(pi*x1*x2) + 20*(x3-.5)^2 + 10*x4 + 5*x5} # MARS MODEL VI - MARS
  else stop("Sorry but the argument for Model= must be I--VI!!")
  mu1 <- mu0 + delta
  y1 <- mu1 + alpha + rnorm(n, sd=sigma)
  
  # OBSERVED OUTCOMES
  y <- y1*trt  + y0*(1-trt)	
  dat <- data.frame(id=1:n, y=y, trt=trt, x1=x1, x2=x2, x3=x3, x4=x4, x5=x5, y1=y1, y0=y0, delta=delta)
  return(dat)
} 

# dat <- rdat(n=100, p=5, rho=0.525, Model="III", sigma=1, digit=2)
# round(dat$delta, digits=4)



# ===========================================
# THE t TEST STATISTIC IN GREEDY SEARCH (GS)
# ===========================================

t2.test <- function(cutpoint, y, trt, x, n0=3, detail=FALSE)
{  
  score <- NA; n <- length(y)
  z <- sign(x <= cutpoint)    
  y0L <- y[trt==0 & z==1]; y0R <- y[trt==0 & z==0]
  y1L <- y[trt==1 & z==1]; y1R <- y[trt==1 & z==0]
  n0L <- length(y0L); n0R <- length(y0R); n1L <- length(y1L);n1R <- length(y1R); 
  if (detail) print(cbind(n0L, n0R, n1L, n1R))
  t0 <- NA
  if (min(n0L, n0R, n1L, n1R) >= n0) {    
    ybar0L <- mean(y0L); ybar0R <- mean(y0R); 
    ybar1L <- mean(y1L); ybar1R <- mean(y1R);
    if (detail) print(c(ybar0L, ybar0R, ybar1L, ybar1R)) 
    DinD <- (ybar1L - ybar0L) - (ybar1R - ybar0R); 
    if (detail) print(DinD)
    mse <- (sum(y^2) - n1L*ybar1L^2 - n1R*ybar1R^2 - n0L*ybar0L^2 - n0R*ybar0R^2)/(n-4)
    if (detail) print(mse) 
    t0 <- DinD/sqrt(mse *(1/n0L + 1/n0R + 1/n1L + 1/n1R))
    score <- t0^2
  }
  return(score)
}

# A REWRITTEN ONE WHERE z IS THE 0(R)-1(L) NODE-ASSIGNING INDICATOR FUNCTION 
t2.interaction <- function(y, trt, z, n0=3, detail=FALSE)
{  
  score <- NA; n <- length(y)
  y0L <- y[trt==0 & z==1]; y0R <- y[trt==0 & z==0]
  y1L <- y[trt==1 & z==1]; y1R <- y[trt==1 & z==0]
  n0L <- length(y0L); n0R <- length(y0R); n1L <- length(y1L);n1R <- length(y1R); 
  if (detail) print(cbind(n0L, n0R, n1L, n1R))
  t0 <- NA
  if (min(n0L, n0R, n1L, n1R) >= n0) {    
    ybar0L <- mean(y0L); ybar0R <- mean(y0R); 
    ybar1L <- mean(y1L); ybar1R <- mean(y1R);
    if (detail) print(c(ybar0L, ybar0R, ybar1L, ybar1R)) 
    DinD <- (ybar1L - ybar0L) - (ybar1R - ybar0R); 
    if (detail) print(DinD)
    mse <- (sum(y^2) - n1L*ybar1L^2 - n1R*ybar1R^2 - n0L*ybar0L^2 - n0R*ybar0R^2)/(n-4)
    if (detail) print(mse) 
    t0 <- DinD/sqrt(mse *(1/n0L + 1/n0R + 1/n1L + 1/n1R))
    score <- t0^2
  }
  return(score)
}




# ===========================================
# THE SMOOTHED t TEST STATISTIC WITH SSS
# ===========================================

t2.SSS <- function(cutpoint, a=10, y, trt, x, detail=F){
  score <- NA; n <- length(y)
  s <- expit(a*(x-cutpoint))
  # SAMPLE SIZES
  n1 <- sum(trt); n0 <- n-n1
  n1L <- sum(trt*s); n1R <- n1-n1L; n0L <- sum((1-trt)*s); n0R <- n0-n0L
  # SUMS OF RESPONSE VALUES
  S1 <- sum(trt*y); S0 <- sum(y) - S1
  S1L <- sum(y*trt*s); S1R <- S1-S1L; S0L <- sum(y*(1-trt)*s); S0R <- S0-S0L
  # MEAN RESPONSES
  ybar1L<-S1L/n1L; ybar1R<-S1R/n1R; ybar0L<-S0L/n0L; ybar0R<-S0R/n0R;
  if (detail) print(c(ybar1L, ybar1R, ybar0L, ybar0R))
  if (detail) print(c(n1L, n1R, n0L, n0R))
  # DIND
  DinD <- (ybar1L - ybar0L) - (ybar1R - ybar0R); 
  if (detail) print(DinD)
  # COMPUTE MSE OR SIGMA SQUARED
  mse <- (sum(y^2) - n1L*ybar1L^2 - n1R*ybar1R^2 - n0L*ybar0L^2 - n0R*ybar0R^2)/(n-4)
  t0 <- DinD/sqrt(mse *(1/n0L + 1/n0R + 1/n1L + 1/n1R))
  score <- t0^2
  return(score)
}


# ==========================
# SSS SPLIT ON ONE VARIABLE
# ==========================

split.SSS <- function(y, trt, x, a=10, 
                      q = 0.01,     					# PERCENTILE 
                      n.intervals=4)					# n.intervals MUST BE >= 1. SETTING IT TO 1 WOULD NULLIFY THE MULTI-START STRATEGY
{
  LB <- quantile(x, probs = q); 
  UB <- quantile(x, probs =1-q); 	
  n.intervals <- round(n.intervals)
  if (n.intervals < 1) stop("n.intervals MUST BE AN INTEGER >= 1.")
  B <- seq(LB, UB, length.out=n.intervals+1)
  Q.max <- -1e15; c.hat.SSS <- NA
  for (b in 1:n.intervals) {
    fit.tmp <-  optimize(t2.SSS, lower=B[b], upper=B[b+1], maximum=TRUE, 
                         a=a, y=y, trt=trt, x=x, detail=FALSE)
    # print(fit.tmp)
    if (!is.nan(fit.tmp$objective) && fit.tmp$objective > Q.max) {
      Q.max <- fit.tmp$objective
      c.hat.SSS <- fit.tmp$maximum
    }
  }				
  return(list(bcut=c.hat.SSS, Q.max=Q.max))
}




# ===========================================================================
# THE power.set() FUNCTION PROVIDES THE POWER SET FOR A CATEGORICAL VARIABLE
# ===========================================================================

power.set <- function(x) {
  if(length(x) == 0) return(vector(mode(x), 0))
  x <- sort(unique(x)); n <- length(x); K <- NULL
  for(m in x) K <- rbind(cbind(K, FALSE), cbind(K, TRUE))
  out <- apply(K, 1, function(x, s) s[x], s = x)
  out <- out[-c(1, length(out))]
  l <- length(out); i <- 1
  out[!sapply(out, length)>=ceiling(n/2+.5)]
}


# THIS FUNCTION as.numeric.factor() CONVERTS FACTOR INTO NUMERIC
as.numeric.factor <- function(x){as.numeric(levels(x))[x]}

# ========================================================================================
# FUNCTION prepare.nominal() PREPARES NOMINAL VARIABLES (WITH MORE THAN TWO LEVELS) BY FIRST 
# MERGING THE LEVELS OF A NOMRINAL VARIABLE WITH FEW OBSERVATIONS AND THEN SORTING LEVELS 
# OF A NOMINAL VARIABLES IN ASCENDING ORDER OF TREATMENT EFFECTS
# ========================================================================================

prepare.nominal <- function(y, trt, x, n.00=0){
  z <- x
  # MERGE LEVELS IF NECESSARY	
  dat.tmp <- data.frame(aggregate(y, by=list(x, trt), FUN=length),
                        aggregate(y, by=list(x, trt), FUN=mean)$x)
  colnames(dat.tmp) <- c("level", "trt", "n", "ybar")
  dat1.tmp <- dat.tmp[dat.tmp$trt==1, ]
  dat0.tmp <- dat.tmp[dat.tmp$trt==0, ]
  dat.tmp <- merge(dat1.tmp, dat0.tmp, by =c("level"), all=TRUE) 
  dat.tmp$diff <- dat.tmp$ybar.x - dat.tmp$ybar.y
  dat.tmp$n.x[is.na(dat.tmp$n.x)] <- 0
  dat.tmp$n.y[is.na(dat.tmp$n.y)] <- 0
  ybar1 <- dat.tmp$ybar.x; ybar0 <- dat.tmp$ybar.y; 
  diff <- dat.tmp$diff
  level.new <- dat.tmp$level 
  Eligible <-  (dat.tmp$n.x > n.00 & dat.tmp$n.y > n.00)
  for (i in 1:NROW(dat.tmp)){
    eligible <- Eligible; eligible[i] <- FALSE   
    level.candidate <- dat.tmp$level[eligible]
    level.replace <- dat.tmp$level[i]
    if (dat.tmp$n.x[i] ==0 && dat.tmp$n.y[i] > n.00) {	
      level.replace <- level.candidate[which.min(abs(dat.tmp$ybar.y[i]-dat.tmp$ybar.y[eligible]))] 
    } else if (dat.tmp$n.y[i] ==0 && dat.tmp$n.x[i] > n.00) {
      level.replace <- level.candidate[which.min(abs(dat.tmp$ybar.x[i]-dat.tmp$ybar.x[eligible]))] 
    } else if (dat.tmp$n.y[i] <= n.00 || dat.tmp$n.x[i] <= n.00) {
      if(is.na(dat.tmp$diff[i])) {
        dat.tmp$diff[i] <- ifelse(is.na(dat.tmp$ybar.x[i]), dat.tmp$ybar.y[i], dat.tmp$ybar.x[i])
        level.replace <- level.candidate[which.min(abs(dat.tmp$diff[i]-dat.tmp$diff[eligible]))]
      } else {
        level.replace <- level.candidate[which.min(abs(dat.tmp$diff[i]-dat.tmp$diff[eligible]))] 
      }
    }
    level.new[i] <- level.replace
    z[z==dat.tmp$level[i]] <- level.replace
  }
  dat.tmp$level.new <- level.new 
  info.pre.merge <- dat.tmp
  
  # ORDINALIZATION 
  tbl <- data.frame(aggregate(y[trt==1], by=list(z[trt==1]), FUN=length), 
                    aggregate(y[trt==1], by=list(z[trt==1]), FUN=mean)$x, 
                    aggregate(y[trt==0], by=list(z[trt==0]), FUN=length)$x, 
                    aggregate(y[trt==0], by=list(z[trt==0]), FUN=mean)$x)
  colnames(tbl) <- c("level", "n.1", "ybar1", "n.0", "ybar0")
  tbl$diff <- tbl$ybar1 - tbl$ybar0
  info.post.merge <- tbl
  levels.sorted <- tbl$level[order(tbl$diff)] 
  z <- as.numeric(ordered(z, levels=levels.sorted))   # z2 IS NUMERICAL 
  return(list(info.pre.merge=info.pre.merge, info.post.merge=info.post.merge, z=z))
}

if (FALSE){
  n <- 20
  dat <- rdat1(n=n, beta=rep(2, 6), sigma=1, cut1=.5, cut2=.5, nominal=TRUE, m=6); head(dat)
  result <- prepare.nominal(dat$y, dat$trt, dat$x5, n.00=0); result
}





# ===============================
# ONE SINGLE SPLIT OF DATA
# ===============================

# WHEN USING FOR RANDOM FORESTS, SET test=NULL. 
# min.node.size= SETS THE MINIMUM NUMBER OF OBSERVATIONS FOR CLAIMING A TERMINAL NODE
# n0= SETS THE MINIMUM NUMBER OF OBSERVATIONS FOR (n11, n10, n01, n00)
# cols.cols.split.var= ASKS FOR THE COLUMNS OF SPLITTING VARIABELS, INCLUDING BOTH CONTINUOUS AND CATEGORICAL ONES
# max.depth= SPECIFIED THE MAXIMUM HEIGHT OF A TREE (ANOTHER WAY OF STOPPING THE GROWTH).
# mtry= SPECIFIES THE NUMBER OF COVARIATES IN THE RANDOMLY SELECTED SUBSETS FOR SPLITTING

partition.IT <- function(dat, test=NULL, name="1", 
                         min.node.size=20, n0=5, max.depth=15, 
                         SSS=TRUE, a=10, q=0.01, n.intervals=3,	
                         col.y=NULL, col.trt=NULL, cols.split.var, cols.nominal=NULL, 
                         mtry=length(cols.split.var))
{   
  call <- match.call(); out <- match.call(expand = F)
  out$info <- out$name.l <- out$name.r <- out$left <- out$right <- out$... <- NULL
  name.l <- paste(name, 1, sep=""); name.r <- paste(name, 2, sep="")
  n <- nrow(dat); 
  if (!is.null(test)) {n.test <- NA; score.test <- NA;}  ########## TEST SAMPLE ########  
  bvar.col <- vname <- best.cut <- NA; max.score <- -1e20;   
  if (!is.null(col.y)) y <- dat[, col.y] else y <- dat$y; 
  if (!is.null(col.trt)) trt <- dat[, col.trt] else trt <- dat$trt;  
  vnames <- colnames(dat)
  # COMPUTE THE TREATMENT EFFECT IN CURRENT NODE
  trt.effect <- NA; n.1 <- sum(trt==1); n.0 <- n-n.1
  # if (min(n.1, n.0)> n0) 
  trt.effect <- mean(y[trt==1]) - mean(y[trt==0])
  # CONTROL THE MAX TREE DEPTH
  depth <- nchar(name) 
  if (depth <= max.depth && n >= min.node.size) {
    m.try <- ifelse(is.null(mtry), length(cols.split.var), mtry)  
    for(j in sample(cols.split.var, size=m.try, replace=FALSE)) {
      x <- dat[,j]; v.name <- vnames[j]; 
      xx <- x
      if (is.element(j, cols.nominal)) xx <- prepare.nominal(y, trt, x, n.00=0)$z   # DEAL WITH CATEGORICAL VARIABLES
      temp <- sort(unique(xx)); 
      if (length(temp)==1){
        score.x <- NA 	# UNARY VARIABLE
      } else if (length(temp)==2) {
        # CHARACTER-VALUES ARE AUTOMATICALLY SORTED IN THE ALPHABETICAL ORDER.   
        z <- sign(xx==temp[1]) 		
        score.x <- t2.interaction(y, trt, z, n0=n0, detail=FALSE)
        cut.x <- temp[1]
      } else {
        if (SSS) {
          # SMOOTH SIGMOID SURROGATE (SSS)
          cut.x <- split.SSS(y, trt, xx, a=a, q=q, n.intervals=n.intervals)$bcut
          score.x <- t2.interaction(y, trt, z=sign(xx<=cut.x), n0=n0, detail=FALSE)
          if (is.element(j, cols.nominal)) cut.x <- paste(sort(unique(x[xx <=cut.x])), collapse=" ")
        } else {
          # GREEDY SEARCH (GS)
          zcut <- temp[-length(temp)]
          score.x <- -1e20;
          for(k in zcut) {
            score <- NA; 
            score <- t2.interaction(y, trt, z=sign(xx<=k), n0=n0, detail=FALSE)
            if (!is.na(score) && score >= score.x) {
              score.x <- score; cut.x <- k
              if (is.element(j, cols.nominal)) cut.x <- paste(sort(unique(x[xx <=cut.x])), collapse=" ")
            } 
          }
        }
      }
      if (!is.na(score.x) && score.x >= max.score) {
        max.score <- score.x; bvar.col <- j; 
        vname <- v.name; best.cut <- cut.x
      }
    }
  }
  if (!is.null(test)) { 
    n.test <- NROW(test); score.test <- NA;
    if (!is.na(best.cut)) {
      if (is.element(bvar.col, cols.nominal)){			
        var.split <- as.character(test[, bvar.col])
        cut1 <- unlist(strsplit(as.character(best.cut), split=" ")) ######## NOMINAL PREDICTOR
        z.test <- sign(is.element(var.split, cut1))
      } else  {z.test <- sign(test[, bvar.col] <= best.cut)}   
      score.test <- t2.interaction(test$y, test$trt, z=z.test, n0=n0/2, detail=FALSE)
      if (!is.na(score.test)){
        out$name.l <- name.l; out$name.r <- name.r
        out$left.test <- test[z.test==1,  ]
        out$right.test <- test[z.test==0,  ]
        if (is.element(bvar.col, cols.nominal)){	######## NOMINAL PREDICTOR
          z.data <- is.element(as.character(dat[,bvar.col]), cut1)
          out$left <- dat[z.data==1, ]
          out$right <- dat[z.data==0, ]
        } else {	
          out$left  <- dat[dat[,bvar.col] <= best.cut,]
          out$right <- dat[dat[,bvar.col] > best.cut, ]
        }
      } else {bvar.col <- NA; vname <- NA; best.cut <- NA;  max.score <- NA}
    }
    out$info <- data.frame(node=name, n=n, n.1=n.1, n.0=n.0, trt.effect=trt.effect,
                           bvar.col = bvar.col, vname=vname, best.cut= best.cut, 
                           score=ifelse(max.score==-1e20, NA, max.score), 
                           score.test, n.test=n.test)
  } else {
    if (!is.na(best.cut)) {
      out$name.l <- name.l; out$name.r <- name.r
      if (is.element(bvar.col, cols.nominal)){	   ######## NOMINAL PREDICTOR
        cut1 <- unlist(strsplit(as.character(best.cut), split=" "))
        z.data <- is.element(as.character(dat[,bvar.col]), cut1)
        out$left <- dat[z.data==1, ]
        out$right <- dat[z.data==0, ]
      } else {	
        out$left  <- dat[dat[,bvar.col] <= best.cut,]
        out$right <- dat[dat[,bvar.col] > best.cut, ]
      }
    }		
    out$info <- cbind(node=name, n = n, n.1=n.1, n.0=n.0, trt.effect=trt.effect,
                      bvar.col = bvar.col, vname=vname, best.cut= best.cut, 
                      score=ifelse(max.score==-1e20, NA, max.score))
  }	
  return(out) 
}


if (FALSE){
  dat <- rdat1(n=n, beta=rep(2, 6), sigma=1, cut1=.5, cut2=.5, nominal=TRUE, m=6)
  head(dat)
  test <- rdat1(n=100, beta=rep(2, 6), sigma=1, cut1=.5, cut2=.5, nominal=TRUE, m=6)
  pt <- partition.IT(dat, test=test, name="1", 
                     min.node.size=20, n0=5, max.depth=15, 
                     SSS=TRUE, a=10, q=0.01, n.intervals=3,				
                     cols.split.var, cols.nominal=7, mtry=length(cols.split.var))
  pt$info
}


# =================================================
# THE grow.IT() FUNCTION CONSTRUCTS A LARGE TREE 
# =================================================

grow.IT <- function(data, test=NULL, 
                    min.node.size=20, n0=5, max.depth=15, 
                    SSS=TRUE, a=10, q=0.01, n.intervals=3,				
                    col.y=NULL, col.trt=NULL, cols.split.var, cols.nominal=NULL, 
                    mtry=length(cols.split.var))
{
  out <- list.nd <- list.test <- temp.list <- temp.test <- temp.name <- NULL
  list.nd <- list(data); 
  if (!is.null(test)) list.test <- list(test)
  name <- 1
  while (length(list.nd)!=0) {    
    for (i in 1:length(list.nd)){
      # print(i)
      if (!is.null(dim(list.nd[[i]])) && nrow(list.nd[[i]]) > 1){ 
        test0 <- NULL
        if (!is.null(test)) test0 <- list.test[[i]]
        split <- partition.IT(list.nd[[i]], test0, name[i], min.node.size=min.node.size, n0=n0, 
                              max.depth=max.depth, SSS=SSS, a=a, q=q, n.intervals=n.intervals,				
                              col.y=col.y, col.trt=col.trt, cols.split.var=cols.split.var, 
                              cols.nominal=cols.nominal, mtry=mtry)
        out <- rbind(out, split$info)
        if (!is.null(split$left) && is.null(test)) {
          temp.list <- temp.test <- c(temp.list, list(split$left, split$right))
          temp.name <- c(temp.name, split$name.l, split$name.r)
          
        } else if (!is.null(split$left) && !is.null(test) && !is.null(split$left.test)) {
          temp.list <- c(temp.list, list(split$left, split$right))
          temp.name <- c(temp.name, split$name.l, split$name.r)
          temp.test <- c(temp.test, list(split$left.test, split$right.test))
        }
      }
    }
    list.nd <- temp.list; list.test <- temp.test; name <- temp.name
    temp.list <- temp.test <- temp.name <- NULL
  }   
  out <- as.data.frame(out)    
  out$node <- as.character(out$node)
  out <- out[order(out$node), ]
  out
}


# EXAMPLE
if (FALSE){
  n <- 400
  dat <- rdat1(n=n, beta=rep(2, 6), sigma=1, cut1=.5, cut2=.5, nominal=TRUE, m=6)
  test <- rdat1(n=200, beta=rep(2, 6), sigma=1, cut1=.5, cut2=.5, nominal=TRUE, m=6)
  
  cols.split.var <- c(1, 2, 3, 4, 7); cols.nominal <- 7
  tree0 <- grow.IT(data=dat, test=test, 
                   min.node.size=20, n0=5, max.depth=15, 
                   SSS=TRUE, a=10, q=0.01, n.intervals=3,				
                   cols.split.var, cols.nominal, mtry=length(cols.split.var))
  tree0 
}



# ===================================================================
# PLOTTING IT TREE STRUCTURE, MODIFIED FROM PETER CALHOUN'S CODES
# ===================================================================

plot.tree <- function(tree, cols.nominal=NULL, 
                      textDepth=3, lines="rectangle", digits=4)
{
  depth<-max(nchar(tree[,1]))
  par(xaxs='i')
  par(mar=c(1,1,1,1))
  par(xpd=TRUE)
  plot(1, type="n", xlab="",ylab="",xlim=c(0,1),ylim=c(0,1), axes=FALSE,xaxs="i",yaxs="i")
  nodes <- tree$node
  nodesBin <- gsub("1", "0", nodes)
  nodesBin <- gsub("2", "1", nodesBin)
  lastObs<-nchar(nodesBin)
  nodesBin<-substr(nodesBin,2,lastObs)
  var <- tree$vname 
  col.var <- tree$bvar.col 
  cut <- as.character(tree$best.cut) 
  size <- tree$n  
  effect <- tree$trt.effect
  
  for(i in 1:length(nodesBin)){
    nChar<-nchar(nodesBin[i])
    if(!is.na(var[i])){
      if(lines=="rectangle"){
        lines(c((1/(2^(nChar+2)))*(4*strtoi(nodesBin[i], base = 2)+1),
                (1/(2^(nChar+2)))*(4*strtoi(nodesBin[i], base = 2)+3)),
              c((depth-nChar)/(depth+1),(depth-nChar)/(depth+1)))
        
        lines(c((1/(2^(nChar+2)))*(4*strtoi(nodesBin[i], base = 2)+1),
                (1/(2^(nChar+2)))*(4*strtoi(nodesBin[i], base = 2)+1)),
              c((depth-nChar)/(depth+1),(depth-nChar-1)/(depth+1)))
        
        lines(c((1/(2^(nChar+2)))*(4*strtoi(nodesBin[i], base = 2)+3),
                (1/(2^(nChar+2)))*(4*strtoi(nodesBin[i], base = 2)+3)),
              c((depth-nChar)/(depth+1),(depth-nChar-1)/(depth+1)))
        
      } else if(lines=="triangle"){
        lines(c((1/(2^(nChar+2)))*(4*strtoi(nodesBin[i], base = 2)+2),
                (1/(2^(nChar+2)))*(4*strtoi(nodesBin[i], base = 2)+1)),
              c((depth-nChar)/(depth+1),(depth-nChar-1)/(depth+1)))
        
        lines(c((1/(2^(nChar+2)))*(4*strtoi(nodesBin[i], base = 2)+2),
                (1/(2^(nChar+2)))*(4*strtoi(nodesBin[i], base = 2)+3)),
              c((depth-nChar)/(depth+1),(depth-nChar-1)/(depth+1)))
      }         
      
      if(nChar <= textDepth){ 
        if (is.element(col.var[i], cols.nominal)){
          cutpoint <- unlist(strsplit(as.character(cut[i]), split=" "))
          cutpoint <- paste("{", paste(cutpoint, collapse=","), "}", sep="")
          text((1/(2^(nChar+2)))*(4*strtoi(nodesBin[i], base = 2)+2), (depth-nChar)/(depth+1)+1/(depth+20), 
               bquote(.(as.character(var[i]))%in%.(cutpoint)),cex=1, col="blue")
        } else {
          cutpoint <- round(as.numeric(cut[i]), digits=digits)
          text((1/(2^(nChar+2)))*(4*strtoi(nodesBin[i], base = 2)+2), (depth-nChar)/(depth+1)+1/(depth+20), 
               bquote(.(as.character(var[i]))<=.(cutpoint)),cex=1, col="blue")
        }
      }
    } else {
      if(nChar <= textDepth){
        effect.i <- round(as.numeric(effect[i]), digits=digits)
        text((1/(2^(nChar+2)))*(4*strtoi(nodesBin[i], base = 2)+2),(depth-nChar)/(depth+1), 
             paste("n=", size[i], sep=""),cex=1, offset=1, col="red")
        text((1/(2^(nChar+2)))*(4*strtoi(nodesBin[i], base = 2)+2), (depth-nChar)/(depth+1)-0.025, 
             paste("d=", effect.i, sep=""),cex=1, offset=1, col="red")
      }
    }
  }
}

# EXAMPLE
# plot.tree(tree0, cols.nominal=7, textDepth=5, lines="rectangle")


# ==========================================================
# FUNCTION de() FINDS ALL THE DESCENDENTS OF NODE x IN tree
# ==========================================================

de <- function(x, tree)
{
  if(length(x) != 1) stop("The length of x in function de must be 1.")    
  y <- tree$node;  de <- NA
  if(sum(match(x, y), na.rm = T) != 0) {
    temp <- 1:length(y)
    start <- match(x, y) + 1    
    end <- length(y)
    if(start <= length(y) & nchar(y[start]) > nchar(x)) {
      temp1 <- temp[temp >= start & nchar(y) <= nchar(x)][1] - 1
      if(!is.na(temp1)) end <- temp1
      de <- y[start:end]
    }}
  de
}

# de(1112, tree0)


# =========================================================================================
# FUNCTION plot.tree.latex() GENERATES LATEX CODES FOR PLOTTING TREE IN PACKAGE pstricks
# =========================================================================================

plot.tree.latex <- function(tree, file="tree-code.tex", 
                            digits=4, cols.nominal=NULL, landscape=FALSE)
{
  n.node <- nrow(tree)
  sink(file=file)
  #  SET UP THE LATEX FILE
  cat("\\documentclass[12pt]{article} \n")
  cat("\\usepackage{pstricks,pst-node,pst-tree,lscape} \n")
  cat("\\pagenumbering{gobble} \n\n")   # SUPPRESS PAGE NUMBERING
  cat("\\begin{document} \n \n")
  
  # DEFINE SYMBOLS
  if (landscape) cat("\\begin{landscape} \n")
  cat("\\begin{figure} \n \\centering \n")
  cat("\\newcommand{\\Tgcircle}{\\Tcircle[linecolor=black, fillcolor=gray]} \n")
  cat("\\newcommand{\\Tgoval}{\\Toval[linecolor=black, fillcolor=gray]} \n")
  cat("\\newcommand{\\Tgframe}[1]{\\Tr{\\psframebox[linecolor=black, fillstyle=solid, fillcolor=orange]{#1}}} \n")
  # OPTION
  cat("\\psset{nodesep=0.7pt, linecolor=black, treesep=1.2cm, levelsep=1.8cm} \n")
  I0 <- i0 <- NULL
  for (i in 1:n.node) {
    node.i <- as.character(tree[i, 1])
    de.i <- de(node.i, tree)
    blanks <- paste(rep(" ", (nchar(node.i)-1)*8), sep="")  # 8 SPACES IN ONE TAB		
    n.i <- tree$n[i]
    trt.effect <- round(as.numeric(tree$trt.effect[i]), digits=digits)
    lable.i <- ifelse(!is.null(tree$trt.effect), trt.effect, " "); 	##### THIS MEASURE MAY BE DIFFERENT FOR OTHER TYPES OF TREES
    if (!is.na(de.i[1])) {	# INTERNAL NODE
      if (nchar(node.i)==1 ||  substr(node.i, nchar(node.i), nchar(node.i))=="2") 
        cat(blanks, "\\pstree{\\Tgcircle{~~}} \n", sep = "")
      else cat(blanks, "\\pstree{\\Tgcircle{~~} \\tlput{\\color{blue}", rule.i, "\\hspace{-.6in}}} \n", sep = "")				
      cat(blanks, "{ \n", sep = "") 
      I0 <- c(I0, i)
      i0 <- c(i0, i + length(de.i))
      # UPDATE THE SPLITTING RULE
      vname.i <- tree$vname[i]; col.var <- tree$bvar.col[i]; 
      cutpoint <- tree$best.cut[i]
      if (is.element(col.var, cols.nominal)){
        cutpoint <- unlist(strsplit(as.character(cutpoint), split=" "))
        cutpoint <- paste("\\{", paste(cutpoint, collapse=","), "\\}", sep="")
        rule.i <- paste("\\texttt{", vname.i, "}", "$\\,\\in\\,$", cutpoint, sep="") 
      } else {
        # cutpoint <- strtrim(as.character(cutpoint), width=digits); 				
        cutpoint <- round(as.numeric(cutpoint), digits=digits)
        rule.i <- paste("\\texttt{", vname.i, "}", "$\\,\\leq\\,", cutpoint, "$", sep="")
      }
    } else if (substr(node.i, nchar(node.i), nchar(node.i))=="1") { # TERMINAL NODE
      cat(blanks, "\\Tgframe{",  lable.i, "} \\nput{d}{\\pssucc}{\\color{black} \\textbf{", n.i, "}}",
          "\\tlput{\\color{blue} ", rule.i, " \\hspace{-.3in}} \n", sep = "")
    } else cat(blanks, "\\Tgframe{",  lable.i, "} \\nput{d}{\\pssucc}{\\color{black} \\textbf{",  n.i, "}} \n", sep = "")
    if (is.element(i, i0)) {
      rep0 <- rep("}", sum(i0==i))
      node.i0 <- as.character(tree[I0[i0==i][1] , 1])
      blanks0 <- paste(rep(" ", (nchar(node.i0)-1)*8), sep="")  		
      cat(blanks0, rep0, "\n", sep = "") 
    }
  }
  cat("\\end{figure} \n\n")
  if (landscape) cat("\\end{landscape} \n\n")
  cat("\\end{document} \n")
  sink()  
}
# plot.tree.latex(tree0, file="tree-code.tex", digits=5, cols.nominal=7)







# ======================================================
# BUILD RANDOM FORESTS OF INTERACTION TREES (RFIT) 
# ======================================================

Build.RFIT <- function(dat, min.node.size=20, n0=5, max.depth=15,
                       SSS=TRUE, a=10, q=0.01, n.intervals=3,   # SSS OPTIONS
                       col.y=NULL, col.trt=NULL, cols.split.var, cols.nominal=NULL,
                       ntree=500, mtry=max(floor(length(cols.split.var)/3), 1),
                       avoid.nul.tree=TRUE)
{
  out <- as.list(NULL); 
  names(dat)[c(col.y, col.trt)] <- c("y", "trt");
  n <- NROW(dat); 
  out$ID.Boots.Samples <- out$TREES <- as.list(1:ntree)
  out$N.Bn <- matrix(0, ntree, n) 
  b <- 1
  while (b <= ntree) {
    print(b)
    # TAKE BOOTSTRAP SAMPLES
    id.b <- sample(1:nrow(dat), size=nrow(dat), replace = T)	
    dat.b <- dat[id.b,]
    tre.b <- grow.IT(data=dat.b, test=NULL, 
                     min.node.size=min.node.size, n0=n0, max.depth=max.depth, 
                     SSS=SSS, a=a, q=q, n.intervals=n.intervals,				
                     col.y=col.y, col.trt=col.trt, cols.split.var=cols.split.var, cols.nominal=cols.nominal, 
                     mtry=mtry)
    if (avoid.nul.tree) {
      if (nrow(tre.b) > 1) {
        out$ID.Boots.Samples[[b]] <- id.b
        out$N.Bn[b, ] <- as.vector(table(factor(id.b, levels=1:n)))	
        out$TREES[[b]] <- tre.b; 
        b <- b +1			
      }
    }
    else {
      out$ID.Boots.Samples[[b]] <- id.b
      out$N.Bn[b, ] <- as.vector(table(factor(id.b, levels=1:n)))	
      out$TREES[[b]] <- tre.b; 
      b <- b + 1		
    }
  }
  Model.Specification <- as.list(NULL)
  Model.Specification$data <- dat; Model.Specification$cols.split.var <- cols.split.var; 
  Model.Specification$cols.nominal <- cols.nominal; Model.Specification$col.y <- col.y; 
  Model.Specification$col.trt <- col.trt;
  Model.Specification$n <- n;  
  out$Model.Specification <- Model.Specification
  return(out)
}

# EXAMPLE 
# NEED TO FIGURE OUT WHERE THE WARNINGS ARE FROM
if (FALSE) {
  fit.rfit <- Build.RFIT(dat, min.node.size=20, n0=5, max.depth=10,
                         SSS=TRUE, a=10, q=0.01, n.intervals=3,   # SSS OPTIONS
                         col.y=col.y, col.trt=col.trt, cols.split.var=cols.split.var, cols.nominal=cols.nominal,
                         ntree=50, mtry = 3, avoid.nul.tree=TRUE)
  names(fit.rfit)
}



# =============================================================
# FUNCTION combine.RFIT() COMBINES TWO Build.RF.IT() OBJECTS
# =============================================================

combine.RFIT <- function (...){
  rflist <- list(...)
  n.rf <- length(rflist)
  rf1 <- rflist[[1]]
  components <- names(rf1)
  n.comp <- length(components) 
  OUT <- as.list(1:n.comp)
  for (j in 1:3) {
    if (j ==3) comp.j <- NULL
    else comp.j <- as.list(NULL)
    for (i in 1:n.rf) {
      rf.i <- rflist[[i]]
      comp.ij <- rf.i[[j]]
      if (j ==3) comp.j <- rbind(comp.j, comp.ij)  # CONCATENATE THE N.Bn MATRICES
      else comp.j <- c(comp.j, comp.ij) 
    }
    OUT[[j]] <- comp.j
  }
  OUT[[n.comp]] <- rf1$Model.Specification			
  names(OUT) <- components 
  return(OUT)
}



# ------------------------------------------------------------------
# THIS senddown() FUNCTION IS WRITTEN FOR THE VIARIABE IMPORTANCE PART
# USING RANDOM FORESTS 
# ------------------------------------------------------------------

send.down <- function(dat.new, tre, cols.nominal=NULL)
{
  # print(tre)
  cut.point <- as.vector(tre$best.cut); 
  split.var <- as.numeric(as.vector(tre$bvar.col)); 
  dat.node <- rep(1, NROW(dat.new))
  for (i in 1:nrow(tre)){
    in.node <- (dat.node)==(tre$node[i]);                      
    if (!is.na(split.var[i])){
      # print(cbind(i, var=split.var[i], cut=cut.point[i]))
      var.split <- dat.new[,split.var[i]]; 
      cut0 <- cut.point[i]
      if (!is.element(split.var[i], cols.nominal)) { 
        cut0 <- as.numeric(cut0)    
        l.nd <- dat.node[in.node & var.split <= cut0] 
        r.nd <- dat.node[in.node & var.split > cut0]
        dat.node[in.node & var.split <= cut0] <- paste(l.nd, 1, sep="")
        dat.node[in.node & var.split >  cut0] <- paste(r.nd, 2, sep="")  
      }
      else {
        var.split <- as.character(var.split)
        cut0 <- unlist(strsplit(as.character(cut0), split=" ")) #####################
        l.nd <- dat.node[in.node & is.element(var.split, cut0 )] 
        r.nd <- dat.node[in.node & !is.element(var.split, cut0)]                  
        dat.node[in.node & is.element(var.split, cut0)] <- paste(l.nd, 1, sep="")  
        dat.node[in.node & !is.element(var.split, cut0)] <- paste(r.nd, 2, sep="")}                   
    }}
  # OBTAIN PREDICTED TREATMETN EFFECT
  terminals <- sort(unique(dat.node))
  pred.effect <- rep(NA, NROW(dat.new))
  effect <- as.numeric(as.vector(tre$trt.effect))
  for (k in terminals) pred.effect[dat.node==k] <- effect[tre$node==k]
  # data.frame(id=1:NROW(dat.new), node=dat.node, pred.effect=pred.effect)	
  return(pred.effect)
}

# -----------------------------------------------------------------------------------
# THe FUNCTION predict.RFIT() RETURNS PREDICYTED ITE AND SE VIA BIAS-CORRECTED IJ
# -----------------------------------------------------------------------------------

predict.RFIT <- function(fit.rfit, newdata, SE=TRUE, bias.EFRON=TRUE)
{
  out <- as.list(NULL); 
  n1 <- NROW(newdata); n <- fit.rfit$Model.Specification$n
  B <- length(fit.rfit$TREES)
  PRED <- matrix(NA, B, n1)
  cols.nominal <- fit.rfit$Model.Specification$cols.nominal
  for (b in 1:B) PRED[b,] <- send.down(dat.new=newdata, fit.rfit$TREES[[b]], cols.nominal=cols.nominal)
  # print(dim(PRED))
  out$predicted.ite <- apply(PRED, 2, mean)
  out$OTR.II <- apply(PRED, 2, FUN=function(x){sum(x>0)})/B
  if (SE) {
    out$se.ite.naive <- apply(PRED, 2, sd)/sqrt(B)    # VERY NATIVE APPROACH
    N.Bn <- fit.rfit$N.Bn
    
    # WAGER, HASTIE, AND EFRON (2016, JMLR)
    PRED.centered <- scale(PRED, center = TRUE, scale = FALSE)
    N.mat <- (N.Bn-1)/B
    dim(PRED.centered)
    V.IJ <- apply((t(N.mat) %*% PRED.centered)^2, 2, sum)
    bias1 <- (n-1)*apply(PRED, 2, var)/B 
    V.IJu <- V.IJ - bias1
    # print(cbind(V.IJ, bias, V.IJu))
    out$se.ite.uncorrected <- sqrt(V.IJ)
    out$se.ite.bias <- bias1
    out$se.ite.corrected <- sqrt(V.IJu)
    
    # EFRON (2014, JASA) - NEEDS MATRIX Z FOR EACH OBS IN NEWDATA
    # HARDER TO COMPUTE, YIELDING SIMILAR RESULTS 
    if (bias.EFRON) {
      V <- bias2 <- V.u <- rep(0, n1)
      for (k in 1:n1) {
        pred.k <- PRED[, k]   # B by 1				
        pred.k.centered <- scale(pred.k, center=TRUE, scale=FALSE)   # B by 1
        Z.k <- apply(N.Bn-1 , 2,  FUN=function(x){x*pred.k.centered})	 # B by n
        V[k] <- sum(apply(Z.k, 2, mean)^2)  # scalar
        bias2[k] <- sum(apply(Z.k, 2, var))/B
        V.u[k] <- V[k] - bias2[k]
      }
      # print(cbind(V, V.IJ, bias1, bias2, V.IJu, V.u))  # V and V.IJ SHOULD BE IDENTICAL; CHECK TO SEE. 
      out$se.ite.bias.Efron <- bias2 
      out$se.ite.corrected.Efron <- sqrt(V.u)
    }	
  }
  out$OTR.I <- sign(out$predicted.ite >0)
  return(out)
}






# ================================================================================
# PREDICTING ITE USING THE REGRESSION APPROACH - METHOD I
# NOTE THAT I HAVE IMPLEMENTED THE SO-CALLED AIPWE METHOD IN ZHANG ET AL. (2012)
# ================================================================================

# FORMULA SHOULD BE OF FORM y ~ x1 + x2 + x3 + x4 + x5 | trt

# install.packages("randomForest")
require(randomForest)
Predict.ITE.SR <- function(formula, dat, test=test, ntree=500, mtry=NULL)
{
  vars <- all.vars(formula); p <- length(vars)-1  
  trt.name <- vars[length(vars)]
  trt <- dat[, trt.name]
  if (length(unique(trt)) !=2) stop("This method can only handle binary treatments!")
  trt <- as.numeric(as.vector(factor(trt, labels=0:1)))
  if (is.null(mtry)) mtry <- max(floor((p-1)/3), 1)    # ORIGINAL SUGGESTION OF RF
  form0 <- as.formula(paste(vars[1], paste(vars[2:p], collapse=" + "), sep=" ~ "))	
  RF1 <- randomForest(formula=form0, data=dat, subset=(trt==1), ntree=ntree, nodesize=5, mtry=mtry) 
  RF0 <- randomForest(formula=form0, data=dat, subset=(trt==0), ntree=ntree, nodesize=5, mtry=mtry)
  p.y1 <- predict(object=RF1, newdata=test, type="response") 
  p.y0 <- predict(object=RF0, newdata=test, type="response")
  pred.ICE <- p.y1-p.y0
  return(pred.ICE)
}






#