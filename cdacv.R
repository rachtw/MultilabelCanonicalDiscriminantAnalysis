args<-commandArgs()
print(args)

distNo<-args[8]
print(length(args))
if (length(args)==9) {
  suffix<-args[9]
  suffix2<-""
} else if (length(args)==10) {
  suffix<-args[9]
  suffix2<-args[10]
} else {
  suffix<-""
  suffix2<-""
}

lm.ridge <- function(formula, data, subset, na.action,
    lambda = 0, model = FALSE, x = FALSE, y = FALSE, contrasts = NULL, ...)
{
    m <- match.call(expand.dots = FALSE)
    m$model <- m$x <- m$y <- m$contrasts <- m$... <- m$lambda <- NULL
    m[[1L]] <- quote(stats::model.frame)
    # m <- eval.parent(m)
      m <- eval(m, parent.frame()) # added
    # print(paste("m=",m))

    Terms <- attr(m, "terms")
    Y <- model.response(m)
    X <- model.matrix(Terms, m, contrasts)
    # print("X=")
    # print(X)
      x <- X # added
    na.action <- attr(m, "na.action")
    xlevels <- .getXlevels(Terms, m)

    z <- .Call(stats:::C_Cdqrls, x, y, tol=1e-07, FALSE)

    n <- nrow(X); p <- ncol(X)
      assign <- attr(X, "assign")  # added
    offset <- model.offset(m)
    if(!is.null(offset)) {
      Y <- Y - offset
    }
    # print("Y=")
    # print(Y)
      y <- Y # added

    if(Inter <- attr(Terms, "intercept"))
    {
      if (dim(X)[1]==1) {
        Xm <- as.matrix(X[, -Inter])
        rownames(Xm) <- colnames(X)[-Inter]
        p <- p - 1
        X <- t(as.matrix(X[, -Inter]))
        colnames(X) <- rownames(Xm) # must assign names after transforming matrix
        rownames(X) <- rownames(Y)
      } else {
        Xcolnames <- colnames(X)[-Inter]
        X <- as.matrix(X[, -Inter])
        colnames(X) <- Xcolnames
        rownames(X) <- rownames(Y)
        Xm <- as.matrix(colMeans(X))
        p <- p - 1
        X <- X - rep(Xm, rep(n, p))
      }
      # print("offset X=")
      # print(X)
      # print("Xm=")
      # print(Xm)
      Ym <- if (is.null(dim(Y))) mean(Y) else colMeans(Y)
      # print("Ym=")
      # print(Ym)
      if (is.null(dim(Y))) Y <- Y - Ym else Y <- Y - matrix(rep(colMeans(Y),each=n),nrow=n)
    } else {
      Ym <- Xm <- NA
    }
    if (is.null(dim(X)) || dim(X)[1]==1) {
      Xscale <- array(0,c(dim(X)[2],1))
    } else {
      Xscale <- as.matrix(drop(rep(1/n, n) %*% X^2)^0.5) # equivalent to colMeans(x^2)^0.5
      # print("Xscale=")
      # print(Xscale)
      X_ <- X/rep(Xscale, rep(n, p))
      # print("X_=")
      # print(X_)
      if (all(is.na(X_))) {
        # don't alter X
      } else {
        X <- X_
      }
    }
    rownames(Xscale) <- colnames(X)
    # print("Xscale=")
    # print(Xscale)
    # print("scaled X=")
    # print(X)
    Xs <- svd(X)              # X = UDV'
    # print("Xs=")    print(Xs)
    rhs <- t(Xs$u) %*% Y      # rhs = U'Y
    # print("rhs=")    print(rhs)
    d <- Xs$d
    # print("d=")    print(d)
    lscoef <-  Xs$v %*% (rhs/d)     # V * (U'Y/D)
    lsfit <- X %*% lscoef
    resid <- Y - lsfit              # Y = V * (U'Y/D) + Resid
    s2 <- sum(resid^2)/(n - p - Inter)
    HKB <- (p-2)*s2/sum(lscoef^2)
    LW <- (p-2)*s2*n/sum(lsfit^2)
    k <- length(lambda)
    dx <- length(d)

    #div <- d^2 + rep(lambda, rep(dx,k))  # D^2
    #a <- drop(d*rhs)/div    # D/(D^2 + lambda * I) U'Y
    #dim(a) <- c(dx, if (is.matrix(y)) ncol(y) else k) # modified
    #coef <- Xs$v %*% a    # V * (D/(D^2 + lambda * I)) U'Y
    #dimnames(coef) <- list(names(Xscale), format(lambda))
    #dimnames(coef)[[1]] <- names(Xscale)
    #dimnames(coef)[[2]] <- colnames(y)

    if (is.null(dim(rhs))) {
      coefarr <- array(0L,dim=c(length(rhs),1,k))
    } else {
      #coefarr <- array(0L,dim=c(dim(rhs)[1],dim(rhs)[2],k))
      coefarr <- array(0L,dim=c(dim(Xs$v)[1],dim(rhs)[2],k))
    }
    GCV <- array(-1L,k)
    for (kk in 1:k) {
      div <- d^2 + lambda[kk] # (D^2 + lambda * I)
      coefarr[,,kk] <- Xs$v %*% (d*rhs/div) # V * (D/(D^2 + lambda * I)) U'Y; rhs = U'Y
      GCV[kk] <- mean((Y - X %*% coefarr[,,kk])^2/(1-sum(d^2/div)/n)^2)
      names(GCV)[kk] <- lambda[kk]
    }
    coef <- coefarr[,,which.min(GCV)] # need checking
    # print("coef=")
    # print(coef)
    div <- d^2 + lambda[which.min(GCV)] # need checking

    # How to scale back?
    # Yridge=H %*% Y
    # Anova uses model.matrix, residual and coefficients; Are they scaled? Anova generate E & H; Are they scaled?
    # candisc calculates raw coefficients using E & H. Features must be "centered" before multiply with raw coefficients
    # So it's reasonable to calculate hat matrix without scaling back first

    coef_<-as.matrix(coef)    
    if (any(Xscale==0)) {
      if (dim(coef_)[1]==dim(Y)[2]) {
        scaledcoef <- coef_
      } else {
        scaledcoef <- t(coef_)
      }
    } else {
      scaledcoef <- t(coef_/rep(Xscale,ncol(coef_)))
    }
    if (Inter) {
        # print("scaledcoef=")
        # print(scaledcoef)
        # print("Xm=")
        # print(Xm)
        inter <- Ym - scaledcoef %*% Xm
        # print("inter=")
        # print(inter)
        if (is.null(dim(scaledcoef))) {
          scaledcoef <- c(inter, scaledcoef)
          names(scaledcoef) <- c("(Intercept)",names(Xscale))
        } else {
          scaledcoef <- cbind(inter, scaledcoef)
          # print("scaledcoef with inter=")
          # print(scaledcoef)
          # print("Xscale=")
          # print(Xscale)
          # print(c("(Intercept)",rownames(Xscale)))
          if (is.null(dim(Xscale))) {
            colnames(scaledcoef) <- c("(Intercept)",names(Xscale))
          } else {
            colnames(scaledcoef) <- c("(Intercept)",rownames(Xscale))
          }
          rownames(scaledcoef) <- colnames(Y)
        }
    }
    #coefficients<-drop(scaledcoef)
    coefficients<-t(scaledcoef)

    ### added (start)
    fitted.values2 <- X %*% coef + Ym
    if(!is.null(offset)) {
      fitted.values2 <- fitted.values2 + offset
    }
    fitted.values <- x %*% coefficients
    residuals <- y - fitted.values
    residuals2 <- y - fitted.values2              # Y = XB + R
    ### added (end)

    #GCV <- colSums((Y - X %*% coef)^2)/(n-colSums(matrix(d^2/div, dx)))^2

    df.residual = nrow(x)-z$rank # added

    res <- list(scales = Xscale,
                Inter = Inter, lambda = lambda, ym = Ym, xm = Xm,
                GCV = GCV, kHKB = HKB, kLW = LW)
    res$inter <- inter
    res$na.action <- na.action
    res$df.residual = df.residual # checked
    res$residuals = residuals # bad... large
    res$residuals2 = residuals2 # worse... larger
    res$coef=drop(coef)
    res$coefficients=coefficients # checked - close
    res$offset = offset
    res$contrasts = contrasts # NULL
    res$call = match.call()
    res$terms = Terms # checked
    res$model = m
    res$X = X
    res$Y = Y
    res$x = x # checked
    res$y = y # checked
    res$assign = assign # checked
    res$fitted.values = fitted.values # checked - close
    res$fitted.values2 = fitted.values2
    res$xlevels <- xlevels
    res$qr <- NULL
    res$rank <- z$rank # checked
    res$d <- d

    class(res) <- "ridgelm"

    res
}

validator.names <- function(model, ...) {
  UseMethod("validator.names")
}

validator.names.default<-function(model, ...){
  validateors <- attr(terms(model), "variables")
  as.character(validateors[3:length(validateors)])
}

# from discproj.R in library(fpc)
tdecomp <- function(m){
  wm <- eigen(m, symmetric=TRUE)
  p <- ncol(m)
  wmd <- wm$values
  sqrtwmd <- try(sqrt(wmd),silent=TRUE)
  if (class(sqrtwmd)!="try-error")
    out <- t(wm$vectors %*% diag(sqrtwmd))
  else
    out <- NA
  out
}

candisc.ridgelm <- function(mod, term, type="2", manova, ndim=rank, E, H, ...) candisc.mlm(mod, term, type="2", manova, ndim=rank, E, H, dfh, ...)

candisc.mlm <- function(mod, term, type="2", manova, ndim=rank, E, H, dfh, ...) {
  if (!inherits(mod, "mlm") && !inherits(mod, "ridgelm")) stop("Not an mlm/ridgelm object")
  #if (missing(manova)) manova <- Anova(mod, type=as.character(type))
  if (!missing(manova)) {
    terms <- manova$terms
    if (missing(term)) term <- terms[1]
    if (missing(E)) E <- manova$SSPE
    if (missing(H)) H <- manova$SSP[[term]]
    dfe <- manova$error.df
    dfh <- manova$df[[term]]
  }
  #dfh <- mod$rank - 1
  dfe <- mod$df.residual
  W <- E / dfe

  Tm <- tdecomp(E)
  eInv <- solve(Tm)
  eHe <- t(eInv) %*% H %*% eInv
  dc <- eigen(eHe, symmetric=TRUE)
  rank <- min(dfh, sum(dc$values>0))
  pct <- 100 * dc$values / sum(dc$values)

  #if (ndim > rank) {
  #  warning(paste("You asked for", ndim, "dimensions, but rank is", rank, ". ndim has been reset to", rank))
  #  ndim <- rank
  #}
  ndim <- rank  # modified

  coeffs.raw <- eInv %*% dc$vectors * sqrt(dfe)
  # should we drop the coeffs corresponding to 0 eigenvalues here or at the end?
  coeffs.raw <- as.matrix(coeffs.raw[,1:ndim])
  rownames(coeffs.raw) <- rownames(H)
  colnames(coeffs.raw) <- cn <- paste('Can', 1:ndim, sep="")

  # These are what SAS calls pooled within-class std. can. coefficients
  coeffs.std <- diag(sqrt(diag(W))) %*% coeffs.raw
  rownames(coeffs.std) <- rownames(H)
  colnames(coeffs.std) <- cn

  data <- model.frame(mod)
  Y <- model.response(data)
  Y <- scale(Y, center=TRUE, scale=FALSE)

  scores <- Y %*% coeffs.raw
  scores <- as.matrix(scores[,1:ndim])
  colnames(scores) <- cn

  # Get the factor(s) corresponding to the term...
  #all.factors <- data[, sapply(data, is.factor), drop=FALSE]
  #factor.names <- unlist(strsplit(term, ":"))
  #factors <- data[factor.names]

  # Canonical means for levels of the factor(s)
  #means <- aggregate(scores, factors, mean)
  #rownames(means) <- do.call(paste,c(means[factor.names],sep=':'))
  #means <- means[, -(1:length(factor.names))]

  # These are what SAS calls total canonical structure coefficients
  # and what I plot as vectors in my canplot macro
  structure <- cor(Y, scores)

  canrsq <- dc$values[1:ndim] / (1+dc$values[1:ndim])

  # make scores into a data frame containing factors in mod
#  scores <- cbind( model.frame(mod)[,-1], scores )
#### FIXME: scores should also include regressors in the model
#  scores <- cbind( all.factors, as.data.frame(scores) )
  #scores <- cbind( model.frame(mod)[validator.names(mod)], as.data.frame(scores) )
  result <- list(
    dfh=dfh,
    dfe=dfe,
    eigenvalues=dc$values, canrsq=canrsq,
    pct=pct, rank=rank, ndim=ndim,
    #means=means,
    #factors=factors, term=term, terms=terms,
    coeffs.raw=coeffs.raw, coeffs.std=coeffs.std,
    structure=structure,
    scores=scores, E=E, H=H
    )
  class(result) <- "candisc"
  result
}

calWr <- function(W, Wr.lambda, regularize) { # generalized
  if (regularize==1) { # no regularization

    Wr <- W

  } else if (regularize==2) {

    if (length(W)>1) {
      Ip <- diag(dim(W)[1])
      D <- diag(W)*diag(dim(W)[1])
      D.eig <- eigen(D)
      D.sqrt <- D.eig$vectors %*% diag(sqrt(D.eig$values)) %*% t(D.eig$vectors)
      Dinv <- try(solve(D.sqrt),silent=TRUE)
      if (class(Dinv)=="try-error") return(NA)
      R <- Dinv %*% W %*% Dinv
      R <- (1-Wr.lambda)*R + Wr.lambda*Ip
      Wr <- D.sqrt %*% R %*% D.sqrt
      colnames(Wr) <- colnames(W)
      rownames(Wr) <- rownames(W)
    } else {
      D <- W
      D.sqrt <- sqrt(W)
      Dinv <- 1/D.sqrt
      R <- Dinv * W * Dinv
      R <- (1-Wr.lambda)*R + Wr.lambda
      Wr <- D.sqrt * R * D.sqrt
    }

  } else if (regularize==3) {
    Ip <- diag(dim(W)[1])
    Wr <- Wr.lambda*W+(1-Wr.lambda)*Ip

  }
  Wr
}

# Choice 1 (Option 2)

calDfe <- function(dummyClass) { # calculate directly generalized one less than df.residual
  numLevels<-0
  totalLength<-0
  reducedClasses<-colnames(dummyClass)
  for (s in reducedClasses) {
    index<-which(dummyClass[,s]>0)
    if (length(index)>1) {
      numLevels<-numLevels+1
    }
  }
  #numLevels <- numLevels + 1
  totalLength <- dim(dummyClass)[1]
  return(totalLength - numLevels)
}

# multi label
calE.W <- function(trainFeatures,dummyClass) # calculate directly generalized
{
  reducedClasses<-colnames(dummyClass)
  if (is.null(dim(trainFeatures))) {
    W <- 0
  } else {
    W <- matrix(0,nrow=dim(trainFeatures)[2],ncol=dim(trainFeatures)[2],dimnames=list(colnames(trainFeatures),colnames(trainFeatures)))
  }

  numLevels<-0
  totalLength<-0
  for (s in reducedClasses) {
    index<-which(dummyClass[,s]>0)
    if (length(index)>1) {
      if (is.null(dim(trainFeatures))) {
        levelIData <- trainFeatures[index]
        covClass <- var(levelIData)*(length(index)-1)
      } else {
        levelIData <- trainFeatures[index,]
        covClass <- cov(levelIData)*(length(index)-1)
      }
      W <- W + covClass
      numLevels<-numLevels+1
    }
  }

  if (is.null(dim(trainFeatures))) {
    totalLength <- length(trainFeatures)
  } else {
    totalLength <- dim(trainFeatures)[1]
  }

  intercept
  #index <- 1:totalLength
  #levelIData <- trainFeatures[index,]
  #covClass <- cov(levelIData)*(length(index)-1)
  #W <- W + covClass
  #numLevels <- numLevels + 1

  W <- W / (totalLength - numLevels)

  return(W)
}

calH.B <- function(trainFeatures,dummyClass) {
  cm <- colMeans(trainFeatures)
  reducedClasses<-colnames(dummyClass)
  B <- matrix(0,nrow=dim(trainFeatures)[2],ncol=dim(trainFeatures)[2],dimnames=list(colnames(trainFeatures),colnames(trainFeatures)))
  numLevels<-0
  totalLength<-0
  for (s in reducedClasses) {
    index<-which(as.vector(dummyClass[,s])==TRUE)
    levelIData <- trainFeatures[index,]
    if (!is.null(dim(levelIData)[1])) {
      diff <- colMeans(levelIData)-cm
      covClass <- diff %*% t(diff)
      B <- B + length(index) * covClass
      numLevels<-numLevels+1
    }
  }
  B <- B / (numLevels - 1)
  B
}

# Choice 2 = Call Anova for mlm
# Choice 2.2 (Option 3)

linearHypothesis.ridgelm <- function(model, hypothesis.matrix, rhs=NULL, V, ...)
  linearHypothesis.mlm(model, hypothesis.matrix, rhs=NULL, V, ...)

linearHypothesis.mlm <- function(mod, hypothesis.matrix, rhs=NULL, V){
  test <- "Wilks"
  mod.coef <- mod$coefficients
  rownames(mod.coef)[1] <- "(Intercept)"
  B <- mod.coef # B = coefficients(mod)
  L <- if (is.null(dim(hypothesis.matrix))) t(hypothesis.matrix) else hypothesis.matrix
  #rank <- sum(eigen(SSPE, only.values=TRUE)$values >= sqrt(.Machine$double.eps))
  r <- ncol(B)
  rhs <- matrix(0, nrow(L), r)
  rownames(rhs) <- rownames(L)
  colnames(rhs) <- colnames(B)
  q <- NROW(L)
  SSPH <- t(L %*% B - rhs) %*% solve(L %*% V %*% t(L)) %*% (L %*% B - rhs) # V is Q in the note
  rval <- list(SSPH=SSPH, df=q, r=r, title="", test=test)
  class(rval) <- "linearHypothesis.mlm"
  rval
}

df.terms <- function(model, term, ...){
  UseMethod("df.terms")
}

df.terms.default <- function(model, term, ...){
  #if (is.aliased(model)) stop("Model has aliased term(s); df ambiguous.")
  if (!missing(term) && 1 == length(term)){
    assign <- model$assign
    which.term <- which(term == labels(terms(model)))
    if (0 == length(which.term)) stop(paste(term, "is not in the model."))
    sum(assign == which.term)
  }
  else {
    terms <- if (missing(term)) labels(terms(model)) else term
    result <- numeric(0)
    for (term in terms) result <- c(result, Recall(model, term))
    names(result) <- terms
    result
  }
}

calH.SSP.DF <- function(mod,idClassMat,allQ,classWeighted=FALSE) { # for mlm; calculated by Anova
  x <- mod$x #added
  if (classWeighted) {
    wts <- allQ[idClassMat[as.integer(rownames(x)),2]]
  } else {
    wts <- rep(1,nrow(x)) # wts <- rep(1,nrow(model.matrix(mod)))
  }
  # V = (X'WX)^{-1}; x is the Y in Y = XB + U
  V<-try(solve(wcrossprod(x, w=wts)),silent=TRUE) # V <- solve(wcrossprod(model.matrix(mod), w=wts))
  if (class(V)=="try-error") {
    if (is.null(mod$lambda)) return(NA)
    # B = (Xr'Xr)^(-1)X'Y
    V <- mrdivide(mod$coefficients,t(mod$x)%*%mod$y)
  }

  fac <- attr(mod$terms, "factors")
  B <- mod$coefficients # B = coefficients(mod)
  p <- nrow(B)  # B = coefficients(mod)
  I.p <- diag(p)
  assign <- mod$assign
  terms <- labels(terms(mod))
  n.terms <- length(terms)

  #assign <- array(0,dim=length(mod$assign))
  #for (i in 1:n.terms) assign[match(paste(terms[i],1,sep=""),colnames(mod$x))]<-1

  #SSPE <- calE.SSPE(mod,idClassMat,allQ)

  SSP <- as.list(rep(0, n.terms))
  df <- rep(0, n.terms)
  names(df) <- names(SSP) <- terms
  for (i in 1:n.terms){
    term<-terms[i]
    which.term <- which(term == terms)
    subs.term <- which(assign == which.term)
    relatives <- NULL
    subs.relatives <- integer(0)
    hyp.matrix.1 <- I.p[subs.relatives,,drop=FALSE]
    hyp.matrix.2 <- I.p[c(subs.relatives, subs.term),,drop=FALSE]
    SSP1 <- if (length(subs.relatives) == 0) 0 else linearHypothesis(mod, hypothesis.matrix=hyp.matrix.1, V=V)$SSPH
    SSP2 <- linearHypothesis(mod, hypothesis.matrix=hyp.matrix.2, V=V)$SSPH
    SSP[[i]]<-(SSP2 - SSP1)

    df[i]<- df.terms(mod, terms[i])
  }
  list(SSP=SSP,dfh=sum(df)-1)
}

calE.SSPE <- function(mod,idClassMat,allQ,classWeighted=FALSE) { # idClassMat, allQ; for mlm; calculated by Anova # equals calSSPE below
  if (classWeighted) {
    wts <- allQ[idClassMat[as.integer(rownames(mod$x)),2]]
  } else {
    wts <- rep(1, nrow(mod$x))
  }
  SSPE <- wcrossprod(mod$residuals, w=wts)
  SSPE
}

calH.SSPS <- function(trainFeatures,dummyClass,SSP) {
  if (proportionalH) {
    nfeatures<-ncol(trainFeatures)
    n<-dim(dummyClass)[1]
    nclass<-ncol(dummyClass)
    lsum<-0
#for (i in 1:nclass){ l1<-length(which(dummyClass[,i]>0)); l2<-length(which(dummyClass[,i]==0)); lsum<-lsum+min(l1,l2) }
#SSPS<-array(0,dim=c(nfeatures,nfeatures),dimnames=list(colnames(trainFeatures),colnames(trainFeatures)))
#for (i in 1:nclass){ l1<-length(which(dummyClass[,i]>0)); l2<-length(which(dummyClass[,i]==0)); SSPS<-SSPS+SSP[[i]]*min(l1,l2)/lsum }
    for (i in 1:nclass){ l1<-length(which(dummyClass[,i]>0)); l2<-length(which(dummyClass[,i]==0)); lsum<-lsum+l1 }
    SSPS<-array(0,dim=c(nfeatures,nfeatures),dimnames=list(colnames(trainFeatures),colnames(trainFeatures)))
    for (i in 1:nclass){ l1<-length(which(dummyClass[,i]>0)); l2<-length(which(dummyClass[,i]==0)); SSPS<-SSPS+SSP[[i]]*l1/lsum }
    SSPS*nclass
  } else {
    nfeatures<-ncol(trainFeatures)
    nclass<-ncol(dummyClass)
    SSPS<-array(0,dim=c(nfeatures,nfeatures),dimnames=list(colnames(trainFeatures),colnames(trainFeatures)))
    for (i in 1:nclass){ SSPS<-SSPS+SSP[[i]] }
    SSPS
  }
}

calSp <- function(mod,SSPE,dummyClass) { # for mlm; estimated by Anova & candisc
 Sp <- SSPE / mod$df.residual
 Sp
}

# Choice 2.3

# calH.SSPS <- function(trainFeatures,dummyClass,SSP) { # 較少
  # nfeatures<-ncol(trainFeatures)
  # nclass<-ncol(dummyClass)
  # SSPS<-array(0,dim=c(nfeatures,nfeatures),dimnames=list(colnames(trainFeatures),colnames(trainFeatures)))
  # for (i in 1:nclass){ SSPS<-SSPS+SSP[[i]] }
  # SSPS
# }

# calSp <- function(mod,SSPE,dummyClass) { # for mlm; estimated by Anova & candisc
  # Sp <- SSPE / mod$df.residual
  # nclass<-ncol(dummyClass)
  # Sp*nclass
# }

# Choice 3

# calSSPE <- function(mod) { # for lm; formula in SAS book # equals calSSPE above
  # return(t(mod$y)%*%mod$y - t(mod$coefficients)%*%t(mod$x)%*%mod$y)
# }

# calSSPT <- function(mod) { # for lm; formula in SAS book
  # return(t(mod$y)%*%mod$y)
# }

# calBt <- function(mod) { # for lm; formula in SAS book
  # Bt <- ((dim(mod$y)[1]-1) * calSSPT(mod) - mod$df.residual * calSSPE(mod)) / (mod$rank - 1)
  # Bt
# }

# Choice 4

# calEHr <- function(mod,lambda) {
  # x <- mod$x[,-1]
  # n <- dim(mod$y)[1]
  # Syy <- cor(mod$y)
  # Sxx <- cor(x)
  # Sxx_r <-(1-lambda)*Sxx+lambda*diag(dim(Sxx)[1])
  # Syx <- cor(mod$y,x)
  # Sxy <- cor(x,mod$y)
  # E <- (n-1)*(Syy - Syx %*% solve(Sxx_r) %*% Sxy)
  # H <- (n-1)*(Syx %*% solve(Sxx_r) %*% Sxy)
  # return(list(E=E,H=H))
# }

# calEH <- function(mod) { # formula by SAS; estimating Sp & Bt (between-class variance-covariance matrix)
  # x <- mod$x[,-1]
  # n <- dim(mod$y)[1]
  # Syy <- cor(mod$y)
  # Sxx <- cor(x)
  # Syx <- cor(mod$y,x)
  # Sxy <- cor(x,mod$y)
  # E <- (n-1)*(Syy - Syx %*% solve(Sxx) %*% Sxy)
  # H <- (n-1)*(Syx %*% solve(Sxx) %*% Sxy)
  # return(list(E=E,H=H))
# }

#system.time(R<-sapply(train_classes,function(x) { exp(-0.5*mahalanobis(validate_scores,trainCentroid[as.character(x),],train_Sp,tol=1e-20))})) # distance to each centroid
#    user  system elapsed
#920.737   2.624 923.830
#system.time(R2 <- do.call(rbind,as.list(by(trainCentroid, as.factor(train_classes), function(x) {exp(-0.5*mahalanobis(validate_scores,as.matrix(x),train_Sp,tol=1e-20))}))))
#R2 <- t(R2)
#   user  system elapsed
#917.802   1.432 919.702
#   user  system elapsed ... with matrix
#903.696   0.372 904.531
#ptm <- proc.time()
#R3 <- array(0,dim=c(dim(validate_scores)[1],length(train_classes)))
#for (i in 1:length(train_classes)) {
#  R3[,i] <- exp(-0.5*mahalanobis(validate_scores,trainCentroid[i,],train_Sp,tol=1e-20))
#}
#proc.time() - ptm
#   user  system elapsed
#916.569   3.636 920.674
#system.time(R4 <- apply(trainCentroid,1,function(x) exp(-0.5*mahalanobis(validate_scores,x,train_Sp,tol=1e-20))))
#   user  system elapsed
#920.310   1.288 922.068
   #user  system elapsed
#905.649   2.100 908.353

checkinverse <- function(m) class(try(solve(m),silent=T))=="matrix"
           
calDistanceAndPredict <- function(id, scores, trainSp, trainQ, trainCentroid, allClasses, predictStruct) {
  trainClasses <- as.integer(rownames(trainCentroid))

  system.time(R <- do.call(rbind,as.list(by(trainCentroid, as.factor(trainClasses), function(x) {exp(-0.5*maha(scores, as.matrix(x), trainSp))}))))
  #   user  system elapsed
  #316.408   0.000 316.570

  R<-t(R)
  colnames(R)<-trainClasses
  rownames(R)<-id

  WS<-(R %*% trainQ)

  P<-sapply(trainClasses,function(x) R[,as.character(x)]*trainQ[as.character(x)])
  P<-t(as.matrix(sapply(1:dim(P)[1],function(x) if (is.na(P[x,1])) {c(rep(NA,length(P[x,])))} else {P[x,]/WS[x]})))
  rownames(P)<-rownames(R)
  colnames(P)<-trainClasses

  P[is.na(P)]<-0
  pClass<-do.call(rbind,as.list(apply(P,1,function(x) { m<-which(x>=max(x,na.rm = TRUE)); if (length(m)==0 || length(m)==dim(P)[2]) list(NA) else names(m); })))
  if (length(which(is.na(pClass)))==length(pClass)) {
    print("return NA")
    return(NA)
  }

  if (length(predictStruct$pClass)!=0) {
    predictStruct$id <- c(predictStruct$id,id)
    ncols <- ncol(predictStruct$pClass) - ncol(pClass)
    if (ncols>0) {
      pClass<-cbind(pClass, array(NA, dim=c(dim(pClass)[1],ncols)))
    } else {
      predictStruct$pClass<-cbind(predictStruct$pClass, array(NA, dim=c(dim(predictStruct$pClass)[1],-ncols)))
    }
    predictStruct$pClass <- rbind(predictStruct$pClass,pClass)
  } else {
    predictStruct$id <- id
    predictStruct$pClass <- pClass
  }

  # pClass is a matrix rather than a vector

  dif <- setdiff(allClasses,trainClasses)
  if (length(dif)>0) {
    for (d in 1:length(dif)) {
      P <- cbind(P,matrix(NA,nrow=dim(P)[1]))
      colnames(P)[dim(P)[2]] <- dif[d]
    }
  }
  P <- P[, order(as.integer(colnames(P)))]
  predictStruct$prob <- rbind(predictStruct$prob,P)
  predictStruct$prob <- apply(predictStruct$prob,2,function(x) replace(x,which(is.na(x)),0))

  predictStruct
}

calNumAnnotations <- function(id, classes, idClassMat) {
  Q<-c()
  totalLength<-0
  for (s in classes) {
    itsct<-intersect(idClassMat[which(idClassMat[,2]==s),1],id)
    if (length(itsct)>0) {
      Q<-c(Q,length(itsct))
      names(Q)[length(Q)]<-as.character(s)
    }
  }
  Q<-Q[order(Q)]
  return(Q)
}

calQ <- function(id, classes, idClassMat) {
  Q<-c()
  totalLength<-0
  for (s in classes) {
    itsct<-intersect(idClassMat[which(idClassMat[,2]==s),1],id)
    if (length(itsct)>0) {
      Q<-c(Q,length(itsct))
      names(Q)[length(Q)]<-as.character(s)
      totalLength<-totalLength+length(itsct)
    }
  }
  Q<-Q/totalLength
  return(Q)
}

calCentroids <- function(id,features,classes,idClassMat) {
  if (is.null(dim(features))) {
    trainCentroid<-array(0,dim=c(length(classes)))
  } else {
    trainCentroid<-array(0,dim=c(length(classes),ncol(features)))
  }
  i<-1
  for (s in classes) {
    itsct<-intersect(idClassMat[which(idClassMat[,2]==s),1],id)
    if (is.null(dim(features))) {
      trainCentroid[i]<-mean(features[as.character(itsct)])
    } else {
      trainCentroid[i,]<-colMeans(features[as.character(itsct),])
    }
    i<-i+1
  }
  if (is.null(dim(features))) {
    names(trainCentroid)<-as.character(classes)
  } else {
    rownames(trainCentroid)<-as.character(classes)
    colnames(trainCentroid)<-colnames(features)
  }
  return(trainCentroid)
}

calSD <- function(id,features,classes,idClassMat) {
  SD<-array(0,dim=c(length(classes),ncol(features)))
  i<-1
  for (s in classes) {
    itsct<-intersect(idClassMat[which(idClassMat[,2]==s),1],id)
    SD[i,]<-apply(features[as.character(itsct),],2,SD)
    i<-i+1
  }
  rownames(SD)<-as.character(classes)
  colnames(SD)<-colnames(features)
  return(SD)
}

outputMeasures <- function(measures,sortBy2="Macro",filePrefixMeasure) {  
  #print("Mean")  
  evaluationArrMean=rbind(by(measures$evaluationArr[,sortBy2],measures$evaluationArr[,c("State","Measure")],function(x) mean(x)))
  #print(evaluationArrMean)
  saveRDS(evaluationArrMean,file=paste(filePrefixMeasure,".evaluationArrMean",sep=""))
  
  #print("OSE")    
  evaluationArrOSE=rbind(by(measures$evaluationArr[,sortBy2],measures$evaluationArr[,c("State","Measure")],function(x) sd(x)/sqrt(cvFold)))
  #print(evaluationArrOSE)
  saveRDS(evaluationArrOSE,file=paste(filePrefixMeasure,".evaluationArrOSE",sep=""))
   
  #print("Mean by class")
  evaluationArrByClassMean=by(measures$evaluationArrByClass[,sortBy2],measures$evaluationArrByClass[,c("State","Class","Measure")],function(x) mean(x))
  for (cl in dimnames(evaluationArrByClassMean)$Class) {
    #print(paste("Class",cl))
    #print(evaluationArrByClassMean[,cl,])
    saveRDS(evaluationArrByClassMean[,cl,],file=paste(filePrefixMeasure,".evaluationArrByClassMean.",cl,sep=""))
  }  
 
  #print("OSE by class")
  evaluationArrByClassOSE=by(measures$evaluationArrByClass[,sortBy2],measures$evaluationArrByClass[,c("State","Class","Measure")],function(x) sd(x)/sqrt(cvFold))
  for (cl in dimnames(evaluationArrByClassOSE)$Class) {
    #print(paste("Class",cl))
    #print(evaluationArrByClassOSE[,cl,])
    saveRDS(evaluationArrByClassOSE[,cl,],file=paste(filePrefixMeasure,".evaluationArrByClassOSE.",cl,sep=""))
  }
  
  #print("Mean of Class SD")
  evaluationSDArrMean=rbind(by(measures$evaluationSDArr[,sortBy2],measures$evaluationSDArr[,c("State","Measure")],function(x) mean(x)))
  #print(evaluationSDArrMean)
  saveRDS(evaluationSDArrMean,file=paste(filePrefixMeasure,".evaluationSDArrMean",sep=""))
  
  #print("OSE of Class SD")
  evaluationSDArrOSE=rbind(by(measures$evaluationSDArr[,sortBy2],measures$evaluationSDArr[,c("State","Measure")],function(x) sd(x)/sqrt(cvFold)))
  #print(evaluationSDArrOSE)
  saveRDS(evaluationSDArrOSE,file=paste(filePrefixMeasure,".evaluationSDArrOSE",sep=""))  
  
  return(list(evaluationArrMean=evaluationArrMean,evaluationArrOSE=evaluationArrOSE,evaluationArrByClassMean=evaluationArrByClassMean,evaluationArrByClassOSE=evaluationArrByClassOSE))
}

assignMeasures <- function(measures,evaluateStruct,idClassMatEvaluate,idFeatureMatEvaluate,allNumAnn,allQ,filePrefixMeasure,cv,sortBy=0,state,nclasses) {
  sink(paste(filePrefixMeasure,".",cvFold,"-",cv,".",state,".txt",sep=""))
  perf <- evaluateMultilabel(evaluateStruct=evaluateStruct,idClassMatEvaluate,idFeatureMatEvaluate,allNumAnn,allQ,filePrefixMeasure,output=TRUE,sortBy=sortBy)
  sink()
 
  # 改進
  n <- length(evaluationColnamesArr)*(cv-1)*3 + length(evaluationColnamesArr)*(state-1)
  measures$evaluationArr$Measure[(n+1):(n+length(evaluationColnamesArr))] <- evaluationColnamesArr
  measures$evaluationArr$CV[(n+1):(n+length(evaluationColnamesArr))] <- cv
  measures$evaluationArr$State[(n+1):(n+length(evaluationColnamesArr))] <- state
  measures$evaluationArr$Perf[(n+1):(n+length(evaluationColnamesArr))] <- perf$evaluationArr
  #measures$evaluationArr<-rbind(measures$evaluationArr,data.frame(Measure=as.factor(evaluationColnamesArr),CV=as.factor(cv),State=as.factor(state),perf$evaluationArr))

  n2 <- length(evaluationColnamesArrByClass)*(cv-1)*3*nclasses + length(evaluationColnamesArrByClass)*(state-1)*nclasses
  for (cl in 1:nclasses) {    
    
    # 改進    
    measures$evaluationArrByClass$Measure[(n2+1):(n2+length(evaluationColnamesArrByClass))] <- evaluationColnamesArrByClass
    measures$evaluationArrByClass$CV[(n2+1):(n2+length(evaluationColnamesArrByClass))] <- cv
    measures$evaluationArrByClass$State[(n2+1):(n2+length(evaluationColnamesArrByClass))] <- state
    measures$evaluationArrByClass$Perf[(n2+1):(n2+length(evaluationColnamesArrByClass))] <- perf$evaluationArrByClass[,,cl]
    #measures$evaluationArrByClass<-rbind(measures$evaluationArrByClass,data.frame(Class=as.factor(cl),Measure=as.factor(evaluationColnamesArrByClass),CV=as.factor(cv),State=as.factor(state),perf$evaluationArrByClass[,,cl]))        
    
    n2 <- n2 + length(evaluationColnamesArrByClass)
  }
  
  # 改進
  n3 <- length(evaluationColnamesArrByClass)*(cv-1)*3*nclasses + length(evaluationColnamesArrByClass)*(state-1)
  measures$evaluationSDArr$Measure[(n3+1):(n3+length(evaluationColnamesArrByClass))] <- evaluationColnamesArrByClass
  measures$evaluationSDArr$CV[(n3+1):(n3+length(evaluationColnamesArrByClass))] <- cv
  measures$evaluationSDArr$State[(n3+1):(n3+length(evaluationColnamesArrByClass))] <- state
  measures$evaluationSDArr$Perf[(n3+1):(n3+length(evaluationColnamesArrByClass))] <- perf$evaluationSDArr
  #measures$evaluationSDArr<-rbind(measures$evaluationSDArr,data.frame(Measure=as.factor(evaluationColnamesArrByClass),CV=as.factor(cv),State=as.factor(state),perf$evaluationSDArr))        
  
  return(measures)
}

assignThreeMeasures <- function(measures,trainStruct,validateStruct,testStruct,idClassMat,idFeatureMat,idClassMatTest,idFeatureMatTest,
                                allNumAnn,allQ,filePrefixMeasure,cv,sortBy=0,nclasses) {

  measures<-assignMeasures(measures=measures,
    evaluateStruct=trainStruct,
    idClassMatEvaluate=idClassMat,
    idFeatureMatEvaluate=idFeatureMat,
    allNumAnn,allQ,
    filePrefixMeasure=paste(filePrefixMeasure,".train",sep=""),
    sortBy=sortBy,cv=cv,state=1,nclasses=nclasses)
    
  measures<-assignMeasures(measures=measures,
    evaluateStruct=validateStruct,
    idClassMatEvaluate=idClassMat,
    idFeatureMatEvaluate=idFeatureMat,
    allNumAnn,allQ,
    filePrefixMeasure=filePrefixMeasure,
    sortBy=sortBy,cv=cv,state=2,nclasses=nclasses)
    
  measures<-assignMeasures(measures=measures,
    evaluateStruct=testStruct,
    idClassMatEvaluate=idClassMatTest,
    idFeatureMatEvaluate=idFeatureMatTest,
    allNumAnn,allQ,
    filePrefixMeasure=paste(filePrefixMeasure,".test",sep=""),
    sortBy=sortBy,cv=cv,state=3,nclasses=nclasses)
    
  measures
}

initialMeasures <- function(measures) {
  n <- length(evaluationColnamesArr)*cvFold*3
  n3 <- length(evaluationColnamesArrByClass)*cvFold*3
  n2 <- n3*nclasses
  
  measures$evaluationArr<-list()
  measures$evaluationArr$Measure<-character(n)
  measures$evaluationArr$CV<-numeric(n)
  measures$evaluationArr$State<-numeric(n)
  measures$evaluationArr$Perf<-matrix(nrow=n,ncol=length(evaluationColnamesArrByClass2))

  measures$evaluationArrByClass<-list()
  measures$evaluationArrByClass$Class<-character(n2)
  measures$evaluationArrByClass$Measure<-character(n2)
  measures$evaluationArrByClass$CV<-numeric(n2)
  measures$evaluationArrByClass$State<-numeric(n2)
  measures$evaluationArrByClass$Perf<-matrix(nrow=n2,ncol=length(evaluationColnamesArrByClass2))
  
  measures$evaluationSDArr<-list()
  measures$evaluationSDArr$Measure<-character(n3)
  measures$evaluationSDArr$CV<-numeric(n3)
  measures$evaluationSDArr$State<-numeric(n3)
  measures$evaluationSDArr$Perf<-matrix(nrow=n3,ncol=length(evaluationColnamesArrByClass2))
}

otherMethods <- function(filePrefixClass,idClassMat,idFeatureMat,idClassMatTest,idFeatureMatTest,allQ,allNumAnn,
                        methodVec=c("binary_svm","binary_lda","binary_logit","ml_nc","ml_lda")) { # "ml_svm"
                        
  allClasses<-sort(unique(idClassMat[,2]))
  nclasses<-length(allClasses)

  methodMeasures <- list()
  i <- 1 
  for (method in methodVec) {
    # 改進
    methodMeasures[[i]]<-list()
    initialMeasures(methodMeasures[[i]])
    i <- i + 1
  }
  testData<-readRDS(testData, file = paste(filePrefixClass,".testData.",1,".Rda",sep=""))

  cvEnd <- 1
  for (cv in cvFold:cvEnd) {
  
    trainData<-readRDS(file = paste(filePrefixClass,".trainData.",1,".",cvFold,"-",cv,".Rda",sep=""))              
    if (cvFold==1) {
      validateData<-trainData
    } else {
      validateData<-readRDS(file = paste(filePrefixClass,".validateData.",1,".",cvFold,"-",cv,".Rda",sep=""))
    }
    
    if (("binary_svm" %in% methodVec || "binary_lda" %in% methodVec || "binary_logit" %in% methodVec)) {
      myTrain <- data.frame(id=trainData$id,trainData$features,trainData$dummyClass)
      for (co in colnames(trainData$dummyClass)) {
        myTrain[,co]=as.factor(myTrain[,co])
      }
      myValidate <- data.frame(id=validateData$id,validateData$features)
      myTest <- data.frame(id=testData$id,testData$features)
    }
            
    ############ Binary SVM #############
    index<-which(methodVec=="binary_svm")
    if (length(index)>0) {
      trainStruct <- list(id=trainData$id, pClass=c(), prob=array(0,dim=c(length(trainData$id),nclasses),dimnames=list(trainData$id,allClasses)))
      validateStruct <- list(id=validateData$id, pClass=c(), prob=array(0,dim=c(length(validateData$id),nclasses),dimnames=list(validateData$id,allClasses)))
      testStruct <- list(id=testData$id, pClass=c(), prob=array(0,dim=c(length(testData$id),nclasses),dimnames=list(testData$id,allClasses)))
      
      for (sName in colnames(trainData$dummyClass)) {
        formulastr<-paste(sName,"~ v1")
        if (trainData$nfeatures>1) {
          formulastr <- paste(formulastr,paste(paste("+",paste("v",2:trainData$nfeatures,sep="")),collapse = ''))
        }
        mySVM <- svm(as.formula(formulastr), data = myTrain, probability=TRUE)
        
        s <- substr(sName,2,nchar(sName))  
        myPredict <- predict(mySVM, myTrain, probability=TRUE)
        trainStruct$prob[,s] <- attr(myPredict,"probabilities")[,1]
        myPredict <- predict(mySVM, myValidate, probability=TRUE)
        validateStruct$prob[,s] <- attr(myPredict,"probabilities")[,1]
        myPredict <- predict(mySVM, myTest, probability=TRUE)
        testStruct$prob[,s] <- attr(myPredict,"probabilities")[,1]
      }      

      filePrefixMeasure=paste(filePrefixClass,".binary_svm",sep="")
      methodMeasures[[index]] <- assignThreeMeasures(measures=methodMeasures[[index]],trainStruct,validateStruct,testStruct,
                                idClassMat,idFeatureMat,idClassMatTest,idFeatureMatTest,
                                allNumAnn,allQ,filePrefixMeasure,
                                cv,nclasses=nclasses)
    }            
    
    ############ Binary LDA #############
    
    # Venables, W. N. and Ripley, B. D. (2002) Modern Applied Statistics with S. Fourth edition. Springer.
    # Ripley, B. D. (1996) Pattern Recognition and Neural Networks. Cambridge University Press.
    
    index<-which(methodVec=="binary_lda")
    if (length(index)>0) {
      trainStruct <- list(id=trainData$id, pClass=c(), prob=array(0,dim=c(length(trainData$id),nclasses),dimnames=list(trainData$id,allClasses)))
      validateStruct <- list(id=validateData$id, pClass=c(), prob=array(0,dim=c(length(validateData$id),nclasses),dimnames=list(validateData$id,allClasses)))
      testStruct <- list(id=testData$id, pClass=c(), prob=array(0,dim=c(length(testData$id),nclasses),dimnames=list(testData$id,allClasses)))
      
      for (sName in colnames(trainData$dummyClass)) {
        formulastr<-paste(sName,"~ v1")
        if (trainData$nfeatures>1) {
          formulastr <- paste(formulastr,paste(paste("+",paste("v",2:trainData$nfeatures,sep="")),collapse = ''))
        }
        myLDA <- lda(as.formula(formulastr), data = myTrain)
        
        s <- substr(sName,2,nchar(sName))  
        trainStruct$prob[,s] <- predict(myLDA, myTrain)$posterior[,1]
        validateStruct$prob[,s] <- predict(myLDA, myValidate)$posterior[,1]
        testStruct$prob[,s] <- predict(myLDA, myTest)$posterior[,1]
      }      
      
      filePrefixMeasure=paste(filePrefixClass,".binary_lda",sep="")
      methodMeasures[[index]] <- assignThreeMeasures(measures=methodMeasures[[index]],trainStruct,validateStruct,testStruct,
                                idClassMat,idFeatureMat,idClassMatTest,idFeatureMatTest,
                                allNumAnn,allQ,filePrefixMeasure,
                                cv,nclasses=nclasses)
      
    }
    
    ############ Binary logistic regression #############
    # Venables, W. N. and Ripley, B. D. (2002) Modern Applied Statistics with S. New York: Springer.
    
    index<-which(methodVec=="binary_logit")
    if (length(index)>0) {        
      trainStruct <- list(id=trainData$id, pClass=c(), prob=array(0,dim=c(length(trainData$id),nclasses),dimnames=list(trainData$id,allClasses)))
      validateStruct <- list(id=validateData$id, pClass=c(), prob=array(0,dim=c(length(validateData$id),nclasses),dimnames=list(validateData$id,allClasses)))
      testStruct <- list(id=testData$id, pClass=c(), prob=array(0,dim=c(length(testData$id),nclasses),dimnames=list(testData$id,allClasses)))
      
      for (sName in colnames(trainData$dummyClass)) {
        formulastr<-paste(sName,"~ v1")
        if (trainData$nfeatures>1) {
          formulastr <- paste(formulastr,paste(paste("+",paste("v",2:trainData$nfeatures,sep="")),collapse = ''))
        }
        myLogit <- glm(as.formula(formulastr), data = myTrain, family = "binomial")
        
        s <- substr(sName,2,nchar(sName))        
        trainStruct$prob[,s] <- predict(myLogit, newdata = myTrain, type="response")
        validateStruct$prob[,s] <- predict(myLogit, newdata = myValidate, type="response")
        testStruct$prob[,s] <- predict(myLogit, newdata = myTest, type="response")
      }
      
      methodMeasures[[index]] <- assignThreeMeasures(measures=methodMeasures[[index]],trainStruct,validateStruct,testStruct,
                                idClassMat,idFeatureMat,idClassMatTest,idFeatureMatTest,
                                allNumAnn,allQ,filePrefixMeasure=paste(filePrefixClass,".binary_logit",sep=""),
                                cv,nclasses=nclasses)
      
    }
    
    ############ Multi-label LDA #############
    
    # Hua Wang, Chris Ding, and Heng Huang. Multi-label Linear Discriminant Analysis, ECCV 2010
    
    index<-which(methodVec=="ml_lda")
    if (length(index)>0) {      
      emptyProb<-matrix(NA, nrow = 0, ncol=nclasses)
  
      X<-trainData$features
      if (is.factor(trainData$dummyClass[,1])) {
        Y<-data.matrix(trainData$dummyClass)-1
      } else {
        Y<-trainData$dummyClass
      }
      nominator<-t(Y) %*% X
      denominator<-colSums(Y)
      Mk<-nominator/denominator
      M<-colSums(nominator)/sum(denominator)
      Xwav<-t(sapply(1:nrow(X),function(i) X[i,]-M))
      CY<-cosine(Y)
      Z<-t(sapply(1:nrow(X),function(i) Y[i,]%*%CY/sum(abs(Y[i,]))))
      W<-diag(colSums(Z))
      Sb<-t(Xwav) %*% Z %*% inv(W) %*% t(Z) %*% Xwav
      L<-diag(rowSums(Z))
      St<-t(Xwav) %*% L %*% Xwav
      Sw<-St-Sb
      ev<-eigen(ginv(Sw) %*% Sb)      
      
      Xnew <- X %*% ev$vectors
      idFeatureMat2 <- cbind(as.integer(rownames(Xnew)),Xnew)
      
      trainData2 <- createData(id=trainData$id,idClassMat,idFeatureMat2,dummyType=1)
      trainCentroids2 <- calCentroids(trainData2$id,Xnew,trainData2$class,idClassMat)
      Q2 <- calQ(trainData2$id,trainData2$class,idClassMat)
           
      trainStruct <- list(id=c(), pClass=c(), prob=emptyProb)
      trainStruct <- calDistanceAndPredict(
              id=trainData2$id,
              scores=as.matrix(Xnew),
              trainSp=diag(dim(Xnew)[2]),
              trainQ=Q2, trainCentroid=trainCentroids2,
              allClasses, predictStruct=trainStruct)
              
      Xnew <- validateData$features %*% ev$vectors
      idFeatureMatVali <- cbind(as.integer(rownames(Xnew)),Xnew)
      idFeatureMat2 <- rbind(idFeatureMat2,idFeatureMatVali)
      
      validateStruct <- list(id=c(), pClass=c(), prob=emptyProb)
      validateStruct <- calDistanceAndPredict(
              id=validateData$id,
              scores=as.matrix(Xnew),
              trainSp=diag(dim(Xnew)[2]),
              trainQ=Q2, trainCentroid=trainCentroids2,
              allClasses, predictStruct=validateStruct)

      Xnew <- testData$features %*% ev$vectors
      idFeatureMat3 <- cbind(as.integer(rownames(Xnew)),Xnew)
      
      testStruct <- list(id=c(), pClass=c(), prob=emptyProb)
      testStruct <- calDistanceAndPredict(
              id=testData$id,
              scores=as.matrix(Xnew),
              trainSp=diag(dim(Xnew)[2]),
              trainQ=Q2, trainCentroid=trainCentroids2,
              allClasses, predictStruct=testStruct)
      
      methodMeasures[[index]] <- assignThreeMeasures(measures=methodMeasures[[index]],trainStruct,validateStruct,testStruct,
                                idClassMat,idFeatureMat2,idClassMatTest,idFeatureMat3,
                                allNumAnn,allQ,filePrefixMeasure=paste(filePrefixClass,".ml_lda",sep=""),
                                cv,nclasses=nclasses)      
    }
      
    ############ Multi-label Nearest Centroid #############

    index<-which(methodVec=="ml_nc")
    if (length(index)>0) {
      emptyProb<-matrix(NA, nrow = 0, ncol=nclasses)
      Q <- readRDS(paste(filePrefixClass,".Q.",cvFold,"-",cv,sep=""))
      trainCentroid <- readRDS(paste(filePrefixClass,".trainCentroid.",cvFold,"-",cv,sep=""))
    
      trainStruct <- list(id=c(), pClass=c(), prob=emptyProb)
      trainStruct <- calDistanceAndPredict(
              id=trainData$id,
              scores=as.matrix(trainData$features),
              trainSp=diag(dim(trainData$features)[2]),
              trainQ=Q, trainCentroid,
              allClasses, predictStruct=trainStruct)

      validateStruct <- list(id=c(), pClass=c(), prob=emptyProb)
      validateStruct <- calDistanceAndPredict(
              id=validateData$id,
              scores=as.matrix(validateData$features),
              trainSp=diag(dim(validateData$features)[2]),
              trainQ=Q, trainCentroid,
              allClasses, predictStruct=validateStruct)
      
      testStruct <- list(id=c(), pClass=c(), prob=emptyProb)
      testStruct <- calDistanceAndPredict(
              id=testData$id,
              scores=as.matrix(testData$features),
              trainSp=diag(dim(testData$features)[2]),
              trainQ=Q, trainCentroid,
              allClasses, predictStruct=testStruct)          
      
      methodMeasures[[index]] <- assignThreeMeasures(measures=methodMeasures[[index]],trainStruct,validateStruct,testStruct,
                                idClassMat,idFeatureMat,idClassMatTest,idFeatureMatTest,
                                allNumAnn,allQ,filePrefixMeasure=paste(filePrefixClass,".ml_nc",sep=""),
                                cv,nclasses=nclasses)     
      
    }
    
    ############ Multi-label SVM ###########      
    index<-which(methodVec=="ml_svm")
    if (length(index)>0) {
      P<-read.csv("testStruct.P")
      colnames(P)<-sapply(colnames(P),function(x) substr(x,2,nchar(x)))
      P<-P[,-1]
      p<-P[,1]
    }    
  }
  
  for (i in 1:length(methodMeasures)) {
    om<-outputMeasures(measures=methodMeasures[[i]],filePrefixMeasure=filePrefixMeasure)
    methodMeasures[[i]]$evaluationArrMean=om$evaluationArrMean
    methodMeasures[[i]]$evaluationArrOSE=om$evaluationArrOSE
    methodMeasures[[i]]$evaluationArrByClassMean=om$evaluationArrByClassMean
    methodMeasures[[i]]$evaluationArrByClassOSE=om$evaluationArrByClassOSE
  }
  
  for (i in 1:length(methodMeasures)) {
    print(methodVec[i])
    print(methodMeasures[[i]]$evaluationArrMean)
  }
}

transformPredictAndEvaluate <- function(filePrefixClass,idClassMat,idFeatureMat,idClassMatTest,idFeatureMatTest,allQ,allNumAnn,
                        classWeightedVec,proportionalHVec, 
                        linearModelTypeVec,regularizeWTypeVec,regularizeWLambdaVec,
                        sortBy,sortBy2,
                        validateOnly,
                        otherMethodMeasures) {

  allClasses<-sort(unique(idClassMat[,2]))
  nclasses<-length(allClasses)
  emptyProb<-matrix(NA, nrow = 0, ncol=nclasses)
  allClasses<-as.character(allClasses)
  regularizeWLambdaVecOptimized <- regularizeWLambdaVec

  print(paste("validateOnly=",validateOnly))

  # currently assume single species
  mlcdaMeasures <- list()
  evaluationArrList <- list()
  n <- (length(evaluationColnamesArr)*length(classWeightedVecName)*length(proportionalHVecName)*length(regularizeWLambdaVecOptimized))        
  #for (sp in 1:nspecies) {
  for (dt in 1:length(dummyTypeVec)) {
    mlcdaMeasures[[dt]] <- list()  
    evaluationArrList[[dt]] <- list()
    for (lt in 1:length(linearModelTypeVec)) {
      mlcdaMeasures[[dt]][[lt]] <- list()  
      evaluationArrList[[dt]][[lt]] <- list()
      for (rt in 1:length(regularizeWTypeVec)) {
        evaluationArrList[[dt]][[lt]][[rt]]$Measure=character(n)
        evaluationArrList[[dt]][[lt]][[rt]]$CV=numeric(n)
        evaluationArrList[[dt]][[lt]][[rt]]$Weighted.class=character(n)
        evaluationArrList[[dt]][[lt]][[rt]]$Weighted.H=as.factor(n)
        evaluationArrList[[dt]][[lt]][[rt]]$Lambda=as.factor(n)
        evaluationArrList[[dt]][[lt]][[rt]]$evaluationArrList=numeric(n)
      }
    }
  }    
  #}
  
  testData<-readRDS(file = paste(filePrefixClass,".testData.Rda",sep=""))
  testids<-testData$id
  
  if (validateOnly) {
    trainvalidate.compare.test.state<-2
  } else {
    trainvalidate.compare.test.state<-1
  }
  #cv <- cvFold
  cvEnd <- 1
  #repeatValidate <- FALSE
  #while (TRUE) {
  for (cv in cvFold:cvEnd) {
    print(paste("cv=",cv))
    #print(paste("trainvalidate.compare.test.state=",trainvalidate.compare.test.state))

    filePrefixOptimized<-paste(filePrefixClass,".",sortBy,".",sortBy2,".",cvFold,"-",cv,sep="")
             
    for (dt in 1:length(dummyTypeVec)) {
      print(paste(cv,dummyTypeVec[dt],":"))
      
      trainData<-readRDS(file = paste(filePrefixClass,".trainData.",dummyTypeVec[dt],".",cvFold,"-",cv,".Rda",sep=""))              
      if (cvFold==1) {
        validateData<-trainData
      } else {
        validateData<-readRDS(file = paste(filePrefixClass,".validateData.",dummyTypeVec[dt],".",cvFold,"-",cv,".Rda",sep=""))
      }
    
      for (lt in 1:(length(linearModelTypeVec))) {
        print(paste(cv,dummyTypeVec[dt],linearModelTypeVec[lt],":"))

        #if (trainvalidate.compare.test.state==1) {
        lmMod<-readRDS(file = paste(filePrefixClass,".lm.",dummyTypeVec[dt],".",linearModelTypeVec[lt],".",cvFold,"-",cv,".Rda",sep=""))

        if (trainvalidate.compare.test.state==3) {
          classWeightedVecOptimized <- bestCVParams[which(bestCVParams[,"Y.type"]==dummyTypeVec[dt] && Regression==linearModelTypeVec[lt]),"cwv"]
          proportionalHVecOptimized <- bestCVParams[which(bestCVParams[,"Y.type"]==dummyTypeVec[dt] && Regression==linearModelTypeVec[lt]),"phv"]
        } else {
          classWeightedVecOptimized <- classWeightedVec
        }
        
        for (cwv in 1:length(classWeightedVecOptimized)) {
          classWeighted <- classWeightedVecOptimized[cwv]
          print(paste(cv,dummyTypeVec[dt],linearModelTypeVec[lt],classWeighted,":"))
          
          for (phv in 1:length(proportionalHVecOptimized)) {
            proportionalH <- proportionalHVecOptimized[phv]
            print(paste(cv,dummyTypeVec[dt],linearModelTypeVec[lt],classWeighted,proportionalH,":"))

            if (is.factor(trainData$dummyClass[,1])) {            
              dummyMatrix<-data.matrix(trainData$dummyClass)-1
            } else {
              dummyMatrix<-trainData$dummyClass
            }
            if (estimateE==3) {
              E <- calE.SSPE(lmMod,idClassMat,allQ,classWeighted=classWeighted)
            } else if (estimateE==2) {
              dfe <- calDfe(dummyMatrix)
              E <- calE.W(trainData$features,dummyMatrix) * dfe # checked
            } else if (estimateE==1) {
              E<-lmMod.anova$SSPE
              dfe<-lmMod.anova$error.df
            }
            if (estimateH==3) {
              SSPDF <- calH.SSP.DF(lmMod,idClassMat,allQ,classWeighted=classWeighted)
              if (any(is.na(SSPDF))) next
              H <- calH.SSPS(trainData$features,dummyMatrix,SSPDF$SSP,proportionalH=proportionalH)
              if (any(is.na(H))) next
            } else if (estimateH==2) {
              H <- calH.B(trainData$features,dummyMatrix) * (ncol(dummyMatrix)-1)  # checked
              if (any(is.na(H))) next
              dfh <- lmMod$rank - 1
            } else if (estimateH==1) {
              H<-lmMod.anova$SSP$g1
              dfh<-lmMod.anova$df
            }

            #resp<-try(Spr <- readRDS(paste(filePrefixClass,".Spr.",fileSuffixParam,sep="")),silent=TRUE)
            #if (class(resp)=="try-error") next

            #if (checkinverse(Spr)) {
            #}
          
            for (rt in 1:length(regularizeWTypeVec)) {
              print(paste(cv,dummyTypeVec[dt],linearModelTypeVec[lt],classWeighted,proportionalH,regularizeWTypeVec[rt],":"))

              if (trainvalidate.compare.test.state==3) {
                regularizeWLambdaVecOptimized <- bestCVParams[which(bestCVParams[,"Y.type"]==dummyTypeVec[dt] 
                  && Regression==linearModelTypeVec[lt] && regularizeWTypeVec[rt]),"la"]
              } else {
                regularizeWLambdaVecOptimized <- regularizeWLambdaVec
              }
        
              for (la in 1:length(regularizeWLambdaVecOptimized)) {
                print(paste(cv,dummyTypeVec[dt],linearModelTypeVec[lt],classWeighted,proportionalH,regularizeWTypeVec[rt],regularizeWLambdaVecOptimized[la],":"))

                if (regularizeWTypeVec[rt]==1) {
                  fileSuffixParam <- paste(cvFold,"-",cv,".",dummyTypeVec[dt],".",linearModelTypeVec[lt],".",classWeighted,".",proportionalH,".",regularizeWTypeVec[rt],sep="")
                } else {
                  fileSuffixParam <- paste(cvFold,"-",cv,".",dummyTypeVec[dt],".",linearModelTypeVec[lt],".",classWeighted,".",proportionalH,".",regularizeWTypeVec[rt],".",regularizeWLambdaVecOptimized[la],sep="")
                }

                if (trainvalidate.compare.test.state==1) {
                  Spr <- calWr(calSp(lmMod,E,dummyMatrix),Wr.lambda=regularizeWLambdaVecOptimized[la],regularize=regularizeWTypeVec[rt])
                  if (any(is.na(Spr))) next
                  saveRDS(Spr, paste(filePrefixClass,".Spr.",fileSuffixParam,sep=""))

                  Er <- lmMod$df.residual * Spr
                  Tm <- tdecomp(Er)
                  if (any(is.na(Tm)) || !checkinverse(Tm)) next
                  if (estimateH==3) {
                    myCDA <- candisc(mod=lmMod, E=Er, H=H, dfh=SSPDF$dfh)
                  } else {
                    myCDA <- candisc(mod=lmMod, E=Er, H=H, dfh=dfh)
                  }

                  trainIDs <- rownames(myCDA$scores)
                  idFeatureMat2 <- cbind(as.integer(rownames(myCDA$scores)),myCDA$scores)

                  trainData2 <- createData(id=trainIDs,idClassMat,idFeatureMat2,dummyType=dummyType)
                  #W <- calE.W(trainData2$features,trainData2$dummyClass)
                  Q2 <- calQ(trainData2$id,trainData2$class,idClassMat)
                  trainCentroid2 <- calCentroids(trainData2$id,trainData2$features,trainData2$class,idClassMat)

                  saveRDS(myCDA$coeffs.raw, paste(filePrefixClass,".coeffs.",fileSuffixParam,sep=""))
                  #saveRDS(W, paste(filePrefixClass,".W.",fileSuffixParam,sep=""))
                  saveRDS(Q2, paste(filePrefixClass,".Q2.",fileSuffixParam,sep=""))
                  saveRDS(trainCentroid2, paste(filePrefixClass,".trainCentroid2.",fileSuffixParam,sep=""))
                }

                coeffs <- readRDS(paste(filePrefixClass,".coeffs.",fileSuffixParam,sep=""))
                #W <- readRDS(paste(filePrefixClass,".W.",fileSuffixParam,sep=""))
                Q2 <- readRDS(paste(filePrefixClass,".Q2.",fileSuffixParam,sep=""))
                trainCentroid2 <- readRDS(paste(filePrefixClass,".trainCentroid2.",fileSuffixParam,sep=""))
                # if (cvFold==1) {
                  # trainSD2 <- calSD(trainData2$id,trainData2$features,trainData2$class,idClassMat)
                  # saveRDS(trainSD2,paste(filePrefixClass,".trainSD2.",fileSuffixParam,sep=""))
                # }

                Spr2 <- diag(dim(validateScores)[2])
                if (checkinverse(Spr2)) {
                  if (cvFold!=1 && trainvalidate.compare.test.state==3) {

                    # training
                    
                    trainStruct <- list(id=c(), pClass=c(), prob=emptyProb)
                    trainStruct <- calDistanceAndPredict(
                            id=trainData2$id,
                            validateScores=myCDA$scores,
                            trainSp=diag(dim(myCDA$scores)[2]),
                            trainQ=Q2, trainCentroid2,
                            allClasses, validateStruct=trainStruct, idClassMat)
                              
                    # validation
                     
                    validateFeaturesScaled <- scale(as.matrix(validateData$features), center=TRUE, scale=FALSE)              
                    validateScores <- (validateFeaturesScaled  %*% as.matrix(coeffs))
                    write.csv(validateScores, paste(filePrefixClass,".validateScores.",fileSuffixParam,sep=""))
                    idFeatureMat2 <- rbind(idFeatureMat2,cbind(as.integer(rownames(validateScores)),validateScores))
                    
                    validateStruct <- list(id=c(), pClass=c(), prob=emptyProb)
                    validateStruct <- calDistanceAndPredict(
                            id=validateData$id,
                            validateScores=validateScores,
                            trainSp=Spr2, trainQ=Q2, trainCentroid2,
                            allClasses, validateStruct=validateStruct, idClassMat)                    
                    
                    # testing
                
                    testFeaturesScaled <- scale(as.matrix(testData$features), center=TRUE, scale=FALSE)
                    testScores <- (testFeaturesScaled  %*% as.matrix(coeffs))
                    idFeatureMat3 <- cbind(as.integer(rownames(testScores)),testScores)

                    testStruct <- list(id=c(), pClass=c(), prob=emptyProb)
                    testStruct <- calDistanceAndPredict(
                            id=testData$id,
                            validateScores=testScores,
                            trainSp=diag(dim(testScores)[2]), trainQ=Q2, trainCentroid2,
                            allClasses, validateStruct=testStruct, idClassMatTest)           
                    
                    mlcdaMeasures[[dt]][[lt]][[rt]] <- assignThreeMeasures(measures=mlcdaMeasures[[dt]][[lt]][[rt]],trainStruct,validateStruct,testStruct,
                                              idClassMat,idFeatureMat2,idClassMatTest,idFeatureMat3,
                                              allNumAnn,allQ,filePrefixMeasure=paste(filePrefixClass,".ml_rcda",sep=""),
                                              cv,nclasses=nclasses)
                                              
                  } else if (trainvalidate.compare.test.state!=3) {
                    validateStruct <- list(id=c(), pClass=c(), prob=emptyProb)
                    validateStruct <- calDistanceAndPredict(
                            id=validateData$id,
                            validateScores=validateScores,
                            trainSp=Spr2, trainQ=Q2, trainCentroid2,
                            allClasses, validateStruct=validateStruct, idClassMat)
                    if (!is.na(validateStruct)) {
                    
                      perf <- evaluateMultilabel(evaluateStruct=validateStruct,idClassMat,idFeatureMat,allNumAnn,allQ,filePrefixEvaluate=paste(filePrefixClass,".",fileSuffixParam,sep=""),output=TRUE,sortBy=sortBy)
                      
                      # 改進
                      evaluationArrList[[dt]][[lt]][[rt]]<- data.frame(Measure=as.factor(evaluationColnamesArr),CV=as.factor(cv),Weighted.class=as.factor(classWeightedVecName[cwv]),Weighted.H=as.factor(proportionalHVecName[phv]),Lambda=as.factor(regularizeWLambdaVecOptimized[la]),perf$evaluationArrList))
                    }
                  }
                }
              }
            }
          }
        }
      }
    }
    print(paste("cv=",cv,"finished"))
  
    index<-which(colnames(evaluationArrList)=="CV")

    if (trainvalidate.compare.test.state==1) {
      saveRDS(evaluationArrList, paste(filePrefixClass,".evaluationArrList",sep=""))
      trainvalidate.compare.test.state <- 2
    } else if (trainvalidate.compare.test.state==3) {
      regularizeWLambdaVecOptimized <- regularizeWLambdaVec
      if (!validateOnly) {
        trainvalidate.compare.test.state <- 1
      } else {
        trainvalidate.compare.test.state <- 2
      }
      cv <- cv - 1
      if (cv==0) break
    }

    # compare validation results
    if (trainvalidate.compare.test.state==2) {
      evaluationArrList[[dt]][[lt]][[rt]] <- readRDS(paste(filePrefixClass,".evaluationArrList[[dt]][[lt]][[rt]]",sep=""))
      
      # Goal 1: Find best cwv, phv, la for this fold & for every sortBy-sortBy2 & for every dt-lt-rt      

      for (dt in 1:length(dummyTypeVec)) {
        for (lt in 1:(length(linearModelTypeVec))) {
          for (rt in 1:length(regularizeWTypeVec)) {
            for (sortBy in sortByVec) {
              for (sortBy2 in sortByVec2) {
                tmp<-do.call(rbind,by(
                  evaluationArrList[[dt]][[lt]][[rt]][which(evaluationArrList[[dt]][[lt]][[rt]][,"CV"]==cv),sortBy2],
                  evaluationArrList[which(evaluationArrList[[dt]][[lt]][[rt]][,"CV"]==cv),c("Weighted.class","Weighted.H","Lambda")],
                  max)
                )
                evaluationArrList[[dt]][[lt]][[rt]]<-evaluationArrList[[dt]][[lt]][[rt]][-which(evaluationArrList[[dt]][[lt]][[rt]][,"CV"]==cv),]
                ### ????? how to find out the index of max
                # 改進                
                bestCVParams<-rbind(bestCVParams,data.frame(Y.type=dummyTypeVec[dt],Regression=linearModelTypeVec[lt],Regularized.W=regularizeWTypeVec[rt],Measure=sortBy,Measure2=sortBy2,cwv=cwv,phv=phv,la=la))
              }
            }
          }
        } 
      }
      trainvalidate.compare.test.state <- 3
    }
  }
   
  # Goal 2: Pick the best sortBy-sortBy2 avaeraging over folds over every dt-lt-rt that wins over most other methods 
  #         Report that particular sortBy-sortBy2 for all dt-lt-rt
  for (dt in 1:length(dummyTypeVec)) {
    for (lt in 1:(length(linearModelTypeVec))) {
      for (rt in 1:length(regularizeWTypeVec)) {  
        for (sortBy in sortByVec) {
          for (sortBy2 in sortByVec2) {
            indexArr<-which(mlcdaMeasures[[dt]][[lt]][[rt]]$evaluationArr[,"Measure"]==sortBy)
            tmp<-do.call(rbind,by(
                          mlcdaMeasures[[dt]][[lt]][[rt]]$evaluationArr[indexArr,sortBy2],
                          mlcdaMeasures[[dt]][[lt]][[rt]]$evaluationArr[indexArr,"CV"],
                          mean)
            ### ????
            
          }
        }
      }
    }
  }
   
  sink(paste(filePrefixClass,".evaluation.txt",sep=""))
  for (lt in 1:length(linearModelTypeVec)) {
    for (rt in 1:length(regularizeWTypeVec)) {
      outputMeasures(mlcdaMeasures[[lt]][[rt]],sortBy2,filePrefixMeasure=paste(filePrefixClass,".ml_rcda.",linearModelTypeVec[lt],".",regularizeWTypeVec[rt],sep=""))   
    }
  }
  sink()
      
  # BUG: return the results of species 1 only
  species<-1
  return(mlcdaMeasures)
}

calEigen <- function(classWeightedVec, proportionalHVec, linearModelTypeVec, regularizeWTypeVec, regularizeWLambdaVec,idClassMat,allQ,filePrefixClass) {

  if (doCalEigen) {
    
    eigenMat <- array(0,dim=8,dimnames=c("cv","Y.type","regression","weighted.class","weighted.H","regularized.W","lambda","eigen"))
    for (cv in 1:cvFold) {
      print(paste("cv=",cv))

      for (dt in 1:length(dummyTypeVec)) {        
        trainData<-readRDS(file = paste(filePrefixClass,".trainData.",dummyTypeVec[dt],".",cvFold,"-",cv,".Rda",sep=""))
        
        for (lt in 1:(length(linearModelTypeVec))) {
          print(paste(cv,linearModelTypeVec[lt],":"))
          lmMod<-readRDS(file = paste(filePrefixClass,".lm.",dummyTypeVec[dt],".",linearModelTypeVec[lt],".",cvFold,"-",cv,".Rda",sep=""))

          for (cwv in 1:length(classWeightedVec)) {
            classWeighted <- classWeightedVec[cwv]
            print(paste("classWeighted=",classWeighted))
            for (phv in 1:length(proportionalHVec)) {
              proportionalH <- proportionalHVec[phv]
              print(paste("proportionalH=",proportionalH))
              
              if (is.factor(trainData$dummyClass[,1])) {            
                dummyMatrix<-data.matrix(trainData$dummyClass)-1
              } else {
                dummyMatrix<-trainData$dummyClass
              }        
              if (estimateE==3) {
                E <- calE.SSPE(lmMod,idClassMat,allQ,classWeighted)
              } else if (estimateE==2) {
                dfe <- calDfe(dummyMatrix)
                E <- calE.W(trainData$features,dummyMatrix) * dfe
              } else if (estimateE==1) {
                E <- lmMod.anova$SSPE
              }
              if (estimateH==3) {
                SSPDF <- calH.SSP.DF(lmMod,idClassMat,allQ,classWeighted)
                if (any(is.na(SSPDF))) next
                H <- calH.SSPS(trainData$features,dummyMatrix,SSPDF$SSP,proportionalH)
                if (any(is.na(H))) next
              } else if (estimateH==2) {
                H <- calH.B(trainData$features,dummyMatrix) * (ncol(dummyMatrix)-1)
                if (any(is.na(H))) next
              } else if (estimateH==1) {
                H <- lmMod.anova$SSP$g1
              }

              for (rt in 1:length(regularizeWTypeVec)) {
                print(paste(cv,linearModelTypeVec[lt],regularizeWTypeVec[rt],":"))

                if (regularizeWTypeVec[rt]==1) {
                  end <- 1
                }
                else {
                  end <- length(regularizeWLambdaVec)
                }

                for (la in 1:end) {
                  print(paste(cv,linearModelTypeVec[lt],regularizeWTypeVec[rt],regularizeWLambdaVec[la],":"))
                  Spr <- calWr(calSp(lmMod,E,dummyMatrix),Wr.lambda=regularizeWLambdaVec[la],regularize=regularizeWTypeVec[rt])
                  if (any(is.na(Spr))) next
                  Er <- lmMod$df.residual * Spr
                  Tm <- tdecomp(Er)
                  if (!any(is.na(Tm)) && checkinverse(Tm)) {
                    eInv <- solve(Tm)
                    eHe <- t(eInv) %*% H %*% eInv
                    dc <- eigen(eHe, symmetric=TRUE)
                    # 改進                    
                    eigenMat<-rbind(eigenMat,c(cv,dummyTypeVecName[dt],linearModelTypeVecName[lt],classWeightedVecName[cwv],proportionalHVecName[phv],regularizeWTypeVecName[rt],regularizeWLambdaVec[la],(dc$values[1]/(1+dc$values[1]))))
                  }
                }
              }
            }          
          }
        }
      }
    }
    saveRDS(eigenMat, file = paste(filePrefixClass,".eigenMat.Rda",sep=""))
  } else {
    eigenMat <- readRDS(file = paste(filePrefixClass,".eigenMat.Rda",sep=""))
  }
  
  if (searhByEigen) {
    eigenMat2 <- rbind(by(eigenMat[,"eigen"],eigenMat[,"cv"],function(x) mean(x)))
    maxEigen=rbind(by(eigenMat2[,"eigen"],eigenMat[,c("Y.type","regression","weighted.class","weighted.H","regularized.W")],function(x) max(x)))
    
    # output
    # filter
  }

  return(list(dummyTypeVec=dummyTypeVec,linearModelTypeVec=linearModelTypeVec,classWeightedVec=classWeightedVec,proportionalHVec=proportionalHVec,regularizeWTypeVec=regularizeWTypeVec,regularizeWLambdaVec=regularizeWLambdaVec))
}

linearModel <- function(trainData,linearModelType,ridgeRegressionLambdaVec,validateData) { #,error.mat,maxLambda.vec) {
  myTrainData <- data.frame(trainData$id,trainData$features,trainData$dummyClass)
  if (is.null(dim(trainData$dummyClass))) {
    colnames(myTrainData)[dim(trainData$features)[2]+2]<-"g1"
  }
  if (is.null(dim(trainData$features))) {
    colnames(myTrainData)[2]<-"v1"
  }
  if (linearModelType==1) {

    lmMod <- lm(as.formula(trainData$formulastr),data=myTrainData, x=TRUE, y=TRUE, model=TRUE) # trainData
    lmModAnova <- try(Anova(lmMod,test="Wilks",type="2"),silent=TRUE)
    lmPred<-predict(lmMod)

    myValidateData <- data.frame(validateData$id,validateData$features,validateData$dummyClass)
    lmPred.vali <- predict(lmMod,myValidateData)

    return(list(lmMod=lmMod,lmModAnova=lmModAnova,error=mean((lmPred-lmMod$y)^2),error.vali=mean((lmPred.vali - validateData$features)^2)))

  } else if (linearModelType==2) {

    # load customed "lm.ridge"
    #Terms <- lmMod$term
    #Y <- lmMod$y
    #y <- Y
    #trainFeatures <- lmMod$x
    #x <- trainFeatures
    #xlevels <- lmMod$xlevels
    #assign <- lmMod$assign
    #offset <- NULL
    #Inter <- 1
    #n <- nrow(lmMod$x); p <- ncol(lmMod$x)

    lmMod <- lm.ridge(as.formula(trainData$formulastr),data=myTrainData, lambda=ridgeRegressionLambdaVec) # trainData$data
    maxLambda <- lmMod$lambda[which.min(lmMod$GCV)]

    #lmPred <- scale(lmMod$x[,2:ncol(lmMod$x)],scale = lmMod$scales)%*% lmMod$coef + lmMod$ym
    lmPred <- lmMod$fitted.values     # = lmMod$x %*% lmMod$coefficients

    missing <- setdiff(colnames(validateData$dummyClass),rownames(lmMod$coefficient)[-1])
    na.matrix <- matrix(NA,length(missing),dim(lmMod$coefficients)[2])
    rownames(na.matrix) <- missing
    coeff <- rbind(lmMod$coefficients,na.matrix)
    coeff <- coeff[order(rownames(coeff)),]

    na.index <- which(is.na(coeff[,1]))
    if (length(na.index)!=0) {
      lmPred.vali <- cbind(1,validateData$dummyClass)[,-na.index] %*% coeff[-na.index,]
    } else {
      lmPred.vali <- cbind(1,validateData$dummyClass) %*% coeff
    }

    return(list(lmMod=lmMod,maxLambda=maxLambda,error=mean((lmPred-lmMod$y)^2),error.vali=mean((lmPred.vali - validateData$features)^2)))

  } #else if (linearModelType==3) { # can't proceed
  #  class.f <- as.factor(lm.data.cv$class)
  #  dummies = model.matrix(~class.f)
  # m.lasso <- lars(as.matrix(sampled_features),dummies,trace=TRUE)
  #  jpeg(filePrefixClass,".lasso.",sep="")
  #  plot(m.lasso)
  #  dev.off()
  #  r <- cv.lars(train_classs,sampled_features)
  #  bestfraction <- r$index[which.min(r$cv)]
    #coef.lasso <- predict(m.lasso,train_classs,s=bestfraction,type="coefficient",mode="fraction")
    #coef.lasso
  #  sampled_features.pred.lasso <- predict(m.lasso,train_classs,s=bestfraction,type="fit",mode="fraction")$fit
  #  error.mat[lt]<-mean((sampled_features.pred.lasso - sampled_features)^2)
  #}
}

removeClassSizeOne<-function(idClassMat2) {
  classes<-idClassMat2[,2]
  classSize<-do.call(rbind,as.list(by(classes,as.factor(classes),function(x) length(x))))
  excludedClasses<-as.integer(rownames(which(classSize==1,arr.ind=TRUE)))
  #print(paste("excludedClasses=",excludedClasses))

  reducedClasses <- sort(setdiff(unique(classes),excludedClasses))
  idClassMatReduced <- idClassMat2[idClassMat2[,2] %in% reducedClasses,]
  #print(paste("excluded_proteins",idClassMat2[!idClassMat2[,2] %in% reducedClasses,1]))

  return(idClassMatReduced)
}

createData<-function(id,idClassMat,idFeatureMat,dummyType=1,removeSizeOne=FALSE,createDummy=TRUE) {
  idClassMat2 <- idClassMat[idClassMat[,1] %in% id,]

  # for testing - single label
  #idClassMat2<-do.call("rbind", as.list(by(idClassMat2,idClassMat2[,1],function(x) x[1,])))

  if (length(idClassMat2)>2) {
    if (removeSizeOne) {
      i<-1
      #print(i)
      idClassMatReduced<-removeClassSizeOne(idClassMat2)
      i<-i+1
      while (dim(idClassMatReduced)[1]!=dim(idClassMat2)[1]) {
        #print(i)
        idClassMat2<-idClassMatReduced
        idClassMatReduced<-removeClassSizeOne(idClassMat2)
        i<-i+1
      }
    } else {
      idClassMatReduced <- idClassMat2
    }

    reducedIDs <- unique(idClassMatReduced[,1])
    reducedClasses <- sort(unique(idClassMatReduced[,2]))

    # for Iris testing
    #reducedClasses <- reducedClasses[reducedClasses!=1]
  } else {
    idClassMatReduced <- idClassMat2
    reducedIDs <- as.integer(idClassMat2[1])
    reducedClasses <- as.integer(idClassMat2[2])
  }

  if (!is.null(idFeatureMat)) {
    reducedFeatures <- idFeatureMat[as.character(reducedIDs),2:ncol(idFeatureMat)]
    if (!is.null(dim(reducedFeatures))) {
      colnames(reducedFeatures) <- paste("v",1:ncol(reducedFeatures),sep="")
    } else {
      names(reducedFeatures) <- paste("v",1:length(reducedFeatures),sep="")
    }
    # data <- data.frame(id=reducedIDs, reducedFeatures)
  } else {
    # data <- data.frame(id=reducedIDs)
  }

  reducedClassNames <- paste("g",reducedClasses,sep="")
  reducedDummyClass <- NULL
  if (createDummy) {  
    reducedDummyClass <- array(0, dim=c(length(reducedIDs),length(reducedClasses)), dimnames=list(as.character(reducedIDs),reducedClassNames))
    if (length(reducedDummyClass)==1) {
      reducedDummyClass[1,1] <- 1
    } else {
      for (s in reducedClasses) {
        reducedDummyClass[as.character(idClassMatReduced[which(idClassMatReduced[,2]==s),1]),paste("g",s,sep="")]<-1
      }
    }
         
    if (missing(dummyType) || dummyType==1) { # binary
      reducedDummyClass<-data.frame(reducedDummyClass)
      for (s in reducedClasses) {
        reducedDummyClass[,paste("g",s,sep="")]<-as.factor(reducedDummyClass[,paste("g",s,sep="")])
      }
    } else if (dummyType==4 || dummyType==5) { # cross entropy & normalized    
      X<-reducedFeatures
      Y<-reducedDummyClass
      nominator<-t(Y) %*% X
      denominator<-colSums(Y)
      Mk<-nominator/denominator
      M<-colSums(nominator)/sum(denominator)
      Xwav<-t(sapply(1:nrow(X),function(i) X[i,]-M))
      CY<-cosine(Y)
      
      if (dummyType==4) {
        Z<-t(sapply(1:nrow(X),function(i) Y[i,]%*%CY/sum(abs(Y[i,]))))
      } else if (dummyType==5) {
        Z<-t(sapply(1:nrow(X),function(i) Y[i,]%*%CY/sum(abs(Y[i,]%*%CY))))
      }
      rownames(Z)<-rownames(reducedDummyClass)
      colnames(Z)<-colnames(reducedDummyClass)
      
      reducedDummyClass<-Z
    } else if (dummyType==6) {
      if (length(reducedDummyClass)==1) {
        reducedDummyClass <- idClassMatReduced[2]
      } else {
        reducedDummyClass <- idClassMatReduced[as.character(reducedIDs),2]
      }
    }
  }

  if (!is.null(idFeatureMat)) {

    if (is.null(dim(reducedFeatures))) {
      formulastr <- paste("cbind(",
                    paste(names(reducedFeatures), collapse=","),
                    ") ~ ")
      nfeatures<-length(reducedFeatures)
    } else {
      formulastr <- paste("cbind(",
                    paste(colnames(reducedFeatures), collapse=","),
                    ") ~ ")
      nfeatures<-ncol(reducedFeatures)
    }

    if (dummyType!=6) {
      formulastr<-paste(formulastr,reducedClassNames[1])
      if (length(reducedClassNames)>1) {
        for (i in 2:length(reducedClassNames)) {
          formulastr<-paste(formulastr,"+ ",reducedClassNames[i])
        }
      }
    } else {
      formulastr<-paste(formulastr,"g1")
    }
  } else {
    formulastr <- ""
    nfeatures<-0
  }

  if (length(reducedIDs)==1) {
    return(list(id=reducedIDs, features=t(reducedFeatures),dummyClass=reducedDummyClass,formulastr=formulastr,class=reducedClasses,nfeatures=nfeatures))
  } else {
    return(list(id=reducedIDs, features=reducedFeatures,dummyClass=reducedDummyClass,formulastr=formulastr,class=reducedClasses,nfeatures=nfeatures))
  }
}

dataToSVM <- function(createddata,filename) {
  multilabel.set<-apply(createddata$dummyClass,1,function(x) sapply(names(x)[which(x==1)],function(y) as.integer(substr(y,2,nchar(y)))))
  multilabel.set.str<-unlist(lapply(multilabel.set,function(z) paste(z,collapse=",")))
  feature.str<-paste(paste(1:ncol(createddata$features),":",createddata$features[1,],sep=""),collapse=" ")    
  fileConn<-file(filename)
  writeLines(paste(multilabel.set.str,feature.str,sep="  "), fileConn)
  close(fileConn)
}

createSaveData <- function(idClassMat,idFeatureMat,filePrefixClass,idClassMatTest,idFeatureMatTest) {
  # order data
  allIDs <- sort(unique(idClassMat[,1]))  # 3680

  if (cvFold==1) {
    sampled1 <- c()
    if (nspecies>=2) {
      sampled2 <- c()
    }
  } else {
    if (!useAncestor) {
      sampled1<-readRDS(file = paste(commonDir,dataNames,".sampled1.Rda",sep=""))
      if (nspecies>=2) {
        sampled2<-readRDS(file = paste(commonDir,dataNames,".sampled2.Rda",sep=""))
      }
    } else {
      sampled1<-readRDS(file = paste(commonDir,dataNames,".sampled1.ancestor.Rda",sep=""))
      if (nspecies>=2) {
        sampled2<-readRDS(file = paste(commonDir,dataNames,".sampled2.ancestor.Rda",sep=""))
      }
    }
  }
  
  for (cv in 1:cvFold) {
    print(paste("cv=",cv))
    if (cvFold==1) {
      sampled <- c()
    } else {
      sampled <- sampled1[(floor(length(sampled1)/cvFold*(cv-1))+1):(if (cv==cvFold) length(sampled1) else floor(length(sampled1)/cvFold*cv))]
      if (nspecies>=2) {
        sampled <- c(sampled,sampled2[(ceiling(length(sampled2)/cvFold*(cv-1))+1):(if (cv==cvFold) length(sampled2) else ceiling(length(sampled2)/cvFold*cv))])
      }
      sampled <- sort(sampled) # 1227
    }
    for (dt in 1:length(dummyTypeVec)) {
      validateData<-createData(sampled,idClassMat,idFeatureMat,dummyType=dummyTypeVec[dt])
      saveRDS(validateData, file = paste(filePrefixClass,".validateData.",dummyTypeVec[dt],".",cvFold,"-",cv,".Rda",sep=""))
     
      trainData<-createData(id=setdiff(allIDs,sampled),idClassMat,idFeatureMat,dummyType=dummyTypeVec[dt])
      saveRDS(trainData, file = paste(filePrefixClass,".trainData.",dummyTypeVec[dt],".",cvFold,"-",cv,".Rda",sep=""))
      
      trainValidateData<-createData(id=allIDs,idClassMat,idFeatureMat,dummyType=dummyTypeVec[dt])
      dataToSVM(createddata=trainValidateData,filename=paste(filePrefixClass,".train.validate.",dummyTypeVec[dt],".",cvFold,"-",cv,".svm",sep=""))      
      
      if (cv==1) {
        testids <- sort(unique(idClassMatTest[,1])) 
        testData<-createData(id=testids,idClassMatTest,idFeatureMatTest,dummyType=dummyTypeVec[dt])  
        saveRDS(testData, file = paste(filePrefixClass,".testData.",dummyTypeVec[dt],".Rda",sep=""))
        dataToSVM(createddata=testData,filename=paste(filePrefixClass,".test.",dummyTypeVec[dt],".svm",sep=""))      
      }
    }

    Q <- calQ(trainData$id,trainData$class,idClassMat)
    saveRDS(Q, paste(filePrefixClass,".Q.",cvFold,"-",cv,sep=""))
    
    trainCentroid <- calCentroids(trainData$id,trainData$features,trainData$class,idClassMat)    
    saveRDS(trainCentroid, paste(filePrefixClass,".trainCentroid.",cvFold,"-",cv,sep=""))
  }  
}

fitLinearModels <- function(idClassMat,idFeatureMat,filePrefixClass)
{
  # measurements
  n <- cvFold*length(dummyTypeVec)*length(linearModelTypeVec)*2
  # 改進
  errorMat <- data.frame(cv=numeric(n), Y.type=character(n), regression=character(n), set=character(n), error=numeric(n))
  maxLambdaArr <- array(0,dim=c(cv,length(dummyTypeVec)))

  for (cv in 1:cvFold) {
    print(paste("cv=",cv))
    for (dt in 1:length(dummyTypeVec)) {
      validateData<-readRDS( file = paste(filePrefixClass,".validateData.",dummyTypeVec[dt],".",cvFold,"-",cv,".Rda",sep=""))    
      trainData<-readRDS(file = paste(filePrefixClass,".trainData.",dummyTypeVec[dt],".",cvFold,"-",cv,".Rda",sep=""))       
      for (lt in 1:(length(linearModelTypeVec))) {
        print(paste("linearModelTypeVec[lt]=",linearModelTypeVec[lt]))
        lmOut<-linearModel(trainData,linearModelType=linearModelTypeVec[lt],ridgeRegressionLambdaVec,validateData)
        # 改進        
        errorMat<-rbind(errorMat,c(cv,dummyTypeVecName[dt],linearModelTypeVecName[lt],"training",lmOut$error))
        errorMat<-rbind(errorMat,c(cv,dummyTypeVecName[dt],linearModelTypeVecName[lt],"validation",lmOut$error.vali))
        saveRDS(lmOut$lmMod, file = paste(filePrefixClass,".lm.",dummyTypeVec[dt],".",linearModelTypeVec[lt],".",cvFold,"-",cv,".Rda",sep=""))
        if (linearModelTypeVec[lt]==1) {
          saveRDS(lmOut$lmModAnova, file = paste(filePrefixClass,".lm.anova.",dummyTypeVec[dt],".",linearModelTypeVec[lt],cvFold,"-",cv,".Rda",sep=""))
        } else {
          maxLambdaArr[cv,dt]<-lmOut$maxLambda
        }
      }
    }
  }
  for (i in 1:(ncol(errorMat)-1)) {
    errorMat[,i]<-as.factor(errorMat[,i])
  }
  print(paste("errorMat=",errorMat))
  print(paste("maxLambdaArr=",maxLambdaArr))
  saveRDS(errorMat, file = paste(filePrefixClass,".lm.error.Rda",sep="")) # can be plotted against # dimensions; paired with validation mse
  saveRDS(maxLambdaArr, file = paste(filePrefixClass,".maxLambdaArr.Rda",sep=""))
  
  # use last fold to plot
  unit<-ceiling(length(id)/100)
  # randIDs<-sample(id)
  randIDs<-id
  numGaps<-floor(length(randIDs)/unit)  
  errorMat2 <- array(0,dim=5,dimnames=c("training.size","Y.type","regression","set","error"))

  for (i in 0:(numGaps-1)) {
    trainSize=i*unit+1
    print(paste("trainSize=",trainSize))
    for (dt in 1:length(dummyTypeVec)) {
      trainData<-createData(randIDs[1:trainSize],idClassMat,idFeatureMat,dummyType=dummyTypeVec[dt])
      for (lt in 1:(length(linearModelTypeVec))) {
        print(paste("linearModelTypeVec[lt]=",linearModelTypeVec[lt]))
        lmOut<-linearModel(trainData,linearModelType=linearModelTypeVec[lt],ridgeRegressionLambdaVec,validateData)        
        # 改進        
        errorMat2<-rbind(errorMat2,c(trainSize,dummyTypeVecName[dt],linearModelTypeVecName[lt],"training",lmOut$error))
        errorMat2<-rbind(errorMat2,c(trainSize,dummyTypeVecName[dt],linearModelTypeVecName[lt],"validation",lmOut$error.vali))        
      }
    }
  }
  for (i in 2:(ncol(errorMat2)-1)) {
    errorMat2[,i]<-as.factor(errorMat2[,i])
  }  
  saveRDS(errorMat2, file = paste(filePrefixClass,".lm.error.trainSize.Rda",sep="")) # can be plotted against train data size; paired with validation mse
}

errorplot.dimension <- function(filePrefixClassArr,xArr) # xArr might hold dimension
{
  errorMat3 <- array(0,dim=6,dimnames=c("dimension","cv","Y.type","regression","set","error"))
  i<-1
  for (filePrefixClass in filePrefixClassArr) {
    errorMat<-readRDS(file = paste(filePrefixClass,".lm.error.Rda",sep=""))
    # 改進    
    errorMat3<-rbind(errorMat3,cbind(xArr[i],errorMat))
    i<-i+1
  }  
  errorMat4=rbind(by(errorMat3[,"error"],errorMat3[,"cv"],function(x) mean(x)))

  g <- ggplot(errorMat4, aes(dimension, error, color=set)) + geom_point()
  g + facet_grid(Y.type ~ regression)
}

errorplot.trainSize.dimension <- function(filePrefixClassArr,xArr) # xArr might hold dimension
{
  errorMat3 <- array(0,dim=6,dimnames=c("dimension","training.size","Y.type","regression","set","error"))
  i<-1
  for (filePrefixClass in filePrefixClassArr) {
    errorMat2<-readRDS(file = paste(filePrefixClass,".lm.error.trainSize.Rda",sep=""))    
    # 改進    
    errorMat3<-rbind(errorMat3,cbind(xArr[i],errorMat2))
    
    # combined plot
    g <- ggplot(errorMat3, aes(training.size, error, color=set)) + geom_point()
    g + facet_grid(Y.type ~ regression)
    
    i<-i+1
  }

}

evaluateMultilabel <- function(evaluateStruct,idClassMat,idFeatureMat,allNumAnn,allQ,filePrefixEvaluate="",byLevel=FALSE,output=FALSE,sortBy=0) {  # evaluateStructMat
  evaluateStructSorted <- evaluateStruct
  evaluateStructSorted$prob <- evaluateStructSorted$prob[order(evaluateStructSorted$id),]
  evaluateStructSorted$pClass <- evaluateStructSorted$pClass[order(evaluateStructSorted$id)]
  evaluateStructSorted$id <- sort(evaluateStructSorted$id)
  nC <- ncol(evaluateStructSorted$prob) # number of classes

  perArr <- array(0,dim=c(length(evaluationColnamesArrByClass2),1),dimnames=list(evaluationColnamesArrByClass2,"total"))
  perArrByClass <- array(0,dim=c(length(evaluationColnamesArrByClass2)-1,nC),dimnames=list(evaluationColnamesArrByClass2[-which(evaluationColnamesArrByClass2=="Micro")],colnames(evaluateStructSorted$prob)))
  perSDArr <- array(0,dim=c(length(evaluationColnamesArrByClass2)-1,1),dimnames=list(evaluationColnamesArrByClass2[-which(evaluationColnamesArrByClass2=="Micro")],"total"))

  evaluationArr <- array(0,dim=c(length(evaluationColnamesArr),length(evaluationColnamesArrByClass2)),dimnames=list(evaluationColnamesArr,evaluationColnamesArrByClass2))
  evaluationArrByClass <- array(0,dim=c(length(evaluationColnamesArrByClass),length(evaluationColnamesArrByClass2)-1,nC),dimnames=list(evaluationColnamesArrByClass,evaluationColnamesArrByClass2[-which(evaluationColnamesArrByClass2=="Micro")],colnames(evaluateStructSorted$prob)))
  evaluationSDArr <- array(0,dim=c(length(evaluationColnamesArrByClass),length(evaluationColnamesArrByClass2)-1),dimnames=list(evaluationColnamesArrByClass,evaluationColnamesArrByClass2[-which(evaluationColnamesArrByClass2=="Micro")]))
  
  if (byLevel) {  
    evaluationArrByLevel <- array(0,dim=c(length(evaluationColnamesArr),length(evaluationColnamesArrByClass2),length(classFrequencyArr)-1),dimnames=list(evaluationColnamesArr,evaluationColnamesArrByClass2,paste(classFrequencyArr[-length(classFrequencyArr)]+1,"-",classFrequencyArr[-1],sep="")))
    evaluationArrByClassByLevel <- array(0,dim=c(length(evaluationColnamesArrByClass),length(evaluationColnamesArrByClass2)-1,nC,length(classFrequencyArr)-1),dimnames=list(evaluationColnamesArrByClass,evaluationColnamesArrByClass2[-which(evaluationColnamesArrByClass2=="Micro")],colnames(evaluateStructSorted$prob),paste(classFrequencyArr[-length(classFrequencyArr)]+1,"-",classFrequencyArr[-1],sep="")))
    evaluationSDArrByLevel <- array(0,dim=c(length(evaluationColnamesArrByClass),length(evaluationColnamesArrByClass2)-1,length(classFrequencyArr)-1),dimnames=list(evaluationColnamesArrByClass,evaluationColnamesArrByClass2[-which(evaluationColnamesArrByClass2=="Micro")],paste(classFrequencyArr[-length(classFrequencyArr)]+1,"-",classFrequencyArr[-1],sep="")))  
  } else {
    evaluationArrByLevel <- c()
    evaluationArrByClassByLevel <- c()
    evaluationSDArrByLevel <- c()
  }
   
    ids <- sort(evaluateStructSorted$id) # evaluateStructMat[,'id']

    # if (species==1) {
      # id <- id[which(id<=dividingID)]
    # } else if (species==2) {
      # id <- id[which(id>dividingID)]
    # }
    # if (length(id)==0) next

    evaluateData <- createData(id=ids,idClassMat,idFeatureMat,dummyType=1,removeSizeOne=FALSE)
    if (is.factor(evaluateData$dummyClass[,1])) {            
      dummyMatrix<-data.matrix(evaluateData$dummyClass)-1
    } else {
      dummyMatrix<-evaluateData$dummyClass
    }    
    
    # evaluate by protein; 2013
    
    if (output || sortBy=="F" || sortBy=="AUC") {
    
      longestCutoffs <- c()
      for (id in ids) { # s is TRUE; others are FALSE
        print(paste("protein=",id))
        positiveIDs <- dummyMatrix[as.character(id),]>0
        numTrue <- sum(positiveIDs)
        if (numTrue!=ncol(dummyMatrix) && numTrue!=0) {
          cutoffs <- prediction(evaluateStructSorted$prob[as.character(id),],positiveIDs)@cutoffs[[1]]
          longestCutoffs <- unique(sort(c(longestCutoffs,cutoffs)))
        }
      }
      cutoffs<-longestCutoffs

      NN <- 0      
      NNvec <- c()
      for (cut in cutoffs) {
        totalPredicted <- 0
        IDpredicted <- c()
        for (id in ids) {
          #print(paste("protein=",id))
          positiveIDs <- dummyMatrix[as.character(id),]>0
          predicted <- as.integer(sum(evaluateStructSorted$prob[as.character(id),]>=cut)>0)      
          if (predicted > 0) {
            totalPredicted <- totalPredicted +1
            IDpredicted <- c(IDpredicted,id)
          }
        }
        if (totalPredicted > NN) {
          NN <- totalPredicted
          NNvec <- IDpredicted
        }      
      }
      
      fVec <- c()
      precisionVec <- c()
      recallVec <- c()
      aucVec <- c()
      TPSumVec <- c()
      FPSumVec <- c() 
      PSumVec <- c()
      NSumVec <- c()
      for (cut in cutoffs) {
        precisionVecProtein <- c()
        recallVecProtein <- c()
        proteinvecArray <- c()
        PSum <- 0
        NSum <- 0
        TPSum <- 0
        FPSum <- 0
        positiveIDs <- c()        
        probVec <- c()
        for (id in NNvec) {
          #print(paste("protein=",id))
          
          probVec <- c(probVec,evaluateStructSorted$prob[as.character(id),])
          predicted <- sum(probVec>=cut)
          
          positiveIDs <- c(positiveIDs,dummyMatrix[as.character(id),]>0)
          # 改進          
          proteinvecArray <- rbind(proteinvecArray,dummyMatrix[as.character(id),]>0)          
          
          numTrue <- sum(positiveIDs)
          #print(numTrue)          
          if (numTrue!=length(probVec) && numTrue!=0) {
            pred <- prediction(probVec,positiveIDs)            

            P <- length(which(positiveIDs==TRUE))
            PSum <- PSum + P
            N <- length(which(positiveIDs==FALSE))
            NSum <- NSum + N
          
            perf <- performance(pred, "tpr")
            TP <- perf@y.values[[1]][length(which(perf@x.values[[1]]>=cut))]*P
            TPSum <- TPSum + TP
            perf <- performance(pred, "fpr")
            FP <- perf@y.values[[1]][length(which(perf@x.values[[1]]>=cut))]*N
            FPSum <- FPSum + FP

            if (TP==0) {
              precision <- 0
            } else {
              precision <- TP/(TP+FP)
            }
            if (predicted>0) {
              niter <- length(positiveIDs)/ncol(dummyMatrix)
              for (i in 1:niter) {
                #print(paste("add",precision))
                precisionVecProtein <- c(precisionVecProtein,precision)
              }
            }

            if (TP==0) {
              recall <- 0
            } else {
              recall <- TP/P
            }
            niter <- length(positiveIDs)/ncol(dummyMatrix)
            for (i in 1:niter) {
              #print(paste("add",recall))
              recallVecProtein <- c(recallVecProtein,recall)
            }
            positiveIDs <- c()
            probVec <- c()
          }
        }
        #print(paste("cut=",cut))
        #print(precisionVecProtein)
        #print(recallVecProtein)
        precision <- mean(precisionVecProtein) # over predicted proteins
        recall <- mean(recallVecProtein) # over NN proteins
        f <- 2*precision*recall/(precision+recall) # over NN proteins
        #print(f)
        rownames(proteinvecArray) <- rownames(evaluateStructSorted$prob)                
        pred <- prediction(evaluateStructSorted$prob,proteinvecArray) # over NN proteins
        auc <- performance(pred, "auc", fpr.stop=1)@y.values[[1]] 
        
        fVec <- c(fVec,f)
        precisionVec <- c(precisionVec,precision)
        recallVec <- c(recallVec,recall)
        aucVec <- c(aucVec,auc)
          
        # For micro AUC; over NN proteins
        PSumVec <- c(PSumVec,PSum/NN)
        NSumVec <- c(NSumVec,NSum/NN)
        TPSumVec <- c(TPSumVec,TPSum/NN)
        FPSumVec <- c(FPSumVec,FPSum/NN)        
      }      
      evaluationArr["F",c("Weighted macro","Macro")] <- max(fVec,na.rm=TRUE)
      evaluationArr["F","Micro"]<-max(2*precisionVec*recallVec/(recallVec+precisionVec),na.rm = TRUE)
      evaluationArr["AUC",c("Weighted macro","Macro")] <- max(aucVec,na.rm=TRUE)
      
      tprVec<-TPSumVec/PSumVec
      #tprVec<-c(tprVec,1)
      fprVec<-FPSumVec/NSumVec
      #fprVec<-c(fprVec,1)   
      evaluationArr["AUC","Micro"] <- abs(trapz(fprVec, tprVec))      
    
      if (sortBy=="F") {
        perArr[,] <- evaluationArr["F",]
      } else if (sortBy=="AUC") {
        perArr[,] <- evaluationArr["AUC",]
      }
    }

    if (output || sortBy=="Max class rank") {
      rankbyprotein<-t(apply(1-evaluateStructSorted$prob,1,function(x) rank(x,ties.method="first")))

      # Yu et al.; IEEE/ACM TCBB 2013
      # average of max class ranks to cover all truths for each protein
      evaluationArr["Max class rank",]=mean(sapply(rownames(rankbyprotein),function(x) max(rankbyprotein[x,][dummyMatrix[x,]==1])))-1
      
      if (sortBy=="Max class rank") {
        perArr[,] <- evaluationArr["Max class rank",1]
      }
    }
    
    if (output || sortBy=="Median rank by class") {
      rankbyclass<-apply(1-evaluateStructSorted$prob,2,function(x) rank(x,ties.method="first"))
      
      # Pandey et al.; BMC 2009
      # average of median protein ranks to uncover all truths for each class
      evaluationArrByClass["Median rank by class","Macro",]<-sapply(colnames(rankbyclass),function(x) median(rankbyclass[,x][dummyMatrix[,paste("g",x,sep="")]==1]))
      evaluationArrByClass["Median rank by class","Weighted macro",]=evaluationArrByClass["Median rank by class","Macro",]

      evaluationArr["Median rank by class",] <- mean(evaluationArrByClass["Median rank by class",2,])
      evaluationSDArr["Median rank by class",] <- sd(evaluationArrByClass["Median rank by class",2,])   

      if (sortBy=="Median rank by class") {
        perArr[,] <- evaluationArr["Median rank by class",]
        perArrByClass[,] <- evaluationArrByClass["Median rank by class",,]        
        perSDArr[,] <- evaluationSDArr["Median rank by class",]
      }
    }
    
    if (output || sortBy=="MSE" || sortBy=="Cross entropy") {
      evaluationArr["MSE",] <- mean((evaluateStructSorted$prob-dummyMatrix)^2)
      
      if (sortBy=="MSE") {
        perArr[,] <- evaluationArr["MSE",]
      }
    }
    
    if (output || sortBy=="MSE" || sortBy=="Cross entropy") {
      entropyprob <- evaluateStructSorted$prob*dummyMatrix
      evaluationArr["Cross entropy",] <- -mean(log2(entropyprob[entropyprob!=0]))
      
      if (sortBy=="Cross entropy") {
        perArr[,] <- evaluationArr["Cross entropy",]
      }
    }    

    if (output || (sortBy!=0 && which(evaluationColnamesArr==sortBy)<=which(evaluationColnamesArr=="F by class"))) {

      classByFrequency <- as.character(evaluateData$class[order(allNumAnn[as.character(evaluateData$class)])])

      longestCutoffs <- c()
      for (s in classByFrequency) { # s is TRUE; others are FALSE
        #print(paste("class=",s))
        class <- paste("g",s,sep="")
        classVec <- dummyMatrix[,class]>0
        cutoffs <- prediction(evaluateStructSorted$prob[evaluateStructSorted$id %in% evaluateData$id,s],classVec)@cutoffs[[1]]
        longestCutoffs <- unique(sort(c(longestCutoffs,cutoffs)))
      }
      cutoffs<-longestCutoffs

      TPSumVec <- numeric(length(cutoffs))
      TPSumVecList <- list(numeric(length(cutoffs)),numeric(length(cutoffs)),numeric(length(cutoffs)),numeric(length(cutoffs)))
      FPSumVec <- numeric(length(cutoffs))
      FPSumVecList <- list(numeric(length(cutoffs)),numeric(length(cutoffs)),numeric(length(cutoffs)),numeric(length(cutoffs)))
      PSumList <- list(0,0,0,0)
      NSumList <- list(0,0,0,0)
      QSumList <- list(0,0,0,0)
      classByLevel <- list(c(),c(),c(),c())

      # micro average
      PSum <- 0
      NSum <- 0
      QSum <- 0

      b <- 1
      nC <- length(classByFrequency)
      for (s in classByFrequency) { # s is TRUE; others are FALSE
        #print(paste("class=",s))
        while (allNumAnn[s]>classFrequencyArr[b]) {
          b <- b + 1
        }
        #print(paste("in slot",b))
        if (b>1) {
          classByLevel[[b-1]] <- c(classByLevel[[b-1]], s)
          class <- paste("g",s,sep="")
          classVec <- dummyMatrix[,class]>0
          P <- length(which(classVec==TRUE)) # 21
          N <- length(which(classVec==FALSE)) # 76

          pred <- prediction(evaluateStructSorted$prob[evaluateStructSorted$id %in% evaluateData$id,s],classVec) 

          if (output || sortBy=="AUC by class") {
            evaluationArrByClass["AUC by class",c("Weighted macro","Macro"),s] <- performance(pred, "auc", fpr.stop=1)@y.values[[1]]
            if (byLevel) {
              evaluationArrByClassByLevel["AUC by class","Weighted macro",s,b-1] <- evaluationArrByClass["AUC by class","Macro",s]*allQ[as.character(s)]
              evaluationArrByClassByLevel["AUC by class","Macro",s,b-1] <- evaluationArrByClass["AUC by class","Macro",s]
            }
            if (sortBy=="AUC by class") {
              perArrByClass[c("Weighted macro","Macro"),s] <- evaluationArrByClass["AUC by class","Macro",s]
            }
          }

          PSum <- PSum + P
          NSum <- NSum + N
          QSum <- QSum + allQ[as.character(s)]
          if (byLevel) {
            PSumList[[b-1]] <- PSumList[[b-1]] + P
            NSumList[[b-1]] <- NSumList[[b-1]] + N
            QSumList[[b-1]] <- QSumList[[b-1]] + allQ[as.character(s)]
          }

          perf <- performance(pred, "tpr")
          TPvec <- c()
          i <- 1
          for (cut in cutoffs) {
            TP <- perf@y.values[[1]][length(which(perf@x.values[[1]]>=cut))]*P
            TPvec <- c(TPvec,TP)
            TPSumVec[i]<- TPSumVec[i]+TP
            if (byLevel) {
              TPSumVecList[[b-1]][i] <- TPSumVecList[[b-1]][i]+TP
            }
            i <- i + 1
          }

          perf <- performance(pred, "fpr")
          FPvec <- c()
          i <- 1
          for (cut in cutoffs) {
            FP <- perf@y.values[[1]][length(which(perf@x.values[[1]]>=cut))]*N
            FPvec <- c(FPvec,FP)
            FPSumVec[i]<- FPSumVec[i]+FP
            if (byLevel) {
              FPSumVecList[[b-1]][i]<- FPSumVecList[[b-1]][i]+FP
            }
            i <- i + 1
          }
          #TN=64

          if (output || sortBy=="Precision by class" || sortBy=="F by class") {
            precisionVec<-TPvec/(TPvec+FPvec)
            precisionVec<-c(precisionVec,0)
            recallVec<-TPvec/P
            recallVec<-c(recallVec,1)
            recallindex<-which(recallVec>=0.1)[1]

            evaluationArrByClass["Precision by class",c("Weighted macro","Macro"),s] <- precisionVec[recallindex]
            if (byLevel) {
              evaluationArrByClassByLevel["Precision by class","Weighted macro",s,b-1] <- evaluationArrByClass["Precision by class","Macro",s]*allQ[as.character(s)]
              evaluationArrByClassByLevel["Precision by class","Macro",s,b-1] <- evaluationArrByClass["Precision by class","Macro",s]
            }

            evaluationArrByClass["F by class",c("Weighted macro","Macro"),s] <- max(2*precisionVec*recallVec/(recallVec+precisionVec),na.rm = TRUE)
            if (byLevel) {
              evaluationArrByClassByLevel["F by class","Weighted macro",s,b-1] <- evaluationArrByClass["F by class","Macro",s]*allQ[as.character(s)]
              evaluationArrByClassByLevel["F by class","Macro",s,b-1] <- evaluationArrByClass["F by class","Macro",s]
            }
            if (sortBy=="Precision by class") {
              perArrByClass[c("Weighted macro","Macro"),s] <- evaluationArrByClass["Precision by class","Macro",s]
            } else if (sortBy=="F by class") {
              perArrByClass[c("Weighted macro","Macro"),s] <- evaluationArrByClass["F by class","Macro",s]
            }
          }
        }
      }

      if (output || sortBy=="Precision by class" || sortBy=="F by class") {
        precisionVec<-TPSumVec/(TPSumVec+FPSumVec)
        #precisionVec<-c(precisionVec,0)
        recallVec<-TPSumVec/PSum
        #recallVec<-c(recallVec,1)
      }
       
      if (output || sortBy=="AUC by class") {       
        tprVec<-TPSumVec/PSum
        #tprVec<-c(tprVec,1)
        fprVec<-FPSumVec/NSum
        #fprVec<-c(fprVec,1)
      
        evaluationArr["AUC by class","Micro"]<-abs(trapz(fprVec, tprVec))

        evaluationArr["AUC by class","Weighted macro"]<-sum(evaluationArrByClass["AUC by class","Macro",]*allQ[classByFrequency])/QSum
        evaluationSDArr["AUC by class","Weighted macro"]<-sd(evaluationArrByClass["AUC by class","Macro",]*allQ[classByFrequency]/QSum)
        evaluationArr["AUC by class","Macro"]<-mean(evaluationArrByClass["AUC by class","Macro",])
        evaluationSDArr["AUC by class","Macro"]<-sd(evaluationArrByClass["AUC by class","Macro",])

        if (sortBy=="AUC by class") {
          perArr["Micro","total"]<-evaluationArr["AUC by class","Micro"]
        }
        
        if (output) {
          jpeg(paste(filePrefixEvaluate,".ROC.jpg",sep=""))
          plot(fprVec,tprVec,type="l",xlab="False Positive Rate", ylab="True Positive Rate",xaxt="n", col="blue",main="ROC Curves")
          axis(1, c(0,0.2,0.4,0.6,0.8,1))
          dev.off()               
        }          
      } 
      if (output || sortBy=="Precision by class") {
        recallindex<-which(recallVec>=0.1)[1]
      
        evaluationArr["Precision by class","Micro"]<-precisionVec[recallindex]

        evaluationArr["Precision by class","Weighted macro"]<-sum(evaluationArrByClass["Precision by class","Macro",]*allQ[classByFrequency])/QSum
        evaluationSDArr["Precision by class","Weighted macro"]<-sd(evaluationArrByClass["Precision by class","Macro",]*allQ[classByFrequency]/QSum)
        evaluationArr["Precision by class","Macro"]<-mean(evaluationArrByClass["Precision by class","Macro",])
        evaluationSDArr["Precision by class","Macro"]<-sd(evaluationArrByClass["Precision by class","Macro",])
              
        if (sortBy=="Precision by class") {
          perArr["Micro","total"]<-evaluationArr["Precision by class","Micro"]
        }
        
        if (output) {
          jpeg(paste(filePrefixEvaluate,".PR.jpg",sep=""))
          plot(recallVec,precisionVec,type="l",xlab="Recall", ylab="Precision",xaxt="n", col="blue",main="PR Curves")
          axis(1, c(0,0.2,0.4,0.6,0.8,1))
          dev.off()         
        }          
      } 
      if (output || sortBy=="F by class") {        
        evaluationArr["F by class","Micro"]<-max(2*precisionVec*recallVec/(recallVec+precisionVec),na.rm = TRUE)

        evaluationArr["F by class","Weighted macro"]<-sum(evaluationArrByClass["F by class","Macro",]*allQ[classByFrequency])/QSum
        evaluationSDArr["F by class","Weighted macro"]<-sd(evaluationArrByClass["F by class","Macro",]*allQ[classByFrequency]/QSum)
        evaluationArr["F by class","Macro"]<-mean(evaluationArrByClass["F by class","Macro",])
        evaluationSDArr["F by class","Macro"]<-sd(evaluationArrByClass["F by class","Macro",])
      
        if (sortBy=="F by class") {
          perArr["Micro","total"]<-evaluationArr["F by class","Micro"]
        }
      }
      
      if (sortBy!=0 && which(evaluationColnamesArr==sortBy)<=which(evaluationColnamesArr=="F by class")) {
        perArr["Weighted macro","total"]<-sum(evaluationArrByClass[sortBy,"Macro",]*allQ[classByFrequency])/QSum
        perSDArr["Weighted macro","total"]<-sd(evaluationArrByClass[sortBy,"Macro",]*allQ[classByFrequency]/QSum)
        perArr["Macro","total"]<-mean(evaluationArrByClass[sortBy,"Macro",])
        perSDArr["Macro","total"]<-sd(evaluationArrByClass[sortBy,"Macro",])
      }
      
      if (output && byLevel) {
        for (b in 1:(length(classFrequencyArr)-1)) {
          if (length(classByLevel[[b]])==0) next
          filePrefixevaluateb <- paste(filePrefixEvaluate,".r",classFrequencyArr[b]+1,"-",classFrequencyArr[b+1],sep="")

          print(paste("range ",classFrequencyArr[b]+1,"-",classFrequencyArr[b+1],sep=""))

          precisionVec<-TPSumVecList[[b]]/(TPSumVecList[[b]]+FPSumVecList[[b]])
          precisionVec<-c(precisionVec,0)
          recallVec<-TPSumVecList[[b]]/PSumList[[b]]
          recallVec<-c(recallVec,1)

          tprVec<-TPSumVecList[[b]]/PSumList[[b]]
          tprVec<-c(tprVec,1)
          fprVec<-FPSumVecList[[b]]/NSumList[[b]]
          fprVec<-c(fprVec,1)

          evaluationArrByLevel["AUC by class","Micro",b]<-trapz(fprVec, tprVec)
          evaluationArrByLevel["AUC by class","Weighted macro",b]<-sum(evaluationArrByClassByLevel["AUC by class","Macro",classByLevel[[b]],b])/QSumList[[b]]
          evaluationSDArrByLevel["AUC by class","Weighted macro",b]<-sd(evaluationArrByClassByLevel["AUC by class","Macro",classByLevel[[b]],b]/QSumList[[b]])
          evaluationArrByLevel["AUC by class","Macro",b]<-mean(evaluationArrByClassByLevel["AUC by class","Macro",classByLevel[[b]],b])
          evaluationSDArrByLevel["AUC by class","Macro",b]<-sd(evaluationArrByClassByLevel["AUC by class","Macro",classByLevel[[b]],b])

          jpeg(paste(filePrefixevaluateb,".ROC.jpg",sep=""))
          plot(fprVec,tprVec,type="l",xlab="False Positive Rate", ylab="True Positive Rate",xaxt="n", col="blue",main="ROC Curves")
          axis(1, c(0,0.2,0.4,0.6,0.8,1))
          dev.off()

          recallindex<-which(recallVec>=0.1)[1]
          evaluationArrByLevel["Precision by class","Micro",b]<-precisionVec[recallindex]
          evaluationArrByLevel["Precision by class","Weighted macro",b]<-sum(evaluationArrByClassByLevel["Precision by class","Weighted macro",classByLevel[[b]],b])/QSumList[[b]]
          evaluationSDArrByLevel["Precision by class","Weighted macro",b]<-sd(evaluationArrByClassByLevel["Precision by class","Weighted macro",classByLevel[[b]],b]/QSumList[[b]])
          evaluationArrByLevel["Precision by class","Macro",b]<-mean(evaluationArrByClassByLevel["Precision by class","Macro",classByLevel[[b]],b])
          evaluationSDArrByLevel["Precision by class","Macro",b]<-sd(evaluationArrByClassByLevel["Precision by class","Macro",classByLevel[[b]],b])

          jpeg(paste(filePrefixevaluateb,".PR.jpg",sep=""))
          plot(recallVec,precisionVec,type="l",xlab="Recall", ylab="Precision",xaxt="n", col="blue",main="PR Curves")
          axis(1, c(0,0.2,0.4,0.6,0.8,1))
          dev.off()

          evaluationArrByLevel["F by class","Micro",b]<-max(2*precisionVec*recallVec/(recallVec+precisionVec),na.rm = TRUE)
          evaluationArrByLevel["F by class","Weighted macro",b]<-sum(evaluationArrByClassByLevel["F by class","Weighted macro",classByLevel[[b]],b])/QSumList[[b]]
          evaluationSDArrByLevel["F by class","Weighted macro",b]<-sd(evaluationArrByClassByLevel["F by class","Weighted macro",classByLevel[[b]],b]/QSumList[[b]])
          evaluationArrByLevel["F by class","Macro",b]<-mean(evaluationArrByClassByLevel["F by class","Macro",classByLevel[[b]],b])
          evaluationSDArrByLevel["F by class","Macro",b]<-sd(evaluationArrByClassByLevel["F by class","Macro",classByLevel[[b]],b])

        }
      }
    }
    
    return(list(
    perArr=perArr,perArrByClass=perArrByClass,perSDArr=perSDArr,
    evaluationArr=evaluationArr,evaluationArrByClass=evaluationArrByClass,evaluationSDArr=evaluationSDArr,
    evaluationArrByLevel=evaluationArrByLevel,evaluationArrByClassByLevel=evaluationArrByClassByLevel,evaluationSDArrByLevel=evaluationSDArrByLevel
    ))
}

if (ANOVA) {
  for (dummyType in dummyTypeVec) {
    trainValidateData<-createData(id=allIDs,idClassMat,idFeatureMat,dummyType=dummyType)
    saveRDS(trainValidateData, file = paste(filePrefixClass,".trainValidateData.",dummyType,".Rda",sep=""))
    #library(biotools) # requires GUI
    #sink(paste(filePrefixClass,".BoxM.txt",sep=""))
    #nosig <- 0
    #for (c in 1:ncol(trainData$dummyClass)) {
    #  tryCatch(tmp<-boxM(trainData$features, trainData$dummyClass[,c]),warning=function(w) {print(w); print(colnames(trainData$dummyClass)[c])})
    #  if (tmp$id.value<0.001) {
    #    nosig <- nosig + 1
    #    print(colnames(trainData$dummyClass)[c])
    #  }
    #}
    #print(paste("covariance not equal: ",nosig))
    #sink()
    myData <- data.frame(trainValidateData$id,trainValidateData$features,trainValidateData$dummyClass)
    lmMod <- lm(as.formula(trainValidateData$formulastr),data=myData, x=TRUE, y=TRUE, model=TRUE) # trainValidateData$data
    lmModAnova <- try(Anova(lmMod,test="Wilks",type="2"),silent=TRUE)
    sink(paste(filePrefixClass,".anova.",dummyType,".Rda",sep=""))
    print(lmModAnova)
    sink()
  }
}

createSaveData(idClassMat,idFeatureMat,filePrefixClass,idClassMatTest,idFeatureMatTest)
otherMethods(filePrefixClass,idClassMat,idFeatureMat,idClassMatTest,idFeatureMatTest,allQ,allNumAnn)

if (doFitLinearModels) {
  fitLinearModels(idClassMat,idFeatureMat,filePrefixClass)
}

parameters <- calEigen(classWeightedVec=classWeightedVec,proportionalHVec=proportionalHVec,
  linearModelTypeVec=linearModelTypeVec,regularizeWTypeVec=regularizeWTypeVec,regularizeWLambdaVec=regularizeWLambdaVecOrigin,
  idClassMat,allQ,filePrefixClass)
if (is.na(parameters$linearModelTypeVec)) {
  print("no result")
  next
}

results<-transformPredictAndEvaluate(idClassMat,idFeatureMat,allQ,
                    classWeightedVec=parameters$classWeightedVec,
                    proportionalHVec=parameters$proportionalHVec,
                    linearModelTypeVec=parameters$linearModelTypeVec,
                    regularizeWTypeVec=parameters$regularizeWTypeVec,
                    regularizeWLambdaVec=parameters$regularizeWLambdaVec,
                    validateOnly,
                    idClassMatTest=idClassMatTest,
                    idFeatureMatTest=idFeatureMatTest,
                    filePrefixClass)    

q()

#clusters <- list()
#nums <- c(114, 400, 715, 4644)
#for (i in 1:4) {
#  clusters[[i]] <- read.csv(paste("clus2",suffix2,".",nums[i],".csv",sep=""))
#}

complexfile<-read.csv(paste(commonDir,dataNames,,".complexs",sep=""),header=FALSE)
className<-"complex"
idClassMat<-complexfile

ecfile<-read.csv(paste(commonDir,dataNames,,".EC",sep=""),header=FALSE)
className<-"EC"
idClassMat<-ecfile

orthologfile<-read.csv(paste(commonDir,dataNames,".orthologs",sep=""),header=FALSE)
className<-"ortholog"
idClassMat<-orthologfile

slimtermfile<-read.csv(paste(commonDir,dataNames,".slim_term.no",sep=""),header=FALSE)
className<-"slimterm"
idClassMat<-slimtermfile

q()

lm.mod<-lm(cbind(Sepal.Length,Sepal.Width,Petal.Length,Petal.Width)~Species, data=iris,y=TRUE,x=TRUE)
lmPred<-predict(lm.mod)
mean((lm.mod$fitted.values-lm.mod$y)^2)
#0.148829
lm.anova <- Anova(lm.mod)
lm.mod.ridge<-lm.ridge(cbind(Sepal.Length,Sepal.Width,Petal.Length,Petal.Width)~Species, data=iris,lambda=c(0,10^seq(10,-2,length=100)))
mean((lm.mod.ridge$fitted.values-lm.mod.ridge$y)^2)
#[1] 0.1488314

distNo<-36
suffix<-""
suffix2<-""
idClassMat<-cbind(1:dim(lm.mod$y)[1],iris$Species)
idFeatureMatOrigin<-cbind(1:dim(lm.mod$y)[1],lm.mod$y)

dataNames<<-"iris3"
className<-"iris"
dividingID<<-150
commonDir<<-""
useAncestor<<-FALSE
filePrefix<<-paste(dataNames,".dist.",distNo,suffix,sep="")

dummyType<<-1 # dummyType<<-5
factorize<<-FALSE

id<-1:150

dummyType=1
[1] "combineParameters-sortBy-validateOnly-TRUE"
[1] "pertotaltotal= 0.998933333333333"
[1] "pertotaltotal2= 0.994266666666667"
  with proportional calH.SSPS => No change

dummyType=1; factorize
[1] "combineParameters-sortBy-validateOnly-TRUE"
[1] "pertotaltotal= 0.998933333333333"
[1] "pertotaltotal2= 0.994266666666667"
  with proportional calH.SSPS => No change

dummyType=5
[1] "combineParameters-sortBy-validateOnly-TRUE"
[1] "pertotaltotal= 0.998933333333333"
[1] "pertotaltotal2= 0.994266666666667"
  with proportional calH.SSPS => No change

dummyType=5

q()

######### Prostate ##########

testProstate <- function() {

  library(lasso2)
  data(Prostate)
  Prostate$gleasonb<-rep(0,nrow(Prostate))
  Prostate$gleasonb[which(Prostate$gleason>=7)]<-1

  lm.mod.prostate<-lm(cbind(lcavol,age,lpsa,lweight,lbph,lcp,pgg45) ~ I(svi==1)+I(gleasonb==1), data=Prostate,y=TRUE, x=TRUE)

    #lmPred.prostate<-predict(lm.mod.prostate)
    #mean((lm.mod.prostate$fitted.values-lm.mod.prostate$y)^2)
      #[1] 67.11472
    #lm.anova.prostate <- Anova(lm.mod.prostate)
    #lm.mod.ridge.prostate<-lm.ridge(cbind(lcavol,age,lpsa,lweight,lbph,lcp,pgg45) ~ I(svi==1)+I(gleasonb==1),data=Prostate, lambda=c(0,10^seq(10,-2,length=100)))
    #mean((lm.mod.ridge.prostate$fitted.values-lm.mod.ridge.prostate$y)^2)
      #[1] 67.19024

  id<-1:dim(lm.mod.prostate$y)[1]
  index<-(1:length(id))[Prostate$svi==0&Prostate$gleasonb==0]
  idClassMat<-cbind(index,array(1,dim=c(length(index))))
  index2<-(1:length(id))[Prostate$svi==1]
  idClassMat<-rbind(idClassMat,cbind(index2,array(2,dim=length(index2))))
  index3<-(1:length(id))[Prostate$gleasonb==1]
  idClassMat<-rbind(idClassMat,cbind(index3,array(3,dim=length(index3))))
  idClassMat<-idClassMat[order(idClassMat[,1],idClassMat[,2]),]
  idClassMatOrigin<-idClassMat

  randomtest<-sort(sample(id,15))
  idClassMatTest <- idClassMat[idClassMat[,1] %in% randomtest,]
  idClassMatTest <- t(sapply(1:dim(idClassMatTest)[1],function(x) c(which(idClassMatTest[x,1]==randomtest)+length(id),idClassMatTest[x,2])))
  rownames(idClassMatTest)<-as.character(idClassMatTest[,1])

  idFeatureMatOrigin<-cbind(id,lm.mod.prostate$y)

  idFeatureMatTest<-idFeatureMatOrigin[idFeatureMatOrigin[,1] %in% randomtest,]
  idFeatureMatTest<-cbind((1:length(randomtest))+length(id),idFeatureMatTest[,-1])
  rownames(idFeatureMatTest)<-as.character(idFeatureMatTest[,1])

  idFeatureMatOrigin<-rbind(idFeatureMatOrigin,idFeatureMatTest)

    # lmData <- createData(id,idClassMat,idFeatureMatOrigin,dummyType=dummyType)
    ## index<-which(colnames(lmData$data)=="g1")
    ## dfh<-2  # checked
    ## dummyClass<-as.matrix(lmData$data[,index:(index+dfh)])
    ## dummyClass <- t(apply(dummyClass,1,function(x) as.numeric(x)))
    ## trainFeatures <- lm.mod.prostate$y

    # data <- data.frame(lmData$id,lmData$features,lmData$dummyClass)
    # lm.mod.prostate.I<-lm(as.formula(lmData$formulastr), data=data,y=TRUE, x=TRUE) # lmData$data
    # lmPred.prostate.I<-predict(lm.mod.prostate.I)
    # mean((lm.mod.prostate.I$fitted.values-lm.mod.prostate.I$y)^2)
      ##[1] 67.11472
    ## lm.anova.prostate.I <- Anova(lm.mod.prostate.I)
    # lm.mod.ridge.prostate.I<-lm.ridge(as.formula(lmData$formulastr), data=data, lambda=c(0,10^seq(10,-2,length=100))) # lmData$data
    # mean((lm.mod.ridge.prostate.I$fitted.values-lm.mod.ridge.prostate.I$y)^2)
      ##[1] 67.23859 #unfactorize; dummyType=1
      ##[1] 67.39024 #unfactorize; dummyType=4, 5

  distNo<-36
  suffix<-""
  suffix2<-""
  className<-"prostate"

  dataNames<<-"prostate"
  dividingID<<-120
  commonDir<<-""
  filePrefix<<-paste(dataNames,".dist.",distNo,suffix,sep="")
  
}

##################################################################

#dummyType=1
combineParameters-sortBy-validateOnly-FALSE-2-TRUE
[1] "pertotaltotal= 0.928545530036135"
[1] "pertotaltotal2= 0.927968876537982"
    with proportional calH.SSPS
"combineParameters-sortBy-validateOnly-FALSE-2-TRUE"
[1] "pertotaltotal= 0.929277985833912"
[1] "pertotaltotal2= 0.924387581857871"

#dummyType=1; factorize
  with proportional calH.SSPS
[1] "pertotaltotal= 0.929277985833912"
[1] "pertotaltotal2= 0.924387581857871"

#dummyType=5
"combineParameters-sortBy-validateOnly-FALSE-2-TRUE"
[1] "pertotaltotal= 0.949402483792861"
[1] "pertotaltotal2= 0.939678200421776"
    with proportional calH.SSPS
"combineParameters-sortBy-validateOnly-FALSE-2-TRUE"
[1] "pertotaltotal= 0.949402483792861"
[1] "pertotaltotal2= 0.940193704600484"

#dummyType=5; factorize
    with proportional calH.SSPS
"combineParameters-sortBy-validateOnly-FALSE-2-TRUE"
[1] "pertotaltotal= 0.88082871202062"
[1] "pertotaltotal2= 0.958990080449895"

prostate2.dist.36.FALSE.prostate.TRUE-3
    with no proportional calH.SSPS; with weights
[1] "pertotaltotal= 0.947359993751465"
[1] "pertotaltotal2= 0.89190424119347"