#################################################################
#                                                               #
#   CERFIT (causal inference random forest interaction trees)   #
#                                                               #
#################################################################
# This code is created for causal inference random forest interaction trees.
# This code can implement individualized treatment effect prediction for binary treatments,
# as well as non-bianry treatments. 
# For observational data, the propensity scores will be first estimated depending on
# number of treatments/levels. For binary, we recommend GBM, CBPS and 
# Random Forst methods. For less than 4 treatments, CBPS ran much fast than GBM,
# If more than 4 discrete treatments, using GBM. For continuous, we recommend CBPS
# or HI methods. 
# Code cannot handle missing data.
# Code for Binary response or survival time to event data is still under development
# Tree converted to constparty object using partykit R package.
# Code is developed by Laura Li and Peter Calhoun  


### Common input variables and description ###
# 
# ntrees - Number of trees grown
# formula - Formula to build CERFIT.  Categorical predictors must be listed as
#           a factor. 
#           e.g., Y ~ x1 + x2 | treatment
# data - Data to grow a tree
# search - Method to search through candidate splits.
#          Options are "exhaustive","sss". *"sss" is still under development.
# method - For observational stuy data, method="observation"; 
#         for randomized study data, method="RCT".
# split - Impurity measure splitting statistic
#         For continuous, split can be "ttest.obs" or "ttest.rand"
#         For binary, split can be "ttest.glm" 
# mtry - Number of variables randomly selected 
# nsplit - Number of cutpoints selected (for "exhaustive")
# nsplit.random - Logical: indicates if process to select cutpts are random (for "exhaustive")
# minsplit - Number of observations required to continue growing tree
# minbucket - Number of observations required in each child node
# maxdepth - Maximum depth of tree
# a - Sigmoid approximation variable (for "sss" which is still under development) 
# sampleMethod - Method to sample learning sample.
#                Options are "bootstrap", "subsample", "subsampleByID", "learning"
# scale.y - Logical, standardize y when creating splits
#           For "sss" to increase stability


### Load Functions ###
#rm(list=ls(all=TRUE))
library(partykit)
library(parallel)
library(pROC)
library(CBPS)
library(randomForest)
library(twang)
library(glmnet)
library(sandwich)


###Finds candidate cutpoints ###
findCutpts <- function(x, minbucket) {
  nX <- length(x)
  #Due to ties, it's possible minbucket cannot be satisfied
  if (sort(x)[minbucket]==sort(x, decreasing=TRUE)[minbucket]) {
    cutpts=NULL} 
  else {
    #Cutpoints considered must satisfy minbucket
    cutpts <- unique(sort(x)[minbucket:(nX-minbucket+1)])

    if(length(cutpts)==1){stop(paste0("Only 1 cutpt??? ", cutpts, x))}
    cutpts <- (cutpts[1:(length(cutpts)-1)]+cutpts[2:length(cutpts)])/2
  }
  return(cutpts)
}

### CERFIT: Extract scores from interaction model###

ttest.obs <- function(y, trt, x, cutpt, propensity,minbucket)#,detail=FALSE)# 
{  
  score <- NA; n <- length(y)
  z <- sign(x<=cutpt)
  if(length(unique(trt))==2) {
    mingroup<- min(lengths(split(y,paste(trt,z))))} else {
      yz<-unique(cbind(y,z))
      mingroup<- min(lengths(split(yz[,1],yz[,2])))
      mingroup<-mingroup/2}## if continuous trt, unique y in each side of childnodes should greater than 2*minbucket
  if(min(sum(z),sum(1-z))<0.1*length(y) | mingroup < minbucket/2 ) {
    return(NA)} else{
      mod1 <- try(lm(y ~ trt + z + propensity+z*trt), silent=TRUE)#, data=dat))#,weights = w) # fitting linear model with interaction with propensity score
      if (inherits(mod1, "try-error" ) | is.nan(mod1$coefficients[5])){
        return(NA)} else {
          t<-try(summary(mod1)$coefficients[5,3],silent=TRUE)
          if (inherits(t, "try-error" )) {
            return(NA)} else {
              SE<-sqrt(diag(vcovHC(mod1,type="HC0")))
              score=try((coef(mod1)[5]/SE[5])^2,silent=TRUE)
              if (inherits(score, "try-error" )) {
                return(t^2)} else {
                  return(score)}
            }
        }
    } 
}


ttest.rand <- function(y, trt, x, cutpt, minbucket)#,detail=FALSE)
{  
  score <- NA; n <- length(y)#;dat=as.data.frame(c(y,trt,z,propensity))
  z <- sign(x<=cutpt)
  if(length(unique(trt))==2) {
    mingroup<- min(lengths(split(y,paste(trt,z))))} else {
      yz<-unique(cbind(y,z))
      mingroup<- min(lengths(split(yz[,1],yz[,2])))
      mingroup<-mingroup/2}
  
  if(min(sum(z),sum(1-z))<0.1*length(y) | mingroup < minbucket/2 ) {
    return(NA)}   else {
      mod1 <- try(lm(y ~ trt + z + z*trt), silent=TRUE)#, data=dat))#,weights = w) # fitting linear model with interaction with propensity
      if (inherits(mod1, "try-error")| is.na(mod1$coefficients[4])){
        return(NA)}  else {
          score=(summary(mod1)$coefficients[4,3])^2
          return(score)}
    } }

pickCutpt <- function(allVars, y, minbucket){ 
  v <- sample.int(NCOL(allVars), 1)
  #If variable selected is a factor, pick a random subset of categories.
  #Do not sort by mean.  Otherwise, will be very likely to select categorical variable
  xTemp <- ordinalize(x=allVars[,v], y,  sortCat=FALSE)
  x <- xTemp$x
  #If all x values the same, do not check optimal split
  if (abs(max(x) - min(x)) > 1e-8) {
    cutpts <- findCutpts(x, minbucket)
    if (!is.null(cutpts) ){
      if (is.factor(allVars[,v])) {return(list(varid=v, cutpt=randCutpt, x=x, cutToLvl=xTemp$cutToLvl))
      } else{ 
        randCutpt <-ifelse(length(cutpts)==1, cutpts, sample(cutpts,1))
        return(list(varid=v, cutpt=randCutpt, x=x, cutToLvl=xTemp$cutToLvl))}
    }
  }
}

splitrule <- function(y, x, trt, cutpts, method, propensity, minbucket){#split,
  if (method=="RCT") {
    #MSE: Q_m(T)=1/n_L*sum((y-ybar)^2)
    stat <- vapply(cutpts, ttest.rand, y=y, x=x,trt=trt, minbucket=minbucket,trtlevels=trtlevels,numeric(1))
  } else {
    stat <- vapply(cutpts, ttest.obs, y=y, x=x,trt=trt, propensity=propensity, minbucket= minbucket, numeric(1))
    #print(stat)
  }
  if (all(is.na(stat))) {return(list(cutoff=NA, stat=NA))
  } else {return(list(cutoff = cutpts[which.max(stat)], stat = max(stat,na.rm=TRUE)))}
}


### Convert factors to numerical value. ###
ordinalize <- function(x, y, sortCat=TRUE){
  if (is.factor(x)) {
    x <- factor(x) #Remove factors not listed
    #One can randomly assign a category a distinct numerical value
    if (!sortCat) {
      cutToLvl <- t(sample.int(length(levels(x))))
      colnames(cutToLvl)=levels(x)
    } else {
      #For binary, sort data by proportion in class 1.  For continuous, sort by means
      if (is.factor(y)) {
        cutToLvl <- prop.table(table(y,x),2)[1,,drop=FALSE]
      } else {cutToLvl <- t(vapply(levels(x), function(z){mean(y[x==z])}, numeric(1)))}
    }
    #Convert lvls to numerical value. Slow method. Make this faster later.
    xTemp <- rep(NA,length(x))
    for (lvls in levels(x)) {
      xTemp[x==lvls] <- cutToLvl[colnames(cutToLvl)==lvls]
    }
  } else {
    xTemp <- x
    cutToLvl <- NULL
  }
  return(list(x=xTemp, cutToLvl=cutToLvl))
}


partition<- function(vars, y, trt, propensity, subset, search, method, split, nsplit, nsplit.random,
                     minsplit, minbucket, a, scale.y, useSearch, useOptim,trtlevels){#, allVars
  if (sum(subset) < minsplit) {return(NULL)}
  
  vars <- vars[subset,,drop=FALSE]
  y <- y[subset]
  trt<- trt[subset]
  #print(unique(trt))
  if(length(trtlevels)>2 & length(trtlevels<10)) {
    propensity<-propensity[subset,]
  } else {
    propensity<-propensity[subset]
  }
  trt.length<-length(trtlevels)
  if(trt.length>2 & trt.length < 10) { ## if less than 10 treatments/levels
    ran<- sample(unique(trt),2)
    vars<-subset(vars,trt==ran[1] | trt==ran[2])
    y<-subset(y,trt==ran[1] | trt==ran[2])
    propensity<-subset(propensity,trt==ran[1] | trt==ran[2])
    trt<-subset(trt,trt==ran[1] | trt==ran[2])
    trt<-ifelse(trt==ran[1],1,0)
    propensity<-propensity[,ran[1]] # need to make sure trt is levels as propensity nameorders
  }
  
  
  if (NROW(vars) < 2*minbucket) {return(NULL)}
  if (length(unique(y))==1) {return(NULL)}  
  
  stats<- cutoff<- breakLeft<-NA
  findStats<-sapply(vars,function(x){
    
    if (search=="exhaustive" && !is.null(nsplit) && nsplit.random) {xTemp <- ordinalize(x, y, sortCat=FALSE) 
    } else {xTemp <- ordinalize(x, y, sortCat=TRUE)}
    x <- xTemp$x
    #If all x values the same, do not check optimal split
    if (abs(max(x) - min(x)) > 1e-8) {
      #The SSS partition deals with problems when there is a very small number of observations
      #Use exhaustive search in this case (or set minsplit >= 5)
      if (search=="sss") { #leave sss here for now
        print("sss not ready")
      } else if (search=="exhaustive") { #current codes only work exhaustive search
        cutpts <- findCutpts(x, minbucket)
        if (is.null(nsplit)) { nsplit<- length(cutpts)} 
        #Take nsplit cutpoints (if applicable)
        if (!is.null(nsplit) && !is.null(cutpts) && length(cutpts) > 1) {
          #If nsplit.random is TRUE, take nsplit cutpts randomly.  Otherwise, take nsplit cutpts equally spread out across cutpts
          if (!nsplit.random & length(cutpts) > nsplit) { #if not random select nsplit cut
            cutpts <- unique(cutpts[seq(1, length(cutpts), length.out=nsplit)])
          } else {
            cutpts <- sort(sample(cutpts, min(c(nsplit, length(cutpts))), replace=FALSE))
          }
        }
        #
        #It is possible (unlikely) no cutpoint can satisfy minbucket
        if (!is.null(cutpts)) {
          mod <- splitrule(y=y, x=x, trt=trt,cutpts=cutpts, method=method, propensity = propensity,minbucket=minbucket)
          if (!is.na(mod$stat)) {
            stats <- mod$stat
            if (is.factor(x)) {
              cutoff<- "factor"
              breakLeft <- rep(NA, length(levels(x)))
              breakLeft[levels(x) %in% colnames(xTemp$cutToLvl)[xTemp$cutToLvl <= mod$cutoff]]=1L
              breakLeft[levels(x) %in% colnames(xTemp$cutToLvl)[xTemp$cutToLvl > mod$cutoff]]=2L
              if (all(is.na(breakLeft)) & length(unique(breakLeft))<=1) {stop("Did not find correct cutpoints")}
            }
            else {cutoff <- mod$cutoff; breakLeft<-NA}
          }
        }
        
      }
      else {stop("Unexpected search")}
    }
    return(c(stats,cutoff,breakLeft))})
  #If randomly picking a subset of categories, do not sort by mean.  Would be more likely to select variables when sorted
  
  #If each candidate variable cannot be split (e.g. cannot satisfy minbucket), return null
  if (all(is.na(findStats[1,]))) {return(NULL)}
  if (findStats[2,which.max(findStats[1,])]=="factor") {
    #Index is used for categorical variable splits
    return(partysplit(varid=as.integer(colnames(findStats)[which.max(findStats[1,])]),
                      index=findStats[3,which.max(findStats[1,])],
                      info=list(stats=findStats[1,])))
  } else {
    
    #Breaks is used for continuous variable splits
    return(partysplit(varid=as.integer(colnames(findStats)[which.max(findStats[1,])]),
                      breaks=findStats[2, which.max(findStats[1,])],
                      info=list(stats=findStats[1,])))
  }
}



### Grow tree by using partition() function several times in a recursive loop ###
growTemp <- function(id=1L, depth=1L, data, response, treatment, Propensity, subset, search, method, split,
                     mtry, nsplit, nsplit.random, minsplit, minbucket, maxdepth,
                     a, scale.y, trtlevels){ 
  if (depth > maxdepth) {return(partynode(id=id))}
  y <- data[[response]]
  trt<-data[[treatment]]
  propensity<-data[[Propensity]]
  varSelected <- sort(sample.int(ncol(data)-4, mtry))
  vars <- data[varSelected]
  colnames(vars) <- varSelected #Have columns represent varid
  
  sp <- partition(vars=vars, y=y,  subset=subset,trt=trt,propensity=propensity,
                  search=search, method=method, split=split, nsplit=nsplit, nsplit.random=nsplit.random,
                  minsplit=minsplit, minbucket=minbucket, a=a, scale.y=scale.y,
                  useSearch=useSearch, useOptim=useOptim,trtlevels=trtlevels)
  
  if (is.null(sp)) {return(partynode(id=id))}
  
  # Split the data
  kidids <- kidids_split(sp, data=data)
  depth <- depth + 1
  
  kids <- vector(mode="list", length=max(kidids, na.rm=TRUE))
  for (kidid in seq_along(kids)) {
    s <- subset
    s[kidids != kidid] <- FALSE
    # Node ID
    if (kidid > 1) {myid <- max(nodeids(kids[[kidid-1]]))
    } else {myid <- id}
    # Start recursion on this daugther node
    kids[[kidid]] <- growTemp(id=as.integer(myid+1), depth=depth, data=data, response=response, treatment=treatment, Propensity=Propensity,
                              subset=s, search=search, method=method, split=split, mtry=mtry, nsplit=nsplit, nsplit.random=nsplit.random,
                              minsplit=minsplit, minbucket=minbucket, maxdepth=maxdepth,
                              a=a, scale.y=scale.y, trtlevels=trtlevels)
  }
  return(partynode(id=as.integer(id), split=sp, kids=kids,
                   info=list(stats=max(info_split(sp)$stats, na.rm=TRUE))))
}

### Grows a tree with formula interface. Converts to constparty object ###
growTree <- function(formula, data, subset=NULL, search=c("exhaustive","sss"),
                     method=c("RCT","observation"),
                     split=c("t.test", "pvalue"),#, "gini", "entropy", "information"),
                     mtry=NULL, nsplit=NULL, nsplit.random=TRUE, minsplit=20, minbucket=round(minsplit/3), maxdepth=20,
                     a=50, useRes, scale.y=FALSE, trtlevels)
{
  search <- match.arg(search,c("exhaustive","sss"))
  method <- match.arg(method,c("RCT","observation"))
  split <- match.arg(split,c("t.test", "pvalue"))
  stopifnot(is.logical(nsplit.random), is.logical(scale.y), is.logical(useRes))#, is.logical(useRpart))
  if (is.numeric(nsplit) && !nsplit.random && nsplit < 5) {"Selecting <5 ordered splits may yield unexpected results"}
  
  response <- all.vars(formula)[1]
  if(grepl("\\|", as.character(formula)[3])){
    treatment <- trimws(strsplit(as.character(formula)[3], "\\|")[[1]][2], which="both")
  } else {stop("Please specify the treatment in formula")} 
  if(useRes) {oriResponse<- data$yo}
  
  Propensity <- "prop"
  Iptw <- "iptw"
  data <- data[c(all.vars(formula)[-1], response, Iptw,Propensity)] #Rearrange data so that response comes last
  
  if (!all(complete.cases(data)) & !is.null(subset)) { paste0("Specifying subset with missing data can yield unexpected results") }
  data <- data[complete.cases(data),] 
  
  
  if (is.null(mtry)){mtry <- length(all.vars(formula[[3]]))-1} #defult mtry=p
  
  #if(is.factor(data[[response]])){data[[response]]=as.numeric(data[[response]]==levels(data[[response]])[1])}
  
  if (is.null(subset)){subset <- rep(TRUE, nrow(data))}
  
  # Grow tree
  nodes <- growTemp(id=1L, depth=1L, data=data, response=response, treatment=treatment, Propensity=Propensity, subset=subset, search=search, method=method, split=split,
                    mtry=mtry, nsplit=nsplit, nsplit.random=nsplit.random, minsplit=minsplit, minbucket=minbucket,
                    maxdepth=maxdepth, a=a, scale.y=scale.y, trtlevels=trtlevels)
  
  # Compute terminal node number for each observation
  fitted <- fitted_node(nodes, data=data)
  
  if(ncol(data[[Propensity]])!=1) {
    daprop<-cbind(data[[treatment]],data[[Propensity]])
    ps<-apply(daprop,1,function(v){x<-v[1];return(v[x+1])})
  }
  ps<- data[[Propensity]]
  # Return rich constparty object
  ret <- party(nodes, data = data,
               fitted = data.frame("(fitted)" = fitted,
                                   "(response)" = data[[response]],
                                   "(treatment)" = data[[treatment]],
                                   "(propensity)"= ps,#data[[Propensity]],
                                   "(iptw)"= data[[Iptw]],
                                   check.names = FALSE),
               terms = terms(formula))
  as.constparty(ret)
  #as.simpleparty(ret)
}

##quantile function to trancate weight
truncquant <- function(prop,q=0.9){
  qmax <- quantile(prop,q)
  qmin <- quantile(prop,1-q)
  prop[prop >= qmax] <- qmax
  prop[prop<=qmin] <- qmin
  return(prop)
}

### Grows a random forest ###
CERFIT <- function( formula, data, ntrees, subset=NULL, search=c("exhaustive","sss"),
                    method=c("RCT","observation"), PropForm=c("randomForest","CBPS","GBM", "HI"),
                    split=c("t.test"),
                    mtry=NULL, nsplit=NULL, nsplit.random=TRUE, minsplit=20, minbucket=round(minsplit/3), maxdepth=30,
                    a=50, sampleMethod=c('bootstrap','subsample','subsampleByID','learning'),
                    useRes=TRUE, scale.y=FALSE)#
{
  sampleMethod <- match.arg(sampleMethod, c('bootstrap','subsample','subsampleByID','learning'))
  
  if (missing(formula)) stop("A formula must be supplied.", call. = FALSE) 
  if (missing(data)) data <- NULL 
  
  if(useRes){
    resformula<- as.formula(paste(all.vars(formula)[1], paste(all.vars(formula)[2:(length(all.vars(formula))-1)], collapse=" + "), sep=" ~ "))
    reslm<-lm(resformula,data)
    eres<-resid(reslm)
    data$yo <- data[[all.vars(formula)[1]]]
    data[[all.vars(formula)[1]]] <- eres
  }
  
  TrT<-data[all.vars(formula)[length(all.vars(formula))]]
  trt.length<-nrow(unique(TrT))
  if (trt.length<2) stop("Only one treatment?", call. = FALSE)
  trt.type <- ifelse(trt.length==2,"binary","multiple")
  trt.type <- ifelse(trt.length>10, "continuous", trt.type )
  trtlevels<-c(1:trt.length)
  print(trt.type)
  
  if(method=="observation"){
    propformula <- as.formula(paste(all.vars(formula)[length(all.vars(formula))], paste(all.vars(formula)[2:(length(all.vars(formula))-1)], collapse=" + "), sep=" ~ "))
    if(trt.type=="continuous"){
      if(PropForm=="CBPS"){
        propfun <- CBPS(propformula, data = data[,all.vars(formula)[-1]],ATT=FALSE,method = "exact")# 
        prop <- propfun$fitted.values
        Iptw <- propfun$weights
      } else if(PropForm=="HI") {
        propfun <- lm(propformula,data=data[all.vars(formula)[-1]])
        prt <- predict(propfun)
        sigm <- summary(propfun)$sigma
        prop <- dnorm(TrT,prt,sigm)
        modhi=lm(TrT~1)
        ps.num=dnorm((TrT-modhi$fitted)/(summary(modhi))$sigma,0,1)
        Iptw=ps.num/prop
      } 
    } else if(trt.type=="binary") {
      if (PropForm=="GBM") {
        propfun<-ps(propformula,data=data[,all.vars(formula)[-1]],interaction.depth = 4, stop.method = "es.max",estimand="ATE",verbose=FALSE,n.trees = 10000)
        prop<-propfun$ps
        Iptw<-get.weights(propfun,stop.method = "es.max",estimand="ATE")
      } else if (PropForm=="CBPS") {
        propfun <- CBPS(propformula, data = data[,all.vars(formula)[-1]],ATT=FALSE,method = "exact")# 
        prop <- propfun$fitted.values
        Iptw <- propfun$weights  
      } else if (PropForm=="randomForest") {
        propfun<- suppressWarnings(randomForest(propformula,data=data[all.vars(formula)[-1]]))
        prop <- propfun$predicted
        Iptw <-sum(Trt)/length(Trt)*Trt/prop+sum(1-Trt)/length(Trt)*(1-Trt)/(1-prop)
        #Iptw <-TrT/prop+(1-TrT)/(1-prop)
        Iptw<-truncquant(Iptw,q=0.9)
      }
    } else if (trt.type=="multiple"){
      if(PropForm=="GBM") {
        data[,all.vars(formula)[length(all.vars(formula))]]<-as.factor(data[,all.vars(formula)[length(all.vars(formula))]])
        propfun<-mnps(propformula,data=data[,all.vars(formula)[-1]],interaction.depth = 4, stop.method = "es.max",estimand="ATE",verbose=FALSE,n.trees = 10000)
        pslist<-propfun$psList
        prop<-matrix(NA,ncol=trt.length,nrow=nrow(data))
        for(i in 1:trt.length){
          prop[,i]<-unlist(pslist[[i]]$ps)
        }
        colnames(prop)<-levels(data[,all.vars(formula)[length(all.vars(formula))]])
        levels(data[,all.vars(formula)[length(all.vars(formula))]])<-c(1:trt.length)
        Iptw<-get.weights(propfun,stop.method = "es.max",estimand="ATE")
      } else if (PropForm=="CBPS" & trt.length<5 ) {
        data[,all.vars(formula)[length(all.vars(formula))]]<-as.factor(data[,all.vars(formula)[length(all.vars(formula))]])
        propfun <- CBPS(propformula, data = data[,all.vars(formula)[-1]],ATT=FALSE,method = "exact")# 
        prop <- propfun$fitted.values
        Iptw <- propfun$weights
        levels(data[,all.vars(formula)[length(all.vars(formula))]])<-c(1:trt.length)
      } 
    }  else stop("Please specify a propensity score method: randomForest or CBPS or GBM", call. = FALSE)
  } else if (method=="RCT") {
    prop <- rep("none",nrow(data)) # for observational no prop need
    Iptw<- rep(1,nrow(data))}
  data[,all.vars(formula)[length(all.vars(formula))]]<-as.numeric(as.character(data[,all.vars(formula)[length(all.vars(formula))]]))
  data$iptw <- Iptw
  data$prop <- prop
  #Construct random forest
  randFor <- lapply(1:ntrees,function(b){
    if(b%%10==0){print(paste0("Tree Number: ",b))}
    print(paste0("Tree Number: ",b))
    obs.b <- switch(sampleMethod,
                    bootstrap = sample.int(nrow(data), size=nrow(data), replace=TRUE, prob=data$iptw), #inverse weighting in boostrapping
                    subsample = sample.int(nrow(data), size=round(nrow(data)*0.632), replace=FALSE,prob=data$iptw), # stratified sampling
                    subsampleByID = {nIds <- length(unique(data[[idVar]]))
                    unlist(lapply(sample(unique(data[[idVar]]), size=round(nIds*0.632), replace=FALSE),
                                  function(x){which(data[[idVar]] == x)}))},
                    learning = 1:nrow(data))
    sample.b <- data[obs.b,]
    tree.b <- growTree(formula=formula, data=sample.b, subset=subset, search=search, method=method, split=split,
                       mtry=mtry, nsplit=nsplit, nsplit.random=nsplit.random, minsplit=minsplit, minbucket=minbucket, maxdepth=maxdepth, a=a,
                       scale.y=scale.y, useRes=useRes, trtlevels=trtlevels)#, useRpart=useRpart, minpvalue=minpvalue, corstr=corstr)
    list(tree=tree.b,cases=sort(unique(obs.b)))
  })
  return(randFor)
}


### Return predicted values from a tree. ###
predictTree <- function(tree, newdata=tree$data, gridval, LB, UB, ntrt,type="response",alpha ){
  da <- fitted(tree)
  colnames(da)[2:5]<-c("y","Trt","prop","ww")
  ufit<-sort(unique(da[["(fitted)"]]))
  nodesNewdata <- predict(tree, newdata=newdata, type="node")
  Y.min<-ifelse(min(da[,2])<0,2*min(da[,2]),min(da[,2])/2)
  Y.max<-ifelse(max(da[,2])<0,max(da[,2])/2,2*max(da[,2]))
  #if (ntrt<=10){ #for binary and multiple trt, ignore gridval
  if(ntrt<=10){
    pred<- lapply(split(da ,list(da[["(fitted)"]],da[,3])), function(da){
      ytemp<-try(weighted.mean(da[,2],da[,5],na.rm=T))
      if(inherits(ytemp,"try-error")) { 
        return(NA)} else {
          return(ytemp)}
    })
    nodepred<- cbind(ufit,t(matrix(unlist(pred), ncol = length(ufit), byrow = TRUE)))
  } else{ 
    pred<- lapply(split(da ,da[["(fitted)"]]), function(da){
      Trt<-da$Trt
      x<-cbind(Trt,Trt^2,Trt^3)
      lambdas <- 10^seq(5, -3, by = -.1)
      fit <- try(cv.glmnet(x, da$y,  family = "gaussian", alpha = alpha, lambda = lambdas,nfolds =10),silent=TRUE)
      if (inherits(fit, "try-error")){
        fit2<-try(cv.glmnet(x, jitter(da$y), family = "gaussian", alpha = alpha, lambda = lambdas, nfolds =10),silent=TRUE)
        if (inherits(fit2, "try-error")) {
          return(NA)} else {fit<-fit2}
      }
      
      Trt<-gridval
      ext<-Trt>max(da[,3])|Trt<min(da[,3]) 
      nd<-cbind(gridval,gridval^2,gridval^3)
      ytemp <- predict(fit, newx = nd, s=fit$lambda.min)
      ytemp[!ext]=ifelse(ytemp[!ext]>Y.max,Y.max,ytemp[!ext])##avoid extrem values
      ytemp[!ext]=ifelse(ytemp[!ext]<Y.min,Y.min,ytemp[!ext])#mean(da[,2])
      ytemp[ext]=ifelse(ytemp[ext]>Y.max,NA,ytemp[ext])
      ytemp[ext]=ifelse(ytemp[ext]<Y.min,NA,ytemp[ext])
      if (type!="opT") {
        return(ytemp)
      }else {
        #top<-gridval[which.max(ytemp)]
        #yop<-max(ytemp)
        lengthout<-5
        B <- seq(LB, UB, length.out=lengthout)
        opY<-Y.min; opTrt <- NA
        pref<-function(Trt){
          trtt<-cbind(Trt,Trt^2,Trt^3)
          yp<-predict(fit, newx = trtt, s=fit$lambda.min)
          return(yp)}
        for (b in 1:(lengthout-1)) {
          fit.tmp <-  suppressWarnings(optimize(pref, lower=B[b], upper=B[b+1], maximum=TRUE))
          if (is.na(fit.tmp$objective)) {
            opY<-opY
            opTrt<-opTrt
          } else {
            if (!is.nan(fit.tmp$objective) && fit.tmp$objective > opY && fit.tmp$objective < Y.max ) {
              opY <- fit.tmp$objective
              opTrt <- fit.tmp$maximum
            } 
          } }				
        return(cbind(top,yop))}
    })
    nodepred<- cbind(ufit,matrix(unlist(pred), ncol = length(pred[[1]]), byrow = TRUE))
  }
  if(type=="opT") ntrt<-2
  predictions<-as.data.frame(cbind(nodesNewdata,matrix(NA,ncol=ntrt,nrow=nrow(newdata))))
  predictions[,2:(ntrt+1)] <- nodepred[match(predictions$nodesNewdata,nodepred[,1]),2:(ntrt+1)]
  return(predictions[,2:(ntrt+1)])
}

### Return predicted values from a random forest ###
predict.CERFIT <- function(cerfit, data,newdata, gridval=NULL, prediction=c("overall","by iter"), type=c("response","ITE","node","opT"), alpha=0.5,useRse=TRUE){ 
  
  #Return prediction using all trees ("overall") or using first i trees ("by iter")
  prediction <- match.arg(prediction, c("overall","by iter"))
  
  type <- match.arg(type, c("response","ITE","node","opT"))
  
  cumMeanNA <- function(x){
    xTemp<-x;
    xTemp[is.na(xTemp)] <- 0
    cumsum(xTemp)/cumsum(!is.na(x))}
  #utrt<- sort(unique(c(fitted(cerfit[[1]]$tree)[,3],fitted(cerfit[[2]]$tree)[,3],fitted(cerfit[[3]]$tree)[,3])))
  formulaTree <- formula(cerfit[[1]]$tree$terms)
  treatment <- all.vars(formulaTree)[length(all.vars(formulaTree))]
  utrt<-sort(unique(data[[treatment]]))
  LB<-min(data[[treatment]])
  UB<-max(data[[treatment]])
  qu<-seq(LB,UB,length.out = 6)
  ## add a statement warning if gridvalue beyond the LB and UB
  ## should add warnings here if gridbalue beyond min or max utrt
  #ntrt <- length(utrt)
  # if grival is null, use the 10th quantile
  if(useRse){
    resformula<-  as.formula(paste(all.vars(formulaTree)[1], paste(all.vars(formulaTree)[2:(length(all.vars(formulaTree))-1)], collapse=" + "), sep=" ~ "))
    reslm<-lm(resformula,data)
    ylmp<-predict(reslm,newdata)
  } else {
    ylmp<-rep(0,nrow(newdata))
  }
  if(length(utrt)<=20){ ## if less than 20 unique treatments/levels using unique treatments
    ntrt=length(utrt)
    gridval<-utrt
  } else if(is.null(gridval)) { # if more than 20, and gridval is null, use percentiles at 5% increment
    gridval<-quantile(utrt, prob = seq(0, 1, length = 21))
    ntrt<-length(gridval)-1 
  } else {
    ntrt<-length(gridval)}
  print(gridval)
  
  if(type!="opT"){
    predictMat<-lapply(lapply(cerfit, "[[" , "tree"), predictTree, newdata=newdata,gridval=gridval,ntrt=ntrt,type=type,LB=LB,UB=UB,alpha=alpha)
    ypre<- do.call(cbind,predictMat)
    #yp<- lapply(1:ntrt,function(i,k) k[,seq(i, by = ntrt, length = NCOL(ypre) / ntrt)],k=ypre)
    ypre<- lapply(1:ntrt,function(i,k) k[,seq(i, NCOL(ypre), by = ntrt)], k=ypre)
    y.pre<- t(matrix(unlist(lapply(ypre,rowMeans,na.rm=TRUE)), ncol=NROW(newdata),byrow = TRUE))
    y.pre<-y.pre+ylmp
    #y.pre: by row observation, each column is the corresponding predition for 1 treatment.
  } else{
    predictMat<-lapply(lapply(cerfit , "[[" , "tree"), predictTree, newdata=newdata,gridval=gridval,ntrt=ntrt,type="opT", alpha=alpha, LB=LB,UB=UB,alpha=alpha)
    ntrt<-2
    ypre<- do.call(cbind,predictMat)
    ypre<- lapply(1:ntrt,function(i,k) k[,seq(i, NCOL(ypre), by = ntrt)], k=ypre)
    topt<-as.matrix(ypre[[1]])
    yopt<-as.matrix(ypre[[2]])
    y.opt<-rowMeans(yopt)+ylmp
    t.opt<-rowMeans(topt)
    y.pre<- cbind(t.opt,y.opt)
  }
  
  yname<-NA
  if (prediction=="overall") {
    if(type=="response") {
      resp <- y.pre
      yname<- paste("y=",gridval,sep="") 
      colnames(resp) <- yname
      return(resp)}
    if(type=="ITE") { #using the first level or smallest value as reference group
      yname<-paste("y",utrt,"-y",utrt[1],sep="") 
      ite<- y.pre-y.pre[,1]
      colnames(ite) <- c(yname)
      return(ite[,-1])
    }
    if(type=="opT") {
      yname<-c("opTreat","opResponse")
      opTY<-y.pre
      colnames(opTY) <- c(yname)
      return(opTY)
    }
  }
  else if(prediction=="by iter"){
    Ypre<-as.list(NA)
    for(i in 1: ntrt){
      Ypre[[i]]<-t(apply(ypre[[i]],1,cumMeanNA))
    }
    cumypre<-t(matrix(unlist(Ypre),ncol=NROW(newdata),byrow = TRUE))
    ntree<-length(cerfit)
    cumypre.l<- lapply(seq(1,(ntrt*ntree),by=ntree),function(i,k) k[,i:(i+ntree-1)], k=cumypre)
    print(cumypre.l)
    if(type=="response"){
      yname<-paste("ycum",utrt,sep="") 
      names(cumypre.l) <- yname  
      return(cumypre.l)}
    if(type=="ITE")  {
      cumite<-as.list(NA)
      for(i in 1:ntrt){
        cumite[[i]]<- cumypre.l[[i]]-cumypre.l[[1]]  
      }
      yname<-paste("ycum",utrt,"-ycum",utrt[1],sep="") 
      names(cumite)<-yname
      print(yname)
      return(cumite[[-1]])
    }}
  
}


### Minimal depth for variable importance ###
MinDepth <- function(cerfit){  # need to given number of levels if observation
  Term<-cerfit[[1]]$tree$terms
  dataTemp<-all.vars(Term[[3]])
  vars<-dataTemp[-length(dataTemp)]
  
  mindepth <- rep(0, length(vars))
  for (t in seq_along(cerfit)) {
    intNodes <- nodeids(cerfit[[t]]$tree)[-nodeids(cerfit[[t]]$tree, terminal = TRUE)]
    varsInTree <- vars[unique(unlist(nodeapply(cerfit[[t]]$tree, ids=intNodes, FUN=function(n){split_node(n)$varid})))]
    varsAtNode <- unlist(nodeapply(cerfit[[t]]$tree, ids=intNodes, FUN=function(n){split_node(n)$varid}))
    #Root node should be counted as 0
    depthAtNode <- table(unlist(lapply(intNodes, function(x) intersect(intNodes, nodeids(cerfit[[t]]$tree, from=x)))))-1
    treeDepth <- depth(cerfit[[t]]$tree)
    
    for (j in seq_along(vars)) {
      if (is.element(vars[j], varsInTree)) { #If variable is in tree
        mindepth[j]=mindepth[j]+min(depthAtNode[varsAtNode==j])
      } else {  #If variable not in tree, set mindepth to maximum depth+1
        mindepth[j]=mindepth[j]+treeDepth+1
      }
    }
  }
  mindepth <- mindepth/length(cerfit)
  names(mindepth) <- vars
  return(mindepth)
}

### need further checking on the stablity of the following functions for parallel computing function###
### Grows tree setup for parallel computing ###
CERFIT_Parallel_temp <- function(b, formula, data, search, method,
                                 split, mtry, nsplit, nsplit.random,
                                 minsplit, minbucket, maxdepth,
                                 a, sampleMethod,
                                 scale.y, useRes) {
  library(partykit)
  
  if (b %% 100 == 0) {print(paste0("Tree Number: ",b))}
  
  if (grepl("\\|", as.character(formula)[3])) {
    idVar <- trimws(strsplit(as.character(formula)[3], "\\|")[[1]][2], which="both")
  } else {idVar <- ""}
  
  obs.b <- switch(sampleMethod,
                  bootstrap = sample.int(nrow(data), size=nrow(data), replace=TRUE),
                  subsample = sample.int(nrow(data), size=round(nrow(data)*0.632), replace=FALSE),
                  subsampleByID = {nIds <- length(unique(data[[idVar]]))
                  unlist(lapply(sample(unique(data[[idVar]]), size=round(nIds*0.632), replace=FALSE),
                                function(x){which(data[[idVar]] == x)}))},
                  learning = 1:nrow(data))
  sample.b <- data[obs.b,]
  tree.b <- growTree(formula=formula, data=sample.b, search=search, method=method, split=split,
                     mtry=mtry, nsplit=nsplit, nsplit.random=nsplit.random, minsplit=minsplit, minbucket=minbucket, maxdepth=maxdepth, a=a,
                     scale.y=scale.y, useRes=useRes, trtlevels=trtlevels )
  list(tree=tree.b,cases=sort(unique(obs.b)))
}


### Grows a random forest using parallel computing: need further check on the function ###
### Note: not as many checks performed.  Recommend using CERFIT() function first and get code to work, ###
### then switch to CERFITparallel() function ###
CERFITparallel <- function(ntrees, formula, data, search, method,
                           split=c("MSE", "ttest", "gini", "entropy", "information"),
                           mtry, nsplit=NULL, nsplit.random=TRUE, minsplit, minbucket, maxdepth,
                           a=50, sampleMethod=c('bootstrap','subsample','subsampleByID','learning'),
                           scale.y=TRUE, useRes=useRes, iseed=111111){
  
  sampleMethod <- match.arg(sampleMethod, c('bootstrap','subsample','subsampleByID','learning'))
  
  #Construct random forest
  num_cores <- detectCores() - 1
  cl <- makeCluster(num_cores)
  clusterExport(cl,c("findCutpts","pickCutpt","splitrule","ttest.obs",
                     "ttest.rand","ordinalize","partition","growTemp",
                     "growTree","CERFIT_Parallel_temp"))
  clusterSetRNGStream(cl=cl, iseed=iseed)
  randFor <- parLapply(cl,1:ntrees, CERFIT_Parallel_temp,
                       formula=formula, data=data, search=search, method=method,
                       split=split, mtry=mtry, nsplit=nsplit, nsplit.random=nsplit.random,
                       minsplit=minsplit, minbucket=minbucket, maxdepth=maxdepth,
                       a=a, sampleMethod=sampleMethod,
                       scale.y=scale.y, useRes=useRes, trtlevels=trtlevels)
  stopCluster(cl)
  return(randFor)
}
