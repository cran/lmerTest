Dev <- function(rho, vec.matr, nll = FALSE) {
### Deviance of a LMM as a function of the variance-covariance
### parameters.
### nll: should the negative log-likelihood rather than the deviance
### (= 2*nll) be returned?  
  sigma<-vec.matr[[1]]
  Lambda<-makeLambda(rho, vec.matr) 
  Ut <- crossprod(Lambda,rho$Zt)
  L <- Cholesky(tcrossprod(Ut), LDL = FALSE, Imult = 1)
  cu <- solve(L, solve(L, Ut %*% rho$y, sys = "P"), sys = "L")
  RZX <- solve(L, solve(L, Ut %*% rho$X, sys = "P"), sys = "L")
  RX <- chol(rho$XtX - crossprod(RZX))
  cb <- solve(t(RX),crossprod(rho$X,rho$y)- crossprod(RZX, cu))
  beta <- solve(RX, cb)
  u <- solve(L,solve(L,cu - RZX %*% beta, sys="Lt"), sys="Pt")
  fitted <- as.vector(crossprod(Ut, u) + rho$X %*% beta)
  ## evaluate using dnorm?
  prss <- sum(c(rho$y - fitted, as.vector(u))^2)
  ## rho$prss <- prss
  n <- length(fitted); p <- ncol(RX)
  ## ML deviance:
  dev <- as.vector(n * log(2 * pi * sigma^2) + prss / sigma^2 +
                   c(2 * determinant(L)$modulus))
  if(rho$REML) ## REML deviance:
    dev <- dev + as.vector(c(2 * determinant(RX)$modulus) -
                           p * log(2 * pi * sigma^2))
  if(nll) ## return negative log-likelihood rather than deviance?
    dev <- dev/2
  return(as.vector(dev))
}




#### Rune's function ################################
getVcov <- function(rho, vec.matr) {
### get the variance-covariance matrix of the fixed effects
### parameters, beta at the values of vec.matr. 
### This implementation use a profiled formulation of the deviance. 
  sigma<-vec.matr[[1]]
  Lambda<-makeLambda(rho, vec.matr)
  Ut <- crossprod(Lambda,rho$Zt)
  rho$L <- Cholesky(tcrossprod(rho$Zt), LDL = FALSE, Imult = 1, super = TRUE)
  L <- update(rho$L, Ut, mult = 1)
  RZX <- solve(L, solve(L, Ut %*% rho$X, sys = "P"), sys = "L")
  RX <- chol(rho$XtX - crossprod(RZX))
  vcov <- sigma^2 * chol2inv(RX)
  return(vcov)
}

#### Rune's function ################################
Ct.rhbc <- function(rho, vec.matr, Lc) {
### Returns the [ind, ind] element of the variance-covariance matrix
### of the fixed effects parameters, beta evaluated at the value of
### vec.matr. The gradient of this wrt. vec.matr are needed to estimate the
### Satterthwaite's degrees of freedom for the t-statistic.
  
  vcov <- getVcov(rho, vec.matr)
  return(as.matrix(Lc %*% as.matrix(vcov) %*% t(Lc)))
}

###########################################################################
#calculates Lambda matrix (see lmer theory)
###########################################################################
makeLambda<-function(rho,vec.matr)
{  
  #if there is correlation between intercept and slope in random term
  if(rho$corr.intsl)
  {
     Lambda<-matrix(nrow=0,ncol=0)
     for(i in 1:length(rho$nlev))
     {
    
       lambda1<-vec.matr[which(rho$param$vec.num==i)+1]
       #the correlation between intercept and slope is present
       if(length(lambda1)>1)
       {
          #one random coefficient
          #ST.full<-c(lambda1[1:2],0,lambda1[length(lambda1)])
          #ST<-matrix(ST.full,nrow=2,ncol=2)
          #multiple random coefficients
          ST<-matrix(0,nrow=rho$param$STdim[i],ncol=rho$param$STdim[i])
          ST[lower.tri(ST, diag=TRUE)]<-lambda1
          Lambda <- bdiag(Lambda,kronecker(ST, Diagonal(rho$nlev[i])))
          # a new one
          #Lambda <- bdiag(Lambda,kronecker(Diagonal(rho$nlev[i]), ST))
          
       }
       else
       {
          Lambda <- bdiag(Lambda,kronecker(lambda1, Diagonal(rho$nlev[i])))
       }
     }
   }
   if(!(rho$corr.intsl))
   { 
     Lambda<-Diagonal(x=rep.int(vec.matr[-1],rho$nlev))
   }
   
   return(Lambda)
}

##########################################################################
# Check if the data is balanced with respect to factors in it ############ 
##########################################################################
isbalanced<-function(data)
{
   nvar<-dim(data)[2]
   data.fac<-data
   var.quant<-NULL
   for(i in 1:nvar)
   {
      if(!is.factor(data[,i]))
         var.quant<-c(var.quant,i)
   }
   data.fac<-data[,-var.quant]
   return(!is.list(replications(~ . , data.fac)))
}

##########################################################################

##########################################################################
# Create rho environmental variable of mixed model ####################### 
##########################################################################
rhoInit<-function(model)
{
   # creating rho
   rho <- new.env(parent = emptyenv()) # create an empty environment
   rho$y <- model@y                   # store arguments and derived values
   rho$X <- model@X
   chol(rho$XtX <- crossprod(model@X))       # check for full column rank

   rho$REML <-  model@dims['REML']
   rho$Zt <- model@Zt
   #rho$nlev <- sapply(model@flist, function(x) length(levels(factor(x))))
   rho$L <- Cholesky(tcrossprod(model@Zt), LDL = FALSE, Imult = 1, super = TRUE)
   ls.str(rho)
 
   # change rho$nlev to suit random coefficients
   rf.model<-ranef(model)
   rho$nlev<-NULL
   nlev.names<-NULL
   for(i in 1:length(rf.model))
   {   
       #nrow(rf.model[[i]])
       nlev.names<-c(nlev.names,rep(names(rf.model[i]),ncol(rf.model[[i]])))
       rho$nlev<-c(rho$nlev,rep(nrow(rf.model[[i]]),ncol(rf.model[[i]])))    
   }
   names(rho$nlev)<-	nlev.names 
   
    rho$s<-summary(model,"lme4")
   
   rho$fixEffs<-fixef(model)



   
   #correlation between intercept and slope is present
   #put all necessary info about correlation in param variable
   param<-NULL
   #add std dev to the vector
   #param$vec.matr<-as.numeric(rho$s@REmat[nrow(rho$s@REmat),4])
   # a new one
   param$vec.matr<-attr(VarCorr(model), "sc")
   param$vec.num<-NULL

   param$STdim<-NULL
   
   for(i in 1:length(model@ST)) 
   {
          
     #correlation between intercept and slope is present
     if(nrow(model@ST[[i]])>1)
     {       
       S<-diag(diag(model@ST[[i]]))       
       T<-model@ST[[i]]
       T[!lower.tri(model@ST[[i]])]<-0
       diag(T)<-1
       lambda1<-T %*% S
       #for one random coefficient
       #param$vec.matr<-c(param$vec.matr,as.vector(lambda1)[-3])
       #param$vec.num<-c(param$vec.num,rep(i,3))       
       # for one random coefficient
       #rho$nlev<-rho$nlev[-i]
       
       # for multiple random coefficients
       param$vec.matr<-c(param$vec.matr,lambda1[lower.tri(lambda1, diag=TRUE)])
       param$vec.num<-c(param$vec.num,rep(i,length(which(as.vector(lower.tri(lambda1, diag=TRUE))==TRUE))))       
       rho$nlev<-rho$nlev[-((i+1):(i+ncol(lambda1)-1))]
       param$STdim<-c(param$STdim,ncol(lambda1))
     }
     else
     {
       param$vec.matr<-c(param$vec.matr,model@ST[[i]])
       param$vec.num<-c(param$vec.num,i)
       param$STdim<-c(param$STdim,1)
     }       
   }
   
   #check if there are correlations between intercepts and slopes
   rho$corr.intsl<-checkCorr(model)
   
   rho$param<-param
   return(rho)

}

       
##############################################################################################
# function to calculate summary of F test with Satterthwaite's approximation of denominator df 
##############################################################################################
calcSatterth<-function(Lc, rho, method.grad)
{
  # F statistics for tested term
  if(is.vector(Lc))
     C.theta.optim<-Ct.rhbc(rho, rho$param$vec.matr, t(Lc))
  else
     C.theta.optim<-Ct.rhbc(rho, rho$param$vec.matr, Lc)
  #invC.theta<-ginv(C.theta.optim)
  invC.theta<-solve(C.theta.optim)
  q<-qr(C.theta.optim)$rank
  F.stat<-(t(Lc %*% rho$fixEffs) %*% invC.theta %*% (Lc %*% rho$fixEffs))/q
  
  
  #df for F statistics for tested term
  svdec<-eigen(C.theta.optim) 
  
  
  PL<-t(svdec$vectors) %*% Lc
  
  nu.m<-NULL
  for(m in 1:length(svdec$values))
  {   
     g<-grad(function(x)  Ct.rhbc(rho,x,t(PL[m,])), rho$param$vec.matr , method = method.grad)
     nu.m<-c(nu.m, 2*(svdec$values[m])^2/(t(g) %*% rho$A %*% g))
  }
  
  E<- sum( (nu.m/(nu.m-2)) * as.numeric(nu.m>2))
  nu.F<-2*E*as.numeric(E>q)/(E-q)
  
  pvalueF<-1-pf(F.stat,qr(Lc)$rank,nu.F)
  return(list(denom = nu.F, Fstat = F.stat, pvalue = pvalueF))

}

       
###########################################################################
# function to calculate F stat and pvalues for a given term
###########################################################################
calcFpvalue<-function(term, L, model, rho, ddf, method.grad="simple")
{

  #L<-calcGeneralSetForHypothesis(DesignMat, rho)
  Lc<-makeContrastType3SAS(model, term, L)
    
  
  if(is.null(Lc))
    return(NULL)
  
  # for running rune's vcov function
  if(is.vector(Lc))
  {
    #Lc<-Lc[which(rho$s.test!=0)]
    Lc<-Lc[rho$nums.Coefs]
  }  
  else
  {
     #Lc<-Lc[,which(rho$s.test!=0)]
     Lc<-Lc[,rho$nums.Coefs]
  }
   
  if(ddf=="Kenward-Roger")
  {
    if(is.vector(Lc))
      res.KR<-KRmodcomp(model,t(as.matrix(Lc)))
    else
      res.KR<-KRmodcomp(model,Lc)
   # return(list(denom = res.KR$stats["df2"], Fstat = res.KR$stats["Fstat"], pvalue =  res.KR$stats["p.value"]))
    return(list(denom = res.KR$test[1,"ddf"], Fstat = res.KR$test[1,"F.scaling"], pvalue =  res.KR$test[1,"p.value"]))
  }
  else
  {
    ## apply satterthwaite's approximation of ddf
    return(calcSatterth(Lc, rho, method.grad))
  }
}


###############################################################################
# function to calculate T test 
###############################################################################
calculateTtest<-function(rho, Lc, nrow.res, method.grad)
{
  #(Lc<-t(popMatrix(m, c("Product"))))
  #resultTtest <- matrix(0, nrow = ncol(Lc), ncol = 3)
  #define Lc contrast matrix for t-test
  #Lc <- diag(rep(1,nrow(rho$s@coefs)))
  resultTtest <- matrix(0, nrow = nrow.res, ncol = 4)
  colnames(resultTtest) <- c("df", "t value", "p-value", "sqrt.varcor")
  #rownames(resultTtest) <- rownames(rho$s@coefs)

  #
  for(i in 1:nrow.res)
  {
     g <- grad(function(x) Ct.rhbc (rho, x, t(Lc[,i])) , rho$param$vec.matr, method = method.grad)   
     #denominator df
     denom<-t(g) %*% rho$A %*% g
     varcor<-Ct.rhbc(rho, rho$param$vec.matr, t(Lc[,i]))
     #df
     resultTtest[i,1] <- 2*(varcor)^2/denom
     #statistics
     resultTtest[i,2] <- (Lc[,i] %*%rho$s@coefs[,1])/sqrt(varcor)
     resultTtest[i,3] <- 2*(1 - pt(abs(resultTtest[i,2]), df = resultTtest[i,1]))
     resultTtest[i,4] <- sqrt(varcor) 
  }
  
  return(resultTtest)
}

###############################################################################
# construct design matrix for F test 
###############################################################################
createDesignMat<-function(model,data)
{
model.term <- terms(model)
fixed.term <- attr(model.term,"term.labels") #parsed.formula$labels[parsed.formula$random==0]
X.design<-NULL
names.design<-NULL
for(i in 1:length(fixed.term))
{
   formula.term<-as.formula(paste("~", fixed.term[i], "- 1"))
   X.design<-cbind(X.design,model.matrix(formula.term, data))
   names.design<-c(names.design,rep(fixed.term[i],ncol(model.matrix(formula.term, data))))
}
X.design<-cbind(rep(1,dim(X.design)[1]),X.design)
names.design<-c("Intercept",names.design)
colnames(X.design)<-names.design
return(X.design)
}


###############################################################################
# initialize anova table for F test 
###############################################################################
initAnovaTable<-function(model, isFixReduce)
{
  anm<-anova(model, ddf="lme4")
  NumDF<-anm[,1]
  p.value<-DenDF<-F.value<-as.numeric(rep("",length(NumDF)))
  anova.table<-cbind(NumDF,DenDF,F.value,p.value)
  colnames(anova.table)<-c("NumDF","DenDF","F.value","Pr(>F)")
  rownames(anova.table)<-rownames(anm)
  if(isFixReduce)
  {
    if(nrow(anova.table)==1)
    {
      anova.table<-c(anova.table[,1:3],0,anova.table[,4])
      anova.table<-matrix(anova.table, nrow=1, ncol=length(anova.table))
      colnames(anova.table)<-c("NumDF","DenDF","F.value","elim.num","Pr(>F)")
      rownames(anova.table)<-rownames(anm)
      return(anova.table)
    }
    elim.num<-rep(0,nrow(anova.table))
    anova.table<-cbind(anova.table[,1:3],elim.num,anova.table[,4])
    colnames(anova.table)[5]<-"Pr(>F)"
  }
  return(anova.table)    
}


###############################################################################
# get terms contained 
###############################################################################
getIndTermsContained<-function(allterms, ind.hoi)
{
  
  terms.hoi.split<-strsplit(allterms[ind.hoi],":")
  ind.terms.contain<-NULL
  #check which of the terms are contained in the highest order terms
  for(i in (1:length(allterms))[-ind.hoi]) 
  {
    isContained<-FALSE
    for(j in 1:length(terms.hoi.split))
    {
      #if the term is contained in some of the highest order interactions then 
      #we cannot test it for significance
      if(length(which(unlist(strsplit(allterms[i],":")) %in% terms.hoi.split[[j]] == FALSE))==0)
      {
        isContained<-TRUE
        break
      }                
    }
    if(isContained)
      ind.terms.contain<-c(ind.terms.contain,i)
    
  }
  # if there are no terms that are contained in the maximum order effects
  # then compare all the terms between each other for the maximum p value
  if(is.null(ind.terms.contain))
    return(NULL)
  return(ind.terms.contain)
}

###############################################################################
# get terms to compare in anova.table
###############################################################################
getTermsToCompare<-function(model)
{
  order.terms<-attr(terms(model),"order")
  allterms<-attr(terms(model),"term.labels")
  ind.hoi<- which(order.terms == max(order.terms))
  ind.terms.contain<-getIndTermsContained(allterms, ind.hoi)
  
  #get the rest of the terms to compare
  allterms.rest<- allterms[-c(ind.terms.contain, ind.hoi)]
  if(length(allterms.rest)==0)
    terms.compare<-allterms[ind.hoi]
  else
  {
    #get highest order terms in the remaining ones
    order.rest<-unlist(lapply(allterms.rest, function(x) length(unlist(strsplit(x,":")))))
    ind.hoi.rest<- which(order.rest == max(order.rest))
    gtc<-getIndTermsContained(allterms.rest, ind.hoi.rest)
    if(!is.null(gtc))
      terms.compare<-c(allterms[ind.hoi],allterms.rest[-getIndTermsContained(allterms.rest, ind.hoi.rest)])
    else
      terms.compare<-c(allterms[ind.hoi],allterms.rest)
  }
  return(terms.compare)
}

###############################################################################
# find NS effect from the model (starting from highest order interactions)
###############################################################################
getNSFixedTerm<-function(model, anova.table, data, alpha)
{
  
  pv.max<-0
  
  terms.compare<-getTermsToCompare(model)
  
  for(tcmp in terms.compare)
  {
    ind<-which(rownames(anova.table)==tcmp)
    if(anova.table[ind,which(colnames(anova.table)=="Pr(>F)")]>=pv.max)
    {
      ns.term<-tcmp
      pv.max<-anova.table[ind,which(colnames(anova.table)=="Pr(>F)")]
    }
  }  
  if(pv.max>=alpha)
    return(ns.term)  
  else
    return(NULL)  
}
  
 
###############################################################################
# eliminate NS effect from the model
############################################################################### 
elimNSFixedTerm<-function(model, anova.table, data, alpha, elim.num, l)
{
  ns.term<-getNSFixedTerm(model, anova.table, data, alpha)
  if(is.null(ns.term))
    return(NULL)
  anova.table[which(rownames(anova.table) == ns.term),which(colnames(anova.table)=="elim.num")]<-elim.num
  fm <- formula(model)
  fm[3] <- paste(fm[3], "-", ns.term)
  mf.final<- as.formula(paste(fm[2],fm[1],fm[3], sep=""))
  #if(!is.null(l))
  #  model<-eval(substitute(lmer(mf.final, data=data, contrasts=l),list(mf.final=mf.final)))
  #else
  #  model<-eval(substitute(lmer(mf.final, data=data),list(mf.final=mf.final)))
  model<-updateModel(model, mf.final, model@dims[["REML"]], l)
  #model<-update(model,formula. = mf.final)
  return(list(model=model, anova.table=anova.table))
}
  
 
  
#################################################################
# find which effect contains effect term
#################################################################
relatives <- function(classes.term, term, names, factors)
{
  # checks if the terms have the same number of covariates (if any)
  checkCovContain<-function(term1, term2)
  {        
    num.numeric<-which(classes.term=="numeric")
    num.numeric.term1<-which((num.numeric %in% which(factors[,term1]!=0))==TRUE)
    num.numeric.term2<-which((num.numeric %in% which(factors[,term2]!=0))==TRUE)
    if((length(num.numeric.term1)>0 && length(num.numeric.term2)>0)||(length(num.numeric.term1)==0 && length(num.numeric.term2)==0))
       return(all(num.numeric.term2 == num.numeric.term1))
    else
       return(FALSE)
  }
  is.relative <- function(term1, term2) 
  {
    return(all(!(factors[,term1]&(!factors[,term2]))) && checkCovContain(term1,term2))
  }
  if(length(names) == 1) return(NULL)
  	 which.term <- which(term==names)
	  (1:length(names))[-which.term][sapply(names[-which.term], 
		  			function(term2) is.relative(term, term2))]
}

############################################################################       
# caclulate the General contrast matrix for the hypothesis (as in SAS)
############################################################################
calcGeneralSetForHypothesis<-function(X.design, rho)
{
  #zero out dependent columns: in order to calculate g2 inverse
  #X.design2<-X.design
  #X.design2[,which(rho$s.test==0)]<-rep(0,nrow(X.design2))
  #X.design2[,rho$nums.zeroCoefs]<-rep(0,nrow(X.design2))
  
  
  xtx<-t(X.design) %*% X.design
  #xtx2<-t(X.design2) %*% X.design2
  
  g2<-matrix(0,ncol=ncol(xtx), nrow=nrow(xtx))
  g2[rho$nums.Coefs,rho$nums.Coefs]<-solve(xtx[rho$nums.Coefs,rho$nums.Coefs])
  #g2<-ginv(xtx2)
  g2[abs(g2)<1e-10]<-0
  
  #check g2:
  #all.equal(xtx %*% g2 %*% xtx, xtx)
  ######all.equal(g2 %*% xtx %*% xtx, xtx)
  #all.equal(g2 %*% xtx %*% g2, g2)

  #general set of estimable function
  L<-g2 %*% xtx
  L[abs(L)<1e-6]<-0
  return(L)
}
     
       
############################################################################       
# type 3 hypothesis SAS
############################################################################
makeContrastType3SAS<-function(model, term, L)
{
  
  eps <- 1e-8
  #apply rule 1 (Goodnight 1976)
  
  #find all effects that contain term effect
  model.term <- terms(model)
  fac <- attr(model.term,"factors")
  names<-attr(model.term,"term.labels")
  classes.term<-attr(model.term,"dataClasses")
  
  cols.eff<-which(colnames(L)==term)
  num.relate<-relatives(classes.term,term,names,fac)
  if(length(num.relate)==0)
    colnums<-setdiff(1:ncol(L),cols.eff)
  if(length(num.relate)>0)
  {
    cols.contain<-NULL
    for(i in 1:length(num.relate))
      cols.contain<-c(cols.contain,which(colnames(L)==names[num.relate[i]]))
    colnums<-setdiff(1:ncol(L),c(cols.eff,cols.contain))   
  }
    
  for(colnum in colnums)
  {
    
    pivots<-which(abs(L[,colnum]) > eps)
    #pivots<-which(L[,colnum]!=0)
    if(length(pivots)>0)
    {
      L[pivots[1],]<-L[pivots[1],]/L[pivots[1],colnum]
      nonzeros<-setdiff(pivots,pivots[1])
      if(length(nonzeros)!=0)
      {
         for(nonzero in nonzeros)
         {
           L[nonzero,]<-L[nonzero,]-L[nonzero,colnum]*L[pivots[1],]
         }
      }
     
      L[pivots[1],]<-rep(0,ncol(L))
    }
  }
    
  nums<-which(apply(L,1,function(y) sum(abs(y)))!=0) 
  L<-L[nums,]
  
  if(is.vector(L))
    return(L)
  
  #orthogonalization
  if(length(cols.eff)>1)
      zero.rows<-which(apply(L[,cols.eff],1,function(y) sum(abs(y)))==0)
  else
      zero.rows<-which(L[,cols.eff]==0)
      
  for(zero.row in zero.rows) 
  {
    w<-L[zero.row,]
    for(i in setdiff(1:nrow(L),zero.row))
    {
      if(sum(abs(L[i,]))!=0)
        L[i,]<-L[i,]-((w %*% L[i,])/(w %*% w)) %*% w
    }
    L[zero.row,]<-rep(0,ncol(L))
  }

  L[abs(L)<1e-6]<-0
  
  nums<-which(apply(L,1,function(y) sum(abs(y)))!=0) 
  L<-L[nums,]
  return(L)
}

############################################################################
#get formula for model 
############################################################################
getFormula<-function(model, withRand=TRUE)
{
  fmodel<-formula(model)
  terms.fm<-attr(terms(fmodel),"term.labels")
  ind.rand.terms<-which(unlist(lapply(terms.fm,function(x) substring.location(x, "|")$first))!=0)
  terms.fm[ind.rand.terms]<-unlist(lapply(terms.fm[ind.rand.terms],function(x) paste("(",x,")",sep="")))
  fm<-paste(fmodel)
  if(withRand)
    fm[3]<-paste(terms.fm,collapse=" + ")
  else
    fm[3]<-paste(terms.fm[-ind.rand.terms],collapse=" + ")
  return(fm)
}


###################################################################
#get the combinatoion of the fixed factors for the lsmeans
###################################################################
getFacCombForLSMEANS<-function(split.eff, data)
{
  if(length(split.eff)==1)
    data.merge<-as.data.frame(levels(data[,split.eff]))
  if(length(split.eff)>=2)
    data.merge<-merge(levels(data[,split.eff[1]]),levels(data[,split.eff[2]]))
  if(length(split.eff)>=3)
  {
    for(i in 3:length(split.eff))
    {
      d.split.eff_i<-as.data.frame(levels(data[,split.eff[i]]))
      names(d.split.eff_i)<-paste("l",i)
      data.merge<-merge(data.merge,d.split.eff_i)
    }
              
  }
  names(data.merge)<-split.eff
  return(as.matrix(data.merge))
}


###################################################################
#checks if all the terms in interaction are covariates
###################################################################
checkAllCov<-function(split.eff, data)
{
  for(spleff in split.eff)
  {
    if(!is.factor(data[,spleff]))
    {
      return(TRUE)  
    }
  }
  return(FALSE)
}

###################################################################
#concatenate levels of the effects to form the rownames
###################################################################
concatLevs<-function(matr, row.names)
{
  
  if(is.vector(matr))
    levs<-paste(names(matr),matr)
  else
  {
    levs<-paste(rownames(matr),matr[,1])
    for(i in 2:ncol(matr))
    {
      levs<-paste(levs,matr[,i])
    }    
  }
  
    
  return(levs)
}
    

#convert facs into numeric
convertFacsToNum<-function(data, begin, end)
{
  
 #convert vars to numeric
 for(i in begin:end)
  data[,i]<-as.numeric(levels(data[,i])[as.integer(data[,i])]) 

  return(data)
}

#convert numeric to facs
convertNumsToFac<-function(data, begin, end)
{
  
 #convert vars to numeric
 for(i in begin:end)
  data[,i]<-as.factor(data[,i]) 

  return(data)
}

###################################################################
#fills the LSMEANS and DIFF summary matrices
###################################################################
fillLSMEANStab<-function(mat, rho, summ.eff, nfacs, alpha, method.grad)
{
   estim.lsmeans<-mat%*%rho$fixEffs
   summ.eff[,nfacs+1]<-estim.lsmeans
   ttest.res<-calculateTtest(rho, t(mat), nrow(mat), method.grad)
   summ.eff[,nfacs+2]<-ttest.res[,4]#stdErrLSMEANS(rho, std.rand, mat)
   #df
   summ.eff[,(nfacs+3)]<-ttest.res[,1]
   #t values
   summ.eff[,(nfacs+4)]<-ttest.res[,2]
   #p values
   summ.eff[,(nfacs+7)]<-ttest.res[,3]
   # CIs
   summ.eff[,nfacs+5]<-estim.lsmeans-abs(qt(alpha/2,ttest.res[,1]))*ttest.res[,4]
   summ.eff[,nfacs+6]<-estim.lsmeans+abs(qt(alpha/2,ttest.res[,1]))*ttest.res[,4]
   return(summ.eff)
}

###################################################################
#round the columns of LSMEANS or DIFFS tables
###################################################################
roundLSMEANStab<-function(summ.eff, nfacs)
{
   summ.eff[,nfacs+1]<-round(summ.eff[,nfacs+1],4)  
   summ.eff[,nfacs+2]<-round(summ.eff[,nfacs+2],4)#stdErrLSMEANS(rho, std.rand, mat)
   #df
   summ.eff[,(nfacs+3)]<-round(summ.eff[,(nfacs+3)],1)
   #t values
   summ.eff[,(nfacs+4)]<-round(summ.eff[,(nfacs+4)],2)
   #p values
   summ.eff[,(nfacs+7)]<-round(summ.eff[,(nfacs+7)],4)
   # CIs
   summ.eff[,nfacs+5]<-round(summ.eff[,nfacs+5],4)
   summ.eff[,nfacs+6]<-round(summ.eff[,nfacs+6],4)
   return(summ.eff)
}

############################################################################
#function to identify the colors of bar according to significance of effects
############################################################################
calc.cols<-function(x)
{
  if(x<0.001) 
    return("red") 
  if(x<0.01) 
    return("orange") 
  if(x<0.05) 
    return("yellow") 
  return("grey")
}

#get names for ploting barplots for the effects
getNamesForPlot<-function(names, ind)
{
  namesForPlot<-unlist(lapply(names, function(y) substring2(y, 1, substring.location(y, " ")$first[1]-1)))
  namesForLevels<-unlist(lapply(names, function(y) substring2(y, substring.location(y, " ")$first[1]+ind, nchar(y))))
  return(list(namesForPlot=namesForPlot, namesForLevels=namesForLevels))
}

#plots for LSMEANS or DIFF of LSMEANS
plotLSMEANS<-function(table, response, which.plot=c("LSMEANS", "DIFF of LSMEANS"))
{
    if(which.plot=="LSMEANS")
      names<-getNamesForPlot(rownames(table),2)
    else
      names<-getNamesForPlot(rownames(table),1)
    namesForPlot<-names$namesForPlot
    namesForLevels<-names$namesForLevels
    un.names<-unique(namesForPlot)
    
    
    for(i in 1:length(un.names))
    {
      inds.eff<-namesForPlot %in% un.names[i]
      split.eff <- unlist(strsplit(un.names[i],":"))
      col.bars<- lapply(table[inds.eff,][,"p-value"], calc.cols)
      #windows()
      #par(mfrow=c(1,1))
      x11()
      layout(matrix(c(rep(1,3),2,rep(1,3),2), 2, 4, byrow = TRUE))
      barplot2(table[inds.eff,"Estimate"],col=unlist(col.bars), ci.l=table[inds.eff,ncol(table)-2], ci.u=table[inds.eff,ncol(table)-1], plot.ci=TRUE, names.arg=namesForLevels[inds.eff], xlab=un.names[i], ylab=response, main=paste(which.plot," and CI plot for", un.names[i]))
      plot.new()
      legend("topright", c("ns","p<0.05", "p<0.01", "p<0.001"), pch=15, col=c("grey","yellow","orange","red"), title="SIGNIFICANCE", bty="n", cex=0.8)
      if(which.plot=="LSMEANS")
      {
        if(length(split.eff)==2)
        {
          par(mfrow=c(1,1))
          #windows()
          interaction.plot(table[inds.eff,split.eff[1]], table[inds.eff,split.eff[2]], table[inds.eff,"Estimate"], xlab=split.eff[1], ylab=response, trace.label=paste(split.eff[2]), main="2-way Interaction plot", col=1:nlevels(table[inds.eff,split.eff[2]]))
        }
      }             
    }
}


#calculate DIFFERENCES OF LSMEANS and STDERR for effect
calcDiffsForEff<-function(facs, fac.comb, split.eff, eff, effs, data, rho, alpha, mat, method.grad)
{
   ###calculating diffs for 2 way interaction
   if(length(split.eff)>=1 && length(split.eff)<=2)
   {   
     if(length(split.eff)==2)
     {            
       fac.comb.names<-concatLevs(fac.comb)
       main.eff<-effs[effs %in% split.eff]
       ### check if the difference between one factor within another is required
       #if(length(main.eff)==1)
       #{
       #  levs.main.eff<-levels(data[,main.eff])
         #fac.comb[,main.eff]
         #ttt<-unlist(lapply(fac.comb[,effs[effs %in% split.eff]], function(x) which(levels(data[,effs[effs %in% split.eff]])==x)))
       #  mat.names.diffs<-NULL        
       #  for(dupls in levs.main.eff)
       #  {
       #    ind.dupls<-which(fac.comb[,main.eff]==dupls)
        #   mat.names.diffs<-cbind(mat.names.diffs,combn(fac.comb.names[ind.dupls],2))
       #  }
       #  mat.nums.diffs<-apply(mat.names.diffs, c(1,2), function(x) which(fac.comb.names==x))
      # }
      # else
      # {
          mat.names.diffs<-combn(fac.comb.names,2)
          mat.nums.diffs<-apply(mat.names.diffs, c(1,2), function(x) which(fac.comb.names==x))
      # }     
     }
     else
     {
       mat.names.diffs<-combn(fac.comb,2)
       mat.nums.diffs<-apply(mat.names.diffs, c(1,2), function(x) which(fac.comb==x))
     }
    
     mat.diffs<-matrix(0, nrow=ncol(mat.nums.diffs), ncol=ncol(mat))
     colnames(mat.diffs)<-colnames(mat)
     for(ind.diffs in 1:ncol(mat.nums.diffs))
     {
       mat.diffs[ind.diffs,]<-mat[mat.nums.diffs[1,ind.diffs],]- mat[mat.nums.diffs[2,ind.diffs],]
     }
     names.combn<-apply(mat.names.diffs, 2, function(x) paste(x[1],x[2],sep="-"))
     rownames(mat.diffs)<-paste(eff, names.combn)
     
     diffs.summ<- matrix(NA, ncol=7, nrow=nrow(mat.diffs))
     colnames(diffs.summ)<-c("Estimate","Standard Error", "DF", "t-value", "Lower CI", "Upper CI", "p-value")
     rownames(diffs.summ)<-rownames(mat.diffs)
          
     diffs.summ<-as.data.frame(fillLSMEANStab(mat.diffs, rho, diffs.summ, 0, alpha, method.grad))
     return(roundLSMEANStab(diffs.summ, 0))
    }
   
}

#calculate LSMEANS and STDERR for effect
calcLsmeansForEff<-function(lsmeans.summ, fac.comb, eff, split.eff, alpha, mat, rho, facs, method.grad)
{
   
   summ.eff<-matrix(NA, ncol=ncol(lsmeans.summ), nrow=nrow(fac.comb))
   colnames(summ.eff)<-colnames(lsmeans.summ)
   #rownames(summ.eff)<-rep(eff, nrow(fac.comb))
   summ.eff[,split.eff]<-fac.comb
   names.arg<-concatLevs(summ.eff[,split.eff])
   summ.eff<-as.data.frame(fillLSMEANStab(mat, rho, summ.eff, length(facs), alpha, method.grad))
   summ.eff<-convertFacsToNum(summ.eff, length(facs)+1, ncol(summ.eff))
   #estim.lsmeans<-mat%*%rho$fixEffs
   #summ.eff[,length(facs)+1]<-round(estim.lsmeans,4)
   #ttest.res<-calculateTtest(rho, t(mat), nrow(mat))
   #summ.eff[,length(facs)+2]<-round(ttest.res[,4],4)#stdErrLSMEANS(rho, std.rand, mat)
   #df
   #summ.eff[,(length(facs)+3)]<-round(ttest.res[,1],1)
   #t values
   #summ.eff[,(length(facs)+4)]<-round(ttest.res[,2],2)
   #p values
   #summ.eff[,(length(facs)+7)]<-round(ttest.res[,3],4)
   # CIs
   #summ.eff[,length(facs)+5]<-round(estim.lsmeans-abs(qt(alpha/2,ttest.res[,1]))*ttest.res[,4],4)
   #summ.eff[,length(facs)+6]<-round(estim.lsmeans+abs(qt(alpha/2,ttest.res[,1]))*ttest.res[,4],4)
   
   summ.eff<-roundLSMEANStab(summ.eff, length(facs))
   
   #summ.eff.data<-as.data.frame(summ.eff, row.names="")
   #summ.eff.data<-convertFacsToNum(summ.eff.data, length(facs)+1)   
   rownames(summ.eff)<-paste(rep(eff, nrow(fac.comb)), names.arg)
   return(summ.eff) 
}


###################################################################
#calculate LSMEANS DIFFS and CI for all effects
###################################################################
calcLSMEANS<-function(model, data, rho, alpha, test.effs = NULL, method.grad="Richardson", lsmeansORdiff=TRUE, l)
{  
 
 #library(gplots)
 ####old code#########################
 #fm<-getFormula(model, withRand=FALSE)
 #if(fm[3]=="")
 #   m<-lm(as.formula(paste(fm[2],fm[1],1, sep="")), data=data)
 #else
 #   m<-lm(as.formula(paste(fm[2],fm[1],fm[3], sep="")), data=data)
 #####################################
 m<-lm(model, data=model$data, contrasts=l)
 #m<-lm(model, data=model$data)
 effs<-attr(terms(m),"term.labels")
 if(!is.null(test.effs))
    effs<-effs[effs %in% test.effs]
 dclass<-attr(terms(m),"dataClasses")
 facs<-names(dclass[which(dclass=="factor")])
 #Get standard deviation of random parameters from model
 std.rand<-as.numeric(rho$s@REmat[,3])

 
 #init lsmeans summary
 if(lsmeansORdiff)
 {
   lsmeans.summ<- matrix(ncol=length(facs)+7,nrow=0)
   colnames(lsmeans.summ)<-c(facs,"Estimate","Standard Error", "DF", "t-value", "Lower CI", "Upper CI", "p-value")
   summ.data<-as.data.frame(lsmeans.summ)
 }
 else
 {
   #init diff summary
   diff.summ<- matrix(ncol=7,nrow=0)
   colnames(diff.summ)<-c("Estimate","Standard Error", "DF", "t-value", "Lower CI", "Upper CI", "p-value")
   summ.data<-as.data.frame(diff.summ)
 }
 
  
 for(eff in effs)
 {
   
   split.eff <- unlist(strsplit(eff,":"))
   if(checkAllCov(split.eff, data))
     next
   mat <- popMatrix(m, split.eff)
   fac.comb <- getFacCombForLSMEANS(split.eff, data)  
   #if(length(split.eff)>=2)
   #  par(mfrow=c(1,1))
   #fill
   if(!lsmeansORdiff)
     summ.data <- rbind(summ.data,   calcDiffsForEff(facs, fac.comb, split.eff, eff, effs, data, rho, alpha, mat, method.grad))
   else
     summ.data <- rbind(summ.data,   calcLsmeansForEff(lsmeans.summ, fac.comb, eff, split.eff, alpha, mat, rho, facs, method.grad))
 }
 return(list(summ.data = summ.data))
}

#check if there are correlations between intercepts and slopes
checkCorr<-function(model)
{
   corr.intsl<-FALSE
   lnST<-length(model@ST)
   for(i in 1:lnST)
   {    
      if(nrow(model@ST[[i]])>1)
         corr.intsl<-TRUE
   } 
   return(corr.intsl) 
}

# get dummy coefficients of the fixed part of the model
getNumsDummyCoefs2<-function(model, data)
{
  ###old code #############################
  #fm<-getFormula(model, withRand=FALSE)
  #if(fm[3]=="")
  #   m<-lm(as.formula(paste(fm[2],fm[1],1, sep="")), data=data)
  #else
  #   m<-lm(as.formula(paste(fm[2],fm[1],fm[3], sep="")), data=data) 
  m<-lm(model, data=model$data)
  
  #get full coefficients
  dc<-dummy.coef(m)
  names.dc<-names(dc)[1]
  for (i in 2:length(dc))
  {
    # if the terms are covariates
    if(is.null(names(dc[[i]])) || names(dc)[i]==names(dc[[i]])[[1]])
      names.dc<-c(names.dc,names(dc)[i])
    else
    {
      #check the presence of covariates in terms
      effsTerm<-unlist(strsplit(names(dc)[i],":"))
      ln.names.dc.i<-length(effsTerm)
      ln.names.dc.i.1<-length(unlist(strsplit(names(dc[[i]])[[1]], ":")))
      # the covariances are present in interaction
      is.interact<-substring.location(names(dc)[i],":")$last!=0
      if((ln.names.dc.i!=ln.names.dc.i.1) && is.interact)
      {
        covs<-paste(effsTerm[1:(ln.names.dc.i-ln.names.dc.i.1)], collapse=":")      
        facsTerm<-unlist(lapply(strsplit(names(dc[[i]]), ":"), function(x) paste(paste(effsTerm[(ln.names.dc.i-ln.names.dc.i.1+1):ln.names.dc.i],x, sep=""),collapse=":")))
        names.dc<-c(names.dc, unlist(lapply(facsTerm, function(x) paste(c(covs, x), collapse=":"))))    
      }
      else
      {
        names.dc<-c(names.dc,unlist(lapply(strsplit(names(dc[[i]]), ":"), function(x) paste(paste(unlist(strsplit(names(dc)[i],":")),x, sep=""),collapse=":"))))  
      }    
    }
  }
  
  fullCoefs<-unlist(dc)
  names(fullCoefs)<-names.dc
  is.zeroCoef<-names(fullCoefs) %in% names(coef(m))
  return(list(nums.zeroCoefs = which(is.zeroCoef==FALSE), nums.Coefs = which(is.zeroCoef==TRUE)))
}


# get dummy coefficients of the fixed part of the model
getNumsDummyCoefs<-function(model, data, l)
{
  ### old code ######################
  #fm<-getFormula(model, withRand=FALSE)
  #if(fm[3]=="")
  #   m<-lm(as.formula(paste(fm[2],fm[1],1, sep="")), data=data)
  #else
  #   m<-lm(as.formula(paste(fm[2],fm[1],fm[3], sep="")), data=data) 
  ####################################
  
  m<-lm(model, data=model$data, contrasts=l)  
  #m<-lm(model, data=model$data)  
  
  #get full coefficients
  dc<-dummy.coef(m)
  zeroCoefs<-which(unlist(dc)==0)
  nonzeroCoefs<-which(unlist(dc)!=0)
  return(list(nums.zeroCoefs = zeroCoefs, nums.Coefs = nonzeroCoefs))
}



##############################################################################################################
## functions for popMatrix for LSMEANS (from doBy package)
##############################################################################################################
.get_xlevels <- function(obj){
  UseMethod(".get_xlevels")
}

.get_xlevels.default <- function(obj){
	obj$xlevels
}


.covariateAve <- function(object, at=NULL, tt=terms(object)){
  tt  <- delete.response(tt)
  att <- attributes(tt)
  rhs.terms <- rownames(att$factors)[rowSums(att$factors)>0]
  rhs.class <- att$dataClass[match(rhs.terms, names(att$dataClass))]
  nums      <- rhs.terms[rhs.class=="numeric"]

  ans  <- lapply(model.frame(object)[,nums, drop=FALSE], mean) 
  
  nn <- match(names(ans), names(at))
  nn <- nn[!is.na(nn)]
  at.num <- at[nn]
  ans[names(at[nn])] <- at.num
  attr(ans, "at.num") <- at.num
  ans
}


.get_vartypes <- function(object){
  tt <- terms(object)
  tt  <- delete.response(tt)
  att <- attributes(tt)
  rhs.terms <- rownames(att$factors)[rowSums(att$factors)>0]
  rhs.class <- att$dataClass[match(rhs.terms, names(att$dataClass))]
  nums      <- rhs.terms[rhs.class=="numeric"]
  fact      <- rhs.terms[rhs.class=="factor"]
  list(numeric=nums, factor=fact)
}


.set_xlevels <- function(xlev, at){
  nam    <- names(xlev)
  nn <- match(nam, names(at))
  nn <- nn[!is.na(nn)]
  at.fact <- at[nn]
  xlev[names(at[nn])]  <- at.fact
  attr(xlev, "at.fact") <- at.fact
  xlev
}



.getX <- function(object, newdata){
  tt <- terms(object)
  Terms  <- delete.response(tt)
  mf  <- model.frame(Terms, newdata, xlev = .get_xlevels(object))
  X   <- model.matrix(Terms, mf, contrasts.arg = .get_contrasts(object))
  attr(X,"assign")<-NULL
  attr(X, "contrasts") <- NULL
  X
}


.get_contrasts <- function(obj){
  UseMethod(".get_contrasts")
}

.get_contrasts.default <- function(obj){
  obj$contrasts
}


popMatrix <- function(object, effect=NULL, at=NULL, only.at=TRUE){
  tt <- terms(object)
  Terms   <- delete.response(tt)
  xlev    <- .get_xlevels(object)
  ccc     <- .covariateAve(object,at)
  vartype <- .get_vartypes(object)


##   cat("INPUT: effect:\n"); str(effect)
##   cat("INPUT: at:\n"); str(at)
##   cat("---------------------------\n")
  xlev   <- .get_xlevels(object)

  if (is.null(effect)){
    at.factor <- at[intersect(vartype$factor, names(at))]
    xxx       <- if(length(at.factor)>0)
      at.factor
  } else {
    xlev   <- .set_xlevels(xlev, at=at)
    at.fact <- names(attr(xlev, "at.fact"))
    effect <- setdiff(effect, at.fact)
    xxx    <- xlev[c(effect,at.fact)]
  }

#  print(ccc)
#  print(xxx)
  

  #print(xxx)
  if (is.null(xxx)){
    ## No 'effect' and no 'at'; just to a global average.
    newdata <- expand.grid(xlev)
    newdata[,names(ccc)] <- ccc   
    mf  <- model.frame(Terms, newdata, xlev = .get_xlevels(object))
    X   <- model.matrix(Terms, mf, contrasts.arg = .get_contrasts(object))
    res <- apply(X,2,mean)
    res <- do.call(rbind, list(res))
    attr(res,"at") <- at[intersect(vartype$numeric, names(at))]
  } else {
    eff.grid  <- expand.grid(xxx)
    eff.grid  <- as.data.frame(lapply(eff.grid, as.character),stringsAsFactors=FALSE)
    #cat("eff.grid:\n"); print(eff.grid)
    res <- list()
    for (ii in 1:nrow(eff.grid)){
      conf  <- eff.grid[ii,,drop=FALSE]
      xlev2 <- .set_xlevels(xlev,  at=conf)
      #cat("xlev2 (which defines the grid):\n"); str(xlev2)
      newdata <- expand.grid(xlev2)
      newdata[,names(ccc)] <- ccc   

      #print(newdata)
      mm   <- .getX(object, newdata)
      X    <- apply(mm,2,mean)
      res[[ii]] <- X
    }

    res <- do.call(rbind, res)
#    print(eff.grid)
    uuu <- at[intersect(vartype$numeric, names(at))]
#    print(uuu)
#    print(vartype)
#    print(at)
#    print(ccc)
    #eff.grid[,names(ccc)] <- at[intersect(vartype$numeric, names(at))]
    eff.grid[,names(ccc)] <- ccc
    attr(res,"grid") <- eff.grid
    attr(res,"at") <- at
  }
  class(res) <- c("popMatrix", "conMatrix","matrix")
  res 
}

#isWellSpecifiedModel<-function(model)
#{
#  if(sum(anova(model)$Df)!=dim(model@X)[2])
#    return(FALSE)
#  return(TRUE)
#}

emptyAnovaLsmeansTAB<-function()
{
  result<-NULL
  anova.table<- matrix(ncol=5,nrow=0)
  colnames(anova.table)<-c("Estimate","Standard Error", "DF", "F-value", "p-value")
  anova.table<-as.data.frame(anova.table)
  result$TAB.fixed<-anova.table
  lsmeans.summ<- matrix(ncol=7,nrow=0)
  colnames(lsmeans.summ)<-c("Estimate","Standard Error", "DF", "t-value", "Lower CI", "Upper CI", "p-value")
  lsmeans.summ<-as.data.frame(lsmeans.summ)
  result$TAB.lsmeans<-lsmeans.summ
  return(result)
}

#function checks if the model with covariates is well defined
#checkModelWithCovsTEST<-function(model)
#{
#  tt<-delete.response(terms(model))
#  num.ord<-which(attr(tt, "order")>1)
#  effs<-attr(tt,"term.labels")
#  for(nums in num.ord)
#  {
#    covs<-unlist(lapply(unlist(strsplit(effs[nums],":")), function(x) names(which(attr(tt,"dataClasses")[x]=="numeric"))))    
#    if(length(covs)!=0)
#    {
#      combs.lo<-unlist(lapply(combn(unlist(strsplit(effs[nums],":")),2, simplify=FALSE), function(x) paste(x,collapse=":")))
#      if(combs.lo %in% effs)    
#    }
#  }
#}


### initialize table for random terms
initRandTable<-function(terms, reduce.random)
{
  if(reduce.random)
  {
    rand.table<-matrix(0,ncol=4, nrow=length(terms))
    colnames(rand.table)<-c("Chi.sq", "Chi.DF", "elim.num", "p.value")
    rownames(rand.table)<-terms
  }
  else
  {
    rand.table<-matrix(0,ncol=3, nrow=length(terms))
    colnames(rand.table)<-c("Chi.sq", "Chi.DF", "p.value")
    rownames(rand.table)<-terms
  }
  return(rand.table)
}

### get names of terms out of rownames of rand.table
getTermsRandtable<-function(names.rand.table)
{
  return(unlist(lapply(names.rand.table, function(x) substring2(x,substring.location(x,"(")$first, nchar(x)))))
}


### fill a row for the random matrix
fillRowRandTable<-function(term, rand.table, rand.terms.upd=NULL, elim.num, reduce.random)
{
  nrow.term<-which(getTermsRandtable(rownames(rand.table))==term$term)
  
  rand.table[nrow.term, "Chi.sq"]<-term$chisq
  rand.table[nrow.term, "Chi.DF"]<-term$chisq.df
  rand.table[nrow.term, "p.value"]<-term$pv
  if(reduce.random)
    rand.table[nrow.term, "elim.num"]<-elim.num 
  if(!is.null(rand.terms.upd) && length(rand.terms.upd)!=0)
  {     
    rand.table.upd<-matrix(0,ncol=4, nrow=length(rand.terms.upd))
    colnames(rand.table.upd)<-c("Chi.sq", "Chi.DF", "elim.num", "p.value")
    #rownames(rand.table.upd)<-paste(paste(rep(" ",max(nchar(rownames(rand.table)))), collapse=""), rand.terms.upd, sep="")
    #rownames(rand.table.upd)<-rand.terms.upd
    nspace.term<-nchar(substring2(rownames(rand.table)[nrow.term],1,substring.location( rownames(rand.table)[nrow.term],"(")$first))
    rownames(rand.table.upd)<-paste(paste(rep(" ", nspace.term + 5), collapse=""), rand.terms.upd, sep="")
    if(nrow.term==nrow(rand.table))
    {
      rand.table<-rbind(rand.table, rand.table.upd)
      #rownames(rand.table)[1]<-term$term
     }
    else
    {
      rnames<-c(rownames(rand.table)[1:nrow.term], rownames(rand.table.upd),rownames(rand.table)[(nrow.term+1):nrow(rand.table)])
      rand.table<-rbind(rand.table[1:nrow.term,], rand.table.upd, rand.table[(nrow.term+1):nrow(rand.table),])  
      rownames(rand.table)<-rnames 
    } 
  }
  return(rand.table)
}

### update table for random terms
updateRandTable<-function(infoForTerm, rand.table, rand.terms.upd=NULL, elim.num=0, reduce.random)
{
  
  if(!is.null(infoForTerm$term))
  {   
    rand.table<-fillRowRandTable(infoForTerm, rand.table, rand.terms.upd, elim.num, reduce.random)    
    return(rand.table)
  }    
  else
  {
    for(iterm in infoForTerm)
      rand.table<-fillRowRandTable(iterm, rand.table, rand.terms.upd, elim.num, reduce.random)  
  }    
  return(rand.table)
}

############################################################################    
#save info (pvalues, chisq val, std, var...) for term
############################################################################    
saveInfoForTerm<-function(term, chisq, chisq.df, pv)
{
  term.info<-NULL 
  term.info$pv<-pv
  term.info$chisq<-chisq
  term.info$chisq.df<-chisq.df
  term.info$term<-term
  return(term.info)
}

### check if there are no random terms in the model
checkPresRandTerms<-function(mf.final)
{
  sub.loc.rand<-substring.location(paste(mf.final)[3], "|")
  if(sub.loc.rand$first==0 && sub.loc.rand$last==0)
    return(FALSE)
  return(TRUE)
}

### compare mixed model versus fixed
compareMixVSFix<-function(model, mf.final, data, name.term, rand.table, alpha, elim.num, reduce.random)
{
  #library(nlme)
  #return(NULL)
  #mframe<-model.frame(mf.final, data=data, na.action=na.pass)
  #if(length(which(names(mframe) %in% names(data)==FALSE))!=0)
  # {
  #   data$response<-mframe[,1]
  #   fm<-paste(mf.final)
  #   fm[2]<-"response"
  #   mf.final<- as.formula(paste(fm[2],fm[1],fm[3], sep=""))
  #   mf.final<-update.formula(mf.final,mf.final)       
  # }
   
  #model.red<-gls(model = mf.final, data=data, method = "REML", na.action=na.omit)
  model.red <- lm(model, data=model$data)
  
  l.fix <- -2*logLik(model)[1]
  l.red <- -2*logLik(model.red, REML=TRUE)[1]
  
  p.chisq <- 1 - pchisq (l.red -l.fix ,1)
  infoForTerm <- saveInfoForTerm(name.term, l.red -l.fix, 1, p.chisq)
  #detach(package:nlme)
  
  if(infoForTerm$pv > alpha)
  {    
    rand.table<-updateRandTable(infoForTerm, rand.table, elim.num=elim.num, reduce.random=reduce.random)
    model.last<-model.red
  }
  else
  {
    rand.table<-updateRandTable(infoForTerm, rand.table, reduce.random=reduce.random)
    model.last<-model    
  }
  return(list(model=model.last, TAB.rand=rand.table))   
}


### check if the correlation between intercept and slopes is present
isCorrInt<-function(term)
{  
  if(substring.location(term,"+ 1 +")$last !=0)
    return(TRUE)
  if(substring.location(term,"(1 +")$last !=0)
    return(TRUE)
  if(substring.location(term,"+ 1 |")$last !=0)
    return(TRUE)
  return(FALSE)  
}

### check if the correlation between slopes is present
isCorrSlope<-function(term)
{  
  if(substring.location(term,"+ 0 +")$last !=0)
    return(TRUE)
  if(substring.location(term,"(0 +")$last !=0)
    return(TRUE)
  if(substring.location(term,"+ 0 |")$last !=0)
    return(TRUE)
  return(FALSE)  
}

# modify (reduce) the random part when there are slopes 
changeSlopePart<-function(term, isCorr)
{
  sub.loc.div<-substring.location(term," |")
  slopepart<-substring2(term,2,sub.loc.div$first)
  grouppart<-substring2(term,sub.loc.div$last, nchar(term))
  parts<-unlist(strsplit(slopepart, "\\+"))
  
  if(isCorr)
  {
    ind.int<-if(length(which(parts==" 1 "))!=0) which(parts==" 1 ") else which(parts=="1 ")    
    new.terms<-c(paste("(",paste(c(parts[-ind.int], " 0 "), collapse="+"),grouppart, sep=""),paste("(1 ",grouppart,sep=""))
  }
  else
  {
    new.terms<-NULL
    ind.int<-if(length(which(parts==" 0 "))!=0) which(parts==" 0 ") else which(parts=="0 ")
    for(part in parts[-ind.int])
      new.terms<-c(new.terms,paste("(",paste(c(part, " 0 "), collapse="+"),grouppart, sep=""))
    new.terms<-c(new.terms,paste("(1 ",grouppart,sep=""))
  }
  return(new.terms)      
}

### get the random terms
getRandTerms<-function(fmodel)
{
  terms.fm<-attr(terms(fmodel),"term.labels")
  ind.rand.terms<-which(unlist(lapply(terms.fm,function(x) substring.location(x, "|")$first))!=0)
  return(unlist(lapply(terms.fm[ind.rand.terms],function(x) paste("(",x,")",sep=""))))
}

### get names of variables of the slope and group part of random term
getGrSlrand<-function(rand.term)
{
  # find the names of variables for slope part sl and for group part gr
  rand.term1<-substring2(rand.term,2,nchar(rand.term)-1)
  splGrSl<-unlist(strsplit(rand.term1, "|", fixed = TRUE))
  gr<-substring2(splGrSl[2],2,nchar(splGrSl[2]))
  sl<-unlist(strsplit(splGrSl[1], "+", fixed = TRUE))
  for(i in 1:length(sl))
  {
    if(i==1)
      sl[i]<-substring2(sl[i],1,nchar(sl[i])-1)
    else
    {
      sl[i]<-substring2(sl[i],2,nchar(sl[i])-1)
    }
  }
  # change 1 to (Intercept) or eliminate 0 in slope part
  sl[which(sl=="1")]<-"(Intercept)"
  if(length(which(sl=="0"))!=0)
    sl<-sl[-which(sl=="0")]
  return(list(gr=gr,sl=sl))
}

### Find Index of the random term in VarCorr matrix
findIndTerm<-function(vcr,GrSl)
{
  indsGr<-which(names(vcr)==GrSl$gr)
  for(indGr in indsGr)
  {
    if(length(which((names(attr(vcr[[indGr]],"stddev"))==GrSl$sl)==FALSE))==0)
      return(indGr)
  }
}

#### check if the term has zero variance or correlation equal to +-1
checkIsZeroVarOrCorr<-function(model, rand.term, isCorr)
{
 
  vcr<-VarCorr(model)
  ind.term<-findIndTerm(vcr,getGrSlrand(rand.term))
  if(isCorr)
  {
    matcorr<-attr(vcr[[ind.term]],"correlation")
    if(length(which(abs(abs(matcorr[lower.tri(matcorr)])-1)<1e-07 || matcorr[lower.tri(matcorr)]=="NaN"))!=0)
      return(TRUE)      
  }
  else
  {
    stddev<-attr(vcr[[ind.term]],"stddev")
    if(abs(stddev)<1e-07)
      return(TRUE)
  }
  return(FALSE)
}

#### eliminate components with zero variance or correlation +-1, NaN
elimZeroVarOrCorr<-function(model, data, l)
{
  stop=FALSE
  while(!stop)
  {
    fmodel<-formula(model)    
    rand.terms<-getRandTerms(fmodel)
    for(rand.term in rand.terms)
    {
      
      isCorr.int<-isCorrInt(rand.term)
      isCorr.slope<-isCorrSlope(rand.term)
      if(checkIsZeroVarOrCorr(model, rand.term, isCorr.int || (isCorr.slope && length(substring.location(rand.term,"+")$first)>1)))
      {
              
        fm<-paste(fmodel)
        if(isCorr.int || (isCorr.slope && length(substring.location(rand.term,"+")$first)>1)) 
        {
          new.terms<-changeSlopePart(rand.term,isCorr.int)
          fm[3]<-paste(fm[3], "-", rand.term, "+" , paste(new.terms,collapse="+"))
          #print(paste("Random term",rand.term, "was eliminated because of having correlation +-1 or NaN", sep=" "))
          warning(paste(" \n Random term", rand.term, "was eliminated because of having correlation +-1 or NaN \n", sep=" "), call. = FALSE, immediate. = TRUE)
        }
        else
        {
          fm[3]<-paste(fm[3], "-", rand.term)
          warning(paste("\n Random term",rand.term, "was eliminated because of standard deviation being equal to 0 \n", sep=" "), call. = FALSE, immediate. = TRUE)
        }
          
        mf.final<- as.formula(paste(fm[2],fm[1],fm[3], sep=""))
        mf.final<-update.formula(mf.final,mf.final)
        is.present.rand<-checkPresRandTerms(mf.final)
        if(!is.present.rand)
        {
          #library(nlme)
          #model=gls(model = mf.final, data=data, method = "REML", na.action=na.omit)
          #detach(package:nlme)
          model <- lm(model, data=model$data)
          return(list(model=model, TAB.rand=NULL))
        }
        #update model
        #if(!is.null(l))
        #  model<-eval(substitute(lmer(mf.final, data=data, REML=model@dims[["REML"]], contrasts=l),list(mf.final=mf.final)))
        #else
        #  model<-eval(substitute(lmer(mf.final, data=data, REML=model@dims[["REML"]]),list(mf.final=mf.final)))
        model<-updateModel(model, mf.final, model@dims[["REML"]], l)
        elimZero<-TRUE
        break       
      }
      else
        elimZero<-FALSE      
    }
    if(!elimZero)
      return(list(model=model, TAB.rand=NULL))
  }
}

### eliminate NS random terms 
elimRandEffs<-function(model, data, alpha, reduce.random, l)
{
  isInitRand<-TRUE
  elim.num<-1
  stop=FALSE
  while(!stop)
  {
    fmodel<-formula(model)    
    rand.terms<-getRandTerms(fmodel)
    
    if(isInitRand)
    {
      rand.table<-initRandTable(rand.terms, reduce.random)
      isInitRand<-FALSE
    }      
    fm<-paste(fmodel)
    pv.max<-0
    infoForTerms<-vector("list", length(rand.terms))
    names(infoForTerms)<-  rand.terms
    
    for(rand.term in rand.terms)
    {
      fm<-paste(fmodel)
      isCorr.int<-isCorrInt(rand.term)
      isCorr.slope<-isCorrSlope(rand.term)
      if(isCorr.int || (isCorr.slope && length(substring.location(rand.term,"+")$first)>1)) 
      {
        new.terms<-changeSlopePart(rand.term,isCorr.int)
        fm[3]<-paste(fm[3], "-", rand.term, "+" , paste(new.terms,collapse="+"))
      }
      else
        fm[3]<-paste(fm[3], "-", rand.term)
      mf.final<- as.formula(paste(fm[2],fm[1],fm[3], sep=""))
      mf.final<-update.formula(mf.final,mf.final)
      is.present.rand<-checkPresRandTerms(mf.final)
      
      # no more random terms in the model
      if(!is.present.rand)
      {
        return(compareMixVSFix(model, mf.final, data, rand.term, rand.table, alpha, elim.num, reduce.random))
        
      } 
      #if(!is.null(l))
      #  model.red<-eval(substitute(lmer(mf.final, data=data, contrasts=l),list(mf.final=mf.final)))
      #else
      #  model.red<-eval(substitute(lmer(mf.final, data=data),list(mf.final=mf.final)))
      model.red<-updateModel(model, mf.final, model@dims[["REML"]], l)
      anova.red<-anova(model, model.red)
      infoForTerms[[rand.term]]<-saveInfoForTerm(rand.term, anova.red$Chisq[2], anova.red[2,6] , anova.red$Pr[2])
            
      if((anova.red$Pr[2] >= pv.max) && reduce.random)
      { 
           pv.max<-anova.red$Pr[2]
           infoForTermElim<-infoForTerms[[rand.term]]
           model.final<-model.red 
           if(anova.red$Pr[2]==1)
             break
      }
    }    
    
    if(!reduce.random)
    {
      rand.table<-updateRandTable(infoForTerms, rand.table, reduce.random=reduce.random)
      model.last<-model
      break
    }
    
    rand.terms.upd<-getRandTerms(formula(model.final))
    
    if(infoForTermElim$pv > alpha)
    {
      rand.table<-updateRandTable(infoForTermElim, rand.table, rand.terms.upd[!rand.terms.upd %in% rand.terms] , elim.num, reduce.random)
      elim.num=elim.num+1      
    }
    else
    {
      rand.table<-updateRandTable(infoForTerms, rand.table, reduce.random=reduce.random)
      model.last<-model
      break
    }
    
    model<-model.final  
  }  
  return(list(model=model.last, TAB.rand=rand.table))
}

#####save results for fixed effects for model with only fixed effects
saveResultsFixModel<-function(result, model)
{
  result$anova.table<-anova(model)
  result$model<-model
  lsmeans.summ<- matrix(ncol=7,nrow=0)
  colnames(lsmeans.summ)<-c("Estimate","Standard Error", "DF", "t-value", "Lower CI", "Upper CI", "p-value")
  result$lsmeans.table<-lsmeans.summ
  result$diffs.lsmeans.table<-lsmeans.summ
  return(result)
}

################# UNUSED function #########################
#code from lme4 package
###########################################################
formatVC <- function(varc, digits = max(3, getOption("digits") - 2))
### "format()" the 'VarCorr' matrix of the random effects -- for show()ing
{
    sc <- unname(attr(varc, "sc"))
    recorr <- lapply(varc, attr, "correlation")
    reStdDev <- c(lapply(varc, attr, "stddev"), list(Residual = sc))
    reLens <- unlist(c(lapply(reStdDev, length)))
    nr <- sum(reLens)
    reMat <- array('', c(nr, 4),
       list(rep.int('', nr),
			c("Groups", "Name", "Variance", "Std.Dev.")))
    reMat[1+cumsum(reLens)-reLens, 1] <- names(reLens)
    reMat[,2] <- c(unlist(lapply(varc, colnames)), "")
    reMat[,3] <- format(unlist(reStdDev)^2, digits = digits)
    reMat[,4] <- format(unlist(reStdDev), digits = digits)
    if (any(reLens > 1)) {
	maxlen <- max(reLens)
	corr <-
	    do.call("rBind",
		    lapply(recorr,
			   function(x, maxlen) {
			       x <- as(x, "matrix")
			       cc <- format(round(x, 3), nsmall = 3)
			       cc[!lower.tri(cc)] <- ""
			       nr <- dim(cc)[1]
			       if (nr >= maxlen) return(cc)
			       cbind(cc, matrix("", nr, maxlen-nr))
			   }, maxlen))
	colnames(corr) <- c("Corr", rep.int("", maxlen - 1))
	cbind(reMat, rBind(corr, rep.int("", ncol(corr))))
    } else reMat
}


#update model
updateModel<-function(model, mf.final, reml, l)
{
  return(suppressWarnings(update(object=model, formula.=mf.final, REML=reml, contrasts=l)))
  
  
#   if(!is.null(l)) 
#   {
#     #l<-l[names(l) %in% strsplit(paste(mf.final[3]), split=" + ", fixed=TRUE)[[1]]]
#     #print(mf.final)
#     #print(l)
#     model<-suppressWarnings(update(model, mf.final, REML=model@dims[["REML"]], contrasts=l))
#   }
#   else
#     model<-suppressWarnings(update(model, mf.final, REML=model@dims[["REML"]]))
#   return(model)
}


