totalAnovaRandLsmeans <- function(model, ddf="Satterthwaite", type = 3, alpha.random = 0.1, alpha.fixed = 0.05, reduce.fixed = TRUE, reduce.random = TRUE, lsmeans.calc = TRUE, difflsmeans.calc=TRUE,  isTotal=FALSE, isAnova=FALSE, isRand=FALSE, isLSMEANS=FALSE, isDiffLSMEANS=FALSE, isTtest=FALSE, test.effs=NULL, method.grad="simple")
{
  #errors in specifying the parameters
  if(!isRand && !(ddf %in% c("Satterthwaite","Kenward-Roger")))
  {
    print ('Error: parameter ddf is wrongly specified')
    stop()
  }
  if(!isRand && !(type %in% c(1,3)))
  {
    print ('Error: parameter type is wrongly specified')
    stop()
  } 
  
  data <- model.frame(model) #summary(model,"lme4")@frame 
  
  
  #update contrasts
  mm <- model.matrix(model)
  l <- attr(mm,"contrasts")
  contr <- l
  ### change contrasts for F tests calculations
  #list of contrasts for factors
  if( isAnova || isTotal )
  {    
    if( length(which(unlist(contr)!="contr.SAS")) > 0 )
    {
      names.facs <- names(contr)
      l <- as.list(rep("contr.SAS",length(names.facs)))
      names(l) <- names(contr)
      #warning(" \nmodel has been refitted with contrasts=contr.SAS \n")
      #model<-update(model,.~., data=data, contrasts=l)
	    model <- updateModel(model, .~., getREML(model), l) 
    }    
  }
  else
  {
    #update model to mer class
	model <- updateModel(model, .~., getREML(model), l)
  }
  
  
  
  
  #not to show the warnings  
  #options(warn=-1) 
  
  result <- NULL
  anova.table <- NULL
  
  result$response <- rownames(attr(terms(model),"factors"))[1]#names(attr(terms(model),"dataClasses")[1])
  
  #model<-lmer(formula=formula, data=data)
  
  #update model
  # change unordered contrasts to contr.SAS
  # change REML to TRUE
  #options(contrasts=c(unordered="contr.SAS", ordered="contr.poly"))
  
  #model<-update(model, REML=TRUE)
  ## deleted because use ML for anova(m1, m2) for random effects
#   if( isRand || isTotal || (ddf=="Kenward-Roger" && (isTotal || isAnova)) )
#   {
#     
#       model<-
#       if (getREML(model) == 1)
#       {
#         model
#       }
#       else
#       {
#         warning("\n model has been refitted with REML=TRUE \n")
#         updateModel(model, .~., reml=TRUE, l)
#       }
#   }
#   
  mf.final <- update.formula(formula(model),formula(model)) 
  
  
  
  data <- data[complete.cases(data),]
  
  
    
  model <- updateModel(model, mf.final, getREML(model), l)
  
   
  
  # save the call of the model              
  result$call <- model@call
  
  #check if there are correlations between intercept and slope
  result$corr.intsl <- checkCorr(model)  
  
  
 
  if(isRand || isTotal)
  {
    result.rand <- elimZeroVarOrCorr(model, data, l)
    model <- result.rand$model      
  }    
  
  
  #save results for fixed effects for model with only fixed effects
  if(class(model) == "lm" | class(model) == "gls")
  {
    result <- saveResultsFixModel(result, model)
    result$rand.table=NULL
    return(result)
  }
  
  
  
 
  #analysis of the random part  
  if(isRand || isTotal)
  {
    if(isRand)
      reduce.random <- FALSE
    result.rand <- elimRandEffs(model, data, alpha.random, reduce.random, l)  
   
    model <- result.rand$model
    #convert rand table to data frame
    rt <- as.data.frame(result.rand$TAB.rand)
    rt$Chi.DF <- as.integer(rt$Chi.DF)
    if(!is.null(rt$elim.num))
      rt$elim.num <- as.integer(rt$elim.num)
    
    result$rand.table <- rt
    if(isRand)
      return(result)
  }
  
  
  
  
  
  #save results for fixed effects for model with only fixed effects
  if(class(model) == "lm" | class(model) == "gls")
    return(saveResultsFixModel(result, model))
  
  
  #perform reduction of fixed effects for model with mixed effects
  stop = FALSE
  is.first.anova <- TRUE
  is.first.sign <- TRUE  
  
  
  
  while(!stop)
  {      
    
    # if there are no fixed terms
    if(nrow(anova(model, ddf="lme4"))==0)
    {
      if(is.null(anova.table))
      {
        if(isLSMEANS || isDiffLSMEANS)
        {
          lsmeans.summ <-  matrix(ncol=7,nrow=0)
          colnames(lsmeans.summ) <- c("Estimate","Standard Error", "DF", "t-value", "Lower CI", "Upper CI", "p-value")
          lsmeans.summ <- as.data.frame(lsmeans.summ)
          if(isLSMEANS)
            result$lsmeans.table <- lsmeans.summ
          if(isDiffLSMEANS)
            result$diffs.lsmeans.table <- lsmeans.summ
          return(result)
        }
        if(isTtest)
        {
          # save lmer outcome in rho environmental variable
          rho <- rhoInit(model)     
          
          # calculate asymptotic covariance matrix A
          h  <-  hessian(function(x) Dev(rho,x), rho$param$vec.matr)
          
          rho$A <- 2*solve(h)
          #rho$A <- 2*ginv(h)
          
          tsummary <- calculateTtest(rho, diag(rep(1,length(rho$fixEffs))),length(rho$fixEffs), method.grad)
          result$ttest <- list(df=tsummary[,"df"], tvalue=tsummary[,"t value"], tpvalue=tsummary[,"p-value"])
        }
        result$model <- model
        result$anova.table <- anova(model, ddf="lme4")
        return(result)        
        
      }          
      break
    }        
    
    
    # save lmer outcome in rho environmental variable
    rho <- rhoInit(model)     
    
    # calculate asymptotic covariance matrix A
    h  <-  hessian(function(x) Dev(rho,x), rho$param$vec.matr)
    #h  <-  myhess(function(x) Dev(rho,x), rho$param$vec.matr)
    ch <- try(chol(h), silent=TRUE)
    if(inherits(ch, "try-error")) {
      message("Model is not identifiable...")
    }
    rho$A <- 2*chol2inv(ch)
    
    eigval <- eigen(h, symmetric=TRUE, only.values=TRUE)$values
    isposA <- TRUE
    if(min(eigval) < sqrt(.Machine$double.eps)) ## tol ~ sqrt(.Machine$double.eps)
      isposA <- FALSE
    
    #rho$A  <-  2*solve(h)

    
    #Check if A is positive-definite
    #isposA <- all(eigen(rho$A)$values>0)      
    if(!isposA)
    {
      print("Asymptotic covariance matrix A is not positive!")
      #previous code: return ERROR
      #print("ERROR: asymptotic covariance matrix A is not positive!")
      #result$model <- model
      #TABs <- emptyAnovaLsmeansTAB()
      #result$anova.table <- TABs$TAB.fixed
      #result$lsmeans.table <- TABs$TAB.lsmeans
      #result$diffs.lsmeans.table <- TABs$TAB.lsmeans
      #return(result)
    }
    
    
    #calculate ttest and p-values for summary
    if(isTtest)
    {
      tsummary <- calculateTtest(rho, diag(rep(1,length(rho$fixEffs))), 
                                 length(rho$fixEffs), method.grad)
      result$ttest <- list(df=tsummary[, "df"], tvalue=tsummary[, "t value"], 
                           tpvalue=tsummary[, "p-value"])
      return(result)
    }
    
    
    #calculate lsmeans of differences of LSMEANS of the final model
    if(isLSMEANS || isDiffLSMEANS)
    {
      if(isLSMEANS)
      {
        lsmeans.tab <- calcLSMEANS(model, data, rho, alpha.fixed, test.effs=test.effs, method.grad=method.grad, lsmeansORdiff=TRUE, l)
        result$lsmeans.table <- lsmeans.tab$summ.data
        result$diffs.lsmeans.table <- NULL
      }
      if(isDiffLSMEANS)
      {
        lsmeans.tab <- calcLSMEANS(model, data, rho, alpha.fixed, test.effs=test.effs, method.grad=method.grad, lsmeansORdiff=FALSE, l)
        result$diffs.lsmeans.table <- lsmeans.tab$summ.data
        result$lsmeans.table <- NULL
      }
      return(result)
    }
    
    
    # Calculate  F-test with Satterthwaite's approximation
    # create X design matrix for fixed effects
    X.design.list <- createDesignMat(model,data)
    X.design <- X.design.list$X.design
    names.design.withLevels <- X.design.list$names.design.withLevels
    
    
    #save full coefficients in rho
    ## old code, worked with dummy.coef
    #nums.dummy.coefs <- getNumsDummyCoefs(model, data, l)
    #rho$nums.zeroCoefs <- nums.dummy.coefs$nums.zeroCoefs
    #rho$nums.Coefs <- nums.dummy.coefs$nums.Coefs
    #fullCoefs <- rep(0, ncol(X.design))
    #fullCoefs[rho$nums.Coefs] <- rho$fixEffs
    #fullCoefs <- setNames(fullCoefs, names.design.withLevels) 
    ###new code with X.design matrix
    fullCoefs <- rep(0, ncol(X.design))
    fullCoefs <- setNames(fullCoefs, names.design.withLevels) 
    names(fullCoefs)[1] <- "(Intercept)"
    fullCoefs[names(rho$fixEffs)] <- rho$fixEffs
    rho$nums.Coefs <- which(names(fullCoefs) %in% names(rho$fixEffs))
      #unique(c(which(names(fullCoefs)=="(Intercept)"),which(names(fullCoefs) %in% names(rho$fixEffs))))
    rho$nums.Coefs <- setNames(rho$nums.Coefs, names(fullCoefs[rho$nums.Coefs]))
    #define the terms that are to be tested
    test.terms <- attr(terms(model),"term.labels")
    
    #initialize anova table
    if(is.first.anova)
    {
      anova.table <- initAnovaTable(model, test.terms, reduce.fixed)
      is.first.anova <- FALSE
      elim.num <- 1
    }
    
    # calculate general set matrix for type 3 hypothesis
    if( type==3 )
      L <- calcGeneralSetForHypothesis(X.design, rho)  
    
    
    # calculate type 1 hypothesis matrices for each term
    # TODO: fix bug with noint sens1 + Homesize
    if( type==1 )
    {
      #X <- model.matrix(model)
      X <- X.design
      p <- ncol(X)
      XtX <- crossprod(X)
      U <- doolittle(XtX)$U
      d <- diag(U)
      for(i in 1:nrow(U))
        if(d[i] > 0) U[i, ] <- U[i, ] / d[i]
      L <- U

    }
 
    #calculate ss F value ddf and p value for each term 
    resultFpvalueSS <- lapply(test.terms, calcFpvalueMAIN, L=L, X.design=X.design, fullCoefs=fullCoefs, model=model, rho=rho, ddf=ddf, method.grad=method.grad, type=type)  
    #fill anova table
    anova.table <- fillAnovaTable(resultFpvalueSS,  anova.table)
    
   
    if(!reduce.fixed)             
      break     
    else
    {
      resNSelim <- elimNSFixedTerm(model, anova.table, data, alpha.fixed, elim.num, l)
      if(is.null(resNSelim))
        break
      else
      {
        model <- resNSelim$model
        mf.final <- update.formula(formula(model),formula(model))
        model <- updateModel(model, mf.final, getREML(model), l)        
        anova.table <- updateAnovaTable(resNSelim)
        elim.num <- elim.num+1
       
      }        
    }      
    
  }
  
  #convert anova table to data frame
  anova.table <- as.data.frame(anova.table)
  anova.table$NumDF <- as.integer(anova.table$NumDF)
  if(!is.null(anova.table$elim.num))
    anova.table$elim.num <- as.integer(anova.table$elim.num)
  
  if(isTotal || isAnova)
  {
    result$anova.table <- anova.table
    if(isAnova)
      return(result)
  }
  
  #if in step function least squares means of diffs of LSMEANS are required
  if(lsmeans.calc)
  {
    lsmeans.tab <- calcLSMEANS(model, data, rho, alpha.fixed, test.effs=test.effs, method.grad=method.grad, lsmeansORdiff=TRUE, l)
    result$lsmeans.table <- lsmeans.tab$summ.data
  }
  else
  {
    result$lsmeans.table <- NULL
  }
  if(difflsmeans.calc)
  {
    lsmeans.tab <- calcLSMEANS(model, data, rho, alpha.fixed, test.effs=test.effs, method.grad=method.grad, lsmeansORdiff=FALSE, l)
    result$diffs.lsmeans.table <- lsmeans.tab$summ.data
  }
  else
  {
    result$diffs.lsmeans.table <- NULL
  }
  
  tsummary <- calculateTtest(rho, diag(rep(1,length(rho$fixEffs))), length(rho$fixEffs), method.grad)
  result$ttest <- list(df=tsummary[,"df"], tvalue=tsummary[,"t value"], tpvalue=tsummary[,"p-value"])
  
  #format anova.table and random.table according to elim.num column
  result$anova.table <- formatElimNumTable(result$anova.table) 
  result$rand.table <- formatElimNumTable(result$rand.table) 
  
  #update final model
  mf.final <- update.formula(formula(model),formula(model))
  model <- updateModel(model, mf.final, getREML(model), contr)
  
  #save model
  if(inherits(model, "merMod"))
    model <- as(model,"merModLmerTest")
  #else if (inherits(model, "mer"))    
  #  model <- as(model,"merLmerTest")
  #if(class(model)=="merLmerTest")
  #  model@t.pval <- result$ttest$tpvalue
  result$model <- model
  return(result)
}



step <- function(model, ddf="Satterthwaite", type=3, alpha.random = 0.1, alpha.fixed = 0.05, reduce.fixed = TRUE, reduce.random = TRUE, lsmeans.calc=TRUE, difflsmeans.calc=TRUE, test.effs=NULL, method.grad="simple", ...)
{  
  result <- totalAnovaRandLsmeans(model=model, ddf=ddf , type=type,  alpha.random=alpha.random, alpha.fixed=alpha.fixed, reduce.fixed=reduce.fixed, reduce.random=reduce.random, lsmeans.calc=lsmeans.calc, difflsmeans.calc=difflsmeans.calc, isTotal=TRUE, isTtest=FALSE, test.effs=test.effs, method.grad=method.grad)
  class(result) <- "step"
  result
}

### UNUSED function
#totalAnalysis.formula <- function(formula, data, ...)
#{
#  model <- lmer(formula=formula, data=data)
#  resAnalysis<-totalAnalysis.default(model, data, ...)
#  resAnalysis$call<-match.call()
#  resAnalysis
#}

print.step <- function(x, ...)
{
  
  if(!is.null(x$rand.table))
  {
    cat("\nRandom effects:\n") 
    x$rand.table[,"p.value"] <- format.pval(x$rand.table[,"p.value"], digits=4, eps=1e-7)
    x$rand.table[,"Chi.sq"] <- round(x$rand.table[,"Chi.sq"],2)
    print(x$rand.table)
    #printCoefmat(x$rand.table, digits=3 , dig.tst=1  ,tst.ind=which(colnames(x$rand.table)=="Chi.DF"), P.values=TRUE, has.Pvalue=TRUE, na.print = "KEEP")
  }
  
  if(nrow(x$anova.table) != 0)
  {
    if(class(x$model) == "lm" | class(x$model) == "gls")
    {
      cat("\nFixed effects:\n")
      print(x$anova.table)
      cat("\nLeast squares means:\n")
      print(x$lsmeans.table) 
      cat("\nFinal model:\n")
      print(x$model)
      return()
    }
    else
    {
      cat("\nFixed effects:\n")
      x$anova.table[,"Pr(>F)"] <- format.pval(x$anova.table[,"Pr(>F)"], digits=4, eps=1e-7)
      x$anova.table[,c("Sum Sq","Mean Sq", "F.value")] <- round(x$anova.table[,c("Sum Sq","Mean Sq", "F.value")],4)
      x$anova.table[,"DenDF"] <- round(x$anova.table[,"DenDF"],2)
      print(x$anova.table) 
       #printCoefmat(x$anova.table, dig.tst=3, tst.ind=3, cs.ind=3, digits=3 ,P.values = TRUE, has.Pvalue=TRUE)
      if(!is.null(x$lsmeans.table))
      {
        cat("\nLeast squares means:\n")
        printCoefmat(x$lsmeans.table, dig.tst=3 ,tst.ind=c(1:(which(colnames(x$lsmeans.table)=="Estimate")-1),which(colnames(x$lsmeans.table)=="DF")), digits=3 ,P.values = TRUE, has.Pvalue=TRUE)
      }
      if(!is.null(x$diffs.lsmeans.table))
      {
        cat("\n Differences of LSMEANS:\n")
        printCoefmat(x$diffs.lsmeans.table, dig.tst=1  ,tst.ind=c(1:(which(colnames(x$diffs.lsmeans.table)=="Estimate")-1),which(colnames(x$diffs.lsmeans.table)=="DF")), digits=3 ,P.values = TRUE, has.Pvalue=TRUE)
      }
      
    }    
  }
  else
    print(x$anova.table)
  cat("\nFinal model:\n")
  print(x$model@call) 
}


plot.step <- function(x, ...)
{
  if(!is.null(x$lsmeans.table) && nrow(x$lsmeans.table)>0)
    plotLSMEANS(x$lsmeans.table, x$response, "LSMEANS")     
  if(!is.null(x$diffs.lsmeans.table) && nrow(x$diffs.lsmeans.table)>0)
    plotLSMEANS(x$diffs.lsmeans.table, x$response, "DIFF of LSMEANS")
}

# lmer <-
#   function(formula, data, family = NULL, REML = TRUE,
#            control = list(), start = NULL, verbose = FALSE, doFit = TRUE,
#            subset, weights, na.action, offset, contrasts = NULL,
#            model = TRUE, x = TRUE, ...)
 lmer <- function(formula, data = NULL, REML = TRUE,
          control = lmerControl(), start = NULL, verbose = 0L,
          subset, weights, na.action, offset, contrasts = NULL,
          devFunOnly = FALSE, ...)
  {
    mc <- match.call()
    mc[[1]] <- quote(lme4::lmer)
    model <- eval.parent(mc)
    if(inherits(model, "merMod"))
      model <- as(model,"merModLmerTest")
    #else if (inherits(model, "mer"))    
    #  model <- as(model,"merLmerTest")   
    #model <- as(model,"merLmerTest")
#     if(class(model) == "lmerMod")
#      model <- as(model,"merModLmerTest")
#     else if (class(model) == "mer")
#      model <- as(model,"merLmerTest")
    #tryCatch(  { result = glm( y~x , family = binomial( link = "logit" ) ) } , error = function(e) { print("test") } )
    #t.pval <- tryCatch( {totalAnovaRandLsmeans(model=model, ddf="Satterthwaite", isTtest=TRUE)$ttest$tpvalue}, error = function(e) { NULL })
    #if(!is.null(t.pval))
    #{
    #  model@t.pval  <- t.pval
    #}
    #else
    #{
    #  model <- as(model,"lmerMod")
    #}
    return(model)
  }




setMethod("anova", signature(object="merModLmerTest"),
          function(object,..., ddf="Satterthwaite", type=3, method.grad="simple")  
          {
            
            mCall <- match.call(expand.dots = TRUE)
            dots <- list(...)
            modp <- if (length(dots))
              sapply(dots, is, "merModLmerTest") | sapply(dots, is, "merMod") | sapply(dots, is, "lm") else logical(0)
            if (any(modp)) {
              return(callNextMethod())
            }
            else
            {
              cnm <- callNextMethod()
              if(!is.null(ddf)&& ddf=="lme4") 
                return(cnm)              
             {
                  table <- cnm
                  
                  
                  an.table <- tryCatch({totalAnovaRandLsmeans(model=object, ddf=ddf, type=type, isAnova=TRUE, reduce.random=FALSE, reduce.fixed=FALSE, method.grad=method.grad)$anova.table}, error = function(e) { NULL })
                  if(!is.null(an.table))
                  {
                    table <- an.table
#                     rnames <- rownames(table)
#                     if(nrow(an.table)>0)
#                     {
#                       
#                       table <- as.data.frame(cbind(table$Df, table$"Sum Sq", table$"Mean Sq", an.table[,"F.value"], an.table[,"DenDF"], an.table[,"Pr(>F)"]))
#                       colnames(table) <- c("Df", "Sum Sq", "Mean Sq", "F value", "Denom", "Pr(>F)")
#                       dimnames(table) <- list(rnames,
#                                               c("Df", "Sum Sq", "Mean Sq", "F value", "Denom", "Pr(>F)"))
#                     }
#                     else                    
#                       table <- an.table
                    
                    attr(table, "heading") <- paste("Analysis of Variance Table of type", type ," with ", ddf, "\napproximation for degrees of freedom")
                  }
                  
                  
                  
                  class(table) <- c("anova", "data.frame")
                  return(table)
                }  
              
            }
            
          })

setMethod("summary", signature(object = "merModLmerTest"),
          function(object, ddf="Satterthwaite", ...)
          {
            
            cl <- callNextMethod()
            if(!is.null(ddf) && ddf=="lme4") return(cl)
            else
            {
              tsum <- tryCatch( {totalAnovaRandLsmeans(model=object, ddf="Satterthwaite", isTtest=TRUE)$ttest}, error = function(e) { NULL })
              coefs.satt <- cbind(cl$coefficients[,1:2, drop=FALSE], tsum$df, tsum$tvalue, tsum$tpvalue)
               # t.pval <- tryCatch( {totalAnovaRandLsmeans(model=object, ddf="Satterthwaite", isTtest=TRUE)$ttest$tpvalue}, error = function(e) { NULL })
               # coefs.satt <- cbind(cl$coefficients, t.pval) 
                cl$coefficients <- coefs.satt
                colnames(cl$coefficients)[3:5] <- c("df","t value","Pr(>|t|)")              
            }   
            
            return(cl)
          }
          
)

#randTAB.default<-function(model, data, ...)
rand <- function(model, ...)
{
  result <- totalAnovaRandLsmeans(model=model, isRand=TRUE, reduce.random=FALSE)  
  res <- list(rand.table=result$rand.table, isCorr = result$corr.intsl)
  class(res) <- "rand"
  res
}

print.rand <- function(x, ...)
{
  
  cat("Analysis of Random effects Table:\n")
  if(!is.null(x$rand.table))
    printCoefmat(x$rand.table, digits=3 , dig.tst=1  ,tst.ind=c(which(colnames(x$rand.table)=="Chi.DF"),which(colnames(x$rand.table)=="elim.num")), P.values=TRUE, has.Pvalue=TRUE)        
}





lsmeans <- function(model, test.effs=NULL, method.grad="simple", ...)
{
  result <- totalAnovaRandLsmeans(model=model, ddf="Satterthwaite", isLSMEANS=TRUE, test.effs=test.effs, reduce.random=FALSE, reduce.fixed=FALSE, method.grad=method.grad)  
  res <- list(lsmeans.table=result$lsmeans.table, response=result$response)
  class(res) <- "lsmeans"
  res 
}

print.lsmeans <- function(x, ...)
{
  
  cat("Least Squares Means table:\n")
  printCoefmat(data.matrix(x$lsmeans.table), dig.tst=1, tst.ind=c(1:(which(colnames(x$lsmeans.table)=="Estimate")-1),which(colnames(x$lsmeans.table)=="DF")), digits=3 , P.values=TRUE, has.Pvalue=TRUE)       
}

plot.lsmeans <- function(x, ...)
{
  
  #plots for LSMEANS
  if(!is.null(x$lsmeans.table) && nrow(x$lsmeans.table)>0)
    plotLSMEANS(x$lsmeans.table, x$response, "LSMEANS")     
}

difflsmeans <- function(model, test.effs=NULL, method.grad="simple", ...)
{
  result <- totalAnovaRandLsmeans(model=model, ddf="Satterthwaite", isDiffLSMEANS=TRUE, test.effs=test.effs, reduce.random=FALSE, reduce.fixed=FALSE, method.grad=method.grad)  
  res <- list(diffs.lsmeans.table=result$diffs.lsmeans.table, response=result$response)
  class(res) <- "difflsmeans"
  res 
}

print.difflsmeans <- function(x, ...)
{
  
  cat("Differences of LSMEANS:\n")
  printCoefmat(data.matrix(x$diffs.lsmeans.table), dig.tst=1, tst.ind=c(1:(which(colnames(x$diffs.lsmeans.table)=="Estimate")-1),which(colnames(x$diffs.lsmeans.table)=="DF")), digits=3 , P.values=TRUE, has.Pvalue=TRUE)
  
}

plot.difflsmeans <- function(x, ...)
{
  
  #plots for DIFF of LSMEANS
  if(!is.null(x$diffs.lsmeans.table) && nrow(x$diffs.lsmeans.table)>0)
    plotLSMEANS(x$diffs.lsmeans.table, x$response, "DIFF of LSMEANS")   
}
