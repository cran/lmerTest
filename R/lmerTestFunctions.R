totalAnovaRandLsmeans <- function(model, ddf="Satterthwaite", type=3, alpha.random = 0.1, alpha.fixed = 0.05, reduce.fixed = TRUE, reduce.random = TRUE, lsmeans.calc = TRUE, difflsmeans.calc=TRUE,  isTotal=FALSE, isAnova=FALSE, isRand=FALSE, isLSMEANS=FALSE, isDiffLSMEANS=FALSE, isTtest=FALSE, test.effs=NULL, method.grad="simple")
{
  #errors in specifying the parameters
  if(!isRand && !(ddf %in% c("Satterthwaite","Kenward-Roger")))
  {
    print ('Error: parameter ddf is wrongly specified')
    stop()
  }
  if(!isRand && type!=3)
  {
    print ('Error: parameter type is wrongly specified')
    stop()
  } 
  
  data <- summary(model,"lme4")@frame
  
   
  
  #update contrasts
  mm <- model.matrix(model)
  l <- attr(mm,"contrasts")
  contr<-l
  ### change contrasts for F tests calculations
  #list of contrasts for factors
  if(isAnova || isTotal)
  {
    
    if(length(which(unlist(contr)!="contr.SAS"))>0)
    {
      names.facs <- names(contr)
      l <- as.list(rep("contr.SAS",length(names.facs)))
      names(l) <- names(contr)
      #warning(" \nmodel has been refitted with contrasts=contr.SAS \n")
      #model<-update(model,.~., data=data, contrasts=l)
      model <- updateModel(model, .~., model@dims[["REML"]], l)
     
    }
    
  }
  else
  {
    #update model to mer class
    model <- updateModel(model, .~., model@dims[["REML"]], l)
  }
  
  
  
  
  #not to show the warnings  
  #options(warn=-1)    
   
 
  
  result<-NULL
  anova.table<-NULL
  
  result$response<-names(attr(terms(model),"dataClasses")[1])
  
  #model<-lmer(formula=formula, data=data)
  
  #update model
  # change unordered contrasts to contr.SAS
  # change REML to TRUE
  #options(contrasts=c(unordered="contr.SAS", ordered="contr.poly"))
  
  #model<-update(model, REML=TRUE)
  if(isRand || isTotal || (ddf=="Kenward-Roger" && (isTotal || isAnova)))
  model<-
        if (model@dims['REML'] == 1)
        {
            model
        }
        else
        {
  	    warning("\n model has been refitted with REML=TRUE \n")
           updateModel(model, .~., reml=TRUE, l)
        }
  mf.final<-update.formula(formula(model),formula(model))
  
  
  
  #update data 
  #eliminate missing values
  ### old code ############
  #if(ncol(get_all_vars(mf.final, data))>ncol(data))
  #{
  #   #mframe<-model.frame(mf.final, data=data, na.action=na.pass)
  #   #data$response<-mframe[,1]
  #   data$response<-data[,1]
  #   data<-data[,-1]
  #   fm<-paste(mf.final)
  #   fm[2]<-"response"
  #   mf.final<- as.formula(paste(fm[2],fm[1],fm[3], sep=""))
  #   mf.final<-update.formula(mf.final,mf.final)
  #}
  
  #data<-get_all_vars(mf.final, data)
  ###############################################################
  
  data<-data[complete.cases(data),]
  
  
   
  
  #if(!is.null(l))
  #   model<-eval(substitute(lmer(mf.final, data=data, REML=model@dims[["REML"]], contrasts=l),list(mf.final=mf.final)))
  #else
  #   model<-eval(substitute(lmer(mf.final, data=data, REML=model@dims[["REML"]]),list(mf.final=mf.final)))
  #if(!is.null(l))
  #   model<-update(model, mf.final, data=data, REML=model@dims[["REML"]], contrasts=l)
  #else
  #   model<-update(model, mf.final, data=data, REML=model@dims[["REML"]])

  model<-updateModel(model, mf.final, model@dims[["REML"]], l)
 
    # check if there are no fixed effects
  #if(length(attr(delete.response(terms(model)),"term.labels"))==0)
 # {
  #  fm<-getFormula(model, withRand=TRUE)
  #  fm<-as.formula(paste(fm[2],fm[1],paste("1+", fm[3], sep=""), sep=""))
  #  model<-eval(substitute(lmer(mf.final, data=data),list(mf.final=mf.final)))
  #}
   
  
   
  # save the call of the model              
  result$call<-model@call
  
  #check if there are correlations between intercept and slope
  result$corr.intsl<-checkCorr(model)
  
  
  
  
  #Reduction of random effects
  #if(isRand || isTotal )
  #{     
    # perform reduction of random effects
    #if(!isRandReduce)
    #  result.rand<-elimRandLmer(model, data, 1)
    #else
    #  result.rand<-elimRandLmer(model, data, alpha.rand)
    
    if(isRand || isTotal)
    {
      result.rand<-elimZeroVarOrCorr(model, data, l)
      model<-result.rand$model      
    }
    
     
    
    #save results for fixed effects for model with only fixed effects
    if(class(model) == "lm" | class(model) == "gls")
    {
      result<-saveResultsFixModel(result, model)
      result$rand.table=NULL
      return(result)
    }
    
  
   
    #if(!isRandReduce)
    #  result.rand<-elimRandEffs(model, data, 1)
  
    #analysis of the random part  
    if(isRand || isTotal)
    {
      if(isRand)
        reduce.random<-FALSE
      #if(reduce.random)
      result.rand<-elimRandEffs(model, data, alpha.random, reduce.random, l)  
      #if(!reduce.random && isRand)
      #  result.rand<-elimRandEffs(model, data, 1)
  
      model<-result.rand$model
      #convert rand table to data frame
      rt<-as.data.frame(result.rand$TAB.rand)
      rt$Chi.DF<-as.integer(rt$Chi.DF)
      if(!is.null(rt$elim.num))
        rt$elim.num<-as.integer(rt$elim.num)
  
      result$rand.table<-rt
      if(isRand)
        return(result)
    }
    
       
    
    
  #}
  
   
  
  #save results for fixed effects for model with only fixed effects
  if(class(model) == "lm" | class(model) == "gls")
    return(saveResultsFixModel(result, model))
  
    
  #perform reduction of fixed effects for model with mixed effects
  stop = FALSE
  is.first.anova<-TRUE
  is.first.sign<-TRUE  
   
  
  
  while(!stop)
  {
  
      
    
      # if there are no fixed terms
      if(nrow(anova(model, ddf="lme4"))==0)
      {
        if(is.null(anova.table))
        {
          
          if(isLSMEANS || isDiffLSMEANS)
          {
             lsmeans.summ<- matrix(ncol=7,nrow=0)
             colnames(lsmeans.summ)<-c("Estimate","Standard Error", "DF", "t-value", "Lower CI", "Upper CI", "p-value")
             lsmeans.summ<-as.data.frame(lsmeans.summ)
             if(isLSMEANS)
               result$lsmeans.table<-lsmeans.summ
             if(isDiffLSMEANS)
               result$diffs.lsmeans.table<-lsmeans.summ
             return(result)
          }
          if(isTtest)
          {
            # save lmer outcome in rho environmental variable
            rho<-rhoInit(model)     
            
            # calculate asymptotic covariance matrix A
            h <- hessian(function(x) Dev(rho,x), rho$param$vec.matr)
            rho$A <- 2*solve(h)
            #rho$A <- 2*ginv(h)
            
            tsummary<-calculateTtest(rho, diag(rep(1,nrow(rho$s@coefs))), nrow(rho$s@coefs), method.grad)
            result$ttest<-list(df=tsummary[,"df"], tvalue=tsummary[,"t value"], tpvalue=tsummary[,"p-value"])
          }
          result$model<-model
          result$anova.table<-anova(model, ddf="lme4")
          return(result)        
          
        }          
        break
      }
        
          
      
      # save lmer outcome in rho environmental variable
      rho<-rhoInit(model)     
      
      # calculate asymptotic covariance matrix A
      h <- hessian(function(x) Dev(rho,x), rho$param$vec.matr)
      rho$A <- 2*solve(h)
      #rho$A <- 2*ginv(h)
      
      #Check if A is positive-definite
      isposA<-all(eigen(rho$A)$values>0)
      
      
      
      
      if(!isposA)
      {
        print("Asymptotic covariance matrix A is not positive!")
        #previous code: return ERROR
        #print("ERROR: asymptotic covariance matrix A is not positive!")
        #result$model<-model
        #TABs<-emptyAnovaLsmeansTAB()
        #result$anova.table<-TABs$TAB.fixed
        #result$lsmeans.table<-TABs$TAB.lsmeans
        #result$diffs.lsmeans.table<-TABs$TAB.lsmeans
        #return(result)
      }
      
      
      #calculate ttest and p-values for summary
      if(isTtest)
      {
        tsummary<-calculateTtest(rho, diag(rep(1,nrow(rho$s@coefs))), nrow(rho$s@coefs), method.grad)
        result$ttest<-list(df=tsummary[,"df"], tvalue=tsummary[,"t value"], tpvalue=tsummary[,"p-value"])
        return(result)
      }
        
      
      #calculate lsmeans of differences of LSMEANS of the final model
      if(isLSMEANS || isDiffLSMEANS)
      {
        if(isLSMEANS)
        {
          lsmeans.tab<-calcLSMEANS(model, data, rho, alpha.fixed, test.effs=test.effs, method.grad=method.grad, lsmeansORdiff=TRUE, l)
          result$lsmeans.table<-lsmeans.tab$summ.data
          result$diffs.lsmeans.table<-NULL
        }
        if(isDiffLSMEANS)
        {
          lsmeans.tab<-calcLSMEANS(model, data, rho, alpha.fixed, test.effs=test.effs, method.grad=method.grad, lsmeansORdiff=FALSE, l)
          result$diffs.lsmeans.table<-lsmeans.tab$summ.data
          result$lsmeans.table<-NULL
        }
        return(result)
      }
      
    
      # Calculate  F-test with Satterthwaite's approximation
      # create X design matrix for fixed effects
      X.design <- createDesignMat(model,data)
       
      
      # define full coefficients for model
      #coefs.real<-fullCoefs(model, data)
      #coefs.real<-getDummyCoefs(model, data)
      
      
      #save full coefficients in rho
      #rho$s.test<-coefs.real
      nums.dummy.coefs <- getNumsDummyCoefs(model, data, l)
      rho$nums.zeroCoefs <- nums.dummy.coefs$nums.zeroCoefs
      rho$nums.Coefs <- nums.dummy.coefs$nums.Coefs
      
     
      
      #define the terms that are to be tested
      test.terms <- attr(terms(model),"term.labels")
      
      #initialize anova table
      if(is.first.anova)
      {
        anova.table<-initAnovaTable(model, reduce.fixed)
        is.first.anova<-FALSE
        elim.num<-1
      }
      
      # calculate general set matrix for type 3 hypothesis
      L <- calcGeneralSetForHypothesis(X.design, rho)
      
        
      for(i in 1:length(test.terms))
      {
                         
         result.fstat <- calcFpvalue(test.terms[i], L, model, rho, ddf, method.grad=method.grad)   
         if(!is.null(result.fstat))
         {
            anova.table[which(rownames(anova.table) == test.terms[i]),2]<-result.fstat$denom#round(result.fstat$denom,2)
            anova.table[which(rownames(anova.table) == test.terms[i]),3]<-result.fstat$Fstat#round(result.fstat$Fstat,2)
            anova.table[which(rownames(anova.table) == test.terms[i]),which(colnames(anova.table)=="Pr(>F)")]<-result.fstat$pvalue#round(result.fstat$pvalue,4)
            #if(!is.first.anova)
            #{
            #  anova.model<-anova(model)
            #  anova.table[which(rownames(anova.table) == test.terms[i]),1:4]<-anova.model[which(rownames(anova.model) == test.terms[i]),]
            #}
              
         }
       
      }
      
      #return("abc")
      
     
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
          #model <- eval(substitute(lmer(mf.final, data=data),list(mf.final=mf.final)))
          
          model<-updateModel(model, mf.final, model@dims[["REML"]], l)
         
          anova.table <- resNSelim$anova.table
          elim.num <- elim.num+1
        }        
      }      
    
  }
  

  
  #convert anova table to data frame
  anova.table<-as.data.frame(anova.table)
  anova.table$NumDF<-as.integer(anova.table$NumDF)
  if(!is.null(anova.table$elim.num))
    anova.table$elim.num<-as.integer(anova.table$elim.num)
  
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
  
  tsummary<-calculateTtest(rho, diag(rep(1,nrow(rho$s@coefs))), nrow(rho$s@coefs), method.grad)
  result$ttest<-list(df=tsummary[,"df"], tvalue=tsummary[,"t value"], tpvalue=tsummary[,"p-value"])
  
  
  #update final model
  mf.final <- update.formula(formula(model),formula(model))
  #if(!is.null(l))
  #  model<-eval(substitute(lmer(mf.final, data=data, REML=model@dims[["REML"]], contrasts=l),list(mf.final=mf.final)))
  #else
  #  model<-eval(substitute(lmer(mf.final, data=data, REML=model@dims[["REML"]]),list(mf.final=mf.final)))
  model<-updateModel(model, mf.final, model@dims[["REML"]], contr)
  
  #save model
  model <- as(model,"merLmerTest")
  model@t.pval<-result$ttest$tpvalue
  result$model<-model
  return(result)
}


# generic functions for total analysis on mixed models
#totalAnalysis <- function(model, data,  alpha.rand = 0.05, alpha.fix = 0.05, isFixReduce = FALSE, isRandReduce = FALSE, test.effs=NULL, plot=FALSE, ...) UseMethod("totalAnalysis")


step <- function(model, ddf="Satterthwaite", type=3, alpha.random = 0.1, alpha.fixed = 0.05, reduce.fixed = TRUE, reduce.random = TRUE, lsmeans.calc=TRUE, difflsmeans.calc=TRUE, test.effs=NULL, method.grad="simple",...)
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
  #cat("Call:\n")
  #print(x$call)
  
  if(!is.null(x$rand.table))
  {
    cat("\nRandom effects:\n")  
    printCoefmat(x$rand.table, digits=3 , dig.tst=1  ,tst.ind=c(which(colnames(x$rand.table)=="Chi.DF"),which(colnames(x$rand.table)=="elim.num")), P.values=TRUE, has.Pvalue=TRUE)
  }
  
  if(nrow(x$anova.table)!=0)
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
      printCoefmat(x$anova.table, dig.tst=1, tst.ind=c(1,2), cs.ind=3, digits=3 ,P.values = TRUE, has.Pvalue=TRUE)
      if(!is.null(x$lsmeans.table))
      {
        cat("\nLeast squares means:\n")
        printCoefmat(x$lsmeans.table, dig.tst=1  ,tst.ind=c(1:(which(colnames(x$lsmeans.table)=="Estimate")-1),which(colnames(x$lsmeans.table)=="DF")), digits=3 ,P.values = TRUE, has.Pvalue=TRUE)
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


plot.step<-function(x, ...)
{
  if(!is.null(x$lsmeans.table) && nrow(x$lsmeans.table)>0)
    plotLSMEANS(x$lsmeans.table, x$response, "LSMEANS")     
  if(!is.null(x$diffs.lsmeans.table) && nrow(x$diffs.lsmeans.table)>0)
    plotLSMEANS(x$diffs.lsmeans.table, x$response, "DIFF of LSMEANS")
}

lmer <-
    function(formula, data, family = NULL, REML = TRUE,
             control = list(), start = NULL, verbose = FALSE, doFit = TRUE,
             subset, weights, na.action, offset, contrasts = NULL,
             model = TRUE, x = TRUE, ...)
    {
      mc<-match.call()
      mc[[1]] <- quote(lme4::lmer)
      model<-eval.parent(mc)
      model<-as(model,"merLmerTest")
      #tryCatch(  { result = glm( y~x , family = binomial( link = "logit" ) ) } , error = function(e) { print("test") } )
      t.pval <- tryCatch( {totalAnovaRandLsmeans(model=model, ddf="Satterthwaite", isTtest=TRUE)$ttest$tpvalue}, error = function(e) { NULL })
      if(!is.null(t.pval))
      {
        model@t.pval <-t.pval
      }
      else
      {
        model<-as(model,"mer")
      }
      return(model)
    }


setMethod("anova", signature(object="merLmerTest"),
    function(object,...)  
    {
      mCall <- match.call(expand.dots = TRUE)
      dots <- list(...)
      modp <- if (length(dots))
      sapply(dots, is, "merLmerTest") | sapply(dots, is, "mer") | sapply(dots, is, "lm") else logical(0)
	    if (any(modp)) {
	      return(callNextMethod())
      }
      else
      {
        cnm <- callNextMethod()
        #if(!is.null(ddf)&& ddf=="lme4") 
        #  return(cnm)
        if(("ddf" %in% names(dots)) && dots$"ddf"=="lme4") 
           return(cnm)
        if("ddf" %in% names(dots))
           ddf <- dots$"ddf"
        else
           ddf <- "Satterthwaite"  
        if("method.grad" %in% names(dots))
          method.grad <- dots$"method.grad"
        else
          method.grad <- "simple"
        {
          table <- cnm
          an.table <- tryCatch({totalAnovaRandLsmeans(model=object, ddf=ddf, type=3, isAnova=TRUE, reduce.random=FALSE, reduce.fixed=FALSE, method.grad=method.grad)$anova.table}, error = function(e) { NULL })
          if(!is.null(an.table))
          {
            rnames<-rownames(table)
            if(nrow(an.table)>0)
            {
              table<-as.data.frame(cbind(table$Df, table$"Sum Sq", table$"Mean Sq", an.table[,"F.value"], an.table[,"DenDF"], an.table[,"Pr(>F)"]))
              colnames(table) <- c("Df", "Sum Sq", "Mean Sq", "F value", "Denom", "Pr(>F)")
              dimnames(table) <- list(rnames,
                                      c("Df", "Sum Sq", "Mean Sq", "F value", "Denom", "Pr(>F)"))
            }
            else
              table<-an.table
            attr(table, "heading") <- paste("Analysis of Variance Table with ",ddf," approximation for degrees of freedom")
          }
          
          class(table) <- c("anova", "data.frame")
          return(table)
        }  
      }
      
    })
setMethod("summary", signature(object = "merLmerTest"),
    function(object, ddf="Satterthwaite", ...)
    {
      cl<-callNextMethod()
      if(!is.null(ddf) && ddf=="lme4") return(cl)
      else
      {
         #coefs.satt <- cbind(cl@coefs,totalAnovaRandLsmeans(model=object, ddf=ddf, isTtest=TRUE)$ttest$tpvalue) 
         coefs.satt <- cbind(cl@coefs,object@t.pval) 
         cl@coefs<-coefs.satt
         
#            coefs.satt <- cbind("df"= result$ttest$df, "p value"= result$ttest$tpvalue)
         colnames(cl@coefs)[4]<-"Pr(>|t|)"
      }
      return(as(cl,"summary.merLmerTest"))
    }
) 


# generic functions for ANOVA
#anovaTAB <- function(model, data, ...) UseMethod("anovaTAB")

#setMethod("anova", signature(object = "mer"),
#    function(object, ddf , type, ...)
#    {
#      cat("This the lmerTest version of anova\n\n")
#            
#      result <- totalAnovaRandLsmeans(model=model, ddf=ddf, type=type, isAnova=TRUE, reduce.random=FALSE, reduce.fixed=FALSE, method.grad=method.grad)  
#      anova.table<-result$anova.table
#        
#      if(nrow(anova.table)!=0)
#      {
#        cat("Analysis of Variance Table:\n")    
#        printCoefmat(anova.table, dig.tst=1, tst.ind=c(1,2), cs.ind=3, digits=3 , P.values=TRUE, has.Pvalue=TRUE)
#      }
#      else
#        print(anova.table)
#     }
#)



#setMethod("anova", signature(object = "mer"),
#    function(object, ..., ddf="Satterthwaite" , method.grad="simple")
#      {
#	  mCall <- match.call(expand.dots = TRUE)
#	  dots <- list(...)
#	  modp <- if (length(dots))
#	      sapply(dots, is, "mer") | sapply(dots, is, "lm") else logical(0)
#	  if (any(modp)) {		# multiple models - form table
#	      opts <- dots[!modp]
#	      mods <- c(list(object), dots[modp])
#	      names(mods) <- sapply(as.list(mCall)[c(FALSE, TRUE, modp)],
#				    as.character)
#	      mods <- mods[order(sapply(lapply(mods, logLik, REML = FALSE),
#					attr, "df"))]
#	      calls <- lapply(mods, slot, "call")
#	      data <- lapply(calls, "[[", "data")
#	      if (any(data != data[[1]]))
#		  stop("all models must be fit to the same data object")
#	      header <- paste("Data:", data[[1]])
#	      subset <- lapply(calls, "[[", "subset")
#	      if (any(subset != subset[[1]]))
#		  stop("all models must use the same subset")
#	      if (!is.null(subset[[1]]))
#		  header <-
#		      c(header, paste("Subset", deparse(subset[[1]]), sep = ": "))
#	      llks <- lapply(mods, logLik, REML = FALSE)
#	      Df <- sapply(llks, attr, "df")
#	      llk <- unlist(llks)
#	      chisq <- 2 * pmax(0, c(NA, diff(llk)))
#	      dfChisq <- c(NA, diff(Df))
#	      val <- data.frame(Df = Df,
#				AIC = sapply(llks, AIC),
#				BIC = sapply(llks, BIC),
#				logLik = llk,
#				"Chisq" = chisq,
#				"Chi Df" = dfChisq,
#				"Pr(>Chisq)" = pchisq(chisq, dfChisq, lower = FALSE),
#				row.names = names(mods), check.names = FALSE)
#	      class(val) <- c("anova", class(val))
#             attr(val, "heading") <-
#                  c(header, "Models:",
#                    paste(rep(names(mods), times = unlist(lapply(lapply(lapply(calls,
#                                           "[[", "formula"), deparse), length))),
#                         unlist(lapply(lapply(calls, "[[", "formula"), deparse)),
#                         sep = ": "))
#	      return(val)
#	  }
#	  else { ## ------ single model ---------------------
#            if (length(object@muEta))
#              stop("single argument anova for GLMMs not yet implemented")
#            if (length(object@V))
#              stop("single argument anova for NLMMs not yet implemented")
#            
#            ### added code for anova with ddf approximations         
#              p <- object@dims[["p"]]
#              ss <- (.Call(mer_update_projection, object)[[2]])^2
#              names(ss) <- names(object@fixef)
#              asgn <- attr(object@X, "assign")
#
#              terms <- terms(object)
#              nmeffects <- attr(terms, "term.labels")
#              if ("(Intercept)" %in% names(ss))
#                nmeffects <- c("(Intercept)", nmeffects)
#              ss <- unlist(lapply(split(ss, asgn), sum))
#              df <- unlist(lapply(split(asgn,  asgn), length))
#              ## dfr <- unlist(lapply(split(dfr, asgn), function(x) x[1]))
#              ms <- ss/df                              
#              
#              ## P <- pf(f, df, dfr, lower.tail = FALSE)
#              ## table <- data.frame(df, ss, ms, dfr, f, P)
#              table <- data.frame(df, ss, ms)
#              dimnames(table) <-
#  	          list(nmeffects,
#  		        ## c("Df", "Sum Sq", "Mean Sq", "Denom", "F value", "Pr(>F)"))
#  		        c("Df", "Sum Sq", "Mean Sq"))
#             
#              if(!(ddf %in% c("Satterthwaite","Kenward-Roger")))
#              {
#                f <- ms/(sigma(object)^2)
#                table[,"F value"]<-f
#                if ("(Intercept)" %in% nmeffects)
#                  table <- table[-match("(Intercept)", nmeffects), ]
#                attr(table, "heading") <- "Analysis of Variance Table"
#              }
#              else
#              {
#                 if ("(Intercept)" %in% nmeffects)
#                  table <- table[-match("(Intercept)", nmeffects), ]
#                 res.anova <- totalAnovaRandLsmeans(model=object, ddf=ddf, type=3, isAnova=TRUE, reduce.random=FALSE, reduce.fixed=FALSE, method.grad=method.grad)$anova.table  
#                 table[,"Denom"] <- res.anova[,"DenDF"]
#                 table[,"F value"]<-res.anova[,"F.value"]
#                 table[,"Pr(>F)"]<-res.anova[,"Pr(>F)"]
#                 dimnames(table) <-
#                 list(rownames(table),
#    	            c("Df", "Sum Sq", "Mean Sq", "Denom", "F value", "Pr(>F)"))
#  		           attr(table, "heading") <- paste("Analysis of Variance Table with ",ddf," approximation for degrees of freedom")
#
#              }
#                 
#              
#             
#           
#              
#              class(table) <- c("anova", "data.frame")
#              table
#        
#	  }
#      })
#


#anovaTAB.default<-function(model, data, ...)
#anova <- function(model, ddf, type=3, method.grad="Richardson",...)
#{
#  result <- totalAnovaRandLsmeans(model=model, ddf=ddf, type=type, isAnova=TRUE, reduce.random=FALSE, reduce.fixed=FALSE, method.grad=method.grad)  
#  res <- list(anova.table=result$anova.table)
#  class(res) <- "anova"
#  res
#}

#print.anova <- function(x, ...)
#{

  #if(nrow(x$anova.table)!=0)
  #{
  #  cat("Analysis of Variance Table:\n")    
  #  printCoefmat(x$anova.table, dig.tst=1, tst.ind=c(1,2), cs.ind=3, digits=3 , P.values=TRUE, has.Pvalue=TRUE)
  #}
  #else
  #  print(x$anova.table)
#}


# generic functions for random effects
#randTAB <- function(model, data, ...) UseMethod("randTAB")


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
  #if(x$isCorr)
  #  print(x$rand.table)
  #else
  if(!is.null(x$rand.table))
    printCoefmat(x$rand.table, digits=3 , dig.tst=1  ,tst.ind=c(which(colnames(x$rand.table)=="Chi.DF"),which(colnames(x$rand.table)=="elim.num")), P.values=TRUE, has.Pvalue=TRUE)        
}



# generic functions for LSMEANS
#lsmeansTAB <- function(model, data, test.effs=NULL, plot=FALSE , ...) UseMethod("lsmeansTAB")


#lsmeansTAB.default<-function(model, data, test.effs=NULL, plot=FALSE , ...)
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






#setGeneric("sigma", function(object, ...) standardGeneric("sigma"))

#setMethod("sigma", signature(object = "mer"),
#          function (object, ...) {
#              dd <- object@dims
#        if(!dd[["useSc"]]) return(1)
#	      object@deviance[[if(dd[["REML"]]) "sigmaREML" else "sigmaML"]]
#          })

#setMethod("summary", signature(object = "summary.mer"), function(object) object)

#setMethod("coef", signature(object = "summary.mer"),
#          function(object, ...) object@coefs)

#setMethod("summary", signature(object = "mer"),
#    function(object, ddf="Satterthwaite", ...)
#      {
#          
#          
#          #calculate df, t-value and p-values
#          coefs.satt<-NULL
#          if(ddf=="Satterthwaite")
#          {
#            #added code: refit model with contr.SAS
#            #options(contrasts=c(unordered="contr.SAS", ordered="contr.poly"))
#            #warning("\n model has been refitted with contrasts=contr.SAS \n")
#            #object<-update(object, .~., object@dims[["REML"]])
# 
#            result <- totalAnovaRandLsmeans(model=object, ddf=ddf, isTtest=TRUE)  
#            coefs.satt <- cbind("df"= result$ttest$df, "p value"= result$ttest$tpvalue)
#          }
#          
#                    
#          REML <- object@dims[["REML"]]
#          fcoef <- fixef(object)
#          vcov <- vcov(object)
#          corF <- vcov@factors$correlation
#          dims <- object@dims
#          coefs <- cbind("Estimate" = fcoef, "Std. Error" = corF@sd) #, DF = DF)
#          llik <- logLik(object, REML)
#          dev <- object@deviance
#          mType <- if((non <- as.logical(length(object@V)))) "NMM" else "LMM"
#          if (gen <- as.logical(length(object@muEta)))
#              mType <- paste("G", mType, sep = '')
#          mName <- switch(mType, LMM = "Linear", NMM = "Nonlinear",
#                          GLMM = "Generalized linear",
#                          GNMM = "Generalized nonlinear")
#	  method <- {
#	      if (mType == "LMM")
#		  if(REML) "REML" else "maximum likelihood"
#	      else
#		  paste("the", if(dims[["nAGQ"]] == 1) "Laplace" else
#			"adaptive Gaussian Hermite",
#			"approximation")
#	  }
#          AICframe <- data.frame(AIC = AIC(llik), BIC = BIC(llik),
#                                 logLik = as.vector(llik),
#                                 deviance = dev[["ML"]],
#                                 REMLdev = dev[["REML"]],
#                                 row.names = "")
#          if (is.na(AICframe$REMLdev)) AICframe$REMLdev <- NULL
#          varcor <- VarCorr(object)
#          REmat <- formatVC(varcor)
#          if (is.na(attr(varcor, "sc")))
#              REmat <- REmat[-nrow(REmat), , drop = FALSE]
#
#          if (nrow(coefs) > 0) {
#              if (!dims[["useSc"]]) {
#                  coefs <- coefs[, 1:2, drop = FALSE]
#                  stat <- coefs[,1]/coefs[,2]
#                  pval <- 2*pnorm(abs(stat), lower = FALSE)
#                  coefs <- cbind(coefs, "z value" = stat, "Pr(>|z|)" = pval)
#              } else {
#                  stat <- coefs[,1]/coefs[,2]
#                  ##pval <- 2*pt(abs(stat), coefs[,3], lower = FALSE)
#                  coefs <- cbind(coefs, "t value" = stat) #, "Pr(>|t|)" = pval)
#                  #added code: inclusion of p-values, df and t-values SAS
#                  coefs <- cbind(coefs, coefs.satt)
#
#              }
#          } ## else : append columns to 0-row matrix ...
#          new("summary.mer",
#              object,
#              methTitle = paste(mName, "mixed model fit by", method),
#              logLik = llik,
#              ngrps = sapply(object@flist, function(x) length(levels(x))),
#              sigma = sigma(object),
#              coefs = coefs,
#              vcov = vcov,
#              REmat = REmat,
#              AICtab= AICframe
#              )
#      })## summary()


### code from lme4 package
## This is modeled a bit after  print.summary.lm :
#printMer <- function(x, digits = max(3, getOption("digits") - 3),
#                     correlation = TRUE, symbolic.cor = FALSE,
#                     signif.stars = getOption("show.signif.stars"), ...)
#{
#    so <- summary(x)
#    REML <- so@dims[["REML"]]
#    llik <- so@logLik
#    dev <- so@deviance
#    dims <- x@dims
#
#    cat(so@methTitle, "\n")
#    if (!is.null(x@call$formula))
#        cat("Formula:", deparse(x@call$formula),"\n")
#    if (!is.null(x@call$data))
#        cat("   Data:", deparse(x@call$data), "\n")
#    if (!is.null(x@call$subset))
#        cat(" Subset:", deparse(x@call$subset),"\n")
#    print(so@AICtab, digits = digits)
#
#    cat("Random effects:\n")
#    print(so@REmat, quote = FALSE, digits = digits, ...)
#
#    ngrps <- so@ngrps
#    cat(sprintf("Number of obs: %d, groups: ", dims[["n"]]))
#    cat(paste(paste(names(ngrps), ngrps, sep = ", "), collapse = "; "))
#    cat("\n")
#    if (is.na(so@sigma))
#  cat("\nEstimated scale (compare to 1):",
#            sqrt(exp(so@deviance[["lr2"]])/so@dims[["n"]]), "\n")
#    if (nrow(so@coefs) > 0) {
#  ### here the code is changed in order to put p-values
#  cat("\nFixed effects:\n")
#  if(ncol(so@coefs)==5)
#  printCoefmat(so@coefs, zap.ind = 3, #, tst.ind = 4
#		     digits = digits, signif.stars = signif.stars, has.Pvalue=TRUE ,P.values = TRUE)
#  else
#  printCoefmat(so@coefs, zap.ind = 3, #, tst.ind = 4
#  	     digits = digits, signif.stars = signif.stars)

#  if(correlation) {
#	    corF <- so@vcov@factors$correlation
#	    if (!is.null(corF)) {
#		p <- ncol(corF)
#		if (p > 1) {
#		    rn <- rownames(so@coefs)
#		    rns <- abbreviate(rn, minlength=11)
#		    cat("\nCorrelation of Fixed Effects:\n")
#		    if (is.logical(symbolic.cor) && symbolic.cor) {
#			corf <- as(corF, "matrix")
#			dimnames(corf) <- list(rns,
#					       abbreviate(rn, minlength=1, strict=TRUE))
#			print(symnum(corf))
#		    }
#		    else {
#			corf <- matrix(format(round(corF@x, 3), nsmall = 3),
#				       ncol = p, dimnames = list(rns,
#					       abbreviate(rn, minlength=6)))
#			corf[!lower.tri(corf)] <- ""
#			print(corf[-1, -p, drop=FALSE], quote = FALSE)
#		    }
#		}
#	    }
#	}
#   }
#    invisible(x)
#}

#setMethod("print", "mer", printMer)
#setMethod("show", "mer", function(object) printMer(object))
#####setMethod("print", "summary.mer", printMer)

#setMethod("deviance", signature(object = "summary.mer"), function(object) object@deviance)
#setMethod("logLik", signature(object = "summary.mer"), function(object) object@logLik)
#setMethod("vcov", signature(object = "summary.mer"), function(object) object@vcov)
#setMethod("summary", signature(object = "summary.mer"), function(object) object)

#setClass("summary.mer",                 # Additional slots in a summary object
#         representation(           
#  		methTitle = "character",
#			logLik= "logLik",
#			ngrps = "integer",
#			sigma = "numeric", # scale, non-negative number
#			coefs = "matrix",
#			vcov = "dpoMatrix",
#			REmat = "matrix",
#			AICtab= "data.frame"),
 #        contains = "mer")