
### fill a row for the random matrix
fillRowRandTable <- function(term, rand.table, elim.num, reduce.random)
{
#   nrow.term <- which(rownames(rand.table)==names(term))
#   um"] <- elim.num  
#   return(rand.table)
  
  #
  nrow.term <- which(gsub(" ","",rownames(rand.table))==names(term$term))
  rand.table[nrow.term, "Chi.sq"] <- term$chisq
  rand.table[nrow.term, "Chi.DF"] <- term$chisq.df
  rand.table[nrow.term, "p.value"] <- term$pv
  if(reduce.random)
  {
    rand.table[nrow.term, "elim.num"] <- elim.num 
    if(term$chisq.df==2 && elim.num!=0)
    {
      rand.table.upd <- matrix(0,ncol=4, nrow=1)
      colnames(rand.table.upd) <- c("Chi.sq", "Chi.DF", "elim.num", "p.value")
     
      rownames(rand.table.upd) <- paste(paste(rep(" ", 
                                            substring.location( names(term$term), term$term[[1]]$gr.part)$first-2), collapse=""), term$term[[1]]$gr.part, collapse="")
      if(nrow.term==nrow(rand.table))
      {
        rand.table <- rbind(rand.table, rand.table.upd)
        #rownames(rand.table)[1] <- term$term
      }
      else
      {
        rnames <- c(rownames(rand.table)[1:nrow.term], rownames(rand.table.upd),rownames(rand.table)[(nrow.term+1):nrow(rand.table)])
        rand.table <- rbind(rand.table[1:nrow.term,], rand.table.upd, rand.table[(nrow.term+1):nrow(rand.table),])  
        rownames(rand.table) <- rnames 
      } 
    }
      
  }
      
  rand.table  
}

### update table for random terms
updateRandTable <- function(infoForTerm, rand.table, 
                            elim.num=0, reduce.random)
{  
  if(!is.null(infoForTerm$term))
  {   
    rand.table <- fillRowRandTable(infoForTerm, rand.table, elim.num, reduce.random)  
  } 
  else
  {
    for(i in 1:length(infoForTerm))
    {
       iterm <- infoForTerm[i]
       rand.table <- fillRowRandTable(iterm[[1]], rand.table, 
                                     elim.num, reduce.random) 
    }
  }  
  rand.table
}



### get the random terms out of a model formula
# getRandTerms <- function(fmodel)
# {
#   terms.fm <- attr(terms(fmodel),"term.labels")
#   ind.rand.terms <- which(unlist(lapply(terms.fm,
#                             function(x) substring.location(x, "|")$first))!=0)
#   return(terms.fm[ind.rand.terms])
# }

### function get the slope part of a random term
findSlopePart <- function(term)
{
 
  
  sub.loc.div <- substring.location(term," |")
  slopepart <- substring2(term,1,sub.loc.div$first)
  grouppart <- substring2(term,sub.loc.div$last, nchar(term))
  parts <- unlist(strsplit(slopepart, split=c("[[:punct:]]")))
  parts <- gsub(" ","", parts , fixed=TRUE)
  parts
}

getSlGrParts <- function(term)
{
  randTerm.split <- unlist(strsplit(term, "\\|")) 
  sl.part <- sapply( strsplit(randTerm.split[1], split=c("[[:punct:]]")), 
                     function(x) gsub(" ","", x , fixed=TRUE))
  gr.part <- gsub(" ","", randTerm.split[2] , fixed=TRUE)
  return(list(sl.part=sl.part, gr.part=gr.part))
}


getRowNamesForRandTable <- function(term)
{
  sl.part <- getSlGrParts(term)$sl.part
  gr.part <- getSlGrParts(term)$gr.part
  if(length(sl.part) > 1 || sl.part !="1")
  {
     sl.part <- sl.part[! sapply(sl.part, function(x) x=="1" || x=="0")]
     listterms <- function(x)
     {
       lst <- list(sl.part=x, gr.part=gr.part)
       lst
     }
     la <- lapply(sl.part, listterms)
     names(la) <- sapply(sl.part, 
              function(x) paste(c(x, gr.part), collapse=":"), USE.NAMES=FALSE)     
  }
  else
  {
     lst <- list(sl.part=sl.part, gr.part=gr.part) 
     la <- list(lst)
     names(la) <- gr.part     
  } 
  
  la
}

### initialize table for random terms
initRandTable <- function(names, reduce.random)
{ 
  
  if(reduce.random)
  {
    rand.table <- matrix(0,ncol=4, nrow=length(names))
    colnames(rand.table) <- c("Chi.sq", "Chi.DF", "elim.num", "p.value")
    rownames(rand.table) <- names
  }
  else
  {
    rand.table <- matrix(0,ncol=3, nrow=length(names))
    colnames(rand.table) <- c("Chi.sq", "Chi.DF", "p.value")
    rownames(rand.table) <- names
  }
  return(rand.table)
}


# return formula with the simplified random structure
fmElimRandTerm <- function(rnm, rand.terms, fm)
{
  nSameGr <- length(rand.terms[sapply(rand.terms, 
                function(x) (getSlGrParts(x)$gr.part==rnm[[1]]$gr.part))])
  term.simplify <-  rand.terms[sapply(rand.terms, 
                     function(x) (getSlGrParts(x)$gr.part==rnm[[1]]$gr.part) && 
                                  (rnm[[1]]$sl.part %in% findSlopePart(x)))] 
  sub.loc.div <- substring.location(term.simplify," |")
  slopepart <- gsub(" ","", substring2(term.simplify, 1, sub.loc.div$first) , 
                    fixed=TRUE)
  mf.term.simpl <- as.formula(paste(fm[2],fm[1],paste(slopepart, "-", 
                                                      rnm[[1]]$sl.part), sep=""))
  mf.term.simpl <- update.formula( mf.term.simpl, mf.term.simpl)
  if(rnm[[1]]$sl.part != "1" && as.character(mf.term.simpl[3])!="1 - 1" 
                                 && nSameGr==1)
  {
   
    fm[3] <- paste(fm[3], "-", paste("(",term.simplify,")", sep=""), "+" , 
                       paste("(",as.character(mf.term.simpl[3]),
                             "|",rnm[[1]]$gr.part, ")", sep=""))
  }
  else
  {
    fm[3] <- paste(fm[3], "-", paste("(",term.simplify,")", sep=""), sep="")
  }
 
  mf.final <-  as.formula(paste(fm[2],fm[1],fm[3], sep=""))
  mf.final <- update.formula( mf.final, mf.final)
  mf.final
}

getRandTermsTable <- function(rand.terms)
{
  rand.terms.table <- NULL
  for(rand.term in rand.terms)
  {
    rand.terms.table <- c(rand.terms.table, getRowNamesForRandTable(rand.term))
  }
  rand.terms.table
}

############################################################################    
#save info (pvalues, chisq val, std, var...) for term
############################################################################    
saveInfoForTerm <- function(term, chisq, chisq.df, pv)
{
  term.info <- NULL 
  term.info$pv <- pv
  term.info$chisq <- chisq
  term.info$chisq.df <- chisq.df
  term.info$term <- term
  return(term.info)
}

############################################################################
#get formula for model 
############################################################################
# getFormula <- function(model, withRand=TRUE)
# {
#   fmodel <- formula(model)
#   terms.fm <- attr(terms.formula(fmodel),"term.labels")
#   ind.rand.terms <- which(unlist(lapply(terms.fm,function(x) substring.location(x, "|")$first))!=0)
#   terms.fm[ind.rand.terms] <- unlist(lapply(terms.fm[ind.rand.terms],function(x) paste("(",x,")",sep="")))
#   fm <- paste(fmodel)
#   if(withRand)
#     fm[3] <- paste(terms.fm,collapse=" + ")
#   else
#     fm[3] <- paste(terms.fm[-ind.rand.terms],collapse=" + ")
#   
#   if(fm[3]=="")
#     fo <- as.formula(paste(fm[2],fm[1],1, sep=""))
#   else
#     fo <- as.formula(paste(fm[2],fm[1],fm[3], sep=""))
#   return(fo)
# }


# refitLM <- function(obj, l="contr.SAS") {
#   #   cl <- match.call()
#   #   bits <- getME(obj,c("X","y","offset"))
#   #   w <- weights(obj)
#   #   m <- if (!all(w==1)) {
#   #     with(bits,lm.fit(X,y,offset=offset, contrasts=l))
#   #   } else {
#   #     with(bits,lm.wfit(X,y,w,offset=offset))
#   #   }
#   #   class(m) <- "lm"
#   #   m$offset <- offset
#   #   m$call <- cl
#   #   m
#   
#   #mm <- model.frame.fixed(obj)
#   mm <- model.frame(obj)
#   colnames(mm)[1] <- "y"
#   fo <- getFormula(obj, withRand=FALSE)# formula(obj,fixed.only=TRUE)
#   if(fo != as.formula(.~.))
#   {
#     inds <-  names(l) %in% attr(terms(fo), "term.labels")
#     #update contrast l
#     l <- l[inds]
#   }
#   fo <- update(fo, y ~ .)
#   lm(fo, data=mm, contrasts = l)
# }


# ### compare mixed model versus fixed
# compareMixVSFix <- function(model, mf.final, data, name.term, rand.table, alpha, elim.num, reduce.random)
# {
#   #library(nlme)
#   #return(NULL)
#   #mframe <- model.frame(mf.final, data=data, na.action=na.pass)
#   #if(length(which(names(mframe) %in% names(data)==FALSE))!=0)
#   # {
#   #   data$response <- mframe[,1]
#   #   fm <- paste(mf.final)
#   #   fm[2] <- "response"
#   #   mf.final <-  as.formula(paste(fm[2],fm[1],fm[3], sep=""))
#   #   mf.final <- update.formula(mf.final,mf.final)       
#   # }
#   
#   #model.red <- gls(model = mf.final, data=data, method = "REML", na.action=na.omit)
#   #model.red  <-  lm(formula(model,fixed.only=TRUE), data=model.frame.fixed(model))#lm(model, data=summary(model,"lme4")@frame)
#   model.red <- refitLM(model)
#   
#   l.fix <- -2*logLik(model)[1]
#   l.red <- -2*logLik(model.red, REML=TRUE)[1]
#   
#   p.chisq <- 1 - pchisq (l.red -l.fix ,1)
#   infoForTerm <- saveInfoForTerm(name.term, l.red -l.fix, 1, p.chisq)
#   #detach(package:nlme)
#   
#   if(infoForTerm$pv > alpha)
#   {    
#     rand.table <- updateRandTable(infoForTerm, rand.table, elim.num=elim.num, reduce.random=reduce.random)
#     model.last <- model.red
#   }
#   else
#   {
#     rand.table <- updateRandTable(infoForTerm, rand.table, reduce.random=reduce.random)
#     model.last <- model    
#   }
#   return(list(model=model.last, TAB.rand=rand.table))   
# }


### check if there are no random terms in the model
# checkPresRandTerms <- function(mf.final)
# {
#   sub.loc.rand <- substring.location(paste(mf.final)[3], "|")
#   if(sub.loc.rand$first==0 && sub.loc.rand$last==0)
#     return(FALSE)
#   return(TRUE)
# }


### eliminate NS random terms 
elimRandEffs <- function(model, data, alpha, reduce.random, l)
{
  isInitRand <- TRUE
  elim.num <- 1
  stop <- FALSE
  while(!stop)
  {
    fmodel <- formula(model)    
    rand.terms <- sapply(getRandTerms(fmodel), function(x) substr(x,2,nchar(x)-1)
                         , USE.NAMES = FALSE) 
    rand.terms.table <- getRandTermsTable(rand.terms)
    
    if(isInitRand)
    {
      rand.table <- initRandTable(names(rand.terms.table), reduce.random)
      isInitRand <- FALSE
    }      
    
    fm <- paste(fmodel)
    pv.max <- 0
    infoForTerms <- vector("list", length(rand.terms.table))
    names(infoForTerms) <- names(rand.terms.table)  
    
    for(i in 1:length(rand.terms.table))
    {
      rnm <- rand.terms.table[i]
      fm <- paste(fmodel)
      mf.final <- fmElimRandTerm(rnm, rand.terms, fm)
      is.present.rand <- checkPresRandTerms(mf.final)
      
      # no more random terms in the model
      if(!is.present.rand)
      {
        return(compareMixVSFix(model, mf.final, data, rnm, rand.table, alpha, elim.num, reduce.random))
        
      } 
      model.red <- updateModel(model, mf.final, getME(model, "is_REML"), l)
      anova.red <- suppressMessages(anova(model, model.red))
      infoForTerms[[names(rnm)]] <- saveInfoForTerm(rnm, anova.red$Chisq[2], anova.red$"Chi Df"[2] , anova.red$'Pr(>Chisq)'[2])
      if((anova.red$'Pr(>Chisq)'[2] >= pv.max) && reduce.random)
      { 
        pv.max <- anova.red$'Pr(>Chisq)'[2]
        infoForTermElim <- infoForTerms[[names(rnm)]]
        model.final <- model.red 
        #if(anova.red$'Pr(>Chisq)'[2]==1)
        #  break
      }
    }
    
    if(!reduce.random)
    {
      rand.table <- updateRandTable(infoForTerms, rand.table, reduce.random=reduce.random)
      model.last <- model
      break
    }
    
    #rand.terms.upd <- getRandTerms(formula(model.final))
    
    if(infoForTermElim$pv > alpha)
    {
      rand.table <- updateRandTable(infoForTermElim, rand.table, elim.num, reduce.random)
      elim.num <- elim.num+1      
    }
    else
    {
      rand.table <- updateRandTable(infoForTerms, rand.table, reduce.random=reduce.random)
      model.last <- model
      break
    }
    
    model <- model.final  
    
  }
  return(list(model=model.last, TAB.rand=rand.table))
}

