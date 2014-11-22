##deviance function
## parameter  -  covariate parameters
devfun5.vars <- function(fm,  reml = TRUE) 
{
  stopifnot(is(fm, "merMod"))
  
  vlist <- sapply(fm@cnms, length)
  
  pp <- fm@pp$copy()
  
  resp <- fm@resp$copy()
  np <- length(pp$theta)
  nf <- length(fixef(fm)) 
  if (!isGLMM(fm)) 
    np <- np + 1L
  n <- nrow(pp$V)
  
  
  ff <- updateModel(fm, .~., getREML(fm), 
                    attr(model.matrix(fm),"contrasts"), 
                    devFunOnly.lmerTest.private = TRUE) 
  
  envff <- environment(ff)
  
  if (isLMM(fm)) {
    ans <- function(vars) {
      stopifnot(is.numeric(vars), length(vars) == np)
      
      sigsq  <- vars[np]
      thpars <- Vv_to_Cv(vars, n = vlist, s = sqrt(sigsq))
      ff(thpars)
      
      #       thpars <- Vv_to_Cv(pars, n = vlist, s = sqrt(sigma2))
      dev <- envff$pp$ldL2() + (envff$resp$wrss() + envff$pp$sqrL(1))/sigsq + n * 
        log(2 * pi * sigsq)      
      if(reml){
        p <- ncol(envff$pp$RX())
        dev <- dev + 2*determinant(envff$pp$RX())$modulus - p * log(2 * pi * sigsq)              
      }
      
      return(dev)     
    }
  }
  
  #attr(ans, "thopt") <- pp$theta
  class(ans) <- "devfun5.vars"
  ans
}




## not calling the C code
vcovJSStheta2.var <- function(fm)
{
  stopifnot(is(fm, "merMod"))
  
  vlist <- sapply(fm@cnms, length)
  
  pp <- fm@pp$copy()
  
  
  resp <- fm@resp$copy()
  np <- length(pp$theta)
  nf <- length(fixef(fm))
  if (!isGLMM(fm)) 
    np <- np + 1L
  
 
  
  ff2 <- updateModel(fm, .~., getREML(fm), 
                     attr(model.matrix(fm),"contrasts"), 
                     devFunOnly.lmerTest.private = TRUE) 
  
  envff2 <- environment(ff2)
  
  if (isLMM(fm)) {
    ans <- function(Lc, vars) {     
      stopifnot(is.numeric(vars), length(vars) == np)      
      sigma2  <- vars[np]
      thpars <- Vv_to_Cv(vars, n = vlist, s = sqrt(sigma2))
      ff2(thpars)
      
      
      #.Call("lmer_Deviance", pp$ptr(), resp$ptr(), thpars[-np], PACKAGE = "lme4")      
      vcov_out <- sigma2 * tcrossprod(envff2$pp$RXi()) 
      
      return(as.matrix(Lc %*% as.matrix(vcov_out) %*% t(Lc)))        
    }
  } 
  
  class(ans) <- "vcovJSStheta2.var"
  ans
}







test <- FALSE

if(test){
  
  
  ## test 
  (fm1 <- lmer(Reaction ~ Days + (1| Subject), sleepstudy))
  vv <- as.data.frame(VarCorr(fm1))[, "vcov"]
  dd <- devfun5.vars(fm1)
  dd(vv)
  
  ## Evaluate gradient:
  require(numDeriv)
  grad(dd, vv)
  
  ## Evaluate Hessian:
  TOL <- 1e-5
  (Hess <- hessian(dd, vv))
  stopifnot(all.equal(2*solve(Hess), 
                      matrix(c(255801, -1146, -1145.93, 11459), ncol=2, nrow=2),
                      tol = TOL))
  
  m <- lmer(Informed.liking ~ Product*Information*Gender 
            + (1|Consumer) + (1|Product:Consumer)  + (1|Information:Consumer), data=ham)
  vv <- as.data.frame(VarCorr(m))[, "vcov"]
  dd <- devfun5.vars(m)
  dd(vv)
  Vv_to_Cv(vv, n = sapply(m@cnms, length), s = getME(m, "sigma") )
  Cv_to_Vv(getME(m, "theta"), n = sapply(m@cnms, length), s = getME(m, "sigma") )
  getME(m, "theta")
  (Hess <- hessian(dd, vv))
  2*solve(Hess)
  
  TOL <- 1e-3
  stopifnot(all.equal(2*solve(Hess), 
                      matrix(c(0.1422, 0.002592, -0.03554, -0.01037,
                               0.002592, 0.007363, -0.00368, -0.00518,
                               -0.03554, -0.00368, 0.05824, 0.002592,
                               -0.01037, -0.00518, 0.002592, 0.02073), ncol=4, nrow=4),
                      tol = TOL))
  
  ## check models with slopes
  (fm2 <- lmer(Reaction ~ Days + (1 + Days| Subject), sleepstudy))
  vv <- as.data.frame(VarCorr(fm2))[, "vcov"] ## does not work
  vv <- Cv_to_Vv(getME(fm2, "theta"), n = sapply(fm2@cnms, length), 
                 s = getME(fm2, "sigma") )
  dd <- devfun5.vars(fm2)
  dd(vv)
  
  (Hess <- hessian(dd, vv))
  
  ## check that asymptotic var cov matrix is as in SAS
  stopifnot(all.equal(2*solve(Hess), 
                      matrix(c(83395, -2688.41, 105.21, -2058.08,
                               -2688.41, 2178.88, -136.11, 324.96,
                               105.21, -136.11, 218.51, -72.2134,
                               -2058.08, 324.96, -72.2134, 5957.61), 
                             ncol=4, nrow=4),
                      tol = TOL))
  
  
  all.equal(Vv_to_Cv(Cv_to_Vv(getME(fm2, "theta"), n = sapply(fm2@cnms, length), 
                              s = getME(fm2, "sigma") ), n = sapply(fm2@cnms, length), 
                     s = getME(fm2, "sigma")),
            getME(fm2, "theta"), check.attributes = FALSE)
  
  
  
  m.carrots <- lmer(Preference ~ sens1*Age*sens2*Homesize
                    +(1 + sens2|Consumer) + (1|product), data=carrots)
  
  
  vv <- Cv_to_Vv(getME(m.carrots, "theta"), n = sapply(m.carrots@cnms, length), 
                 s = getME(m.carrots, "sigma") )
  dd <- devfun5.vars(m.carrots)
  dd(vv)
  
  (Hess <- hessian(dd, vv))
  
  stopifnot(all.equal(2*solve(Hess), matrix(c(0.001697, 0.000021, 1.367E-6, 
                                              -7.88E-7, -0.00018,
                                              0.000021, 0.000029, 7.057E-7,
                                              1.393E-7, 4.733E-8, 
                                              1.367E-6, 7.057E-7, 1.996E-6, 
                                              1.185E-7, -0.00001,
                                              -7.88E-7, 1.393E-7, 1.185E-7,
                                              0.000586, -0.00002,
                                              -0.00018, 4.733E-8, -0.00001,
                                              -0.00002, 0.002137), 
                                            nrow = 5, ncol = 5),
                      tol = TOL))
  
  ## check the asymptotic var cor for model with slopes as factors
  load(system.file("testdata","bugSummaryData.RData", package="lmerTest"))
  lmer4 <- lmer(cog_Mid ~ (1|item) + (vowel - 1|speaker) + Language*vowel*sex, 
                data = data1.frame)

  all.equal(Vv_to_Cv(Cv_to_Vv(getME(lmer4, "theta"), n = sapply(lmer4@cnms, length), 
                              s = getME(lmer4, "sigma") ), 
                     n = sapply(lmer4@cnms, length), 
                     s = getME(lmer4, "sigma")),
            getME(lmer4, "theta"), check.attributes = FALSE)
  
}
