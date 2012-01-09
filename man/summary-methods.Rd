\name{summary-methods}
\docType{methods}
\alias{summary,merModLmerTest-method}
\title{  Methods for Function \code{summary} in Package \pkg{lmerTest}}
\description{
 Methods for function \code{summary} in package \pkg{lmerTest} 
}
\section{Methods}{
\describe{



\item{\code{signature(object = "merModLmerTest", \dots)}}{
 summary of the results of linear mixed effects model fitting of object. Returns a similar to summary function of "merMod" class but with additional columns with the names \code{"df"}, \code{"t value"} and \code{"Pr(>t)"} representing degrees of freedom, \code{t}-statistics and \code{p}-values respectively  calculated based on Satterthwaite's approximation.
      \code{\link{summary}}
}
}}
\keyword{methods}


\examples{
(fm1 <- lmer(Reaction ~ Days + (Days | Subject), sleepstudy))

# will give you an additional column with p values for the t test
summary(fm1)

#will give the summary of lme4 package
summary(fm1, ddf="lme4")
}