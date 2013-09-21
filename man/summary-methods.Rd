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
 summary of the results of linear mixed effects model fitting of object. Returns the same as summary function of "merMod" class but  with an additional column \code{t.pval}, which gives \code{p}-values for the \code{t} test of fixed terms with Satterthwaite approximation for degrees of freedom.
      \code{\link{summary}}
}
}}
\keyword{methods}


\examples{
(fm1 <- lmer(Reaction ~ Days + (Days | Subject), sleepstudy))

# will give you an additional column with p values for the t test
summary(fm1)
}