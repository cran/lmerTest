\name{anova-methods}
\docType{methods}
\alias{anova,ANY-method}
\alias{anova,merModLmerTest-method}
\title{Methods for function \code{anova} in package \pkg{lmerTest}}
\description{
Methods for Function \code{anova} in Package \pkg{lmerTest}
}
%\section{Methods}{
%\describe{

%\item{\code{signature(object = "merModLmerTest")}}{
%   returns the data frame with F test of type 3 hypothesis for each term in formula and corresponding p-values. For multiple arguments, model comparison statistics. The definition of type 3 hypothesis comes from SAS theory proc mixed. The approximations for degrees of freedom by default are Satterthwaite's (ddf="Satterthwaite"). If ddf="Kenward-Roger", then the Kenward-Roger's approximations are calculated using \code{KRmodcomp} function from \code{pbkrtest} package. If ddf="lme4" then the anova table that comes from \code{lme4} package is returned. 
%method.grad  is an argument for the approximation method for the grad function, which is used in calculation of denominator degrees of freedom. method.grad could be "simple" or "Richardson". method.grad="simple" is the default one and the fastest one, method.grad="Richardson" gives more accurate results.
%\code{\link{anova}}
%}

\usage{\S4method{anova}{merModLmerTest}(object, \dots , ddf="Satterthwaite", type=3, method.grad="simple")}
\arguments{
\item{object}{object of class "merModLmerTest"}
\item{\dots}{object of class "merModLmerTest". Then the model comparison statistisc will be calculated}
\item{ddf}{By default the Satterthwaite's approximation to degrees of freedom is calculated. If ddf="Kenward-Roger", then the Kenward-Roger's approximation is calculated using \code{KRmodcomp} function from \pkg{pbkrtest} package. If ddf="lme4" then the anova table that comes from \pkg{lme4} package is returned.}
\item{type}{type of hypothesis to be tested. Could be type=3 or type=1 (The definition comes from SAS theory) }
\item{method.grad}{"simple" is the default one and the fastest one. "Richardson" gives more accurate results. This argument only applies to ddf="Satterthwaite"}
}


%}}

\references{
SAS Technical Report R-101 1978 Tests of Hypotheses in Fixed-Effects Linear Models \emph{Copyright (C)}  (SAS Institute Inc., Cary, NC, USA)

Goodnight, J.H. 1976 General Linear Models Procedure (S.A.S. Institute, Inc.)

Schaalje G.B., McBride J.B., Fellingham G.W. 2002 Adequacy of approximations to distributions of test Statistics in complex mixed linear models 
}

\keyword{methods}
\examples{
#import lmerTest package
library(lmerTest)

#import ham data from lmerTest package
data(ham)

m.ham <- lmer(Informed.liking ~ Product*Information*Gender 
+ (1|Consumer), data=ham)

# type 3 anova table with  denominator degrees of freedom
# calculated based on Satterthwaite's approximation
anova(m.ham)

# type 1 anova table with denominator degrees of freedom
# calculated based on Satterthwaite's approximation
anova(m.ham, type=1)

# type3 anova table with additional F statistics and denominator degrees of freedom
# calculated based on Kenward-Roger's approximation
anova(m.ham, ddf="Kenward-Roger")

# anova table, that is returned by lme4 package
anova(m.ham, ddf="lme4")
}

