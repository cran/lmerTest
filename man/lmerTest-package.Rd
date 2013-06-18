\name{lmerTest-package}
\alias{lmerTest-package}
\alias{lmerTest}
\docType{package}
\title{
 The package performs different kinds of tests on lmer objects, such as F tests of type 3 hypothesis for the fixed part, 
log-likelihood ratio tests for the random part, least squares means (population means) and differences of least squares means for the factors of the fixed part with corresponding plots. The package also provides with a function step, that preforms backward elimination of non-significant effects, starting from the random effects, and then fixed ones. 
}
\description{
 The package provides anova function, that gives data frame similar to what gives lme4 package but with p-values calculated from F statistics of type 3 hypothesis. There are two options for denominator degrees of freedom of F statistics: "Satterthwaite" and "Kenward-Roger". The calculation of Kenward-Roger's approximation is based on function from pbkrtest package, the calculation of Satterthwaite's approximation is based on SAS proc mixed theory (see reference).The type 3 hypothesis (marginal) is calculated according to SAS theory (SAS Institute Inc., 1978). The package also provides summary function, which gives the same as lme4 package summary function but with p-values added for the t-test (based on Satterthwaite approximation for denominator degrees of freedom).
The tests on random effects are performed using log-likelihood ratio tests with one degree of freedom, which means, testing one effect in a time.
}

\details{
\tabular{ll}{
Package: \tab lmerTest\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2012-01-10\cr
License: \tab GPL\cr
}

The calculation of statistics for the fixed part was developed according to SAS Proc Mixed Theory (see reference).   
If the elimination of the random part is required, the model with REML=FALSE is updated automatically to REML=TRUE

}
\author{
Alexandra Kuznetsova <alku@imm.dtu.dk>, Per Bruun Brockhoff, Rune Haubo Bojesen Christensen 
}
\references{
Tests of Hypotheses in Fixed-Effects Linear Models Copyright 1978 by SAS Institute Inc., Cary, NC, USA.
J.H. Goodnight. General Linear Model Procedure S.A.S. Institute, Inc.
Approximations to Distributions of Test Statistics in Complex Mixed Linear Models Using SAS Proc MIXED
}
%%~~ Optionally other standard keywords, one per line, from file KEYWORDS in the R documentation ~~
%%~~ directory ~~
%%\keyword{ package }
%%\seealso{
%%~~ Optional links to other man pages, e.g. ~~
%%~~ \code{\link[<pkg>:<pkg>-package]{<pkg>}} ~~
%%}
\examples{
#import lme4 package and lmerTest package
library(lmerTest)

#import ham data from lmerTest package
data(ham)

# an object of class merLmerTest
m <- lmer(Informed.liking ~ Gender+Information+Product +(1|Consumer), data=ham)

# gives summary of lmer object. The same as of class summary.mer but with 
# additional p-values calculated based on Satterthwate's approximations
summary(m)

# anova table the same as of class mer but with additional F statistics and 
# p-values calculated based on Satterthwaite's approximations
anova(m)

# anova table the same as of class mer but with additional F statistics and  
# p-values calculated based on Kenward-Roger's approximations
anova(m, ddf="Kenward-Roger")

# anova table the same as of class mer 
anova(m, ddf="lme4")

# backward elimination of non-significant effects of model m
st<-step(m)

plot(st)

}

