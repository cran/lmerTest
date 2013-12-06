\name{lsmeans}
\alias{lsmeans}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
 Calculates Least Squares Means and Confidence Intervals for the factors of a fixed part of mixed effects model of lmer object. 
}
\description{
 Produces a data frame which resembles to what SAS software gives in proc mixed statement. The approximation of degrees of freedom is Satterthwate's.
}
\usage{
lsmeans(model, test.effs=NULL, method.grad="simple",...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{model}{
   linear mixed effects model (lmer object).
}


\item{test.effs}{
    
     charachter vector specyfying the names of terms to be tested. If NULL all the terms are tested.
}

\item{method.grad}{
   approximation method for the grad function, which is used in calculation of denominator degrees of freedom. Could be "simple" or "Richardson". "simple" is the default one (and the faster one).
}

  \item{\dots}{
     other potential arguments.
}
}

\value{
 Produces Least Squares Means (population means) table with p-values and Confidence intervals. 
 }
\references{
 doBy package, gplots package
 }
\author{
Alexandra Kuznetsova, Per Bruun Brockhoff, Rune Haubo Bojesen Christensen
}

\note{
  For construction of the contrast matrix popMatrix function from doBy package was used. 
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{step}}, \code{\link{rand}}, \code{\link{difflsmeans}}
}
\examples{


#import lme4 package and lmerTest package
library(lmerTest)

#specify lmer model
m1 <- lmer(Informed.liking ~ Gender*Information +(1|Consumer), data=ham)

#calculate least squares means for interaction Gender:Information
lsmeans(m1, test.effs="Gender:Information")

#import TVbo data from lmerTest package
data(TVbo) 

m <- lmer(Coloursaturation ~ TVset*Picture + (1|Assessor), data=TVbo)
plot(lsmeans(m))
lsmeans(m, test.effs="TVset")
}
% R documentation directory.
%%\keyword{ ~kwd1 }
%%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
