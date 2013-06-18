\name{merLmerTest-class}
\title{Mixed Model Representations}
\docType{class}
% Classes
\alias{merLmerTest-class}
%\alias{coef.mer-class}
\alias{summary.merLmerTest-class}
%\alias{ranef.mer-class}
% Methods
\alias{anova,merLmerTest-method}

\alias{summary,merLmerTest-method}
\alias{summary,summary.merLmerTest-method}

\description{The \code{merLmerTest} \emph{contains} mer class of lme4 package and overloads \code{anova} and \code{summary} functions.  The \code{summary.merLmerTest} class represents the summary of objects of class \code{merLmerTest}.


  The \code{"summary.merLmerTest"} class \emph{contains} the \code{"summary.mer"}
}

\section{Objects from the Class}{Objects can be created via the
  \code{\link{lmer}}  functions.
} 
 
\section{Methods}{
anova \code{signature(object = "merLmerTest", \dots)}: returns the data frame with F test of type 3 hypothesis for each term in formula and corresponding p-values. For multiple arguments, model comparison statistics. The definition of type 3 hypothesis comes from SAS theory proc mixed. The approximations for degrees of freedom by default are Satterthwaite's (ddf="Satterthwaite"). If ddf="Kenward-Roger", then the Kenward-Roger's approximations are calculated using \code{KRmodcomp} function from \code{pbkrtest} package. If ddf="lme4" then the anova table that comes from \code{lme4} package is returned. 
method.grad  is an argument for the approximation method for the grad function, which is used in calculation of denominator degrees of freedom. method.grad could be "simple" or "Richardson". method.grad="simple" is the default one and the fastest one, method.grad="Richardson" gives more accurate results.
      \code{\link{anova}}
}
 
\seealso{
  \code{\link{lmer}()}
}
\examples{

(m <- lmer(Reaction ~ Days + (1|Subject) + (0+Days|Subject),
             data = sleepstudy))

# anova table with additional F statistics and denominator degrees of freedom
# calculated based on Satterthwaite's approximation
anova(m)

# anova table with additional F statistics and denominator degrees of freedom
# calculated based on Kenward-Roger's approximation
anova(m, ddf="Kenward-Roger")

# anova table, that is returned by lme4 package
anova(m, ddf="lme4")

# summary of merLmerTest object. Returns the same as mer object but with an 
#additional column of p values for the t test. 
summary(m)

}
\keyword{classes}
