# AdaptiveConfidenceIntervals


  This adaptive.ci function computes 95% confidence intervals limits
  for any single coefficient in a linear model having fixed effects.
  The observations are adaptively weighted by default.

  This program is written in the base R language so no other packages
  are needed.

  Author: T. W. O'Gorman            email:  twogorman@gmail.com

  Language: R  version 3.2.2

  Operating System:  Windows 7 Professional

  Packages needed: None

  Before attempting to use this function, the user should read   

    Reference: O'Gorman, T. W. Reducing the width of confidence intervals
               for the difference between two population means by inverting
               adaptive tests.  2016.

  The function is called by:

    adaptive.ci(ci.df, depvar, complete, reduced, civar,
               equalwts, nblocks, s1, s2, s3, details)

  The function will return a vector that contains the lower and upper limits.

  The first five arguments are required; the last six are optional.
  The arguments are:


    1)  ci.df is a R data frame that includes all of the variables that
        will be used in the analysis.

    2)  depvar is the name of the dependent variable in the model.

    3)  complete is a character string that specifies the full model
        including the confidence interval variable.

    4)  reduced is a character string that species the model that does
        not include the confidence interval variable.

    5)  civar is the character string that specifies the confidence
        inverval variable.

    6)  if equalwts = 1 the observations are given equal weights,
        if equalwts = 0 the observations are given adaptive weights,
        which is the default.

    7)  The number of blocks used in the approximation. The default is
        is 2 blocks, so that the last permutation test will use 8000
        permutations.

    8)  s1 is one of three random number seeds. It can be any integer
        in the range of 1 to 30000.

    9)  s2 is another random number seed in the range of 1 to 30000.

   10)  s3 is the last random number seed in the range of 1 to 30000.

   11)  if details = 1 (default) the details of the search will be printed.
        if details = 0 only the limits will be returned.


  Notes:

    1) The first five arguments are required. If you want to find the
       limits using adaptive weighting, and 2 blocks will give sufficient
       accuracy, and the default random number seeds are acceptable,
       then you only need to enter the first five arguments.
    2) The data frame cannot contain missing values for any variables
       used in the complete model.
    3) This function calls the adonetailp, adaptiveweights, rootcdffast,
       cdfhat, and shufflewh functions.

  Examples:

    If blood pressure data is used to create a data frame (dfbp) that has
    blood pressure (bp), the age (age), and a treatment indicator (group),
    then we could find an adaptive 95% confidence interval for the
    group effect by using this code:

    source("adaptive.ci.r")
    depvar   <- c("bp");
    complete <- c("bp~group");
    reduced  <- c("bp~1);
    civar    <- c("group");
    bplimits <- adaptive.ci(dfbp,depvar,complete,reduced,civar)

    The vector bplimits will contain the lower and upper limits.

    We could expand this example if we needed to include age as a
    covariate we would use:

    source("adaptive.ci.r")
    depvar   <- c("bp");
    complete <- c("bp~age+group")
    reduced  <- c("bp~age)
    civar    <- c("group");
    bplimits <- adaptive.ci(dfbp,depvar,complete,reduced,civar)

  These functions were carefully checked on May 24, 2016, and I believe
  that these functions are correct.  However, the author is not
  responsible for any errors that may still exist in the code.

  Please report any problems with this code to T. W. O'Gorman via 
  email at twogorman@gmail.com

