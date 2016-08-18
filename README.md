
  README file for:

  Program:  adaptive.ci.r
  Revision date: August 18, 2016


  The adaptive.ci function computes 95% confidence intervals limits
  for any single coefficient in a linear model having fixed effects.
  The observations are adaptively weighted by default.

  This program is written in the base R language so no other packages
  are needed.

  Author: T. W. O'Gorman            email:  twogorman@gmail.com

  Language: R  version 3.2.2

  Operating System:  Windows 7 Professional

  Packages needed: None
  

    Reference: O'Gorman, T. W. Reducing the width of confidence intervals
               for the difference between two population means by inverting
               adaptive tests. Statistical Methods in Medical Research,
               in press, 2016.

  The function is called by:

    adaptive.ci(ci.df, depvar, complete, reduced, civar,
               equalwts, nblocks, s1, s2, s3, details)

  The function will return a vector that contains the lower and upper limits.

  The first five arguments are required; the last six are optional.

   The function arguments are:
 
     1)  ci.df is an R data frame that includes all of the variables that
         will be used in the analysis.
     2)  depvar is the name of the dependent variable in the model.
     3)  complete is a character string that specifies the full model
         including the confidence interval variable. This model uses
         the same syntax as the lm() function.
     4)  reduced is a character string that species the reduced model 
         that does not include the confidence interval variable.
     5)  civar is the character string that specifies the confidence
         inverval variable.
     6)  if equalwts = 1 the observations are given equal weights,
         if equalwts = 0 the observations are given adaptive weights,
         which is the default.
     7)  The number of blocks used in the approximation. The default is
         is 2 blocks, so that the last permutation test will use 8000
         permutations.
     8)  s1 is one of three random number seeds. It can be any integer
         in the range of 1 to 30268.
     9)  s2 is another random number seed in the range of 1 to 30306.
    10)  s3 is the last random number seed in the range of 1 to 30322.
    11)  if details = 1 (default) the details of the search will be printed.
         if details = 0 only the limits will be returned by the function,
         but no output will be printed.
         if details = 2 the weights given to the unpermuted data in the last
         block will be printed, in addition to the details = 1 output.
 
   Notes:
 
     1) The first five arguments are required. If you want to find the
        limits using adaptive weighting, and 2 blocks will give sufficient
        accuracy, and the default random number seeds are acceptable,
        then you only need to enter the first five arguments.
     2) The data frame cannot contain missing values for any variables
        used in the complete model.
     3) This function calls the adonetailp, adaptiveweights, rootcdf,
        cdfhat, and shufflewh functions.
     4) This function is written in base R.  No packages are required. 
 
   Examples:
 
     If blood pressure data is used to create a data frame (bp.df) that has
     blood pressure (bp), age (age), and a treatment indicator (group),
     then we could find an adaptive 95% confidence interval for the
     treatment effect by using this code:
 
       source("adaptive.ci.r")
       bplimits <- adaptive.ci(ci.df=bp.df, depvar=c("bp"),
                      complete=c("bp~group"), reduced=c("bp~1),
                      civar=c("group") )
 
     The vector bplimits will contain the lower and upper limits.
 
     We could expand this example if we needed to include age as a
     covariate. If we wanted to use 3 blocks and we wanted to specify
     random number seeds we would use:
 
       source("adaptive.ci.r")
       bplimits <- adaptive.ci(ci.df=bp.df, depvar=c("bp"),
                      complete=c("bp~age + group"), reduced=c("bp~age),
                      civar=c("group"), nblocks=3,
                      s1 = 3682, s2 = 27812, s3 = 12973 )
 
   These R functions were carefully checked and I believe
   that the functions are correct.  However, the author is not
   responsible for any errors that may still exist in the code.
 
   Please report any issues concerning this code to T. W. O'Gorman via 
   email at twogorman@gmail.com
 
