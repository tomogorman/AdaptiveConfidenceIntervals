
  README file for:

  Program:  adaptive.ci.r
  Revision date: October 11, 2016


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

       adaptive.ci <- function(ci.df, model, equalwts=0, nblocks=2,
         s1=7832, s2=25933, s3=19857, details=1)
 
   The function arguments are:
 
     1)  ci.df is an R data frame that includes all of the variables that
         will be used in the analysis.
     2)  model is a character string that specifies the full model
         including the confidence interval variable. This model uses
         the same syntax as the lm() function.
         Note: The confidence interval will be computed for the last
               variable in the model character string.
     3)  if equalwts = 1 the observations are given equal weights,
         if equalwts = 0 the observations are given adaptive weights,
         which is the default.
     4)  The number of blocks used in the approximation. The default is
         is 2 blocks, so that the last permutation test will use 8000
         permutations.
     5)  s1 is one of three random number seeds. It can be any integer
         in the range of 1 to 30268.
     6)  s2 is another random number seed in the range of 1 to 30306.
     7)  s3 is the last random number seed in the range of 1 to 30322.
     8)  if details = 1 (default) the details of the search will be printed.
         if details = 0 the limits will be returned by the function, but
         no output will be printed.
         if details = 2 the weights given to the unpermuted data at the last
         block will be printed, in addition to the details = 1 output.

  Notes:

     1) The first two arguments are required. If you want to find the
        limits using adaptive weighting, and 2 blocks will give sufficient
        accuracy, and the default random number seeds are acceptable,
        then you only need to enter the first two arguments.
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
       bplimits <- adaptive.ci(ci.df=bp.df, model=c("bp~group") )
 
     The vector bplimits will contain the lower and upper limits.
 
     We could expand this example if we needed to include age as a
     covariate. If we wanted to use 3 blocks and we wanted to specify
     the three random number seeds we could use:
 
       source("adaptive.ci.r")
       bplimits <- adaptive.ci(ci.df=bp.df, complete=c("bp ~ age + group"),
                      nblocks=3, s1 = 3682, s2 = 27812, s3 = 12973 )
 
   Note that the group variable was specified as the last variable in the
   model because we wanted the confidence interval for the group effect.

 
   These R functions were carefully checked and I believe
   that the functions are correct.  However, the author is not
   responsible for any errors that may still exist in the code.
 
   Please report any issues concerning this code to T. W. O'Gorman via 
   email at twogorman@gmail.com
 
