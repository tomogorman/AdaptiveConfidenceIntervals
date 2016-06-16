#
#  This adaptive.ci function computes 95% confidence intervals limits
#  for any single coefficient in a linear model having fixed effects.
#  The observations are adaptively weighted by default. 
#
#  Reference: O'Gorman, T. W. Reducing the width of confidence intervals
#             for the difference between two population means by inverting
#             adaptive tests.  2016.
#
#  The function arguments are:
#
#    1)  ci.df is a R data frame that includes all of the variables that
#        will be used in the analysis.
#    2)  depvar is the name of the dependent variable in the model.
#    3)  complete is a character string that specifies the full model
#        including the confidence interval variable.
#    4)  reduced is a character string that species the model that does
#        not include the confidence interval variable.
#    5)  civar is the character string that specifies the confidence
#        inverval variable.
#    6)  if equalwts = 1 the observations are given equal weights,
#        if equalwts = 0 the observations are given adaptive weights,
#        which is the default.
#    7)  The number of blocks used in the approximation. The default is
#        is 2 blocks, so that the last permutation test will use 8000
#        permutations.
#    8)  s1 is one of three random number seeds. It can be any integer
#        in the range of 1 to 30000.
#    9)  s2 is another random number seed in the range of 1 to 30000.
#   10)  s3 is the last random number seed in the range of 1 to 30000.
#   11)  if details = 1 (default) the details of the search will be printed.
#        if details = 0 only the limits will be returned.
#
#  Notes:
#
#    1) The first five arguments are required. If you want to find the
#       limits using adaptive weighting, and 2 blocks will give sufficient
#       accuracy, and the default random number seeds are acceptable,
#       then you only need to enter the first five arguments.
#    2) The data frame cannot contain missing values for any variables
#       used in the complete model.
#    3) This function calls the adonetailp, adaptiveweights, rootcdffast,
#       cdfhat, and shufflewh functions.
#
#  Examples:
#
#    If blood pressure data is used to create a data frame (dfbp) that has
#    blood pressure (bp), the age (age), and a treatment indicator (group),
#    then we could find an adaptive 95% confidence interval by this code:
#
#    source("adaptive.ci.r")
#    depvar   <- c("bp");
#    complete <- c("bp~group");
#    reduced  <- c("bp~1);
#    civar    <- c("group");
#    bplimits <- adaptive.ci(dfbp,depvar,complete,reduced,civar)
#
#    The vector bplimits will contain the lower and upper limits.
#
#    We could expand this example if we needed to include age as a
#    covariate. We would use:
#
#    source("adaptive.ci.r")
#    depvar   <- c("bp");
#    complete <- c("bp~age+group")
#    reduced  <- c("bp~age)
#    civar    <- c("group");
#    bplimits <- adaptive.ci(dfbp,depvar,complete,reduced,civar)
#
#  These functions were carefully checked on May 24, 2016, and I believe
#  that these functions are correct.  However, the author is not
#  responsible for any errors that may still exist in the code.
#
#  Please report any problems with this code to T. W. O'Gorman via 
#  email at twogorman@gmail.com
#

adaptive.ci <- function(ci.df, depvar, complete, reduced, civar,
               equalwts=0, nblocks=2, s1=7832, s2=25933, s3=19857, details=1)                {

depvar  <- gsub(" ","",depvar)
reduced <- gsub(" ","",reduced)
civar   <- gsub(" ","",civar)
n <- length(ci.df[,depvar])
limits <- double(2)
                  
# The value of r specifies the number of permutations for the first block.

r <- c(4000)

if(details == 1) {
cat("\n")
cat("  Input to adaptive.ci function:","\n\n")
cat("    Data set:",dataset,"\n")
cat("    Depandent variable: ",depvar, "\n")
cat("    Complete model: ",complete,"\n")
cat("    Reduced model : ",reduced,"\n")
cat("    Confidence Interval for = : ",civar, "\n")
cat("    Maximum number of blocks = ", nblocks,"\n")
cat("    Random seeds = ",s1,s2, s3,"\n")
cat("    Number of permutations for first blocks = ",r,"\n\n")
                  }

#     The next three lines reset the random number seeds.

    s1 <- (171*s1) %% 30269
    s2 <- (172*s2) %% 30307
    s3 <- (170*s3) %% 30323

alpha   <- 0.05
alphad2 <- alpha/2
probit  <- function(p) qnorm(p)
probitad2 <- probit(alphad2)

compmodel <- lm(as.formula(complete), data = ci.df)
beta <- summary(compmodel)$coefficients[civar,1]
se <- summary(compmodel)$coefficients[civar,2]
tunperm <- summary(compmodel)$coefficients[civar,3]
ndf <- df.residual(compmodel)

t01  <- qt(0.01, ndf)
t025 <- qt(0.025, ndf)
t04  <- qt(0.04, ndf)
t96  <- qt(0.96, ndf)
t975 <- qt(0.975, ndf)
t99  <- qt(0.99, ndf)

e01 <- beta + t01*se
e025<- beta + t025*se
e04 <- beta + t04*se
e96 <- beta + t96*se
e975<- beta + t975*se
e99 <- beta + t99*se

if(details==1) cat("\n"," Traditional 95% Confidence Interval = ("
     ,round(e025,5),", ", round(e975,5),")", "\n","\n")

ciadj.df <- ci.df

if(details == 1) {
if(equalwts ==1)     {
      cat("Equal Weights Used to Compute Limits; No adaptation","\n")
    } else {
      cat("Adaptive Weights Used to Compute Limits","\n")
                      }
                 }

for(lefttail in 1:2 )                  {

if( (lefttail == 1) && (details ==1) ) cat("\n\n"," Begin Left Tail","\n\n\n")
if( (lefttail != 1) && (details ==1) ) cat("\n\n"," Begin Right Tail","\n\n\n")

  if(lefttail == 2) {
                          e01 <- e99
                          e04 <- e96
                          }

repeat {

ciadj.df[,depvar] <- ci.df[,depvar] - e01*ci.df[,civar]

plist <- adonetailp(ciadj.df, depvar, complete, reduced, civar, r,
                    s1, s2, s3, n, lefttail, equalwts)
  p01 <- plist[[1]]
  s1  <- plist[[2]]
  s2  <- plist[[3]]
  s3  <- plist[[4]]

ciadj.df[,depvar] <- ci.df[,depvar] - e04*ci.df[,civar]

plist <- adonetailp(ciadj.df, depvar, complete, reduced, civar, r,
                    s1, s2, s3, n, lefttail, equalwts)
  p04 <- plist[[1]]
  s1  <- plist[[2]]
  s2  <- plist[[3]]
  s3  <- plist[[4]]

if( ( lefttail == 1) && (details ==1) ) {
cat("left tail lower estimate = ", e01, " p-value = ", p01, "\n")
cat("left tail upper estimate = ", e04, " p-value = ", p04, "\n\n")
                        }
if( ( lefttail != 1) && (details ==1) ) {
cat("right tail upper estimate = ", e01, " p-value = ", p01, "\n")
cat("right tail lower estimate = ", e04, " p-value = ", p04, "\n\n")
                        }
if( (p01 < alphad2) && (alphad2 < p04) ) break

shift <- (e04 - e01)/2

if(p04 < alphad2) {
  e01 <- e01 + shift
  e04 <- e04 + shift
                       }

if(p01 > alphad2) {
  e01 <- e01 - shift
  e04 <- e04 - shift
                       }

    }

slopel <- (probit(p04) - probit(p01))/(e04 - e01)

l<- double(nblocks)
p<- double(nblocks)

l[1] <- e04-(probit(p04) - probitad2)/slopel
if(details==1) cat(" First Interpolated Estimate = ", l[1], "\n\n")
if(details ==1) {
  cat("          Estimate                                   Estimate","\n")
  cat(" Block  Used in Test    R    p-value     SE(p)   Updated after test",
      "\n\n")
                     }
rnext <- r
for (i in 2:nblocks) {
  rnext <- 2*rnext
  se <- sqrt(0.025*(1-0.025)/rnext)
  im1 <- i - 1
  yadj <- ci.df[,depvar] - l[im1]*ci.df[,civar]
  ciadj.df[,depvar] <- yadj
  plist <- adonetailp(ciadj.df, depvar, complete, reduced, civar, rnext,
                      s1, s2, s3, n, lefttail, equalwts )
  p[im1] <- plist[[1]]
  s1  <- plist[[2]]
  s2  <- plist[[3]]
  s3  <- plist[[4]]
#cat("seeds after reestimation calls = ", s1, s2, s3,"\n","\n") 

  l[i]  <- (l[im1] + (l[im1] - (probit(p[im1]) - probitad2)/slopel))/2
if(details == 1) { 
  cat(format(i,width=4),
      format(round(l[im1],5), nsmall=5, width=13),
      format(rnext,width= 7),
      format(round(p[im1],5), nsmall=5, width=10),
      format(round(se,5),     nsmall=5, width=10),
      format(round(l[i],5),   nsmall=5, width=12),"\n\n")
                     }
                         }

  if(i==nblocks) limits[lefttail] <- l[i]
                                         }

return(limits)
                                                                 }



#
#  The adonetailp function produces one-tailed p-values
#  based on t test statistics from an adaptively weighted model.
#  If the variable lefttail=1 then the left tail p-value will
#  be computed, otherwise the right tail p-value will be computed.
#
#  If the reduced model contains only an intercept term, the weights
#  will be permuted, rather than computed from the permuted data.
#
#  These functions were carefully checked on May 17, 2016.
#

adonetailp <- function(dfadonetail, depvar, complete, reduced, indvar, r,
                       s1, s2, s3, n, lefttail, equalwts) {
localwt             <- double(n)
dfadonetail$w2      <- double(n)
dfadonetail$w2perm  <- double(n)
  
redu <- lm(as.formula(reduced), data = dfadonetail)
yhat <- predict(redu)
yresidual <- residuals(redu)

if(equalwts == 1) {
          dfadonetail$w2 <- rep(1,n)
  } else {
          dfadonetail    <- adaptiveweights(dfadonetail,reduced)
          }

compu <- lm(as.formula(complete), data = dfadonetail, weights = dfadonetail$w2)

tunperm <- summary(compu)$coefficients[indvar,3]

e <- 0

simple <- paste(depvar,c("~1"),sep="")

if(reduced==simple) {
  countnum <- double(n)
  ynew     <- double(n)
  yresshuf <- double(n)
  dfadonetail$w2perm   <- double(n)
  countnum <- c(1:n)
  for (k in 1:r) {
    permlist <- shufflewh(countnum,s1,s2,s3,n)
    s1 <- permlist[[2]]
    s2 <- permlist[[3]]
    s3 <- permlist[[4]]

    permnums <- permlist[[1]]

    yresshuf <- yresidual[permnums]
    dfadonetail$w2perm <- dfadonetail$w2[permnums]

#    for (i in 1:n) {
#      index <- permnums[i]
#      yresshuf[i] <- yresidual[index]
#      dfadonetail$w2perm[i] <- dfadonetail$w2[index]
#                    }

    dfadonetail[,depvar] <- yhat + yresshuf

 compw <-lm(as.formula(complete), data = dfadonetail,
            weights = dfadonetail$w2perm)
      tperm <- summary(compw)$coefficients[indvar,3]
      if( (lefttail == 1) & (tperm >= tunperm) ) e  <-  e + 1
      if( (lefttail != 1) & (tperm <= tunperm) ) e  <-  e + 1
                    }
                           }

if(reduced != simple)         {
for (k in 1:r) {
  sresidlist <- shufflewh(yresidual,s1,s2,s3,n)
  s1 <- sresidlist[[2]]
  s2 <- sresidlist[[3]]
  s3 <- sresidlist[[4]]
  dfadonetail[,depvar] <- yhat+sresidlist[[1]]
if(equalwts == 1) {
          dfadonetail$w2 <- rep(1,n)
  } else {
          dfadonetail    <- adaptiveweights(dfadonetail,reduced)
          }
  compw  <- lm(as.formula(complete), data = dfadonetail,
                 weights = dfadonetail$w2)

  tperm  <- summary(compw)$coefficients[indvar,3]
      if( (lefttail == 1) & (tperm >= tunperm) ) e  <-  e + 1
      if( (lefttail != 1) & (tperm <= tunperm) ) e  <-  e + 1
                   }
                                  }
p <- (e+1)/(r+1)
plist <- list(p,s1,s2,s3)
return(plist)
}



adaptiveweights <- function(dfweights,reduced) {
red <- lm(as.formula(reduced), data=dfweights)
resid <- residuals(red)
#                               compute traditional quantiles

q25 <- quantile(resid,0.25,type=6)
q75 <- quantile(resid,0.75,type=6)
iqr <- q75 - q25
sigmat <- iqr/1.349

#                               compute bandwidth (h)

n <- length(resid)
h <- 1.587*sigmat*n^(-0.33333333333)

minr  <- min(resid)
maxr  <- max(resid)
range <- maxr - minr
lower <- minr - range/2
upper <- maxr + range/2
tol   <- 0.0000001*iqr
cdf25 <- rootcdffast(resid,h,0.25,lower,maxr,tol)
cdf50 <- rootcdffast(resid,h,0.50,minr, maxr,tol)
cdf75 <- rootcdffast(resid,h,0.75,minr,upper,tol)
sigma <- (cdf75-cdf25)/1.349

#                               compute adaptive weights

s <- (resid-cdf50)/sigma
tolerance <- 100*tol
w  <- double(n)
dfweights$w2 <- double(n)

residdh <- resid/h

for (i in 1:n) {
  phi  <- pnorm(residdh[i] - residdh)
  fhat <- mean(phi)
  z    <- qnorm(fhat)
  if( abs(s[i]) >= tolerance ) w[i] <- z/s[i] else w[i] <- 1 
                }
  dfweights$w2 <- w*w
                
return(dfweights)
#print(dfweights)

  }


rootcdffast <- function(x,h,p,xlow,xhigh,tolerance) {
#
#  This root finding function uses bisection for the first
#  few iterations, then switches to the false position method for the
#  remaining iterations.  Convergence is achieved when the cdf
#  is within a small tolerance around the desired percentile.
#
  flow  <-  cdfhat(x,h,xlow)
  fhigh  <-  cdfhat(x,h,xhigh)  
  if ( (flow >p) | (fhigh < p) ) {
    cat("Stop in rootcdffast flow,fhigh=",flow,fhigh);q() }
  for (i in seq(1:60) ) {
    if( i < 8 ) xmiddle <- (xlow+xhigh)/2
             else xmiddle <- ((fhigh-p)*xlow-(flow-p)*xhigh)/(fhigh-flow)
    fmiddle <- cdfhat(x,h,xmiddle)
    if( fmiddle == p ) break
    if( fmiddle < p  ) {xlow   <-  xmiddle;  flow  <-  fmiddle}
    if( fmiddle > p  ) {xhigh  <-  xmiddle; fhigh  <-  fmiddle}
    if(  (abs(fmiddle - p)) <= tolerance ) break
    if( i==60 ){print (" stop in rootcdffast with over 60 iterations");q()}
                          }
  return(xmiddle)
  }


cdfhat <- function(xvector,h,xpoint){
#
#   This function computes the smooth estimate of the c.d.f. at xpoint.
#
  phi <- pnorm((xpoint-xvector)/h)
  cdf <- mean(phi)
  return(cdf)
  }


shufflewh <- function(y, s1, s2, s3, n) {

#     Compute n-1 uniform random numbers using the Wichmann-Hill method
#     and then perform Durstenfeld shuffle on y.
#     Return seeds in list yandseeds.

k <- c(1:n)
yold <- double(n)
yold <- y
for (i in n:2) {
    s1 <- (171*s1) %% 30269
    s2 <- (172*s2) %% 30307
    s3 <- (170*s3) %% 30323
    u  <- (s1/30269.0 + s2/30307.0 + s3/30323.0) %% 1.000000
  itrade      <- floor(i*u + 1) 
  ktemp       <- k[itrade]
  k[itrade]   <- k[i]
  k[i]        <- ktemp
                }
for (i in 1:n)  {
  iplace <- k[i]
  y[iplace] <- yold[i]
                 }
yandseeds <- list(y, s1, s2, s3)
return(yandseeds)
                                         }


