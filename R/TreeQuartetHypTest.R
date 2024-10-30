#' Power divergence statistic of Cressie & Read
#'
#' Computes any of the family of power-divergence statistics for multinomial data of \insertCite{CR84;textual}{MSCquartets},
#' to compare observed and expected counts. Includes Likelihood Ratio and Chi-squared statistics as special cases.
#'
#' @references
#' \insertRef{CR84}{MSCquartets}
#'
#' @param obs  observation vector
#' @param expd  expected vector
#' @param lambda  statistic parameter (e.g., 0=Likelihood Ratio, 1=Chi-squared)
#'
#' @return value of statistic
#'
#' @examples
#' obs=c(10,20,30)
#' expd=c(20,20,20)
#' powerDivStat(obs,expd,0)
#'
#' @export
powerDivStat <- function(obs,
                         expd,
                         lambda) {
  if (is.vector(obs) == FALSE ||
      is.vector(expd) == FALSE ||
      length(obs) != length(expd) ||
      min(obs) < 0 ||
      min(expd) < 0 ||
      is.vector(lambda) == FALSE || length(lambda) != 1) {
    # Return an error and stop if arguments are misspecified. If any expectations are 0 then the function fails.
    stop(
      "Invalid arguments: obs and expd must be non-negative vectors of the same length. lambda must be real."
    )
  }
  else{
    if (all.equal(sum(obs), sum(expd)) != TRUE)
    {
      diff = sum(obs) - sum(expd)
      warning(paste0(
        "sum(obs)-sum(expd) = ",
        diff,
        ", which is not within tolerance of 0."
      ))
    }
    if (all(obs == expd) == TRUE) {
      # If observations and expectations are equal, there is no need for an algorithm for computing the power-divergence statistic since the power-divergence statistic must be 0.
      stat  <- 0
    } else {
      if (lambda %in% c(0, 1,-1 / 2,-1,-2, 2 / 3) == FALSE) {
        warning(
          "Non-standard choice of lambda. Standard choices are lambda=0,1,-1/2,-1,-2,2/3."
        )
      }
      if (lambda == 0) {
        # If lambda=0, compute the likelihood ratio statistic.
        stat <- sum(ifelse(obs == 0, 0, 2 * obs * log(obs / expd)))
      } else if (lambda == -1) {
        # If lambda=-1, compute the modified likelihood ratio statistic.
        if (min(obs) == 0) {
          # If one of the observations is 0, then the modified likelihood statistic is undefined and an error is returned, advising the user to choose a different lambda.
          stop("Some observations are 0. Use any lambda other than lambda=-1.")
        } else {
          stat <- sum(2 * expd * log(expd / obs))
        }
      } else {
        # If lambda is anything other than 0 or -1, then the function in the limit is not required.
        stat <-
          sum((2 / (lambda * (lambda + 1))) * obs * ((obs / expd) ^ lambda -
                                                       1))
      }
    }
    return(stat)
  }
}

#########################################################

# Error function
#
#@param x argument
#@return function value
#
erf <- function(x) {
  2 * pnorm(x * sqrt(2)) - 1
}

#########################################################

# Inverse error function
#
#@param x argument
#@return function value
#
erf.inv <- function(x) {
  qnorm((x + 1) / 2) / sqrt(2)
}

#########################################################

#' Modified Struve function
#'
#' This function is used in computing the probability density for
#' Model T1. The code is closely based on the \code{I0L0} function implemented in Python for the package
#' RandomFieldUtils, which was previously on CRAN up to 12/2022).
#'
#'@param x  function argument
#'
#'@return value of negative modified Struve function function
#'
#'@export
M0 <- function(x) {
  if (x<0) {
    r <- NA
  } else if (x>=0&&x<16) {
    g2 <- c(0.52468736791485599138,-0.35612460699650586196,0.20487202864009927687,-0.10418640520402693629,0.4634211095548429228*10^(-1),-0.1790587192403498630*10^(-1),0.597968695481143177*10^(-2),-0.171777547693565429*10^(-2),0.42204654469171422*10^(-3),-0.8796178522094125*10^(-4),0.1535434234869223*10^(-4),-0.219780769584743*10^(-5),0.24820683936666*10^(-6),-0.2032706035607*10^(-7),0.90984198421*10^(-9),0.2561793929*10^(-10),-0.710609790*10^(-11),0.32716960*10^(-12),0.2300215*10^(-13),-0.292109*10^(-14),-0.3566*10^(-16),0.1832*10^(-16),-0.10*10^(-18),-0.11*10^(-18))
    r <- -0.5*g2[1]
    ac <- acos((6*x-40)/(x+40))
    for (i in 2:24) {
      r <- r-g2[i]*cos((i-1)*ac)
    }
  } else {
    g3 <- c(2.00326510241160643125,0.195206851576492081*10^(-2),0.38239523569908328*10^(-3),0.7534280817054436*10^(-4),0.1495957655897078*10^(-4),0.299940531210557*10^(-5),0.60769604822459*10^(-6),0.12399495544506*10^(-6),0.2523262552649*10^(-7),0.504634857332*10^(-8),0.97913236230*10^(-9),0.18389115241*10^(-9),0.3376309278*10^(-10),0.611179703*10^(-11),0.108472972*10^(-11),0.18861271*10^(-12),0.3280345*10^(-13),0.565647*10^(-14),0.93300*10^(-15),0.15881*10^(-15),0.2791*10^(-16),0.389*10^(-17),0.70*10^(-18),0.16*10^(-18))
    r <- 0.5*g3[1]
    x2 <- x^2
    ac <- acos((800-x2)/(288+x2))
    for (i in 2:24) {
      r <- r+g3[i]*cos((i-1)*ac)
    }
    r <- -r*2/(pi*x)
  }
  return(r)
}


#########################################################

#' Probability density function for Model T1
#'
#' Value of probability density function for Model T1 of \insertCite{MAR19;textual}{MSCquartets}, Proposition 5.2.
#'
#'@references
#'\insertRef{MAR19}{MSCquartets}
#'
#'@param x  statistic value (e.g., likelihood ratio statistic, or other power divergence statistic)
#'@param mu0  parameter of density function
#'
#'@return value of density function
#'
#'@seealso \code{\link{T3density}}
#'
#'@export
T1density <-
  function(x,
           mu0) {
    (1 / 4) * exp(-(1 / 2) * x) * (sqrt(2 / (pi * x)) * (1 + erf(mu0 / sqrt(2))) -
                                     exp(-(1 / 2) * mu0 ^ 2) * Vectorize(M0)(mu0 * sqrt(x)))
  }

#########################################################

#' Probability density function for Model T3
#'
#' Value of probability density function for Model T3 of \insertCite{MAR19;textual}{MSCquartets}, Proposition 4.2.
#'
#' @references
#' \insertRef{MAR19}{MSCquartets}
#'
#' @param x  statistic value (e.g., likelihood ratio statistic, or other power divergence statistic)
#' @param mu0 parameter of density function
#' @param alpha0 parameter of density function
#' @param beta0 parameter of density function
#'
#' @return  value of density function
#'
#' @seealso \code{\link{T1density}}
#'
#' @export
T3density <- function(x,
                      mu0,
                      alpha0,
                      beta0) {
  (1 / (2 * sqrt(2 * pi * x))) * (exp(-(1 / 2) * x) * (1 - erf((1 / sqrt(
    2
  )) *
    (
      sqrt(x) * tan(beta0) - mu0
    ))) + exp(-(1 / 2) * (sqrt(x) - mu0 * cos(alpha0)) ^ 2) *
    (1 - erf((1 / sqrt(
      2
    )) * (
      sqrt(x) * tan(beta0) + mu0 * sin(alpha0)
    ))) +
    exp(-(1 / 2) * (sqrt(x) + mu0 * cos(alpha0)) ^ 2) * (1 - erf((1 / sqrt(
      2
    )) *
      (
        sqrt(x) * tan(alpha0) + mu0 * sin(alpha0)
      ))))
}

#########################################################

#' Hypothesis test for quartet counts fitting a tree under the MSC
#'
#' Test the hypothesis H_0= T1 or T3 model of \insertCite{MAR19;textual}{MSCquartets}, vs. H_1 = everything else.
#' T1 is for a specific species quartet topology, and T3 for any species quartet topology.
#'
#' @details
#' This function implements two of
#' the versions of the test given by \insertCite{MAR19;textual}{MSCquartets} as well as parametric
#' bootstrapping, with other procedures for when some expected counts are small.
#' When the topology and/or the internal
#' quartet branch length is not specified by the null hypothesis these are more accurate tests than,
#' say, a Chi-square with one degree of freedom, which is not theoretically
#' justified near the singularities and boundaries of the models.
#'
#' If \code{method="MLtest"}, this uses the test by that name described in Section 7 of \insertCite{MAR19;textual}{MSCquartets}.
#' For both the T1 and T3 models the test is slightly anticonservative over a small range of true internal edges of the quartet species tree.
#' Although the test generally performs well in practice, it lacks a uniform asymptotic guarantee over
#' the full parameter space for either T1 or T3.
#'
#' If \code{method="conservative"}, a conservative test described by \insertCite{MAR19;textual}{MSCquartets} is used. For model T3 this
#' uses the Chi-square distribution with 1 degree of freedom
#' (the "least favorable" approach), while for model T1
#' it uses the Minimum Adjusted Bonferroni, based on precomputed values from simulations with n=1e+6.
#' These conservative tests are asymptotically guaranteed to reject the null
#' hypothesis at most at a specified level, but at the expense of increased type II errors.
#'
#' If \code{method="bootstrap"}, then parametric bootstrapping is performed, based on parameter estimates of the quartet topology
#' and internal edge length. The bootstrap sample size is given by the \code{bootstrap} argument.
#'
#'
#' When some or all expected topology counts are small, the methods \code{"MLest"} and \code{"conservative"} are not appropriate.
#' The argument \code{smallsample} determines whether a precomputed bootstrap of 1e+8 samples, or actual boostrapping with the specified size,
#' is used when the total sample is small (<30).
#' The argument \code{smallcounts} determines whether bootstrapping or a faster approximate method is used when only some counts are small.
#' The approximate approach returns a precomputed p-value, found by replacing the largest observed count
#' with 1e+6 and performing 1e+8 bootstraps for the model T3. When n >30 and
#' some expected counts are small, the quartet tree error probability is small and the bootstrap p-value is
#' approximately independent of the choice of T3 or T1 and of the largest observed count.
#'
#' For model T1, the first entry of \code{obs} is treated as the count of gene quartets concordant with the species tree.
#'
#' The returned p-value should be taken with caution when there is a small sample size, e.g. less than 30 gene trees.
#' The returned value of \code{$edgelength} is a consistent estimator, but not the MLE, of the internal
#' edge length in coalescent units. Although consistent, the MLE for this length is biased.
#' Our consistent estimator is still biased, but with less bias than the MLE. See \insertCite{MAR19;textual}{MSCquartets}
#' for more discussion on dealing with the bias of parameter estimates in the
#' presence of boundaries and/or singularities of parameter spaces.
#'
#' @references
#' \insertRef{MAR19}{MSCquartets}
#'
#' @param obs  vector of 3 counts of resolved quartet frequencies
#' @param model  \code{"T1"} or \code{"T3"}, for the models of \insertCite{MAR19;textual}{MSCquartets}
#' @param lambda  parameter for power-divergence statistic (e.g., 0 for likelihood ratio statistic, 1 for Chi-squared statistic)
#' @param method \code{"MLtest"},\code{"conservative"}, or \code{"bootstrap"}
#' @param smallsample \code{"precomputed"} or \code{"bootstrap"}, method of obtaining p-value when sample is small (<30)
#' @param smallcounts \code{"precomputed"} or \code{"bootstrap"}, method of obtaining p-value when some (but not all) counts are small
#' @param bootstraps  number of samples for bootstrapping
#'
#' @return \code{output} where \code{output$p.value} is the p-value and \code{output$edgelength} is a consistent estimator of the
#' internal edge length in coalescent units, possibly \code{Inf}.
#'
#' @examples
#'  quartetTreeTest(c(17,72,11),"T3")
#'  quartetTreeTest(c(17,72,11),"T1")
#'  quartetTreeTest(c(72,11,17),"T1")
#'  quartetTreeTest(c(11,17,72),"T1")
#'
#' @seealso \code{\link{quartetTreeTestInd}}
#'
#' @importFrom stats pbinom rmultinom integrate pnorm uniroot
#'
#' @export
quartetTreeTest <- function (obs, model = "T3", lambda = 0, method = "MLest", smallsample = "precomputed",
                             smallcounts = "precomputed", bootstraps = 10^4)
{
  T1bias <- function(mu02, mu0) {
    (1/sqrt(2 * pi)) * exp(-(1/2) * mu02^2) + (1/2) * mu02 *
      (1 + erf(mu02/sqrt(2))) - mu0
  }
  integrand <- function(y, mu02, alpha02, beta02) {
    (1/sqrt(2 * pi)) * y * (exp(-(1/2) * (y - mu02)^2) *
                              erf(y * (1/tan(beta02))/sqrt(2)) + exp(-(1/2) * (y +
                                                                                 mu02 * sin(alpha02))^2) * (erf((y * (1/tan(alpha02)) +
                                                                                                                   mu02 * cos(alpha02))/sqrt(2)) + erf((y * tan(alpha02 +
                                                                                                                                                                  beta02) - mu02 * cos(alpha02))/sqrt(2))))
  }
  T3bias <- function(mu02, mu0, n) {
    phi02 <- min((1/(4 * (n + mu02^2))) * (4 * n + 3 * mu02^2 -
                                             mu02 * sqrt(8 * n + 9 * mu02^2)), 1)
    alpha02 <- atan(1/sqrt(3 * (3 - 2 * phi02)))
    beta02 <- (1/2) * (pi/2 - alpha02)
    bias <- integrate(integrand, lower = 0, upper = Inf,
                      mu02 = mu02, alpha02 = alpha02, beta02 = beta02)$value -
      mu0
    return(bias)
  }
  if (is.numeric(obs) == FALSE || is.vector(obs) == FALSE ||
      length(obs) != 3 || min(obs) < 0 || sum(obs) < 1 || all.equal(sum(obs),
                                                                    round(sum(obs))) != TRUE || model %in% c("T1", "T3") ==
      FALSE || is.numeric(lambda) == FALSE || is.vector(lambda) ==
      FALSE || length(lambda) != 1 || smallsample %in% c("bootstrap",
                                                         "precomputed") == FALSE || smallcounts %in% c("bootstrap",
                                                                                                       "precomputed") == FALSE || is.numeric(bootstraps) ==
      FALSE || is.vector(bootstraps) == FALSE || length(bootstraps) !=
      1 || bootstraps <= 0 || bootstraps%%1 != 0 || method %in%
      c("MLest", "conservative", "bootstrap") == FALSE) {
    stop("Invalid arguments: obs must be a numeric non-negative vector of length 3 summing to a positive integer;\n      model must be \"T1\" or \"T3\";\n      lambda must be a real number;\n      method must be \"MLest\", \"conservative\", or \"bootstrap\";\n      smallsample must be \"bootstrap\" or \"precomputed\";\n      smallcounts must be \"bootstrap\" or \"approximate\";\n      bootstraps must be a positive integer.")
  }
  else {
    obs <- as.numeric(obs)
    n <- round(sum(obs))
    if (model == "T3") {
      obs <- sort(obs, decreasing = TRUE)
      expd <- c(obs[1], (1/2) * (n - obs[1]), (1/2) * (n -
                                                         obs[1]))
    }
    else {
      expd <- c(max(obs[1], n/3), (1/2) * (n - max(obs[1],
                                                   n/3)), (1/2) * (n - max(obs[1], n/3)))
    }
    phi0 <- 3 * (n - expd[1])/(2 * n)
    mu0 <- max(0, sqrt(2 * n) * (1 - phi0)/sqrt(phi0 * (3 -
                                                          2 * phi0)))
    if (all.equal(min(expd), 0) != TRUE) {
      if (model == "T1") {
        if (mu0 <= 1/sqrt(2 * pi)) {
          mu0unbiased <- 0
        }
        else if (sign(T1bias(0, mu0)) == sign(T1bias(mu0,
                                                     mu0))) {
          mu0unbiased <- mu0
        }
        else {
          mu0unbiased <- uniroot(T1bias, c(0, mu0), mu0 = mu0)$root
        }
        phi0unbiased <- (1/(4 * (n + mu0unbiased^2))) *
          (4 * n + 3 * mu0unbiased^2 - mu0unbiased *
             sqrt(8 * n + 9 * mu0unbiased^2))
      }
      else {
        if (mu0 <= 3 * sqrt(3)/(2 * sqrt(2 * pi))) {
          mu0unbiased <- 0
        }
        else if (sign(T3bias(0, mu0, n)) == sign(T3bias(mu0,
                                                        mu0, n))) {
          mu0unbiased <- mu0
        }
        else {
          mu0unbiased <- uniroot(T3bias, c(0, mu0), mu0 = mu0,
                                 n = n)$root
        }
        phi0unbiased <- (1/(4 * (n + mu0unbiased^2))) *
          (4 * n + 3 * mu0unbiased^2 - mu0unbiased *
             sqrt(8 * n + 9 * mu0unbiased^2))
        if (all.equal(sum(abs(obs - expd)), 0) != TRUE) {
          alpha0unbiased <- atan(1/sqrt(3 * (3 - 2 *
                                               phi0unbiased)))
          beta0unbiased <- (1/2) * (pi/2 - alpha0unbiased)
        }
      }
    }
    else {
      phi0unbiased <- phi0
    }
    t <- -log(phi0unbiased)
    if (all.equal(sum(abs(obs - expd)), 0) == TRUE) {
      p <- 1
    }
    else {
      if (n < 30) {
        if (smallsample == "precomputed" && (all.equal(sum(abs(obs -
                                                               trunc(obs))), 0) == TRUE || all.equal(sum(abs(obs -
                                                                                                             trunc(obs) - rep(1/3, 3))), 0) == TRUE || all.equal(sum(abs(obs -
                                                                                                                                                                         trunc(obs) - rep(2/3, 3))), 0) == TRUE)) {
          if (model == "T1") {
            obs <- c(obs[1], sort(obs[2:3], decreasing = TRUE))
          }
          if (all.equal(sum(abs(obs - trunc(obs))), 0) ==
              TRUE) {
            if (model == "T1") {
              pvaluetable <- precompT1int
            }
            else if (model == "T3") {
              pvaluetable <- precompT3int
            }
          }
          else if (all.equal(sum(abs(obs - trunc(obs) -
                                     rep(1/3, 3))), 0) == TRUE) {
            if (model == "T1") {
              pvaluetable <- precompT1onethird
            }
            else if (model == "T3") {
              pvaluetable <- precompT3onethird
            }
          }
          else if (all.equal(sum(abs(obs - trunc(obs) -
                                     rep(2/3, 3))), 0) == TRUE) {
            if (model == "T1") {
              pvaluetable <- precompT1twothirds
            }
            else if (model == "T3") {
              pvaluetable <- precompT3twothirds
            }
          }
          r <- which(colSums(abs(t(pvaluetable[, 1:3]) -
                                   obs) < .Machine$double.eps^0.5) == 3)
          p <- as.numeric(pvaluetable[r, 4])
        }
        else {
          if (all.equal(sum(abs(obs - trunc(obs))), 0) ==
              TRUE || all.equal(sum(abs(obs - trunc(obs) -
                                        rep(1/3, 3))), 0) == TRUE || all.equal(sum(abs(obs -
                                                                                       trunc(obs) - rep(2/3, 3))), 0) == TRUE) {
            warning(paste0("Consider choosing smallsample=",
                           "\"precomputed\"", " for accurate p-values from 10^8 bootstraps for any n<30."))
          }
          stat <- powerDivStat(obs, expd, lambda)
          if (bootstraps == 0) {
            bootstraps <- 10^4
          }
          count <- 0
          sims <- rmultinom(bootstraps, n, prob = c(1 -
                                                      (2/3) * phi0unbiased, (1/3) * phi0unbiased,
                                                    (1/3) * phi0unbiased))
          if (model == "T3") {
            sims <- apply(sims, 2, sort, decreasing = TRUE)
          }
          for (i in 1:bootstraps) {
            if (model == "T3") {
              expected <- c(sims[1, i], (1/2) * (n -
                                                   sims[1, i]), (1/2) * (n - sims[1, i]))
            }
            else {
              expected <- c(max(sims[, i], n/3), (1/2) *
                              (n - max(sims[, i], n/3)), (1/2) * (n -
                                                                    max(sims[, i], n/3)))
            }
            simstat <- powerDivStat(sims[, i], expected,
                                    lambda)
            count <- ifelse(simstat >= stat, count +
                              1, count)
          }
          p <- count/bootstraps
        }
      }
      else {
        stat <- powerDivStat(obs, expd, lambda)
        if (method == "bootstrap" || (method != "bootstrap" &&
                                      min(expd) < 5 && smallcounts == "bootstrap") ||
            (method != "bootstrap" && min(expd) < 5 &&
             smallcounts != "bootstrap" && ((all.equal(sum(abs(obs *
                                                               3 - round(obs * 3))), 0) != TRUE) || (all.equal(sum(abs(obs *
                                                                                                                       3 - round(obs * 3))), 0) == TRUE && all.equal(var(obs - trunc(obs)), 0) != TRUE)))) {
          warning("Bootstrapping has been selected. p-values are approximate. Bootstrapping is only necessary when some expectations are small and approximate bootstrap p-values not available.")
          if (method == "bootstrap" && min(expd) < 5 &&
              smallcounts == "precomputed") {
            warning("method is prioritized over smallcounts. Bootstrap has been chosen.")
          }
          if (method != "bootstrap" && min(expd) < 5 &&
              smallcounts != "bootstrap" && ((all.equal(sum(abs(obs *
                                                                3 - round(obs * 3))), 0) != TRUE) || (all.equal(sum(abs(obs *
                                                                                                                        3 - round(obs * 3))), 0) == TRUE && all.equal(var(obs - trunc(obs)), 0) != TRUE))) {
            warning("Approximate bootstrap p-values not available when two expected counts are below 5 and observed counts are not all integers, all integers + 1/3 or all integers + 2/3.")
          }
          if (bootstraps == 0) {
            bootstraps <- 10^4
          }
          count <- 0
          sims <- rmultinom(bootstraps, n, prob = c(1 -
                                                      (2/3) * phi0unbiased, (1/3) * phi0unbiased,
                                                    (1/3) * phi0unbiased))
          if (model == "T3") {
            sims <- apply(sims, 2, sort, decreasing = TRUE)
          }
          for (i in 1:bootstraps) {
            if (model == "T3") {
              expected <- c(sims[1, i], (1/2) * (n -
                                                   sims[1, i]), (1/2) * (n - sims[1, i]))
            }
            else {
              expected <- c(max(sims[, i], n/3), (1/2) *
                              (n - max(sims[, i], n/3)), (1/2) * (n -
                                                                    max(sims[, i], n/3)))
            }
            simstat <- powerDivStat(sims[, i], expected,
                                    lambda)
            count <- ifelse(simstat >= stat, count +
                              1, count)
          }
          p <- count/bootstraps
        }
        else if (method != "bootstrap" && min(expd) <
                 5 && smallcounts == "precomputed" && (all.equal(sum(abs(obs *
                                                                         3 - round(obs * 3))), 0) == TRUE) && all.equal(var(obs - trunc(obs)), 0) == TRUE) {
          sortobs23times3 = sort(round(3 * obs[2:3]))
          if (sortobs23times3[1]%%3 == 0) {
            if (all.equal(sortobs23times3, c(0, 27)) ==
                TRUE) {
              p <- 0.000919
            }
            else if (all.equal(sortobs23times3, c(3,
                                                  24)) == TRUE) {
              p <- 0.0209
            }
            else if (all.equal(sortobs23times3, c(6,
                                                  21)) == TRUE) {
              p <- 0.103
            }
            else if (all.equal(sortobs23times3, c(9,
                                                  18)) == TRUE) {
              p <- 0.359
            }
            else if (all.equal(sortobs23times3, c(12,
                                                  15)) == TRUE) {
              p <- 0.79
            }
            else if (all.equal(sortobs23times3, c(0,
                                                  24)) == TRUE) {
              p <- 0.00191
            }
            else if (all.equal(sortobs23times3, c(3,
                                                  21)) == TRUE) {
              p <- 0.0397
            }
            else if (all.equal(sortobs23times3, c(6,
                                                  18)) == TRUE) {
              p <- 0.175
            }
            else if (all.equal(sortobs23times3, c(9,
                                                  15)) == TRUE) {
              p <- 0.523
            }
            else if (all.equal(sortobs23times3, c(0,
                                                  21)) == TRUE) {
              p <- 0.00416
            }
            else if (all.equal(sortobs23times3, c(3,
                                                  18)) == TRUE) {
              p <- 0.0763
            }
            else if (all.equal(sortobs23times3, c(6,
                                                  15)) == TRUE) {
              p <- 0.297
            }
            else if (all.equal(sortobs23times3, c(9,
                                                  12)) == TRUE) {
              p <- 0.768
            }
            else if (all.equal(sortobs23times3, c(0,
                                                  18)) == TRUE) {
              p <- 0.00871
            }
            else if (all.equal(sortobs23times3, c(3,
                                                  15)) == TRUE) {
              p <- 0.129
            }
            else if (all.equal(sortobs23times3, c(6,
                                                  12)) == TRUE) {
              p <- 0.476
            }
            else if (all.equal(sortobs23times3, c(0,
                                                  15)) == TRUE) {
              p <- 0.0183
            }
            else if (all.equal(sortobs23times3, c(3,
                                                  12)) == TRUE) {
              p <- 0.24
            }
            else if (all.equal(sortobs23times3, c(6,
                                                  9)) == TRUE) {
              p <- 0.737
            }
            else if (all.equal(sortobs23times3, c(0,
                                                  12)) == TRUE) {
              p <- 0.0393
            }
            else if (all.equal(sortobs23times3, c(3,
                                                  9)) == TRUE) {
              p <- 0.439
            }
            else if (all.equal(sortobs23times3, c(0,
                                                  9)) == TRUE) {
              p <- 0.086
            }
            else if (all.equal(sortobs23times3, c(3,
                                                  6)) == TRUE) {
              p <- 0.681
            }
            else if (all.equal(sortobs23times3, c(0,
                                                  6)) == TRUE) {
              p <- 0.197
            }
            else if (all.equal(sortobs23times3, c(0,
                                                  3)) == TRUE) {
              p <- 0.478
            }
          }
          else if (sortobs23times3[1]%%3 == 1) {
            if (all.equal(sortobs23times3, c(1, 28)) ==
                TRUE) {
              p <- 0.00222
            }
            else if (all.equal(sortobs23times3, c(4,
                                                  25)) == TRUE) {
              p <- 0.0237
            }
            else if (all.equal(sortobs23times3, c(7,
                                                  22)) == TRUE) {
              p <- 0.119
            }
            else if (all.equal(sortobs23times3, c(10,
                                                  19)) == TRUE) {
              p <- 0.358
            }
            else if (all.equal(sortobs23times3, c(13,
                                                  16)) == TRUE) {
              p <- 0.777
            }
            else if (all.equal(sortobs23times3, c(1,
                                                  25)) == TRUE) {
              p <- 0.00452
            }
            else if (all.equal(sortobs23times3, c(4,
                                                  22)) == TRUE) {
              p <- 0.0436
            }
            else if (all.equal(sortobs23times3, c(7,
                                                  19)) == TRUE) {
              p <- 0.203
            }
            else if (all.equal(sortobs23times3, c(10,
                                                  16)) == TRUE) {
              p <- 0.52
            }
            else if (all.equal(sortobs23times3, c(1,
                                                  22)) == TRUE) {
              p <- 0.0094
            }
            else if (all.equal(sortobs23times3, c(4,
                                                  19)) == TRUE) {
              p <- 0.0794
            }
            else if (all.equal(sortobs23times3, c(7,
                                                  16)) == TRUE) {
              p <- 0.295
            }
            else if (all.equal(sortobs23times3, c(10,
                                                  13)) == TRUE) {
              p <- 0.754
            }
            else if (all.equal(sortobs23times3, c(1,
                                                  19)) == TRUE) {
              p <- 0.0194
            }
            else if (all.equal(sortobs23times3, c(4,
                                                  16)) == TRUE) {
              p <- 0.143
            }
            else if (all.equal(sortobs23times3, c(7,
                                                  13)) == TRUE) {
              p <- 0.47
            }
            else if (all.equal(sortobs23times3, c(1,
                                                  16)) == TRUE) {
              p <- 0.0404
            }
            else if (all.equal(sortobs23times3, c(4,
                                                  13)) == TRUE) {
              p <- 0.233
            }
            else if (all.equal(sortobs23times3, c(7,
                                                  10)) == TRUE) {
              p <- 0.72
            }
            else if (all.equal(sortobs23times3, c(1,
                                                  13)) == TRUE) {
              p <- 0.085
            }
            else if (all.equal(sortobs23times3, c(4,
                                                  10)) == TRUE) {
              p <- 0.422
            }
            else if (all.equal(sortobs23times3, c(1,
                                                  10)) == TRUE) {
              p <- 0.113
            }
            else if (all.equal(sortobs23times3, c(4,
                                                  7)) == TRUE) {
              p <- 0.665
            }
            else if (all.equal(sortobs23times3, c(1,
                                                  7)) == TRUE) {
              p <- 0.237
            }
            else if (all.equal(sortobs23times3, c(1,
                                                  4)) == TRUE) {
              p <- 0.532
            }
          }
          else {
            if (all.equal(sortobs23times3, c(2, 26)) ==
                TRUE) {
              p <- 0.00827
            }
            else if (all.equal(sortobs23times3, c(5,
                                                  23)) == TRUE) {
              p <- 0.0436
            }
            else if (all.equal(sortobs23times3, c(8,
                                                  20)) == TRUE) {
              p <- 0.199
            }
            else if (all.equal(sortobs23times3, c(11,
                                                  17)) == TRUE) {
              p <- 0.516
            }
            else if (all.equal(sortobs23times3, c(2,
                                                  23)) == TRUE) {
              p <- 0.0164
            }
            else if (all.equal(sortobs23times3, c(5,
                                                  20)) == TRUE) {
              p <- 0.0777
            }
            else if (all.equal(sortobs23times3, c(8,
                                                  17)) == TRUE) {
              p <- 0.298
            }
            else if (all.equal(sortobs23times3, c(11,
                                                  14)) == TRUE) {
              p <- 0.739
            }
            else if (all.equal(sortobs23times3, c(2,
                                                  20)) == TRUE) {
              p <- 0.0237
            }
            else if (all.equal(sortobs23times3, c(5,
                                                  17)) == TRUE) {
              p <- 0.146
            }
            else if (all.equal(sortobs23times3, c(8,
                                                  14)) == TRUE) {
              p <- 0.466
            }
            else if (all.equal(sortobs23times3, c(2,
                                                  17)) == TRUE) {
              p <- 0.0468
            }
            else if (all.equal(sortobs23times3, c(5,
                                                  14)) == TRUE) {
              p <- 0.238
            }
            else if (all.equal(sortobs23times3, c(8,
                                                  11)) == TRUE) {
              p <- 0.703
            }
            else if (all.equal(sortobs23times3, c(2,
                                                  14)) == TRUE) {
              p <- 0.0926
            }
            else if (all.equal(sortobs23times3, c(5,
                                                  11)) == TRUE) {
              p <- 0.409
            }
            else if (all.equal(sortobs23times3, c(2,
                                                  11)) == TRUE) {
              p <- 0.185
            }
            else if (all.equal(sortobs23times3, c(5,
                                                  8)) == TRUE) {
              p <- 0.645
            }
            else if (all.equal(sortobs23times3, c(2,
                                                  8)) == TRUE) {
              p <- 0.377
            }
            else if (all.equal(sortobs23times3, c(2,
                                                  5)) == TRUE) {
              p <- 0.525
            }
          }
        }
        else if (method != "bootstrap" && min(expd) >=
                 5) {
          if (method == "MLest") {
            if (model == "T3") {
              p <- ifelse(stat == 0, 1, min(1, max(0,
                                                   integrate(T3density, lower = stat, upper = Inf,
                                                             mu0 = mu0unbiased, alpha0 = alpha0unbiased,
                                                             beta0 = beta0unbiased)$value)))
            }
            else {
              p <- ifelse(stat == 0, 1, min(1, max(0,
                                                   integrate(T1density, lower = stat, upper = Inf,
                                                             mu0 = mu0unbiased)$value)))
            }
          }
          else if (method == "conservative") {
            if (model == "T3") {
              p <- pchisq(stat, 1, lower.tail = FALSE)
            }
            else {
              if (mu0 >= 0 && mu0 <= 0.125) {
                beta <- 0.5
                q <- pT1beta50
              }
              else if (mu0 > 0.125 && mu0 <= 1.75) {
                beta <- 0.45
                q <- pT1beta45
              }
              else if (mu0 > 1.75 && mu0 <= 2.75) {
                beta <- 0.3
                q <- pT1beta30
              }
              else if (mu0 > 2.75 && mu0 <= 3.75) {
                beta <- 0.2
                q <- pT1beta20
              }
              else {
                beta <- 0.05
                q <- pT1beta5
              }
              mu0lower <- max(0, mu0 + sqrt(2) * erf.inv(2 *
                                                           beta - 1))
              p1 <- ifelse(stat == 0, 1, min(1, max(0,
                                                    integrate(T1density, lower = stat, upper = Inf,
                                                              mu0 = mu0lower)$value)))
              p <- (length(q[q <= p1]))/length(q)
            }
          }
        }
      }
    }
    output <- list(p, t)
    names(output) <- c("p.value", "edgelength")
    return(output)
  }
}




#########################################################

#' Multiple independent hypothesis tests for quartet counts fitting a species tree under the MSC
#'
#' Perform a tree hypothesis test for all quartet counts in an input table, as if the counts for different choices of 4 taxa
#' are independent.
#'
#' @details This function assumes all quartets are resolved.  The test performed and the arguments
#' are described more fully in \code{quartetTreeTest}.
#'
#' @references
#' \insertRef{MAR19}{MSCquartets}
#'
#' @param rqt  table of resolved quartet counts, as produced by \code{quartetTableResolved}, or \code{quartetStarTestInd}
#' @param model  \code{"T1"} for a specific species tree topology, or \code{"T3"} for any species tree topology, with these
#' models explained more fully by \insertCite{MAR19;textual}{MSCquartets}
#' @param lambda power divergence statistic parameter (e.g., 0 for likelihood ratio statistic, 1 for Chi-squared statistic)
#' @param method \code{"MLest"}, \code{"conservative"}, or \code{"bootstrap"}; see \code{quartetTreeTest} for explanation
#' @param smallsample \code{"precomputed"} or \code{"bootstrap"}, method of obtaining p-value when sample is small (<30)
#' @param smallcounts \code{"precomputed"} or \code{"bootstrap"}, method of obtaining p-value when some counts are small, so
#' the chosen \code{method} is inappropriate
#' @param bootstraps  number of samples for bootstrapping
#' @param speciestree  species tree, in Newick as text, to determine quartet for T1 test; required for \code{model="T1"},
#' ignored for \code{model="T3"}
#' @return
#'   if \code{model="T3"}, a copy of \code{rqt} with a new column \code{"p_T3"} appended with p-values for each quartet;
#'   if \code{model="T1"}, a copy of \code{rqt} with 2 columns appended: \code{"p_T1"} with p-values, and \code{"qindex"}
#'   giving index of quartet consistent with specified species tree,
#'   i.e., 1 if 12|34 on species tree, 2 if 13|24, 3 if 14|23
#'
#' @seealso \code{\link{quartetTreeTest}}, \code{\link{quartetTestPlot}}, \code{\link{quartetStarTestInd}}, \code{\link{quartetTableResolved}}
#'
#' @examples
#' gtrees=read.tree(file=system.file("extdata","dataGeneTreeSample",package="MSCquartets"))
#' tnames=c("t1","t2","t3","t4","t5","t6")
#' QT=quartetTable(gtrees,tnames)
#' RQT=quartetTableResolved(QT)
#' stree="(((t5,t6),t4),((t1,t2),t3));"
#' pTable3=quartetTreeTestInd(RQT,"T3")
#' quartetTablePrint(pTable3[1:6,])
#' stree="((((t5,t6),t4),t7),((t8,t9),((t1,t2),t3)));"
#' pTable1=quartetTreeTestInd(RQT,"T1",speciestree=stree)
#' quartetTablePrint(pTable1[1:6,])
#'
#' @importFrom ape unroot read.tree cophenetic.phylo
#' @export
quartetTreeTestInd <- function (rqt,
                                model = "T3",
                                lambda = 0,
                                method = "MLest",
                                smallsample = "precomputed",
                                smallcounts = "precomputed",
                                bootstraps = 10 ^ 4,
                                speciestree = NULL)
{
  if (model %in% c("T1", "T3") ==
      FALSE || is.numeric(lambda) == FALSE || is.vector(lambda) ==
      FALSE ||
      length(lambda) != 1 || smallsample %in% c("bootstrap", "precomputed") == FALSE ||
      smallcounts %in% c("bootstrap", "precomputed") == FALSE ||
      is.numeric(bootstraps) ==
      FALSE ||
      is.vector(bootstraps) == FALSE || length(bootstraps) !=
      1 || bootstraps <= 0 || bootstraps %% 1 != 0 || method %in%
      c("MLest", "conservative", "bootstrap") == FALSE) {
    stop(
      "Invalid arguments: lambda must be real; smallcounts must be \"bootstrap\" or \"approximate\";\n      bootstraps must be a positive integer; method must be \"MLest\",\"conservative\", or \"bootstrap\"."
    )
  }
  colnames = colnames(rqt)
  taxanames = setdiff(
    colnames,
    c(
      "12|34",
      "13|24",
      "14|23",
      "1234",
      "p_cut",
      "cutindex",
      "p_star",
      "p_T3",
      "p_T1"
    )
  )
  M = dim(rqt)[1]
  n = length(taxanames)
  if (model == "T3") {
    pTable = cbind(rqt, p_T3 = 0)
    message("Applying hypothesis test for model T3 to ", M, " quartets.")
    for (m in 1:M) {
      obs = unname(rqt[m, c("12|34", "13|24", "14|23")])
      pvec = quartetTreeTest(
        obs,
        model,
        lambda = lambda,
        smallsample = smallsample,
        smallcounts = smallcounts,
        bootstraps = bootstraps,
        method = method
      )
      pTable[m, "p_T3"] = pvec$p.value
    }
  }
  else {
    if (model != "T1") {
      stop("Invalid model: must use 'T1' or 'T3'")
    }
    else {
      if (is.null(speciestree)) {
        stop("Species tree topology must be supplied for model T1")
      }
      else {
        stree = unroot(read.tree(text = speciestree))
        nedges = dim(stree$edge)[1]
        stree$edge.length = rep(1, nedges)
        D = cophenetic.phylo(stree)
        D = D[order(rownames(D)), order(colnames(D))]
        pTable = cbind(rqt, p_T1 = 0, qindex = 0)
        message("Applying hypothesis test for model T1 to ", M, " quartets.")
        for (m in 1:M) {
          qnames = which(rqt[m, 1:n] == 1)
          a = D[qnames[1], qnames[2]] + D[qnames[3], qnames[4]]
          b = D[qnames[1], qnames[3]] + D[qnames[2], qnames[4]]
          c = D[qnames[1], qnames[4]] + D[qnames[2], qnames[3]]
          squartet = which.min(c(a, b, c))
          qcounts = rqt[m, c("12|34", "13|24", "14|23")]
          temp = qcounts[1]
          qcounts[1] = qcounts[squartet]
          qcounts[squartet] = temp
          pTable[m, "qindex"] = squartet
          pvec = quartetTreeTest(
            unname(qcounts),
            model,
            lambda = lambda,
            smallcounts = smallcounts,
            bootstraps = bootstraps
          )
          pTable[m, "p_T1"] = pvec$p.value
        }
      }
    }
  }
  return(pTable)
}


#############################################################

#' Produce simplex plot with  results of quartet hypothesis tests
#'
#' Plot a 2-d probability simplex, with points for all quartet count vectors. Colors
#' indicate rejection or failure to reject for tests at specified levels.
#'
#' @details The first argument of this function is a table of quartets and p-values. The
#' plot may show results of either the T1, T3, or 2-cut
#' test, with or without a star tree test (depending on whether a \code{"p_star"} column is in the table and/or \code{beta =1}).
#' The p-values must be computed by previous calls to
#' \code{quartetTreeTestInd} (for \code{"T1"} or \code{"T3"} p-values)
#' and \code{quartetStarTestInd} (for \code{"star"} p-values). The \code{NANUQ} and \code{NANUQdist}
#' functions include calls to these tree test functions.
#'
#' @param pTable table of quartets and p-values, as produced by \code{quartetTreeTestInd},
#' \code{quartetStarTestInd}, or \code{NANUQ}
#' @param test  model to use, for tree null hypothesis; options are \code{"T1"}, \code{"T3"}, \code{"cut"}, \code{"NANUQ"}
#' @param alpha level for tree test with null hypothesis given by \code{test}
#' @param beta level for test with null hypothesis star tree;
#' test results plotted only if \code{beta<1} and \code{"p_star"} column present in \code{pTable}
#' @param cex scaling factor for size of plotted symbols
#' @return NULL
#'
#' @seealso \code{\link{quartetTreeTestInd}}, \code{\link{quartetStarTestInd}},
#' \code{\link{NANUQ}}, \code{\link{NANUQdist}}
#'
#' @examples
#' gtrees=read.tree(file=system.file("extdata","dataGeneTreeSample",package="MSCquartets"))
#' tnames=c("t1","t2","t3","t4","t5","t6")
#' QT=quartetTable(gtrees,tnames[1:6])
#' RQT=quartetTableResolved(QT)
#' stree="(((t5,t6),t4),((t1,t2),t3));"
#' pTable=quartetTreeTestInd(RQT,"T1",speciestree=stree)
#' pTable=quartetStarTestInd(pTable)
#' quartetTestPlot(pTable, "T1", alpha=.05, beta=.95)
#'
#' @export
quartetTestPlot <- function(pTable,
                            test,
                            alpha = .05,
                            beta = 1,
                            cex = 1) {
  if (!(is.numeric(alpha) && is.numeric(beta))) {
    stop("Test levels alpha and beta must be numeric.")
  }

  lineWidth = 1.5
  yellowColor = "goldenrod3" #pick specific colors
  orangeColor = "darkgreen"
  redColor = "red"
  blueColor = "blue"

  counts = c("12|34", "13|24", "14|23")
  M = dim(pTable)[1]

  if ((test == "T1") & (beta == 1)) {
    #  T1 only
    p_tree.small = (pTable[, "p_T1"] < alpha)
    p_star.small = rep(TRUE, M)
    titletext = bquote("T1 Model," ~ alpha * "=" * .(alpha))
    legtext = c("reject tree", "fail to reject tree")
    legpch = c(2, 1)
    legcol = c(redColor, blueColor)
    model = "T1"
    #  for model 1 need to adjust point to reflect species tree
    for (m in 1:M) {
      qind = pTable[m, "qindex"]
      if (qind == 3) {
        pTable[m, counts] = pTable[m, counts[c(3, 1, 2)]]
      } else {
        if (qind == 2) {
          pTable[m, counts] = pTable[m, counts[c(2, 3, 1)]]
        }
      }
    }
  } else {
    if ((test == "T1") & (beta < 1)) {
      #  T1  and star
      p_tree.small = (pTable[, "p_T1"] < alpha)
      p_star.small = (pTable[, "p_star"] <= beta)
      titletext = bquote("T1 Model," ~ alpha * "=" * .(alpha) * "," ~ beta *
                           "=" * .(beta))
      legtext = c(
        "reject tree & star",
        "fail to reject tree/reject star",
        "fail to reject tree & star",
        "reject tree/fail to reject star"
      )
      legpch = c(2, 1, 0, 4)
      legcol = c(redColor, blueColor, yellowColor, orangeColor)
      model = "T1"
      #  for model 1 need to adjust point to reflect species tree
      for (m in 1:M) {
        qind = pTable[m, "qindex"]
        if (qind == 3) {
          pTable[m, counts] = pTable[m, counts[c(3, 1, 2)]]
        } else {
          if (qind == 2) {
            pTable[m, counts] = pTable[m, counts[c(2, 3, 1)]]
          }
        }
      }
    } else{
      if ((test == "T3") & (beta == 1)) {
        #T3 only
        p_tree.small = (pTable[, "p_T3"] < alpha)
        p_star.small = rep(TRUE, M)
        titletext = bquote("T3 Model," ~ alpha * "=" * .(alpha))
        legtext = c("reject tree", "fail to reject tree")
        legpch = c(2, 1)
        legcol = c(redColor, blueColor)
        model = "T3"
      } else {
        if ((test == "NANUQ") | ((test == "T3") & (beta < 1))) {
          #T3 and star
          p_tree.small = (pTable[, "p_T3"] < alpha)
          p_star.small = (pTable[, "p_star"] <= beta)

          if (test == "NANUQ") {
            titletext = bquote("NANUQ," ~ alpha * "=" * .(alpha) * "," ~ beta *
                                 "=" * .(beta))
          } else {
            titletext = bquote("T3 Model," ~ alpha * "=" * .(alpha) * "," ~ beta *
                                 "=" * .(beta))
          }
          legtext = c(
            "reject tree & star",
            "fail to reject tree/reject star",
            "fail to reject tree & star",
            "reject tree/fail to reject star"
          )
          legpch = c(2, 1, 0, 4)
          legcol = c(redColor, blueColor, yellowColor, orangeColor)
          model = "T3"
        } else {
          if ((test == "cut") & (beta == 1)) {
            #  cut only
            p_tree.small = (pTable[, "p_cut"] < alpha)
            p_star.small = rep(TRUE, M)
            titletext = bquote("2-cut Model," ~ alpha * "=" * .(alpha))
            legtext = c("reject cut", "fail to reject cut")
            legpch = c(2, 1)
            legcol = c(redColor, blueColor)
            model = "cut"
          } else {
            if ((test == "cut") & (beta < 1)) {
              #  cut  and star
              p_tree.small = (pTable[, "p_cut"] < alpha)
              p_star.small = (pTable[, "p_star"] <= beta)
              titletext = bquote("2-cut Model," ~ alpha * "=" * .(alpha) * "," ~ beta *
                                   "=" * .(beta))
              legtext = c(
                "reject cut & star",
                "fail to reject cut/reject star",
                "fail to reject cut & star",
                "reject cut/fail to reject star"
              )
              legpch = c(2, 1, 0, 4)
              legcol = c(redColor, blueColor, yellowColor, orangeColor)
              model = "cut"
            } else
            stop("Invalid test name")
          }
        }
      }
    }
  }


  simplexPrepare(model, maintitle = "Quartet Hypothesis Test", titletext =
                   titletext) # draw outline of simplex and lines for model

  blues = (!p_tree.small) & (p_star.small)
  reds = (p_tree.small) & (p_star.small)
  yellows = (!p_tree.small) & (!p_star.small)
  oranges = (p_tree.small) & (!p_star.small)

  bluepoints = pTable[which(blues), counts, drop = FALSE] # points to plot blue
  nblue = dim(bluepoints)[1]
  redpoints = pTable[which(reds), counts, drop = FALSE] #      red
  nred = dim(redpoints)[1]
  yellowpoints = pTable[which(yellows), counts, drop = FALSE] #      yellow
  nyellow = dim(yellowpoints)[1]
  orangepoints = pTable[which(oranges), counts, drop = FALSE] #        orange
  norange = dim(orangepoints)[1]

  if (nblue > 0) {
    for (i in 1:nblue) {
      simplexPoint(
        bluepoints[i,],
        type = "o",
        pch = 1,
        col = blueColor,
        lwd = lineWidth,
        cex = cex
      )
    }
  }
  if (nred > 0) {
    for (i in 1:dim(redpoints)[1]) {
      simplexPoint(
        redpoints[i,],
        type = "o",
        pch = 2,
        col = redColor,
        lwd = lineWidth,
        cex = cex
      )
    }
  }
  if (nyellow > 0) {
    for (i in 1:dim(yellowpoints)[1]) {
      simplexPoint(
        yellowpoints[i,],
        type = "o",
        pch = 0,
        col = yellowColor,
        lwd = lineWidth,
        cex = cex
      )
    }
  }
  if (norange > 0) {
    message("Some points (green) not rejected as star, but rejected as tree.")
    for (i in 1:dim(orangepoints)[1]) {
      simplexPoint(
        orangepoints[i,],
        type = "o",
        pch = 4,
        col = orangeColor,
        lwd = lineWidth,
        cex = cex
      )
    }
  }
  legend(
    'topleft',
    legend = legtext,
    pch = legpch,
    col = legcol,
    cex = 0.8,
    inset = 0.01,
    pt.cex = 1.2,
    lty = 0,
    lwd = 2
  )
}


#########################################

#' Bayesian posterior probability of error in 4-taxon unrooted species tree topology estimate
#'
#' From a gene quartet count concordance factor (qcCF), computes Bayesian posterior probabilities
#' of the three 4-taxon species tree topologies and the Bayesian posterior probability that the
#' assumed topology is incorrect, under the assumption that the counts arise from the MSC on some species tree.
#'
#' @details The Jeffreys prior is used for internal branch length, along with the uniform prior
#' on the resolved topology.
#'
#' @references
#' \insertRef{MAR19}{MSCquartets}
#'
#'@param obs vector of counts for 3 topologies
#'@param model \code{"T3"} or \code{"T1"}, for the models of \insertCite{MAR19;textual}{MSCquartets} describing an unspecified species
#'tree topology (\code{"T3"}), or the topology whose count is the first entry of \code{obs} (\code{"T1"})
#'
#'@return \code{(error.prob, top.probs)} where \code{error.prob} is the species tree error probability
#'and \code{top.probs} is a vector of the three species tree topology probabilities in the order of \code{obs};
#'for model \code{"T1"} the species tree used is the one
#'corresponding to the first count; for model \code{"T3"} the species
#'tree is the one corresponding to the largest count
#'
#'
#'@examples
#' obs <- c(28,32,30)
#' quartetTreeErrorProb(obs,model="T1")
#' quartetTreeErrorProb(obs,model="T3")
#'
#'@importFrom zipfR Ibeta
#'
#'@export
quartetTreeErrorProb <- function(obs,
                                 model = "T3") {
  if (is.vector(obs) == FALSE ||
      length(obs) != 3 ||
      min(obs) < 0 ||
      model %in% c("T1", "T3") == FALSE) {
    stop(
      'Invalid arguments: obs must be an non-negative vector of length 3, model must be "T1" or "T3".'
    )
  } else {
    # If all observations are equal, then all tree topologies have equal probabilities and quartet tree
    # error probability is 2/3.
    if (all.equal(obs[1], obs[2]) == TRUE &&
        all.equal(obs[1], obs[3]) == TRUE &&
        all.equal(obs[2], obs[3]) == TRUE) {
      treeerrorprob <- 2 / 3
      probs <- c(1 / 3, 1 / 3, 1 / 3)
    } else {
      # Compute posterior probabilities for each topology. Prior probabilities for each tree topology
      # are 1/3 and the Jeffreys prior is used for phi0.
      logp1 <-
        obs[1] * log(2) + zipfR::Ibeta(2 / 3, obs[2] + obs[3] + 1 / 2, obs[1] +
                                         1 / 2, log = TRUE)
      logp2 <-
        obs[2] * log(2) + zipfR::Ibeta(2 / 3, obs[1] + obs[3] + 1 / 2, obs[2] +
                                         1 / 2, log = TRUE)
      logp3 <-
        obs[3] * log(2) + zipfR::Ibeta(2 / 3, obs[1] + obs[2] + 1 / 2, obs[3] +
                                         1 / 2, log = TRUE)

      b1 <- exp(logp2 - logp1)
      b2 <- exp(logp3 - logp1)

      p1 <- 1 / (1 + b1 + b2)
      p2 <- b1 / (1 + b1 + b2)
      p3 <- b2 / (1 + b1 + b2)

      probs <- c(p1, p2, p3)

      if (model == "T1") {
        treeerrorprob <- sum(probs[2:3])
      } else {
        treeerrorprob <- sum(probs[setdiff(1:3, which.max(obs))])
      }
    }

    output <- list(treeerrorprob, probs)
    names(output) <-
      c("tree.error.probability", "tree.probabilities")
    return(output)
  }
}

#############################################################################

#' Apply Holm-Bonferroni method to adjust for multiple tests
#'
#' Apply the Holm-Bonferroni method to
#' adjust for multiple hypothesis tests performed on quartets from a data set of gene trees.
#'
#' @param pTable a table of quartets with p-values, as computed by
#' \code{quartetTreeTestInd} or \code{quartetStarTestInd}
#' @param model one of \code{"T1"}, \code{"T3"}, or \code{"star"}, where \code{pTable} contains a column \code{p_model} of p-values
#' @param alpha test level, for rejection of adjusted p-values less than or equal to \code{alpha}
#' @return the same table, with rows reordered, and 2 new columns of 1) adjusted p-values,
#' and 2) "Y" or "N" for indicating "reject" or "fail to reject"
#'
#' @seealso \code{\link{quartetTreeTestInd}}, \code{\link{quartetStarTestInd}}
#'
#'
#' @details When p-values are computed for each quartet using
#' \code{quartetTreeTestInd} or \code{quartetStarTestInd},
#' multiple comparisons are being done for one dataset. The
#' Holm-Bonferroni method \insertCite{HolmBonf79}{MSCquartets} adjusts these p-values upward,
#' controlling the familywise error rate. The probability
#' of at least one false discovery (rejection of the null hypothesis)
#' is no more than the significance level.
#'
#' The Holm-Bonferroni method does not require that test hypotheses are independent, which
#' is important for its application to quartet counts presumed to have arisen on a single
#' species tree.
#'
#' This can be a low power test (often failing to reject when the null hypothesis is false).
#' In particular for the T1 and T3 tests, it does not consider the relationships between edge
#' lengths for different sets of four taxa.
#'
#' Warning: Output of this function should not be used as input for other
#' MSCquartets functions, as it reorders the quartets in the table.
#'
#' @references \insertRef{HolmBonf79}{MSCquartets}
#'
#' @examples
#' gtrees=read.tree(file=system.file("extdata","dataGeneTreeSample",package="MSCquartets"))
#' taxanames=taxonNames(gtrees)
#' QT=quartetTable(gtrees,taxanames[1:6])
#' RQT=quartetTableResolved(QT)
#' pTable=quartetTreeTestInd(RQT,"T3")
#' pTable[1:10,]
#' HBpTable=HolmBonferroni(pTable,"T3",.05)
#' HBpTable[1:10,]
#'
#' @export
HolmBonferroni <- function(pTable,
                           model,
                           alpha = .05) {
  if (!(model %in% c("T1", "T3", "star")))  {
    stop('Argument model must be one of "T1", "T3", or "star".')
  }
  columnname = paste0("p_", model)
  if (!(columnname %in% colnames(pTable))) {
    stop(c('Argument pTable has no column ', columnname, '.'))
  }
  m = nrow(pTable)
  pTable <-
    pTable[order(pTable[, columnname]),] #sort rows of table by columname

  multipliers <- m:1 # descending integers
  HBps <-
    pmin(1, pTable[, columnname] * multipliers) #adjust p-values
  FirstNonFail = min(c(which(HBps > alpha), m + 1))
  Reject = c(rep("Y", FirstNonFail - 1), rep("N", m - FirstNonFail + 1))
  newcols = cbind(HBps, Reject)
  pname = paste0("HB", columnname)
  colnames(newcols) = c(pname, paste0("Reject(", pname, ")"))
  pTable = cbind.data.frame(pTable, newcols)
  message("Holm-Bonferroni method leads to rejection for ",
          (FirstNonFail - 1),
          " tests.")
  return(pTable)
}
