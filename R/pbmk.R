#' @title Decorrelated bootstrapped-based Mann-Kendall Test
#'
#' @description Lancombe et al. (2012) suggest the use of effective prewhitening (Hamed, 2009) to initially decorrelate the data after which the bootstrap based Mann-Kendall test is applied (Yue and Pilon, 2004).
#'
#' @importFrom boot tsboot
#'
#' @usage pbmk(x, ci=0.95, nsim=1000, pw="Hamed")
#'
#' @param  x  - Time series data vector
#'
#' @param  ci - Confidence Interval
#'
#' @param  nsim - Number of Simulations
#'
#' @param pw -  Unbiased pre-whitening suggested by Hamed 2009
#'
#' @return  Z  - Original Mann- Kendall Z-statistic
#'
#' Slp  - Original sen's slope
#'
#' S.orig  - Original Mann-Kendall 'S'- statistic
#'
#' Tau  - Original Mann-Kendall's Tau
#'
#' Zpw  - Bias Corrected Prewhitened  Mann- Kendall Z-statistic
#'
#' Slppw  - Bias Corrected Prewhitened  sen's slope
#'
#' Spw  - Bias Corrected Prewhitened  Mann-Kendall 'S'- statistic
#'
#' Taupw  - Bias Corrected Prewhitened  Mann-Kendall's Tau
#'
#' pval - Bootstrapped P-Value
#'
#' @references Hamed, K. H. (2009). Enhancing the effectiveness of prewhitening in trend analysis of hydrologic data. Journal of Hydrology, 368: 143-155.
#'
#' @references Kendall, M. (1975). Multivariate analysis. Charles Griffin. Londres. 0-85264-234-2.
#'
#' @references Kundzewicz, Z. W. and Robson, A. J. (2004). Change detection in hydrological records - a review of the methodology. Hydrological Sciences Journal, 49(1): 7-19.
#'
#' @references Lancombe, G., McCartney, M., and Forkuor, G. (2012). Drying climate in Ghana over the period 1960-2005: evidence from the resampling-based Mann-Kendall test at local and regional levels. Hydrological Sciences Journal, 57(8): 1594-1609. doi:10.1080/02626667.2012.728291
#'
#' @references Mann, H. B. (1945). Nonparametric Tests Against Trend. Econometrica, 13(3), 245?259. <doi:10.1017/CBO9781107415324.004>
#
#' @references van Giersbergen, N. P. A. (2005). On the effect of deterministic terms on the bias in stable AR models. Economic Letters, 89: 75-82.
#'
#' @references Yue, S. and Pilon, P. (2004). A comparison of the power of the t test, Mann-Kendall and bootstrap tests for trend detection, Hydrological Sciences Journal, 49(1): 21-37.
#'
#' @details The block bootstrap is used along with the non-parametric Mann-Kendall trend test.  A test statistic falling in the tails of the simulated empirical distribution, the results is likely significant.
#'
#' @examples x<-c(Nile[1:10])
#' pbmk(x)
#'
#' @export
#'
pbmk <- function(x, ci=0.95, nsim=1000, pw="Hamed") {
  # Initialize the test Parameters

  # Time-Series Vector
  x = x
  # Confidance Interval
  ci = ci
  #Number of Simulations
  nsim=nsim
  # Mann-Kendall Tau
  Tau = NULL
  # Modified Z-Statistic after Pre-Whitening
  Z = NULL
  # Modified P-value after Pre-Whitening
  pval = NULL
  # Initialize Mann-Kendall 'S'- Statistic - prewhitened
  S = NULL
  # Initialize Mann-Kendall 'S'- Statistic
  S.orig = NULL
  # Sen's slope estimate
  slp = NULL

  # To test whether the data is in vector format

  if (is.vector(x) == FALSE) {
    stop("Input data must be a vector")
  }

  nx<-length(x)

  #Specify minimum input vector length
  if (nx < 4) {
    stop("Input vector must contain at least four values")
  }

  # To test whether the data values are finite numbers and attempting to eliminate non-finite numbers
  if (any(is.finite(x) == FALSE)) {
    x[-c(which(is.finite(x) == FALSE))] -> x
    warning("The input vector contains non-finite numbers. An attempt was made to remove them")
  }

  if (is.null(pw) == FALSE) {
    #Calculate the lag 1 autocorrelation coefficient and the intercept
    zx<-cbind(head(x,n=nx-1),matrix(data=1, nrow=(nx-1),ncol=1),tail(seq(1:nx),n=(nx-1)))
    y<-tail(x,n=nx-1)
    zTrans<-t(zx)
    zTransz<-zTrans%*%zx
    zTranszInv<-solve(zTransz)
    zTranszInvzTrans<-zTranszInv%*%zTrans
    params<-zTranszInvzTrans%*%y
    ACFlag1<-params[1]

    #Correct for bias in the lag-1 acf using eq. 24 of Hamed (2009)
    ACFlag1BC<-((nx*ACFlag1)+2)/(nx-4)

    # Calculating pre-whitened Series

    a=1:(nx-1)
    b=2:nx
    xn<-(x[b]-(x[a]*ACFlag1BC))

    #Bootstrapped using Mann Kendall
    MK.orig <- mkttest(x)
    Z <- round(MK.orig["Z-Value"], digits = 7)
    slp <- round(MK.orig["Sen's slope"], digits = 7)
    Tau <- round(MK.orig["Tau"], digits = 7)
    S.orig <- MK.orig["S"]
    MKpw <- mkttest(xn)
    Zpw <- round(MKpw["Z-Value"], digits = 7)
    slpPW <- round(MKpw["Sen's slope"], digits = 7)
    TauPW <- round(MKpw["Tau"], digits = 7)
    Spw <- MKpw["S"]
    MKS <- function(xn) mkttest(xn)[["S"]]
    boot.out.MKS <- tsboot(xn, MKS, R=nsim, l=1, sim="fixed")
    loc <- suppressWarnings(max(which(sort(boot.out.MKS$t) < Spw)))
    if (loc == -Inf) {
      loc <- 1
    }
    pval <- loc/nsim

    cat(paste("Original Z-Value = ", Z,
              "Original Sen's Slope = ", slp,
              "Original S = ", S.orig,
              "Original Kendall's Tau = ", Tau,
              "Bias Corrected Prewhitened Z-Value = ", Zpw,
              "Bias Corrected Prewhitened Sen's Slope = ", slpPW,
              "Bias Corrected Prewhitened S = ", Spw,
              "Bias Corrected Prewhitened Kendall's Tau = ", TauPW,
              "Bootstrapped P-Value =", pval ,sep="\n"))
  } else {
    #Bootstrapped using Mann Kendall
    MK.orig <- mkttest(x)
    Z <- round(MK.orig["Z-Value"], digits = 7)
    slp <- round(MK.orig["Sen's slope"], digits = 7)
    Tau <- round(MK.orig["Tau"], digits = 7)
    S.orig <- MK.orig["S"]
    MKS1 <- function(x) mkttest(x)[["S"]]
    boot.out.MKS1 <- tsboot(x, MKS1, R=nsim, l=1, sim="fixed")
    loc <- suppressWarnings(max(which(sort(boot.out.MKS1$t) < S.orig)))

    if (loc == -Inf) {
      loc <- 1
    }
    pval <- loc/nsim

    cat(paste("Z-Value = ", Z,
              "Sen's Slope = ", slp,
              "S = ", S.orig,
              "Kendall's Tau = ", Tau,
              "Bootstrapped P-Value =", pval ,sep="\n"))
  }
}



