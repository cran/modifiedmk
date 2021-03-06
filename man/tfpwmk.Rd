% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tfpwmk.R
\name{tfpwmk}
\alias{tfpwmk}
\title{Mann-Kendall Trend Test Applied to Trend-Free Prewhitened Time Series Data in Presence of Serial Correlation Using Yue et al. (2002) Approach}
\usage{
tfpwmk(x)
}
\arguments{
\item{x}{- Time series data vector}
}
\value{
Z-Value  - Z statistic after trend-free prewhitening (TFPW)

Sen's Slope  - Sen's slope for TFPW series

Old Sen's Slope  - Sen's slope for original data series (x)

P-value  - P-value after trend-free prewhitening

S  - Mann-Kendall S statistic

Var(s) - Variance of S

Tau  - Mann-Kendall's Tau
}
\description{
When the time series data are not random and influenced by autocorrelation, the trend component is removed from the data and is prewhitened prior to the application of the trend test.
}
\details{
The linear trend component is removed from the original data and then prewhitened using the lag-1 serial correlation coefficient. The prewhitening data are then tested with Mann-Kendall trend test.
}
\examples{
x<-c(Nile)
tfpwmk(x)

}
\references{
Kendall, M. (1975). Rank Correlation Methods. Griffin, London, 202 pp.

Kulkarni, A. and H. von Storch. 1995. Monte carlo experiments on the effects of serial correlation on the MannKendall test of trends. Meteorologische Zeitschrift N.F, 4(2): 82-85.

Mann, H. B. (1945). Nonparametric Tests Against Trend. Econometrica, 13(3): 245-259.

Salas, J.D. (1980). Applied modeling of hydrologic times series. Water Resources Publication, 484 pp.

Sen, P. K. (1968). Estimates of the Regression Coefficient Based on Kendall’s Tau. Journal of the American Statistical Association, 63(324): 1379. <doi:10.2307/2285891>

von Storch, V. H. (1995). Misuses of statistical analysis in climate research, In: Analysis of Climate Variability: Applications of Statistical Techniques, ed. von H. V. Storch and A. Navarra A. Springer-Verlag, Berlin: 11-26.

Yue, S., Pilon, P., Phinney, B., and Cavadias, G. (2002). The influence of autocorrelation on the ability to detect trend in hydrological series. Hydrological Processes, 16(9): 1807–1829. <doi:10.1002/hyp.1095>
}
