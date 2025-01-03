% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Release_kinetics.R
\name{rk_Elovich}
\alias{rk_Elovich}
\title{Fitting of release kinetics data to the simplified Elovich equation}
\usage{
rk_Elovich(Kt = Kt, Time = Time,...)
}
\arguments{
\item{Kt}{Cumulative amount of element released (mass or mole of element per unit mass of soil, e.g., mg/kg or mmol/kg) after time ‘t’}

\item{Time}{The time intervals used in the kinetics study (user can choose the unit as per convenience)}

\item{...}{Any other argument that can be passed to base plot}
}
\value{
\itemize{
\item R2: Coefficient of determination (more the value, better the fit)
\item SEE: Standard error of estimate (less the value, better the fit)
\item a: Equation constant
\item b: Equation constant (overall release rate coefficient)
}
}
\description{
The function fits the cumulative amount of an element released from soil with a given extractant over a certain time period to the simplified Elovich equation. It generates the coefficient of determination (R-squared) and standard error of estimate to show the goodness of fit. The parameters of simplified Elovich equation, i.e., ‘a’ and ‘b’ can be obtained from this function.
}
\details{
The simplified Elovich equation, expressed as Kt= a + b × ln(t) (Havlin et al. 1985), where ‘Kt’ is the cumulative amount of the element released over time ‘t’, and ‘a’ and ‘b’ are constants. The constant ‘b’ can be considered as the overall release rate coefficient.
}
\examples{
with(data = df_release, rk_Elovich(Kt = Kt, Time = Time))
}
\references{
Havlin, J.L., Westfall, D.G., Olsen, S.R., 1985. Mathematical models for potassium release kinetics in calcareous soil. Soil Science Society of America Journal 49, 371–376. https://doi.org/10.2136/sssaj1985.03615995004900020020x
}
