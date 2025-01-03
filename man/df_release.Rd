% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Release_kinetics.R
\docType{data}
\name{df_release}
\alias{df_release}
\title{Example data frame for release kinetics of an element from soil}
\format{
Write the following notations on the spreadsheet:
Kt	for Cumulative amount of element released (mass or mole of element per unit mass of soil, e.g., mg/kg or mmol/kg) after time ‘t’
Time	The time intervals used in the kinetics study (user can choose the unit as per convenience)
}
\usage{
df_release
}
\description{
User is advised to prepare the data as suggested in the example to study the release kinetics of an element from soil using different equations, viz. power function, simplified Elovich, parabolic diffusion, zero order, first order, and second order.
}
\details{
The rate of release of any nutrient or pollutant element from soil solids to soil solution is very important for availability to plants (Barber 1984), with possible beneficial or harmful effects. The rate of release of any plant nutrient or pollutant element in soil can be assessed by release kinetics studies. In these studies, soils are equilibrated with a given solution for different time intervals, and the amount of the element in question is measured in the solution after each time interval. The cumulative amounts of the element released after different time intervals are plotted as the dependent variable against the time intervals as the independent variable. For every equation, a coefficient of determination (R-squared) and a standard error of estimate will be obtained. In addition, all equations except the second order, two constants, namely ‘a’ and ‘b’ will be obtained; whereas, the maximum desorbable amount of the element and the half the time required to release the maximum amount are obtained from the second order equation.
}
\keyword{datasets}
