#' @title Example data frame for release kinetics of an element from soil
#' @description User is advised to prepare the data as suggested in the example to study the release kinetics of an element from soil using different equations, viz. power function, simplified Elovich, parabolic diffusion, zero order, first order, and second order.
#' @usage df_release
#' @format Write the following notations on the spreadsheet:
#' Kt	for Cumulative amount of element released (mass or mole of element per unit mass of soil, e.g., mg/kg or mmol/kg) after time ‘t’
#' Time	The time intervals used in the kinetics study (user can choose the unit as per convenience)
#' @details The rate of release of any nutrient or pollutant element from soil solids to soil solution is very important for availability to plants (Barber 1984), with possible beneficial or harmful effects. The rate of release of any plant nutrient or pollutant element in soil can be assessed by release kinetics studies. In these studies, soils are equilibrated with a given solution for different time intervals, and the amount of the element in question is measured in the solution after each time interval. The cumulative amounts of the element released after different time intervals are plotted as the dependent variable against the time intervals as the independent variable. For every equation, a coefficient of determination (R-squared) and a standard error of estimate will be obtained. In addition, all equations except the second order, two constants, namely ‘a’ and ‘b’ will be obtained; whereas, the maximum desorbable amount of the element and the half the time required to release the maximum amount are obtained from the second order equation.
#'
#' @export
df_release = read.table(text = "Time	Kt
1	34.4
7	57.8
19	81
31	102
43	121
67	140
91	157
121	176
163	193
223	210
307	228
403	246", header = TRUE)

#' @title Fitting of release kinetics data to the power function equation
#'
#' @description The function fits the cumulative amount of an element released from soil with a given extractant over a certain time period to the power function equation. It generates the coefficient of determination (R-squared) and standard error of estimate to show the goodness of fit. The parameters of power function equation, i.e., ‘a’ and ‘b’ can be obtained from this function.
#' @importFrom graphics abline legend curve
#' @importFrom stats coef nls lm predict var sigma
#'
#' @usage rk_power(Kt = Kt, Time = Time,...)
#'
#' @param Kt	Cumulative amount of element released (mass or mole of element per unit mass of soil, e.g., mg/kg or mmol/kg) after time ‘t’
#' @param Time	The time intervals used in the kinetics study (user can choose the unit as per convenience)
#' @param ... Any other argument that can be passed to base plot
#' @return
#' \itemize{
#'   \item R2: Coefficient of determination (more the value, better the fit)
#'   \item SEE: Standard error of estimate (less the value, better the fit)
#'   \item a: Equation constant
#'   \item b: Equation constant (overall release rate coefficient)
#' }
#'
#' @references Havlin, J.L., Westfall, D.G., Olsen, S.R., 1985. Mathematical models for potassium release kinetics in calcareous soil. Soil Science Society of America Journal 49, 371–376. https://doi.org/10.2136/sssaj1985.03615995004900020020x
#' @details The power function equation is expressed as Kt = a×t^b (Havlin et al. 1985), where ‘Kt’ is the cumulative amount of the element released over time ‘t’, and ‘a’ and ‘b’ are constants. The constant ‘b’ can be considered as the overall release rate coefficient. Though power function does not suggest any particular release mechanism, its parameters can be successfully used to show the effect of any particular treatment on nutrient or pollutant release behaviour in a given soil, or comparing different soils for their ability to release a particular nutrient or pollutant under a given condition.
#' @export
#' @examples
#' with(data = df_release, rk_power(Kt = Kt, Time = Time))

rk_power <- function(Kt = Kt, Time = Time,...){
  
  fit1 <- lm(log(Kt) ~ log(Time))
  
  #Power function equation
  m <- nls(Kt ~ a*(Time^b),
           start = list(a = exp(coef(fit1)[1]), b = coef(fit1)[2])) # power formula: y = a*x^b
  summary(m)
  
  myplot <- function(x, y, xlab=deparse(substitute(x)), ylab=deparse(substitute(y)),...){
    plot(x, y, xlab=xlab, ylab=ylab, ...)
  }
  myplot(x = Time, y = Kt, ...)
  curve(predict(m, newdata = data.frame(Time = x)), add = TRUE, col = "red")
  
  #Get the R2 value
  R2 <- format(var(predict(m))/var(Kt), digits=3)
  ##Printing of the equation and R2 on the plot
  ##Round the coefficients for better output
  cf <- round(coef(m), 3)
  
  ##Make the regression equation
  eq <- as.expression(bquote(y == .(cf[1]) * x^.(cf[2])))

##Printing of the equation and R2 on the plot
legend("topleft", legend = c(eq, as.expression(bquote(R^2 == .(R2)))), bty = "n")

##Printing of the equation and R2 on the plot
#legend("topleft", legend = c(eq, as.expression(bquote(R^2 == .(Rsq)))), bty = "n")
a <- coef(m)[[1]]
b <- coef(m)[[2]]
SEE <- sigma(m)
return(list(a = a, b = b, R2 = R2, SEE = SEE))
}


#' @title Fitting of release kinetics data to the simplified Elovich equation
#'
#' @description The function fits the cumulative amount of an element released from soil with a given extractant over a certain time period to the simplified Elovich equation. It generates the coefficient of determination (R-squared) and standard error of estimate to show the goodness of fit. The parameters of simplified Elovich equation, i.e., ‘a’ and ‘b’ can be obtained from this function.
#' @importFrom graphics abline legend
#' @importFrom stats coef lm
#'
#' @usage rk_Elovich(Kt = Kt, Time = Time,...)
#'
#' @param Kt	Cumulative amount of element released (mass or mole of element per unit mass of soil, e.g., mg/kg or mmol/kg) after time ‘t’
#' @param Time	The time intervals used in the kinetics study (user can choose the unit as per convenience)
#' @param ... Any other argument that can be passed to base plot
#' @return
#' \itemize{
#'   \item R2: Coefficient of determination (more the value, better the fit)
#'   \item SEE: Standard error of estimate (less the value, better the fit)
#'   \item a: Equation constant
#'   \item b: Equation constant (overall release rate coefficient)
#' }
#'
#' @references Havlin, J.L., Westfall, D.G., Olsen, S.R., 1985. Mathematical models for potassium release kinetics in calcareous soil. Soil Science Society of America Journal 49, 371–376. https://doi.org/10.2136/sssaj1985.03615995004900020020x
#' @details The simplified Elovich equation, expressed as Kt= a + b × ln(t) (Havlin et al. 1985), where ‘Kt’ is the cumulative amount of the element released over time ‘t’, and ‘a’ and ‘b’ are constants. The constant ‘b’ can be considered as the overall release rate coefficient.
#' @export
#' @examples
#' with(data = df_release, rk_Elovich(Kt = Kt, Time = Time))

rk_Elovich <- function(Kt = Kt, Time = Time,...){

  #Simplified Elovich
  fit2 <- lm(Kt ~ log(Time))
  summary(fit2)

  #R2 calculation
  R2 <- var(predict(fit2))/var(Kt)

  myplot <- function(x, y, xlab=deparse(substitute(x)), ylab=deparse(substitute(y)),...){
    plot(x, y, xlab=xlab, ylab=ylab, ...)
  }
  myplot(x = log(Time), y = Kt, ...)
  abline(lm(Kt ~ log(Time)), col = "blue")

  ##Printing of the equation and R2 on the plot
  ##Round the coefficients for better output
  cf <- round(coef(fit2), 3)

  ##Make the regression equation
  eq <- paste0("y = ", cf[1],
               ifelse(sign(cf[2])==1, " + ", " - "), abs(cf[2]), " x ")

  #Get the R2 value
  R2 = format(summary(fit2)$r.squared,digits=3)

  ##Printing of the equation and R2 on the plot
  legend("topleft", legend = c(eq, as.expression(bquote(R^2 == .(R2)))), bty = "n")

  a <- coef(fit2)[[1]]
  b <- coef(fit2)[[2]]
  SEE <- sigma(fit2)

  return(list(a = a, b = b, R2 = R2, SEE = SEE))
}

#' @title Fitting of release kinetics data to the parabolic diffusion equation
#'
#' @description The function fits the cumulative amount of an element released from soil with a given extractant over a certain time period to the parabolic diffusion equation. It generates the coefficient of determination (R-squared) and standard error of estimate to show the goodness of fit. The parameters of parabolic diffusion equation, i.e., ‘a’ and ‘b’ can be obtained from this function.
#' @importFrom graphics abline legend
#' @importFrom stats coef lm
#'
#' @usage rk_pd(Kt = Kt, Time = Time,...)
#'
#' @param Kt	Cumulative amount of element released (mass or mole of element per unit mass of soil, e.g., mg/kg or mmol/kg) after time ‘t’
#' @param Time	The time intervals used in the kinetics study (user can choose the unit as per convenience)
#' @param ... Any other argument that can be passed to base plot
#' @return
#' \itemize{
#'   \item R2: Coefficient of determination (more the value, better the fit)
#'   \item SEE: Standard error of estimate (less the value, better the fit)
#'   \item a: Equation constant
#'   \item b: Equation constant (overall release rate coefficient)
#' }
#'
#' @references Havlin, J.L., Westfall, D.G., Olsen, S.R., 1985. Mathematical models for potassium release kinetics in calcareous soil. Soil Science Society of America Journal 49, 371–376. https://doi.org/10.2136/sssaj1985.03615995004900020020x
#' @details The parabolic diffusion equation is expressed as Kt = a + b × √t (Havlin et al. 1985), where ‘Kt’ is the cumulative amount of the element released over time ‘t’, and ‘a’ and ‘b’ are constants. The constants ‘a’ and ‘b’ can be considered as the initial release rate constant and the overall release rate constant, respectively. Though empirical in nature, a good fit to parabolic equation indicates that diffusion is the rate-limiting step in the desorption process.
#' @export
#' @examples
#' with(data = df_release, rk_pd(Kt = Kt, Time = Time))

rk_pd <- function(Kt = Kt, Time = Time,...){

  #Parabolic diffusion
  time0.5 <- (Time)^0.5
  fit3 <- lm(Kt ~ time0.5)
  summary(fit3)

  #R2 calculation
  R2 <- var(predict(fit3))/var(Kt)

  myplot <- function(x, y, xlab=deparse(substitute(x)), ylab=deparse(substitute(y)),...){
    plot(x, y, xlab=xlab, ylab=ylab, ...)
  }
  myplot(x = time0.5, y = Kt, ...)
  abline(lm(Kt ~ time0.5), col = "blue")

  ##Printing of the equation and R2 on the plot
  ##Round the coefficients for better output
  cf <- round(coef(fit3), 3)

  ##Make the regression equation
  eq <- paste0("y = ", cf[1],
               ifelse(sign(cf[2])==1, " + ", " - "), abs(cf[2]), " x ")

  #Get the R2 value
  R2 = format(summary(fit3)$r.squared,digits=3)

  ##Printing of the equation and R2 on the plot
  legend("topleft", legend = c(eq, as.expression(bquote(R^2 == .(R2)))), bty = "n")

  a <- coef(fit3)[[1]]
  b <- coef(fit3)[[2]]
  SEE <- sigma(fit3)

  return(list(a = a, b = b, R2 = R2, SEE = SEE))
}


#' @title Fitting of release kinetics data to the zero-order equation
#'
#' @description The function fits the cumulative amount of an element released from soil with a given extractant over a certain time period to the zero-order equation. It generates the coefficient of determination (R-squared) and standard error of estimate to show the goodness of fit. The parameters of zero-order equation, i.e., ‘a’ and ‘b’ can be obtained from this function.
#' @importFrom graphics abline legend
#' @importFrom stats coef lm
#'
#' @usage rk_zero(Kt = Kt, Time = Time,...)
#'
#' @param Kt	Cumulative amount of element released (mass or mole of element per unit mass of soil, e.g., mg/kg or mmol/kg) after time ‘t’
#' @param Time	The time intervals used in the kinetics study (user can choose the unit as per convenience)
#' @param ... Any other argument that can be passed to base plot
#' @return
#' \itemize{
#'   \item R2: Coefficient of determination (more the value, better the fit)
#'   \item SEE: Standard error of estimate (less the value, better the fit)
#'   \item a: Equation constant
#'   \item b: Equation constant (overall release rate coefficient)
#' }
#'
#' @references Martin, H.W., Sparks, D.L., 1983. Kinetics of nonexchangeable potassium release from two coastal plain soils. Soil Science Society of America Journal 47, 883–887. https://doi.org/10.2136/sssaj1983.03615995004700050008x
#' @details The zero-order equation is expressed as (Km – Kt ) = a – b × t (Martin and Sparks 1983), where ‘Km’ is the total amount of the element released over the entire study period, ‘Kt’ is the cumulative amount of the element released over time ‘t’, and ‘a’ and ‘b’ are constants. The constant ‘b’ can be considered as the overall release rate coefficient.
#' @export
#' @examples
#' with(data = df_release, rk_zero(Kt = Kt, Time = Time, ylab = "Kmax - Kt"))

rk_zero <- function(Kt = Kt, Time = Time,...){

  #Zero order
  Kmax_Kt <- max(Kt, na.rm = T) - Kt
  fit4 <- lm(Kmax_Kt ~ Time)
  summary(fit4)

  #R2 calculation
  #R2 <- format(var(predict(fit4))/var(Kt),digits=3)

  myplot <- function(x, y, xlab=deparse(substitute(x)), ylab=deparse(substitute(y)),...){
    plot(x, y, xlab=xlab, ylab=ylab, ...)
  }
  myplot(x = Time, y = Kmax_Kt, ...)
  abline(lm(Kmax_Kt ~ Time), col = "blue")

  ##Printing of the equation and R2 on the plot
  ##Round the coefficients for better output
  cf <- round(coef(fit4), 3)

  ##Make the regression equation
  eq <- paste0("y = ", cf[1],
               ifelse(sign(cf[2])==1, " + ", " - "), abs(cf[2]), " x ")

  #Get the R2 value
  R2 = format(summary(fit4)$r.squared,digits=3)

  ##Printing of the equation and R2 on the plot
  legend("topright", legend = c(eq, as.expression(bquote(R^2 == .(R2)))), bty = "n")

  a <- coef(fit4)[[1]]
  b <- coef(fit4)[[2]]
  SEE <- sigma(fit4)

  return(list(a = a, b = b, R2 = R2, SEE = SEE))
}


#' @title Fitting of release kinetics data to the first-order equation
#'
#' @description The function fits the cumulative amount of an element released from soil with a given extractant over a certain time period to the first-order equation. It generates the coefficient of determination (R-squared) and standard error of estimate to show the goodness of fit. The parameters of first-order equation, i.e., ‘a’ and ‘b’ can be obtained from this function.
#' @importFrom graphics abline legend
#' @importFrom stats coef lm
#'
#' @usage rk_first(Kt = Kt, Time = Time,...)
#'
#' @param Kt	Cumulative amount of element released (mass or mole of element per unit mass of soil, e.g., mg/kg or mmol/kg) after time ‘t’
#' @param Time	The time intervals used in the kinetics study (user can choose the unit as per convenience)
#' @param ... Any other argument that can be passed to base plot
#' @return
#' \itemize{
#'   \item R2: Coefficient of determination (more the value, better the fit)
#'   \item SEE: Standard error of estimate (less the value, better the fit)
#'   \item a: Equation constant
#'   \item b: Equation constant (overall release rate coefficient)
#' }
#'
#' @references Martin, H.W., Sparks, D.L., 1983. Kinetics of nonexchangeable potassium release from two coastal plain soils. Soil Science Society of America Journal 47, 883–887. https://doi.org/10.2136/sssaj1983.03615995004700050008x
#' @details The first-order equation is expressed as ln(Km – Kt ) = a – b × t (Martin and Sparks 1983), where ‘Km’ is the total amount of the element released over the entire study period, ‘Kt’ is the cumulative amount of the element released over time ‘t’, and ‘a’ and ‘b’ are constants. The constant ‘b’ can be considered as the overall release rate coefficient.  The first-order equation assumes that the amount of the element present in the desorbable form on soil solids or exchange sites determines its release rate.
#' @export
#' @examples
#' with(data = df_release, rk_first(Kt = Kt, Time = Time, ylab = "ln(Kmax - Kt)", xlab = "Time"))

rk_first <- function(Kt = Kt, Time = Time,...){

  #First order equation
  Kmax_Kt <- max(Kt, na.rm = T) - Kt
  df1 <- cbind.data.frame(ln_Kmax_Kt = log(Kmax_Kt), time = Time)

  #Replace NaN & Inf with NA
  df1[is.na(df1) | df1=="-Inf" | df1=="Inf"] = NA

  fit5 <- lm(ln_Kmax_Kt ~ time, data = df1)
  summary(fit5)

  #R2 calculation
  #R2 <- format(var(predict(fit4))/var(Kt),digits=3)

  myplot <- function(x, y, xlab=deparse(substitute(x)), ylab=deparse(substitute(y)),...){
    plot(x, y, xlab=xlab, ylab=ylab, ...)
  }
  myplot(x = df1$time, y = df1$ln_Kmax_Kt, ...)
  abline(lm(ln_Kmax_Kt ~ time, data = df1), col = "blue")

  ##Printing of the equation and R2 on the plot
  ##Round the coefficients for better output
  cf <- round(coef(fit5), 3)

  ##Make the regression equation
  eq <- paste0("y = ", cf[1],
               ifelse(sign(cf[2])==1, " + ", " - "), abs(cf[2]), " x ")

  #Get the R2 value
  R2 = format(summary(fit5)$r.squared,digits=3)

  ##Printing of the equation and R2 on the plot
  legend("topright", legend = c(eq, as.expression(bquote(R^2 == .(R2)))), bty = "n")

  a <- coef(fit5)[[1]]
  b <- coef(fit5)[[2]]
  SEE <- sigma(fit5)

  return(list(a = a, b = b, R2 = R2, SEE = SEE))
}


#' @title Fitting of release kinetics data to the second-order equation
#'
#' @description The function fits the cumulative amount of an element released from soil with a given extractant over a certain time period to the second-order equation. It generates the coefficient of determination (R-squared) and standard error of estimate to show the goodness of fit. The parameters of second-order equation, i.e., ‘Kmax’ and ‘t1/2’ can be obtained from this function.
#' @importFrom graphics abline legend
#' @importFrom stats coef lm
#'
#' @usage rk_second(Kt = Kt, Time = Time,...)
#'
#' @param Kt	Cumulative amount of element released (mass or mole of element per unit mass of soil, e.g., mg/kg or mmol/kg) after time ‘t’
#' @param Time	The time intervals used in the kinetics study (user can choose the unit as per convenience)
#' @param ... Any other argument that can be passed to base plot
#' @return
#' \itemize{
#'   \item R2: Coefficient of determination (more the value, better the fit)
#'   \item SEE: Standard error of estimate (less the value, better the fit)
#'   \item a: Equation constant
#'   \item b: Equation constant (overall release rate coefficient)
#' }
#'
#' @references
#' Grimme, H. 1980. The effect of field strength on quantity of K desorbed from soils by electro-ultrafiltration. Zeitschrift für Pflanzenernährung und Bodenkunde 143, 98–106. https://doi.org/10.1002/jpln.19801430113
#'
#' Lü, X., Xu, J., Ma, W., Lu, Y., 2007. Comparison of seven kinetic equations for k release and application of kinetic parameters. Pedosphere 17, 124–129. https://doi.org/10.1016/S1002-0160(07)60017-4
#'
#' @details The second-order equation is expressed as t ÷ Kt = (t ÷ Kmax) + (t1/2 ÷ Kmax) (Grimme 1980; Lü et al. 2007), where ‘Kmax’ is the maximum desorbable amount of the element, ‘Kt’ is the cumulative amount of the element released over time ‘t’, and ‘t1/2’ is the half the time required to release the maximum amount.
#' @export
#' @examples
#' with(data = df_release, rk_second(Kt = Kt, Time = Time, ylab = "t/Kt"))

rk_second <- function(Kt = Kt, Time = Time,...){

  #Second order equation
  t_kt <- Time/Kt
  fit6 <- lm(t_kt ~ Time)
  summary(fit6)

  #R2 calculation
  #R2 <- format(var(predict(fit4))/var(Kt),digits=3)

  myplot <- function(x, y, xlab=deparse(substitute(x)), ylab=deparse(substitute(y)),...){
    plot(x, y, xlab=xlab, ylab=ylab, ...)
  }
  myplot(x = Time, y = t_kt, ...)
  abline(lm(t_kt ~ Time), col = "blue")

  ##Printing of the equation and R2 on the plot
  ##Round the coefficients for better output
  cf <- round(coef(fit6), 3)

  ##Make the regression equation
  eq <- paste0("y = ", cf[1],
               ifelse(sign(cf[2])==1, " + ", " - "), abs(cf[2]), " x ")

  #Get the R2 value
  R2 = format(summary(fit6)$r.squared,digits=3)

  ##Printing of the equation and R2 on the plot
  legend("topleft", legend = c(eq, as.expression(bquote(R^2 == .(R2)))), bty = "n")

  a <- coef(fit6)[[1]]
  b <- coef(fit6)[[2]]
  SEE <- sigma(fit6)

  return(list(a = a, b = b, R2 = R2, SEE = SEE))
}
