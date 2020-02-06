#' easy filter
#'
#' easier filter than the convolution filter from Silvia Dee. The function is the filter_fuction.R
#'
#' @param dt time step length in years 
#' @param d18O d18O input that comes back filtered
#' @param tau0 tau0 is the mean transit time through the karst in years
#' @param temp temperatur that fits to the data --> later important for butterworth lowpass filter
#'
#' @return list with filtered-d18O and filter function h
#' @export
#'
#' @examples
easy_sensor_wm3 <- function(dt, d18O, tau0, temp){
  
  #source('~/201903_FirstFilter/functions_used.R')
  tau0 <- tau0/dt
  
  h_s=1E-4						# ensure kernel approaches 0
  tau_s=-tau0*log(tau0*h_s)
  tau=seq(from = 1.0, to = ceiling(tau_s), by =1.0)
  
  timeseries_len <- length(d18O)
  if (length(tau)>= (0.5*timeseries_len)){
    writeLines("Warning: mean residence time (tau0) is too close to half of the length of the timeseries.")
  }
  
  h <- (1/tau0)*exp(-tau/tau0) #kernel length as long as time serie
  hint <- DescTools::AUC(seq(0:(length(h)-1)), h)   # obtain normalization constant using Simpson's rule
  h <- h/hint# normalize transit time distribution
  
  d18OK <- filter_function3(d18O,h, tau0,dt)
  
  avg_period <- 10.0
  Tm <- mean(temp)  # extract mean
  # lowpass filter the  temperature
  temp <-c(temp, rep.int(Tm, length(d18OK)-length(temp)))
  Tl <- bwf_filter(temp-Tm,1.0/avg_period,3) + Tm
  # apply thermal fractionation (Eq 11 in [4]) WHYYY???
  d18Oc <- (d18OK + 1000)/1.03086*exp(2780/Tl**2 - 0.00289)-1000.0
  
  Results_SpeleoFilter <- list("d18OK" = d18OK, "h" = h, "d18Oc" = d18Oc)
  
  return (Results_SpeleoFilter)
}
