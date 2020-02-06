MeanSpectrum2 <- function (specList, iRemoveLowest = 1, weights = rep(1, length(specList))) 
{
  print("Average spectra")
  if (length(weights) != length(specList)) 
    stop("specList and weights have a different number of elements")
  weights <- weights/sum(weights)
  specList <- lapply(specList, remove.lowestFreq, iRemove = iRemoveLowest)
  freqRef <- seq(from = min(unlist(lapply(specList, get.fstart.existing))), 
                 to = max(unlist(lapply(specList, get.fend.existing))), 
                 by = min(unlist(lapply(specList, get.df))))
  specList.interpolated <- list()
  print("I am here!")
  for (i in 1:length(specList)) specList.interpolated[[i]] <- SpecInterpolate(freqRef, 
                                                                              specList[[i]])
  for (i in 1:length(specList)) specList.interpolated[[i]]$spec <- specList.interpolated[[i]]$spec * 
    weights[i]
  NSpectra <- length(specList.interpolated)
  result <- list(freq = specList.interpolated[[1]]$freq, spec = rep(0, 
                                                                    length(specList.interpolated[[1]]$spec)))
  specMatrix <- matrix(NA, NSpectra, length(specList.interpolated[[1]]$spec))
  dofMatrix <- matrix(NA, NSpectra, length(specList.interpolated[[1]]$spec))
  print("I am even here!")
  for (i in 1:length(specList.interpolated)) {
    if (sum((result$freq - specList.interpolated[[i]]$freq)^2) > 
        0.1) 
      stop("Different spectra length or resolutions")
    specMatrix[i, ] <- specList.interpolated[[i]]$spec
    dofMatrix[i, ] <- specList.interpolated[[i]]$dof
  }
  result$spec <- colSums(specMatrix, na.rm = TRUE)
  nRecord = colSums(!is.na(specMatrix))
  result$dof <- colSums(dofMatrix, na.rm = TRUE)
  class(result) <- "spec"
  return(list(spec = AddConfInterval(result), nRecord = nRecord))
}