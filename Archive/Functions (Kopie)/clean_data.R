clear_data <- function(data, type){
  
  max_value=array(c(1e7, 45, 1e-2, 2))
  min_value=array(c(0, -20, 0, -15))
  
  index_id = c()

  for (ii in 1:(length(data)-1)){
    if(is.na(data[ii])){
      index_id=c(index_id, ii)
    } else if(data[ii] > max_value[type]){
      index_id=c(index_id, ii)
    } else if (data[ii] < min_value[type]){
      index_id=c(index_id, ii)
    }
  }
  
  gapfiller = mean(data[-1*index_id])
  
  if(length(index_id)==0){
    return(data)
    stop()
  }
  
  for (ii in index_id){
    data[ii] = gapfiller
  }
  return(data)
}

  
  