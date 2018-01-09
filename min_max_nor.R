min_max_nor <- function(m){
  m_size <- dim(m)
  new_m = matrix(0,m_size[1],m_size[2])
  max_m = max(m)
  min_m = min(m)
  for (i in (1:length(m))){
    new_m[i] = (m[i]-min_m)/(max_m-min_m)
  }
  return(new_m)
}

