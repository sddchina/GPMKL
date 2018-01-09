z_score <- function(in_m){
  m_size <- dim(in_m)
  new_m = matrix(0,m_size[1],m_size[2])
  for (i in 1:m_size[2]){
    x = in_m[,i]
    new_m[,i] = (x - mean(x) ) /sd(x)
  }
  return(new_m)
}
