##### netglm ######
####  processing ###


makeSenderMatrix <- function(x, var) {
  out <- matrix(x[,var],length(x[,var]),length(x[,var]))
  rownames(out)<- colnames(out) <- x$ID
  return(out)
}
makeReceiverMatrix <- function(x, var) {
  out <- matrix(x[,var],length(x[,var]),length(x[,var]), byrow = T)
  rownames(out)<- colnames(out) <- x$ID
  return(out)
}
makeSimMatrix <- function(x, var){
  out <- abs(outer(x[,var],x[,var], "-"))*-1
  rownames(out)<- colnames(out) <- x$ID
  return(out)
}
makeSameMatrix <- function(x, var){
  out <- (outer(x[,var],x[,var], "=="))*1
  rownames(out)<- colnames(out) <- x$ID
  return(out)
}