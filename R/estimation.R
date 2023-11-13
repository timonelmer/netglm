##### netglm ######
#### estimation ###
stars <- function(p.value){
  ifelse(p.value > 0.90 | p.value < 0.10,
         ifelse(p.value > 0.95 | p.value < 0.05, 
                ifelse(p.value > 0.99 | p.value < 0.01, 
                       ifelse(p.value > 0.999 | p.value < 0.001, "***", "**"), "*"), "x"), "")
}




