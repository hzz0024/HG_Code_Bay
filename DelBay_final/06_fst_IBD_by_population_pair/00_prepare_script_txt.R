x1 <- c("HC_18", "ARN_18", "COH_18", "SR_18", "NB_18")
x2 <- c("HC_19", "ARN_19", "COH_19", "SR_19", "NB_19") 
x3 <- c("HC_21", "ARN_21", "COH_21", "SR_21", "NB_21") #c("HC_21", "ARN_21", "BEN_21", "BS_21", "COH_21", "SR_21", "NAN_21", "NB_21")
x4 <- c("HC_21spat", "ARN_21spat", "NAN_21spat")
x5 <- c("HC_18", "HC_19", "HC_21", "HC_21spat")
x6 <- c("ARN_18", "ARN_19", "ARN_21", "ARN_21spat")
x7 <- c("NAN_21", "NAN_21spat")
x8 <- c("Sur_19", "Ref_19")
x9 <- c("Sur_20", "Ref_20")
x10 <- c("A_Sur20", "B_Sur20", "Ref_20", "Ref_19")

pop_vector <- function (x) {
  n <- 2
  comb <- t(combn(x,n))
  my_vec <- c()
  for(i in 1:dim(comb)[1]) {              
    my_vec <- c(my_vec, paste("pop", comb[i,][1], comb[i,][2], sep="_"))    # Appending new value to vector
  }
  print(noquote(my_vec))
} 

pop_vector(x1)
pop_vector(x2)
pop_vector(x3)
pop_vector(x4)
pop_vector(x5)
pop_vector(x6)
pop_vector(x7)
pop_vector(x8)
pop_vector(x9)
pop_vector(x10)
