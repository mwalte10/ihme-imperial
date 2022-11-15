set.seed(1)
library(data.table)

#eta
a_obs <- rnorm(100, mean = 2, sd = 0.5)
b_obs <- rnorm(100, mean = -5, sd = 0.1)
c_obs <- rnorm(100, mean = -5, sd = 1)

input <- data.frame(a = a_obs, b = b_obs, c = c_obs)
eta <- function(
    int = list(a = 2, b = -5, c = -5)
                ){
  y = int$a^2 + int$b^2 + int$c^2
  return(y)
}


#x_all
Xall <- function(){
  row = sample(100, 1)
  return(list(a = input[row,1], b = input[row,2], c = input[row,3]))
}

Xset <- function(a, b, c){
  test <- input[sample(100,1),a]
  names(test) = colnames(input)[a]
  test = data.frame(test)
  colnames(test) <- colnames(input)[a]
  return(test)
}

k <- 3
Nv <- 100
Ni <- 100
No <- 100



# Initialize Shapley value for all players
Sh <- rep(0, k)
Sh2 <- rep(0, k)

# Estimate Var[Y] 
Y<-NULL

for (i in 1:Nv){
  X<-Xall()
  # Evaluate the objective function
  val = eta(X)
  Y <- c(Y, val)
}


EY = mean(Y)
VarY = var(Y)

# Generate all k! permutations
perms = gtools::permutations(k, k, 1:k)

# Estimate Shapley effects
m = nrow(perms)


for (p in 1:m){
  pi <- perms[p,]
  prevC = 0
  for (j in 1:k){
    if (j == k){    
      Chat = VarY
      del <- Chat - prevC
    } else {
      cVar <- NULL
      Sj = pi[c(1:j)] # set of the 1st-jth elements in pi 
      Sjc = pi[-c(1:j)] # set of the (j+1)th-kth elements in pi
      for (l in 1:No){
        Y<-NULL
        xjc <- Xset(Sjc, NULL, NULL) # sampled values of the inputs in Sjc
        for (h in 1:Ni) {
          # sample values of inputs in Sj conditional on xjc
          xj <- Xset(Sj, Sjc, xjc) 
          x <- c(xj, xjc)
          pi_s <- sort(pi,index.return=TRUE)$ix
          val <- eta(x[pi_s])
          Y <- c(Y,val)
        }
        cVar <- c(cVar,var(Y))
      }
      
      Chat = mean(cVar)
      del <- Chat - prevC
    }
    
    Sh[pi[j]] = Sh[pi[j]] + del
    Sh2[pi[j]] = Sh2[pi[j]] + del^2
    
    prevC = Chat
  }
}

Sh = Sh/m
Sh2 = Sh2/m

list(Shapley = Sh, SEShapley = sqrt((Sh2 - Sh^2)/m), VarY = VarY, EY = EY)

##Sum of Sh equals to the total variance (VarY)




