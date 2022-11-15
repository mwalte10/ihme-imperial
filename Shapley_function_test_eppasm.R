set.seed(1)
library(data.table)

input <- list(
              ##input into child model
              mat_prev = list(rep(NA, 53 * 10)),
              ## To be split into timing of infections later using existing proportions
              mtct_rate = rep(0.05, 1000),
              cotrim = rnorm(1000, mean = 0.33, 0.1),
              numberonart_scalar =  rnorm(1000, mean = 0.96, sd = 0.082))

#eta
eta <- function(
    int = list(mat_prev, mtct_rate, cotrim, numberonart_scalar)
                ){
  y = eppasm_mod(int$mat_prev, int$mtct_rate, int$cotrim, int$numberonart_scalar)
  return(y)
}


#x_all
Xall <- function(){
  row = sample(1000, 1)
  return(list(mat_prev = input$mat_prev[row], 
              mtct_rate = input$mtct_rate[row], 
              cotrim = input$cotrim[row],
              numberonart_scalar = input$numberonart_scalar[row]
              ))
}

Xset <- function(a, b, c){
  test <- input[sample(1000,1),a]
  names(test) = names(input)[a]
  test = data.frame(test)
  names(test) <- names(input)[a]
  return(test)
}

k <- 4
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




