library("foreign")
library(ggplot2)
library(Ryacas) # need for limit of a function
library(sandwich) #Loading gmm required package: sandwich 
library(gmm) # for using GMM method
library(plyr) # The plyr library has a function round_any that is pretty generic to do all kinds of rounding

# a function for Generating a lag/lead variables
shift<-function(x,shift_by){
  stopifnot(is.numeric(shift_by))
  stopifnot(is.numeric(x))
  
  if (length(shift_by)>1)
    return(sapply(shift_by,shift, x=x))
  
  out<-NULL
  abs_shift_by=abs(shift_by)
  if (shift_by > 0 )
    out<-c(tail(x,-abs_shift_by),rep(NA,abs_shift_by))
  else if (shift_by < 0 )
    out<-c(rep(NA,abs_shift_by), head(x,-abs_shift_by))
  else
    out<-x
  out
}
#********************************************************************************

data_marriage <- read.dta("5-1-LabourForce_Married_NK(wifeOnly)(Stata12)(MuEst).dta") 

# I need consumption & leisure of single indivs in the 2ng stage of algorithm
data_single   <- read.dta("3_LabourForce_Singles(Stata12).dta") 

#********************************************************************************
num_type <- 4
#alpha_j <- 0.8
#alpha_i <- 0.9
#********************************************************************************

# number of married people from each type (as initial guess of the measure of both single females and single males)
m <- vector()
f <- vector()
for (i in 1:num_type) { 
  m[i] <- sum(data_marriage$type_spouse==i )
  f[i] <- sum(data_marriage$type ==i )
}
########################################################################################
#                                      Algorithm
########################################################################################
sigma_theta <- 1

# Step 1: 
#Provide an initial guess of the measure of both single females and single males
#........
TauD_guess <- vector()
for (i in 1:num_type) {
  TauD_guess[i] <- m[i]
}
Ln_TauD_guess <- log(TauD_guess)

TauS_guess <- vector()
for (j in 1:num_type) {
  TauS_guess[j] <- f[j]
}
Ln_TauS_guess <- log(TauS_guess)

# step 2 and 3 
# until convergance
#........
consumption_M_male <- function(alpha,BP,Y,w_i,w_j)
{alpha * (1-BP) *(Y+w_i+w_j)}
consumption_M_female <- function(alpha,BP,Y,w_i,w_j)
{alpha * BP *(Y+w_i+w_j)}
Utility <- function(alpha,consumption,leisure)
{alpha * log( consumption) + (1-alpha) * log( leisure)}

########################################################################################
# define an empty 3X3 matrix for Mu
Mu <- matrix(0, nrow = num_type, ncol = num_type)

for (ii in 1:num_type) { #start for ii 
  for (jj in 1:num_type) {  #start for jj 
    # select match <i,j> as subsample
    
     #ii <- 1
     #jj <- 1
    data_marriage_sub <- subset(data_marriage, ( (data_marriage$type==jj & data_marriage$type_spouse==ii )))
    
    # if there is any match of this type, we can find Pareto weight of this match
    if (nrow(data_marriage_sub)!=0){   #start if
      
      loop_check <- 0
      print("Start")
      while(loop_check==0)  {  #start while
        Mu_temp <- matrix(0, nrow = num_type, ncol = num_type)
        for (i in 1:num_type) {  #start for i 
          for (j in 1:num_type) {   #start for j
            
            #i <- 1
            #j <- 3
            
            data_marriage_sub <- subset(data_marriage, ( (data_marriage$type==j & data_marriage$type_spouse==i )))
            
            if (nrow(data_marriage_sub)!=0){   #start if
              
              mean_Y <- mean(data_marriage_sub$Y)
              mean_W <- mean(data_marriage_sub$W)
              mean_W_spouse <- mean(data_marriage_sub$W_spouse)
              mean_U_single <- mean(data_marriage_sub$U_single)
              mean_U_single_spouse <- mean(data_marriage_sub$U_single_spouse)

              alpha_j <- mean(data_marriage_sub$alpha_j)
              alpha_i <- mean(data_marriage_sub$alpha_i)
              
              C_j  <- function (Mu) alpha_j *(1-Mu)* (mean_Y+mean_W+mean_W_spouse)
              l_j  <- function (Mu) ((1-alpha_j)*Mu* (mean_Y+mean_W+mean_W_spouse))/(2*mean_W)
              
              C_i  <- function (Mu) alpha_i *(1-Mu)* (mean_Y+mean_W+mean_W_spouse)
              l_i  <- function (Mu) ((1-alpha_i)*Mu* (mean_Y+mean_W+mean_W_spouse))/(2*mean_W_spouse)

              
              U_married_i  <- function (Mu) log(C_i(Mu)^alpha_i * l_i(Mu)^(1-alpha_i))
              U_married_j  <- function (Mu) log(C_j(Mu)^alpha_j * l_j(Mu)^(1-alpha_j))
              
              A <-sigma_theta*(Ln_TauS_guess[ii]-Ln_TauD_guess[jj])
              B <- log(mean_U_single)
              C <- log(mean_U_single_spouse)
              D1 <- alpha_j * log( alpha_j * (mean_Y+mean_W+mean_W_spouse)) 
              D2 <- (1-alpha_j) * log((1-alpha_j)*(mean_Y+mean_W+mean_W_spouse))/(2*mean_W)
              D3 <- alpha_i * log( alpha_i * (mean_Y+mean_W+mean_W_spouse)) 
              D4 <- (1-alpha_i) * log((1-alpha_i)*(mean_Y+mean_W+mean_W_spouse))/(2*mean_W_spouse)
              D <- D1+D2-D3-D4
              Marriage_market_clear <- function (Mu) alpha_i*log(Mu/(1-Mu)) -A -B +C +D
              curve(Marriage_market_clear(x), 0, 1)
              uni <- uniroot(Marriage_market_clear,c(0,1))$root
              abline(h = 0, lty = 3)
              points(uni, 0, pch = 16, cex = 2)
              Mu_temp[i,j] <- uni
              
            }  # end if
            else {Mu_temp[i,j]<- 0}
          }   #end for j
        }#end for i
        print("end")
        print(i)
        print(j)
        # step 3: updating
        #.......
        
        mean_Y <- mean(data_marriage_sub$Y)
        mean_W <- mean(data_marriage_sub$W)
        mean_W_spouse <- mean(data_marriage_sub$W_spouse)
        C_i  <- function (Mu) alpha_i *(1-Mu)* (mean_Y+mean_W+mean_W_spouse)
        l_i  <- function (Mu) ((1-alpha_i)*Mu* (mean_Y+mean_W+mean_W_spouse))/(2*mean_W_spouse)
        C_j  <- function (Mu) alpha_j *(1-Mu)* (mean_Y+mean_W+mean_W_spouse)
        l_j  <- function (Mu) ((1-alpha_j)*Mu* (mean_Y+mean_W+mean_W_spouse))/(2*mean_W)
        U_married_i  <- function (Mu) log(C_i(Mu)^alpha_i * l_i(Mu)^(1-alpha_i))
        U_married_j  <- function (Mu) log(C_j(Mu)^alpha_j * l_j(Mu)^(1-alpha_j))
        
        sum_ii <- 0
        for (h in 1:num_type) {sum_ii <- sum_ii + exp(U_married_i(Mu_temp[i,h])/sigma_theta)}
        sum_jj <- 0
        for (h in 1:num_type) {sum_jj <- sum_jj + exp(U_married_j(Mu_temp[h,j])/sigma_theta)}
        
        # for (k in 1:8) {mean_U_single_ii <- aggregate(data_marriage_sub$U_single ~ data_marriage_sub$sex, FUN=mean)}
        # for (k in 1:8) {mean_U_single_jj <- aggregate(data_marriage_sub$U_single ~ data_marriage_sub$sex, FUN=mean)}
        data_marriage_sub   <- subset(data_marriage, ( (data_marriage$type==jj & data_marriage$type_spouse==ii & data_marriage$sex=="Female")))
        #data_marriage_sub_m <- subset(data_marriage_sub, (sex=="Male")  )
        #data_marriage_sub_f <- subset(data_marriage_sub, (sex=="Female")  )
        mean_U_single_ii    <- log(mean(data_marriage_sub$U_single_spouse))
        mean_U_single_jj    <- log(mean(data_marriage_sub$U_single))
        
        
        TauD_iijj      <- f[jj]*(exp(U_married_j(Mu_temp[ii,jj])/sigma_theta))/(exp(mean_U_single_jj/sigma_theta)+sum_jj)
        Ln_TauD_guess_new <- log(TauD_iijj)-((U_married_j(Mu_temp[ii,jj])-mean_U_single_jj)/sigma_theta)
        Ln_TauD_guess_new <- ifelse(Ln_TauD_guess_new<0 , 0.1, Ln_TauD_guess_new)
        TauS_iijj      <- m[ii]*(exp(U_married_i(Mu_temp[ii,jj])/sigma_theta))/(exp(mean_U_single_ii/sigma_theta)+sum_ii)
        Ln_TauS_guess_new <- log(TauS_iijj)-((U_married_i(Mu_temp[ii,jj])-mean_U_single_ii)/sigma_theta)
        Ln_TauS_guess_new <- ifelse(Ln_TauS_guess_new<0 , 0.1, Ln_TauS_guess_new)
        
        #TauD_guess_new <- exp(Ln_TauD_guess_new)
        #TauS_guess_new <- exp(Ln_TauS_guess_new)
        #if (Ln_TauD_guess_new>0 & Ln_TauS_guess_new>0){
        Sys.sleep(0.1)
        print("Ln_TauD_guess:")
        print(Ln_TauD_guess[jj])
        print("Ln_TauS_guess:")
        print(Ln_TauS_guess[ii])
        
        
        epsilonS <- abs(Ln_TauS_guess_new - Ln_TauS_guess[ii])
        epsilonD <- abs(Ln_TauD_guess_new - Ln_TauD_guess[jj])
        
        epsilon  <-  epsilonS + epsilonD
        epsilon[is.na(epsilon)] <- 0
        
        loop_check <- ifelse(epsilon<.01,1,0)
        
        Ln_TauD_guess[jj] <- Ln_TauD_guess_new
        Ln_TauS_guess[ii] <- Ln_TauS_guess_new
        
        #} else {
        #    Ln_TauD_guess[jj] <- ifelse(Ln_TauD_guess_new>0 , Ln_TauD_guess_new, 2*Ln_TauD_guess[jj] )
        #  Ln_TauS_guess[ii] <- ifelse(Ln_TauS_guess_new>0 , Ln_TauS_guess_new, 2*Ln_TauS_guess[ii])
        #}
        
        
        Sys.sleep(0.1)
        print("epsilon:")
        print(epsilon)
        print("Ln_TauD_guess_new:")
        print(Ln_TauD_guess_new)
        print("Ln_TauD_guess:")
        print(Ln_TauD_guess[jj])
        print("Ln_TauS_guess:")
        print(Ln_TauS_guess[ii])
        print("Next")
        
        
        # back to step 2 if loop_check==0   
        
      }  #end while(loop_check==0)
    }  # end if
    else {Mu_temp[ii,jj]<- 0}  # if (nrow(data_marriage_sub)==0)
    Mu[ii,jj] <- Mu_temp[ii,jj]
    Sys.sleep(0.1)
    print(ii)
    print(jj)
    print(epsilon)
    print(Mu[ii,jj])
    print("++++++++++++++++++++++++++++++++++++++++++")
  }#end for jj 
}#end for ii 

write.matrix(Mu, file = "D:/Desktop/Education/2nd year pape/Estimation Parameters/3 Pareto weights/Mu.dta", sep = " ", blocksize)