rm(list=ls(all=TRUE))

args=commandArgs();
####use this to divide to different jobs to parallel run ####
start=as.numeric(args[4]);
end=as.numeric(args[5]);

##############################################################
#Code for computing causal estimates for Data Generation Scenario I as mentioned in Section 3.1 
##############################################################
library(bindata)
library(softmaxreg)
library(e1071)
library(twang)

# Generate the data using the sample size and the seed value and return the
# data frame with the covariates and the outcome varialble and the top two regimens
# along with the G_ computation values
Generate_Data <- function(n,seed_value,a,b,c,d,e,f,g,h,j,k){
  #Covariate Generation
  set.seed(seed_value)
  X_1 <- runif(n,0,1)
  X_2 <- runif(n,0,1)
  X_3 <- runif(n,0,1)
  X_4 <- runif(n,0,1)
  X_5 <- runif(n,0,1)
  X_6 <- runif(n,0,1)
  X_7 <- runif(n,0,1)
  X_8 <- runif(n,0,1)
  X_9 <- runif(n,0,1)
  X_10 <- runif(n,0,1)
  X_11 <- runif(n,0,1)
  X_12 <- runif(n,0,1)
  
  #Medication Generation
  pA_1 <- plogis(0.7*X_1 - 0.87*X_3 + X_5 + X_7 + 0.8*X_2 - 1.6*X_4)
  pA_2 <- plogis(0.9*X_1 - 0.66*X_3 + 0.8*X_5 + X_7 - 0.8*X_6)
  pA_3 <- plogis(1.4*X_1 - 0.45*X_3 + 0.4*X_5 - X_7 - X_8)
  pA_4 <- plogis(1.8*X_1 - 0.95*X_3 - 0.1*X_5 - X_7 - X_12)
  
  m <- diag(4)
  m[1,2] <- 0.15
  m[2,1] <- 0.15
  m[4,3] <- 0.25
  m[3,4] <- 0.25
  
  #Considering positive binary correlations between A_1, A_2 as 0.15 and A_3, A_4 as 0.25
  A <- matrix(nrow=n,ncol=4)
  for (i in 1:n){
    prob <- c(pA_1[i],pA_2[i],pA_3[i],pA_4[i])
    A[i,] <- rmvbin(1, margprob = prob, bincorr = m)
  }
  
  #Separating the array obtained above into separate medications
  Gap <- as.data.frame(A)
  A_1 <- Gap[,1]
  A_2 <- Gap[,2]
  A_3 <- Gap[,3]
  A_4 <- Gap[,4]
  
  #Generating the outcome model
  pY <- plogis(-1.8 + 1.2*X_1*A_3 + 2.15*X_3*A_1 - 0.5*X_5*X_9 + 3.89*X_7*X_1 
               + X_10*A_2 - 1.5*X_11*A_4 + 4.2*(1-A_1)*A_3 + 15*(1-A_2)*A_4)
  Y <- rbinom(n, 1, prob=pY)
  
  #Creating a binary regimen value for each regimens and a value for Site (which denotes 
  #the class of each observation) for each observation generated in the 
  #data which will be later required for the multi-class classification
  Regimen <- matrix(0, nrow= n, ncol = 16)
  Site <- c()
  for(i in 1:n)
    Site[i] <- 0
  for(i in 1:n){
    Sum <- 8*(A_1[i]) + 4*(A_2[i]) + 2*(A_3[i]) + A_4[i] + 1
    Site[i] <- Sum - 1
    Regimen[i,Sum] <- 1
  }
  
  #Indicating the regimens of our interest i.e. R_1 = [1,1,0,0] and R_2 = [1,1,1,1]
  R_1 <- Regimen[,j]
  R_2 <- Regimen[,k]
  
  #Evaluating G Computation 2 for the simulated data
  data1 <- data.frame(X_1,X_2,X_3,X_4,X_5,X_6,X_7,X_8,X_9,X_10,X_11,X_12,A_1,A_2,A_3,A_4,Y)
  model <- glm(Y~.,data1,family = "binomial")
  Cov <- data1[,1:12]
  
  #For regimen R_1: [1,1,0,0]
  data2 <- data.frame(Cov,A_1=a,A_2=b,A_3=c,A_4=d,Y)
  G_comp1 <- predict.glm(model,type="response",newdata=data2)
  
  #For regimen R_2: [1,1,1,1]
  data3 <- data.frame(Cov,A_1=e,A_2=f,A_3=g,A_4=h,Y)
  G_comp2 <- predict.glm(model,type="response",newdata=data3)
  
  
  #Checking for confounding in our simulation data
  data4 <- data.frame(A_1,A_2,A_3,A_4,Y)
  model1 <- glm(Y~.,data = data4, family = "binomial")
  data5 <- data.frame(A_1=a,A_2=b,A_3=c,A_4=d,Y)
  Conf_1 <- predict.glm(model1,type="response",newdata = data5)
  data6 <- data.frame(A_1=e,A_2=f,A_3=g,A_4=h,Y)
  Conf_2 <- predict.glm(model1,type="response",newdata = data6)
  
  
  #Evaluating G Computation 1 for the simulated data
  data7 <- data.frame(X_1,X_2,X_3,X_4,X_5,X_6,X_7,X_8,X_9,X_10,X_11,X_12,Y,R_1)
  data8 <- data7[which(data7$R_1 == 1),]
  data9 <- data8[,1:13]
  model2 <- glm(Y~.,data = data9,family="binomial")
  data10 <- data7[,1:13]
  G_comp3 <- predict.glm(model2,type="response",newdata = data10)
  
  data11 <- data.frame(X_1,X_2,X_3,X_4,X_5,X_6,X_7,X_8,X_9,X_10,X_11,X_12,Y,R_2)
  data12 <- data11[which(data11$R_2 == 1),]
  data13 <- data12[,1:13]
  model3 <- glm(Y~.,data = data13,family="binomial")
  data14 <- data11[,1:13]
  G_comp4 <- predict.glm(model3,type="response",newdata = data14)
  
  
  data11 <- cbind(Site,X_1,X_2,X_3,X_4,X_5,X_6,X_7,X_8,X_9,X_10,X_11,X_12,A_1,A_2,A_3,
                  A_4,Y,R_1,R_2,G_comp1,G_comp2,Conf_1,Conf_2,G_comp3,G_comp4)
  return(data11)
}

#Function to evaluate the true values of the causal estimates or the probability of treatment success 
True_Val <- function(n,m,o){
  True_Value_Regimen_1 <- 0
  True_Value_Regimen_2 <- 0
  for(l in m:(9+m)){
    set.seed(l)
    C_1 <- runif(n,0,1)
    C_2 <- runif(n,0,1)
    C_3 <- runif(n,0,1)
    C_4 <- runif(n,0,1)
    C_5 <- runif(n,0,1)
    C_6 <- runif(n,0,1)
    C_7 <- runif(n,0,1)
    C_8 <- runif(n,0,1)
    C_9 <- runif(n,0,1)
    C_10 <- runif(n,0,1)
    C_11 <- runif(n,0,1)
    C_12 <- runif(n,0,1)
    
    pY1 <- plogis(-1.8 + 2.15*C_3 - 0.5*C_5*C_9 + 3.89*C_7*C_1 + C_10)
    Y1 <- rbinom(o, 1, prob=pY1)
    Expected_Y1 <- sum(Y1)/o
    True_Value_Regimen_1 <- True_Value_Regimen_1 + Expected_Y1/10
    
    pY2 <- plogis(-1.8 + 1.2*C_1 + 2.15*C_3 - 0.5*C_5*C_9 + 3.89*C_7*C_1 + C_10 - 1.5*C_11)
    Y2 <- rbinom(o, 1, prob=pY2)
    Expected_Y2 <- sum(Y2)/o
    True_Value_Regimen_2 <- True_Value_Regimen_2 + Expected_Y2/10
  }
  True_Value <- data.frame(True_Value_Regimen_1,True_Value_Regimen_2)
  return(True_Value)
}

#Function to evalute the causal estimate using the IPTW
IPTW_function <- function(Y,R,pi){
  new_data <- data.frame(Y,R,pi)
  new_data1 <- new_data[which(new_data$R==1),]
  w <- 1/new_data1$pi
  model <- lm(Y~1,weights = w,data=new_data1)
  kappak <- summary(model)$coefficient[1]
  return(kappak)
}

#Function to evalute the causal estimate using the PSA(I)
Propensity_Score_Function1 <- function(Y,R,pi){
  data <- data.frame(R,pi,Y)
  model <- glm(Y~.,data = data, family = "binomial")
  data1 <- data.frame(R=1,pi)
  Prop_Adj <- predict.glm(model,type="response",newdata = data1)
  PAS_Val <- mean(Prop_Adj)
  return(PAS_Val)
}

#Function to evalute the causal estimate using the PSA(II)
Propensity_Score_Function2 <- function(Y,medication,pi,a,b,c,d){
  data <- data.frame(medication,pi,Y)
  model <- glm(Y~.,data = data, family = "binomial")
  data1 <- data.frame(A_1=a,A_2=b,A_3=c,A_4=d,pi)
  Prop_Adj <- predict.glm(model,type="response",newdata = data1)
  PAS_Val <- mean(Prop_Adj)
  return(PAS_Val)
}

#Function to evaluate the causal estimate using TMLE
#Use of G_comp1 and G_comp2 gives rise to TMLE(I) and TMLE(II)
TMLE_function <- function(Y,Q,pi,Regimen_1){
  weight <- 1/pi
  data <- data.frame(weight,Q)
  data1 <- data.frame(weight,Q,Y,Regimen_1)
  data2 <- data1[which(data1$Regimen_1 == 1),]
  #Update Step in the regression
  model <- glm(Y ~ offset(qlogis(Q)),data = data2, weight=weight,family = "binomial")
  #Predict Step
  prediction <- predict.glm(model, type = "response",newdata=data) #Q_star
  TMLE_Val <- mean(prediction)
  return(TMLE_Val)
}

#Function to predict the GPS for a model obtained by Softmax Regression
predict.softmax1 <- function(object, newdata, ...){
  yMat = convertClass2Matrix(as.factor(as.character(object$data$y)))
  param = object$weights
  type = object$type
  funName = object$funName
  nObs = dim(newdata)[1]
  idx = c(1:nObs)
  curFitList = lapply(idx, FUN = forwardPropParallel, newdata, W = param$W, B = param$B, funName)
  yFitMat = convertList2Matrix(curFitList)
  colnames(yFitMat) = colnames(yMat)
  return (yFitMat)
}

Softmax_Regression_Value <- function(data_obt,n,a,b,c,d,e,f,g,h){
  R_1 <- data_obt[,19]
  R_2 <- data_obt[,20]
  Y <- data_obt[,18]
  G_comp1 <- data_obt[,25]
  G_comp2 <- data_obt[,26]
  G_comp3 <- data_obt[,21]
  G_comp4 <- data_obt[,22]
  Medication <- data_obt[,14:17]
  
  data2 <- data_obt[,1:13]
  dat2 <- as.data.frame(data2)
  x2 = dat2[,2:13]
  y2 = dat2$Site
  softmax_model_L1 = softmaxReg(x2, y2, hidden = NULL, maxit = 100, type = "class",
                                algorithm = "rmsprop", L2 = TRUE, penalty = 1e-4, batch = 20)
  yFitMat <- softmax_model_L1$fitted.values
  yFitMat2 <- data.frame(yFitMat)
  
  #Obtain GPS using softmax regression
  pi1 <- yFitMat2[,5]      #Refers to j=13
  pi2 <- yFitMat2[,8]      #Refers to k=16
  
  #Causal Estimates for Regimen 1
  IPTW_1 <- IPTW_function(Y,R_1,pi1)
  
  Prop_Adjusted_1 <- Propensity_Score_Function1(Y,R_1,pi1)
  Prop_Adjusted_2 <- Propensity_Score_Function2(Y,Medication,pi1,a,b,c,d)
  
  TMLE_Compute_1 <- TMLE_function(Y,G_comp1,pi1,R_1)
  TMLE_Compute_2 <- TMLE_function(Y,G_comp3,pi1,R_1)
  
  
  
  #Causal Estimates for Regimen 2
  IPTW_11 <- IPTW_function(Y,R_2,pi2)
  
  Prop_Adjusted_11 <- Propensity_Score_Function1(Y,R_2,pi2)
  Prop_Adjusted_12 <- Propensity_Score_Function2(Y,Medication,pi2,e,f,g,h)
  
  TMLE_Compute_11 <- TMLE_function(Y,G_comp2,pi2,R_2)
  TMLE_Compute_12 <- TMLE_function(Y,G_comp4,pi2,R_2)
  
  #Combine the estimate values in a row
  Solution_softmax <- cbind(IPTW_1,Prop_Adjusted_1,Prop_Adjusted_2,TMLE_Compute_1,TMLE_Compute_2
                            ,IPTW_11,Prop_Adjusted_11,Prop_Adjusted_12,TMLE_Compute_11,TMLE_Compute_12)
  return(Solution_softmax)
}

SVM_Value <- function(data_obt,n,a,b,c,d,e,f,g,h){
  R_1 <- data_obt[,19]
  R_2 <- data_obt[,20]
  Y <- data_obt[,18]
  G_comp1 <- data_obt[,25]
  G_comp2 <- data_obt[,26]
  G_comp3 <- data_obt[,21]
  G_comp4 <- data_obt[,22]
  Medication <- data_obt[,14:17]
  ClassVar <- data_obt[,1]
  Class <- data.frame(as.factor(ClassVar))
  colnames(Class) <- "Classes"
  Cov <- data_obt[,2:13]
  
  MainData <- data.frame(Cov,Class)
  attach(MainData)
  model <- svm(Classes ~ ., data = MainData)
  x <- subset(MainData, select = -Classes)
  y <- Classes
  model <- svm(x, y, probability = TRUE)
  pred <- predict(model, x, decision.values = TRUE, probability = TRUE)
  Prob_Obtained <- attr(pred, "probabilities")
  Prob_Obtained <- data.frame(Prob_Obtained)
  
  #Obtain GPS from Prob_Obtained
  pi1 <- Prob_Obtained$X12
  pi2 <- Prob_Obtained$X15
  
  
  
  #Causal Estimates for Regimen 1
  IPTW_1 <- IPTW_function(Y,R_1,pi1)
  
  Prop_Adjusted_1 <- Propensity_Score_Function1(Y,R_1,pi1)
  Prop_Adjusted_2 <- Propensity_Score_Function2(Y,Medication,pi1,a,b,c,d)
  
  TMLE_Compute_1 <- TMLE_function(Y,G_comp1,pi1,R_1)
  TMLE_Compute_2 <- TMLE_function(Y,G_comp3,pi1,R_1)
  
  
  
  #Causal Estimates for Regimen 2
  IPTW_11 <- IPTW_function(Y,R_2,pi2)
  
  Prop_Adjusted_11 <- Propensity_Score_Function1(Y,R_2,pi2)
  Prop_Adjusted_12 <- Propensity_Score_Function2(Y,Medication,pi2,e,f,g,h)
  
  TMLE_Compute_11 <- TMLE_function(Y,G_comp2,pi2,R_2)
  TMLE_Compute_12 <- TMLE_function(Y,G_comp4,pi2,R_2)
  
  #Combine the estimate values in a row
  Solution_SVM <- cbind(IPTW_1,Prop_Adjusted_1,Prop_Adjusted_2,TMLE_Compute_1,TMLE_Compute_2
                        ,IPTW_11,Prop_Adjusted_11,Prop_Adjusted_12,TMLE_Compute_11,TMLE_Compute_12)
  return(Solution_SVM)
}

GBM_Value <- function(data_obta,n,a,b,c,d,e,f,g,h){
  R_1 <- data_obta[,19]
  R_2 <- data_obta[,20]
  Y <- data_obta[,18]
  G_comp1 <- data_obta[,25]
  G_comp2 <- data_obta[,26]
  G_comp3 <- data_obta[,21]
  G_comp4 <- data_obta[,22]
  Medication <- data_obta[,14:17]
  Cov <- data_obta[,2:13]
  data_req <- data.frame(data_obta)
  
  #Creating a GBM for Regimen 1
  ps.1 <- ps(R_1 ~ X_1+X_2+X_3+X_4+X_5+X_6+X_7+X_8+X_9+X_10+X_11+X_12,
             data = data_req,
             # choosing stopping criterion as the Kolmogovov's Statistic
             stop.method = stop.methods[("ks.stat.max")],
             plots=NULL,
             pdf.plots=FALSE,
             #gbm options,
             n.trees = 10000,
             interaction.depth = 2,
             shrinkage = 0.01,
             perm.test.iters = 0,
             verbose = FALSE)
  
  #Creating a GBM for Regimen 2
  ps.2 <- ps(R_2 ~ X_1+X_2+X_3+X_4+X_5+X_6+X_7+X_8+X_9+X_10+X_11+X_12,
             data = data_req,
             stop.method = stop.methods[("ks.stat.max")],
             plots=NULL,
             pdf.plots=FALSE,
             n.trees = 10000,
             interaction.depth = 2,
             shrinkage = 0.01,
             perm.test.iters = 0,
             verbose = FALSE)
  
  #Obtain GPS using GBM
  pi1 <- ps.1$ps
  pi2 <- ps.2$ps
  
  colnames(pi1) <- c("weight")
  colnames(pi2) <- c("weight")
  
  #Causal Estimates for Regimen 1
  new_data13 <- data.frame(Y,R_1,pi1)
  new_data14 <- new_data13[which(new_data13$R_1==1),]
  w <- 1/new_data14$weight
  model <- lm(Y~1,weights = w,data=new_data14)
  IPTW_1 <- summary(model)$coefficient[1]
  
  Prop_Adjusted_1 <- Propensity_Score_Function1(Y,R_1,pi1)
  Prop_Adjusted_2 <- Propensity_Score_Function2(Y,Medication,pi1,a,b,c,d)
  
  TMLE_Compute_1 <- TMLE_function(Y,G_comp1,pi1,R_1)
  TMLE_Compute_2 <- TMLE_function(Y,G_comp3,pi1,R_1)
  
  
  
  #Causal Estimates for Regimen 2
  new_data23 <- data.frame(Y,R_2,pi2)
  new_data24 <- new_data23[which(new_data23$R_2==1),]
  w <- 1/new_data24$weight
  model <- lm(Y~1,weights = w,data=new_data24)
  IPTW_11 <- summary(model)$coefficient[1]
  
  Prop_Adjusted_11 <- Propensity_Score_Function1(Y,R_2,pi2)
  Prop_Adjusted_12 <- Propensity_Score_Function2(Y,Medication,pi2,e,f,g,h)
  
  TMLE_Compute_11 <- TMLE_function(Y,G_comp2,pi2,R_2)
  TMLE_Compute_12 <- TMLE_function(Y,G_comp4,pi2,R_2)
  
  #Combine the estimate values in a row
  Solution_Twang <- cbind(IPTW_1,Prop_Adjusted_1,Prop_Adjusted_2,TMLE_Compute_1,TMLE_Compute_2
                          ,IPTW_11,Prop_Adjusted_11,Prop_Adjusted_12,TMLE_Compute_11,TMLE_Compute_12)
  return(Solution_Twang)
}

#For the Data Generation Scenario 1, a=1,b=1,c=0,d=0,e=1,f=1,g=1,h=1,j=13,k=16

func <- function(n,Seed_Val,a,b,c,d,e,f,g,h,j,k){
  #Simulate Data with a given seed value
  Simulate_Data <- Generate_Data(n,Seed_Val,a,b,c,d,e,f,g,h,j,k)
  #Estimate the probability of treatment successes using SVM, Softmax Regression and GBM
  est1 <- Softmax_Regression_Value(Simulate_Data,n,a,b,c,d,e,f,g,h)
  est2 <- SVM_Value(Simulate_Data,n,a,b,c,d,e,f,g,h)
  est3 <- GBM_Value(Simulate_Data,n,a,b,c,d,e,f,g,h)
  est <- c(est1,est2,est3)
  return(est)
}

#n denotes the number of observation in each study. As mentioned in Section 3.3
#simulations were carried out by setting the values of n as 500 and 1000.
n <- 500
for(p in start:end){
  Val <- func(n,p,1,1,0,0,1,1,1,1,13,16)
  outpath=paste("/home/arman13/Medication/Output/Error_Re/Firsta/Sim_",start,"_",end,".txt",sep="")
  write.table(Val, outpath, row.names = F, col.names=F, sep="\t")
}
