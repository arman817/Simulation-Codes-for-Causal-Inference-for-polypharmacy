rm(list=ls(all=TRUE))

args=commandArgs();
####use this to divide to different jobs to parallel run ####
start=as.numeric(args[4]);
end=as.numeric(args[5]);

#This code computes the causal estimates for the data using 
#Softmax Regression and SVM for the entire data
##############################################################
#Code for computing causal estimates for Data Generation Scenario II as mentioned in Section 2 of the Supplementary Materials
##############################################################
library(bindata)
library(softmaxreg)
library(twang)
library(e1071)


# Generate the data using the sample size and the seed value and return the
# data frame with the covariates and the outcome varialble and the top regimen
# along with the G_ computation values
Generate_Data <- function(n,seed_value,a,b,c,d,e,f,g,h){
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
  X_9 <- rnorm(n,0,1)
  X_10 <- runif(n,0,1)
  X_11 <- runif(n,0,1)
  X_12 <- rnorm(n,0,1)
  X_13 <- runif(n,0,1)
  X_14 <- runif(n,0,1)
  
  #Medication Generation
  pA_1 <- plogis(0.5*X_1 + 2.3*X_2 + 0.74*X_3 - 1.8*X_4)
  pA_2 <- plogis(0.4*X_1 + 2.5*X_2 + 0.87*X_3 - 1.5*X_5)
  pA_3 <- plogis(0.8*X_1 - 0.4*X_2 + 0.34*X_3 + 0.8*X_6)
  pA_4 <- plogis(0.7*X_1 - 0.65*X_2 + 0.55*X_3 + 1.1*X_7)
  pA_5 <- plogis(1.9*X_1 - 0.33*X_2 + 1.67*X_3 - 0.6*X_8)
  pA_6 <- plogis(1.5*X_1 - 0.53*X_2 + 1.88*X_3 - 0.3*X_9)
  pA_7 <- plogis(-0.85*X_1 + 1.35*X_2 - 0.99*X_3 + 1.4*X_7)
  pA_8 <- plogis(-0.65*X_1 + 1.03*X_2 - 0.875*X_3 + 1.1*X_6)
  
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
  
  #Considering positive binary correlations between A_5, A_6 as 0.15 and A_7, A_8 as 0.25
  B <- matrix(nrow=n,ncol=4)
  for (i in 1:n){
    prob <- c(pA_5[i],pA_6[i],pA_7[i],pA_8[i])
    B[i,] <- rmvbin(1, margprob = prob, bincorr = m)
  }
  
  #Separating the array obtained above into separate medications
  Gap <- as.data.frame(A)
  A_1 <- Gap[,1]
  A_2 <- Gap[,2]
  A_3 <- Gap[,3]
  A_4 <- Gap[,4]
  
  Gap <- as.data.frame(B)
  A_5 <- Gap[,1]
  A_6 <- Gap[,2]
  A_7 <- Gap[,3]
  A_8 <- Gap[,4]
  
  #Generating the outcome model
  pY <- plogis(-1.3 + 0.8*X_1 + 0.89*X_2*A_2 + 1.5*X_11*A_8
               - 1.8*X_13 + 0.9*X_14*A_1*A_8 + 0.4*A_1*A_2*A_5 + 7.2*(1-A_2)*X_9*X_12 - 15*(1-A_7)*A_6
               - 12*X_9*(1-A_7) - X_10*(1-A_8)*(1-A_7) + 5*(1-A_8)*(1-A_1) + 5*(1-A_1)*X_12*X_14)
  Y <- rbinom(n, 1, prob=pY)
  
  #Creating a binary regimen value for each regimens and a value for Site (which denotes 
  #the class of each observation) for each observation generated in the 
  #data which will be later required for the multi-class classification
  
  Regimen <- matrix(0, nrow= n, ncol = 256)
  Site <- c()
  for(i in 1:n)
    Site[i] <- 0
  for(i in 1:n){
    Sum <- (128*(A_1[i])+ 64*(A_2[i]) + 32*(A_3[i]) + 16*(A_4[i]) + 8*(A_5[i])
            + 4*(A_6[i]) + 2*(A_7[i]) + A_8[i] + 1)
    Site[i] <- Sum - 1
    Regimen[i,Sum] <- 1
  }
  Count_Variable <- matrix(0,nrow = 256,ncol = 1)
  for(i in 1:256){
    Count_Variable[i,1] <- sum(Regimen[,i])
  }
  Freq <- NA
  for(i in 1:n){
    Freq[i] <- Count_Variable[Site[i]+1,1]
  }
  
  #Indicating the regimens of our interest i.e. R_1 = [1,1,1,1,1,1,1,1]
  R_1 <- Regimen[,256]
  
  #Evaluating G Computation 2 for the simulated data
  data1 <- data.frame(X_1,X_2,X_3,X_4,X_5,X_6,X_7,X_8,X_9,X_10,X_11,X_12,X_13,X_14,
                      A_1,A_2,A_3,A_4,A_5,A_6,A_7,A_8,Y)
  model <- glm(Y~.,data1,family = "binomial")
  Cov <- data1[,1:14]
  data2 <- data.frame(Cov,A_1=a,A_2=b,A_3=c,A_4=d,A_5=e,A_6=f,A_7=g,A_8=h,Y)
  G_comp1 <- predict.glm(model,type="response",newdata=data2)
  
  
  #Checking for confounding in our simulation data
  data3 <- data.frame(A_1,A_2,A_3,A_4,A_5,A_6,A_7,A_8,Y)
  model1 <- glm(Y~.,data = data3, family = "binomial")
  data4 <- data.frame(A_1=a,A_2=b,A_3=c,A_4=d,A_5=e,A_6=f,A_7=g,A_8=h,Y)
  Conf_1 <- predict.glm(model1,type="response",newdata = data4)
  
  #Evaluating G Computation 1 for the simulated data
  data5 <- data.frame(X_1,X_2,X_3,X_4,X_5,X_6,X_7,X_8,X_9,X_10,X_11,X_12,X_13,X_14,Y,R_1)
  data6 <- data5[which(data5$R_1 == 1),]
  data7 <- data6[,1:15]
  model2 <- glm(Y~.,data = data7,family="binomial")
  data8 <- data5[,1:15]
  G_comp2 <- predict.glm(model2,type="response",newdata = data8)
  
  
  data9 <- cbind(Site,X_1,X_2,X_3,X_4,X_5,X_6,X_7,X_8,X_9,X_10,X_11,X_12,X_13,X_14,
                 A_1,A_2,A_3,A_4,A_5,A_6,A_7,A_8,Y,R_1,G_comp1,Conf_1,Freq,G_comp2)
  return(data9)
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
Propensity_Score_Function2 <- function(Y,medication,pi,a,b,c,d,e,f,g,h){
  data <- data.frame(medication,pi,Y)
  model <- glm(Y~.,data = data, family = "binomial")
  data1 <- data.frame(A_1=a,A_2=b,A_3=c,A_4=d,A_5=e,A_6=f,A_7=g,A_8=h,pi)
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


SVM_Value <- function(data_obt,Original_Data,n){
  Original_Data <- data.frame(Original_Data)
  x2 <- Original_Data[,2:15]
  R_1 <- Original_Data[,25]
  Y <- Original_Data[,24]
  G_comp1 <- Original_Data[,29]
  G_comp2 <- Original_Data[,26]
  Medication <- Original_Data[,16:23]
  
  
  ClassVar <- data_obt[,1]
  Class <- data.frame(as.factor(ClassVar))
  colnames(Class) <- "Classes"
  Cov <- data_obt[,2:15]
  
  MainData <- data.frame(Cov,Class)
  attach(MainData)
  model <- svm(Classes ~ ., data = MainData)
  x <- subset(MainData, select = -Classes)
  y <- Classes
  model <- svm(x, y, probability = TRUE)
  pred <- predict(model, x2, decision.values = TRUE, probability = TRUE)
  Prob_Obtained <- attr(pred, "probabilities")
  Prob_Obtained <- data.frame(Prob_Obtained)
  
  #Obtain GPS from Prob_Obtained
  pi1 <- Prob_Obtained$X255
  
  #Causal Estimates for Regimen 1
  IPTW_1 <- IPTW_function(Y,R_1,pi1)
  
  Prop_Adjusted_1 <- Propensity_Score_Function1(Y,R_1,pi1)
  Prop_Adjusted_2 <- Propensity_Score_Function2(Y,Medication,pi1,1,1,1,1,1,1,1,1)
  
  TMLE_Compute_1 <- TMLE_function(Y,G_comp1,pi1,R_1)
  TMLE_Compute_2 <- TMLE_function(Y,G_comp2,pi1,R_1)
  
  #Combine the estimate values in a row
  Solution_SVM <- cbind(IPTW_1,Prop_Adjusted_1,Prop_Adjusted_2,TMLE_Compute_1,TMLE_Compute_2)
  return(Solution_SVM)
}

Softmax_Regression_Value <- function(data_obt,Original_Data,n){
  Original_Data <- data.frame(Original_Data)
  Xhi <- Original_Data[,2:15]
  R_1 <- Original_Data[,25]
  Y <- Original_Data[,24]
  G_comp1 <- Original_Data[,29]
  G_comp2 <- Original_Data[,26]
  Medication <- Original_Data[,16:23]
  
  data2 <- data_obt[,1:15]
  dat2 <- as.data.frame(data2)
  x2 = dat2[,2:15]
  y2 = dat2$Site
  softmax_model_L1 = softmaxReg(x2, y2, hidden = NULL, maxit = 100, type = "class",
                                algorithm = "rmsprop", L2 = TRUE, penalty = 1e-4, batch = 40)
  #summary(softmax_model2)
  yFitMat <- predict.softmax1(softmax_model_L1,Xhi)
  yFitMat <- data.frame(yFitMat)
  
  #Obtain GPS using softmax regression
  pi1 <- yFitMat$X255
  
  #Causal Estimates for Regimen 1
  IPTW_1 <- IPTW_function(Y,R_1,pi1)
  
  Prop_Adjusted_1 <- Propensity_Score_Function1(Y,R_1,pi1)
  Prop_Adjusted_2 <- Propensity_Score_Function2(Y,Medication,pi1,1,1,1,1,1,1,1,1)
  
  TMLE_Compute_1 <- TMLE_function(Y,G_comp1,pi1,R_1)
  TMLE_Compute_2 <- TMLE_function(Y,G_comp2,pi1,R_1)
  
  #Combine the estimate values in a row
  Solution_softmax <- cbind(IPTW_1,Prop_Adjusted_1,Prop_Adjusted_2,TMLE_Compute_1,TMLE_Compute_2)
  return(Solution_softmax)
}

#For the Data Generation Scenario 2, a=1,b=1,c=1,d=1,e=1,f=1,g=1,h=1
func <- function(n,Seed_Val,a,b,c,d,e,f,g,h){
  Simulate_Data <- Generate_Data(n,Seed_Val,a,b,c,d,e,f,g,h)
  est1 <- Softmax_Regression_Value(Simulate_Data,Simulate_Data,n)
  est2 <- SVM_Value(Simulate_Data,Simulate_Data,n)
  est <- c(est1,est2)
  return(est)
}

#n denotes the number of observation in each study. As mentioned in Section 2 of the Supplementary Materials
#simulations were carried out by setting the values of n as 500.
n <- 500
for(p in start:end){
  Val <- func(n,p,1,1,1,1,1,1,1,1)
  outpath=paste("/home/arman13/Medication/Output/Error_Re/Second/Sim_",start,"_",end,".txt",sep="")
  write.table(Val, outpath, row.names = F, col.names=F, sep="\t")
}
