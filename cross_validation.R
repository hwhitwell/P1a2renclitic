#packages
{
  library(ROCR)
  library(Hmisc)
  library(Amelia)
}

#functions
{
  crossvalidation <- function(data,column_of_case_control,column_of_ModelIndex,number_of_categoricals,spectoset,save_plot=0,CA125.OR.Rule=0)
  {
    Res <- as.data.frame(data) #Convert numeric columns to numeric
    for (i in 2:(ncol(Res)-number_of_categoricals)){
      Res[,i] <- as.numeric(as.character(Res[,i]))
    }
    
    if(number_of_categoricals>0){ #Convert categorical columns to categorical
      for (i in (ncol(Res)-number_of_categoricals):ncol(Res)){
        levels(Res[,i]) <- c("N","Y")
      }
    }
    
    if(OrRuleUsed==1){ #Convert OrRule column to categorical
      Res$CA125.OR.Rule <- factor(Res$CA125.OR.Rule,labels=c("N","Y"))
    }
    
    Data <- Res[,c(column_of_case_control,column_of_ModelIndex)] #Extract case/control and relevant columns
    
    colnames(Data)[1]="output"
    Result <- matrix(data=NA,nrow=nrow(Data),ncol=2) #Matrix to store prediction values
    Result[,1] <- Data[,1]
    for (i in 1:nrow(Data)){ #Generate prediction value based on leave one out.
      Data1 <- Data[-i,]
      mylogit <- glm(output~.,data=Data1,family="binomial")
      prGLM <- predict(mylogit,newdata=Data[i,],type="response",se=TRUE)$fit
      Result[i,2] <- prGLM
    }
    
    auc=somers2(Result[,2],Result[,1])[1] #Calculates AUC
    pred <- prediction(Result[,2],Result[,1]) #Converts to format for performance calculations
    perf<-performance(pred,"tpr","fpr") #Calculates performance calculations
    ind=which((1-perf@x.values[[1]])>=spectoset)
    num=ind[length(ind)] #Index of the value closest to the set specificity
    sens=perf@y.values[[1]][num] #sensitivity at the set specificity
    ind=which(perf@y.values[[1]]==sens)
    spec=max(1-perf@x.values[[1]][ind]) 
    res=list(auc=auc,spec=spec,sens=sens)
    par(mar=c(2,2,2,2)) #because margin error kept occuring
    temp <- tryCatch(plot(perf),error=function(e) 1) #if the plot fails, it will still report values
    if (is.numeric(temp)==T){
      print("Plot failed")
    } else {
      if(save_plot==1) {png(paste(file_location, substr(file_name,1,nchar(file_name)-4),"_sens",spectoset,"_Indices",length(column_of_ModelIndex),"_ROC.png",sep=""))}
      plot(perf); abline(a=0,b=1,col="red",lty=2)
      if(save_plot==1) {dev.off()}
    }
    par(mar=c(5.1,4.1,4.1,2.1)) #return margin to original size
    return(res)
  }
}

#variables to set
{
  file_location <- "20170408_linear/"
  file_name <- "IndexValues_2.csv"
  column_of_case_control <- 2
  column_of_ModelIndex <- c(7) #vector of columns with the model indexes
  number_of_categoricals <- 0
  spectoset <- 0.95 #Set the sensitivity for logistic regression testing.
  save_plot <- 0 #(1=yes, 0=No)
  OrRuleUsed <- 0 #(1=yes, other =no)
  first_columns <- 3
  #SampleNumber <- 7
  
  number_of_indicies <- 1
  
  data=read.csv(paste(file_location,file_name,sep="")) #Read data
  #data$X0 <- rep(1:9,nrow(data)/9)
  #data <-subset(data,Group==1|Group==0)
  
  column_of_ModelIndex <- combn(3:ncol(data),number_of_indicies)
}

#program
{
  #data <- subset(data,X0==SampleNumber | X0.1==0)
  Results <- data.frame(auc=NA,spec=NA,sens=NA)
  for(i in 1:ncol(column_of_ModelIndex)){
    Results <- rbind(Results,
                     as.data.frame(crossvalidation(data,column_of_case_control,column_of_ModelIndex[,i],number_of_categoricals,spectoset,save_plot,OrRuleUsed)))
  }
  Results <- Results[-1,]
  Index <- vector()
  for(i in 1:ncol(column_of_ModelIndex)){
    Index <- c(Index,paste(colnames(data)[column_of_ModelIndex[,i]],collapse="+"))
  }
  Results$Index <- Index
  write.csv(Results,paste(file_location,"CrossValSingleMarkerResults_Thresh_2",".csv",sep=""))
}