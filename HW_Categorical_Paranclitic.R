source("Functions/functions.R")

#Variables to set and data
{
  file_location <- "Data/"
  result_folder <- "20170615_Density/"
  file_name <- "Synthetic_WeightBaseCaseControl.csv"
  first_columns=4 #column number of first parameter -1
  number_of_parameters=30 #Excluding categorical variables - these should be in columns at the end
  column_of_case_control=1
  case_number <- 1 #code for test cases.
  control_number <- 0 #code for control cases
  contour_number <- 2 #code for network controls
  auc_case <- 3 #Cases fro making the weights. This should be contour_number +1. NA = not using
  number_of_categoricals <- 9
  number_of_indexes <- 17
  bestmarker <- c("X1","X2","X3") #names of top 3 single markers
  Best_col <- 5 #Best marker
  Second_col <- 6 #second best marker
  Third_col <- 7 #third best marker
  grid_size <- 16 #size of the density matrix grid (e.g. 10 = 10x10 grid)
  IndicesInModel <- 1 #the number of indices to include in the model
  senstoset <- 0.95 #Set the specificity for logistic regression testing.
  TotalConnectionsPlot <- 1 #Plot graphs of total connections for each threshold (1=yes)
  
  data=read.csv(paste(file_location,file_name,sep=""),na.strings=c("NA","#VALUE!",""))
  data <- data[,-35] #Removing a categorical that was artificially set to 100% accuracy.
  data <- data[which(data$V3==9),] #Use only samples closest to diagnosis.
  
  #Examply code for imputing missing values using Amelia algorithm.
  # #Impute Missing Values
  # set.seed=123
  # 
  # data_imp <- amelia(data[,(first_columns+1):ncol(data)],m=5)
  # data_imp <- (data_imp$imputations$imp1+data_imp$imputations$imp2+
  #                data_imp$imputations$imp3 +
  #                data_imp$imputations$imp4 +
  #                data_imp$imputations$imp5)/5
  # data[,(first_columns+1):ncol(data)] <- data_imp
  # 
  # #Generate Levels and Data sets
  # data <- data[-which(data$Group=="BD Control"),];data <- droplevels(data)
  # data$Group <- factor(data$Group,labels=c(1,0))
  # data$Group <- as.numeric(as.character(data$Group))
  # data$Group[sample(which(data$Group==0),340)] <- 2
  # data$Group[sample(which(data$Group==1),20)] <- 3
}

#Program 
{  
  st=0.5 #starting threshold
  d=0.05 #increments
  loops=1
  #number of iterations
  
  #Script
  {
    aucii=matrix(0,loops,5)
    for(aucii_index in 1:loops){
      
      threshold = (aucii_index-1)*d+st;
      
      List_of_results <- parenclitic(data,result_folder,first_columns,number_of_parameters,column_of_case_control,case_number,control_number,contour_number,threshold,number_of_categoricals,number_of_indexes,Best_col,Second_col,grid_size,IndicesInModel,senstoset,TotalConnectionsPlot,Third_col,auc_case,bestmarker)
      
      aucii[aucii_index,1]=List_of_results$thres
      aucii[aucii_index,2]=List_of_results$auc
      aucii[aucii_index,3]=List_of_results$specificity
      aucii[aucii_index,4]=List_of_results$sensitivity
      aucii[aucii_index,5]=List_of_results$formula
    }
    
    colnames(aucii) <- c("thres","auc","specificity","sensitivity","formula")
    
    write.csv(aucii,file=paste(result_folder,"aucii.csv",sep=""),row.names=F)
  }
}

#AUC and Sensitivity Graph
{
  data_thres_auc <- read.csv(paste(result_folder,"aucii.csv",sep=""))
  
  sens <- ggplot(data_thres_auc, aes(x=thres, y=sensitivity)) +
    geom_point(shape=1, size=2, colour="red") +
    xlab("Threshold") +
    ylab("Sensitivity") +
    ggtitle("PC Controls \n Late cases") +
    theme(title=element_text(size=8))
  
  auc <- ggplot(data_thres_auc, aes(x=thres, y=auc)) +
    geom_point(shape=1, size=2, colour="red") +
    xlab("Threshold") +
    ylab("AUC") +
    ggtitle("PC Controls \n Late cases") +
    theme(title=element_text(size=8))
  
  plot_grid(sens,auc)  
  
  ggsave(paste(result_folder,"Sens_AUC.png",sep=""),width=8,height=3)
}