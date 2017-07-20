source("Functions/functions.R")

#Variables to set and data
{ #Braces to collapse code section
  file_location <- "Data/" #Path to data files.
  result_folder <- "20170615_Density/" #Where to store output files. This folder will be created if it doesn't already exist.
  file_name <- "Synthetic_WeightBaseCaseControl.csv"
  first_columns <- 4 #column number of first parameter -1
  number_of_continuous <- 30 #Number of continuous data columns. These should be after the sample meta data in your .csv file.
  number_of_categoricals <- 9 #Number of categorical variables. These should be the final columns in your .csv file.
  column_of_case_control <- 1 #Column that cases/controls are defined in.
  case_number <- 1 #code for test cases.
  control_number <- 0 #code for control cases
  base_number <- 2 #code for network controls. These are used to building the underlying density distributions.
  auc_case <- 3 #Cases for making the weights. This should be base_number +1. NA = not using
  number_of_indexes <- 17 #Parenclitic function contains 17 indicies. This should be left at 17 unless the function "parenclitic" is edited.
  bestmarker <- c("X1","X2","X3") #names of top 3 single markers
  Best_col <- 5 #Best marker
  Second_col <- 6 #second best marker
  Third_col <- 7 #third best marker
  grid_size <- 16 #size of the density matrix grid (e.g. 10 = 10x10 grid)
  IndicesInModel <- 1 #the number of indices to include in the model
  spectoset <- 0.95 #Set the specificity for logistic regression testing.
  TotalConnectionsPlot <- 1 #Plot graphs of total connections for each threshold (1=yes)
  
  #Threshold settings. These are required for "Program" section below.
  st=0.5 #Starting threshold
  d=0.05 #Increment. If scanning over multiple thresholds.
  loops=1 #Number of iterations. If performing analysis at one threshold (st), set this to 1.
  
  #Read in the data.
  data=read.csv(paste(file_location,file_name,sep=""),na.strings=c("NA","#VALUE!",""))
  
  #####FOR EXAMPLE DATA#####
  data <- data[,-35] #Removing a categorical that was artificially set to 100% accuracy.
  data <- data[which(data$V3==9),] #Use only samples closest to diagnosis.
  ##########################
  
  #Example code for imputing missing values using Amelia algorithm.
  # #Impute Missing Values
  # set.seed=123
  # 
  # data_imp <- amelia(data[,(first_columns+1):ncol(data)],m=5)
  # data_imp <- (data_imp$imputations$imp1+data_imp$imputations$imp2+
  #                data_imp$imputations$imp3 +
  #                data_imp$imputations$imp4 +
  #                data_imp$imputations$imp5)/5
  # data[,(first_columns+1):ncol(data)] <- data_imp
}

#Program 
{#This loops over thresholds (st) loops number of times, calling "parenclitic" each time. 
  #The threshold is increased by an increment "d" each time.
  #The best model for each threshold, the AUC, Spec and Sens are returned and written in "aucii.csv".
  #If the file "aucii.csv" is not written, the data is stored in the global environment in the matrix aucii.
  aucii=matrix(0,loops,5)
  for(aucii_index in 1:loops){
    
    threshold = (aucii_index-1)*d+st;
    
    List_of_results <- parenclitic(data,result_folder,first_columns,number_of_continuous,column_of_case_control,case_number,control_number,base_number,threshold,number_of_categoricals,number_of_indexes,Best_col,Second_col,grid_size,IndicesInModel,spectoset,TotalConnectionsPlot,Third_col,auc_case,bestmarker)
    
    aucii[aucii_index,1]=List_of_results$thres
    aucii[aucii_index,2]=List_of_results$auc
    aucii[aucii_index,3]=List_of_results$specificity
    aucii[aucii_index,4]=List_of_results$sensitivity
    aucii[aucii_index,5]=List_of_results$formula
  }
  
  colnames(aucii) <- c("thres","auc","specificity","sensitivity","formula")
  
  if(file.exists(paste0(result_folder,"aucii.csv"))==T){
    temp <- as.matrix(read.csv(paste0(result_folder,"aucii.csv")))
    aucii <- rbind(temp,aucii)
  }
  write.csv(aucii,file=paste(result_folder,"aucii.csv",sep=""),row.names=F) #Contains results.
}

#AUC and Sensitivity Graph - if a number of threshold have been tried, this generates a graph for visual selection of the optimal threshold.
{
  data_thres_auc <- read.csv(paste(result_folder,"aucii.csv",sep=""))
  
  sens <- ggplot(data_thres_auc, aes(x=thres, y=sensitivity)) +
    geom_point(shape=1, size=2, colour="red") +
    xlab("Threshold") +
    ylab("Sensitivity") +
    ylim(c(0,1)) +
    ggtitle("Sensitivity") +
    theme(title=element_text(size=8))
  
  auc <- ggplot(data_thres_auc, aes(x=thres, y=auc)) +
    geom_point(shape=1, size=2, colour="red") +
    xlab("Threshold") +
    ylab("AUC") +
    ylim(c(0,1)) +
    ggtitle("AUC") +
    theme(title=element_text(size=8))
  
  plot_grid(sens,auc)  
  
  ggsave(paste(result_folder,"Sens_AUC.png",sep=""),width=8,height=3)
}