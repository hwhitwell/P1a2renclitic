#Packages
{
  library(MASS) 
  library(igraph)
  library(lubridate)
  library(plotrix)
  library(ROCR)
  library(Hmisc)
  library(pROC)
  library(ggplot2)
  library(cowplot)
  library(emdbook)
  library(compiler)
  library(reshape)
  library(Amelia)
}

#Functions
{
  #This function calculates the distance for categorical variables. Base controls are subsetted by the categorical value
  #found at the intersection of 'women_index' and 'cat_column' in 'samples'. con_samples is generated and contains the
  #subsetted values of the continuous data (cat_column). The 1D kernal density estimation contains 'grid_size' number of 
  #bins. 
  cat_cal <- function(grid_size,women_index,cat_column,con_column,base_controls,samples)
  {
    con_samples <- base_controls[which(base_controls[,cat_column]==samples[women_index,cat_column]),con_column]
    density <- density(con_samples,n=grid_size) #1D kernal estimation of density
    
    max_x <- max(density$x) #Minumum value of marker
    min_x <- min(density$x) #Maximum value of marker
    
    x_diff <- abs(samples[women_index,con_column]-density$x)
    
    if(samples[women_index,con_column]<min_x |
       samples[women_index,con_column]>max_x){ #Check to see if sample is outside the density plot
      
      #outside calculation
      if(samples[women_index,con_column]<min_x){ #x=marker value, y=density
        nearest_y <- density$y[1]
      } else {
        nearest_y <- density$y[grid_size]
      }
      
      internal_dist <- abs(max(density$y)-nearest_y)
      xx=internal_dist
      yy=max(density$y)
      
      if(abs(-xx^2*yy*pi)>0.000001){
        if((-xx^2*yy*pi)>(-1/exp(1))){
          sig2=-xx^2/(2*lambertW(-xx^2*yy*pi,b=-1))
        } else sig2=xx^2/2
        
        x_mu <- abs(max(density$x)-samples[women_index,con_column])
        
        temp <- (1.0-exp(-0.5*(x_mu^2/sig2))) 
      } else {
        temp = 1 
      }
      
    } else {#inside calculation
      temp <- density$y[which(x_diff==min(x_diff))]
      temp <- sum(density$y[which(density$y>temp)])
      temp <- temp*(max_x-min_x)
      temp <- temp/grid_size}
    
    return(temp)
  }
  
 #For adding significance * to plots
  func <- function(x)
  {
    if(is.na(x)==T){
      return(" ")
    }
    if(x<0.001){
      return("***")
    } else if(x<0.01) {
      return("**")
    } else if(x<0.05){
      return("*")
    } else return(" ")
  }
  
#Calculation of AUC, Spec, Sens for data with 1st column outcome and other column - predictors
auc_spec_sens <- function(data,Spec)
  {colnames(data)[1]="output"
  mylogit=glm(output~.,data=data,family="binomial") #ERROR OCCURING HERE
  temp <- tryCatch(predict(mylogit,type="response",se=TRUE),error=function(e) 1) #This puts result to 0 if prediction fails
  if (is.numeric(temp)==T){
    res <- list(auc=0,spec=0,sens=0)
  } else{
    prGLM <- predict(mylogit,type="response",se=TRUE)
    Tmp=cbind(data[,1],prGLM$fit)
    auc=somers2(Tmp[,2],Tmp[,1])[1]
    pred1<-prediction(Tmp[,2],Tmp[,1])
    perf1<-performance(pred1,"tpr","fpr")
    ind=which((1-perf1@x.values[[1]])>=Spec)
    num=ind[length(ind)]
    sens=perf1@y.values[[1]][num]
    ind=which(perf1@y.values[[1]]==sens)
    spec=max(1-perf1@x.values[[1]][ind])
    res=list(auc=auc,spec=spec,sens=sens)
  }
  return(res)
  }
  
#looking for parameters for best linear combination: Data - datafile, y - number of column with outcome, x - columns with predictors, x_max - maximal number of variables in the final model, type - "AUC" or "Sens" (in this case "Spec" value will be used)
best_comb_linear <- function(Data, y, x, num_ind, type, Spec)
  {sens_max=0
  spec_max=0
  auc_max=0
  k<-0; ind_max<-rep(0, num_ind)
  i<-num_ind	#looking through all combinations of "i" variables
  tmp=combn(x,i)
  for (j in 1:length(tmp[1,]))
  {data=data.frame(Data[,c(y,tmp[,j])])
  ROC_calc=auc_spec_sens(data,Spec)
  
  if (type=="AUC")
    if (ROC_calc$auc>auc_max)
    {auc_max=ROC_calc$auc
    sens_max=ROC_calc$sens
    spec_max=ROC_calc$spec
    ind_max=tmp[,j]
    }
  if (type=="Sens")
    if (ROC_calc$sens>sens_max)
    {auc_max=ROC_calc$auc
    sens_max=ROC_calc$sens
    spec_max=ROC_calc$spec
    ind_max=tmp[,j]
    }		   
  }
  
  res<-list(auc_max=auc_max, sens_max=sens_max, spec_max=spec_max, ind_max=ind_max)
  }
  
#Calculates the "density" value for points outside the kernal density estimation.
outside <- function(ContourMatrix=contour_matrix[[1]][[2]][[1]],x,y)
  {#x and y co-ordinate of the most dense point
    mu_x <-  ContourMatrix$x[which(ContourMatrix$z==max(ContourMatrix$z),arr.ind=T)[1]]
    mu_y <-  ContourMatrix$y[which(ContourMatrix$z==max(ContourMatrix$z),arr.ind=T)[2]]
    #Gradient and intercep for the positive and negative diagonals through the density matrix
    pos_m <- (range(ContourMatrix$y)[1]-range(ContourMatrix$y)[2])/(range(ContourMatrix$x)[1]-range(ContourMatrix$x)[2])
    pos_c <- range(ContourMatrix$y)[1] - range(ContourMatrix$x)[1] * pos_m
    neg_m <- (range(ContourMatrix$y)[1]-range(ContourMatrix$y)[2])/(range(ContourMatrix$x)[2]-range(ContourMatrix$x)[1])
    neg_c <- range(ContourMatrix$y)[1] - range(ContourMatrix$x)[2] * neg_m
    
    #Determin the correct edge
    if(y>(pos_m * x + pos_c)){
      if(y>(neg_m * x + neg_c)){
        sector <- "upper"
      } else sector <- "left"
    } else if(y>(neg_m * x + neg_c)){
      sector <- "right"
    } else sector <- "lower"
    
    #Determine the x and y co-ord of the nearest edge point
    if(sector=="upper"){
      nearest_y <- max(ContourMatrix$y)
      samp_mu <- y - mu_y
      near_mu <- nearest_y - mu_y
      x_exact <- mu_x + ((x-mu_x) / (samp_mu/near_mu))
      x_diff <- abs(x_exact - ContourMatrix$x)
      nearest_x <- min(ContourMatrix$x[which(x_diff==min(x_diff))])
    }
    
    if(sector=="lower"){
      nearest_y <- min(ContourMatrix$y)
      samp_mu <- y - mu_y
      near_mu <- nearest_y - mu_y
      x_exact <- mu_x + ((x-mu_x) / (samp_mu/near_mu))
      x_diff <- abs(x_exact - ContourMatrix$x)
      nearest_x <- min(ContourMatrix$x[which(x_diff==min(x_diff))])
    }
    
    if(sector=="left"){
      nearest_x <- min(ContourMatrix$x)
      samp_mu <- x - mu_x
      near_mu <- nearest_x - mu_x
      y_exact <- mu_y + ((y-mu_y) / (samp_mu/near_mu))
      y_diff <- abs(y_exact - ContourMatrix$y)
      nearest_y <- min(ContourMatrix$y[which(y_diff==min(y_diff))])
    }
    
    if(sector=="right"){
      nearest_x <- max(ContourMatrix$x)
      samp_mu <- x - mu_x
      near_mu <- nearest_x - mu_x
      y_exact <- mu_y + ((y-mu_y) / (samp_mu/near_mu))
      y_diff <- abs(y_exact - ContourMatrix$y)
      nearest_y <- min(ContourMatrix$y[which(y_diff==min(y_diff))])
    }
    
    #Calculate distance
    internal_dist <- sqrt((mu_x-nearest_x)^2 + (mu_y-nearest_y)^2)
    xx=internal_dist
    yy=ContourMatrix$z[which(ContourMatrix$x==nearest_x),which(ContourMatrix$y==nearest_y)]
    
    if(abs(-xx^2*yy*pi)>0.000001){
      if((-xx^2*yy*pi)>(-1/exp(1))){
        sig2=-xx^2/(2*lambertW(-xx^2*yy*pi,b=-1))
      } else sig2=xx^2/2
      
      x_mu <- sqrt((x-mu_x)^2+(y-mu_y)^2)
      
      dens_val <- (1.0-exp(-0.5*(x_mu^2/sig2))) 
    } else {
      dens_val=1 
    }
    
    return(dens_val)
  }; outside <- cmpfun(outside)
  
#Calculates density estimation for each covariate pair in case/controls.
contour_calculation <- function(contour_matrix,samples,grid_size)
  {
    max_x <- max(contour_matrix$x)
    max_y <- max(contour_matrix$y)
    min_x <- min(contour_matrix$x)
    min_y <- min(contour_matrix$y)
    x_diff <- abs(samples[1,1]-contour_matrix$x)
    y_diff <- abs(samples[1,2]-contour_matrix$y)
    
    if(samples[1,1]>max_x | samples[1,2]>max_y | samples[1,1]<min_x | samples[1,2]<min_y){temp <- outside(contour_matrix,samples[1,2],samples[1,1])
    } else {
      temp <- min(contour_matrix$z[which(x_diff==min(x_diff)),which(y_diff==min(y_diff))])
      temp <- sum(contour_matrix$z[which(contour_matrix$z>=temp)])
      temp <- temp*(max_y-min_y)*(max_x-min_x)
      temp <- temp/(grid_size*grid_size)}
    
    return(temp)
  }
  
#Calls all other functions, performs network analysis and generates results files
  {
    if(dir.exists(result_folder)==FALSE){dir.create(result_folder)}
    
    time1 <- now()
    samples <- subset(data,data[,column_of_case_control]==control_number|data[,column_of_case_control]==case_number)
    base_controls <- subset(data,data[column_of_case_control]==contour_number)
    
    
    print("Generating Contour Matrix")
    contour_matrix <- list(0)
    contour_matrix <- list(replicate(number_of_parameters,contour_matrix))
    contour_matrix <- list(replicate(number_of_parameters,contour_matrix))
    
    for(i in 1:number_of_parameters){
      for(j in 1:number_of_parameters){
        #If interquartile range = 0, then use 1.06*sigma*n^1/5. Otherwise use bandwidth.nrd
        if(
          quantile(base_controls[,c(first_columns+i)],0.25)==quantile(base_controls[,c(first_columns+i)],0.75) |
          quantile(base_controls[,c(first_columns+j)],0.25)==quantile(base_controls[,c(first_columns+j)],0.74)
        ){
          bandwidth <- apply(base_controls[,c(first_columns+i,first_columns+j)],2,function(x)1.06*sd(x)*(length(x)^(1/5)))
          contour_matrix[[1]][[i]][[j]] <- kde2d(base_controls[,first_columns+i],base_controls[,first_columns+j], h=bandwidth, n=grid_size)
        } else {
          contour_matrix[[1]][[i]][[j]] <- kde2d(base_controls[,first_columns+i],base_controls[,first_columns+j], n=grid_size)
        }        
      }
    }
    print(paste("Contour Matrix Completed",round(seconds(interval(time1,now())),2)))
    
    Res<-matrix(0,1,2+number_of_indexes+number_of_categoricals) #This is to hold the topological indexes.
    Res2<-matrix(0,1,2+number_of_parameters+number_of_categoricals) #This is to hold the number of connections for each parameter.
    
    #Generate the network for each sample.
    print(paste("Generating Paranclitic Networks, Threshold=",threshold))
    pdf(paste(result_folder,"Contour_",threshold,"_GridSize_",grid_size,".pdf",sep=""))
    
    layout(matrix(c(1:25),1,1))
    par(mar=c(1,1,1,1)*1)
    
    for(women_index in 1:nrow(samples)){
      network=matrix(0,number_of_parameters,number_of_parameters+number_of_categoricals)
      for(i in 1:number_of_parameters){
        for(j in 1:number_of_parameters){
          if(j!=i){
            network[i,j] <- contour_calculation(contour_matrix[[1]][[i]][[j]],samples[women_index,c(first_columns+i,first_columns+j)],grid_size)
          }
        }
      }
      if(number_of_categoricals>0){
        for(i in (number_of_parameters+1):(number_of_parameters+number_of_categoricals)){
          for(j in 1:number_of_parameters){
            network[j,i] <- cat_cal(grid_size,women_index,(i+first_columns),(j+first_columns),base_controls,samples)
          }
        } 
      }
      
      #To make network square
      if(number_of_categoricals>0){
        network_extension <- t(network[,(number_of_parameters+1):(number_of_parameters+number_of_categoricals)])
        network_extension <- cbind(network_extension,matrix(0,number_of_categoricals,number_of_categoricals))
        network <- rbind(network,network_extension)
      }
      
      #Generate weigth subset, determin AUC, generaet matrix and apply to the network values.
      if(is.na(auc_case)==F){
        auc_samples <- subset(data,data[,column_of_case_control]==auc_case | data[,column_of_case_control]==contour_number)
        ParameterAUC <- apply(auc_samples[,(first_columns+1):(first_columns+number_of_parameters+number_of_categoricals)],2,function(x) auc(auc_samples[,column_of_case_control],x))
        AUCMatrix <- matrix(ParameterAUC,(number_of_parameters+number_of_categoricals),1) %*% matrix(ParameterAUC,1,(number_of_parameters+number_of_categoricals))
        network <- network*AUCMatrix
      }
      
      
      distance <- network[which(network>threshold)]; #For making distance indexes (13 and 14)
      #network <- round(network,digits=4)
      
      network <- ifelse(network-threshold>0.0,1,0)
      
      for(i in 1:number_of_parameters){
        for(j in 1:number_of_parameters){
          network[i,j]=max(network[i,j],network[j,i])
        }
      }
      
      m=as.matrix(network)
      inet=graph.adjacency(m,mode="undirected",weighted=TRUE)
      
      if(samples[women_index,column_of_case_control]>0) { V(inet)$color="red"
      } else {
        V(inet)$color="green"
      }
      
      V(inet)$color[c(Best_col,Second_col, Third_col) - first_columns] <- "blue"
      
      if(number_of_categoricals>0){
        V(inet)$color[(number_of_parameters+1):(number_of_parameters+number_of_categoricals)] <- "yellow"
      }
      
      plot.igraph(inet,layout=layout.fruchterman.reingold,vertex.size=10,edge.color="black",main=paste(women_index,"--",samples[women_index,1],"--",samples[women_index,"Name"],"_",samples[women_index,"Early.Late"]))
      box(lty='1373',col='black')
      print(paste(women_index, round(seconds(interval(time1,now())),2)))
      
      #Calculate The Results for Logistic Regression
      Results<-matrix(0,1,2+number_of_indexes+number_of_categoricals)
      Results[1, 1] = as.character(samples[women_index, 1])
      Results[1, 2] = samples[women_index, column_of_case_control]
      
      Results[1, 3] = max(degree(inet)) # Index 1) max degree.
      Results[1,4] = mean(degree(inet)) # Index 2) mean degree.
      
      shortest_paths=0
      shortest_paths = as.vector(shortest.paths(inet))
      shortest_paths = shortest_paths[shortest_paths != 0]
      shortest_paths = shortest_paths[shortest_paths != Inf]
      
      if(length(shortest_paths) > 0)
      {
        Results[1,5] = mean(1/shortest_paths)*length(shortest_paths)/((number_of_parameters)*((number_of_parameters)-1))
      } else {
        Results[1,5] = 0 # Index 3) inetwork efficiency.
      }
      
      var_betweenness = betweenness(inet)
      if(max(var_betweenness) > 0)
        Results[1,6] = max(var_betweenness) # Index 4) Max betweeness.
      
      Results[1,7] = mean(var_betweenness) # Index 5) Mean betweeness.
      
      var_closeness = closeness(inet)
      if(max(var_closeness) > 0)
        Results[1,8] = max(var_closeness) # Index 6) Max closeness.
      
      Results[1,9] = mean(var_closeness) # Index 7) Mean closenss.
      
      var_page_rank = page.rank(inet)$vector # Google PageRank Score.
      if(max(var_page_rank) > 0)
        Results[1, 10] = max(var_page_rank) # Index 8) Max Google PageRank Score.
      
      Results[1, 11] = mean(var_page_rank)# Index 9) MeanGoogle PageRank Score.
      
      var_authority_score = 0; # Kleinberg's centrality score.
      #var_authority_score = authority.score(inet)$vector
      #if(max(var_authority_score) > 0)
      #Results[1, 12] = max(var_authority_score) # Index 10) max Kleinberg's centrality score.
      
      #Results[1, 13] = mean(var_authority_score)# Index 11) Kleinberg's centrality score.
      
      Results[1, 14] = sum(degree(inet)) # Index12) total number of connections.
      
      if(length(distance)>0){
        Results[1,15] = max(distance)} #Index 13) maximum distance from central population (above threshold).
      
      if(length(distance)>0){
        Results[1,16] = mean(distance)} #Index 14) mean distance from central population (above threshold).
      
      Results[1,17] = sum(inet[(Best_col-first_columns),]) #Index 15) Total number of connections to best marker
      
      Results[1,18] = sum(inet[(Second_col-first_columns),]) #Index 16) Total number of connections to second best marker
      
      Results[1,19] = sum(inet[(Third_col-first_columns),]) #Index 17) Total number of connections to third best marker
      
      #Categorical variables
      if(number_of_categoricals>0){
        for (i in 1:number_of_categoricals){
          Results[1,(2+number_of_indexes+i)] <- as.character(samples[women_index,(first_columns+number_of_parameters+i)])
        }
      }
      
      Res=rbind(Res,Results[1,])
      
      #Calculate The Number of Connections Between Each Marker
      Results2<-matrix(0,1,2+number_of_parameters+number_of_categoricals)
      Results2[1, 1] = as.character(samples[women_index, 1])
      Results2[1, 2] = samples[women_index, column_of_case_control]
      
      Results2[1,3:ncol(Results2)] = apply(inet[],1,sum)
      Res2=rbind(Res2,Results2[1,])
    }
    
    dev.off()
    print(paste("Networks complete",round(seconds(interval(time1,now())),2)))
    
    #Generate Result Plots
    Res <- Res[-1,]
    
    pdf(paste(result_folder,"Results_",threshold,"_GridSize_",grid_size,".pdf",sep=""))
    
    layout(matrix(c(1:20),2,2))
    par(mar=c(1.4,1.4,1.4,1.4)*3)
    
    if(number_of_categoricals>0){
      top_labels=c("0","0","Max degree","Mean degree","inet efficiency","Max betweenness","Mean betweenness","Max closeness","Mean closeness", "Max Page.Rank", "Mean Page.Rank", "Max Keinberg's centrality", "Mean Max Keinberg's centrality","Total degree","Max Distance", "Mean Distance",paste("Connections to",bestmarker[1]),paste("Connections to", bestmarker[2]),paste("Connections to",bestmarker[3]),colnames(samples[(ncol(samples)-number_of_categoricals+1):ncol(samples)]))
    } else {
      top_labels=c("0","0","Max degree","Mean degree","inet efficiency","Max betweenness","Mean betweenness","Max closeness","Mean closeness", "Max Page.Rank", "Mean Page.Rank", "Max Max Keinberg's centrality", "Mean Max Keinberg's centrality","Total degree","Max Distance", "Mean Distance",paste("Connections to",bestmarker[1]),paste("Connections to", bestmarker[2]),paste("Connections to",bestmarker[3]))
    }
    
    
    for (i in 3:(ncol(Results))){
      for (j in 3:(ncol(Results))){
        if(j>i){
          plot(Res[,i],Res[,j],col=as.factor(Res[,2]),xlab=top_labels[i],ylab=top_labels[j],cex=1.4,pch=1)
        }
      }
    }
    
    dev.off()
    colnames(Res) <- top_labels
    write.csv(Res, paste(result_folder,"IndexValues_",threshold,"_GridSize_",grid_size,".csv",sep=""), row.names = FALSE)
    
    print(paste("Result Plots Completed",round(seconds(interval(time1,now())),2)))
    
    #Genereate Total Connections Plot
    if(TotalConnectionsPlot==1){
      colnames(Res2) <- c("PromiseID","CaseControl",colnames(samples)[(first_columns+1):(first_columns+number_of_parameters+number_of_categoricals)])
      Res2 <- Res2[-1,]
      write.csv(Res2, paste(result_folder,"Connections_",threshold,".csv",sep=""),row.names=F)
      
      Res2 <- melt(as.data.frame(Res2),id.vars=c("PromiseID","CaseControl"),variable_name="Marker")
      Res2$value <- as.numeric(as.character(Res2$value))
      
      markers <- unique(Res2$Marker)
      results <- matrix(NA,length(markers),3)
      for(i in 1:length(markers)){
        temp <- subset(Res2,Marker==markers[i])
        temp2 <- tryCatch(t.test(data=temp, value~CaseControl),error=function(e) "NA")
        if(temp2=="NA"){
          results[i,] <- c(mean(temp[which(temp$CaseControl==0),"value"]),mean(temp[which(temp$CaseControl==0),"value"]),1)}
        else{
          results[i,] <- c(temp2$estimate[1],temp2$estimate[2],temp2$p.value)
        }
      }
      
      colnames(results) <- c("Control","Case","PValue")
      results <- as.data.frame(results)
      results$Marker <- markers
      results <- results[,c("Marker","Control","Case","PValue")]
      results$sig <- lapply(results$PValue,func) #significance level
      results$Adjusted <- results$PValue*length(markers) #bonferoni correction
      results$sig.adjusted <- lapply(results$Adjusted,func) #new significance level
      
      ggplot(Res2, aes(x=Marker,y=value)) +
        geom_boxplot(aes(fill=as.factor(CaseControl))) +
        theme(axis.text.x=element_text(angle=90,hjust=1)) +
        annotate(geom="text",x=1:nrow(results),y=-1,label=unlist(results$sig),size=4) +
        annotate(geom="text",x=1:nrow(results),y=-2,label=unlist(results$sig.adjusted),size=4)+
        annotate(geom="text",x=c(1.3,1.5),y=c(-1,-2),label=c("Sig","Bonf"),size=4)
      
      ggsave(paste(result_folder,"Connections_",threshold,".png",sep=""),width=20,height=10)
    }
    
    
    #Generating Linear Regressions and Best Models
    Res <- as.data.frame(Res)
    for (i in 2:(ncol(Res)-number_of_categoricals)){
      Res[,i] <- as.numeric(as.character(Res[,i]))
    }
    
    if(number_of_categoricals>0){
      for (i in (ncol(Res)-number_of_categoricals):ncol(Res)){
        levels(Res[,i]) <- c("N","Y")
      }
    }
    
    number_of_indices<-IndicesInModel #NUMBER OF INDICES YOU WANT TO CHOOSE   
    comb<- best_comb_linear(Res, 2, 3:ncol(Res), number_of_indices, "Sens", spectoset)
    A<- "Res[,2]~"
    for (i in 1:number_of_indices)
      A<-paste(A, "+Res[,",comb$ind_max[i],"]")
    print(A)
    
    mylogit<- glm(as.formula(A), family=binomial(link="logit"), na.action=na.pass)
    prGLM <- predict(mylogit,type="response",se=TRUE)
    Tmp=cbind(Res[,2],prGLM$fit)
    
    statit=somers2(Tmp[,2],Tmp[,1])
    
    
    
    print(statit[1])
    
    pred1<-prediction(Tmp[,2],Tmp[,1])
    perf1<-performance(pred1,"tpr","fpr")
    auc1 <- performance(pred1,"auc")
    
    ind=which((1-perf1@x.values[[1]])>=spectoset)
    num=ind[length(ind)]
    spec=1-perf1@x.values[[1]][num]
    sens=perf1@y.values[[1]][num]
    roc1 <- roc(Tmp[,1],Tmp[,2])
    
    Res_ret<-list(thres=threshold,auc=statit[1],specificity=spec,sensitivity=sens,formula=A, roc_parencl=roc1)
    
    print(paste("Completed",round(seconds(interval(time1,now())),2)))
    
    return(Res_ret)
    
  }
}
