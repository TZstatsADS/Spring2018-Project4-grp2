data_cleaning_MS <- function(data_address){
  train <- read.csv(paste(data_address,'/data_train.csv',sep=""))
  test <- read.csv(paste(data_address,'/data_test.csv',sep=""))
  
  sub_cleaning <- function(data){
    data <- data[,-1]
    C <- unique(data$V2[which(data$V1 == 'C')])
    V <- unique(data$V2[which(data$V1 == 'V')])
    
    df <- data.frame(matrix(0, ncol = length(V), nrow = length(C)))
    rownames(df) <- paste('x',C,sep='')
    colnames(df) <- paste('x',V,sep='')
    
    for(i in which(data$V1=='C')){
      x <- data$V2[i]
      while(data$V1[i+1] == 'V' & i < nrow(data)){
        i <- i+1
        df[which(C==x),which(V==data[i,2])] <- 1
      }
    }
    df <- df[,order(colnames(df))]
  }
  
  transform <- function(train.d,test.d){
    test.dd <- as.data.frame(matrix(0, nrow=nrow(test.d), ncol = ncol(train.d)))
    rownames(test.dd) <- rownames(test.d)
    colnames(test.dd) <- colnames(train.d)
    
    rn <- rownames(train.d)
    cn <- colnames(train.d)
    cn.test <- colnames(test.d)
    
    res.vec <- c()
    for(i in cn.test){
      res.vec <- c(res.vec, which(i == cn))
    }
    
    for(i in 1:nrow(test.d)){
      test.dd[i, res.vec] <- test.d[i, ]
    }
    return(test.dd)
    
  }
  df_train <- sub_cleaning(train)
  test <- sub_cleaning(test)
  df_test <- transform(df_train,test)
  return(list(df_test=df_test,df_train=df_train))
}



data_cleaning_Movie <- function(data_address)
{
  DataTrain <- read.csv(paste(data_address, "/data_train.csv", sep = ""))
  DataTest <- read.csv(paste(data_address, "/data_test.csv", sep = ""))
  
  DataTrain <- DataTrain[,-1]
  DataTest <- DataTest[,-1]
  
  UserNameList <- names(table(DataTrain$User))
  UserNameNum <- length(UserNameList)
  MovieNameList <- names(table(DataTrain$Movie))
  MovieNameNum <- length(MovieNameList)
  
  DataTrain$Movie <- as.numeric(DataTrain$Movie)
  
  #transform Train Data
  TrainMatrix <- dcast(DataTrain, User ~ ...)
  MovieNameList <- names(TrainMatrix)[-1]
  UserNameList <- TrainMatrix[,1]
  TrainMatrix <- data.frame(TrainMatrix)
  TrainMatrix <- TrainMatrix[,-1]
  row.names(TrainMatrix) <- UserNameList
  colnames(TrainMatrix) <- MovieNameList
  
  save(TrainMatrix, file = "Movie_Train.RData")
  
  #transform Test Data.User:same as train data, Movie: included in train data
  TestMatrix <- dcast(DataTest, User ~ ...)
  #MovieNameList <- names(TestMatrix)[-1]
  UserNameList <- TestMatrix[,1]
  TestMatrix <- data.frame(TestMatrix)
  TestMatrix <- TestMatrix[,-1]
  row.names(TestMatrix) <- UserNameList
  colnames(TestMatrix) <- names(table(DataTest$Movie))
  
  Diff <- setdiff(MovieNameList, names(table(DataTest$Movie)))
  for (i in Diff)
    TestMatrix[,i] <- rep(NA, nrow(TestMatrix))
  
  Order <- order(as.numeric(colnames(TestMatrix)))
  TestMatrix <- TestMatrix[,Order]
  
  save(TestMatrix, file = "Movie_Test.RData")
  
  return(list(df_test = TestMatrix, df_train = TrainMatrix))
}