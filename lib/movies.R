find.neighbor <- function(corr.vec, threshold = NULL,n = NULL,  method = "threshold"){
  if (method == "threshold") {
    index  <- which(corr.vec > threshold | corr.vec < -threshold)
    return(index)
  }
  else{
    if(method == "best.n"){
      sort  <- sort(corr.vec, decreasing = T, index.return = T)
      index <- head(sort$ix, n)
      return(index)
    }
    else{
      if(method == "combine"){
        sort         <- sort(abs(corr.vec), decreasing = T, index.return = T)
        breakpoint   <- which.min(sort$x > threshold) -1
        index        <- head(sort$ix[1:breakpoint], n)
        return(index)
      }
    }
  }
}

###function for prediction
weight0_1mat <- function(vec){
  ind <-ifelse(is.na(vec), 0, 1)
  return(ind)
}

###prediction function
predic.func  <- function(user.id, test.data = TestMatrix, train.data = TrainMatrix, Weight.mat = Spearman_mat,
                         arg = list(threshold = NULL,n = NULL,  method = "threshold")){
  
  active.u   <- as.numeric(train.data[paste(user.id), ])
  user.mean  <- mean(active.u, na.rm = T)
  weight.all <- as.numeric(Weight.mat[paste(user.id), ])
  
  neighbor.index <- find.neighbor(weight.all, threshold = arg$threshold, n = arg$n, method = arg$method)
  neighbor.weight<- weight.all[neighbor.index] 
  
  neighbor.mat   <- train.data[neighbor.index,]
  
  item.index     <- which(!is.na(test.data[paste(user.id), ]))
  item.list      <- colnames(test.data)[item.index]
  
  other.mean     <- as.numeric(apply(neighbor.mat, 1, mean, na.rm = T))  
  other.rate     <- neighbor.mat[, item.list]                           
  
  val.mat              <- apply((other.rate - other.mean)*neighbor.weight, 2, sum, na.rm = T ) 
  mat.01               <- apply((other.rate - other.mean)*neighbor.weight, 2, weight0_1mat )
  sum.neighbor         <- as.numeric(neighbor.weight%*%mat.01)
  value                <- val.mat/sum.neighbor
  value[which(is.na(value))] <- rep(0, length(which(is.na(value))))
  
  prediction                         <- user.mean + value
  prediction[which(prediction > 6)]  <- rep(6, length(which(prediction > 6)))  
  prediction[which(prediction < 1)]  <- rep(1, length(which(prediction < 1) ))
  prediction           <- data.frame(Moive = item.list, User = rep(user.id, length(item.list)), prediction = prediction )
  return(prediction)
}
