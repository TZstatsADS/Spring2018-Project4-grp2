simRank <- function(data,C =0.8,k=5){
  cal_sim <- function(ui,uj,method,C){
    calsum <- 0
    if(method == 'item'){
      if(ui == uj){
        return(1)
      }else{
        if(sum(data[ui,]) == 0 | sum(data[uj,]) == 0){
          return(0)
        }
        else{
          useri <- which(data[ui,]>0)
          userj <- which(data[uj,]>0)
          for(k in useri){
            for(h in userj){
              calsum <- calsum + simitem[k,h]
            }
          }
          result <- calsum*C/((sum(data[ui,])*(sum(data[uj,]))))
        }
        
      }
    }else{
      if(ui == uj){
        return(1)
      }else{
        if(sum(data[,ui]) == 0 | sum(data[,uj]) == 0){
          return(0)
        }
        else{
          itemi <- which(data[,ui]>0)
          itemj <- which(data[,uj]>0)
          for(k in itemi){
            for(h in itemj){
              calsum <- calsum + simuser[k,h]
            }
          }
          
          result <- calsum*C/((sum(data[,ui])*(sum(data[,uj]))))
          
        }
        
      }
      
    }
    return(result)      
  }


  simuser <- diag(rep(1,nrow(data)))
  simitem <- diag(rep(1,ncol(data)))
  
  for(t in 1:k){
    newsimuser <- diag(rep(1,nrow(data)))
    newsimitem <- diag(rep(1,ncol(data)))
    for(i in 1:nrow(data)){
      for(j in 1:nrow(data)){
        newsimuser[i,j] <- cal_sim(i,j,'item',C)
      }
      
    }
    
    for(i in 1:ncol(data)){
      for(j in 1:ncol(data)){
        newsimitem[i,j] <- cal_sim(i,j,'user',C)
      }
      
    }
    simuser <- newsimuser
    simitem <- newsimitem
    print(k)
  }
  
  return(simuser)
}


vectorsim <- function(vec1,vec2){
  if(is.null(vec1) == TRUE){
    return(0)
  }else{
    inter <- intersect(which(!is.na(vec1)), which(!is.na(vec2)))
    if (length(inter) == 0) {
      return(0)
    }
    else{
      return(cosine(as.numeric(vec1[inter]),as.numeric(vec2[inter])))
    }
  }
}


spearman <- function(X){
  file=deparse(substitute(X))
  X[is.na(X)] = 0
  X = t(X)
  w = cor(X,use="everything",method="spearman")
  return(w)
}



neighbor.threshold <- function(weight, threshold){
  neighborhood <- list()
  coverage <- rep(0, ncol(weight))
  for(i in 1:nrow(weight)){
    x <- which(abs(weight[i,]) >= threshold)
    neighborhood[[i]] <- x
    coverage[x] <- 1
  }
  coverage <- length(which(coverage != 0))/ncol(weight)
  return(list(neighbor.index = neighborhood, coverage = coverage))
}


neighbor.bestn <- function(weight, n){
  neighborhood <- list()
  coverage <- rep(0, ncol(weight))
  for(i in 1:nrow(weight)){
    x <- order(weight[i,], decreasing =T)[1:n]
    neighborhood[[i]] <- x
    coverage[x] <- 1
  }
  coverage <- length(which(coverage != 0))/ncol(weight)
  return(list(neighbor.index = neighborhood, coverage = coverage))
}


neighbor.both <- function(weight, n, threshold){
  neighborhood <- list()
  coverage <- rep(0, ncol(weight))
  for(i in 1:nrow(weight)){
    y <- which(abs(weight[i,]) >= threshold)
    if(length(y) >=n){
      x <- order(weight[i,y], decreasing = T)[1:n]
    }else{
      x <- which(abs(weight[i,]) >= threshold) 
    }
    neighborhood[[i]] <- x
    coverage[x] <- 1
  }
  coverage <- length(which(coverage != 0))/ncol(weight)
  return(list(neighbor.index = neighborhood, coverage = coverage))
}

prediction.ms <- function(train, test, weight, top.neighbor){
  pred.matrix <- matrix(0, nrow = nrow(train), ncol = ncol(train))
  r.a <- apply(test, 1, mean, na.rm = T)
  r.u <- apply(train, 1 , mean, na.rm=T)
  test.col <- colnames(test)
  test.row <- rownames(test)
  for(i in 1:nrow(train)){
    w.u.i <- weight[i, top.neighbor[[i]]]
    r.u.i <- train[top.neighbor[[i]], ]
    if(length(top.neighbor) == 0){
      pred.matrix[i,] <- r.a[i]
    } else if(length(top.neighbor) == 1){
      pred.matrix[i,] <- r.a[i] + (r.u.i - r.u[i]) * w.u.i / sum(w.u.i, na.rm = T)
    } else{
      pred.matrix[i,] <- r.a[i] + apply((r.u.i - r.u[i]) * w.u.i, 2, sum, na.rm = T) / sum(w.u.i, na.rm = T)
    }
  }
  colnames(pred.matrix) <- colnames(train)
  rownames(pred.matrix) <- rownames(train)
  pred <- pred.matrix[test.row, test.col]
  return(pred)
}


get_weight <- function(TrainMatrix,TestMatrix,method){
  
  Weight.mat <- matrix(NA, nrow = nrow(TestMatrix), ncol = nrow(TrainMatrix))
  colnames(Weight.mat) <- rownames(TrainMatrix)
  rownames(Weight.mat) <- rownames(TestMatrix)
  
  if(method=='spearman'){
    for (i in 1:nrow(TestMatrix)){
      Weight.mat[i,] <- as.numeric(apply(TrainMatrix, 1, cor, y = as.numeric(TrainMatrix[i,]), method = "spearman", use = "pairwise.complete.obs"))
    }
    replace.na  <- function(vec){
      index <- which(is.na(vec))
      vec[index] <- 0
      return(vec)
    }
    Weight.mat <- apply(Weight.mat, 1, replace.na)
    return(Weight.mat)
    
  }

  if(method=='vectorsim'){
    for (i in 1:nrow(TestMatrix)){
      Weight.mat[i,] <- sapply(TrainMatrix,1,vectorsim,vec2=TextMatrix[i,])
    }

    return(Weight.mat)
    }

  
  
  if(method=='simrank'){
    replace.na  <- function(vec){
      index <- which(is.na(vec))
      vec[index] <- 0
      return(vec)
    }
    
    for(i in 1:nrow(TrainMatrix)){
      TrainMatrix[i,] <- replace.na(TrainMatrix[i,])
    }
    
    Weight.mat <- simRank(TrainMatrix,C =0.8,k=3)
    return(Weight.mat)
    
  }
}
