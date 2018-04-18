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




get_weight <- function(TrainMatrix,method){
  
  Weight.mat <- matrix(NA, nrow = nrow(TrainMatrix), ncol = nrow(TrainMatrix))
  colnames(Weight.mat) <- rownames(TrainMatrix)
  rownames(Weight.mat) <- rownames(TrainMatrix)
  
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
    for(i in 1:nrow(TrainMatrix)){
      Weight.mat[i,] <- apply(TrainMatrix, 1, vectorsim, vec1 = TrainMatrix[i,])
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
