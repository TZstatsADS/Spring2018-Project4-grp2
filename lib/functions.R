simRank <- function(data,user1,user2,C,k){
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
  }
  
  return(simuser[user1,user2])
}


vectorsim <- function(vec1,vec2){
  if(is.null(vec1) == TRUE){
    return(0)
  }else{
    n1 <- norm(vec1, type = "2")
    n2 <- norm(vec2, type = "2")
    ip <- t(vec1) %*% vec2
    vs <- ip/(n1*n2)
    return(vs[1])
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

