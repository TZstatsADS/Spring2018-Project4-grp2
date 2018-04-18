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
  r.a <- apply(train, 1, mean, na.rm = T)
  r.u <- apply(train, 1 , mean, na.rm=T)
  name.c <- colnames(test)
  name.r <- rownames(test)
  for(i in 1:nrow(train)){
    w.u.i <- weight[i, top.neighbor[[i]]]
    r.u.i <- train[top.neighbor[[i]], ]
    if(length(top.neighbor) == 0){
      pred.matrix[i,] <- r.a[i]
    } else if(length(top.neighbor) == 1){
      pred.matrix[i,] <- r.a[i] + (r.u.i - r.u[top.neighbor[[i]]]) * w.u.i / sum(w.u.i, na.rm = T)
    } else{
      pred.matrix[i,] <- r.a[i] + apply((r.u.i - r.u[top.neighbor[[i]]]) * w.u.i, 2, sum, na.rm = T) / sum(w.u.i, na.rm = T)
    }
  }
  colnames(pred.matrix) <- colnames(train)
  rownames(pred.matrix) <- rownames(train)
  pred <- pred.matrix[name.r, name.c]
  return(pred)
}


rank_score <- function(pred,test){
  d <- 0.02
  rank_mat_pred <- ncol(pred)+1-t(apply(pred,1,function(x){return(rank(x,ties.method = 'first'))}))
  rank_mat_test <- ncol(test)+1-t(apply(test,1,function(x){return(rank(x,ties.method = 'first'))}))
  vec <- as.data.frame(matrix(0, nrow = nrow(test), ncol = ncol(test)))
  for(i in 1:nrow(test)){
    vec[i, which(test[i,] <= 0)] <- 0
    vec[i, which(test[i,] > 0) ] <- test[i, which(test[i,] > 0)] - d
  }
  R_a <- apply(1/(2^((rank_mat_pred-1)/4)) * vec,1,sum)
  R_a_max <- apply(1/(2^((rank_mat_test-1)/4)) * vec,1,sum)
  return(100*sum(R_a)/sum(R_a_max))
}


# spearman CV - threshold
spearman_threshold <- function(weight,train,test){
  for(i in seq(0.1, 0.8, 0.1)){
    x.neighbor <- neighbor.threshold(weight, i)
    x.pred <- prediction.ms(train, test, weight, x.neighbor[[1]])
    x.coverage <- x.neighbor[[2]]
    x.rankscore <- rank_score(x.pred, test)
    
    #save(x.neighbor, file = paste('spearman.ms.neighbor.RData', i, sep=''))
    #save(x.pred, file = paste('spearman.ms.pred.RData', i, sep=''))
    #save(x.rankscore, file = paste('spearman.ms.rankscore.RData', i, sep=''))
    
    print(x.coverage)
    print(x.rankscore)
    
  }
  
  neighbor.threshold.apply <- function(weight, threshold){
    neighborhood <- list()
    i <- 1
    neighborhood[[i]] <-  apply(weight, 1, f, threshold = threshold)
    i <- i+1
    
    return(list(neighbor.index = neighborhood, coverage = 0))
  }
}



#spearman - best n
spearman_best_n <- function(weight,train,test){
  for(i in c(110,125,150,200)){
    
    x.neighbor <- neighbor.bestn(weight, i)
    x.pred <- prediction.ms(train, test, weight, x.neighbor[[1]])
    x.coverage <- x.neighbor[[2]]
    x.rankscore <- rank_score(x.pred, test)
    
    #save(x.neighbor, file = paste('spearman.bestn.ms.neighbor.RData',i, '.RData', sep=''))
    #save(x.pred, file = paste('spearman.bestn.ms.pred',i, '.RData', sep=''))
    #save(x.rankscore, file = paste('spearman.bestn.ms.rankscore.RData',i, '.RData', sep=''))
    
    print(i)
    print(x.coverage)
    print(x.rankscore)
    
  }
}


#spearman - combined
spearman_combined <- function(weight,train,test){
  
  for(n in c(10,30, 50,75,150)){
    for(t in c(0.1, 0.3, 0.5)){
      x.neighbor <- neighbor.both(weight, n, t)
      x.pred <- prediction.ms(train, test, weight, x.neighbor[[1]])
      x.coverage <- x.neighbor[[2]]
      x.rankscore <- rank_score(x.pred, test)
      
      #save(x.neighbor, file = paste('spearman.both.ms.neighbor.RData',n,t, '.RData', sep=''))
      #save(x.pred, file = paste('spearman.both.ms.pred',n,t, '.RData', sep=''))
      #save(x.rankscore, file = paste('spearman.both.ms.rankscore.RData',n,t, '.RData', sep=''))
      
      print(n)
      print(t)
      print(x.coverage)
      print(x.rankscore)
      print(mae())
    }
  }
}



#vectorsim threshold
vectorsim_threshold <- function(weight,train,test){
  for(i in seq(0.1, 0.2, 0.1)){
    x.neighbor <- neighbor.threshold(weight, i)
    cat("done!1\n")
    x.pred <- prediction.ms(train, test, weight, x.neighbor[[1]])
    cat("done!2\n")
    #x.coverage <- x.neighbor[[2]]
    x.rankscore <- rank_score(x.pred, test)
    cat("done!3\n")
    #save(x.neighbor, file = paste('vector_sim.ms.neighbor.RData', i, sep=''))
    #save(x.pred, file = paste('vector_sim.ms.pred.RData', i, sep=''))
    #save(x.rankscore, file = paste('vector_sim.ms.rankscore.RData', i, sep=''))
    cat("done!4\n")
    
    #print(x.coverage)
    print(x.rankscore)
  }
}



#vectorsim best n
vectorsim_best_n <- function(weight,train,test){
  for(i in c(30,60,90,120,150,180)){##枚举bestN
    
    x.neighbor <- neighbor.bestn(weight, i)
    x.pred <- prediction.ms(train, test, weight, x.neighbor[[1]])
    x.coverage <- x.neighbor[[2]]
    x.rankscore <- rank_score(x.pred, test)
    
    #save(x.neighbor, file = paste('spearman.bestn.ms.neighbor.RData',i, '.RData', sep=''))
    #save(x.pred, file = paste('spearman.bestn.ms.pred',i, '.RData', sep=''))
    #save(x.rankscore, file = paste('spearman.bestn.ms.rankscore.RData',i, '.RData', sep=''))
    print(i)
    print(x.coverage)##输出coverage
    print(x.rankscore)#输出rankscore
  }
}


#vectorsim combined








#SimRank
sim_weight <- read.csv("/Users/jiangyiran/Desktop/Stat\ 2nd/ADS/Project4/simrank.csv")
qiguai <- sim_weight
sim_weight <- sim_weight[,2:4152]
sim_weight <- as.matrix(sim_weight)
dim(sim_weight)

save(sim_weight,file = "/Users/jiangyiran/Documents/GitHub/Spring2018-Project4-group-2/data/MS_simrank.RData")

# simrank CV - threshold
simrank_threshold <- function(weight,train,test){
  for(i in seq(0.01, 0.06, 0.01)){
    x.neighbor <- neighbor.threshold(weight, i)
    print(i)
    x.pred <- prediction.ms(train, test, weight, x.neighbor[[1]])
    print(i)
    x.coverage <- x.neighbor[[2]]
    print(i)
    x.rankscore <- rank_score(x.pred, test)
    #save(x.neighbor, file = paste('simrank.ms.neighbor.RData', i, sep=''))
    #save(x.pred, file = paste('simrank.ms.pred.RData', i, sep=''))
    #save(x.rankscore, file = paste('simrank.ms.rankscore.RData', i, sep=''))
    
    print(x.coverage)
    print(x.rankscore)
    
  }
  
  neighbor.threshold.apply <- function(weight, threshold){
    neighborhood <- list()
    i <- 1
    neighborhood[[i]] <-  apply(weight, 1, f, threshold = threshold)
    i <- i+1
    
    return(list(neighbor.index = neighborhood, coverage = 0))
  }
}





#simrank - best n
simrank_best_n <- function(weight,train,test){
  for(i in c(110,125,150,200)){
    print(i)
    x.neighbor <- neighbor.bestn(weight, i)
    print(i)
    x.pred <- prediction.ms(train, test, weight, x.neighbor[[1]])
    print(i)
    x.coverage <- x.neighbor[[2]]
    print(i)
    x.rankscore <- rank_score(x.pred, test)
    print(i)
    
    #save(x.neighbor, file = paste('simrank.bestn.ms.neighbor.RData',i, '.RData', sep=''))
    #save(x.pred, file = paste('simrank.bestn.ms.pred',i, '.RData', sep=''))
    #save(x.rankscore, file = paste('simrank.bestn.ms.rankscore.RData',i, '.RData', sep=''))
    
    print(i)
    print(x.coverage)
    print(x.rankscore)
    
  }
}


#simrank - combined
simrank_combined <- function(weight,train,test){
  for(n in c(10,30, 50,75,150)){
    for(t in c(0.01, 0.03, 0.05)){
      x.neighbor <- neighbor.both(weight, n, t)
      x.pred <- prediction.ms(train, test, weight, x.neighbor[[1]])
      x.coverage <- x.neighbor[[2]]
      x.rankscore <- rank_score(x.pred, test)
      
      #save(x.neighbor, file = paste('simrank.both.ms.neighbor.RData',n,t, '.RData', sep=''))
      #save(x.pred, file = paste('simrank.both.ms.pred',n,t, '.RData', sep=''))
      #save(x.rankscore, file = paste('simrank.both.ms.rankscore.RData',n,t, '.RData', sep=''))
      
      print(n)
      print(t)
      print(x.coverage)
      print(x.rankscore)
    }
  }
}




