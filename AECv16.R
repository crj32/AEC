AEC <- function(mydata,maxiter=80,clusteralg='spectral',mode='random',
                krange=seq(2,15),RKmin=2, RKmax=10,FASP=FALSE,kernel='density',
                stabilitysearch=TRUE,weightedgraph=FALSE,graphcluster='walktrap',
                FASPk=900,maxK=30,NN=3,NN2=5,attackmethod='random',
                usespectral=TRUE,diffusion=TRUE,diffusion_iter=3,
                specclusteralg='km',showheatmap=FALSE,
                showheatmap2=FALSE,KNN_p=10,forcek=FALSE,myk=2){
  
  # algorithm based on work by Mok et al. (2012)
  # Mok, Pik-Yin, et al. "A robust adaptive clustering analysis method for 
  # automatic identification of clusters." Pattern Recognition 45.8 (2012): 3017-3033.
  
  # maxiter=80
  # clusteralg='spectral'
  # mode='random'
  # krange=seq(2,10)
  # RKmin=2
  # RKmax=10
  # FASP=FALSE
  # kernel='density'
  # stabilitysearch=TRUE
  # weightedgraph=FALSE
  # graphcluster='walktrap'
  # FASPk=900
  # maxK=30
  # NN=3
  # NN2=10
  # attackmethod='random'
  # usespectral=TRUE
  # diffusion=FALSE
  # diffusion_iter=4
  # specclusteralg='km'
  # showheatmap=FALSE
  # showheatmap2=FALSE
  # KNN_p=10
  # forcek=FALSE
  # myk=2
  
  message('***AEC***')
  message(paste('K mode:',mode))
  message(paste('clusteralg: '),clusteralg)
  message(paste('kernel: '),kernel)
  message(paste('stability search: '),as.character(stabilitysearch))
  if (stabilitysearch){
    message(paste('attack method: '),attackmethod)
  }
  message('*********')
  
  ## error handling
  if (class(mydata) == 'list' && clusteralg == 'km'){
    stop('K-means method cannot be used with multi-view data in AEC. Stopping.')
  }
  if (FASP == TRUE && FASPk >= ncol(mydata)){
    stop('FASPk must be less than the number of columns/ samples in data. Stopping.')
  }
  if (forcek == TRUE && stabilitysearch == TRUE){
    stabilitysearch = FALSE
    usespectral = TRUE
  }
  if (stabilitysearch == TRUE && weightedgraph == TRUE){
    message('If stabilitysearch is TRUE the weightedgraph should be set to FALSE. Changed automatically.')
    weightedgraph <- FALSE
  }
  
  ##
  if (class(mydata)=='list'){
    viewn <- length(mydata)
    message(paste('detected views: '),viewn)
  }else{
    message(paste('views: '),'1')
    # if colnames == nothing, then give them names to avoid errors
    # fix for matrices
    if (is.null(colnames(mydata))){
      colnames(mydata) <- paste('col_index',seq(1,ncol(mydata)),sep='')
    }
    if (is.null(row.names(mydata))){
      row.names(mydata) <- paste('row_index',seq(1,nrow(mydata)),sep='')
    }
    ## if we are running FASP
    if (FASP){ # only for one view, centers parameter user defined
      message('calculating k centroids...')
      orig <- mydata
      cs <- kmeans(t(mydata),centers=FASPk)
      csx <- cs$centers
      cas <- cs$cluster 
      mydata <- data.frame(t(csx))
      message('done.')
    }
  }
  # ok
  
  ## spec pre-calculations
  ## if using spectral clustering pre-calculate the kernel, combine data if necessary
  if (clusteralg == 'spectral'){
    message('pre-calculations for spectral clustering...')
    ## multi-view
    if (class(mydata) == 'list'){
      ## integrate data sources
      kernellist <- list()
      for (view in seq(1,length(mydata))){
        ## calculate kernel
        message(paste('calculating kernel',view))
        if (kernel == 'stsc'){
          kerneli <- rbfkernel_b(mydata[[view]],K=NN)
        }else if (kernel == 'density'){ 
          kerneli <- CNN_kernel_mine_b(mydata[[view]],NN=NN,NN2=NN2)
        }else if (kernel == 'density2'){
          kerneli <- CNN_kernel_V2(mydata[[view]],NN=NN,NN2=NN2)
        }
        message('done.')
        colnames(kerneli) <- colnames(mydata[[1]])
        row.names(kerneli) <- colnames(mydata[[1]])
        kernellist[[view]] <- kerneli
      }
      message('combining kernels...')
      A <- Reduce('+', kernellist) # construct A by linear combination
      for (col in seq(1,ncol(A))){
        KNNs <- head(rev(sort(A[,col])),(KNN_p+1)) # find the KNNs (10 default)
        tokeep <- names(KNNs)
        A[!(names(A[,col])%in%tokeep),col] <- 0
      }
      A <- A/rowSums(A) # row normalise A
      Qt <- A
      im <- matrix(ncol=ncol(A),nrow=ncol(A))
      im[is.na(im)] <- 0
      diag(im) <- 1
      for (t in seq(1,diffusion_iter)){ # (5)
        Qt <- A%*%Qt%*%t(A)+im
      }
      A2 <- t(Qt)
      message('done.')
      dv <- 1/sqrt(rowSums(A2))
      l <- dv * A2 %*% diag(dv)
      decomp <- eigen(l)
      xi <- decomp$vectors
    }else{
      ## single view
      if (kernel == 'density'){
        A2 <- CNN_kernel_mine_b(mydata,NN=NN,NN2=NN2)
      }else if (kernel == 'stsc'){
        A2 <- rbfkernel_b(mydata,K=NN)
      }else if (kernel == 'density2'){
        A2 <- CNN_kernel_V2(mydata,NN=NN,NN2=NN2)
      }
      ## add names back on
      colnames(A2) <- colnames(mydata)
      row.names(A2) <- colnames(mydata)
      ##
      if (diffusion == TRUE){
        A <- A2
        for (col in seq(1,ncol(A))){
          KNNs <- head(rev(sort(A[,col])),(KNN_p+1)) # find the KNNs (10 default)
          tokeep <- names(KNNs)
          #print(tokeep)
          A[!(names(A[,col])%in%tokeep),col] <- 0
        }
        #print(A)
        A <- A/rowSums(A) # row normalise A
        Qt <- A
        im <- matrix(ncol=ncol(A),nrow=ncol(A))
        im[is.na(im)] <- 0
        diag(im) <- 1
        for (t in seq(1,diffusion_iter)){ # (default=4)
          Qt <- A%*%Qt%*%t(A)+im
        }
        A2 <- t(Qt)
        #l <- A
      }
      ##
      dv <- 1/sqrt(rowSums(A2))
      l <- dv * A2 %*% diag(dv)
      decomp <- eigen(l)
      xi <- decomp$vectors
      message('done.')
    }
  }
  
  ## for deciding K using stability
  if (forcek == FALSE){
    if (mode == 'conseq'){
      maxiter <- length(krange)
    }
    
    ## create judgement matrix (J)
    if (class(mydata) == 'list'){
      rowNum <- nrow(t(mydata[[1]]))
      J <- matrix(0, rowNum, rowNum)
    }else{
      rowNum <- nrow(t(mydata))
      J <- matrix(0, rowNum, rowNum)
    }
    diffiter <- maxiter-1
    message('building judgement matrix...')
    for (iii in seq(1,maxiter)){
      #message(iii)
      if (mode == 'random'){
        #print(RKmin)
        #print(RKmax)
        clus <- floor(runif(1, min=RKmin, max=RKmax))
        #print(clus)
        # 2 and 60 for single cell, 2 and 10 for complex shapes/ cancer
      }else if (mode == 'conseq'){
        clus <- krange[iii]
      }
      S <- matrix(0, rowNum, rowNum)
      ## cluster
      if (clusteralg == 'km'){ # can't use this if views > 1
        # k means clustering
        res.km <- kmeans(t(mydata), centers=clus)
        cluster <- res.km$cluster
      }else if (clusteralg == 'spectral'){
        # spectral clustering
        xin <- xi[,1:clus]
        yi <- xin/sqrt(rowSums(xin^2))
        yi[which(!is.finite(yi))] <- 0
        if (specclusteralg == 'km'){
          res.km <- suppressWarnings(kmeans(yi, centers=clus))
          #print(clus)
          cluster <- res.km$cluster
        }else{
          gmm <- ClusterR::GMM(yi, clus, verbose = F, seed_mode = "random_spread") # use random spread          
          pr <- ClusterR::predict_GMM(yi, gmm$centroids, gmm$covariance_matrices, gmm$weights)
          rm(.Random.seed, envir=globalenv())
          names(pr)[3] <- 'cluster'
          if (0 %in% pr$cluster){
            pr$cluster <- pr$cluster+1
          }
          cluster <- pr$cluster
        }
      }
      ## build connectivity matrix
      for (j in 1:clus) {
        X <- rep(0, rowNum)
        X[which(cluster == j)] <- 1
        S <- S + X %*% t(X)
      }
      ## add to judgement matrix
      J = J + S
    }
    message('done.')
    
    ## keep track of the names for cluster assignments
    rownames(J) <- rownames(t(mydata))
    colnames(J) <- rownames(t(mydata))
    J_new <- J
    ids <- colnames(J)
    xxx <- data.frame(J_new)
    
    ## graph clustering to determine K
    ## see https://igraph.org/r/doc/communities.html
    ## can just finish here w/o stability search
    if (stabilitysearch == FALSE){
      message('clustering graph without stability search...')
      if (weightedgraph == TRUE){
        G <- igraph::graph_from_adjacency_matrix(J_new, weighted = TRUE, mode = "undirected", diag = FALSE)
      }else if (weightedgraph == FALSE){
        JJ <- J_new
        JJ[JJ>0] <- 1
        G <- igraph::graph_from_adjacency_matrix(JJ, weighted = NULL, mode = "undirected", diag = FALSE)
      }
      if (graphcluster == 'walktrap'){
        done_orig <- igraph::cluster_walktrap(G)
      }else if (graphcluster == 'greedy'){
        done_orig <- igraph::cluster_fast_greedy(G)
      }
      #print(done_orig$membership) ## cluster assignments
      #print(table(done_orig$membership))
      message('done.')
    }
  }
  
  if (usespectral == TRUE && stabilitysearch == FALSE){
    ## for either forcek TRUE or stability search FALSE
    message('performing final spectral clustering...')
    if (forcek){
      optkk <- myk
      done_orig <- list()
    }else{
      res <- done_orig$membership
      optkk <- max(res)
    }
    # spectral clustering for final assignments
    xin <- xi[,1:optkk]
    yi <- xin/sqrt(rowSums(xin^2))
    yi[which(!is.finite(yi))] <- 0
    #res.km <- suppressWarnings(kmeans(yi, centers=optkk))
    #done_orig$membership <- res.km$cluster
    gmm <- ClusterR::GMM(yi, optkk, verbose = F, seed_mode = "random_spread") # use random spread          
    pr <- ClusterR::predict_GMM(yi, gmm$centroids, gmm$covariance_matrices, gmm$weights)
    names(pr)[3] <- 'cluster'
    if (0 %in% pr$cluster){
      pr$cluster <- pr$cluster+1
    }
    #cluster <- res.km$cluster
    done_orig$membership <- pr$cluster
    message('done.')
  }
  
  ## if we are doing the stability search
  if (stabilitysearch == TRUE && forcek == FALSE){
    assignments <- list()
    message('iterative clustering of degrading graph to determine K stability...')
    kv <- c()
    for (z in seq(1,diffiter)){ # for each substraction
      ## decrease matrix by 1
      if (attackmethod == 'conseq'){
        J_new[J_new>0] <- J_new[J_new>0]-1
      }else if (attackmethod == 'random'){
        vvv <- floor(runif(1, min=5, max=round(maxiter/1.25))) # 1.25
        #message(vvv)
        J_new[J_new>0] <- J_new[J_new>0]-vvv
      }
      if (weightedgraph == FALSE){
        # temp
        Jn2 <- J_new
        Jn2[Jn2>0] <- 1 ## works better with this ON
        ## make new graph from decreased matrix
        G <- igraph::graph_from_adjacency_matrix(Jn2, weighted = NULL, mode = "undirected", diag = FALSE)
      }else if (weightedgraph == TRUE){
        ## make new graph from decreased matrix
        G <- igraph::graph_from_adjacency_matrix(J_new, weighted = TRUE, mode = "undirected", diag = FALSE)
      }
      if (graphcluster == 'walktrap'){
        done <- igraph::cluster_walktrap(G)
      }else if (graphcluster == 'greedy'){
        done <- igraph::cluster_fast_greedy(G)
      }
      if (attackmethod == 'random'){
        J_new <- J
      }
      kv[z] <- max(done$membership)
      assignments[[z]] <- done$membership
      ## if we have reached maxK then stop search
      if (length(unique(done$membership)) > maxK && attackmethod == 'conseq'){ # maxk reached
        message('stopping search.')
        break
      }
    }
    message('done.')
  }
  if (clusteralg == 'km'){
    A2 <- 'not found'
  }
  
  ## results processing and plotting
  if (stabilitysearch == TRUE && forcek == FALSE){
    if (length(kv) == 1){
      stop('Error the graph clustering algorithm failed, try different parameters. Stopping.')
    }
    kvs <- sort(kv)
    removefromplot <- 0
    rr <- (length(kv)-removefromplot):length(kv)
    #df <- data.frame(iter=seq(1,length(kv)-removefromplot-1),K=kvs[-rr])
    #p1 <- ggplot(df, aes(x=iter, y=K)) + geom_point() 
    #print(p1)
    ##
    optk <- as.numeric(names(sort(table(kv),decreasing=TRUE)[1]))
    message(paste('optimal k:',optk))
    
    ## if spectral
    if (usespectral){
      message('performing final spectral clustering...')
      optkk <- optk
      # spectral clustering for final assignments
      xin <- xi[,1:optkk]
      yi <- xin/sqrt(rowSums(xin^2))
      yi[which(!is.finite(yi))] <- 0
      if (specclusteralg == 'km'){
        res.km <- suppressWarnings(kmeans(yi, centers=optkk))
        res <- res.km$cluster
        cluster <- res.km$cluster
      }else{
        gmm <- ClusterR::GMM(yi, optkk, verbose = F, seed_mode = "random_spread") # use random spread          
        pr <- ClusterR::predict_GMM(yi, gmm$centroids, gmm$covariance_matrices, gmm$weights)
        names(pr)[3] <- 'cluster'
        if (0 %in% pr$cluster){
          pr$cluster <- pr$cluster+1
        }
        res <- pr$cluster
      }
      message('done.')
    }else{
      ## get the normal assignments out
      optindex <- which(kv==optk)[1] ## selecting index 1 here
      res <- assignments[[optindex]]
    }
    
    ## do plot of stability
    kvs <- sort(kv)
    
    if (attackmethod == 'conseq'){
      ## old plot
      removefromplot <- 0
      rr <- (length(kv)-removefromplot):length(kv) # work out what this is doing
      #df <- data.frame(iter=seq(1,length(kv)-removefromplot-1),K=kvs[-rr])
      #p1 <- ggplot2::ggplot(df, ggplot2::aes(x=iter, y=K)) + ggplot2::geom_point() + 
      #  ggplot2::scale_y_continuous(breaks=seq(0,maxK,by=2))+xlab('Iteration')
      #print(p1)
      ## new plot
      #print(kvs)
      vv <- data.frame(table(kvs[-rr]))
      #print(vv)
      p1 <- ggplot2::qplot(vv$Var1,vv$Freq,group=1,geom=c("point","line"))+ggplot2::ylab('Freq')+
        ggplot2::xlab('K')+ggplot2::theme_bw()+ggplot2::geom_line(colour="skyblue",size = 1.5)++ggplot2::geom_point(colour="black")+
        ggplot2::theme(panel.grid.major = ggplot2::element_blank(), panel.grid.minor = ggplot2::element_blank())
      print(p1)
      ## output table
    }else if (attackmethod == 'random'){
      xm <- seq(2,maxK)
      #print(kvs)
      vv <- data.frame(table(kvs))
      xd <- setdiff(xm,vv$kvs)
      xm2 <- data.frame(kvs=xd,Freq=0)
      xm2$kvs <- as.factor(xm2$kvs)
      vv <- rbind(vv,xm2)
      vv$kvs <- as.numeric(as.character(vv$kvs))
      #print(vv)
      p1 <- ggplot2::qplot(vv$kvs,vv$Freq,group=1,geom=c("point","line"))+ggplot2::ylab('Freq')+
        ggplot2::xlab('K')+ggplot2::theme_bw()+ggplot2::geom_line(colour="skyblue",size = 1.5)+ggplot2::geom_point(colour="black")+
        ggplot2::theme(panel.grid.major = ggplot2::element_blank(), panel.grid.minor = ggplot2::element_blank())+
                         ggplot2::scale_x_continuous(breaks=seq(2,maxK))
      print(p1)
    }
    
    st <- as.data.frame(table(kvs))
    colnames(st) <- c('K','frequency')
    
    ## if FASP
    if (FASP){
      casn <- cas
      casn <- casn[seq_along(cas)] <- res[as.numeric(casn[seq_along(cas)])]
      names(casn) <- names(cas)
      ## get results list ready
      mylist <- list('assignments'=casn,'centroid_assignments'=res,'stability_assessment'=st,'affinity_matrix'=A2,
                     'judgement_matrix'=J,optk=max(res))
    }else{
      ## get results list ready
      mylist <- list('assignments'=res,'stability_assessment'=st,'affinity_matrix'=A2,'judgement_matrix'=J,optk=max(res))
    }
  }else if (stabilitysearch == FALSE){
    res <- done_orig$membership
    
    if (forcek == FALSE){
      ## for forcek == FALSE
      if (FASP){
        casn <- cas
        casn <- casn[seq_along(cas)] <- res[as.numeric(casn[seq_along(cas)])]
        names(casn) <- names(cas)
        ## get results
        mylist <- list('assignments'=casn,'centroid_assignments'=res,'affinity_matrix'=A2,'judgement_matrix'=J,optk=max(res))
        message(paste('optimal k:',max(res)))
      }else{
        ## get results
        mylist <- list('assignments'=res,'affinity_matrix'=A2,'judgement_matrix'=J,optk=max(res))
        message(paste('optimal k:',max(res)))
      }
      ## get results
      mylist <- list('assignments'=res,'affinity_matrix'=A2,'judgement_matrix'=J,optk=max(res))
    }else{
      ## else forcek == TRUE
      if (FASP){
        casn <- cas
        casn <- casn[seq_along(cas)] <- res[as.numeric(casn[seq_along(cas)])]
        names(casn) <- names(cas)
        ## get results
        mylist <- list('assignments'=casn,'centroid_assignments'=res,'affinity_matrix'=A2,optk=max(res))
        message(paste('optimal k:',max(res)))
      }else{
        ## get results
        mylist <- list('assignments'=res,'affinity_matrix'=A2,optk=max(res))
        message(paste('optimal k:',max(res)))
      }
      ## get results
      mylist <- list('assignments'=res,'affinity_matrix'=A2,optk=max(res))
    }
    
  }
  
  ##
  if (showheatmap){
    hmr <- cchc(J,maxiter=maxiter)
  }
  if (showheatmap2){
    hms <- smhc(A2,res)
  }
  
  message('finished.')
  ## return results
  return(mylist)
}

### fast CNN kernel
CNN_kernel_mine_b <- function(mat, NN = 3, NN2 = 7) {
  n <- ncol(mat)
  ## need N nearest neighbour distance per sample (kn), 
  ## and names of NN2 nearest neigbours (nbs)
  nbs <- list()
  dm <- Rfast::Dist(t(mat))
  dimnames(dm) <- list(colnames(mat), colnames(mat))
  kn <- c()
  ## find kth nearest neighbour for each sample 1...N
  for (i in seq_len(n)) {
    # sort the vector to retrieve the N nearest neighbour and the names of the NN2 nearest neighbours
    sortedvec <- sort.int(dm[i, ])
    # append the NNth nearest neighbour distance
    kn <- c(kn, sortedvec[NN + 1])
    # append the names of the NN2 nearest neighbours
    nbs[[i]] <- names(sortedvec[2:(NN2+1)])
    names(nbs)[[i]] <- names(sortedvec)[1]
  }
  ## make the symmetrical matrix of kth nearest neighbours distances
  sigmamatrix <- kn %o% kn
  ## calculate the kernel using the local statistics (scale and density) of each sample pair
  out <- matrix(nrow = n, ncol = n)  # don't overwrite functions
  # calculate the numerator beforehand
  upper <- -dm^2
  for (i in 2:n) {
    for (j in 1:(i - 1)) {
      # shared nearest neighbours between ith and jth sample
      cnns <- length(intersect(nbs[[i]], nbs[[j]]))
      upperval <- upper[i, j]
      # retrieve sigma
      localsigma <- sigmamatrix[i, j]
      # calculate local affinity
      out[i, j] <- exp(upperval / (localsigma * (cnns + 1)))
    }
  }
  # reflect, diag = 1
  out <- pmax(out, t(out), na.rm = TRUE)
  diag(out) <- 1
  ## return kernel
  return(out)
}

### fast CNN kernel V2
CNN_kernel_V2 <- function(mat, NN = 3, NN2 = 7) {
  n <- ncol(mat)
  ## upgraded CNN kernel requires NN distances for global sigma and K CNN's
  nbs <- list()
  dm <- Rfast::Dist(t(mat))
  dimnames(dm) <- list(colnames(mat), colnames(mat))
  kn <- c()
  ## find kth nearest neighbour for each sample 1...N
  for (i in seq_len(n)) {
    # sort the vector to retrieve the N nearest neighbour and the names of the NN2 nearest neighbours
    sortedvec <- sort.int(dm[i, ])
    # append the NNth nearest neighbour distance
    kn <- c(kn, sortedvec[NN + 1])
    # append the names of the NN2 nearest neighbours
    nbs[[i]] <- names(sortedvec[2:(NN2+1)])
    names(nbs)[[i]] <- names(sortedvec)[1]
  }
  ## make the symmetrical matrix of kth nearest neighbours distances
  globalsigma <- mean(kn)*mean(kn) # sigma^2
  ## calculate the kernel using the local statistics (scale and density) of each sample pair
  out <- matrix(nrow = n, ncol = n)  # don't overwrite functions
  # calculate the numerator beforehand
  upper <- -dm^2
  for (i in 2:n) {
    for (j in 1:(i - 1)) {
      # shared nearest neighbours between ith and jth sample
      cnns <- length(intersect(nbs[[i]], nbs[[j]]))
      upperval <- upper[i, j]
      # calculate local affinity
      out[i, j] <- exp(upperval / (globalsigma * (cnns + 1)))
    }
  }
  # reflect, diag = 1
  out <- pmax(out, t(out), na.rm = TRUE)
  diag(out) <- 1
  ## return kernel
  return(out)
}

### fast self tuning kernel
rbfkernel_b <- function (mat, K = 3, sigma = 1) { # calculate gaussian kernel with local sigma
  n <- ncol(mat)
  NN <- K # nearest neighbours (2-3)
  dm <- Rfast::Dist(t(mat))
  kn <- c() # find kth nearest neighbour for each sample 1...N
  for (i in seq_len(n)) {
    sortedvec <- as.numeric(sort.int(dm[i, ]))
    sortedvec <- sortedvec[!sortedvec == 0]
    kn <- c(kn, sortedvec[NN])
  }
  sigmamatrix <- kn %o% kn # make the symmetrical matrix of kth nearest neighbours distances
  upper <- -dm^2 # calculate the numerator beforehand
  out <- matrix(nrow = n, ncol = n)
  for (i in 2:n) {
    for (j in 1:(i - 1)) {
      lowerval <- sigmamatrix[i, j] # retrieve sigma
      upperval <- upper[i, j]
      out[i, j] <- exp(upperval / (lowerval * sigma)) # calculate local affinity
    }
  }
  # reflect, diag = 1
  out <- pmax(out, t(out), na.rm = TRUE)
  diag(out) <- 1
  ## return kernel
  return(out)
}

### display judgement matrix in heatmap
cchc <- function(j,maxiter=maxiter){
  # set up colours
  n <- 10
  seq = rev(seq(0,255,by=255/(n)))
  palRGB = cbind(seq,seq,255)
  mypal <-rgb(palRGB,maxColorValue=255)
  # plot consensus matrix heatmap, do not cluster rows and columns
  j2 <- j/maxiter # normalise matrix
  heatmap(j2, symm=FALSE, scale='none', col=mypal, cexRow=0.1, cexCol=0.1)
  #print(h)
  #return(h)
}

### display similarity matrix in heatmap
smhc <- function(W, group, fsize = 1.5) {
  ind <- sort(as.vector(group),index.return=TRUE)
  ind <- ind$ix
  diag(W) <- 0
  W <- W / rowSums(W)
  W <- W + t(W)
  cols <- c('#FFFFFF',colorRampPalette(RColorBrewer::brewer.pal(9,'Blues'))(100))
  image(1:ncol(W),1:nrow(W),W[ind,ind],col=cols,xlab = 'Samples',ylab='Samples',
        cex.axis=fsize,cex.lab=fsize)
  title(main = paste('K =',max(group)), cex.main = fsize, font.main = 1)
}
