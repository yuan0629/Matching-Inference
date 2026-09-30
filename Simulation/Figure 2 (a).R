#setwd("")

Iteration <- 500
ncores <- 4
d1 <- 100
K <- 5
d2 <- 150*K
rreal <- r <- 2
Time <- 400*K
sigma <- 0.2
nu <- 1/d2

Mmaxlist <- c(5,10,15,20,30,40,50,60,80)
mlist <- c(50,20,20,20,20,20,20,20,20)
etalist <- c(0.10,0.28,0.32,0.34,0.36,0.36,0.38,0.42,0.43)
seedlist <- 1167:1666

lambda <- 0.1
admm_rho <- 0.02
max_iterations <- 100
absolute_tolerance <- 10^(-4)
relative_tolerance <- 10^(-3)

set.seed(1124)
M <- matrix(runif(d1*d2,-20,20),d1,d2)
decomp <- svd(M)
U0 <- decomp$u[,1:rreal]%*%diag(sqrt(decomp$d[1:rreal]))
V0 <- decomp$v[,1:rreal]%*%diag(sqrt(decomp$d[1:rreal]))
Mbase <- U0%*%t(V0)
Mbase <- Mbase/max(abs(Mbase))

soft_threshold <- function(x,threshold) sign(x)*pmax(abs(x)-threshold,0)

convex_estimator <- function(code,y,bound) {
  n <- length(code)
  rows <- unlist(lapply(code,function(x) x[,1]))
  cols <- unlist(lapply(code,function(x) x[,2]))
  index <- rows+(cols-1)*d1
  sumY <- numeric(d1*d2)
  observed_sum <- rowsum(y,index,reorder=TRUE)
  sumY[as.integer(rownames(observed_sum))] <- observed_sum[,1]
  sumY <- matrix(sumY,d1,d2)
  counts <- matrix(tabulate(index,nbins=d1*d2),d1,d2)
  B <- C <- U <- V <- matrix(0,d1,d2)

  for (iter in 1:max_iterations) {
    Z <- (2*sumY/n+admm_rho*(B-U+C-V))/(2*counts/n+2*admm_rho)
    decomp <- svd(Z+U,nu=min(d1,d2),nv=min(d1,d2))
    singular_values <- soft_threshold(decomp$d,lambda/admm_rho)
    keep <- which(singular_values>0)
    Bnew <- if (length(keep)==0) matrix(0,d1,d2) else
      decomp$u[,keep,drop=FALSE]%*%diag(singular_values[keep],length(keep))%*%
      t(decomp$v[,keep,drop=FALSE])
    Cnew <- pmin(pmax(Z+V,-bound),bound)
    U <- U+Z-Bnew
    V <- V+Z-Cnew
    primal <- sqrt(sum((Z-Bnew)^2)+sum((Z-Cnew)^2))
    dual <- admm_rho*sqrt(sum((Bnew-B)^2)+sum((Cnew-C)^2))
    eps_primal <- sqrt(2*d1*d2)*absolute_tolerance+
      relative_tolerance*max(sqrt(sum(Z^2)),sqrt(sum(Bnew^2)+sum(Cnew^2)))
    eps_dual <- sqrt(2*d1*d2)*absolute_tolerance+
      relative_tolerance*admm_rho*sqrt(sum(U^2)+sum(V^2))
    B <- Bnew
    C <- Cnew
    if (primal<=eps_primal && dual<=eps_dual) break
  }
  C
}

double_enhanced_estimator <- function(code,y,bound) {
  Time1 <- floor(length(code)/2)
  code1 <- code[1:Time1]
  n1 <- sum(vapply(code1,nrow,integer(1)))
  initial <- convex_estimator(code1,y[1:n1],bound)
  decomp <- svd(initial,nu=r,nv=r)
  Uhat <- decomp$u[,1:r,drop=FALSE]
  Vhat <- decomp$v[,1:r,drop=FALSE]

  code2 <- do.call(rbind,code[(Time1+1):length(code)])
  y2 <- y[n1+seq_len(nrow(code2))]
  positions_by_row <- split(seq_len(nrow(code2)),factor(code2[,1],levels=1:d1))
  positions_by_col <- split(seq_len(nrow(code2)),factor(code2[,2],levels=1:d2))

  Utilde <- matrix(NA_real_,d1,r)
  for (i in 1:d1) {
    positions <- positions_by_row[[i]]
    Utilde[i,] <- qr.solve(Vhat[code2[positions,2],,drop=FALSE],y2[positions])
  }
  Vtilde <- matrix(NA_real_,d2,r)
  for (j in 1:d2) {
    positions <- positions_by_col[[j]]
    rows_j <- code2[positions,1]
    y_j <- y2[positions]
    Zj <- Uhat[rows_j,,drop=FALSE]
    Vtilde[j,] <- if (length(y_j)<r || qr(Zj)$rank<r)
      solve(crossprod(Zj)+10^(-8)*diag(r),crossprod(Zj,y_j)) else
      qr.solve(Zj,y_j)
  }
  decomp <- svd(Utilde,nu=r,nv=r)
  decomp$u[,1:r,drop=FALSE]%*%t(decomp$v[,1:r,drop=FALSE])%*%t(Vtilde)
}

our_estimator <- function(code,y,M,m,eta) {
  N0 <- Time/(2*m)
  M0 <- matrix(0,d1,d2)
  for (i in 1:N0) {
    index <- code[[i]]
    M0[index] <- M0[index]+y[[i]]/(N0*nu)
  }
  U <- svd(M0)$u[,1:r]
  V <- svd(M0)$v[,1:r]

  update_G <- function(U,V,batch) {
    code_batch <- do.call(rbind,code[batch])
    U_batch <- U[code_batch[,1],,drop=FALSE]
    V_batch <- V[code_batch[,2],,drop=FALSE]
    D <- matrix(0,nrow(code_batch),r^2)
    for (l in 1:r) for (k in 1:r)
      D[,(l-1)*r+k] <- U_batch[,k]*V_batch[,l]
    matrix(solve(crossprod(D))%*%
             crossprod(D,unlist(y[batch],use.names=FALSE)),r,r)
  }

  G <- update_G(U,V,(N0+1):(2*N0))
  decomp_G <- svd(G)
  LG <- decomp_G$u
  Lambda <- diag(decomp_G$d)
  RG <- decomp_G$v

  for (p in 1:(m-1)) {
    Mcurrent <- U%*%G%*%t(V)
    sumXY <- matrix(0,d1,d2)
    for (i in (N0*(2*p)+1):(N0*(2*p+1))) {
      index <- code[[i]]
      sumXY[index] <- sumXY[index]+Mcurrent[index]-y[[i]]
    }
    Um <- U%*%LG-eta/(nu*N0)*sumXY%*%V%*%RG%*%solve(Lambda)
    Vm <- V%*%RG-eta/(nu*N0)*t(sumXY)%*%U%*%LG%*%solve(Lambda)
    U <- svd(Um)$u[,1:r]
    V <- svd(Vm)$u[,1:r]
    G <- update_G(U,V,(N0*(2*p+1)+1):(N0*2*(p+1)))
    decomp_G <- svd(G)
    LG <- decomp_G$u
    Lambda <- diag(decomp_G$d)
    RG <- decomp_G$v
  }
  max(abs(U%*%G%*%t(V)-M))
}

one_replication <- function(iter) {
  set.seed(seedlist[iter])
  cols <- lapply(1:Time,function(t) sample(d2,d1,replace=FALSE))
  noise <- lapply(1:Time,function(t) rnorm(d1,0,sigma))
  code <- lapply(cols,function(x) cbind(1:d1,x))
  error_ours <- error_tang <- numeric(length(Mmaxlist))

  for (cc in seq_along(Mmaxlist)) {
    M <- Mmaxlist[cc]*Mbase
    y <- lapply(1:Time,function(t) M[code[[t]]]+noise[[t]])
    error_ours[cc] <- our_estimator(code,y,M,mlist[cc],etalist[cc])
    Tang <- double_enhanced_estimator(code,unlist(y,use.names=FALSE),Mmaxlist[cc])
    error_tang[cc] <- max(abs(Tang-M))
  }
  list(ours=error_ours,tang=error_tang)
}

cluster <- parallel::makeCluster(ncores)
parallel::clusterExport(cluster,
  c("seedlist","d1","d2","r","Time","sigma","nu","Mmaxlist","mlist",
    "etalist","Mbase","lambda","admm_rho","max_iterations",
    "absolute_tolerance","relative_tolerance","soft_threshold",
    "convex_estimator","double_enhanced_estimator","our_estimator",
    "one_replication"),envir=environment())
result <- parallel::parLapply(cluster,seq_len(Iteration),one_replication)
parallel::stopCluster(cluster)

errormax_ours <- do.call(rbind,lapply(result,function(x) x$ours))
errormax_tang <- do.call(rbind,lapply(result,function(x) x$tang))
colnames(errormax_ours) <- colnames(errormax_tang) <- paste0("Mmax=",Mmaxlist)

x <- Mmaxlist
error_ours <- colMeans(errormax_ours)
error_tang <- colMeans(errormax_tang)
yrange <- range(error_ours,error_tang)
ypad <- 0.06*diff(yrange)

cairo_pdf("estimation oto.pdf",6,4,family="Times New Roman")
par(mar=c(4.4,4.6,0.8,0.8),mgp=c(2.6,0.8,0),tcl=-0.25)
plot(x,error_ours,type="o",pch=16,lwd=1.8,col="#0072B2",xaxt="n",ylab="",
     xlab=expression(paste("||",M,"||")[max]),
     ylim=c(yrange[1]-ypad,yrange[2]+ypad),cex.lab=1.3,cex.axis=1.15)
mtext(expression(paste("||",widehat(M)^{plain(init)}-M,"||")[max]),
      side=2,line=2.9,cex=1.2,family="Cambria Math")
axis(1,at=x,cex.axis=1.15)
lines(x,error_tang,type="o",pch=17,lty=2,lwd=1.8,col="#E69F00")
legend("topleft",c("GD on Grassmannians","Tang et al."),
       col=c("#0072B2","#E69F00"),pch=c(16,17),lty=c(1,2),lwd=1.8,bty="n")
dev.off()
