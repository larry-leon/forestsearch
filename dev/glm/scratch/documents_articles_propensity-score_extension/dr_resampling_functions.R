
# This is used for setting x=exp(z)/(1+exp(z))=1 for z=inf
# where x will be NA
pi.inf<-function(x){
xmod<-ifelse(is.na(x),1.0,x)
xmod
}

########################################################
# The function rescue.resample.binomial
# is used for estimating logistic model
# in resampling stage (or estimation stage 
# however we use the glm() R function) when 
# our Newton-Raphson function fails to converge.
# The function uses various optimization functions.
########################################################

rescue.resample.QLE<-function(b.start,Y,X,q.star,wt=rep(1,length(Y)),family,checking,draw,abstol=0.005){
round<-1
b.zero<-rep(0,length(b.start))
#nlm(f, p, hessian = FALSE, typsize=rep(1, length(p)), fscale=1,
#    print.level = 0, ndigit=12, gradtol = 1e-6,
#    stepmax = max(1000 * sqrt(sum((p/typsize)^2)), 1000),
#    steptol = 1e-6, iterlim = 100, check.analyticals = TRUE, ...)
#ndigit=12, gradtol = 1e-6

if(checking) cat("Going to full iteration (initial round)","\n")

logit.star<-nlm(f=QLE.star,p=c(b.start),iterlim=100,X=X,Y=Y,wt=wt,q.star=q.star,family=family,check.analyticals=FALSE)
Q1<-Q.soln<-logit.star$minimum
converged<-(logit.star$code<=2 & Q.soln<=abstol)
bstar<-logit.star$estimate
b1<-bstar

if(converged & checking){ 
cat("NLM first round converged Q(soln), convergence, converged=",c(Q.soln,logit.star$code,converged),"\n")
grad.nlm<-c(logit.star$gradient)
mine<-del.norm(b=bstar,Y=Y,X=X,q.star=q.star,wt=wt,family=family)
cat("Checking Gradient (max|nlm-mine|)",c(max(abs(grad.nlm-mine))),"\n")
#cat("Estimate=",c(bstar),"\n")
}

# If first round does not converge go to 2nd round
# Begin rounds -- loop 1

if(!converged){

round<-2.0
b.start2<-bstar
if(Q.soln>5) b.start2<-b.start

if(checking) cat("Going to second round Q(soln) round 1 =",c(Q.soln),"\n")

logit.star<-optim(par=c(b.start2),fn=QLE.star,X=X,Y=Y,wt=wt,q.star=q.star,family=family,
control=list(abstol=abstol, maxit=300), method="Nelder-Mead")

bstar<-logit.star$par
Q2<-Q.soln<-logit.star$value
converged<-(logit.star$convergence==0 & Q.soln<=abstol)
b2<-bstar
}

#logit.star<-nlminb(start=c(b.start),objective=Q.logistic.star,X=X,Y=Y,q.star=q.star,wt=wt,gradient=del.norm,
#control=list(abstol=0.005,iter.max=100))
# Evaluate Convergence
#converged<-(logit.star$convergence==0)
#bstar<-logit.star$par
#Q.soln<-logit.star$objective

# If this last guy did not converge then try Nelder-Mead -- 3rd round

if(!converged){
round<-3.0
b.start3<-bstar
if(Q.soln>=5) b.start3<-b.zero

if(checking) cat("Going to third round Q(soln) round 2 =",c(Q.soln),"\n")

logit.star<-optim(par=c(b.start3),fn=QLE.star,X=X,Y=Y,wt=wt,q.star=q.star,gr=del.norm,family=family,
control=list(abstol=abstol, maxit=400), method="BFGS")

bstar<-logit.star$par
Q3<-Q.soln<-logit.star$value
converged<-(logit.star$convergence==0 & Q.soln<=abstol)
b3<-bstar
}

# If this last guy still did not converge  -- Try different optimization routine 
# Let's lower our convergence standards (it's almost last call at the bar)

if(!converged){
round<-4.0
if(checking) cat("Warning -- Going to 4th round Q(soln) round 3 =",c(Q.soln),"\n")

#b.start4<-bstar
#if(Q.soln>=10) b.start4<-b.zero

b.start4<-(b1*(Q1<=1)+b2*(Q2<=1)+b3*(Q3<=1))/((Q1<=1)+(Q2<=1)+(Q3<=1)) # Take average of previous estimators
if(any(is.na(b.start4))) b.start4<-b.start

#logit.star<-nlm(f=QLE.star.withgrad,p=c(b.start4), X=X,Y=Y,wt=wt,q.star=q.star,check.analytical=FALSE,family=family)
#bstar<-logit.star$estimate
#Q.soln<-logit.star$minimum
#converged<-(logit.star$code<=2 & Q.soln<=abstol)
#b4<-bstar

logit.star<-nlminb(start=c(b.start4),objective=QLE.star,X=X,Y=Y,wt=wt,q.star=q.star,
gradient=del.norm,
family=family,control=list(abs.tol=abstol, iter.max=300))
bstar<-logit.star$par
Q4<-Q.soln<-logit.star$objective
converged<-(logit.star$convergence==0 & Q.soln<=abstol)
b4<-bstar
}

if(!converged){
round<-5
if(checking) cat("Warning** -- Going to 5th round Q(soln) round 4 =",c(Q.soln),"\n")

logit.star<-nlminb(start=c(b.zero),objective=QLE.star,X=X,Y=Y,wt=wt,q.star=q.star,
gradient=del.norm,
family=family,control=list(abs.tol=abstol, iter.max=300))

bstar<-logit.star$par
Q5<-Q.soln<-logit.star$objective
converged<-(logit.star$convergence==0 & Q.soln<=abstol)
b5<-bstar
} # End round 5

if(converged & checking){
cat("draw, q*, convergence,Q(soln)",c(draw,sqrt(sum(q.star^2)),converged,c(Q.soln)),"\n")
cat("Converged at round=",c(round),"\n")
}

if(!converged & checking) cat("Warning: Did not converge using multiple attempts","\n")

if(converged) return(list(converged=converged,bstar=bstar))

if(!converged){ 
b.combo<-(b1*(Q1<=1)+b2*(Q2<=1)+b3*(Q3<=1)+b4*(Q4<=1)+b5*(Q5<=1))/((Q1<=1)+(Q2<=1)+(Q3<=1)+(Q4<=1)+(Q5<=1))
if(any(is.na(b.combo))) b.combo<-rep(0,length(b.start))
return(list(converged=converged,bstar=NULL,b.combo=b.combo))
}

}


QLE.NR<-function(Y,X,b.start=rep(0,ncol(X)),q.star=0,wt=rep(1,length(Y)),
family,checking=FALSE,abstol=0.005){
n.iter<-100
beta<-b.start
iter<-0
div.crit<-1e-6
diff<-1
beta.old<-beta+1
norm.beta<-1000
n<-length(Y)
while((norm.beta>abstol & iter<n.iter) | diff>div.crit){
diff<-max(abs((beta-beta.old)/beta.old))
beta.old<-beta
xbeta<-X%*%beta

if(family=="binomial") m.b<-pi.inf(exp(xbeta)/(1+exp(xbeta)))
if(family=="gaussian") m.b<-xbeta

resid<-as.vector((Y-m.b)*wt)
temp<-X*resid # (n x p) matrix
S.beta<-apply(temp,2,mean)

Qstar.beta<-as.vector(S.beta-c(q.star))
norm.beta<-sqrt(sum(Qstar.beta^2))

if(family=="binomial") v.b<-m.b*(1-m.b)
if(family=="gaussian") v.b<-rep(1,length(Y))

info<-t(X)%*%(diag(c(v.b))*c(wt))%*%X
info<-info/n
det.info<-det(info)
if(checking) cat("iter, det, ||Q(beta)||",c(iter,det.info,norm.beta),"\n")
info.qr<-qr(info)
nc<-ncol(info.qr$qr)
if(info.qr$rank!=nc | is.na(norm.beta)){
return(NULL)
}
else{ 
beta<-beta+qr.solve(info,Qstar.beta)
iter<-iter+1
}
}
if(norm.beta<abstol & diff<div.crit){
return(list(beta=beta,norm=norm.beta,iter=iter))
}
else{
return(NULL)
}
}

rmultnorm <- function(n, mu, vmat,vsqrt=NULL,tol = 1e-07)
       {
             if(is.null(vsqrt)) p <- ncol(vmat)
             if(is.null(vmat)) p<-ncol(vsqrt)
               if(length(mu)!=p)
                       stop("mu vector is the wrong length")
                       
             if(is.null(vsqrt)){  if(max(abs(vmat - t(vmat))) > tol)
                       stop("vmat not symmetric")}
                       
               if(is.null(vsqrt)){
               vs <- svd(vmat)
               vsqrt <- t(vs$v %*% (t(vs$u) * sqrt(vs$d)))
               }
               if(all(diag(vsqrt)<=1e-8)) stop("Zero variance in Mult. Normal generation")
               ans <- matrix(rnorm(n * p), nrow = n) %*% vsqrt
               ans <- sweep(ans, 2, mu, "+")
               # dimnames(ans) <- list(NULL, dimnames(vmat)[[2]])
               ans
       }


hist.trim<-function(x,trim=0.05){
x<-na.exclude(x)
quant.x<-c(quantile(x,c(trim,1-trim)))
id.x<-(x>=quant.x[1] & x<=quant.x[2])
x<-x[id.x]
hist(x)
}



Ecdf<-function(x){
xx<-sort(unique(x))
ecdf.x<-unlist(lapply(as.list(xx),sum.toleft,Y=x))/length(x)
out<-list(x=xx,ecdf=ecdf.x)
return(out)
}


sum.toleft<-function(x,Y){
sum((Y<=x))
}

Ecdf.x<-function(x,y,trim=0.0){
y<-na.exclude(y)
y.quant<-quantile(y,c(trim,1-trim))
id.y<-(y>=y.quant[1] & y<=y.quant[2])
y<-y[id.y]
ecdf.y<-unlist(lapply(as.list(x),sum.toleft,Y=y))/length(y)
return(ecdf.y)
}

Ecdf.pair<-function(X1,X2,X3=NULL,trim1=0.0,trim2=0.0,trim3=0.0,xlabel="x",ylabel="F(x)"){
z<-na.exclude(X1)
z.quant<-quantile(z,c(trim1,1-trim1))
id<-(z>=z.quant[1] & z<=z.quant[2])
X1<-z[id]
z<-na.exclude(X2)
z.quant<-quantile(z,c(trim2,1-trim2))
id<-(z>=z.quant[1] & z<=z.quant[2])
X2<-z[id]

if(is.null(X3)){
xx<-sort(c(X1,X2))
F.1<-Ecdf.x(x=xx,y=X1)
F.2<-Ecdf.x(x=xx,y=X2)
plot(xx,F.1,xlab=xlabel,ylab=ylabel,type="l",lty=1)
lines(xx,F.2,lty=1,col="grey")
}

if(!is.null(X3)){
z<-na.exclude(X3)
z.quant<-quantile(z,c(trim3,1-trim3))
id<-(z>=z.quant[1] & z<=z.quant[2])
X3<-z[id]
xx<-sort(c(X1,X2,X3))
F.1<-Ecdf.x(x=xx,y=X1)
F.2<-Ecdf.x(x=xx,y=X2)
F.3<-Ecdf.x(x=xx,y=X3)
plot(xx,F.1,xlab=xlabel,ylab=ylabel,type="l",lty=1)
lines(xx,F.2,lty=1,col="grey")
lines(xx,F.3,lty=1,col="blue")
}

}



var.trim<-function(x,trim=0.05){
x<-na.exclude(x)
x.quant<-quantile(x,c(trim,1-trim))
id.x<-(x>=x.quant[1] & x<=x.quant[2])
x.trim<-x[id.x]
sig2<-var(x.trim)
return(sig2)
}

x.trim<-function(x,trim){
x<-na.exclude(x)
x.quant<-quantile(x,c(trim,1-trim))
id.x<-(x>=x.quant[1] & x<=x.quant[2])
x.trim<-x[id.x]
return(x.trim)
}

logit<-function(x){
log(x/(1-x))
}


logit.adjust<-function(alpha,eta,const){
prb<-exp(alpha+eta)/(1+exp(alpha+eta))
term<-mean(prb)-const
return(mean(prb)-const)
}



del.norm<-function(b,Y,X,q.star,wt=rep(1,length(Y)),family){
xbeta<-X%*%b
if(family=="binomial") m.b<-pi.inf(exp(xbeta)/(1+exp(xbeta)))
if(family=="gaussian") m.b<-xbeta
resid<-as.vector((Y-m.b)*wt)
temp<-X*resid # (n x p) matrix
S.beta<-apply(temp,2,mean)
Qstar.beta<-as.vector(S.beta-c(q.star))
g.beta<-sum(Qstar.beta^2)
if(family=="binomial") v.b<-c(m.b*(1-m.b))
if(family=="gaussian") v.b<-rep(1,length(Y))
a<--v.b*X # (n x p)
n<-length(Y)
dQ<--t(X)%*%diag(c(v.b))%*%X
dQ<-dQ/n
dQstar<-c((dQ%*%Qstar.beta)/sqrt(g.beta)) # (p x 1)

return(c(dQstar))
}

QLE.star<-function(b,Y,X,q.star,family,wt=rep(1,length(Y))){
xbeta<-X%*%b
if(family=="binomial") m.b<-pi.inf(exp(xbeta)/(1+exp(xbeta)))
if(family=="gaussian") m.b<-xbeta
resid<-as.vector((Y-m.b)*wt)
temp<-X*resid # (n x p) matrix
S.beta<-apply(temp,2,mean)
Qstar.beta<-as.vector(S.beta-c(q.star))
norm.beta<-sqrt(sum(Qstar.beta^2))
return(norm.beta)
}



QLE.star.withgrad<-function(b,Y,X,q.star,family,wt=rep(1,length(Y))){
xbeta<-X%*%b
if(family=="binomial") m.b<-pi.inf(exp(xbeta)/(1+exp(xbeta)))
if(family=="gaussian") m.b<-xbeta
resid<-as.vector((Y-m.b)*wt)
temp<-X*resid # (n x p) matrix
S.beta<-apply(temp,2,mean)
Qstar.beta<-as.vector(S.beta-c(q.star))
norm.beta<-sqrt(sum(Qstar.beta^2))
n<-length(Y)
if(family=="binomial") v.b<-c(m.b*(1-m.b))
if(family=="gaussian") v.b<-rep(1,n)
g.beta<-sum(Qstar.beta^2)
dQ<--t(X)%*%diag(c(v.b))%*%X
dQ<-dQ/n
dQstar<-c((dQ%*%Qstar.beta)/sqrt(g.beta)) # (p x 1)
attr(norm.beta,"gradient")<-c(dQstar)
return(norm.beta)
}


QLE.star2<-function(b,Y,X,q.star,family,wt=rep(1,length(Y))){
xbeta<-X%*%b
if(family=="binomial") m.b<-pi.inf(exp(xbeta)/(1+exp(xbeta)))
if(family=="gaussian") m.b<-xbeta
resid<-as.vector((Y-m.b)*wt)
temp<-X*resid # (n x p) matrix
S.beta<-apply(temp,2,mean)
S.beta<-(S.beta-c(q.star))
return(S.beta)
}


#################################################################
#                 Gehan Estimating Functions
#################################################################

simplex <- function(b,del=1,ftol=1e-6,itmax=1000,Y,X,q.star,family,wt=rep(1,length(Y))) {
  iter <- 0
  np <- length(b)
  p0 <- rep(b,np)
  dim(p0) <- c(np,np)
  p0 <- t(p0)
  diag(p0) <- diag(p0)+del
  p0 <- rbind(b,p0)
  y <- rep(0,np+1)
  for (i in 1:(np+1)){ 
y[i] <- QLE.star2(b=p0[i,],Y=Y,X=X,q.star=q.star,wt=wt,family=family)
}
  psum <- apply(p0,2,sum)
  while (iter <= itmax) {
    o <- order(y) 
    p0 <- as.matrix(p0[o,])  
    y <- y[o]
    ilo <- 1
    ihi <- np+1
    inhi <- np
    rtol <- 2*abs(y[ihi]-y[ilo])/(abs(y[ihi])+abs(y[ilo]))
    if (rtol < ftol){ 
return(list(b=p0[ilo,],value=y[ilo],comp=c(iter=iter,error=0)))
}
    if (iter >= itmax) 
{return(list(b=p0[ilo,],value=y[ilo],comp=c(iter=iter,error=1)))
}
    iter <- iter+2
# new point chosen by reflecting the worst current through the plane
# of the others

    z <- smptry(p0=p0,y=y,psum=psum,ihi=ihi,fac=-1,Y=Y,X=X,q.star=q.star,wt=wt)

    if (z[[1]] <= y[ilo]) { # new point is best--try going further

      z <- smptry(z[[4]],z[[2]],z[[3]],ihi,2,Y=Y,X=X,q.star=q.star,wt=wt)

      y <- z[[2]]; psum <- z[[3]]; p0 <- z[[4]]

    } else if (z[[1]] >= y[inhi]) {
      ysave <- z[[2]][ihi] #new point is still worst, try smaller step

      z <- smptry(p0=z[[4]],y=z[[2]],psum=z[[3]],ihi=ihi,fac=0.5,Y=Y,X=X,q.star=q.star,wt=wt)

      y <- z[[2]]; psum <- z[[3]]; p0 <- z[[4]]

      if (z[[1]] >= ysave) { # still bad, shrink simplex
    for (i in (1:(np+1))[-ilo]) {
      psum <- (p0[i,]+p0[ilo,])/2
      p0[i,] <- psum
      y[i] <- QLE.star2(b=psum,Y=Y,X=X,q.star=q.star,wt=wt,family=family)
    }
    iter <- iter+np
    psum <- apply(p0,2,sum)
      }
    } else {
      y <- z[[2]]; psum <- z[[3]]; p0 <- z[[4]]
      iter <- iter-1
    }
  }
}




smptry <- function(p0,y,psum,ihi,fac,Y,X,q.star,wt,family) {
  ndim <- ncol(p0)
  fac1 <- (1-fac)/ndim
  fac2 <- fac1-fac
  ptry <- psum*fac1-p0[ihi,]*fac2
  ytry <- QLE.star2(ptry,Y=Y,X=X,q.star=q.star,wt=wt,family=family)
  if (ytry < y[ihi]) {
    y[ihi] <- ytry
    psum <- psum-p0[ihi,]+ptry
    p0[ihi,] <- ptry
  }
  list(ytry,y,psum,p0)
}


anyMissing<-function(x){
any(is.na(x))
}

scale.tau<-function(y, center = median(y), weights = NULL, init.scale = NULL, tuning = 1.95, na.rm = TRUE)
{
    if(!missing(weights)) {
        if(length(weights) != length(y))
            stop("length of weights must be the same as length of y")
        if(anyMissing(weights))
            stop("missing values not allowed in weights")
        if(min(weights) < 0)
            stop("negative weights not allowed")
    }
    else weights <- rep(1, length(y))
    if(anyMissing(y)) {
        if(!na.rm)
            return(NA)
        nas <- is.na(y)
        y <- y[!nas]
        weights <- weights[!nas]
    }
    y <- y - center
    if(is.null(init.scale))
        init.scale <- mad(y, center = 0, low = T)
    if(init.scale <= 0)
        return(0)
    y <- y/init.scale
    if(tuning == 1.95)
        return(1.048 * init.scale * sqrt(sum(weights * pmin(y^2, 3.8025))/sum(weights)))
    #
    # compute constant for arbitrary tuning parameter
    ff <- function(x)
    pmin(x^2, tuning^2) * dnorm(x)
    xs <- seq(-10, 10, length = 1000)
    h <- xs[2] - xs[1]
    assign("tuning", tuning, frame = 1)
    consistency <- h * sum(c(17/48, 59/48, 43/48, 49/48, rep(1, length = 992), 49/48, 43/48, 59/48, 17/48) * ff(xs))
    init.scale/sqrt(consistency) * sqrt(mean(pmin(y^2, tuning^2)))
}


Qn<-function (x, constant = 2.2219, finite.corr = missing(constant)) 
{
x<-na.exclude(x)
    n <- length(x)
    if (n == 0) 
        return(NA)
    else if (n == 1) 
        return(0)
    r <- constant * .C("Qn0", as.double(x), n, res = double(1), 
        PACKAGE = "robustbase")$res
    if (finite.corr) 
        (if (n <= 9) {
            if (n == 2) 
                0.399
            else if (n == 3) 
                0.994
            else if (n == 4) 
                0.512
            else if (n == 5) 
                0.844
            else if (n == 6) 
                0.611
            else if (n == 7) 
                0.857
            else if (n == 8) 
                0.669
            else if (n == 9) 
                0.872
        }
        else {
            if (n%%2) 
                n/(n + 1.4)
            else n/(n + 3.8)
        }) * r
    else r
}

mad<-function (x, center = median(x), constant = 1.4826, na.rm = TRUE, 
    low = FALSE, high = FALSE) 
{
    if (na.rm) 
        x <- x[!is.na(x)]
    n <- length(x)
    constant * if ((low || high) && n%%2 == 0) {
        if (low && high) 
            stop("'low' and 'high' cannot be both TRUE")
        n2 <- n%/%2 + as.integer(high)
        sort(abs(x - center), partial = n2)[n2]
    }
    else median(abs(x - center))
}

scale.a<-function(y, center = median(y), weights = NULL, init.scale = NULL, tuning = 
    3.85, na.rm = TRUE)
{
    if(!missing(weights)) {
        if(length(weights) != length(y))
            stop("length of weights must be the same as length of y"
                )
        if(anyMissing(weights))
            stop("missing values not allowed in weights")
        if(min(weights) < 0)
            stop("negative weights not allowed")
    }
    else weights <- rep(1, length(y))
    if(anyMissing(y)) {
        if(!na.rm)
            return(NA)
        nas <- is.na(y)
        y <- y[!nas]
        weights <- weights[!nas]
    }
    y <- y - center
    if(is.null(init.scale))
        init.scale <- mad(y, center = 0, low = T)
    if(init.scale <= 0)
        return(0)
    y <- y/init.scale
    if(tuning == 3.85) {
        num <- ifelse(abs(y) > 3.85, 0, y * y * (14.8225 - y * y)^
            4)
        den <- ifelse(abs(y) > 3.85, 0, 219.7065 - 88.935 * y * y +
            5 * y^4)
        return((0.9471 * sqrt(sum(weights)) * init.scale * sqrt(sum(
            weights * num)))/abs(sum(weights * den)))
    }
    #
    # compute constants for arbitrary tuning parameter
    fh <- function(x)
    x * x * (tuning * tuning - x * x)^4 * dnorm(x)
    xs <- seq( - tuning, tuning, length = 600)
    h <- xs[2] - xs[1]
    assign("tuning", tuning, frame = 1)
    h.int2 <- h * sum(c(17/48, 59/48, 43/48, 49/48, rep(1, length = 592),
        49/48, 43/48, 59/48, 17/48) * fh(xs))
    fg <- function(x)
    (tuning^4 - 6 * tuning * tuning * x * x + 5 * x^4) * dnorm(x)
    g.int <- h * sum(c(17/48, 59/48, 43/48, 49/48, rep(1, length = 592),
        49/48, 43/48, 59/48, 17/48) * fg(xs))
    num <- ifelse(abs(y) > tuning, 0, y * y * (tuning * tuning - y * y)^
        4)
    den <- ifelse(abs(y) > tuning, 0, tuning^4 - 6 * tuning * tuning *
        y * y + 5 * y^4)
    (g.int/sqrt(h.int2) * sqrt(sum(weights)) * init.scale * sqrt(sum(
        weights * num)))/abs(sum(weights * den))
}
