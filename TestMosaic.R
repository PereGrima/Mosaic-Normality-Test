#########################################################
# Supplementary material that accompanies the article:  #
#               Mosaic Normality Test                   #
#########################################################

# The ADGofTest package must be installed
library(ADGofTest)

# You can change the values of the following parameters:

n=49         # n: mosaic size (nxn), with 1<n<50. If n is even, then n=n+1
color=c(0.01,0.05,0.10)  # P-value cut points for gray scale
label=F      # label: T=Shows p-value in each cell; F=don't show
bell=0       # 1: Draw the distribution in each box. 0: Do not draw the distribution

sample="NOR" # sample: SAM=from file, NOR=Random Normal(0,1), UNI=Random Uniform(0,1), EXP=Random Exponential(1)

# If sample="SAM" The sample has to be in an ASCCI file with all the values in a column without header and named sample.dat)
mu0=0        # Mean hipotized when data are from a file (sample="SAM"). For Quetelet Data: mean = 40.332 
sigma0=1     # Standard deviation hipotized when data are from a file (sample="SAM"). For Quetelet Data: sd = 2.070

# If sample !="SAM" The data are randomly generated
size=100     # Size for random samples
#set.seed(1)  # You can change the seed (or eliminate it) for random samples

test="AD"    # test: Type of test, AD=Anderson-Darling, KS=Kolmogorov-Smirnov

#########################################################
# The following functions are needed

#### Density function for the skewed power exponential distribution
dsepd=function(x,mu=0,sigma=1,alpha=0.5,p=2){
  A=2*p^(1/p)*((1-alpha)^2-alpha^2)*gamma(2/p)/gamma(1/p)
  A2=(2*p^(1/p))^2*((1-alpha)^3+alpha^3)*gamma(3/p)/gamma(1/p)
  B=sqrt(A2-A^2)
  y=A+B*(x-mu)/sigma
  kp=1/(gamma(1+1/p)*(2*p^(1/p)))
  kp*exp(-(abs(y/(2*ifelse(x<mu-sigma*A/B,alpha,1-alpha)))^p)/p)*B/sigma
}

#### Distribution function for the skewed power exponential distribution
psepd=function(x,mu=0,sigma=1,alpha=0.5,p=2){
  A=2*p^(1/p)*((1-alpha)^2-alpha^2)*gamma(2/p)/gamma(1/p)
  A2=(2*p^(1/p))^2*((1-alpha)^3+alpha^3)*gamma(3/p)/gamma(1/p)
  B=sqrt(A2-A^2)
  x=A+B*(x-mu)/sigma
  alpaux=ifelse(x<0,alpha,1-alpha)
  alpha+sign(x)*alpaux*pgamma(((abs(x)/(2*alpaux))^p)/p,shape=1/p)
}
#
sam=switch(sample,
  SAM = list(values=unlist(read.table("sample.dat",header=F)),mu0=mu0,sigma0=sigma0),
  NOR = list(values=rnorm(size,0,1),mu0=0,sigma0=1), 
  UNI = list(values=runif(size,0,1),mu0=0.5,sigma0=1/sqrt(12)), 
  EXP = list(values=rexp(size,1),mu0=1,sigma0=1)
)
#### Mosaic function
mosaic=function(sam,n,mu,sigma,label=F){
  siz=length(sam)    
  if (n<=1 | n>=50) {
    stop("n must be a value between 2 and 49\n")
  }
  pval=switch(test,
              AD = function(samp,mu=mu,sigma=sigma,alpha=alpha,p=p) ad.test(samp,distr.fun=psepd,mu=mu,sigma=sigma,alpha=alpha,p=p)$p.value,
              KS = function(samp,mu=mu,sigma=sigma,alpha=alpha,p=p) ks.test(samp,psepd,mu=mu,sigma=sigma,alpha=alpha,p=p)$p.value
  )
  if (n%%2==0) n=n+1
  cex1=c(1.2,1,0.8,0.7,0.6)
  old.par=par(mfrow = c(n,n), xaxt="n", yaxt="n", mar=c(0,0,0,0), cex=cex1[findInterval(n,c(1,5,10,15,25,9999))], xaxs="i", yaxs="i")
  s=3.5   
  ymax=c(0,1/sigma)  
  j=(log10(50)/log10(2))^(1/(n/2-0.5))
  alfa=(0:(n-1))/(n-1)
  alfa[alfa==0]=1e-12
  alfa[alfa==1]=1-1e-12
  p=c(1,2,50)
  if (n>3) p=c(1,2^(1/(j^(((n-3)/2):1))),2,2^(j^((1:(n/2-0.5)))))

  graf=function(param,label){
    x=seq(from=mu-s*sigma,to=mu+s*sigma,length.out=101)
    plot(0,0,xlim=c(mu-s*sigma,mu+s*sigma),ylim=ymax,type="n",axes=F)

    pvalue=round(pval(samp=sam,mu=mu,sigma=sigma,alpha=param[1],p=param[2]),2)
    
    id=findInterval(pvalue,c(0,color,1))
    colo=gray(c(0.9,0.85,0.4,0.2,0))[id]
    if(id>1) rect(mu-s*sigma,0,mu+s*sigma,ymax[2],col=colo)
    if (bell==1) curve(dsepd(x,mu=mu,sigma=sigma,alpha=param[1],p=param[2]),col="black",add=T)
    if (param[1]==0.5 & param[2]==2) box(lwd=2,col="white") else box()
    if (label) text(mu,ymax[2]*0.5,round(pvalue,2),col=ifelse(id>2,"white","black"))
  }
  apply(expand.grid(alfa,p),1,graf,label=label)
  par(old.par)
}
mosaic(sam$values,n,sam$mu0,sam$sigma0,label)
