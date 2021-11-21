## simulate data ###

sig = 0.1
sig_a = 0.3
sig_b = 0.1

Ninfo = Ngene
eff = 1/Ninfo

lam = 3
theta = 10

AB1.list = lapply(1:n,function(i){
    matrix(cbind(rnorm(Ngene,0,sig_a) + rnorm(Ngene,0,sig_gene), rnorm(Ngene,0,sig_b)), nrow=Ngene)
})

## survival data
Y = rep(50,n)
x = runif(n,2,6)
beta = 1

for(i in 1:n){
    #print(i)
    AB1 = AB1.list[[i]][1:Ninfo,,drop=FALSE]
    Atrue = sum(AB1[,1])
    Btrue = sum(AB1[,2])
    
    logr = log(runif(1))
    ff <- function(t){(lam/theta^lam)*t^(lam-1)*exp(eff*(Atrue + Btrue*t) + beta*x[i])}
    f = function(x){ integrate(ff,0,x)$value + logr }
    tryCatch({Y[i] = uniroot(f,c(0,50))$root},error=function(e) e )
}

C = rweibull(n,2,3)
U = pmin(Y,C) 
Delta = as.numeric(Y<C)
table(Delta)
data.id = data.frame(ID=1:n,fstat=Delta,ftime=U, x=x)


## longitudinal data
data = lapply(1:n,function(i){
    
    ss=seq(0,data.id$ftime[i],by=len)
    years = rep(ss,Ngene)
    item = rep(paste("gene",1:Ngene,sep=""),each=length(ss))
    ABtmp = AB1.list[[i]][rep(1:Ngene,each=length(ss)),,drop=FALSE]
    
    data.frame(ID=i,item=item,years=years,AB=ABtmp)
})

data = do.call(rbind,data)
value = data[,'AB.1'] + data[,'AB.2']*data[,'years'] + rnorm(nrow(data),0,sig)
data = cbind(data,value)

SS = data.id[match(data$ID,data.id$ID),'ftime']
data = data[data$years<SS,]

LongData = data
SurvData = data.id
