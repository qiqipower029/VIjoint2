## test VBJM method ##
rm(list=ls())
#require(Rcpp)
#sourceCpp("VBJM.cpp")
library(VBJM)
source("VBJM_help.R")

args=commandArgs(trailingOnly = TRUE) 
print(args)
ffs=as.integer(args[1])
#ffs = 1
print(ffs)
set.seed(ffs)

## specify time trend
flex_time_fun <- function(x=NULL){
    #xx = matrix(c(x,x^2,x^3), ncol = 3)
    #colnames(xx) = c("year_1","year_2","year_3")
    xx = matrix(x, ncol = 1)
    colnames(xx) = c("year_l")
    xx
}
ran_time_ind = 1 ## random time-trend effects

### simulate data ####
n = 2000
Ngene = 1
len = 0.2 ## time interval between visits
sig_gene = 0 # 0 indicates no correlation among genes; 0.1 indicates CS covariance matrix
source("simu_data_base.R")


### run VBJM
time_VBJM = system.time({
    init_list = VBJM_init(LongData=LongData, ID_name = "ID", item_name = "item",
                          value_name = "value",  time_name = "years",  
                          fix_cov = NULL, random_cov = NULL, 
                          FUN = flex_time_fun, ran_time_ind=ran_time_ind,
                          SurvData = SurvData, surv_time_name = "ftime", 
                          surv_status_name = "fstat", surv_cov = "x",
                          n_points = 5)
    fitVBJM = VBJM(init_list$data.list,  init_list$para.list, 100, 1e-5)
})

res_VBJM = VBJM_get_summary(init_list=init_list, res=fitVBJM)
round(res_VBJM,3)

### run JM

time_JM = system.time({
    fitLME = lme(value ~ years ,  random = ~ years | ID,
                 data = LongData, control=lmeControl(opt='optim'))
    fitSURV = coxph(Surv(ftime, fstat) ~  x, data = SurvData, x = TRUE)
    fitJOINT = jointModel(fitLME, fitSURV, timeVar = "years")
})


res_JM = summary(fitJOINT)

### run mvJM




### save results 

filename=paste("results/simu1_result",ffs,".rda",sep="")
save(time_VBJM, res_VBJM, res_JM, time_JM, file=filename)

quit(save="no")

