# require(Rcpp)
# sourceCpp("VBJM.cpp")
require(survival)
require(JM)
require(Matrix)
require(statmod)
require(pracma)
require(splines)


# load("simu1.RData")
# ID_name = "ID"
# item_name = "item"
# value_name = "value"
# time_name = "years"
# fix_cov = NULL
# random_cov = NULL
# FUN = flex_time_fun
# ran_time_ind=1
# surv_time_name = "ftime"
# surv_status_name = "fstat"
# surv_cov = "x"
# n_points = 5

# knots = c(0,10)
# flex_time_fun <- function(x=NULL){
#     # xx = matrix(c(x,x^2,x^3), ncol = 3)
#     # colnames(xx) = c("year_1","year_2","year_3")
#     xx = ns(x,3, Boundary.knots=knots)[,1:3,drop=FALSE]
#     colnames(xx) = c("year_1","year_2","year_3")
#     # xx = matrix(x, ncol = 1)
#     # colnames(xx) = c("year_l")
#     xx
# }
# ran_time_ind = 1:3
# FUN = flex_time_fun
# 
# LongData2 = LongData
# LongData2$item ="gene2"
# LongData = rbind(LongData, LongData2)


VBJM_init <- function(LongData=NULL, ID_name = "ID", item_name = "item",
                      value_name = "value",  time_name = "years",  
                      fix_cov = NULL, random_cov = NULL, 
                      FUN = NULL, ran_time_ind=1,
                      SurvData = NULL,surv_time_name = "ftime", 
                      surv_status_name = "fstat", surv_cov = "x",
                      n_points = 5){
    
    data.list = list()
    para.list = list()
    
    flex_time = FUN(1)
    fix_est_name = c("intercept", fix_cov, colnames(flex_time))
    rand_est_name = c("intercept",random_cov, colnames(flex_time)[ran_time_ind])
    surv_est_name = surv_cov
    
    ## run LME to initiate the parameters in Longitudinal submodel
    ## Y_i
    uni_ID = unique(LongData[,ID_name])
    marker.name = sort(unique(LongData[,item_name]))
    
    Y.list = lapply(uni_ID, function(i){
        data.tmp = LongData[LongData[,ID_name]==i,]
        lapply(marker.name,function(x){
            matrix(data.tmp[data.tmp[,item_name]==x,value_name],ncol=1)
            # matrix(data.tmp$value[data.tmp$item==x],ncol=1)
        })
    })
    
    Y.list = do.call(rbind, Y.list)
    
    ## X_i, Z_i
    X.list = lapply(uni_ID, function(i){
        data.tmp =  LongData[LongData[,ID_name]==i,]
        lapply(marker.name,function(x){
            if(!is.null(fix_cov)){
                as.matrix(cbind(1, data.tmp[data.tmp[,item_name]==x,fix_cov],
                                FUN(data.tmp[data.tmp[,item_name]==x,time_name]) ))  
            }else{
                as.matrix(cbind(1, FUN(data.tmp[data.tmp[,item_name]==x,time_name])))
            }
            
            #cbind(1,data.tmp$years[data.tmp$item==x])
        })
    })
    X.list = do.call(rbind, X.list)
    
    
    Z.list = lapply(uni_ID, function(i){
        data.tmp =  LongData[LongData[,ID_name]==i,]
        lapply(marker.name,function(x){
            if(!is.null(random_cov)){
                as.matrix(cbind(1, data.tmp[data.tmp[,item_name]==x,random_cov],
                                FUN(data.tmp[data.tmp[,item_name]==x,time_name])[,ran_time_ind,drop=FALSE]))
            }else{
                as.matrix(cbind(1,  FUN(data.tmp[data.tmp[,item_name]==x,time_name])[,ran_time_ind,drop=FALSE]))
            }
            
            #cbind(1,data.tmp$years[data.tmp$item==x])
        })
    })
    Z.list = do.call(rbind, Z.list)
    
    beta.list = list()
    mu.list = list()
    V.list = list()
    Sigma.list = list()
    sig.vec = rep(NA, length(marker.name))
    alpha.vec = rep(NA,length(marker.name))
    
    for(i in seq_along(marker.name)){
        # i = 1
        # print(i)
        fitLME = init_LME(Y.list[,i], X.list[,i], Z.list[,i], 100, 1e-4)
        beta.list[[i]] = fitLME$beta
        mu.list[[i]] = fitLME$mu
        sig.vec[i] = fitLME$sig2
        Sigma.list[[i]] = fitLME$Sigma
        V.list[[i]] = fitLME$V
        alpha.vec[i] = 0
    }
    
    ## initiate the parameters in Survival submodel
    formula_surv = paste("Surv(", surv_time_name,",", surv_status_name,")" ,"~",
                         paste(surv_cov,collapse ="+"),
                         sep="")
    
    fitSURV = survreg(as.formula(formula_surv), data = SurvData)
    theta = exp(fitSURV$coefficients[1])
    lambda = 1/fitSURV$scale
    gamma = -fitSURV$coefficients[2]/fitSURV$scale
    
    
    Sigma = as.matrix(bdiag(Sigma.list))
    
    V.list = do.call(cbind, V.list)
    V.list = lapply(1:nrow(V.list), function(i){
        as.matrix(bdiag(  V.list[i,] ))
    })
    
    mu.list = do.call(cbind, mu.list)
    
    para.list[["mu"]] = mu.list
    para.list[["V"]] = V.list
    para.list[["Sigma"]] = Sigma
    para.list[["sig2"]] = sig.vec
    
    para.list[["beta"]] = beta.list
    para.list[["weib"]] = c(lambda, theta)
    para.list[["gamma"]] = gamma
    para.list[["alpha"]] = alpha.vec
    
    
    ## prepare extra data for running the algorithm
    
    ## X_i(T_i),  z_i(T_i)
    X_T.list = lapply(uni_ID, function(i){
        T_i = SurvData[SurvData[,ID_name] == i,surv_time_name]
        #T_i = SurvData$ftime[SurvData$ID==i]
        data.tmp =  LongData[LongData[,ID_name]==i,]
        
        lapply(marker.name,function(x){
            vv = data.tmp[data.tmp[,item_name]==x,fix_cov,drop=FALSE]
            matrix(c(1, as.numeric(vv[1,]), as.numeric(FUN(T_i))), ncol=1)
            #matrix(c(1,T_i),ncol=1)
        })
    })
    X_T.list = do.call(rbind, X_T.list)
    
    Z_T.list = lapply(uni_ID, function(i){
        T_i = SurvData[SurvData[,ID_name] == i,surv_time_name]
        #T_i = SurvData$ftime[SurvData$ID==i]
        data.tmp =  LongData[LongData[,ID_name]==i,]
        lapply(marker.name,function(x){
            vv = data.tmp[data.tmp[,item_name]==x,random_cov,drop=FALSE]
            matrix(c(1,as.numeric(vv[1,]) ,as.numeric(FUN(T_i))[ran_time_ind]), ncol=1)
            #matrix(c(1,T_i),ncol=1)
        })
    })
    Z_T.list = do.call(rbind, Z_T.list)
    
    ## X_i(t) , z_i(t) 
    ## this depends on the number of legendre Gaussian quadrature points
    Gauss.point  = gauss.quad(n_points)
    # \int_0^{T_i} f(t)dt
    # t_node = Gauss.point$nodes *(Ti/2) + Ti/2
    # w_node = Gauss.point$weights
    # Ti/2 * sum(w_node * f(t_node))
    
    X_t.list = lapply(uni_ID, function(i){
        Ti = SurvData[SurvData[,ID_name] == i,surv_time_name]
        t_node = Gauss.point$nodes *(Ti/2) + Ti/2
        data.tmp =  LongData[LongData[,ID_name]==i,]
        
        lapply(marker.name,function(x){
            if(!is.null(fix_cov)){
                vv = data.tmp[data.tmp[,item_name]==x,fix_cov,drop=FALSE]
                as.matrix(cbind(1, vv[1,], FUN(t_node) ))
            }else{
                as.matrix(cbind(1,  FUN(t_node) ))
            } 
        })
    })
    
    X_t.list = do.call(rbind, X_t.list)
    
    Z_t.list = lapply(uni_ID, function(i){
        Ti = SurvData[SurvData[,ID_name] == i,surv_time_name]
        t_node = Gauss.point$nodes *(Ti/2) + Ti/2
        data.tmp =  LongData[LongData[,ID_name]==i,]
        
        lapply(marker.name,function(x){
            
            if(!is.null(random_cov)){
                vv = data.tmp[data.tmp[,item_name]==x,random_cov,drop=FALSE]
                as.matrix(cbind(1,  vv[1,], FUN(t_node)[,ran_time_ind,drop=FALSE] ))
            }else{
                as.matrix(cbind(1,  FUN(t_node)[,ran_time_ind,drop=FALSE] ))
            }
        })
    })
    
    Z_t.list = do.call(rbind, Z_t.list)
    
    
    w_node.list = lapply(uni_ID, function(i){
        Ti = SurvData[SurvData[,ID_name] == i,surv_time_name]
        Gauss.point$weights*Ti/2
    })
    t_node.list = lapply(uni_ID, function(i){
        Ti = SurvData[SurvData[,ID_name] == i,surv_time_name]
        Gauss.point$nodes *(Ti/2) + Ti/2
    })
    
    ## covariates in survival submodel
    W = matrix(SurvData[,surv_cov],ncol=1)
    
    data.list[["Y"]] = Y.list
    data.list[["X"]] = X.list
    data.list[["X_t"]] = X_t.list
    data.list[["X_T"]] = X_T.list
    data.list[["Z"]] = Z.list
    data.list[["Z_t"]] = Z_t.list
    data.list[["Z_T"]] = Z_T.list
    data.list[["W"]] = W
    data.list[["GQ_w"]] = w_node.list
    data.list[["GQ_t"]] = t_node.list
    data.list[["ftime"]] = SurvData$ftime
    data.list[["fstat"]] = SurvData$fstat
    
    list(data.list=data.list, para.list=para.list,
         marker.name=marker.name, fix_est_name=fix_est_name,
         rand_est_name=rand_est_name, surv_est_name=surv_est_name
    )
}



VBJM_get_summary <- function(init_list=NULL, res=NULL){
    
    marker.name = init_list$marker.name
    fix_est_name = init_list$fix_est_name
    #rand_est_name = init_list$rand_est_name
    surv_est_name = init_list$surv_est_name
    
    beta.list = lapply(seq_along(marker.name), function(i){
        beta = as.numeric(res$beta[[i]])
        coef_name = paste(marker.name[i],"_fix_",fix_est_name,sep="")
        names(beta) = coef_name
        beta
    })
    
    gamma = as.numeric(res$gamma)
    names(gamma) = paste("Surv_gamma_",surv_est_name,sep="")
    
    alpha = as.numeric(res$alpha)
    names(alpha) = paste(marker.name,"_alpha", sep="")
    
    weib = as.numeric(res$weib)
    names(weib) = c("Weibull_shape","Weibull_scale")
    
    para =c(do.call(c, beta.list),  gamma, alpha, weib)
    
    #res$H = -res$H
    #diag(res$H) = diag(res$H) + 1e-6
    #cov = solve(res$H)
    cov = -pinv(res$H)
    se = round(sqrt(diag(cov)),4)[1:length(para)]
    
    res_summary = data.frame(Estimate=para, SE=se,
                             para-1.96*se, para+1.96*se)
    
    se_weib = se[c(length(se)-1, length(se))]
    se_log_weib = sqrt(se_weib^2 / weib^2)
    
    ci_weib_1 = exp(log(weib[1]) + c(-1.96, 1.96) *se_log_weib[1])
    ci_weib_2 = exp(log(weib[2]) + c(-1.96, 1.96) *se_log_weib[2])
    res_summary[c(length(se)-1, length(se)),3:4] = rbind(ci_weib_1, ci_weib_2)
    
    colnames(res_summary)[3:4] = c("95%CI_lower","95%CI_upper")
    res_summary
}

