require(survival)
require(JM)
require(Matrix)
require(statmod)
require(pracma)
require(splines)


prep_data <- function(LongData=NULL, SurvData = NULL, 
                      control_list=NULL, marker.name=NULL){
    
    ### control_list 
    ID_name = control_list$ID_name
    item_name = control_list$item_name
    value_name = control_list$value_name
    time_name = control_list$time_name
    fix_cov = control_list$fix_cov
    random_cov = control_list$random_cov
    FUN = control_list$FUN
    ran_time_ind=control_list$ran_time_ind
    surv_time_name = control_list$surv_time_name
    surv_status_name = control_list$surv_status_name
    surv_cov =  control_list$surv_cov
    n_points =  control_list$n_points
    
    ###
    data.list = list()
    para.list = list()
    
    flex_time = FUN(1)
    fix_est_name = c("intercept", fix_cov, colnames(flex_time))
    rand_est_name = c("intercept",random_cov, colnames(flex_time)[ran_time_ind])
    surv_est_name = surv_cov
    
    ## run LME to initiate the parameters in Longitudinal submodel
    ## Y_i
    uni_ID = SurvData[,ID_name]
    if(is.null(marker.name)){
        marker.name = unique(LongData[,item_name])
    }
    
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
                time_mat.tmp = FUN(data.tmp[data.tmp[,item_name]==x,time_name])
                cov.tmp = SurvData[SurvData[,ID_name]==i, fix_cov, drop=FALSE]
                cov.tmp = as.matrix(cov.tmp[rep(1,nrow(time_mat.tmp)),])
                as.matrix(cbind(1, cov.tmp, time_mat.tmp ))  
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
                time_mat.tmp = FUN(data.tmp[data.tmp[,item_name]==x,time_name])
                cov.tmp = SurvData[SurvData[,ID_name]==i, random_cov, drop=FALSE]
                cov.tmp = as.matrix(cov.tmp[rep(1,nrow(time_mat.tmp)),])
                as.matrix(cbind(1, cov.tmp, time_mat.tmp[,ran_time_ind,drop=FALSE]))
            }else{
                as.matrix(cbind(1,  FUN(data.tmp[data.tmp[,item_name]==x,time_name])[,ran_time_ind,drop=FALSE]))
            }
            
            #cbind(1,data.tmp$years[data.tmp$item==x])
        })
    })
    Z.list = do.call(rbind, Z.list)
    
    
    ## X_i(T_i),  z_i(T_i)
    X_T.list = lapply(uni_ID, function(i){
        T_i = SurvData[SurvData[,ID_name] == i,surv_time_name]
        #T_i = SurvData$ftime[SurvData$ID==i]
        data.tmp =  LongData[LongData[,ID_name]==i,]
        
        lapply(marker.name,function(x){
            vv = SurvData[SurvData[,ID_name]==i, fix_cov, drop=FALSE]
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
            vv = SurvData[SurvData[,ID_name]==i, random_cov, drop=FALSE]
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
        time_mat.tmp = FUN(t_node)
        
        lapply(marker.name,function(x){
            if(!is.null(fix_cov)){
                cov.tmp = SurvData[SurvData[,ID_name]==i, fix_cov, drop=FALSE]
                cov.tmp = as.matrix(cov.tmp[rep(1,nrow(time_mat.tmp)),])
                as.matrix(cbind(1, cov.tmp, time_mat.tmp ))
            }else{
                as.matrix(cbind(1,  time_mat.tmp ))
            } 
        })
    })
    
    X_t.list = do.call(rbind, X_t.list)
    
    Z_t.list = lapply(uni_ID, function(i){
        Ti = SurvData[SurvData[,ID_name] == i,surv_time_name]
        t_node = Gauss.point$nodes *(Ti/2) + Ti/2
        data.tmp =  LongData[LongData[,ID_name]==i,]
        time_mat.tmp = FUN(t_node)
        
        lapply(marker.name,function(x){
            
            if(!is.null(random_cov)){
                cov.tmp = SurvData[SurvData[,ID_name]==i, random_cov, drop=FALSE]
                cov.tmp = as.matrix(cov.tmp[rep(1,nrow(time_mat.tmp)),])
                as.matrix(cbind(1, cov.tmp, time_mat.tmp[,ran_time_ind,drop=FALSE] ))
            }else{
                as.matrix(cbind(1,  time_mat.tmp[,ran_time_ind,drop=FALSE] ))
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
    W = as.matrix(SurvData[,surv_cov,drop=FALSE])
    
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
    data.list[["ftime"]] = SurvData[,surv_time_name]
    data.list[["fstat"]] = SurvData[,surv_status_name]
    
    list(data.list=data.list, uni_ID=uni_ID, marker.name=marker.name,
         fix_est_name=fix_est_name,rand_est_name=rand_est_name,
         surv_est_name=surv_est_name)
}


VBJM_init <- function(LongData=NULL, SurvData = NULL, control_list=NULL){
    
    ###
    data.list = prep_data(LongData=LongData, SurvData = SurvData, 
                          control_list=control_list)
    #uni_ID = data.list$uni_ID
    marker.name = data.list$marker.name
    fix_est_name = data.list$fix_est_name
    rand_est_name = data.list$rand_est_name
    surv_est_name = data.list$surv_est_name
    data.list = data.list$data.list
    
    ## run LME to initiate the parameters in Longitudinal submodel
    para.list = list()
    beta.list = list()
    mu.list = list()
    V.list = list()
    Sigma.list = list()
    sig.vec = rep(NA, length(marker.name))
    alpha.vec = rep(NA,length(marker.name))
    
    for(i in seq_along(marker.name)){
        # i = 1
        # print(i)
        fitLME = init_LME(data.list$Y[,i], data.list$X[,i], 
                          data.list$Z[,i], 100, 1e-4)
        beta.list[[i]] = fitLME$beta
        mu.list[[i]] = fitLME$mu
        sig.vec[i] = fitLME$sig2
        Sigma.list[[i]] = fitLME$Sigma
        V.list[[i]] = fitLME$V
        alpha.vec[i] = 0
    }
    
    Sigma = as.matrix(bdiag(Sigma.list))
    
    V.list = do.call(cbind, V.list)
    V.list = lapply(1:nrow(V.list), function(i){
        as.matrix(bdiag(  V.list[i,] ))
    })
    
    mu.list = do.call(cbind, mu.list)
    
    ## initiate the parameters in Survival submodel

    fitSURV = survreg(Surv(data.list[["ftime"]],data.list[["fstat"]] ) ~ data.list[["W"]])
    
    theta = exp(fitSURV$coefficients[1])
    lambda = 1/fitSURV$scale
    gamma = -fitSURV$coefficients[2:(1+length(surv_est_name))]/fitSURV$scale
    
    ###
    para.list[["mu"]] = mu.list
    para.list[["V"]] = V.list
    para.list[["Sigma"]] = Sigma
    para.list[["sig2"]] = sig.vec
    
    para.list[["beta"]] = beta.list
    para.list[["weib"]] = c(lambda, theta)
    para.list[["gamma"]] = gamma
    para.list[["alpha"]] = alpha.vec
    
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


VBJM_test_data <- function(init_list=NULL, res=NULL, t_break=1, LongData_test=NULL, 
                           SurvData_test = NULL, control_list=NULL){
    # ### control_list 
    # ID_name = control_list$ID_name
    # item_name = control_list$item_name
    # value_name = control_list$value_name
    time_name = control_list$time_name
    # fix_cov = control_list$fix_cov
    # random_cov = control_list$random_cov
    # FUN = control_list$FUN
    # ran_time_ind=control_list$ran_time_ind
    surv_time_name = control_list$surv_time_name
    surv_status_name = control_list$surv_status_name
    # surv_cov =  control_list$surv_cov
    # n_points =  control_list$n_points
    ##
    ## SurvData = SurvData_test[SurvData_test[,surv_time_name] > t_break, ]
    SurvData = SurvData_test
    LongData = LongData_test[LongData_test[,time_name]<=t_break, ]
    SurvData[,surv_time_name] = t_break
    SurvData[,surv_status_name] = 0
    
    ## prepare  data for running the algorithm

    marker.name = init_list$marker.name

    data.list = prep_data(LongData=LongData, SurvData = SurvData, 
                          control_list=control_list,marker.name=marker.name)
    
    fix_est_name = data.list$fix_est_name
    rand_est_name = data.list$rand_est_name
    surv_est_name = data.list$surv_est_name
    uni_ID = data.list$uni_ID
    data.list = data.list$data.list
    
    ## initiate mu_i and v_i and copy other parameters 
    para.list = list()
    mu.list = list()
    V.list = list()
    mu_mat = matrix(0,nrow=length(rand_est_name), ncol=1)
    V_mat = matrix(0,nrow=length(rand_est_name), ncol=length(rand_est_name))
    diag(V_mat) = 1
    for(i in seq_along(marker.name)){
        mu.list[[i]] = lapply(seq_along(uni_ID), function(id){
            mu_mat
        })
        V.list[[i]] = lapply(seq_along(uni_ID), function(id){
            V_mat
        })
    }
    
    V.list = do.call(cbind, V.list)
    V.list = lapply(1:nrow(V.list), function(i){
        as.matrix(bdiag( V.list[i,] ))
    })
    
    mu.list = do.call(cbind, mu.list)
    
    para.list[["mu"]] = mu.list
    para.list[["V"]] = V.list
    para.list[["Sigma"]] = res$Sigma
    para.list[["sig2"]] = res$sig2
    para.list[["beta"]] = res$beta
    para.list[["weib"]] = res$weib
    para.list[["gamma"]] = res$gamma
    para.list[["alpha"]] = res$alpha
    
    list(data.list=data.list, para.list=para.list, uni_ID=uni_ID,
         marker.name=marker.name, fix_est_name=fix_est_name,
         rand_est_name=rand_est_name, surv_est_name=surv_est_name
    )
}


VBJM_pred <- function(init_list=NULL, res=NULL, t_break=1, tau=1, n_MC=1,
                      LongData_test=NULL, SurvData_test = NULL, control_list=NULL){
    
    testdata_list = VBJM_test_data(init_list=init_list, res=res, t_break=t_break,
                                   LongData_test = LongData_test, 
                                   SurvData_test = SurvData_test,
                                   control_list=control_list)
    
    predVBJM = VBJM_raneff(testdata_list$data.list,  testdata_list$para.list,  1e-5)
    
    # ### control_list 
    ID_name = control_list$ID_name
    fix_cov = control_list$fix_cov
    random_cov = control_list$random_cov
    FUN = control_list$FUN
    ran_time_ind=control_list$ran_time_ind
    surv_time_name = control_list$surv_time_name
    surv_status_name = control_list$surv_status_name
    surv_cov =  control_list$surv_cov
    
    ## calculate risk
    
    Ht <- function(tseq=NULL,bi=NULL,id=NULL){
        #print(bi)
        # print(id)
        time_mat.tmp = FUN(tseq)
        Xbeta = lapply(seq_along(testdata_list$marker.name),function(ii){
            if(!is.null(fix_cov)){
                cov.tmp = SurvData_test[SurvData_test[,ID_name]==id, fix_cov, drop=FALSE]
                cov.tmp = as.matrix(cov.tmp[rep(1,nrow(time_mat.tmp)),])
                as.matrix(cbind(1, cov.tmp, time_mat.tmp )) %*% predVBJM$beta[[ii]]
            }else{
                as.matrix(cbind(1,  time_mat.tmp )) %*% predVBJM$beta[[ii]]
            } 
        })
        
        Zbi = lapply(seq_along(testdata_list$marker.name),function(ii){
            if(!is.null(random_cov)){
                cov.tmp = SurvData_test[SurvData_test[,ID_name]==id, random_cov, drop=FALSE]
                cov.tmp = as.matrix(cov.tmp[rep(1,nrow(time_mat.tmp)),])
                as.matrix(cbind(1, cov.tmp, time_mat.tmp[,ran_time_ind,drop=FALSE] )) %*% bi[[ii]]
            }else{
                as.matrix(cbind(1,  time_mat.tmp[,ran_time_ind,drop=FALSE] )) %*% bi[[ii]]
            } 
        })
        
        XB = lapply(seq_along(testdata_list$marker.name),function(ii){
            predVBJM$alpha[ii] *(Xbeta[[ii]] + Zbi[[ii]])
        })
        
        XB = Reduce("+", XB)
        Wr = as.matrix(SurvData_test[SurvData_test[,ID_name]==id, surv_cov, drop=FALSE]) %*% predVBJM$gamma
        
        XB = XB + as.numeric(Wr)
        h0t = (lam/theta^lam)*tseq^(lam-1)
        h0t*exp(XB)
    }
    
    lam = predVBJM$weib[1]
    theta = predVBJM$weib[2]
    
    H = rep(NA, length(testdata_list$uni_ID))
    for(i in seq_along(testdata_list$uni_ID)){
        id = testdata_list$uni_ID[i]
        bi = predVBJM$mu[i,]
        H[i] = integrate(Ht,t_break,t_break+tau,bi=bi,id=id)$value
    }
    risk = 1 - exp(-H)
    
    if(n_MC > 1){
        n_randeff = length(testdata_list$rand_est_name)
        n_marker = length(testdata_list$marker.name)
        
        HH = lapply(seq_along(testdata_list$uni_ID), function(i){
            id = testdata_list$uni_ID[i]
            mu = predVBJM$mu[i,]
            V = predVBJM$V[[i]]
            b_mat = mvrnorm(n_MC, mu = do.call(c,mu), Sigma = V)
            apply(b_mat,1,function(bb){
                #print(bb)
                bi = split(bb, rep(1:n_marker,each=n_randeff))
                integrate(Ht,t_break,t_break+tau,bi=bi,id=id)$value 
            })
            
            
        })
        risk_MC = do.call(rbind, HH)
        risk_MC = 1 - exp(-risk_MC)
        #risk = cbind(risk, apply(risk_MC,1,mean),  risk_MC)
        risk = cbind(risk, apply(risk_MC,1,mean))
        colnames(risk) = c('risk','risk_MC')
    }
    
    matid = match(testdata_list$uni_ID, SurvData_test[,ID_name])
    df = cbind(SurvData_test[matid,c(ID_name,surv_time_name,surv_status_name)],risk)
    df$risk[df[,surv_time_name] <= t_break] = NA
    df
    
}
