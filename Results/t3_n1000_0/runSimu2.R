## test VBJM method ##
rm(list=ls())
#require(Rcpp)
#sourceCpp("VBJM.cpp")
library(VBJM)
library(joineRML)
source("VBJM_help.R")

args=commandArgs(trailingOnly = TRUE) 
print(args)
ffs=as.integer(args[1])
print(ffs)
set.seed(ffs)

## specify time trend
flex_time_fun <- function(x=NULL){
   # xx = matrix(c(x,x^2,x^3), ncol = 3)
    xx = ns(x,3, Boundary.knots=knots)[,1:3,drop=FALSE]
    colnames(xx) = c("year_1","year_2","year_3")
    # xx = matrix(x, ncol = 1)
    # colnames(xx) = c("year_l")
    xx
}
ran_time_ind = 1:3 ## random time-trend effects

### simu data ####
n = 1000
Ngene = 1
len = 0.2 ## time interval between visits
sig_gene = 0 # 0 indicates no correlation among genes; 0.1 indicates CS covariance matrix
source("simu_data_base.R")

knots = c(0, max(SurvData$ftime)+0.1)

### Make control list
control_list = c(
  ID_name = "ID",
  item_name = "item",
  value_name = "value",
  time_name = "years",
  fix_cov = NULL,
  random_cov = NULL,
  FUN = flex_time_fun,
  ran_time_ind = ran_time_ind,
  surv_time_name = "ftime",
  surv_status_name = "fstat",
  surv_cov =  "x",
  n_points =  5
)

#### VBJM ####

##
time_VBJM = system.time({
  
    init_list = VBJM_init(LongData = LongData, 
                        SurvData = SurvData, 
                        control_list = control_list)
    fitVBJM = VBJM(init_list$data.list,  init_list$para.list, 100, 1e-5)
})

res_VBJM = VBJM_get_summary(init_list=init_list, res=fitVBJM)
round(res_VBJM, 3)

### JM 
## pdDiag() assumes uncorrelated random effects
time_JM = system.time({
    fitLME = lme(value ~ ns(years, 3) ,  
                 random = list(ID = pdDiag(form = ~ ns(years, 3))),
                 data = LongData, control=lmeControl(opt='optim'))
    fitSURV = coxph(Surv(ftime, fstat) ~  x, data = SurvData, x = TRUE)
    fitJOINT = jointModel(fitLME, fitSURV, timeVar = "years")
})

res_JM = summary(fitJOINT)


### run mvJM
time_joineRML = system.time({
    simu.fit = 
        mjoint(
            formLongFixed = list(
                "gene1" = gene1 ~ years + years^2 + years^3
            ),
            formLongRandom = list(
                "effect1" = ~ years|ID,
                "effect2" = ~ years^2|ID,
                "effect3" = ~ years^3|ID
            ),
            formSurv = Surv(ftime, fstat) ~ x,
            data = data.mjoint,
            timeVar = "years"
        )
})

res_joineRML = summary(simu.fit)



### save results 

filename=paste("simu2_result",ffs,".rda",sep="")
save(time_VBJM,res_VBJM, res_JM,time_JM, res_joineRML, time_joineRML, file=filename)

quit(save="no")


