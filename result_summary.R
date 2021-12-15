#------------------------Summarizing Results------------------------#
#------------------------Jieqi Tu------------------------#


# Create data frames to store the result at each iteration
n.sim = 100
VBJM_est_gamma = data.frame(
  estimate = numeric(n.sim),
  se = numeric(n.sim),
  ci_low = numeric(n.sim),
  ci_high = numeric(n.sim),
  cover = numeric(n.sim)
)

VBJM_est_alpha1 = data.frame(
  estimate = numeric(n.sim),
  se = numeric(n.sim),
  ci_low = numeric(n.sim),
  ci_high = numeric(n.sim),
  cover = numeric(n.sim)
)

VBJM_est_beta1 = data.frame(
  estimate = numeric(n.sim),
  se = numeric(n.sim),
  ci_low = numeric(n.sim),
  ci_high = numeric(n.sim)
)

VBJM_time = data.frame(
  user = numeric(n.sim),
  system = numeric(n.sim),
  elapsed = numeric(n.sim)
)

JM_est_gamma = data.frame(
  estimate = rep(NA, n.sim),
  se = rep(NA, n.sim),
  ci_low = rep(NA, n.sim),
  ci_high = rep(NA, n.sim),
  cover = rep(NA, n.sim)
)

JM_est_alpha1 = data.frame(
  estimate = rep(NA, n.sim),
  se = rep(NA, n.sim),
  ci_low = rep(NA, n.sim),
  ci_high = rep(NA, n.sim),
  cover = rep(NA, n.sim)
)

JM_est_beta1 = data.frame(
  estimate = rep(NA, n.sim),
  se = rep(NA, n.sim),
  ci_low = rep(NA, n.sim),
  ci_high = rep(NA, n.sim)
)

JM_time = data.frame(
  user = numeric(n.sim),
  system = numeric(n.sim),
  elapsed = numeric(n.sim)
)

joineRML_est_gamma = data.frame(
  estimate = rep(NA, n.sim),
  se = rep(NA, n.sim),
  ci_low = rep(NA, n.sim),
  ci_high = rep(NA, n.sim),
  cover = rep(NA, n.sim)
)

joineRML_est_alpha1 = data.frame(
  estimate = rep(NA, n.sim),
  se = rep(NA, n.sim),
  ci_low = rep(NA, n.sim),
  ci_high = rep(NA, n.sim),
  cover = rep(NA, n.sim)
)

joineRML_est_beta1 = data.frame(
  estimate = rep(NA, n.sim),
  se = rep(NA, n.sim),
  ci_low = rep(NA, n.sim),
  ci_high = rep(NA, n.sim)
)

joineRML_time = data.frame(
  user = numeric(n.sim),
  system = numeric(n.sim),
  elapsed = numeric(n.sim)
)

# Read in result files
for (i in 1:100) {
  load(file = paste0("./Results/t1_n500_0/Results/simu1_result", i, ".rda"))
  
  #------------------- Store VBJM result -------------------#
  for (j in 1:4) {
    VBJM_est_gamma[i, j] = res_VBJM['Surv_gamma_x', j]
    VBJM_est_alpha1[i, j] = res_VBJM['gene1_alpha', j]
    VBJM_est_beta1[i, j] = res_VBJM['gene1_fix_year_l', j]
  }
  if ((VBJM_est_gamma[i, 3] <= 1) && (VBJM_est_gamma[i, 4] >= 1)) {
    VBJM_est_gamma[i, 5] = 1
  } else VBJM_est_gamma[i, 5] = 0
  
  if ((VBJM_est_alpha1[i, 3] <= 1) && (VBJM_est_alpha1[i, 4] >= 1)) {
    VBJM_est_alpha1[i, 5] = 1
  } else VBJM_est_alpha1[i, 5] = 0
  
  for (k in 1:3) {
    VBJM_time[i, k] = time_VBJM[[k]]
  }
  
  
  
  #------------------- Store JM result -------------------#
  for (j in 1:2) {
    JM_est_alpha1[i, j] = res_JM$`CoefTable-Event`['Assoct', j]
    JM_est_gamma[i, j] = res_JM$`CoefTable-Event`['x', j]
    JM_est_beta1[i, j] = res_JM$`CoefTable-Long`['years', j]
  }
  
  JM_est_alpha1[i,3] = JM_est_alpha1[i,1] - qnorm(0.975)*JM_est_alpha1[i,2]
  JM_est_alpha1[i,4] = JM_est_alpha1[i,1] + qnorm(0.975)*JM_est_alpha1[i,2]
  JM_est_gamma[i,3] = JM_est_gamma[i,1] - qnorm(0.975)*JM_est_gamma[i,2]
  JM_est_gamma[i,4] = JM_est_gamma[i,1] + qnorm(0.975)*JM_est_gamma[i,2]
  JM_est_beta1[i,3] = JM_est_beta1[i,1] - qnorm(0.975)*JM_est_beta1[i,2]
  JM_est_beta1[i,4] = JM_est_beta1[i,1] + qnorm(0.975)*JM_est_beta1[i,2]
  
  if ((JM_est_gamma[i, 3] <= 1) && (JM_est_gamma[i, 4] >= 1)) {
    JM_est_gamma[i, 5] = 1
  } else JM_est_gamma[i, 5] = 0
  
  if ((JM_est_alpha1[i, 3] <= 1) && (JM_est_alpha1[i, 4] >= 1)) {
    JM_est_alpha1[i, 5] = 1
  } else JM_est_alpha1[i, 5] = 0
  
  for (k in 1:3) {
    JM_time[i, k] = time_JM[[k]]
  }
  
  #------------------- Store joineRML result -------------------#
  for (j in 1:2) {
    joineRML_est_alpha1[i, j] = res_joineRML$coefs.surv['gamma_1', j]
    joineRML_est_gamma[i, j] = res_joineRML$coefs.surv['x', j]
    joineRML_est_beta1[i, j] = res_joineRML$coefs.long['years_1', j]
  }
  
  joineRML_est_alpha1[i,3] = joineRML_est_alpha1[i,1] - qnorm(0.975)*joineRML_est_alpha1[i,2]
  joineRML_est_alpha1[i,4] = joineRML_est_alpha1[i,1] + qnorm(0.975)*joineRML_est_alpha1[i,2]
  joineRML_est_gamma[i,3] = joineRML_est_gamma[i,1] - qnorm(0.975)*joineRML_est_gamma[i,2]
  joineRML_est_gamma[i,4] = joineRML_est_gamma[i,1] + qnorm(0.975)*joineRML_est_gamma[i,2]
  joineRML_est_beta1[i,3] = joineRML_est_beta1[i,1] - qnorm(0.975)*joineRML_est_beta1[i,2]
  joineRML_est_beta1[i,4] = joineRML_est_beta1[i,1] + qnorm(0.975)*joineRML_est_beta1[i,2]
  
  # calculate coverage rate of the confidence interval
  if ((joineRML_est_gamma[i, 3] <= 1) && (joineRML_est_gamma[i, 4] >= 1)) {
    joineRML_est_gamma[i, 5] = 1
  } else joineRML_est_gamma[i, 5] = 0
  
  if ((joineRML_est_alpha1[i, 3] <= 1) && (joineRML_est_alpha1[i, 4] >= 1)) {
    joineRML_est_alpha1[i, 5] = 1
  } else joineRML_est_alpha1[i, 5] = 0
  
  for (k in 1:3) {
    joineRML_time[i, k] = time_joineRML[[k]]
  }
}

