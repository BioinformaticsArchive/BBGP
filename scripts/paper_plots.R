ind_all_common=getCommonIndices(c("n1_CMH_r2","n1_CMH_r3","n1_CMH_r4","n1_CMH_r5"),c(1424684,1560423,1631971,1692933),"n1_BBGP_t6_r4",1939941)
source("getPRcurve.R")

# Figure S2:

getPRcurve(c("n1_GP_t6_r2","n1_CMH_r2"),c(1939941,1424684),"n1_BBGP_t6_r2",1939941,"true_n1",100,c("BBGP","GP","CMH"),c("black","blue","red"),c(2,1.5,1),"colored_fv_nofv_PR_curve_t6_r2",ind_all_common)
getPRcurve(c("n1_GP_t6_r3","n1_CMH_r3"),c(1939941,1560423),"n1_BBGP_t6_r3",1939941,"true_n1",100,c("BBGP","GP","CMH"),c("black","blue","red"),c(2,1.5,1),"colored_fv_nofv_PR_curve_t6_r3",ind_all_common)
getPRcurve(c("n1_GP_t6_r4","n1_CMH_r4"),c(1939941,1631971),"n1_BBGP_t6_r4",1939941,"true_n1",100,c("BBGP","GP","CMH"),c("black","blue","red"),c(2,1.5,1),"colored_fv_nofv_PR_curve_t6_r4",ind_all_common)
getPRcurve(c("n1_GP_t6_r5","n1_CMH_r5"),c(1939941,1692933),"n1_BBGP_t6_r5",1939941,"true_n1",100,c("BBGP","GP","CMH"),c("black","blue","red"),c(2,1.5,1),"colored_fv_nofv_PR_curve_t6_r5",ind_all_common)


# Figure 5:

#getPRcurve(c("n1_CMH_r5","n1_CMH_r4","n1_CMH_r3","n1_CMH_r2"),c(1692933,1631971,1560423,1424684),"n1_BBGP_t3_r5",1939941,"true_n1",100,c("redundant","CMH, r=5","CMH, r=4","CMH, r=3","CMH, r=2"),c("black","black","black","black","black"),c(2,1.5,1,0.5),"uncolored_CMH_reps_PR_curve_t3",ind_all_common)
getPRcurve(c("n1_BBGP_t3_r4","n1_BBGP_t3_r3","n1_BBGP_t3_r2"),c(1939941,1939941,1939941),"n1_BBGP_t3_r5",1939941,"true_n1",100,c("BBGP, r=5","BBGP, r=4","BBGP, r=3","BBGP, r=2"),c("black","black","black","black"),c(2,1.5,1,0.5),"uncolored_BBGP_reps_PR_curve_t3",ind_all_common)
getPRcurve(c("n1_BBGP_t6_r4","n1_BBGP_t6_r3","n1_BBGP_t6_r2"),c(1939941,1939941,1939941),"n1_BBGP_t6_r5",1939941,"true_n1",100,c("BBGP, r=5","BBGP, r=4","BBGP, r=3","BBGP, r=2"),c("black","black","black","black"),c(2,1.5,1,0.5),"uncolored_BBGP_reps_PR_curve_t6",ind_all_common)
getPRcurve(c("n1_BBGP_t9_r4","n1_BBGP_t9_r3","n1_BBGP_t9_r2"),c(1939941,1939941,1939941),"n1_BBGP_t9_r5",1939941,"true_n1",100,c("BBGP, r=5","BBGP, r=4","BBGP, r=3","BBGP, r=2"),c("black","black","black","black"),c(2,1.5,1,0.5),"uncolored_BBGP_reps_PR_curve_t9",ind_all_common)
getPRcurve(c("n1_GP_t6_r4","n1_GP_t6_r3","n1_GP_t6_r2"),c(1939941,1939941,1939941),"n1_GP_t6_r5",1939941,"true_n1",100,c("GP, r=5","GP, r=4","GP, r=3","GP, r=2"),c("black","black","black","black"),c(2,1.5,1,0.5),"uncolored_GP_reps_PR_curve_t6",ind_all_common)

# Figure S3:

ind_all_common=getCommonIndices(c("n1_CMH_r4"),c(1631971),"n1_BBGP_t9_r4",1939941)
getPRcurve(c("n1_BBGP_t6_r4","n1_BBGP_t3_r4","n1_CMH_r4"),c(1939941,1939941,1631971),"n1_BBGP_t9_r4",1939941,"true_n1",100,c("BBGP, t=9","BBGP, t=6","BBGP, t=3","CMH"),c("black","blue","green","red"),c(2,1.5,1,0.5),"colored_PR_curve_exp1",ind_all_common)

ind_all_common=getCommonIndices(c("n2_CMH_r4"),c(1639083),"n2_BBGP_t9_r4",1939941)
getPRcurve(c("n2_BBGP_t6_r4","n2_BBGP_t3_r4","n2_CMH_r4"),c(1939941,1939941,1639083),"n2_BBGP_t9_r4",1939941,"true_n2",100,c("BBGP, t=9","BBGP, t=6","BBGP, t=3","CMH"),c("black","blue","green","red"),c(2,1.5,1,0.5),"colored_PR_curve_exp2",ind_all_common)

ind_all_common=getCommonIndices(c("n3_CMH_r4"),c(1638970),"n3_BBGP_t9_r4",1939941)
getPRcurve(c("n3_BBGP_t6_r4","n3_BBGP_t3_r4","n3_CMH_r4"),c(1939941,1939941,1638970),"n3_BBGP_t9_r4",1939941,"true_n3",100,c("BBGP, t=9","BBGP, t=6","BBGP, t=3","CMH"),c("black","blue","green","red"),c(2,1.5,1,0.5),"colored_PR_curve_exp3",ind_all_common)
