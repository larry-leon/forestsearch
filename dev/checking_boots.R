check_boot <- subset(fs_bc$results, is.na(H_biasadj_1))

check15 <- subset(check_boot, boot_id == 15)

check15



getbadboot <- bootstrap_reproduce_aboot(this_boot = 15, fs.est = fs, cox.formula.boot = cox.formula.boot)


getbadboot[c("events_H_0","events_H_1")]

