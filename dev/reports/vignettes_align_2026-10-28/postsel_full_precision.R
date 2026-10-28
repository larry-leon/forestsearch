# §2 record (not a gate): re-run the post-selection vignette's data + fit chunks
# (purled from the committed-to-be qmd; installed build) and print the S3-comparable
# quantities at full precision, beside stage 0's S3 row (ladder_quantities.csv).
# usage: Rscript postsel_full_precision.R <survival-post-selection.qmd> <ladder_quantities.csv> <outdir>
a <- commandArgs(TRUE)
code <- knitr::purl(a[1], output = file.path(a[3], "postsel_purled.R"), quiet = TRUE)
ex <- new.env()
lines <- readLines(code)
stop_at <- grep("^## ----subgroup", lines)[1]
eval(parse(text = lines[seq_len(stop_at - 1L)]), envir = ex)
fs <- ex$fs; g <- fs$mr_inference; fb <- g$field; fc <- fb$complement; js <- fb$joint_s
now <- c(corrected_H = g$debiased$est, field_lower_H = fb$lower_1s,
         field_s_upper_Hc = fc$upper_1s_s, bonf_lower_H = js$bonf_lower_H,
         bonf_upper_Hc = js$bonf_upper_Hc, naive_H = g$naive$est,
         corrected_Hc = g$complement$debiased$est)
s3 <- read.csv(a[2], stringsAsFactors = FALSE); s3 <- s3[s3$rung == "S3", ]
ref <- sapply(names(now), function(k) as.numeric(s3[[k]]))
cat("region:", paste(fs$sg.harm, collapse = " & "), " | S3:", s3$region, "\n")
cat("draws:", g$settings$draws, " seed:", g$settings$seed, " n_family:", g$n_family, " (S3", s3$n_family, ")\n")
cat("field R_out / R_in:", fb$R_out, "/", fb$R_in, " (S3", s3$field_R_out, "/", s3$field_R_in, ")\n")
out <- data.frame(quantity = names(now), this_fit = sprintf("%.6f", now),
                  S3 = sprintf("%.6f", ref), diff = sprintf("%+.6f", now - ref))
print(out, row.names = FALSE)
write.csv(out, file.path(a[3], "postsel_vs_S3.csv"), row.names = FALSE)
