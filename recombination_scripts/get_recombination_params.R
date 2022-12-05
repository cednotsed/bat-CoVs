rm(list = ls())
setwd("../Desktop/git_repos/bat-CoVs/")
require(tidyverse)
require(data.table)

base_dir <- "results/recombination_out"
em_files <- list.files(base_dir, "em.txt")

for (em_file in em_files) {
  prefix <- gsub(".em.txt", "", em_file)
  emsim_file <- gsub(".em.", ".emsim.", em_file)
  emsim <- fread(str_glue("{base_dir}/{emsim_file}")) %>%
    dplyr::rename(rtheta = "R/theta") %>%
    mutate(r_m = rtheta * delta * nu) 
  
  em <- fread(str_glue("{base_dir}/{em_file}"))[1:3, ]
  colnames(em) <- c("param", "post_mean", "post_var", "shape", "rate")
  
  # Calculate r/m using emsim = 100
  r_m_mean <- em[param == "R/theta"]$post_mean * 
    1 / em[param == "1/delta"]$post_mean * 
    em[param == "nu"]$post_mean
  
  r_m_mean <- round(r_m_mean, 2)
  
  r_m_CI <- quantile(emsim$r_m, probs = c(0.025,0.975))
  r_m_CI <- round(r_m_CI, 2)
  r_m_CI <- paste0(r_m_CI, collapse = "-")
  
  
  rtheta_CI <- qgamma(c(0.025,0.975),
                      shape = em[param == "R/theta"]$shape,
                      rate = em[param == "R/theta"]$rate)
  rtheta_CI <- round(rtheta_CI, 3)[order(rtheta_CI)]
  rtheta_CI <- paste0(rtheta_CI, collapse = "-")
  rtheta_mean <- em[param == "R/theta"]$post_mean
  rtheta_mean <- round(rtheta_mean, 3)
  
  
  delta_CI <- 1 / qgamma(c(0.025,0.975),
                         shape = em[param == "1/delta"]$shape,
                         rate = em[param == "1/delta"]$rate)
  delta_CI <- round(delta_CI, 0)[order(delta_CI)]
  delta_CI <- paste0(delta_CI, collapse = "-")
  delta_mean <- 1 / em[param == "1/delta"]$post_mean
  delta_mean <- round(delta_mean, 3)
  
  
  nu_CI <- qgamma(c(0.025,0.975),
                  shape = em[param == "nu"]$shape,
                  rate = em[param == "nu"]$rate)
  nu_CI <- round(nu_CI, 3)[order(nu_CI)]
  nu_CI <- paste0(nu_CI, collapse = "-")
  nu_mean <- em[param == "nu"]$post_mean
  nu_mean <- round(nu_mean, 3)
  
  res <- str_glue("Final results:\nR_theta = {rtheta_mean} [{rtheta_CI}]\ndelta = {delta_mean} [{delta_CI}]\nnu = {nu_mean} [{nu_CI}]")
  bs_res <- str_glue("Boostrapped r/m = {r_m_mean} [{r_m_CI}]")
  fileConn <- file(str_glue("{base_dir}/{prefix}_params.txt"))
  writeLines(c(res, bs_res), fileConn)
  close(fileConn)
  

}
  


