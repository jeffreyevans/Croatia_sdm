setwd("Z:/croatia/birds") 
f <- list.files(getwd(), "txt$", recursive = TRUE, full.names = TRUE)
v <- data.frame(species=unlist(lapply(strsplit(f,"/"), \(x) x[4])), 
                pcc=NA, kappa=NA, log_loss=NA, 
                brier=NA, auc=NA, p_thresh=NA)
  for(i in 1:length(f)){
    d <- readLines(f[i])
    n <- c(  
      as.numeric(strsplit(d[grepl("Percent Correctly Classified: ", d)], ":")[[1]][2]),
      as.numeric(strsplit(d[grepl("Kappa: ", d)], ":")[[1]][2]),
      as.numeric(strsplit(d[grepl("Log loss: ", d)], ":")[[1]][2]),
      as.numeric(strsplit(d[grepl("Brier score: ", d)], ":")[[1]][2]),
      as.numeric(strsplit(d[grepl("Area Under ROC", d)], ":")[[1]][2]),
      as.numeric(strsplit(d[grepl("Probability classification threshold mean: ", d)], ":")[[1]][2]))
      v[i,][2:7] <- n
  }
write.csv(v, "birds_validation_summary.csv")

# spp.dir <- basename(list.dirs())[-1]
# vf <- unlist(lapply(strsplit(list.files(getwd(), "txt$", recursive = TRUE, full.names = TRUE),"/"),
#             \(x) x[4]))		
# spp.dir[which(!spp.dir %in% vf)]
# # c("Anas penelope", "Calidris alpina", "Calidris pugnax")