setwd("Z:/croatia/birds") 
zdir <- file.path("Z:/croatia") 

spp.scores <- read.csv("bird_scores.csv")
  spp.scores$priority[which(is.na(spp.scores$priority))] <- 0
  
mdls <- basename(list.dirs())[-1]
  for(s in 0:4) {
    spp.sub <- spp.scores[which(spp.scores$priority == s),]$species 
    sub.mdls <- spp.sub[which(spp.sub %in% mdls)]
	cat("Compressing", length(sub.mdls),  "rank", s, "bird models", "\n") 
    zip(zipfile=file.path(zdir, paste0("birds_score", s)), 
	    files=file.path(getwd(), sub.mdls))
    mdls <- mdls[-which(mdls %in% spp.sub)]
  } 
if(length(mdls) > 0) {
 cat("Compressing", length(mdls),  "rank", "other bird models", "\n") 
  zip(zipfile=file.path(zdir, paste0("birds_score_other")), 
	  files=file.path(getwd(), mdls))
}
