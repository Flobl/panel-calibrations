### Prepare file info cheet, similar to Antoine
### AS Oct 9 2019 ############################

info <- as.data.frame(matrix(nrow=70,ncol=6))
colnames(info) <- c("prefix", "spectrum_number", 
                    "cleaning","type","group","spec_num") ### cleaning needed? keep just in case

info$prefix <- "gr062519" ### date of measurements
info$spectrum_number <- c(paste0(rep("000",9),as.character(seq(0,9))),
  paste0(rep("00",nrow(info)-10),as.character(seq(10,nrow(info)-1))))
info$cleaning <- "after"
info$type <- c(rep("ref",5), rep("stray",5), rep("target",5),
  rep(c(rep("ref",5), rep("target",5)),5),rep("ref",5))
info$group <- c(rep("1",15), rep("2",10),rep("3",10), rep("4",10), rep("5",10), 
                rep("6",10), rep("7",5))
info$spec_num <- c(seq(1,15), rep(seq(1,10),5), seq(1,5))
info$rep <- rep(seq(1,5),14)

write.csv(info,"./processing_chain/docu/2019-06-25_docu/2019-06-25-Calibration-99AA02-0318-9712.csv", 
          row.names = F)


