create_gc_plot <- function(gc, id){
  gc <- as.double(gc)
  at <- 100-gc
  slices <- c(gc, at)
  lbls <- c("GC", "AT")
  lbls <- paste(lbls, slices)
  lbls <- paste(lbls,"%",sep="")
  file_name <- paste0("data/figures/",id,".png")
  png(file_name)
  pie(slices, labels = lbls, init.angle = 90, col=c("red","blue"),main=id)
  dev.off()

}

sequences <- read.table(snakemake@input[[1]], header=F)
ids <- read.table(snakemake@input[[2]], header=F)

for (i in 1:nrow(sequences)){
  gc <- sequences[i,][1]
  id <- ids[i,]
  create_gc_plot(gc, id)
}
