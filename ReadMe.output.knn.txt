# Readme about knn error model run

- 1. knn.error.models.allData.Rds
    6/7/2023
        all data set.
    cd <- clean.counts(mcd)
# check the final dimensions of the read count matrix
dim(cd)

rm(mcd)
gc()

system.time(
knn <- knn.error.models(cd, k = ncol(cd)/50,
		n.cores = 9, #running parallel::detectCores() to see how many cores on your machine
		 min.count.threshold = 1,#this has to be 1 for umi counts (check the help page on scde website)
		 min.nonfailed = 50,
		save.model.plots=T,max.model.plots = 4, verbose=4)
)

user    system   elapsed
78666.753   807.402 32231.602
