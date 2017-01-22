library(survival)

#method_prefix <- "0.5_4_4_1_4_1"
defaultEncoding <- "UTF8"
#setwd("C:/Users/ewhuang3/Documents/tcm_project")

change.files <- function(filename){
    library_data <- filename
    OV <- read.table(file=library_data,header=TRUE,sep="\t",row.names=NULL)

    exp <- survdiff(Surv(OV$time, OV$death) ~ factor(OV$cluster))

    nclst <- length(exp$n)
    pv <- (1 - pchisq(exp$chisq, nclst - 1))

    p_thresh<-0.01
    if (pv < p_thresh) {
        pv_new<-prettyNum(pv, digits=3, width=4, format="fg")
        fitMeta <- survfit(Surv(OV$time, OV$death) ~ (OV$cluster))

        # PNG filename. TODO.
        underscore_indices<-which(strsplit(filename, "")[[1]]=="_")
        second_idx<-underscore_indices[2]

        # if (grepl('synergy', filename)) {
        slash_indices<-which(strsplit(filename, "")[[1]]=="/")
        last_underscore_idx<-tail(slash_indices, n=1)
        # }
        png(paste('./results/survival_plots', paste(substr(filename, second_idx,
            nchar(filename) - 4), "png", sep="."), sep=''))

        plot(fitMeta, col=rainbow(nclst), xlim=c(0, 48), xlab='Time (months)',
            ylab="Probablity of survival")

        cluster_names<-c()
        for (name in names(exp$n)) {
            cluster_names<-c(cluster_names, substr(name, which(strsplit(name,"")
                [[1]]=="=")[1] + 1, nchar(name)))
        }

        legend_labels<-paste0(exp$n, c(', ', ', ') , cluster_names)

        legend(25, 1.0, legend_labels, lty=c(1,1), lwd=c(2.5,2.5), col=rainbow(
            nclst), title="Number of cluster patients")

        # Title of plot.
        title(main=paste(substr(filename, last_underscore_idx+1, nchar(
            filename) - 4), ', p=', pv_new, sep=""))

        # Save plot.
        dev.off()
    }
}

args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 1) {
    files <- list.files(path=paste("./data/patient_dataframes_", args[1],
        "_features", sep=''), full.names = TRUE)
} else {
    files <- list.files(path=paste("./data/patient_dataframes_", args[1],
        "_features_", args[2], sep=''), full.names = TRUE)    
}
invisible(lapply(files, change.files))
