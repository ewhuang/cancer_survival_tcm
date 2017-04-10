library(survival)

defaultEncoding <- "UTF8"

change.files <- function(filename){
    library_data <- filename
    OV <- read.table(file=library_data,header=TRUE,sep="\t",row.names=NULL)
    exp <- survdiff(Surv(OV$time, OV$death) ~ factor(OV$cluster))
    nclst <- length(exp$n)
    pv <- (1 - pchisq(exp$chisq, nclst - 1))
    # Different p-value thresholds for the full clustering, since for the
    # occurrence clustering, we have a lot of pairs of elements to iterate
    # through.
    if (grepl('full', filename) || grepl('seq', filename)) {
        p_thresh<-1.0
    } else {
        p_thresh<-0.01
    }
    # TODO: Cannot print filename for plotting threshold sensitivity.
    print(filename)
    print(exp$chisq)
    print(pv)
    if (pv < p_thresh) {
        pv_new<-prettyNum(pv, digits=3, width=4, format="fg")
        fitMeta <- survfit(Surv(OV$time, OV$death) ~ (OV$cluster))
        # PNG filename. TODO.
        second_underscore<-which(strsplit(filename, "")[[1]]=="_")[2]
        if (grepl('synergy', filename)) {
            untreat_med<-summary(fitMeta)$table["OV$cluster=drug_only", 'median']
            treat_med<-summary(fitMeta)$table["OV$cluster=both", "median"]
        } else if (grepl('treatment', filename)) {
            untreat_med<-summary(fitMeta)$table["OV$cluster=not_treated", "median"]
            treat_med<-summary(fitMeta)$table["OV$cluster=treated", "median"]
        }
        # TODO: Only write out plots that have clusters where the treatment
        # cluster is better than the untreated cluster.
        if (grepl('full', filename) || grepl('seq', filename) || (!is.na(untreat_med) && !is.na(treat_med) && treat_med > untreat_med)) {
            # if (grepl('synergy', filename)) {
            slash_indices<-which(strsplit(filename, "")[[1]]=="/")
            last_underscore_idx<-tail(slash_indices, n=1)
            # }
            png(paste('./results/survival_plots', paste(substr(filename,
                second_underscore, nchar(filename) - 4), "png", sep="."),
            sep=''))

            plot(fitMeta, col=rainbow(nclst), xlim=c(0, 48),
                xlab='Time (months)', ylab="Probablity of survival")

            cluster_names<-c()
            for (name in names(exp$n)) {
                cluster_names<-c(cluster_names, substr(name, which(
                    strsplit(name,"")[[1]]=="=")[1] + 1, nchar(name)))
            }

            legend_labels<-paste0(exp$n, c(', ', ', ') , cluster_names)

            legend(25, 1.0, legend_labels, lty=c(1,1), lwd=c(2.5,2.5),
                col=rainbow(nclst), title="Number of cluster patients")

            # Title of plot.
            title(main=paste(substr(filename, last_underscore_idx+1, nchar(
                filename) - 4), ', p=', pv_new, sep=""))

            # Save plot.
            dev.off()
        }
    }
}

args <- commandArgs(trailingOnly = TRUE)
files <- list.files(args[1], full.names = TRUE)
invisible(lapply(files, change.files))
