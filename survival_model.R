library(survival)

#method_prefix <- "0.5_4_4_1_4_1"
defaultEncoding <- "UTF8"
#setwd("C:/Users/ewhuang3/Documents/tcm_project")

change.files <- function(filename){
    print(filename)
    library_data<-filename
    OV <- read.table(file=library_data,header=TRUE,sep="\t",row.names=NULL)

    exp<-survdiff(Surv(OV$time, OV$death) ~ factor(OV$cluster))

    nclst <- length(exp$n)
    pv<-1 - pchisq(exp$chisq, nclst - 1)
    pv_new<-prettyNum(pv, digits=3, width=4, format="fg")

    fitMeta <- survfit(Surv(OV$time, OV$death) ~ (OV$cluster))

    png(paste('./results/survival_plots/', paste(substr(filename, 27,
        nchar(filename) - 4), "png", sep="."), sep=''))

    plot(fitMeta,col= rainbow(nclst),xlim=c(0,1000),xlab='Time(days)',
        ylab="Probablity of death event")
    legend(800,0.8, exp$n, lty=c(1,1), lwd=c(2.5,2.5), col=rainbow(nclst))
    title(main=paste('p=', pv_new,sep=""))
    dev.off()
}

files <- list.files(path="./data/patient_dataframes", full.names = TRUE)
lapply(files, change.files)
