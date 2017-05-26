library(survival)
library(GGally)
library(ggplot2)

defaultEncoding <- "UTF8"

change.files <- function(filename){
    options(warn=-1)
    library_data <- filename
    OV <- read.table(file=library_data,header=TRUE,sep="\t",row.names=NULL)
    exp <- survdiff(Surv(OV$time, OV$death) ~ factor(OV$cluster))
    nclst <- length(exp$n)
    pval <- (1 - pchisq(exp$chisq, nclst - 1))

    print(filename)
    print(exp$chisq)
    print(pval)

    fitMeta <- survfit(Surv(OV$time, OV$death) ~ (OV$cluster))
    print(fitMeta, print.rmean=T)
    third_slash<-which(strsplit(filename, "")[[1]]=="/")[3]

    pl <- ggsurv(fitMeta) + guides(linetype = FALSE) + scale_colour_manual(
        # SWITCH COLORS HERE. CHANGE THEM TO FIT WHAT YOU NEED. TODO.
        values=c('0'='#00BFC4', '1'='#F8766D'),
        # CHANGE CLUSTER SIZE HERE. TODO.
        name='Cluster Size', labels = c('21','69'), guide=guide_legend(reverse=TRUE)) + labs(x='Time (months)',
        # name='Cluster Size') + labs(x='Time (months)',
        y='Probability of survival') + theme(aspect.ratio=0.66,
        text=element_text(size=12), legend.justification=c(1,1),
        legend.position=c(1,1), legend.background = element_rect(fill=alpha(
            'white', 0.8))) + xlim(0,40)

    ggsave(filename=paste('./results/survival_plots_seq', substr(filename,
        third_slash, nchar(filename) - 4), '.pdf', sep=''), plot=pl, width=5,
    height=3.3)
}

args <- commandArgs(trailingOnly = TRUE)
change.files(args[1])