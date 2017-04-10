library(survival)

defaultEncoding <- "UTF8"

change.files <- function(filename){
    library_data <- filename
    OV <- read.table(file=library_data,header=TRUE,sep="\t",row.names=NULL)

    setEPS(width = 9, height=6)
    fitMeta <- survfit(Surv(OV$time, OV$death) ~ (OV$cluster))
    postscript('./results/survival_plots_seq/paper_plot_2.eps')
    # postscript('./results/survival_plots_seq/paper_plot_1.eps')
    # postscript('./results/survival_plots_seq/paper_plot_no_prosnet_1.eps')
    # postscript('./results/survival_plots_seq/paper_plot_no_prosnet_2.eps')

    print(fitMeta, print.rmean=T)

    # GGplot?
    
    plot(fitMeta, col=c(4,3), xlim=c(0, 50), lty=c(2,1),
        xlab='Time (months)', ylab="Probablity of survival")

    legend_labels<-fitMeta$n

    legend("topright", 1.0, legend_labels, lty=c(2,1), lwd=c(2.5,2.5),
        col=c(4,3), title="Cluster Size")

    # Title of plot.
    # title('Survival Functions for Squamous-Cell Carcinoma Patients')
    # title('Survival Functions for Non-SQ NSCLC Patients')

    # Save plot.
    dev.off()
}

change.files('./data/patient_dataframes_seq/prosnet_2_100.txt')
# change.files('./data/patient_dataframes_seq/prosnet_1_100.txt')
# change.files('./data/patient_dataframes_seq/without_prosnet_1.txt')
# change.files('./data/patient_dataframes_seq/without_prosnet_2.txt')