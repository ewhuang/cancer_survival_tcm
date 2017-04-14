library(survival)
library(GGally)
library(grid)
library(ggplot2)

defaultEncoding <- "UTF8"

change.files <- function(filename){
    library_data <- filename
    OV <- read.table(file=library_data,header=TRUE,sep="\t",row.names=NULL)

    # setpdf(width = 9, height=6)
    fitMeta <- survfit(Surv(OV$time, OV$death) ~ (OV$cluster))
    # postscript('./results/survival_plots_seq/paper_plot_2.pdf')
    # postscript('./results/survival_plots_seq/paper_plot_1.pdf')
    # postscript('./results/survival_plots_seq/paper_plot_no_prosnet_1.pdf')
    # postscript('./results/survival_plots_seq/paper_plot_no_prosnet_2.pdf')

    print(fitMeta, print.rmean=T)

    # GGplot?

    pl <- ggsurv(fitMeta)
    pl <- pl + ggplot2::guides(linetype = FALSE) + ggplot2::scale_colour_discrete(
      name   = 'Cluster Size',
      labels = fitMeta$n
    ) + labs(x='Time (months)', y='Probability of survival') + theme(
    aspect.ratio=0.66, text=element_text(size=12), legend.justification=c(1,1),
    legend.position=c(1,1), legend.background = element_rect(fill=alpha('white',
        0.8))) + xlim(0,40)


    ggsave(filename='./results/survival_plots_seq/paper_plot_2.pdf',
    # ggsave(filename='./results/survival_plots_seq/paper_plot_1.pdf',
    # ggsave(filename='./results/survival_plots_seq/paper_plot_no_prosnet_1.pdf',
    # ggsave(filename='./results/survival_plots_seq/paper_plot_no_prosnet_2.pdf',
    # ggsave(filename='./results/survival_plots_seq/paper_plot_vkps_1.pdf',
    # ggsave(filename='./results/survival_plots_seq/paper_plot_vkps_2.pdf',
    # ggsave(filename='./results/survival_plots_seq/paper_plot_mean_1.pdf',
    # ggsave(filename='./results/survival_plots_seq/paper_plot_mean_2.pdf',
        plot=pl, width=5, height=3.3)
    # ggsurv(fitMeta, col=c(4,3), xlim=c(0, 50), lty=c(2,1),
    #     xlab='Time (months)', ylab="Probablity of survival")

    # legend_labels<-fitMeta$n

    # legend("topright", 1.0, legend_labels, lty=c(2,1), lwd=c(2.5,2.5),
    #     col=c(4,3), title="Cluster Size")

    # Title of plot.
    # title('Survival Functions for Squamous-Cell Carcinoma Patients')
    # title('Survival Functions for Non-SQ NSCLC Patients')

    # Save plot.
    # dev.off()
}

change.files('./data/patient_dataframes_seq/prosnet_2_100.txt')
# change.files('./data/patient_dataframes_seq/prosnet_1_100.txt')
# change.files('./data/patient_dataframes_seq/without_prosnet_1.txt')
# change.files('./data/patient_dataframes_seq/without_prosnet_2.txt')
# change.files('./data/patient_dataframes_seq/vkps_1.txt')
# change.files('./data/patient_dataframes_seq/vkps_2.txt')
# change.files('./data/patient_dataframes_seq/mean_1.txt')
# change.files('./data/patient_dataframes_seq/mean_2.txt')