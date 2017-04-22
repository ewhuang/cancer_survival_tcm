library(survival)
# library(MASS)
# library(risksetROC)
library(survAUC)

defaultEncoding <- "UTF8"

set.seed(9305)

change.files <- function(filename){
    options(warn=-1)

    library_data <- filename
    OV <- read.table(file=library_data,header=TRUE,sep="\t",row.names=NULL,
        check.names = FALSE)
    # Shuffle the data.
    OV<-OV[sample(nrow(OV)),]
    # Create three folds.
    folds <- cut(seq(1,nrow(OV)),breaks=3,labels=FALSE)
    sum<-0

    res.cox<-coxph(Surv(time, death) ~ . - patient_id - death - time, data=OV)
    test.ph <- cox.zph(res.cox)
    print(test.ph)
    # print(test.ph$table)
    # print(test.ph[is.element(test.ph$p, NaN)])
    exit()

    for(i in 1:3) {
        # Segment data by fold using which() function.
        testIndices<-which(folds==i,arr.ind=TRUE)
        testData<-OV[testIndices, ]
        trainData<-OV[-testIndices, ]
        # Train Cox model.
        train.fit<-coxph(Surv(time, death) ~ . - patient_id - death - time,
            data=trainData, singular.ok=TRUE, na.action=na.omit)
        lp <- predict(train.fit)
        lpnew <- predict(train.fit, newdata=testData)
        Surv.rsp <- Surv(trainData$time, trainData$death)
        Surv.rsp.new <- Surv(testData$time, testData$death)
        # Evaluate every 5 months.
        times <- seq(0, 50, 5)

        # Compute AUC.
        AUC_hc <- AUC.hc(Surv.rsp, Surv.rsp.new, lpnew, times)
        print(times)
        sum<-sum + AUC_hc$iauc
        print(AUC_hc$iauc)
    }
    print(sum / 3.0)
}

args <- commandArgs(trailingOnly = TRUE)
change.files(args[1])
