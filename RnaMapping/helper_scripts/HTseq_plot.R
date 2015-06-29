library(ggplot2)
library(reshape)

# SCRIPT to collect the counting statistics from HTSeq

reader <- function(df){
    data = read.table(df, sep ="\t", 
                      fill =TRUE, stringsAsFactors = FALSE)
    indx <- grepl('__', data[,1])
    out <- rbind('new' =c('__feature_found' , sum(data[,2][!indx])),
                 data[indx,])
    row.names(out) = gsub('__', '', out[,1])
    out[,1] <- NULL
    return(out)
}


HTSeq.logs <- function(path=NULL){
    if(is.null(path)){
        path=getwd()
    }
    logs = list.files(path=path, pattern= "count")
    logsF = file.path(path, logs)
    data = do.call("cbind", lapply(logsF, function(df) reader(df)))
    data1 = apply(data, 2, as.numeric)
    #data1 = apply(data1, 2, function(x) x/sum(x))
    row.names(data1) = row.names(data)
    colnames(data1) = gsub(".count", "", logs)
    return(data1)
}

plot.gg2 <- function(x, main){
    my.theme = theme_bw()+
        theme(axis.text.x = element_text(angle= 90, vjust=0.5, hjust=1),
              legend.position ="none",
              strip.text = element_text(hjust = 0.5, vjust = 0.5),
              plot.title = element_text(size=20, hjust=0.5))
#plot.margin= unit(c(0.5, 0.5, -0.5, 0.5), "lines")
    gg = ggplot(x)+aes(x=X2, y =value)+
        geom_bar(stat="identity",aes(fill =X2))+
            facet_wrap(~X1, scales="free_y", ncol=3)+
                xlab("")+ylab("")+ggtitle(main)+
                    my.theme
    return(gg)
}




data <- HTSeq.logs()
write.table(data, file = "HTSeq_Count_Stats.txt", sep="\t", col.names=T, row.names=T, quote=F)
datM <- melt(data)
datM$X1 <- factor(datM$X1, levels = c("feature_found", "no_feature", "ambiguous",
                               "alignment_not_unique", "too_low_aQual", "not_aligned"))
pdf("HTSeq_Count_Stats.pdf", height=6, width=11)
plot.gg2(datM, "HTSeq counting statistics")
dev.off()

