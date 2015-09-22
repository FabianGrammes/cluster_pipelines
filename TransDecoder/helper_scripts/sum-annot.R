options(echo=TRUE) # if you want see commands in output file
args <- commandArgs(trailingOnly = TRUE)
print(args)

infile <- args[1]
outfile <- args[2]
#-------------------------------------------------------------------------------

library(ggplot2)
library(gridExtra)
library(grid)

dd <- read.table(file = infile, header = TRUE, sep ='\t',
                 stringsAsFactors = FALSE)

#-------------------------------------------------------------------------------
names.tab <- sapply( dd$gene_name, function(x) unlist(strsplit(x, '\\|'))[2], USE.NAMES = FALSE)
names.tab <- as.data.frame(table(names.tab))


g.nt <- ggplot(names.tab)+aes(x = names.tab, y = Freq)+
    geom_bar(stat="identity", fill = 'steelblue2', colour = 'black' )+
    labs(list(title = 'Hit classification', x='', y = '') )+
    theme_bw()+
    theme(legend.position = 'none',
          axis.text.x = element_text( angle = 90, hjust = 1, vjust = 0.5))
#-------------------------------------------------------------------------------
g.ht = ggplot(dd)+aes(x=as.numeric(coverage))+
    geom_histogram(binwidth = 0.05, fill = 'steelblue2', colour = 'black')+
    labs(list(title = 'Histogram: Hit coverage', x='Hit coverage', y = 'count') )+
    theme_bw()
#-------------------------------------------------------------------------------



org.tab <- as.data.frame(table(dd$origin))
org.tab <- head(org.tab[order(org.tab$Freq, decreasing = TRUE),],10)
org.tab$Var1 <- factor(org.tab$Var1, levels =org.tab$Var1, ordered = TRUE)


g.org <- ggplot(head(org.tab,10))+aes(x = Var1, y = Freq)+
    geom_bar(stat="identity", fill = 'steelblue2', colour = 'black' )+
    labs(list(title = 'Top 10 hit sources', y ='count', x = '') )+
    theme_bw()+
    theme(legend.position = 'none',
          axis.text.x = element_text( angle = 90, hjust = 1, vjust = 0.5))

#-------------------------------------------------------------------------------
pdf(outfile, width = 10, height = 5)
lay <- grid.layout(1, 8 )
pushViewport(viewport(layout = lay))
print(g.ht, vp = viewport(layout.pos.col = 1:3, layout.pos.row = 1))
print(g.org, vp = viewport(layout.pos.col = 4:6, layout.pos.row = 1))
print(g.nt, vp = viewport(layout.pos.col = 7:8, layout.pos.row = 1))
dev.off()

cat('-------------------------------------------------------------------------------\n')
names.tab
cat('-------------------------------------------------------------------------------\n')
org.tab
cat('-------------------------------------------------------------------------------\n')
table(grepl('Uncharacterized protein|uncharacterized protein', dd$gene_name))
cat('-------------------------------------------------------------------------------\n')
nrow(dd)
