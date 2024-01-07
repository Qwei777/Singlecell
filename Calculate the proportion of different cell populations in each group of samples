unique(`eye.markers`$Sample)
eye
samples=list.files('eye')
samples 
table(eye$orig.ident)#View the number of cells in each group
prop.table(table(Idents(eye)))
table(Idents(eye), eye$orig.ident)
Cellratio <- prop.table(table(Idents(eye), eye$orig.ident), margin = 2)#Calculate the proportion of different cell populations in each group of samples
Cellratio
Cellratio <- as.data.frame(Cellratio)
colourCount = length(unique(Cellratio$Var1))
library(ggplot2)
##pie plot
##ggplot(Cellratio, aes(x=Var2, y=Freq, fill=Var1)) +
  #geom_bar(stat="identity", width=1) +
  #coord_polar("y", start=0)
ggplot(Cellratio,aes(x=Var2, y=Freq, fill=Var1)) +
  geom_bar(stat="identity", width=1, color="white") +
  coord_polar("y", start=0) +theme_void()+scale_fill_manual(values = c("#DC83783","#BD886E","#638CFE","#0ABF57","#9EA3237"))#+ scale_fill_brewer(palette="Set2")# remove background, grid, numeric labels
