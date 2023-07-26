###data prepare###
load("Data/explore.Rdata")
keep = rowAlls(exp>0)
table(keep)
exp_filter = exp[keep,]
exp_filter["PIR",]
exp_filter[1:4,1:4]
dim(exp_filter)
Group <- c('Wt','Wt','Wt','Tg','Tg','Tg')
Group=factor(Group,levels = c('Wt','Tg'))
Group
exp_cpm <- cpm(exp_filter)
exp_cpm[1:6,1:6]
mythe <- theme_bw() + theme(panel.grid.major=element_blank(),
                            panel.grid.minor=element_blank())
exp <- log10(as.matrix(exp_cpm)+1)
exp[1:6,1:6]

###PCA###
data1=as.data.frame(colnames(exp_filter))
data1$group=Group
data1$condition=Group
data1$cluster=Group
colnames(data1)[1]="sample_id"
data2=log2(exp_filter+1)
data2=as.data.frame(data2)
pca = prcomp(t(data2))
data1 <- data.frame("PC1" = pca$x[,1], 
                    "PC2" = pca$x[,2]) %>% 
  rownames_to_column("sample_id") %>% 
  left_join(., data1, by = "sample_id") 
rownames(data1)=data1$sample_id
ggplot(data=data1,aes(x = PC1, y = PC2)) + 
  geom_hline(yintercept = 0, linetype = 2, size = 0.75*0.47) +
  geom_vline(xintercept = 0, linetype = 2, size = 0.75*0.47) +
  theme_set(theme_set(theme_bw(base_size=20))) + theme_bw() +
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank()) +
  geom_point(aes(colour = group,fill = group)) +
  geom_label_repel(aes(label = rownames(data1),
                       fill = factor(group)), color = 'white',
                   size = 3.5) 

###DEG####

comp <- unlist(strsplit("Tg_vs_Wt",split = "_vs_"))
table(Group)
design <- model.matrix(~0+Group)
rownames(design) <- colnames(exp_filter)
colnames(design) <- levels(factor(comp))
DEG <- DGEList(counts=exp_filter, 
               group=Group)
DEG <- calcNormFactors(DEG)

DEG <- estimateGLMCommonDisp(DEG,design)
DEG <- estimateGLMTrendedDisp(DEG, design)
DEG <- estimateGLMTagwiseDisp(DEG, design)

fit <- glmFit(DEG, design)
lrt <- glmLRT(fit, contrast=c(1,-1)) 
DEG_edgeR <- as.data.frame(topTags(lrt, n=nrow(DEG)))
head(DEG_edgeR)

fc_cutoff <- 2
pvalue <- 0.05

DEG_edgeR$regulated <- "stable"

loc_up <- intersect(which( DEG_edgeR$logFC > log2(fc_cutoff) ),
                    which( DEG_edgeR$PValue < pvalue) )

loc_down <- intersect(which(DEG_edgeR$logFC < (-log2(fc_cutoff))),
                      which(DEG_edgeR$PValue<pvalue))

DEG_edgeR$regulated[loc_up] <- "up"
DEG_edgeR$regulated[loc_down] <- "down"
