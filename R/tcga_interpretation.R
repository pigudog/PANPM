library(survminer)
library(survival)
library(readxl)
library(dplyr)
TCGA_CESC_OS <- read_excel("dl_output/result_tcga_OS_final_12char.xlsx")
TCGA_CESC_OS$case_id = substr(TCGA_CESC_OS$case_id,1,12)
head(TCGA_CESC_OS)

load("./rdata/cli_combined_selected.rda")
# Identify IDs in FUDAN_CESC_OS that are not in combined_Data
missing_ids <- TCGA_CESC_OS %>%
  filter(!case_id %in% os_pfs_combined$ID)

# Print the missing IDs
print(missing_ids)


# Remove rows with missing IDs from FUDAN_CESC_OS
TCGA_CESC_OS <- TCGA_CESC_OS %>%
  filter(case_id %in% os_pfs_combined$ID)
TCGA_CESC_OS$status = abs(TCGA_CESC_OS$censorship-1)
TCGA_CESC_OS$Group=ifelse(TCGA_CESC_OS$risk > median(TCGA_CESC_OS$risk),'High','Low')
head(TCGA_CESC_OS)


maf_data@clinical.data$Tumor_Sample_Barcode
load("./rdata/coad_snp_tcga.Rdata")
maf.coad<-data
load("./rdata/read_snp_tcga.Rdata")
maf.read<-data
library(maftools)
maf_data = rbind(maf.coad,maf.read)
maf_data <- read.maf(maf_data)
# 筛选 MAF 数据，只保留相应的样本
# 映射一下
# 创建一个包含 TCGA_CESC_OS 中所有 case_id 的数据框
# 所以首先顺序是TCGA的顺序
case_id_mapping <- data.frame(
  case_id = TCGA_CESC_OS$case_id,
  stringsAsFactors = FALSE
)

# 为每个 case_id 找到对应的完整条形码
# 将maf数据一一对应上去
case_id_mapping$full_barcode <- maf_data@clinical.data$Tumor_Sample_Barcode[match(case_id_mapping$case_id,
                                                                                  substr(maf_data@clinical.data$Tumor_Sample_Barcode, 1, 12))]

# 检查映射是否成功
print(case_id_mapping)

# 筛选 MAF 数据，保留相应的样本
# 只取full_barcode不为na的，也就是cli和maf的交集
maf_filtered <- subsetMaf(maf = maf_data,
                          tsb = case_id_mapping$full_barcode[!is.na(case_id_mapping$full_barcode)])

# 获取有效的分组信息
# 该函数用于查找 case_id_mapping$case_id 中每个元素在 TCGA_CESC_OS$case_id 中的位置，返回一个索引向量，表示匹配的元素位置。
# 所以拿到的是对应于case_id_mapping的Group信息
# 注意valid_group的顺序跟tcga——cli是一致的
valid_groups <- TCGA_CESC_OS$Group[match(case_id_mapping$case_id, TCGA_CESC_OS$case_id)]
# 去掉maf里面没有的使id跟
valid_groups <- valid_groups[!is.na(case_id_mapping$full_barcode)]  # 确保分组信息不为 NA
# names(valid_groups) = case_id_mapping$full_barcode[!is.na(case_id_mapping$full_barcode)]
identical(substr(maf_filtered@clinical.data$Tumor_Sample_Barcode,1,12),case_id_mapping[!is.na(case_id_mapping$full_barcode),]$case_id)
maf_filtered@clinical.data$Group = valid_groups
# # 假设 maf_filtered@clinical.data 是一个数据框
# maf_filtered@clinical.data <- maf_filtered@clinical.data %>%
#   arrange(Group)
# maf_filtered@data$Tumor_Sample_Barcode = factor(maf_filtered@data$Tumor_Sample_Barcode,levels =maf_filtered@clinical.data$Tumor_Sample_Barcode )

# 分组绘制 oncoplot
oncoplot(maf = maf_filtered, top = 20, clinicalFeatures  = "Group",sortByAnnotation = TRUE)
# oncoplot(maf = maf_filtered, top = 20)
# 改变颜色
ov_palette =c('#A499CC',
              '#E0A7C8',
              '#E069A6',
              "#f1707d",
              "#AFC2D9",
              "#6894B9",
              "#79B99D",
              "#F5D2A8",
              "#D2EBC8")
col = RColorBrewer::brewer.pal(n = 8, name = 'Paired')[1:9]
col = ov_palette[1:9]
names(col) = c('Frame_Shift_Del','Missense_Mutation', 'Nonsense_Mutation',

               'Multi_Hit', 'Frame_Shift_Ins',

                'Splice_Site', 'In_Frame_Del')
oncoplot(maf = maf_filtered, top = 20, clinicalFeatures  = "Group",
         sortByAnnotation = TRUE,
         colors = col)
oncoplot(maf = maf_filtered, top = 20, clinicalFeatures = "Group",
         sortByAnnotation = TRUE,
         colors = col,bgCol = "white"
         # bgCol = "#DFDFDF",borderCol = "#DFDFDF"
         )  # 将缺失值颜色设置为白色


#设置不同样本分组的颜色

# fabcolors = RColorBrewer::brewer.pal(n = 8,name = 'Spectral')

# names(fabcolors) = c("M0", "M1", "M2", "M3", "M4", "M5", "M6", "M7")
#
# fabcolors = list(FAB_classification = fabcolors)

# oncoplot(maf = maf_data,
#
#          colors = col,
#
#          clinicalFeatures = 'FAB_classification',
#
#          sortByAnnotation = TRUE,
#
#          annotationColor = fabcolors)


################################################################################
cnv = data.table::fread("./GDCdata/Gistic2_CopyNumber_Gistic2_all_thresholded.by_genes")
cnv <- cnv %>%
  column_to_rownames("sample")   # Set the new column as row names
cnv[cnv>0] <- "gain"
cnv[cnv<0] <- "loss"
cnv[cnv == 0] <- "neutral"
colnames(cnv) = stringr::str_sub(colnames(cnv),1,12)
library(survminer)
library(survival)
library(readxl)
library(dplyr)
TCGA_CESC_OS <- read_excel("dl_output/result_tcga_OS_final_12char.xlsx")
TCGA_CESC_OS$case_id = substr(TCGA_CESC_OS$case_id,1,12)
head(TCGA_CESC_OS)

load("./rdata/cli_combined_selected.rda")
# Identify IDs in FUDAN_CESC_OS that are not in combined_Data
missing_ids <- TCGA_CESC_OS %>%
  filter(!case_id %in% os_pfs_combined$ID)

# Print the missing IDs
print(missing_ids)


# Remove rows with missing IDs from FUDAN_CESC_OS
TCGA_CESC_OS <- TCGA_CESC_OS %>%
  filter(case_id %in% os_pfs_combined$ID)
TCGA_CESC_OS$status = abs(TCGA_CESC_OS$censorship-1)
TCGA_CESC_OS$Group=ifelse(TCGA_CESC_OS$risk > median(TCGA_CESC_OS$risk),'High','Low')
head(TCGA_CESC_OS)

# 首先根据 TCGA_CESC_OS 中的 case_id 过滤 cnv 的列名
filtered_cnv <- cnv[, colnames(cnv) %in% TCGA_CESC_OS$case_id]

# 将 case_id 添加到 filtered_cnv 中，以便后续合并
filtered_cnv_df <- as.data.frame(t(filtered_cnv))
filtered_cnv_df$case_id <- rownames(filtered_cnv_df)

# 将 TCGA_CESC_OS 中的 Group 信息添加到 filtered_cnv_df
filtered_cnv_df <- merge(filtered_cnv_df, TCGA_CESC_OS[, c("case_id", "Group")], by.x = "case_id", by.y = "case_id")

# 分成 High 和 Low
filtered_cnv_high <- filtered_cnv_df[filtered_cnv_df$Group == "High", ]
filtered_cnv_low <- filtered_cnv_df[filtered_cnv_df$Group == "Low", ]

# # 处理 High 组
# filtered_cnv_high$Group = NULL
# filtered_cnv_high$case_id = NULL
# df3_high <- apply(filtered_cnv_high, 2, table)
# # df3_high <- as.data.frame(t(do.call(rbind,df3_high)))
# df3_high <- as.data.frame(df3_high)
# df3_high["mutation", ] <- df3_high["gain", ] + df3_high["loss", ]
# df3_high_transposed <- as.data.frame(t(df3_high))
# df3_high_transposed$gene <- rownames(df3_high_transposed)
#
# # 提取 gain 和 loss 前 5 个基因
# top_5_gain_high <- head(df3_high_transposed[order(df3_high_transposed$gain, decreasing = TRUE), ], 5)
# top_5_loss_high <- head(df3_high_transposed[order(df3_high_transposed$loss, decreasing = TRUE), ], 5)
#
# # 合并
# top_genes_high <- unique(rbind(top_5_gain_high, top_5_loss_high))
#
# # 绘制 High 组的图
# ggplot(top_genes_high) +
#   geom_linerange(aes(x = gene, ymin = 0, ymax = ifelse(gain > loss, gain, loss)),
#                  linewidth = 2, color = "grey") +
#   geom_point(aes(x = gene, y = gain), color = "#f1707d", size = 4) +
#   geom_point(aes(x = gene, y = loss), color = "#6894B9", size = 4) +
#   labs(x = NULL, y = "CNV frequency") +
#   theme_bw() +
#   theme(
#     axis.text.x = element_text(angle = 45, hjust = 1),
#     panel.grid.major = element_blank(),
#     panel.grid.minor = element_blank()
#   ) +
#   ggtitle("High Group")




# 处理 High 组
filtered_cnv_high$Group = NULL
filtered_cnv_high$case_id = NULL
df3_high <- apply(filtered_cnv_high, 2, table)
df3_high <- as.data.frame(t(do.call(rbind,df3_high)))
# df3_high <- as.data.frame(df3_high)
df3_high["mutation", ] <- df3_high["gain", ] + df3_high["loss", ]
df3_high_transposed <- as.data.frame(t(df3_high))
df3_high_transposed$gene <- rownames(df3_high_transposed)

# 提取 mutation 列并按降序排序
top_genes_high <- df3_high_transposed[order(df3_high_transposed$mutation, decreasing = TRUE), ]
top_10_genes_high <- rownames(head(top_genes_high, 10))
df4_high <- df3_high_transposed[df3_high_transposed$gene %in% top_10_genes_high, ]

# 绘制 High 组的图
ggplot(df4_high) +
  geom_linerange(aes(x = gene, ymin = 0, ymax = ifelse(gain > loss, gain, loss)),
                 linewidth = 2, color = "grey") +
  geom_point(aes(x = gene, y = gain), color = "#f1707d", size = 4) +
  geom_point(aes(x = gene, y = loss), color = "#6894B9", size = 4) +
  labs(x = NULL, y = "CNV frequency") +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.major = element_blank(),  # 删除主网格线
    panel.grid.minor = element_blank()   # 删除次网格线
  ) +
  ggtitle("High Group")

# 处理 Low 组
# 处理 Low 组
filtered_cnv_low$Group = NULL
filtered_cnv_low$case_id = NULL
filtered_cnv_low = as.data.frame(filtered_cnv_low)
df3_low <- apply(filtered_cnv_low, 2, table)
# df3_low <- as.data.frame(df3_low)
df3_low <- as.data.frame(t(do.call(rbind,df3_low)))
df3_low["mutation", ] <- df3_low["gain", ] + df3_low["loss", ]
df3_low_transposed <- as.data.frame(t(df3_low))
df3_low_transposed$gene <- rownames(df3_low_transposed)

# 提取 mutation 列并按降序排序

top_genes_low <- df3_low_transposed[order(df3_low_transposed$mutation, decreasing = TRUE), ]
top_10_genes_low <- rownames(head(top_genes_low, 10))
df4_low <- df3_low_transposed[df3_low_transposed$gene %in% top_10_genes_low, ]

# 绘制 Low 组的图
ggplot(df4_low) +
  geom_linerange(aes(x = gene, ymin = 0, ymax = ifelse(gain > loss, gain, loss)),
                 linewidth = 2, color = "grey") +
  geom_point(aes(x = gene, y = gain), color = "#f1707d", size = 4) +
  geom_point(aes(x = gene, y = loss), color = "#6894B9", size = 4) +
  labs(x = NULL, y = "CNV frequency") +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.major = element_blank(),  # 删除主网格线
    panel.grid.minor = element_blank()   # 删除次网格线
  ) +
  ggtitle("Low Group")
# 我觉得我现在的mutation这个定义不合适，我想gain的取top5，loss的取top5
# 合并 High 和 Low 组的数据
df4_high$Group <- "High"
df4_low$Group <- "Low"
combined_top_genes <- rbind(df4_high, df4_low)

# 绘制 High 和 Low 组的图
p1 = ggplot(df4_high) +
  geom_linerange(aes(x = gene, ymin = 0, ymax = ifelse(gain > loss, gain, loss)),
                 linewidth = 2, color = "grey") +
  geom_point(aes(x = gene, y = gain), color = "#f1707d", size = 4) +
  geom_point(aes(x = gene, y = loss), color = "#6894B9", size = 4) +
  labs(x = NULL, y = "CNV frequency") +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  ggtitle("CNV Frequency by Group") +
  facet_wrap(~ Group, scales = "free_y");p1  # 按组分面，y轴自适应
p2 = ggplot(df4_low) +
  geom_linerange(aes(x = gene, ymin = 0, ymax = ifelse(gain > loss, gain, loss)),
                 linewidth = 2, color = "grey") +
  geom_point(aes(x = gene, y = gain), color = "#f1707d", size = 4) +
  geom_point(aes(x = gene, y = loss), color = "#6894B9", size = 4) +
  labs(x = NULL, y = "CNV frequency") +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  ggtitle("CNV Frequency by Group") +
  facet_wrap(~ Group, scales = "free_y");p2
p1+p2


################################################################################
# 通路富集
library(survminer)
library(survival)
library(readxl)
library(dplyr)
TCGA_CESC_OS <- read_excel("dl_output/result_tcga_OS_final_12char.xlsx")
TCGA_CESC_OS$case_id = substr(TCGA_CESC_OS$case_id,1,12)
head(TCGA_CESC_OS)

load("./rdata/cli_combined_selected.rda")
# Identify IDs in FUDAN_CESC_OS that are not in combined_Data
missing_ids <- TCGA_CESC_OS %>%
  filter(!case_id %in% os_pfs_combined$ID)

# Print the missing IDs
print(missing_ids)


# Remove rows with missing IDs from FUDAN_CESC_OS
TCGA_CESC_OS <- TCGA_CESC_OS %>%
  filter(case_id %in% os_pfs_combined$ID)
TCGA_CESC_OS$status = abs(TCGA_CESC_OS$censorship-1)
TCGA_CESC_OS$Group=ifelse(TCGA_CESC_OS$risk > median(TCGA_CESC_OS$risk),'High','Low')
head(TCGA_CESC_OS)


mRNA = data.table::fread("./GDCdata/HiSeqV2")
library(dplyr)
library(tidyverse)
mRNA <- mRNA %>%
  column_to_rownames("sample")
colnames(mRNA) = substr(colnames(mRNA),1,12)



# 首先根据 TCGA_CESC_OS 中的 case_id 过滤 cnv 的列名
filtered_mRNA <- mRNA[, colnames(mRNA) %in% TCGA_CESC_OS$case_id]

# 将 case_id 添加到 filtered_cnv 中，以便后续合并
filtered_mRNA_df <- as.data.frame(t(filtered_mRNA))
filtered_mRNA_df$case_id <- rownames(filtered_mRNA_df)

# 将 TCGA_CESC_OS 中的 Group 信息添加到 filtered_cnv_df
filtered_mRNA_df <- merge(filtered_mRNA_df, TCGA_CESC_OS[, c("case_id", "Group")], by.x = "case_id", by.y = "case_id")
# head(filtered_mRNA_df[,c(20520:20532)])
rownames(filtered_mRNA_df) = filtered_mRNA_df$case_id
expression_matrix <- filtered_mRNA_df  # 去掉 Group 列

expression_matrix$case_id =NULL
expression_matrix$Group = NULL
expression_matrix = t(expression_matrix)
#数据的归一化处理，两种不同的代码
normalize<-t(apply(expression_matrix,1,function(x) x-(mean(x))))
#normalize<-sweep(exp,1,apply(exp,1,mean,na.rm=T))
# write.table(normalize,file="exp.txt",sep="\t",
#             row.names=T,col.names=TRUE,quote=FALSE)

# expression_matrix =
# rownames(expression_matrix) <- rownames(filtered_mRNA_df)
library(GSVA)
library(pheatmap)
library(clusterProfiler)
library(org.Hs.eg.db)
library(ggplot2)
# devtools::install_github("rcastelo/GSVA")
#选择和肿瘤相关的基因集
hsets <- read.gmt("~/guo_pathology_squamous_cervix/pathway/hallmark_cancersea.gmt")
# 假设 hsets_df 是你的数据框
hsets_list <- split(hsets$gene, hsets$term)
# 去掉前14个元素
hsets_list <- hsets_list[-(1:14)]

PANoptosis_data = data.frame(ID = colnames(expression_matrix))
library(GSEABase)
# gene_set = hsets[marker]
# gene_set = as.vector(hsets[marker]) # 确保是字符向量
# gene_set_collection = GeneSetCollection(GeneSet(gene_set))
# gene_set = as.character(gene_set)
# gene_set=list(gene_set)
# names(gene_set)=marker
# marker = names(gene_set)
#names(gene_set)='ALG5-Epi'
# library(genefilter)
ssgsea_params <- GSVA::gsvaParam( as.matrix(expression_matrix), hsets_list,kcdf = "Gaussian")
ssgsea_data = gsva(ssgsea_params, verbose = TRUE)
ssgsea_df <- as.data.frame(t(as.data.frame(ssgsea_data)))
ssgsea_df$Group <- filtered_mRNA_df$Group  # 添加分组信息

library(dplyr)

# 计算每个通路在 High 和 Low 组的中位数
median_scores <- ssgsea_df %>%
  group_by(Group) %>%
  summarise(across(everything(), median, na.rm = TRUE))

# 转换为数据框并去掉 Group 列
heatmap_data <- as.data.frame(median_scores)
rownames(heatmap_data) <- heatmap_data$Group  # 将 Group 列设置为行名
heatmap_data <- heatmap_data[,-1]  # 去掉 Group 列

p_values <- sapply(1:ncol(heatmap_data), function(i) {
  t.test(ssgsea_df[ssgsea_df$Group == "High", i],
         ssgsea_df[ssgsea_df$Group == "Low", i])$p.value
})
# p_values <- sapply(1:ncol(heatmap_data), function(i) {
#   wilcox.test(ssgsea_df[ssgsea_df$Group == "High", i],
#               ssgsea_df[ssgsea_df$Group == "Low", i])$p.value
# })

# 转换为显著性标记
significant_labels <- ifelse(p_values < 0.05, "*", "")
significant_df <- data.frame(Significance = significant_labels)
rownames(significant_df) <- colnames(heatmap_data)  # 设置行名

library(pheatmap)

# 创建红蓝配色的调色板
color_palette <- colorRampPalette(c("#3474ac", "white", "#e40414"))(100)

# 绘制热图
pheatmap(heatmap_data,
         cluster_rows = FALSE,
         cluster_cols = TRUE,
         scale = "row",
         fontsize = 10,
         show_rownames = TRUE,
         show_colnames = TRUE,
         main = "Median ssGSEA Scores by Group",
         color = color_palette
         # annotation_col = significant_df
)

pheatmap(heatmap_data,
         cluster_rows = FALSE,
         cluster_cols = TRUE,
         scale = "row",
         fontsize = 10,
         show_rownames = TRUE,
         show_colnames = TRUE,
         main = "Median ssGSEA Scores by Group",
         annotation_col = significant_df
         )  # 添加显著性注释



################################################################################
# TIDE
TIDE <- read.csv("./immune_check/exp_TIDE.csv",row.names = 1)

identical(filtered_mRNA_df$case_id,rownames(TIDE))
# 先确保 filtered_mRNA_df$case_id 和 TIDE 的行名是相同的（即元素一致但顺序可能不同）
if (all(filtered_mRNA_df$case_id %in% rownames(TIDE))) {
  # 按照 filtered_mRNA_df$case_id 的顺序重排 TIDE 的行名
  TIDE_sorted <- TIDE[match(filtered_mRNA_df$case_id, rownames(TIDE)), ]
} else {
  stop("case_id 和 TIDE 行名不匹配")
}

# 检查是否排序成功
identical(filtered_mRNA_df$case_id, rownames(TIDE_sorted))  # 应该返回 TRUE

meta_TIDE <- cbind(filtered_mRNA_df$Group,TIDE)
# 设置合并后的第一列列名为 "Group"
colnames(meta_TIDE)[1] <- "cluster"
write.csv(meta_TIDE,"./immune_check/meta_TIDE.csv")

res = meta_TIDE
# 统计学检验一下
table(res$Responder,res$cluster)
f = fisher.test(table(res$cluster,res$Responder))
label = paste("fisher.test p value =",round(f$p.value,3))
label

# 网站中的柱状图
library(ggplot2)
library(dplyr)
res = arrange(res,desc(TIDE))
p1 = ggplot(res,aes(x = 1:nrow(res),
                    y = TIDE,
                    fill = Responder))+
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c("#e04030","#6cb8d2"))+
  xlab("patient")+
  ylab("TIDE value")+
  annotate("text", x = 40, y = -2, label = label,size = 10)+
  theme_bw()+
  theme(legend.position = "none") # 把P1的图注去掉了

########免疫反应与亚型
library(dplyr)
dat=count(res,cluster,Responder)
dat=dat%>%group_by(cluster)%>%
  summarise(Responder=Responder,n=n/sum(n))
dat$Responder=factor(dat$Responder,levels=c("False","True"))
dat

library(ggplot2)
p2=ggplot(data=dat)+
  geom_bar(aes(x=cluster,y=n,
               fill=Responder),
           stat="identity")+
  scale_fill_manual(values=c("#e04030","#6cb8d2"))+
  geom_label(aes(x=cluster,y=n,
                 label=scales::percent(n),
                 fill=Responder),
             color="white",
             size=4,label.size=0,
             show.legend = FALSE,
             position=position_fill(vjust=0.5))+
  ylab("Percentage")+
  theme_minimal()+
  guides(fill = guide_legend(title = "Responder"))  # 仅保留一个图例

library(patchwork)
p1+p2+plot_layout(widths=c(3,2),guides="collect")
ggsave('./immune_check/immune_responder.pdf',width = 20,height = 8)


res$cluster <- factor(res$cluster,levels = c("C1","C2"))
png("TIDE_web.png",width = 1200,height = 1200,res = 300)
ggplot(data=res,aes(x=cluster,y=TIDE,colour = cluster))+ #fill参数不要设置，会不好看
  geom_violin(#color = 'grey',
    alpha = 0.8, #alpha = 0.8 参数控制着小提琴图的透明度。
    scale = 'width',#小提琴宽度
    #linewidth = 1, #外轮廓粗细
    trim = TRUE)+ # trim = TRUE 参数控制着小提琴图的形状。
  geom_boxplot(mapping=aes(x=cluster,y=TIDE,colour=cluster,fill=cluster), #箱线图
               alpha = 0.5,
               size=1.5,
               width = 0.3)+
  geom_jitter(mapping=aes(x=cluster,y=TIDE,colour=cluster), #散点
              alpha = 0.3,size=3)+
  scale_fill_manual(limits=c("C1","C2"),
                    values =c("#e04030","#6cb8d2"))+
  scale_color_manual(limits=c("C1","C2"),
                     values=c("#e04030","#6cb8d2"))+ #颜色
  geom_signif(mapping=aes(x=cluster,y=TIDE), # 不同组别的显著性
              comparisons = list(c("C1","C2")), # 哪些组进行比较
              map_signif_level=T, # T显示显著性，F显示p value
              tip_length=c(0,0),#把向下的帽子去掉，分组数乘以2
              y_position = c(3.5), # 设置显著性线的位置高度
              size=1, # 修改线的粗细
              textsize = 4, # 修改显著性标记的大小
              test = "wilcox.test", # 检验的类型,可以更改
              color = "black")+ # 设置显著性线的颜色
  theme_bw()+ #设置白色背景
  guides(fill = guide_legend(title = "cluster"),  # 设置填充图例的标题
         color = guide_legend(title = "cluster"))+  # 设置颜色图例的标题
  labs(title = "",  # 设置标题
       x="",y= "TIDE value") # 添加标题，x轴，y轴标签
dev.off()


# library(tinyarray)
# res$cluster <- factor(res$cluster,levels = c("C1","C2"))
# colnames(res)
# dat <- t(res[,c(12,14:16,18:22)])
# head(dat)[1:4,1:4]
# draw_boxplot(dat,res$cluster)+
#   facet_wrap(~rows,scales ="free") +
#   scale_fill_manual(values = c("C1" = "#e04030", "C2" = "#6cb8d2"))
# ggsave("together.png",width = 12,height = 12)

#######################################
# 免疫检查点

genes <- c("PDCD1", "CD274", "CTLA4", "LAG3", "HAVCR2", "TIGIT", "VTCN1")
data = expression_matrix[genes,]
colnames(data) = stringr::str_sub(colnames(data),1,12)
data = as.data.frame(t(data))

# metabo = rt[,7:ncol(rt)]
library(limma)
library(corrplot)
#相关性矩阵
M=cor(data)
res1=cor.mtest(data, conf.level=0.95)      #显著性标记的矩阵

#绘制相关性图形
pdf(file="./immune_check/geneCor.pdf", width=7, height=7)
corrplot(M,
         order="original",
         method = "circle",
         type = "upper",
         tl.cex=1, pch=T,
         tl.col="black",
         p.mat = res1$p,
         insig = "label_sig",
         pch.cex = 1.6,
         sig.level=0.05,
         number.cex = 1,
         col=colorRampPalette(c("#3474ac", "white", "#e40414"))(50))
dev.off()
rownames(TCGA_CESC_OS)= TCGA_CESC_OS$case_id
x = TCGA_CESC_OS[rownames(data),]$risk
outTab=data.frame()
for(i in c(1:(ncol(data)))){
  # if(i==geneName){next}
  y=as.numeric(data[,i])
  corT=cor.test(x, y, method = 'spearman')
  cor=corT$estimate
  pvalue=corT$p.value
  if(pvalue<1){
    outTab=rbind(outTab, cbind(Query=colnames(data)[i], Gene=i, cor, pvalue))
  }
}

#输出相关性结果文件
write.table(file="./immune_check/corResult.txt", outTab, sep="\t", quote=F, row.names=F)

#相关性矩阵
data=data[,outTab$Query]
data$riskscore = rt$riskScore
M=cor(data)
res1=cor.mtest(data, conf.level=0.95)
#绘制相关性图形

pdf(file="./immune_check/corpot.pdf",width=5,height=5)
# corrplot(M,
#          order="original",
#          method = "color",
#          number.cex = 0.7,
#          type = "upper",
#          addCoef.col = "black",
#          diag = TRUE,
#          tl.col="black",
#          col=colorRampPalette(c("blue", "white", "red"))(50))
corrplot(M,
         order="original",
         method = "circle",
         type = "upper",
         tl.cex=1, pch=T,
         tl.col="black",
         p.mat = res1$p,
         insig = "label_sig",
         pch.cex = 1.6,
         sig.level=0.05,
         # number.cex = 1,
         col=colorRampPalette(c( "#3474ac", "white", "#e40414"))(50))
dev.off()

################################################################################
# immune check
library(dplyr)
library(tidyverse)
load("./rdata/immune_methods.rda")
frac2$RMSE_CIBERSORT = NULL
frac2$`P-value_CIBERSORT`=NULL
frac2$Correlation_CIBERSORT = NULL
frac2 = frac2 %>%
  column_to_rownames(var = "ID")
frac2 = t(frac2)
rownames(frac2) <- gsub("_CIBERSORT$", "", rownames(frac2))
# frac1$Group = "TIMER"
# frac2$Group = "CiberSort"
# frac3$Group = "MCPCounter"
# frac4$Group = "xCell"
# frac5$Group = "IPS"
# frac6$Group = "epic"
# frac7$Group = "ESTIMATE"
# frac8$Group = "ABIS"
# frac9$Group = "ConsensusTME"
# frac10$Group = "quanTIseq"
frac_list = list(frac1,
                 frac2,
                 frac3,
                 frac4,
                 frac5,
                 frac6,
                 frac7,
                 frac8,
                 frac9,
                 frac10)
names(frac_list) = c("TIMER",
                     "CiberSort",
                     "MCPCounter",
                     "xCell",
                     "IPS",
                     "epic",
                     "ESTIMATE",
                     "ABIS",
                     "ConsensusTME",
                     "quanTIseq")

library(dplyr)
library(survminer)
library(survival)
library(readxl)
library(dplyr)
TCGA_CESC_OS <- read_excel("dl_output/result_tcga_OS_final_12char.xlsx")
TCGA_CESC_OS$case_id = substr(TCGA_CESC_OS$case_id,1,12)
head(TCGA_CESC_OS)
load("./rdata/cli_combined_selected.rda")
# Identify IDs in FUDAN_CESC_OS that are not in combined_Data
missing_ids <- TCGA_CESC_OS %>%
  filter(!case_id %in% os_pfs_combined$ID)

# Print the missing IDs
print(missing_ids)


# Remove rows with missing IDs from FUDAN_CESC_OS
TCGA_CESC_OS <- TCGA_CESC_OS %>%
  filter(case_id %in% os_pfs_combined$ID)
TCGA_CESC_OS$status = abs(TCGA_CESC_OS$censorship-1)
TCGA_CESC_OS$Group=ifelse(TCGA_CESC_OS$risk > median(TCGA_CESC_OS$risk),'High','Low')
head(TCGA_CESC_OS)

# 获取 TCGA_CESC_OS 中的 case_id 前12个字符
case_ids <- substr(TCGA_CESC_OS$case_id, 1, 12)


# 处理 frac_list
frac_list_filtered <- lapply(frac_list, function(df) {


  # 提取列名的前12个字符
  colnames(df) <- substr(colnames(df), 1, 12)
  # 去重列名，保留第一项
  unique_colnames <- colnames(df)[!duplicated(colnames(df))]
  print(length(unique_colnames))
  df <- df[, unique_colnames, drop = FALSE]

  # 只保留在 TCGA_CESC_OS 中存在的 case_id 列
  df_filtered <- df[, colnames(df) %in% case_ids, drop = FALSE]

  return(df_filtered)
})

# 获取 frac_list_filtered 中所有列名
filtered_case_ids <- unique(unlist(lapply(frac_list_filtered, colnames)))

# 筛选 TCGA_CESC_OS 中的行，只保留在 filtered_case_ids 中的 case_id
TCGA_CESC_OS_filtered <- TCGA_CESC_OS[TCGA_CESC_OS$case_id %in% filtered_case_ids, ]

# 确保 case_id 也只保留前12个字符
TCGA_CESC_OS_filtered$case_id <- substr(TCGA_CESC_OS_filtered$case_id, 1, 12)

# 清空环境
# rm(list=ls())

# 加载必要的库
library(tcltk)
library(tidyverse)
library(ggsci)
library(scales)
library(ggtext)
library(reshape2)
library(ggplot2)
library(ggpubr)

# 加载数据
# 确保 TCGA_CESC_OS_filtered 和 frac_list_filtered 已经定义
# TCGA_CESC_OS_filtered 应包含 'risk' 列
# frac_list_filtered 是一个包含免疫细胞比例的数据框列表

# 获取风险得分
gene_exp = TCGA_CESC_OS_filtered$risk
gene_exp = as.data.frame(gene_exp)
rownames(gene_exp) = TCGA_CESC_OS_filtered$case_id
colnames(gene_exp) = "riskScore"

# 结果表格初始化
outTab = data.frame()

# 计算相关性
for (software in names(frac_list_filtered)) {
  # software = names(frac_list_filtered)[7]
  print(software)
  data = frac_list_filtered[[software]]

  # 确保列名与 TCGA_CESC_OS_filtered 的 case_id 对应
  common_ids = intersect(colnames(data), rownames(gene_exp))
  data = t(data)
  data_filtered = data[common_ids, , drop = FALSE]
  gene_exp_filtered = gene_exp[common_ids, , drop = FALSE]

  for (i in colnames(data_filtered)) {
    # i = colnames(data_filtered)[4]
    y = as.numeric(data_filtered[, i])
    if (any(is.na(y))) { next }
    if (sd(y) < 0.0001) { next }
    corT = cor.test(gene_exp_filtered$riskScore, y, method = "spearman")
    cor = corT$estimate
    pvalue = corT$p.value
    if (pvalue < 0.05) {
      outTab = rbind(outTab, cbind(immune = i, cor, pvalue, Software = software))
    }
  }
}

# 输出相关性结果
write.table(file = "csv/immune_corResult.txt", outTab, sep = "\t", quote = FALSE, row.names = FALSE)

ov_palette =c('#A499CC',
              '#E0A7C8',
              '#E069A6',
              "#f1707d",
              "#AFC2D9",
              "#6894B9",
              "#79B99D",
              "#F5D2A8",
              "#D2EBC8")
# 绘制气泡图
corResult = read.table("csv/immune_corResult.txt", head = TRUE, sep = "\t")
corResult$Software = factor(corResult$Software, levels = unique(corResult$Software))
b = corResult[order(corResult$Software),]
b$immune = factor(b$immune, levels = rev(unique(b$immune)))
colslabels = rep(hue_pal()(length(levels(b$Software))), table(b$Software))
library(ggplot2)
ggplot(data = b, aes(x = cor, y = immune, color = Software)) +
  labs(x = "Correlation coefficient", y = "Immune cell") +
  geom_point(size = 4.1)+
  scale_color_manual(values = ov_palette)+
  geom_linerange(aes(xmin = 0, xmax = cor), size = 3, alpha = 0.5) +  # 添加线范围
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )

pdf(file = "correlation.pdf", width = 7, height = 10)
ggplot(data = b, aes(x = cor, y = immune, color = Software)) +
  labs(x = "Correlation coefficient", y = "Immune cell") +
  geom_point(size = 4.1)
# theme(panel.background = element_rect(fill = "white", size = 1, color = "black"),
#       panel.grid = element_line(color = "grey75", size = 0.5),
#       axis.ticks = element_line(size = 0.5),
#       axis.text.y = ggtext::element_markdown(colour = rev(colslabels)))
dev.off()


################################################################################
# oncopredict
library(dplyr)
library(survminer)
library(survival)
library(readxl)
library(dplyr)
TCGA_CESC_OS <- read_excel("dl_output/TCGA-CESC-OS.xlsx")
TCGA_CESC_OS$case_id = substr(TCGA_CESC_OS$case_id,1,12)
head(TCGA_CESC_OS)
load("./rdata/cli_combined_selected.rda")
# Identify IDs in FUDAN_CESC_OS that are not in combined_Data
missing_ids <- TCGA_CESC_OS %>%
  filter(!case_id %in% os_pfs_combined$ID)

# Print the missing IDs
print(missing_ids)


# Remove rows with missing IDs from FUDAN_CESC_OS
TCGA_CESC_OS <- TCGA_CESC_OS %>%
  filter(case_id %in% os_pfs_combined$ID)
TCGA_CESC_OS$status = abs(TCGA_CESC_OS$censorship-1)
TCGA_CESC_OS$Group=ifelse(TCGA_CESC_OS$risk > median(TCGA_CESC_OS$risk),'High','Low')
head(TCGA_CESC_OS)

mRNA = data.table::fread("./GDCdata/TCGA-CESC/HiSeqV2_PANCAN")
library(dplyr)
library(tidyverse)
mRNA <- mRNA %>%
  column_to_rownames("sample")
colnames(mRNA) = substr(colnames(mRNA),1,12)



# 首先根据 TCGA_CESC_OS 中的 case_id 过滤 cnv 的列名
filtered_mRNA <- mRNA[, colnames(mRNA) %in% TCGA_CESC_OS$case_id]

# 将 case_id 添加到 filtered_cnv 中，以便后续合并
filtered_mRNA_df <- as.data.frame(t(filtered_mRNA))
filtered_mRNA_df$case_id <- rownames(filtered_mRNA_df)

# 将 TCGA_CESC_OS 中的 Group 信息添加到 filtered_cnv_df
filtered_mRNA_df <- merge(filtered_mRNA_df, TCGA_CESC_OS[, c("case_id", "Group")], by.x = "case_id", by.y = "case_id")
# head(filtered_mRNA_df[,c(20520:20532)])
rownames(filtered_mRNA_df) = filtered_mRNA_df$case_id
expression_matrix <- filtered_mRNA_df  # 去掉 Group 列

expression_matrix$case_id =NULL
expression_matrix$Group = NULL
expression_matrix = t(expression_matrix)
library(oncoPredict)
GDSC2_Expr=readRDS(file='./oncopredict/GDSC2_Expr (RMA Normalized and Log Transformed).rds')
GDSC2_Res=readRDS(file = './oncopredict/GDSC2_Res.rds')
GDSC2_Res=exp(GDSC2_Res)

# oncopredict
calcPhenotype(trainingExprData = GDSC2_Expr,    #train
              trainingPtype = GDSC2_Res,        #train
              testExprData = expression_matrix,              #test
              batchCorrect = 'eb',              #
              powerTransformPhenotype = TRUE,
              removeLowVaryingGenes = 0.2,      #
              minNumSamples = 10,               #
              printOutput = TRUE,               #
              removeLowVaringGenesFrom = 'rawData')

# load data
rm(list = ls())
gc()
#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("limma")

#install.packages("ggplot2")
#install.packages("ggpubr")



library(limma)
library(ggplot2)
library(ggpubr)

pFilter=0.001                      #

drugFile="./calcPhenotype_Output/DrugPredictions.csv"     #
# load("./output_risk.rda")
# rt$Group=ifelse(rt$riskScore > median(rt$riskScore),'High','Low')
# rt$risk=factor(rt$Group,levels = c('Low','High'))
#

senstivity=read.csv(drugFile, header=T, sep=",", check.names=F, row.names=1)
# colnames(senstivity)=gsub("(.*)\\_(\\d+)", "\\1", colnames(senstivity))
rownames(senstivity) = stringr::str_sub(rownames(senstivity),1,12)
#
# sameSample=intersect(row.names(rt), row.names(senstivity))
# risk = rt
identical(filtered_mRNA_df$case_id,rownames(senstivity))
rownames(TCGA_CESC_OS)= TCGA_CESC_OS$case_id
risk = TCGA_CESC_OS[filtered_mRNA_df$case_id,]$Group
# risk=filtered_mRNA_df$risk
# senstivity=senstivity[sameSample,,drop=F]
rt=cbind(risk, senstivity)

#
# rt$risk=factor(rt$risk, levels=c("low", "high"))
type=levels(factor(rt[,"risk"]))
comp=combn(type, 2)
my_comparisons=list()
for(i in 1:ncol(comp)){my_comparisons[[i]]<-comp[,i]}

#��ҩ�����ѭ��, ��������ͼ
for(drug in colnames(rt)[2:ncol(rt)]){
  rt1=rt[,c(drug, "risk")]
  colnames(rt1)=c("Drug", "Risk")
  rt1=na.omit(rt1)
  rt1$Drug=log2(rt1$Drug+1)
  #�������
  test=wilcox.test(Drug ~ Risk, data=rt1)
  diffPvalue=test$p.value
  #������������ҩ���������ͼ
  if(diffPvalue<0.05){
    boxplot=ggboxplot(rt1, x="Risk", y="Drug", fill="Risk",
                      xlab="Risk",
                      ylab=paste0(drug, " senstivity"),
                      legend.title="Risk",
                      palette=c("#AFC2D9", "#E0A7C8")
    )+
      stat_compare_means(comparisons=my_comparisons)
    #���ͼ��
    pdf(file=paste0("./drug/drugSenstivity.", drug, ".pdf"), width=5, height=4.5)
    print(boxplot)
    dev.off()
  }
}

############################################################################
