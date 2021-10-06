# 2.1 安装CRAN来源常用包
site="https://mirrors.tuna.tsinghua.edu.cn/CRAN"
# 依赖包列表：参数解析、数据变换、绘图和开发包安装、安装依赖、ggplot主题
package_list = c("ggplot2","dplyr")
# 判断R包加载是否成功来决定是否安装后再加载
for(p in package_list){
	if(!suppressWarnings(suppressMessages(require(p, character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)))){
		install.packages(p, repos=site)
		suppressWarnings(suppressMessages(library(p, character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)))
  }
}


# Set ggplot2 drawing parameter, such as axis line and text size, lengend and title size, and so on.
main_theme = theme(panel.background=element_blank(), panel.grid=element_blank(),
	axis.line.x=element_line(size=.5, colour="black"), axis.line.y=element_line(size=.5, colour="black"),
	axis.ticks=element_line(color="black"), axis.text=element_text(color="black", size=7),
	legend.position="right", legend.background=element_blank(), legend.key=element_blank(), legend.text= element_text(size=7),
	text=element_text(family="sans", size=7))

# 实验差异比较结果
x = read.table("result/compare/KO-WT.txt", header=T, row.names= 1, sep="\t", stringsAsFactors = F)
# 只提取前14列
x = x[,1:7]
x = na.omit(x)

# 添加物种注释
taxonomy = read.table("result/taxonomy.txt", sep = "\t", row.names=1, header=T, stringsAsFactors = F)
taxonomy = taxonomy[rownames(x),]
x = cbind(x, taxonomy)

# P值求负对数
x$neglogp = -log10(x$PValue)

x$otu=rownames(x)
x = arrange(x, Kingdom, Phylum, Class, Order, Family, Genus, Species)
x$otu = factor(x$otu, levels=x$otu)   # set x order
x$num = 1:dim(x)[1]

# 读取高丰度门，用于着色
per= read.delim("result/tax/sum_c.txt", sep = "\t", row.names=1, header=T)
mean = rowMeans(per)
per = as.data.frame(mean[order(mean, decreasing = T)])
top_tax=head(rownames(per), n=10)

# 将低丰度的门变为Low Abundance
x$tax = x$Class

# 将门中 proteobacteria 替换为纲
# x[x$tax %in% "Proteobacteria",]$tax = x[x$tax %in% "Proteobacteria",]$Class # no level can get value

if (length(unique(x$tax)) > length(top_tax)){
	x[!(x$tax %in% top_tax),]$tax = "Low Abundance" # no level can get value
}

# 调置标签顺序
label = unique(x$tax)
label = label[!(label %in% "Low Abundance")] # 删除低丰度
x$tax = factor(x$tax, levels = c(label, "Low Abundance"))
# 计算标签中位数
temp = x[x$tax %in% label, c("tax","num")]
mat_mean = aggregate(temp[,-1], by=temp[1], FUN=mean) # mean


# 调整过大的负对数
if (max(x$neglogp)>20){
  x[x$neglogp>20,]$neglogp = 20
}
colnames(x)[2]="log2CPM"
# Manhattan plot
FDR = min(x$neglogp[x$level!="NotSig"])
p = ggplot(x, aes(x=num, y=neglogp, color=tax, size=log2CPM, shape=level)) +
  geom_point(alpha=.7) +
  geom_hline(yintercept=FDR, linetype=2, color="lightgrey") +
  scale_shape_manual(values=c(25, 17, 20))+
  scale_size(breaks=c(5, 10, 15)) +
  labs(x="Phylum of Features", y="-log10(P)") +
  labs(title=paste(gsub(".txt", "", basename("result/compare/KO-WT.txt")), sep=" ")) +
  main_theme +
  # theme(legend.position="top") +
  scale_x_continuous(breaks=mat_mean$x, labels=mat_mean$tax) +
  theme(axis.text.x=element_text(angle=45, hjust=1, vjust=1)) # + ylim(0,20)
  # ylim保持多组y轴范围一致，一般在最终出版时才使用
p
ggsave(file=paste("result/compare/KO-WT.manhattan.c.legend.pdf", sep=""), p, width = 183, height = 149, useDingbats=F, units = "mm")
# ggsave(file=paste("result/compare/KO-WT.txt_manhattan.png", sep=""), p, width = 183, height = 149)

