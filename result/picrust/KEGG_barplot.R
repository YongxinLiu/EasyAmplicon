# 读取文件
KEGG <- read.table("L2.txt",header = TRUE,sep = "\t")

library(ggplot2)
library(reshape2)

# KEGG第一层分类中名字特别长，需要自动换行(替换空格为换行\n)
swr = function(string, nwrap = 12){
  paste(strwrap(string,width = nwrap),collapse = "\n")
}
swr = Vectorize(swr)
KEGG$L1 <- swr(KEGG$L1)

# 绘制L2级和丰度的柱状图，按L1着色并分面
p <- ggplot(KEGG,aes(Abundance,L2)) +
  geom_bar(aes(fill = L1),stat = "identity",width = 0.6) +
  xlab("Relative abundance (%)") + 
  ylab("KEGG Pathway") +
  theme(panel.background = element_rect(fill = "white",colour='black'),
        panel.grid.major = element_line(color = "grey",linetype = "dotted",size = 0.3),
        panel.grid.minor = element_line(color = "grey",linetype = "dotted",size = 0.3),
        axis.ticks.length = unit(0.4,"lines"), 
        axis.ticks = element_line(color='black'),
        axis.line = element_line(colour = "black"), 
        axis.title.x=element_text(colour='black', size=8,face = "bold"),
        axis.title.y=element_text(colour='black', size=8),
        axis.text.x=element_text(colour='black',size=8),
        axis.text.y = element_text(color = "black",size = 8),
        legend.position = "none",
        strip.text.y = element_text(angle = 0,size = 12,face = "bold")) +
  facet_grid(L1~.,space = "free_y",scales = "free_y")
# 预览
p
# 保存矢量图
ggsave("L2.pdf", p, width = 183, height = 240, units = "mm")
