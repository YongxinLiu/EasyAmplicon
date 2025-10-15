[TOC]

# 易扩增子 EasyAmplicon2

    # 作者 Authors: 刘永鑫(Yong-Xin Liu), 陈同(Tong Chen)等
    # 版本 Version: v2.0
    # 更新 Update: 2025-5-30
    # 系统要求 System requirement: WIndows 10+ / Mac OS 10.12+ / Ubuntu 20.04+
    # 引文 Reference: Liu, et al. 2023. EasyAmplicon: An easy-to-use, open-source, reproducible, and community-based pipeline for amplicon data analysis in microbiome research. iMeta 2: e83. https://doi.org/10.1002/imt2.83
    # Yousuf, et al. 2024. Unveiling microbial communities with EasyAmplicon: A user-centric guide to perform amplicon sequencing data analysis. iMetaOmics 1: e42. https://doi.org/10.1002/imo2.42

    # 设置工作(work directory, wd)和软件数据库(database, db)目录
    # 添加环境变量，并进入工作目录 Add environmental variables and enter work directory
    # **每次打开Rstudio必须运行下面4行 Run it**，可选替换${db}为EasyMicrobiome安装位置
    wd=/d/EasyAmplicon2
    db=/d/EasyMicrobiome
    PATH=$PATH:${db}/win
    cd ${wd}


# R语言多样性和物种组成分析

## 1. Alpha多样性

### 1.1 Alpha多样性箱线图

    # 查看帮助,后续所有R代码均可用此方法查看帮助信息
    cd ${wd}/PacBio
    Rscript ${db}/script/alpha_div_box.R -h
    # 完整参数，多样性指数可选richness chao1 ACE shannon simpson invsimpson
    Rscript ${db}/script/alpha_div_box.R --input result/alpha/vegan.txt \
    --metadata result/metadata.txt \
    --alpha_index richness,chao1,ACE,shannon,simpson,invsimpson \
    --group Group \
    --out_prefix result/alpha/vegan
    #--alpha_index控制绘制多样性指数类型，可绘制单个
    

### 1.2 稀释曲线
    
    # PacBio data
    Rscript ${db}/script/alpha_rare_curve.R \
      --input result/alpha/alpha_rare.txt --design result/metadata.txt \
      --group Group --output result/alpha/ \
      --width 120 --height 78


### 1.3 多样性维恩图

    # PacBio data
    # 支持2-4组比较
    Rscript ${db}/script/venn.R \
      --input result/alpha/otu_group_exist.txt \
      --groups All,feces,plaque,saliva \
      --output result/alpha/venn.pdf
    # 在线绘制韦恩图推荐在线网站http://www.ehbio.com/test/venn


## 2. Beta多样性

### 2.1 距离矩阵热图pheatmap
    
    # PacBio data
    # 以bray_curtis为例，-f输入文件,-h是否聚类TRUE/FALSE,-u/v为宽高英寸
    bash ${db}/script/sp_pheatmap.sh \
      -f result/beta/bray_curtis.txt \
      -H 'TRUE' -u 6 -v 5
    # 添加分组注释，如2，4列的基因型和地点
    mkdir -p temp
    cut -f 1-2 result/metadata.txt > temp/group.txt
    # -P添加行注释文件，-Q添加列注释
    bash ${db}/script/sp_pheatmap.sh \
      -f result/beta/bray_curtis.txt \
      -H 'TRUE' -u 6.9 -v 5.6 \
      -P temp/group.txt -Q temp/group.txt
    # 距离矩阵与相关类似，可尝试corrplot或ggcorrplot绘制更多样式
    # - [绘图相关系数矩阵corrplot](http://mp.weixin.qq.com/s/H4_2_vb2w_njxPziDzV4HQ)
    # - [相关矩阵可视化ggcorrplot](http://mp.weixin.qq.com/s/AEfPqWO3S0mRnDZ_Ws9fnw)


### 2.2 主坐标分析PCoA

    # PacBio data
    # 输入文件，选择分组，输出文件
    Rscript ${db}/script/beta_PCoA.R \
      --input result/otutab.txt \
      --metadata result/metadata.txt \
      --group Group \
      --output result/beta/PCoa2.pdf


### 2.3 限制性主坐标分析CPCoA

    # PacBio data
    Rscript ${db}/script/beta_cpcoa.R \
      --input result/beta/bray_curtis.txt --design result/metadata.txt \
      --group Group --output result/beta/bray_curtis.cpcoa.pdf \
      --width 89 --height 59
    # 添加样本标签 --label TRUE
    Rscript ${db}/script/beta_cpcoa.R \
      --input result/beta/bray_curtis.txt --design result/metadata.txt \
      --group Group --label TRUE --width 89 --height 59 \
      --output result/beta/bray_curtis.cpcoa.label.pdf


## 3. 物种组成Taxonomy

### 3.1 堆叠柱状图Stackplot

    # PacBio data
    # 以门(p)水平为例，结果包括output.sample/group.pdf两个文件
    Rscript ${db}/script/tax_stackplot.R \
      --input result/tax/sum_p.txt --design result/metadata.txt \
      --group Group -t 10 --color manual1 --legend 7 --width 89 --height 59 \
      --output result/tax/sum_p.stackplot
    # 补充-t 筛选丰度前多少物种及-s x轴排列位置-s "feces,plaque,saliva" 
    # 修改颜色--color ggplot, manual1(30), Paired(12) or Set3(12)

    # 批量绘制输入包括p/c/o/f/g共5级
    for i in p c o f g s; do
    Rscript ${db}/script/tax_stackplot.R \
      --input result/tax/sum_${i}.txt --design result/metadata.txt \
      --group Group -t 10 --output result/tax/sum_${i}.stackplot \
      --legend 8 --width 89 --height 59; done
      
    # 组间比较堆叠柱状图
    Rscript ${db}/script/tax_stack_compare.R \
      --input result/tax/data_practice53_16S.txt \
      --compare Illumina-PacBio \
      --output result/tax/
      
    # 连线堆叠柱状图
    Rscript ${db}/script/microeco_alluvial.R \
      --otu result/otutab2.txt \
      --metadata result/metadata.txt \
      --taxonomy result/taxonomy.txt \
      --output result/tax/


### 3.2 弦/圈图circlize

    # PacBio data
    # 以纲(class,c)为例，绘制前5组
    i=c
    Rscript ${db}/script/tax_circlize.R \
      --input result/tax/sum_${i}.txt --design result/metadata.txt \
      --group Group --legend 5
    # 结果位于当前目录circlize.pdf(随机颜色)，circlize_legend.pdf(指定颜色+图例)
    # 移动并改名与分类级一致
    mv circlize.pdf result/tax/sum_${i}.circlize.pdf
    mv circlize_legend.pdf result/tax/sum_${i}.circlize_legend.pdf
    

### 3.3 气泡图

    # PacBio data
    # 以属为例（genus，g），绘制丰度前15的属水平丰度气泡图；输入物种丰度表和样本metadata文件，输出分组物种丰度气泡图
    i=g
    Rscript ${db}/script/tax_bubble.R \
    -i result/tax/sum_${i}.txt \
    -g result/metadata.txt \
    -c Group -w 7 -e 4 -s 15 -n 15 \
    -o result/tax/sum_g.bubble3.pdf
    

### 3.4 甜甜圈图(Donut)和雷达图(Ladar plot)展示相对丰度组成

    # PacBio data
    Rscript ${db}/script/Donut_plot.R \
      --otu_table result/otutab.txt \
      --metadata result/metadata.txt \
      --taxonomy result/taxonomy.txt \
      --output_dir result/tax/
    

### 3.5 核心微生物(至少在80%样本中存在)

    # PacBio data
    Rscript ${db}/script/core_ASVS_and_other.R \
      --otu_table result/otutab.txt \
      --metadata result/metadata.txt \
      --taxonomy result/taxonomy.txt \
      --output_dir result/tax/
    
    # 核心微生物散点图
    Rscript ${db}/script/Core_Abundance_ScatterPlot.R \
      --otu_table result/otutab.txt \
      --metadata result/metadata.txt \
      --output result/tax/


# 4、差异比较

## 1. R语言差异分析

### 1.1 差异比较

    # PacBio data
    mkdir -p result/compare/
    # 输入特征表、元数据；指定分组列名、比较组和丰度
    # 选择方法 wilcox/t.test/edgeR、pvalue和fdr和输出目录
    # 这里选择默认的wilcox
    compare="saliva-plaque"
    Rscript ${db}/script/compare.R \
      --input result/otutab.txt --design result/metadata.txt \
      --group Group --compare ${compare} --threshold 0.01 \
      --pvalue 0.05 --fdr 0.2 \
      --output result/compare/
    # 并筛选丰度为前20%的ASV（--threshold 0.05 ）用来画热图
    compare="feces-saliva"
    Rscript ${db}/script/compare.R \
      --input result/otutab.txt --design result/metadata.txt \
      --group Group --compare ${compare} --threshold 0.2 \
      --pvalue 0.05 --fdr 0.2 \
      --output result/compare/
    # 常见错误：Error in file(file, ifelse(append, "a", "w")) : 无法打开链结 Calls: write.table -> file
    # 解决方法：输出目录不存在，创建目录即可
    

### 1.2 火山图

    # PacBio data
    Rscript ${db}/script/volcano2.R \
      --input result/compare/saliva_plaque2.txt \
      --group saliva-plaque \
      --output result/compare/
      
    # 多组比较火山图
    Rscript ${db}/script/multigroup_compare_volcano.R \
      --input result/otutab2.txt \
      --metadata result/metadata.txt \
      --output result/compare/


### 1.3 热图

    # PacBio data
    # 输入compare.R的结果，筛选列数，指定元数据和分组、物种注释，图大小英寸和字号
    compare="saliva_plaque2"
    bash ${db}/script/compare_heatmap.sh \
       -i result/compare/${compare}.txt -l 7 \
       -d result/metadata.txt -A Group \
       -t result/taxonomy.txt \
       -w 12 -h 20 -s 14 \
       -o result/compare/${compare}
      
    # 多样本分组比较热图
    Rscript ${db}/script/multisample_compare_heatmap.R \
      --input result/compare/data5.txt \
      --output result/compare/
    

### 1.4 三元图
    
    # PacBio data
    Rscript ${db}/script/Ternary_plot.R   \
      --input result/tax/sum_p.txt   \
      --metadata result/metadata.txt  \
      --group Group   \
      --taxlevel Phylum \
      --output result/compare/ternary_p.pdf   \
      --topn 10
    
    # Error: package or namespace load failed for ‘ggtern’: .onLoad failed in loadNamespace() for 'ggtern', details: call: NULL error: <ggplot2::element_line> object properties are invalid: - @lineend must be <character> or <NULL>, not S3<arrow>
    # 遇到以上报错是因为ggplot2和ggtern软件包不匹配，需要将ggplot2和ggtern匹配到正确的版本
    # remove.packages("ggplot2")
    # install.packages("https://cran.r-project.org/src/contrib/Archive/ggplot2/ggplot2_3.4.4.tar.gz", repos = NULL, type = "source")
    # library(ggplot2)
    # remove.packages("ggtern")
    # install.packages("https://cran.r-project.org/src/contrib/Archive/ggtern/ggtern_3.4.2.tar.gz",repos = NULL, type = "source")
    # library(ggtern)


## 1.5. STAMP差异分析图

### 1.5.1 生成输入文件(备选)

    # PacBio data
    Rscript ${db}/script/format2stamp.R -h
    mkdir -p result/stamp
    Rscript ${db}/script/format2stamp.R --input result/otutab.txt \
      --taxonomy result/taxonomy.txt --threshold 0.01 \
      --output result/stamp/tax
    # 可选Rmd文档见result/format2stamp.Rmd
    
    
### 1.5.2 R语言绘制STAMP图

    # 请到https://github.com/YongxinLiu/EasyMicrobiome/tree/master/script下载最新的代码脚本compare_stamp.R

    # PacBio data
    compare="feces-plaque"  
    Rscript ${db}/script/compare_stamp.R \
      --input result/tax/sum_g2.txt --metadata result/metadata.txt \
      --group Group --compare ${compare} --threshold 0.1 \
      --method "wilcox" --pvalue 0.2 --fdr "none" \
      --width 150 --height 120 \
      --output result/tax/stamp_${compare}
    
    # 报错library(patchwork)载入报错，是因为ggplot2版本不匹配，这里需要对ggplot2软件包进行更新  
    # install.packages("ggplot2")
      

### Lefse分析及组间比较柱状图

    Lefse分析在线网站：https://www.bioincloud.tech/standalone-task-ui/lefse
    Reference: Gao, Yunyun, Guoxing Zhang, Shunyao Jiang, and Yong‐Xin Liu. 2024. “
    Wekemo Bioincloud: A User‐friendly Platform for Meta‐omics Data Analyses.” i
    Meta e175. https://doi.org/10.1002/imt2.17
 
 
### 2.环境因子关联分析

### Mantel检验相关性分析热图
 
    # PacBio data
    Rscript ${db}/script/env_mantel_heatmap.R \
      --input result/tax/otutab2.txt \
      --env result/tax/env_amplicon.txt \
      --output result/tax/

### Redundancy analysis (RDA) 冗余分析

    # PacBio data
    # 软件按照可能需要一些时间
    Rscript ${db}/script/RDA_microeco.R \
      --input result/tax/otutab.txt \
      --metadata result/tax/metadata.txt \
      --tax result/tax/taxonomy.txt \
      --phylo result/tax/otus.tree \
      --output result/tax/


### 3.系统发育树

### 3.1 无进化距离系统发育树

    # PacBio data
    Rscript ${db}/script/phylogenetic_tree01.R \
      --input result/tax/otus.nwk \
      --anno result/tax/annotation2.txt \
      --output result/tax/
      
### 3.2 有进化距离系统发育树
    # PacBio data
    Rscript ${db}/script/phylogenetic_tree02.R \
      --input result/tax/otus.nwk \
      --anno result/tax/annotation3.txt \
      --output result/tax/
    

### 4.网络分析

### 4.1 多组Spearman相关性网络比较分析 (Multi-group Spearman correlation network)
    
    # PacBio data
    Rscript ${db}/script/Spearman_network01.R \
      --input result/tax/otutab_amplicon.txt \
      --group result/tax/sample_amplicon.txt \
      --tax result/tax/taxonomy_amplicon.txt \
      --output result/tax/


### 5.随机森林模型

### 5.1 运行随机森林模型

    # PacBio data
    # 注意：这是第一版本的随机森林模型，其中的一些参数之间内置在代码中，使用时需要根据具体样本情况进行调整
    Rscript ${db}/script/random_forest01.R \
      --input result/tax/PacBio_data.txt \
      --group result/tax/PacBio_metadata.txt \
      --output result/tax/RF_model/
      
### 5.2 绘制柱状图
    
    # PacBio data
    Rscript ${db}/script/RF_plot01.R \
      --input result/tax/RF_model/Species_imp_rf21.txt \
      --tax result/tax/taxonomy_amplicon.txt \
      --optimal 4 \
      --output result/tax/RF_model/
    

### 6.功能差异分析

### 6.1 数据处理及差异分析
    
    # PacBio data
    Rscript ${db}/script/function_data_process.R \
      --input result/tax/KEGG.PathwayL2.raw.txt \
      --group result/tax/metadata_amplicon.txt \
      --output result/tax/pathway/
    
### 6.2 热图结合柱状图展示差异

    # PacBio data
    Rscript ${db}/script/function_diff_plot01.R \
      --input result/tax/pathway/Difference_pathway21.txt \
      --pathway result/tax/pathway/pathway_count_data.txt \
      --statistic result/tax/pathway/pathway_percent_abundance2.txt \
      --output result/tax/pathway/


### 注意事项

1.绘图代码中的部分参数需要根据特定研究进行修改；
2.可到EasyMicrobiome(https://github.com/YongxinLiu/EasyMicrobiome)文件夹下获取相应的R语言脚本，可自行对代码进行调整;
3.如有疑问欢迎提问交流(https://github.com/YongxinLiu/EasyAmplicon/pulls)。



