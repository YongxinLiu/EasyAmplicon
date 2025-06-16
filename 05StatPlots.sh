[TOC]

# 易扩增子EasyAmplicon2

    # 作者 Authors: 刘永鑫(Yong-Xin Liu), 陈同(Tong Chen)等
    # 版本 Version: v2.0
    # 更新 Update: 2025-5-30
    # 系统要求 System requirement: WIndows 10+ / Mac OS 10.12+ / Ubuntu 20.04+
    # 引文 Reference: Liu, et al. 2023. EasyAmplicon: An easy-to-use, open-source, reproducible, and community-based pipeline for amplicon data analysis in microbiome research. iMeta 2: e83. https://doi.org/10.1002/imt2.83
    # Yousuf, et al. 2024. Unveiling microbial communities with EasyAmplicon: A user-centric guide to perform amplicon sequencing data analysis. iMetaOmics 1: e42. https://doi.org/10.1002/imo2.42


```
    # 设置工作(work directory, wd)和软件数据库(database, db)目录
    # 添加环境变量，并进入工作目录 Add environmental variables and enter work directory
    # **每次打开Rstudio必须运行下面4行 Run it**，可选替换${db}为EasyMicrobiome安装位置
    wd=/c/amplicon2
    db=/c/EasyMicrobiome
    PATH=$PATH:${db}/win
    cd ${wd}
    
```

# R语言多样性和物种组成分析

## 1. Alpha多样性

### 1.1 Alpha多样性箱线图
    
```
    # PacBio data
    # 查看帮助,后续所有R代码均可用此方法查看帮助信息
    Rscript ${db}/script/alpha_div_box.R -h
    # 完整参数，多样性指数可选richness chao1 ACE shannon simpson invsimpson
    Rscript ${db}/script/alpha_div_box.R --input result2/alpha/vegan.txt \
    --metadata result2/metadata.txt \
    --alpha_index richness,chao1,ACE,shannon,simpson,invsimpson \
    --group Group \
    --out_prefix result2/alpha/vegan
    #--alpha_index控制绘制多样性指数类型，可绘制单个
    
    # Illumina data
    # 查看帮助,后续所有R代码均可用此方法查看帮助信息
    Rscript ${db}/script/alpha_div_box.R -h
    # 完整参数，多样性指数可选richness chao1 ACE shannon simpson invsimpson
    Rscript ${db}/script/alpha_div_box.R --input result_illumina/alpha/vegan.txt \
    --metadata result_illumina/metadata.txt \
    --alpha_index richness,chao1,ACE,shannon,simpson,invsimpson \
    --group Group \
    --out_prefix result_illumina/alpha/vegan
    #--alpha_index控制绘制多样性指数类型，可绘制单个
    
    # Nanopore data
    # 查看帮助,后续所有R代码均可用此方法查看帮助信息
    Rscript ${db}/script/alpha_div_box.R -h
    # 完整参数，多样性指数可选richness chao1 ACE shannon simpson invsimpson
    Rscript ${db}/script/alpha_div_box.R --input result_nanopore/alpha/vegan.txt \
    --metadata result_nanopore/metadata.txt \
    --alpha_index richness,chao1,ACE,shannon,simpson,invsimpson \
    --group Group \
    --out_prefix result_nanopore/alpha/vegan
    #--alpha_index控制绘制多样性指数类型，可绘制单个
    
```

### 1.2 稀释曲线

```
    # PacBio data
    Rscript ${db}/script/alpha_rare_curve.R \
      --input result2/alpha/alpha_rare.txt --design result2/metadata.txt \
      --group Group --output result2/alpha/ \
      --width 120 --height 78
      
    # Illumina data
    Rscript ${db}/script/alpha_rare_curve.R \
      --input result_illumina/alpha/alpha_rare.txt --design result_illumina/metadata.txt \
      --group Group --output result_illumina/alpha/ \
      --width 120 --height 78
      
    # Nanopore data
    Rscript ${db}/script/alpha_rare_curve.R \
      --input result_nanopore/alpha/alpha_rare.txt --design result_nanopore/metadata.txt \
      --group Group --output result_nanopore/alpha/ \
      --width 120 --height 78
      
```

### 1.3 多样性维恩图

```
    # PacBio data
    # 支持2-4组比较
    Rscript ${db}/script/venn.R \
      --input result2/alpha/otu_group_exist.txt \
      --groups All,feces,plaque,saliva \
      --output result2/alpha/venn.pdf
    # 在线绘制韦恩图推荐在线网站http://www.ehbio.com/test/venn
    
    # Illumina data
    Rscript ${db}/script/venn.R \
      --input result_illumina/alpha/otu_group_exist.txt \
      --groups All,feces,plaque,saliva \
      --output result_illumina/alpha/venn.pdf
    
    # Nanopore data
    Rscript ${db}/script/venn.R \
      --input result_nanopore/alpha/otu_group_exist.txt \
      --groups East,North,South,Outside \
      --output result_nanopore/alpha/venn.pdf
    
```

## 2. Beta多样性

### 2.1 距离矩阵热图pheatmap
    
```
    # PacBio data
    # 以bray_curtis为例，-f输入文件,-h是否聚类TRUE/FALSE,-u/v为宽高英寸
    bash ${db}/script/sp_pheatmap.sh \
      -f result2/beta/bray_curtis.txt \
      -H 'TRUE' -u 6 -v 5
    # 添加分组注释，如2，4列的基因型和地点
    cut -f 1-2 result2/metadata.txt > temp2/group.txt
    # -P添加行注释文件，-Q添加列注释
    bash ${db}/script/sp_pheatmap.sh \
      -f result2/beta/bray_curtis.txt \
      -H 'TRUE' -u 6.9 -v 5.6 \
      -P temp2/group.txt -Q temp2/group.txt
    # 距离矩阵与相关类似，可尝试corrplot或ggcorrplot绘制更多样式
    # - [绘图相关系数矩阵corrplot](http://mp.weixin.qq.com/s/H4_2_vb2w_njxPziDzV4HQ)
    # - [相关矩阵可视化ggcorrplot](http://mp.weixin.qq.com/s/AEfPqWO3S0mRnDZ_Ws9fnw)

    # Illumina data
    # 以bray_curtis为例，-f输入文件,-h是否聚类TRUE/FALSE,-u/v为宽高英寸
    bash ${db}/script/sp_pheatmap.sh \
      -f result_illumina/beta/bray_curtis.txt \
      -H 'TRUE' -u 6 -v 5
    # 添加分组注释，如2，4列的基因型和地点
    cut -f 1-2 result_illumina/metadata.txt > temp_illumina/group.txt
    # -P添加行注释文件，-Q添加列注释
    bash ${db}/script/sp_pheatmap.sh \
      -f result_illumina/beta/bray_curtis.txt \
      -H 'TRUE' -u 6.9 -v 5.6 \
      -P temp_illumina/group.txt -Q temp_illumina/group.txt
    
    # Nanopore data
    # 以bray_curtis为例，-f输入文件,-h是否聚类TRUE/FALSE,-u/v为宽高英寸
    bash ${db}/script/sp_pheatmap.sh \
      -f result_nanopore/beta/bray_curtis.txt \
      -H 'TRUE' -u 6 -v 5
    # 添加分组注释，如2，4列的基因型和地点
    cut -f 1-2 result_nanopore/metadata.txt > temp_nanopore/group.txt
    # -P添加行注释文件，-Q添加列注释
    bash ${db}/script/sp_pheatmap.sh \
      -f result_nanopore/beta/bray_curtis.txt \
      -H 'TRUE' -u 6.9 -v 5.6 \
      -P temp_nanopore/group.txt -Q temp_nanopore/group.txt
 
```


### 2.2 主坐标分析PCoA

```
    # PacBio data
    # 输入文件，选择分组，输出文件
    Rscript ${db}/script/beta_PCoA.R \
      --input result2/otutab.txt \
      --metadata result2/metadata.txt \
      --group Group \
      --output result2/beta/PCoa2.pdf
      
    # Illumina data
    # 输入文件，选择分组，输出文件
    Rscript ${db}/script/beta_PCoA.R \
      --input result_illumina/otutab.txt \
      --metadata result_illumina/metadata.txt \
      --group Group \
      --output result_illumina/beta/PCoa2.pdf
      
    # Nanopore data
    # 输入文件，选择分组，输出文件
    Rscript ${db}/script/beta_PCoA.R \
      --input result_nanopore/otutab.txt \
      --metadata result_nanopore/metadata.txt \
      --group Group \
      --output result_nanopore/beta/PCoa2.pdf
      
```

### 2.3 限制性主坐标分析CPCoA

```
    # PacBio data
    Rscript ${db}/script/beta_cpcoa.R \
      --input result2/beta/bray_curtis.txt --design result2/metadata.txt \
      --group Group --output result2/beta/bray_curtis.cpcoa.pdf \
      --width 89 --height 59
    # 添加样本标签 --label TRUE
    Rscript ${db}/script/beta_cpcoa.R \
      --input result2/beta/bray_curtis.txt --design result2/metadata.txt \
      --group Group --label TRUE --width 89 --height 59 \
      --output result2/beta/bray_curtis.cpcoa.label.pdf
      
    # Illumina data
    Rscript ${db}/script/beta_cpcoa.R \
      --input result_illumina/beta/bray_curtis.txt --design result_illumina/metadata.txt \
      --group Group --output result_illumina/beta/bray_curtis.cpcoa.pdf \
      --width 89 --height 59
    # 添加样本标签 --label TRUE
    Rscript ${db}/script/beta_cpcoa.R \
      --input result_illumina/beta/bray_curtis.txt --design result_illumina/metadata.txt \
      --group Group --label TRUE --width 89 --height 59 \
      --output result_illumina/beta/bray_curtis.cpcoa.label.pdf
      
    # Nanopore data
    Rscript ${db}/script/beta_cpcoa.R \
      --input result_nanopore/beta/bray_curtis.txt --design result_nanopore/metadata.txt \
      --group Group --output result_nanopore/beta/bray_curtis.cpcoa.pdf \
      --width 89 --height 59
    # 添加样本标签 --label TRUE
    Rscript ${db}/script/beta_cpcoa.R \
      --input result_nanopore/beta/bray_curtis.txt --design result_nanopore/metadata.txt \
      --group Group --label TRUE --width 89 --height 59 \
      --output result_nanopore/beta/bray_curtis.cpcoa.label.pdf
      
```

## 3. 物种组成Taxonomy

### 3.1 堆叠柱状图Stackplot

```
    # PacBio data
    # 以门(p)水平为例，结果包括output.sample/group.pdf两个文件
    Rscript ${db}/script/tax_stackplot.R \
      --input result2/tax/sum_p.txt --design result2/metadata.txt \
      --group Group -t 10 --color manual1 --legend 7 --width 89 --height 59 \
      --output result2/tax/sum_p.stackplot
    # 补充-t 筛选丰度前多少物种及-s x轴排列位置-s "feces,plaque,saliva" 
    # 修改颜色--color ggplot, manual1(30), Paired(12) or Set3(12)

    # 批量绘制输入包括p/c/o/f/g共5级
    for i in p c o f g s; do
    Rscript ${db}/script/tax_stackplot.R \
      --input result2/tax/sum_${i}.txt --design result2/metadata.txt \
      --group Group -t 10 --output result2/tax/sum_${i}.stackplot \
      --legend 8 --width 89 --height 59; done
      
    # 组间比较堆叠柱状图
    Rscript ${db}/script/tax_stack_compare.R \
      --input result2/tax/data_practice53_16S.txt \
      --compare Illumina-PacBio \
      --output result2/tax/
      
    # 连线堆叠柱状图
    Rscript ${db}/script/microeco_alluvial.R \
      --otu result2/otutab2.txt \
      --metadata result2/metadata.txt \
      --taxonomy result2/taxonomy.txt \
      --output result2/tax/
      
    # Illumina data
    # 以门(p)水平为例，结果包括output.sample/group.pdf两个文件
    Rscript ${db}/script/tax_stackplot.R \
      --input result_illumina/tax/sum_p.txt --design result_illumina/metadata.txt \
      --group Group -t 10 --color manual1 --legend 7 --width 89 --height 59 \
      --output result_illumina/tax/sum_p.stackplot
    # 补充-t 筛选丰度前多少物种及-s x轴排列位置-s "feces,plaque,saliva" 
    # 修改颜色--color ggplot, manual1(30), Paired(12) or Set3(12)

    # 批量绘制输入包括p/c/o/f/g共5级
    for i in p c o f g s; do
    Rscript ${db}/script/tax_stackplot.R \
      --input result_illumina/tax/sum_${i}.txt --design result_illumina/metadata.txt \
      --group Group -t 10 --output result_illumina/tax/sum_${i}.stackplot \
      --legend 8 --width 89 --height 59; done
      
    # 组间比较堆叠柱状图
    Rscript ${db}/script/tax_stack_compare.R \
      --input result_illumina/tax/data_illumina_pacbio2.txt \
      --compare Illumina-PacBio \
      --output result_illumina/tax/
      
    # 连线堆叠柱状图
    Rscript ${db}/script/microeco_alluvial.R \
      --otu result_illumina/otutab2.txt \
      --metadata result_illumina/metadata.txt \
      --taxonomy result_illumina/taxonomy.txt \
      --output result_illumina/tax/
      
    # Nanopoare data
    # 以门(p)水平为例，结果包括output.sample/group.pdf两个文件
    Rscript ${db}/script/tax_stackplot.R \
      --input result_nanopore/tax/sum_p.txt --design result_nanopore/metadata.txt \
      --group Group -t 10 --color manual1 --legend 7 --width 89 --height 59 \
      --output result_nanopore/tax/sum_p.stackplot
    # 补充-t 筛选丰度前多少物种及-s x轴排列位置-s "feces,plaque,saliva" 
    # 修改颜色--color ggplot, manual1(30), Paired(12) or Set3(12)

    # 批量绘制输入包括p/c/o/f/g共5级
    for i in p c o f g s; do
    Rscript ${db}/script/tax_stackplot.R \
      --input result_nanopore/tax/sum_${i}.txt --design result_nanopore/metadata.txt \
      --group Group -t 10 --output result_nanopore/tax/sum_${i}.stackplot \
      --legend 8 --width 89 --height 59; done
      
    # 连线堆叠柱状图
    Rscript ${db}/script/microeco_alluvial.R \
      --otu result_nanopore/otutab2.txt \
      --metadata result_nanopore/metadata.txt \
      --taxonomy result_nanopore/taxonomy.txt \
      --output result_nanopore/tax/
      
```

### 3.2 弦/圈图circlize

```
    # PacBio data
    # 以纲(class,c)为例，绘制前5组
    i=c
    Rscript ${db}/script/tax_circlize.R \
      --input result2/tax/sum_${i}.txt --design result2/metadata.txt \
      --group Group --legend 5
    # 结果位于当前目录circlize.pdf(随机颜色)，circlize_legend.pdf(指定颜色+图例)
    # 移动并改名与分类级一致
    mv circlize.pdf result2/tax/sum_${i}.circlize.pdf
    mv circlize_legend.pdf result2/tax/sum_${i}.circlize_legend.pdf
    
    # Illumina data
    # 以纲(class,c)为例，绘制前5组
    i=c
    Rscript ${db}/script/tax_circlize.R \
      --input result_illumina/tax/sum_${i}.txt --design result_illumina/metadata.txt \
      --group Group --legend 5
    # 结果位于当前目录circlize.pdf(随机颜色)，circlize_legend.pdf(指定颜色+图例)
    # 移动并改名与分类级一致
    mv circlize.pdf result_illumina/tax/sum_${i}.circlize.pdf
    mv circlize_legend.pdf result_illumina/tax/sum_${i}.circlize_legend.pdf
    
    # Nanopore data
    # 以纲(class,c)为例，绘制前5组
    i=c
    Rscript ${db}/script/tax_circlize.R \
      --input result_nanopore/tax/sum_${i}.txt --design result_nanopore/metadata.txt \
      --group Group --legend 5
    # 结果位于当前目录circlize.pdf(随机颜色)，circlize_legend.pdf(指定颜色+图例)
    # 移动并改名与分类级一致
    mv circlize.pdf result_nanopore/tax/sum_${i}.circlize.pdf
    mv circlize_legend.pdf result_nanopore/tax/sum_${i}.circlize_legend.pdf
    
    
```

### 3.3 气泡图

```
    # PacBio data
    # 以属为例（genus，g），绘制丰度前15的属水平丰度气泡图；输入物种丰度表和样本metadata文件，输出分组物种丰度气泡图
    i=g2
    Rscript ${db}/script/tax_bubble.R \
    -i result2/tax/sum_${i}.txt \
    -g result2/metadata.txt \
    -c Group -w 7 -e 4 -s 15 -n 15 \
    -o result2/tax/sum_g.bubble3.pdf
    
    # Illumina data
    # 以属为例（genus，g），绘制丰度前15的属水平丰度气泡图；输入物种丰度表和样本metadata文件，输出分组物种丰度气泡图
    i=g
    Rscript ${db}/script/tax_bubble.R \
    -i result_illumina/tax/sum_${i}.txt \
    -g result_illumina/metadata.txt \
    -c Group -w 7 -e 4 -s 15 -n 15 \
    -o result_illumina/tax/sum_g.bubble3.pdf
    
    # Nanopore data
    # 以属为例（genus，g），绘制丰度前15的属水平丰度气泡图；输入物种丰度表和样本metadata文件，输出分组物种丰度气泡图
    i=g
    Rscript ${db}/script/tax_bubble.R \
    -i result_nanopore/tax/sum_${i}.txt \
    -g result_nanopore/metadata.txt \
    -c Group -w 7 -e 4 -s 15 -n 15 \
    -o result_nanopore/tax/sum_g.bubble3.pdf
    
```

### 3.4 甜甜圈图(Donut)和雷达图(Ladar plot)展示相对丰度组成

```
    # PacBio data
    Rscript ${db}/script/Donut_plot.R \
      --otu_table result2/otutab2.txt \
      --metadata result2/metadata.txt \
      --taxonomy result2/taxonomy.txt \
      --output_dir result2/tax/
    
    # Illumina data
    Rscript ${db}/script/Donut_plot.R \
      --otu_table result_illumina/otutab2.txt \
      --metadata result_illumina/metadata.txt \
      --taxonomy result_illumina/taxonomy.txt \
      --output_dir result_illumina/tax/
    
    # Nanopore data
    Rscript ${db}/script/Donut_plot.R \
      --otu_table result_nanopore/otutab2.txt \
      --metadata result_nanopore/metadata.txt \
      --taxonomy result_nanopore/taxonomy.txt \
      --output_dir result_nanopore/tax/
```

### 3.5 核心微生物(至少在80%样本中存在)

```
    # PacBio data
    Rscript ${db}/script/core_ASVS_and_other.R \
      --otu_table result2/otutab2.txt \
      --metadata result2/metadata.txt \
      --taxonomy result2/taxonomy.txt \
      --output_dir result2/tax/
    
    # 核心微生物散点图
    Rscript ${db}/script/Core_Abundance_ScatterPlot.R \
      --otu_table result2/otutab2.txt \
      --metadata result2/metadata.txt \
      --output result2/tax/
      
    # Illumina data
    Rscript ${db}/script/core_ASVS_and_other.R \
      --otu_table result_illumina/otutab2.txt \
      --metadata result_illumina/metadata.txt \
      --taxonomy result_illumina/taxonomy.txt \
      --output_dir result_illumina/tax/
      
    # 核心微生物散点图
    Rscript ${db}/script/Core_Abundance_ScatterPlot.R \
      --otu_table result_illumina/otutab2.txt \
      --metadata result_illumina/metadata.txt \
      --output result_illumina/tax/

    # Nanopore data
    Rscript ${db}/script/core_ASVS_and_other.R \
      --otu_table result_nanopore/otutab2.txt \
      --metadata result_nanopore/metadata.txt \
      --taxonomy result_nanopore/taxonomy.txt \
      --output_dir result_nanopore/tax/
      
    # 核心微生物散点图
    Rscript ${db}/script/Core_Abundance_ScatterPlot.R \
      --otu_table result_nanopore/otutab2.txt \
      --metadata result_nanopore/metadata.txt \
      --output result_nanopore/tax/
```


# 4、差异比较

## 1. R语言差异分析

### 1.1 差异比较

```
    # PacBio data
    mkdir -p result2/compare/
    # 输入特征表、元数据；指定分组列名、比较组和丰度
    # 选择方法 wilcox/t.test/edgeR、pvalue和fdr和输出目录
    # 这里选择默认的wilcox
    compare="saliva-plaque"
    Rscript ${db}/script/compare.R \
      --input result2/otutab.txt --design result2/metadata.txt \
      --group Group --compare ${compare} --threshold 0.01 \
      --pvalue 0.05 --fdr 0.2 \
      --output result2/compare/
    # 并筛选丰度为前20%的ASV（--threshold 0.05 ）用来画热图
    compare="feces-saliva"
    Rscript ${db}/script/compare.R \
      --input result2/otutab.txt --design result2/metadata.txt \
      --group Group --compare ${compare} --threshold 0.2 \
      --pvalue 0.05 --fdr 0.2 \
      --output result2/compare/
    # 常见错误：Error in file(file, ifelse(append, "a", "w")) : 无法打开链结 Calls: write.table -> file
    # 解决方法：输出目录不存在，创建目录即可
    
    # Illumina data
    mkdir -p result_illumina/compare/
    # 输入特征表、元数据；指定分组列名、比较组和丰度
    # 选择方法 wilcox/t.test/edgeR、pvalue和fdr和输出目录
    # 这里选择默认的wilcox
    compare="saliva-plaque"
    Rscript ${db}/script/compare.R \
      --input result_illumina/otutab.txt --design result_illumina/metadata.txt \
      --group Group --compare ${compare} --threshold 0.01 \
      --pvalue 0.05 --fdr 0.2 \
      --output result_illumina/compare/
    # 并筛选丰度为前20%的ASV（--threshold 0.05 ）用来画热图
    compare="feces-saliva"
    Rscript ${db}/script/compare.R \
      --input result_illumina/otutab.txt --design result_illumina/metadata.txt \
      --group Group --compare ${compare} --threshold 0.2 \
      --pvalue 0.05 --fdr 0.2 \
      --output result_illumina/compare/
    # 常见错误：Error in file(file, ifelse(append, "a", "w")) : 无法打开链结 Calls: write.table -> file
    # 解决方法：输出目录不存在，创建目录即可
    
    # Nanopore data
    mkdir -p result_nanopore/compare/
    # 输入特征表、元数据；指定分组列名、比较组和丰度
    # 选择方法 wilcox/t.test/edgeR、pvalue和fdr和输出目录
    # 这里选择默认的wilcox
    compare="East-Outside"
    Rscript ${db}/script/compare.R \
      --input result_nanopore/otutab.txt --design result_nanopore/metadata.txt \
      --group Group --compare ${compare} --threshold 0.01 \
      --pvalue 0.05 --fdr 0.2 \
      --output result_nanopore/compare/
    # 并筛选丰度为前20%的ASV（--threshold 0.05 ）用来画热图
    compare="East-Outside"
    Rscript ${db}/script/compare.R \
      --input result_nanopore/otutab.txt --design result_nanopore/metadata.txt \
      --group Group --compare ${compare} --threshold 0.2 \
      --pvalue 0.05 --fdr 0.2 \
      --output result_nanopore/compare/
    # 常见错误：Error in file(file, ifelse(append, "a", "w")) : 无法打开链结 Calls: write.table -> file
    # 解决方法：输出目录不存在，创建目录即可
    

```

### 1.2 火山图
    
``` 
    # PacBio data
    Rscript ${db}/script/volcano2.R \
      --input result2/compare/saliva_plaque2.txt \
      --group saliva-plaque \
      --output result2/compare/
      
    # 多组比较火山图
    Rscript ${db}/script/multigroup_compare_volcano.R \
      --input result2/otutab2.txt \
      --metadata result2/metadata.txt \
      --output result2/compare/
    
    # Illumina data
    Rscript ${db}/script/volcano2.R \
      --input result_illumina/compare/saliva_plaque2.txt \
      --group saliva-plaque \
      --output result_illumina/compare/
      
    # 多组比较火山图
    Rscript ${db}/script/multigroup_compare_volcano.R \
      --input result_illumina/otutab2.txt \
      --metadata result_illumina/metadata.txt \
      --output result_illumina/compare/
      
    # Nanopore data
    Rscript ${db}/script/volcano2.R \
      --input result_nanopore/compare/East-Outside2.txt \
      --group East-Outside \
      --output result_nanopore/compare/
      
    # 多组比较火山图
    Rscript ${db}/script/multigroup_compare_volcano.R \
      --input result_nanopore/otutab2.txt \
      --metadata result_nanopore/metadata.txt \
      --output result_nanopore/compare/

```

### 1.3 热图

```
    # PacBio data
    # 输入compare.R的结果，筛选列数，指定元数据和分组、物种注释，图大小英寸和字号
    compare="saliva_plaque2"
    bash ${db}/script/compare_heatmap.sh \
       -i result2/compare/${compare}.txt -l 7 \
       -d result2/metadata.txt -A Group \
       -t result2/taxonomy.txt \
       -w 12 -h 20 -s 14 \
       -o result2/compare/${compare}
      
    # 多样本分组比较热图
    Rscript ${db}/script/multisample_compare_heatmap.R \
      --input result2/compare/data5.txt \
      --output result2/compare/
    
    # Illumina data
    # 输入compare.R的结果，筛选列数，指定元数据和分组、物种注释，图大小英寸和字号
    compare="saliva_plaque2"
    bash ${db}/script/compare_heatmap.sh \
       -i result_illumina/compare/${compare}.txt -l 7 \
       -d result_illumina/metadata.txt -A Group \
       -t result_illumina/taxonomy.txt \
       -w 12 -h 20 -s 14 \
       -o result_illumina/compare/${compare}
      
    # 多样本分组比较热图
    Rscript ${db}/script/multisample_compare_heatmap.R \
      --input result_illumina/compare/data5.txt \
      --output result_illumina/compare/
      
    # Nanopore data
    # 输入compare.R的结果，筛选列数，指定元数据和分组、物种注释，图大小英寸和字号
    compare="East-Outside2"
    bash ${db}/script/compare_heatmap.sh \
       -i result_nanopore/compare/${compare}.txt -l 7 \
       -d result_nanopore/metadata.txt -A Group \
       -t result_nanopore/taxonomy.txt \
       -w 12 -h 20 -s 14 \
       -o result_nanopore/compare/${compare}
      
    # 多样本分组比较热图
    Rscript ${db}/script/multisample_compare_heatmap.R \
      --input result_nanopore/compare/data5.txt \
      --output result_nanopore/compare/
    
```

### 1.4 三元图

```
    # PacBio data
    Rscript ${db}/script/Ternary_plot.R   \
    --input result2/tax/sum_p.txt   \
    --metadata result2/metadata.txt  \
    --group Group   \
    --taxlevel Phylum \
    --output result2/compare/ternary_p.pdf   \
    --topn 10
    
    # Illumina data
    Rscript ${db}/script/Ternary_plot.R   \
    --input result_illumina/tax/sum_p.txt   \
    --metadata result_illumina/metadata.txt  \
    --group Group   \
    --taxlevel Phylum \
    --output result_illumina/compare/ternary_p.pdf   \
    --topn 10
    
    # Nanopore data
    Rscript ${db}/script/Ternary_plot.R   \
    --input result_nanopore/tax/sum_p.txt   \
    --metadata result_nanopore/metadata.txt  \
    --group Group   \
    --taxlevel Phylum \
    --output result_nanopore/compare/ternary_p.pdf   \
    --topn 10

```



## 1.5. STAMP差异分析图

```
### 1.5.1 生成输入文件(备选)

    # PacBio data
    Rscript ${db}/script/format2stamp.R -h
    mkdir -p result2/stamp
    Rscript ${db}/script/format2stamp.R --input result2/otutab.txt \
      --taxonomy result2/taxonomy.txt --threshold 0.01 \
      --output result2/stamp/tax
    # 可选Rmd文档见result/format2stamp.Rmd
    
    # Illumina data
    Rscript ${db}/script/format2stamp.R -h
    mkdir -p result_illumina/stamp
    Rscript ${db}/script/format2stamp.R --input result_illumina/otutab.txt \
      --taxonomy result_illumina/taxonomy.txt --threshold 0.01 \
      --output result_illumina/stamp/tax
    # 可选Rmd文档见result/format2stamp.Rmd
    
    # Nanopore data
    Rscript ${db}/script/format2stamp.R -h
    mkdir -p result_nanopore/stamp
    Rscript ${db}/script/format2stamp.R --input result_nanopore/otutab.txt \
      --taxonomy result_nanopore/taxonomy.txt --threshold 0.01 \
      --output result_nanopore/stamp/tax
    # 可选Rmd文档见result/format2stamp.Rmd
    
    
### 1.5.2 R语言绘制STAMP图

    # 此处对统计还有问题，应该用wilcox.test，不能用t.test，还需要后续调整代码

    # PacBio data
    Rscript ${db}/script/stamp2.R \
      --input result2/tax/sum_g2.txt \
      --group result2/metadata.txt \
      --compare feces-plaque \
      --output result2/tax/
      
    # Illumina data
    Rscript ${db}/script/stamp2.R \
      --input result_illumina/tax/sum_g.txt \
      --group result_illumina/metadata.txt \
      --compare feces-plaque \
      --output result_illumina/tax/
    
    # Illumina data
    Rscript ${db}/script/stamp2.R \
      --input result_nanopore/tax/sum_g.txt \
      --group result_nanopore/metadata.txt \
      --compare East-Outside \
      --output result_nanopore/tax/

```

### Lefse分析及组间比较柱状图

    Lefse分析在线网站：https://www.bioincloud.tech/standalone-task-ui/lefse
    Reference: Gao, Yunyun, Guoxing Zhang, Shunyao Jiang, and Yong‐Xin Liu. 2024. “
    Wekemo Bioincloud: A User‐friendly Platform for Meta‐omics Data Analyses.” i
    Meta e175. https://doi.org/10.1002/imt2.17
 
 
### 2.环境因子关联分析
 
```
### Mantel检验相关性分析热图
 
    # Nanopore data
    Rscript ${db}/script/env_mantel_heatmap.R \
      --input result2/tax/otutab2.txt \
      --env result2/tax/env_amplicon.txt \
      --output result2/tax/

```

```
### Redundancy analysis (RDA) 冗余分析

    # Nanopore data
    Rscript ${db}/script/RDA_microeco.R \
      --input result2/tax/otutab.txt \
      --metadata result2/tax/metadata.txt \
      --tax result2/tax/taxonomy.txt \
      --phylo result2/tax/otus.tree \
      --output result2/tax/

``` 

### 3.系统发育树

```
### 3.1 无进化距离系统发育树

    # PacBio data
    Rscript ${db}/script/phylogenetic_tree01.R \
      --input result2/tax/otus.nwk \
      --anno result2/tax/annotation2.txt \
      --output result2/tax/
      
### 3.2 有进化距离系统发育树
    # PacBio data
    Rscript ${db}/script/phylogenetic_tree02.R \
      --input result2/tax/otus.nwk \
      --anno result2/tax/annotation3.txt \
      --output result2/tax/
    
```


### 4.网络分析

```
### 4.1 多组Spearman相关性网络比较分析 (Multi-group Spearman correlation network)
    
    # PacBio data
    Rscript ${db}/script/Spearman_network01.R \
      --input result2/tax/otutab_amplicon.txt \
      --group result2/tax/sample_amplicon.txt \
      --tax result2/tax/taxonomy_amplicon.txt \
      --output result2/tax/

```

### 5.随机森林模型

```
### 5.1 运行随机森林模型

    # PacBio data
    # 注意：这是第一版本的随机森林模型，其中的一些参数之间内置在代码中，使用时需要根据具体样本情况进行调整
    Rscript ${db}/script/random_forest01.R \
      --input result2/tax/PacBio_data.txt \
      --group result2/tax/PacBio_metadata.txt \
      --output result2/tax/RF_model/
      
### 5.2 绘制柱状图
    
    # PacBio data
    Rscript ${db}/script/RF_plot01.R \
      --input result2/tax/RF_model/Species_imp_rf21.txt \
      --tax result2/tax/taxonomy_amplicon.txt \
      --optimal 4 \
      --output result2/tax/RF_model/
    
``` 

### 6.功能差异分析

```
### 6.1 数据处理及差异分析
    
    # PacBio data
    Rscript ${db}/script/function_data_process.R \
      --input result2/tax/KEGG.PathwayL2.raw.txt \
      --group result2/tax/metadata_amplicon.txt \
      --output result2/tax/pathway/
    
### 6.2 热图结合柱状图展示差异

    # PacBio data
    Rscript ${db}/script/function_diff_plot01.R \
      --input result2/tax/pathway/Difference_pathway21.txt \
      --pathway result2/tax/pathway/pathway_count_data.txt \
      --statistic result2/tax/pathway/pathway_percent_abundance2.txt \
      --output result2/tax/pathway/

```

### 注意事项

# 1.绘图代码中的部分参数需要根据特定研究进行修改；
# 2.可到EasyMicrobiome(https://github.com/YongxinLiu/EasyMicrobiome)文件夹下获取相应的R语言脚本，可自行对代码进行调整;
# 3.如有疑问欢迎提问交流(https://github.com/YongxinLiu/EasyAmplicon/pulls)。



