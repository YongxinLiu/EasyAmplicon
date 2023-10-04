[TOC]

# QIIME2 2023.7分析流程

## 0. 软件安装(附录1)

    # 仅限在Linux/Mac系统中运行，
    # Windows用户使用内置Linux子系统(支持右键粘贴)、Linux服务器环境下安装并运行
    # 详细教程参阅官网https://docs.qiime2.org/2023.7/，或QIIME2教程专辑 https://mp.weixin.qq.com/s/farGisfX3fVL_5WgXS8lIg
    # 安装Windows子系统，https://mp.weixin.qq.com/s/0PfA0bqdvrEbo62zPVq4kQ

## 1. 准备工作

    # 设置工作目录，如服务器为~/amplicon/qiime2，Win子系统如下：
    wd=/mnt/c/amplicon/qiime2/
    # 进入工作目录
    mkdir -p ${wd}
    cd ${wd}
    # 激活QIIME2工作环境，旧版conda使用source替换conda运行
    conda activate qiime2-2023.7
    
    # 准备样本元数据metadata.txt、原始数据seq/*.fq.gz
    
    ## 此处代码基于metadata.txt从公共数据下载测序数据，按GSA的CRA(批次)和CRR(样品)编号下载数据
    mkdir -p seq
    # # 公网下载
    # cd seq
    # awk '{system("wget -c ftp://download.big.ac.cn/gsa/"$5"/"$6"/"$6"_f1.fq.gz -O "$1"_1.fq.gz")}' \
    #     <(tail -n+2 ../metadata.txt)
    #  awk '{system("wget -c ftp://download.big.ac.cn/gsa/"$5"/"$6"/"$6"_r2.fq.gz -O "$1"_2.fq.gz")}' \
    #     <(tail -n+2 ../metadata.txt)
    # cd .. && ls -lsh seq
    # 从其他地方链接(不额外占用空间)
    ln /mnt/c/amplicon/seq/* seq/
    ln /mnt/c/amplicon/result/metadata.txt ./
    # 根据metadata生成manifest文件
    awk 'NR==1{print "sample-id\tforward-absolute-filepath\treverse-absolute-filepath"} \
      NR>1{print $1"\t$PWD/seq/"$1"_1.fq.gz\t$PWD/seq/"$1"_2.fq.gz"}' \
      metadata.txt > manifest
    head -n3 manifest
    
    # 数据导入qiime2，格式为双端33格式
    qiime tools import \
      --type 'SampleData[PairedEndSequencesWithQuality]' \
      --input-path manifest \
      --output-path demux.qza \
      --input-format PairedEndFastqManifestPhred33V2
    # fq文件1G用时7m，fq.gz压缩格式仅34s，测试数据9s


## 2. 生成特征表和代表序列

### 方法1. DADA2(慢，依赖R和Python包多容易报错)

    # 支持多线程加速，90万条PE250数据，0/96p, 34m；24p, 44m；8p, 77m；1p, 462m
    # 27万，8p, 9m; 4p, 11m;
    time qiime dada2 denoise-paired \
      --i-demultiplexed-seqs demux.qza \
      --p-n-threads 4 \
      --p-trim-left-f 29 --p-trim-left-r 18 \
      --p-trunc-len-f 0 --p-trunc-len-r 0 \
      --o-table dada2-table.qza \
      --o-representative-sequences dada2-rep-seqs.qza \
      --o-denoising-stats denoising-stats.qza
    # 确定使用dada2结果并导入主流程
    cp dada2-table.qza table.qza
    cp dada2-rep-seqs.qza rep-seqs.qza

### 方法2. 外部导入特征表和代表序列(常用)

    # 上传我们生成的OTU表otutab.txt和代表序列otus.fa
    # 转换文本为Biom1.0，注意biom --version 2.1.5/8可以，2.1.7报错
    biom convert -i otutab.txt -o otutab.biom \
      --table-type="OTU table" --to-json
    # 导入特征表，9s
    qiime tools import --input-path otutab.biom \
      --type 'FeatureTable[Frequency]' --input-format BIOMV100Format \
      --output-path table.qza
    # 导入代表序列，8s
    qiime tools import --input-path otus.fa \
      --type 'FeatureData[Sequence]' \
      --output-path rep-seqs.qza

### 特征表和代表序列统计

    qiime feature-table summarize \
      --i-table table.qza \
      --o-visualization table.qzv \
      --m-sample-metadata-file metadata.txt
    qiime feature-table tabulate-seqs \
      --i-data rep-seqs.qza \
      --o-visualization rep-seqs.qzv
    # 下载qzv在线查看，dada2只有1千多个ASV


## 3. Alpha和beta多样性分析

### 构建进化树用于多样性分析 53s

    qiime phylogeny align-to-tree-mafft-fasttree \
      --i-sequences rep-seqs.qza \
      --o-alignment aligned-rep-seqs.qza \
      --o-masked-alignment masked-aligned-rep-seqs.qza \
      --o-tree unrooted-tree.qza \
      --o-rooted-tree rooted-tree.qza

### 计算核心多样性 

    # 13s，采样深度通常选择最小值，来自table.qzv
    qiime diversity core-metrics-phylogenetic \
      --i-phylogeny rooted-tree.qza \
      --i-table table.qza \
      --p-sampling-depth 7439 \
      --m-metadata-file metadata.txt \
      --output-dir core-metrics-results

### Alpha多样性组间显著性分析和可视化

    # 7s, 可选的alpha指数有 faith_pd、shannon、observed_features、evenness
    index=observed_features
    qiime diversity alpha-group-significance \
      --i-alpha-diversity core-metrics-results/${index}_vector.qza \
      --m-metadata-file metadata.txt \
      --o-visualization core-metrics-results/${index}-group-significance.qzv

### Alpha多样性稀疏曲线

    # 25s, max-depth选最大值，来自table.qzv
    qiime diversity alpha-rarefaction \
      --i-table table.qza \
      --i-phylogeny rooted-tree.qza \
      --p-max-depth 10298 \
      --m-metadata-file metadata.txt \
      --o-visualization alpha-rarefaction.qzv
    # 结果有observed_otus, shannon, 和faith_pd三种指数可选

### Beta多样性组间显著性分析和可视化

    # 可选的beta指数有 unweighted_unifrac、bray_curtis、weighted_unifrac和jaccard
    # 7s, 指定分组是减少计算量，置换检验较耗时
    distance=weighted_unifrac
    column=Group
    qiime diversity beta-group-significance \
      --i-distance-matrix core-metrics-results/${distance}_distance_matrix.qza \
      --m-metadata-file metadata.txt \
      --m-metadata-column ${column} \
      --o-visualization core-metrics-results/${distance}-${column}-significance.qzv \
      --p-pairwise


## 4. 物种组成分析

    # 物种注释，数据库见附录，可选silva-138-99-nb-classifier.qza 或 2022.10.backbone.full-length.nb.qza，详见附录或官方教程
    # 此步计算时间5m，内存使用6GB
    time qiime feature-classifier classify-sklearn \
      --i-classifier 2022.10.backbone.full-length.nb.qza \
      --i-reads rep-seqs.qza \
      --o-classification taxonomy.qza
    # 可视化物种注释
    qiime metadata tabulate \
      --m-input-file taxonomy.qza \
      --o-visualization taxonomy.qzv
    # 堆叠柱状图展示
    qiime taxa barplot \
      --i-table table.qza \
      --i-taxonomy taxonomy.qza \
      --m-metadata-file metadata.txt \
      --o-visualization taxa-bar-plots.qzv


## 5. 差异分析ancom

    # 格式化特征表，添加伪计数，4s
    qiime composition add-pseudocount \
      --i-table table.qza \
      --o-composition-table comp-table.qza
    
    # 计算差异特征，指定分组类型比较，1m
    column=Group
    time qiime composition ancom \
      --i-table comp-table.qza \
      --m-metadata-file metadata.txt \
      --m-metadata-column ${column} \
      --o-visualization ancom-${column}.qzv
    
    # 按属水平合并，并统计
    ## 按属水平合并，6s
    qiime taxa collapse \
      --i-table table.qza \
      --i-taxonomy taxonomy.qza \
      --p-level 6 \
      --o-collapsed-table table-l6.qza
    # 格式化特征表，添加伪计数，6s
    qiime composition add-pseudocount \
      --i-table table-l6.qza \
      --o-composition-table comp-table-l6.qza
    # 计算差异属，指定分组类型比较，16s
    qiime composition ancom \
      --i-table comp-table-l6.qza \
      --m-metadata-file metadata.txt \
      --m-metadata-column ${column} \
      --o-visualization ancom-l6-${column}.qzv


# 附录

## 1. qiime2 2023.7安装

### 安装Conda

    # 下载、安装和启动conda
    wget -c https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
    bash Miniconda3-latest-Linux-x86_64.sh -b -f
    ~/miniconda3/condabin/conda init
    source ~/.bashrc

### 方法1. Conda在线安装QIIME

    # 附软件在线安装和打包代码
    n=qiime2-2023.7
    # 下载软件列表
    wget -c https://data.qiime2.org/distro/core/${n}-py38-linux-conda.yml
    # 备用链接
    wget -c http://www.imeta.science/db/conda/${n}-py38-linux-conda.yml
    # 新环境安装，可在不同电脑服务器上安装成功后打包分发
    conda env create -n ${n} --file ${n}-py38-linux-conda.yml
    # 环境打包(可选，1.2G)
    conda pack -n ${n} -o ${n}.tar.gz

### 方法2. 本地安装QIIME

    n=qiime2-2023.7
    # 安装包下载链接,或百度云下载链接：https://pan.baidu.com/s/1Ikd_47HHODOqC3Rcx6eJ6Q?pwd=0315
    wget -c ftp://download.nmdc.cn/tools/conda/${n}.tar.gz
    # 新环境安装
    mkdir -p ~/miniconda3/envs/${n}
    tar -xzf ${n}.tar.gz -C ~/miniconda3/envs/${n}
    # 激活并初始化环境
    conda activate ${n}
    conda unpack

## 2. 物种注释数据训练集

### Silva 138 99% OTUs full-length sequences

    # 官网下载
    wget -c https://data.qiime2.org/2023.7/common/silva-138-99-nb-classifier.qza
    # 备用链接
    wget -c ftp://download.nmdc.cn/tools/amplicon/silva/silva-138-99-nb-classifier.qza

### Greengenes2 2022.10 full length sequences

    # 官网下载
    wget -c http://ftp.microbio.me/greengenes_release/2022.10/2022.10.backbone.full-length.nb.qza
    # 备用链接
    wget -c ftp://download.nmdc.cn/tools/amplicon/GreenGenes/2022.10.backbone.full-length.nb.qza
    

## 3. 物种注释数据训练集

    wd=/mnt/c/amplicon/qiime2
    mkdir -p $wd
    cd $wd
    # 下载数据库文件(greengenes, 320M)
    # wget -c ftp://greengenes.microbio.me/greengenes_release/gg_13_5/gg_13_8_otus.tar.gz
    # 国内备用链接
    # wget -c ftp://download.nmdc.cn/tools/amplicon/GreenGenes/gg_13_8_otus.tar.gz
    # 国内备份核心99数据库(60M)
    wget -c ftp://download.nmdc.cn/tools/amplicon/GreenGenes/gg_13_8_otus_99.tar.gz
    mv gg_13_8_otus_99.tar.gz gg_13_8_otus.tar.gz
    # 解压
    tar -zxvf gg_13_8_otus.tar.gz
    
    # 使用rep_set文件中的99_otus.fasta数据和taxonomy中的99_OTU_taxonomy.txt数据作为参考物种注释
    # 导入参考序列，50s
    qiime tools import \
      --type 'FeatureData[Sequence]' \
      --input-path gg_13_8_otus/rep_set/99_otus.fasta \
      --output-path 99_otus.qza
    # Fontconfig error: Cannot load default config file 不影响结果
    
    # 导入物种分类信息，6s
    qiime tools import \
      --type 'FeatureData[Taxonomy]' \
      --input-format HeaderlessTSVTaxonomyFormat \
      --input-path gg_13_8_otus/taxonomy/99_otu_taxonomy.txt \
      --output-path ref-taxonomy.qza
    
    # Train the classifier（训练分类器）——全长，30m
    time qiime feature-classifier fit-classifier-naive-bayes \
      --i-reference-reads 99_otus.qza \
      --i-reference-taxonomy ref-taxonomy.qza \
      --o-classifier classifier_gg_13_8_99.qza

    # 引物提取参考序列的扩增区段 Extract reference reads
    # 常用Greengenes 13_8 99% OTUs from 341F CCTACGGGNGGCWGCAG/805R GACTACHVGGGTATCTAATCC region of sequences（分类器描述），提供测序的引物序列，截取对应的区域进行比对，达到分类的目的。
    # 本次使用引物799F-1193R，请根据实际替换, 8m
    time qiime feature-classifier extract-reads \
      --i-sequences 99_otus.qza \
      --p-f-primer AACMGGATTAGATACCCKG \
      --p-r-primer ACGTCATCCCCACCTTCC \
      --o-reads ref-seqs.qza
    # Train the classifier（训练分类器）
    # 基于筛选的指定区段，生成实验特异的分类器，7 min
    time qiime feature-classifier fit-classifier-naive-bayes \
      --i-reference-reads ref-seqs.qza \
      --i-reference-taxonomy ref-taxonomy.qza \
      --o-classifier classifier_gg_13_8_99_V5-V7.qza
    
    # 常见问题1：scikit-learn版本不兼容，重新fit-classifier-naive-bayes构建训练集即可
    Plugin error from feature-classifier:
      The scikit-learn version (0.21.2) used to generate this artifact does not match the current version of scikit-learn installed (0.22.1). Please retrain your classifier for your current deployment to prevent data-corruption errors.
    Debug info has been saved to /tmp/qiime2-q2cli-err-5ngzk2hm.log
    
    # 常见问题2: dada2运行报错，环境配置不完全，运行conda unpack初始化
    Plugin error from dada2:
    An error was encountered while running DADA2 in R (return code 255), please inspect stdout and stderr to learn more.
    Debug info has been saved to /tmp/qiime2-q2cli-err-utwt1cmu.log
    
## 引文Citation

    使用此脚本，请引用下文：
    If used this script, please cited:
    
    Evan Bolyen, Jai Ram Rideout, Matthew R. Dillon, Nicholas A. Bokulich, Christian C. Abnet, Gabriel A. Al-Ghalith, Harriet Alexander, Eric J. Alm, Manimozhiyan Arumugam, Francesco Asnicar, Yang Bai, Jordan E. Bisanz, Kyle Bittinger, Asker Brejnrod, Colin J. Brislawn, C. Titus Brown, Benjamin J. Callahan, Andrés Mauricio Caraballo-Rodríguez, John Chase, Emily K. Cope, Ricardo Da Silva, Christian Diener, Pieter C. Dorrestein, Gavin M. Douglas, Daniel M. Durall, Claire Duvallet, Christian F. Edwardson, Madeleine Ernst, Mehrbod Estaki, Jennifer Fouquier, Julia M. Gauglitz, Sean M. Gibbons, Deanna L. Gibson, Antonio Gonzalez, Kestrel Gorlick, Jiarong Guo, Benjamin Hillmann, Susan Holmes, Hannes Holste, Curtis Huttenhower, Gavin A. Huttley, Stefan Janssen, Alan K. Jarmusch, Lingjing Jiang, Benjamin D. Kaehler, Kyo Bin Kang, Christopher R. Keefe, Paul Keim, Scott T. Kelley, Dan Knights, Irina Koester, Tomasz Kosciolek, Jorden Kreps, Morgan G. I. Langille, Joslynn Lee, Ruth Ley, **Yong-Xin Liu**, Erikka Loftfield, Catherine Lozupone, Massoud Maher, Clarisse Marotz, Bryan D. Martin, Daniel McDonald, Lauren J. McIver, Alexey V. Melnik, Jessica L. Metcalf, Sydney C. Morgan, Jamie T. Morton, Ahmad Turan Naimey, Jose A. Navas-Molina, Louis Felix Nothias, Stephanie B. Orchanian, Talima Pearson, Samuel L. Peoples, Daniel Petras, Mary Lai Preuss, Elmar Pruesse, Lasse Buur Rasmussen, Adam Rivers, Michael S. Robeson, Patrick Rosenthal, Nicola Segata, Michael Shaffer, Arron Shiffer, Rashmi Sinha, Se Jin Song, John R. Spear, Austin D. Swafford, Luke R. Thompson, Pedro J. Torres, Pauline Trinh, Anupriya Tripathi, Peter J. Turnbaugh, Sabah Ul-Hasan, Justin J. J. van der Hooft, Fernando Vargas, Yoshiki Vázquez-Baeza, Emily Vogtmann, Max von Hippel, William Walters, Yunhu Wan, Mingxun Wang, Jonathan Warren, Kyle C. Weber, Charles H. D. Williamson, Amy D. Willis, Zhenjiang Zech Xu, Jesse R. Zaneveld, Yilong Zhang, Qiyun Zhu, Rob Knight, J. Gregory Caporaso. 2019. Reproducible, interactive, scalable and extensible microbiome data science using QIIME 2. **Nature Biotechnology** 37: 852-857. https://doi.org/10.1038/s41587-019-0209-9
    
    Copyright 2016-2023 Yong-Xin Liu <liuyongxin@caas.cn> 测试整理
