[TOC]

# EasyAmplicon2 base Usearch、Vsearch

    # 设置工作(work directory, wd)和软件数据库(database, db)目录
    # 添加环境变量，并进入工作目录 Add environmental variables and enter work directory
    # **每次打开Rstudio必须运行下面4行 Run it**，可选替换${db}为EasyMicrobiome安装位置
    wd=/d/EasyAmplicon_paper_materials
    db=/c/EasyMicrobiome
    PATH=$PATH:${db}/win
    cd ${wd}



## 1. 起始文件 start files

    # 1. 分析流程pipeline.sh
    # 2. 样本元信息metadata.txt，保存于result目录
    # 3. 测序数据fastq文件保存于seq目录，通常以`.fq.gz`结尾，每个样品一对文件
    # 4. 进入对应工作目录，创建临时文件存储目录，分析结束可删除
    cd Pacbio/Nanopore
    mkdir -p seq result temp 

### 1.1. 元数据/实验设计 metadata

    # 准备样本元数据result/metadata.txt
    # csvtk统计表行(样本数，不含表头)列数，-t设置列分隔为制表符，默认为;
    csvtk -t stat result/metadata_raw.txt
    # 元数据至少3列，首列为样本ID(SampleID)，结尾列为描述(Description)
    # cat查看文件，-A显示符号，"|"为管道符实现命令连用，head显示文件头，-n3控制范围前3行
    cat -A result/metadata_raw.txt | head -n3
    # windows用户结尾有^M，运行sed命令去除，再用cat -A检查
    sed 's/\r//' result/metadata_raw.txt > result/metadata.txt
    sed 's/\r//' result/metadata.txt
    cat -A result/metadata.txt | head -n3

### 1.2. 测序数据 sequencing data

    # 公司返回的三代（Pacbio或者Nanopore）全长扩增子测序结果，通常为一个样品只有一个fq/fastq.gz格式压缩文件
    # 文件名与样品名务必对应：不一致时手工修改，批量改名见"常见问题6"
    # 如果测序数据是.gz的压缩文件，有时需要使用gunzip解压后使用，vsearch通常可直接读取压缩文件
    # gunzip seq/*.gz
    # zless按页查看压缩文件，空格翻页、q退出；head默认查看前10行，-n指定行
    ls -sh seq/
    head seq/feces1.fastq
    # 每行太长，指定查看每行的1-60个字符
    head seq/feces1.fastq | cut -c 1-60
    # 统计测序数据，依赖seqkit程序
    seqkit stat seq/feces1.fastq
    # 批量统计测序数据并汇总表
    seqkit stat seq/*.fastq > result/seqkit.txt
    head result/seqkit.txt

### 1.3. 流程和数据库 pipeline & database

    # 数据库第一次使用必须解压，以后可跳过此段

    # usearchs可用16S/18S/ITS数据库：
    SILVA和UNITE，本地文件位置 ${db}/usearch/
    # 解压ITS UNITE数据库，需自行从官网或网盘db/amplicon/usearch中下载
    # gunzip -c ${db}/usearch/utax_reference_dataset_all_29.11.2022.fasta.gz >${db}/usearch/unite.fa
    seqkit stat ${db}/usearch/unite.fa # 32.6万
    # Greengene数据库用于功能注释: ftp://greengenes.microbio.me/greengenes_release/gg_13_5/gg_13_8_otus.tar.gz
    # 默认解压会删除原文件，-c指定输出至屏幕，> 写入新文件(可改名)
    gunzip -c ${db}/gg/97_otus.fasta.gz > ${db}/gg/97_otus.fa
    seqkit stat ${db}/gg/97_otus.fa


## 2. 序列合并和重命名 reads merge and rename

      
### 2.1 文件改名 Reads rename

    # # 单个序列改名示例
    # i=feces1
    # gunzip -c seq/${i}_1.fq.gz > seq/${i}.fq
    # usearch -fastx_relabel seq/${i}.fq -fastqout temp/${i}.merged.fq -prefix ${i}.
    # 
    # # 批量改名，需要有单端fastq文件，且解压(usearch不支持压缩格式)
    gunzip seq/*.gz
    time for i in `tail -n+2 result/metadata.txt|cut -f1`;do
     usearch -fastx_relabel seq/${i}.fastq -fastqout temp/${i}.merged.fastq -prefix ${i}.
     done &
    # # vsearch大数据方法参考“常见问题2”

### 2.2 改名后序列整合 integrate renamed reads

    #合并所有样品至同一文件
    cat temp/*.merged.fastq > temp/all.fq
    #查看文件大小223M，软件不同版本结果略有差异
    ls -lsh temp/all.fq
    # 查看序列名，“.”之前是否为样本名，样本名绝不允许有点 (".")
    # 样本名有点 (.) 的一个显著特征是生成的特征表会很大，特征表里面列很多，导致后面分析出现内存不足。
    # 后面分析获得特征表后要看一眼有没有问题，遇到内存不足问题，也要回头来排查。
    head -n 6 temp/all.fq|cut -c1-60


## 3. 切除引物与质控 Cut primers and quality filter

    # 用cutadapt过滤序列，去除接头："前向引物...反向引物",12s
    cutadapt -g "AGRGTTYGATYMTGGCTCAG...RGYTACCTTGTTACGACTT" \
    --error-rate=0.1 \
    -j 10 \
    --discard-untrimmed \
    -o temp/allfilter.fastq temp/all.fq

    # 细菌16S片段长度通常约为1500bp，为避免过长或过短序列干扰，使用vsearch进行片段长度筛选。筛选长度在1200~1800 bp之间的序列，输出为fasta格式
    vsearch --fastx_filter temp/allfilter.fastq --fastq_minlen 1200 --fastq_maxlen 1800 --fastaout temp/filtered.fa --fastq_qmax 93

    # 查看文件了解fa文件格式
    head temp/filtered.fa


## 4. 去冗余挑选OTU/ASV Dereplicate and cluster/denoise

### 4.1 序列去冗余 Dereplicate

    #去冗余：完全一致的序列合并，统计size信息，输出唯一序列，-minsize表示保留序列的最小条数,-sizeout输出丰度,--fasta_width 0输出FASTA时不换行（单行序列）减少文件体积 --relabel必须加序列前缀更规范, 1s
    time vsearch --derep_fulllength temp/filtered.fa --fasta_width 0 --sizeout --relabel Uni_  --output temp/uniques.fa --minuniquesize 2 --threads 8
     
    #查看uniques.fa文件
    ls -lsh temp/uniques.fa
    # Uni_1;size=165  - 去冗余后序列的名字 Uni_1；该序列在所有样品测序数据中出现 165 次

### 4.2 聚类OTU/去噪ASV Cluster or denoise

    #有两种方法：推荐unoise3去噪获得单碱基精度ASV，备选传统的97%聚类OTU(属水平精度)
    #usearch两种特征挑选方法均自带de novo去嵌合体
    #-minsize二次过滤，控制OTU/ASV数量至1-5千，方便下游统计分析

    #方法1. 97%聚类OTU，适合大数据/ASV规律不明显/reviewer要求
    #结果耗时3s, 产生264 OTUs, 去除4chimeras
    #usearch -cluster_otus temp/uniques.fa -minsize 2 \
      -otus temp/otus.fa \
      -relabel OTU_

    #方法2. ASV去噪 Denoise: predict biological sequences and filter chimeras
    #5s, 1234 good, 0 chimeras, 序列百万条可能需要几天/几周
    usearch -unoise3 temp/uniques.fa -minsize 2 \
      -zotus temp/zotus.fa
    ##三代测序数据size数会出现普遍较小的情况，根据实际情况考虑是否调整 -minsize 
    #修改序列名：Zotu改为ASV方便识别
    sed 's/Zotu/ASV_/g' temp/zotus.fa > temp/otus.fa
    head -n 2 temp/otus.fa

    #方法3. 数据过大无法使用usearch时，备选vsearch方法见"常见问题3"

### 4.3 基于参考去嵌合 Reference-based chimera detect

    # 不推荐，容易引起假阴性，因为参考数据库无丰度信息
    # 而de novo时要求亲本丰度为嵌合体16倍以上防止假阴性
    # 因为已知序列不会被去除，数据库选择越大越合理，假阴性率最低
    mkdir -p result/raw

    # 方法1. vsearch+silva去嵌合
    # vsearch --uchime_ref temp/otus.fa \
      -db ${db}/usearch/SILVA_modified.fasta \
      --nonchimeras result/raw/otus.fa
    # 1m5s, 17 chimeras
    # Win vsearch结果添加了windows换行符^M需删除，mac不要执行此命令
    # sed -i 's/\r//g' result/raw/otus.fa
    # 方法2. 不去嵌合
    cp -f temp/otus.fa result/raw/otus.fa


## 5. 特征表构建和筛选 Feature table create and filter

    # OTU和ASV统称为特征(Feature)，它们的区别是：
    # OTU通常按97%聚类后挑选最高丰度或中心的代表性序列；
    # ASV是基于序列进行去噪(排除或校正错误序列，并挑选丰度较高的可信序列)作为代表性序列

### 5.1 生成特征表

    # 方法1. usearch生成特征表，小样本(<30)快；但大样本受限且多线程效率低，78.3%, 4核27s
    # time usearch -otutab temp/filtered.fa \
    #   -otus result/raw/otus.fa \
    #   -threads 4 \
    #   -otutabout result/raw/otutab.txt

    # 方法2. vsearch生成特征表
    # id(1)：100%相似度比对Matching unique query sequences: 6946 of 32215 (21.56%) Writing OTU table (classic) 100% real time 6m56.208s
    # id(0.97)推荐：(更高数据使用率，更快)26294 of 32215 (81.62%),3m15s
    time vsearch --usearch_global temp/filtered.fa \
      --db result/raw/otus.fa \
      --id 1 --threads 4 \
    	--otutabout result/raw/otutab.txt 
    	
    	
    time vsearch --usearch_global temp/filtered.fa \
      --db result/raw/otus.fa \
      --id 0.97 --threads 4 \
    	--otutabout result/raw/otutab.txt

    # vsearch结果windows用户删除换行符^M校正为标准Linux格式
    sed -i 's/\r//' result/raw/otutab.txt
    head -n6 result/raw/otutab.txt | cut -f 1-6 |cat -A
    # csvtk统计表行列
    # 这里一定看好列数，是不是等于你的样品数；如果不等，一般是样品命名存在问题，具体看上面解释
    csvtk -t stat result/raw/otutab.txt
    
### 5.2 物种注释，且/或去除质体和非细菌 Remove plastid and non-Bacteria

    # 物种注释-去除质体和非细菌/古菌并统计比例(可选)
    # SILVA数据库(SILVA_138.2_SSURef_NR99_tax_silva.fasta)更好注释真核、质体序列
    # 修改SILVA数据库(SILVA_138.2_SSURef_NR99_tax_silva.fasta)格式以适应注释代码，输出数据库文件SILVA_modified.fasta，包含种和株系的信息
    python database_silva.py
    # 置信阈值通常0.6/0.8，vserch最低0.1/usearch可选0输出最相似物种注释用于观察潜在分类,1m58s
    vsearch --sintax result/raw/otus.fa \
      --db ${db}/usearch/SILVA_modified.fasta \
      --sintax_cutoff 0.1 \
      --tabbedout result/raw/otus.sintax 
    head result/raw/otus.sintax | cat -A
    sed -i 's/\r//' result/raw/otus.sintax


    # 方法1. 原始特征表行数
    wc -l result/raw/otutab.txt
    #R脚本选择细菌古菌(真核)、去除叶绿体、线粒体并统计比例；输出筛选并排序的OTU表
    #输入为OTU表result/raw/otutab.txt和物种注释result/raw/otus.sintax
    #输出筛选并排序的特征表result/otutab.txt和
    #统计污染比例文件result/raw/otutab_nonBac.txt和过滤细节otus.sintax.discard
    #真菌ITS数据，请改用otutab_filter_nonFungi.R脚本，只筛选真菌
    # Rscript ${db}/script/otutab_filter_nonBac.R -h # 显示参数说明
    Rscript ${db}/script/otutab_filter_nonBac.R \
      --input result/raw/otutab.txt \
      --taxonomy result/raw/otus.sintax \
      --output result/otutab.txt\
      --stat result/raw/otutab_nonBac.stat \
      --discard result/raw/otus.sintax.discard
    # 筛选后特征表行数
    wc -l result/otutab.txt
    #过滤特征表对应序列
    cut -f 1 result/otutab.txt | tail -n+2 > result/otutab.id
    usearch -fastx_getseqs result/raw/otus.fa \
        -labels result/otutab.id -fastaout result/otus.fa
    #过滤特征表对应序列注释
    awk 'NR==FNR{a[$1]=$0}NR>FNR{print a[$1]}'\
        result/raw/otus.sintax result/otutab.id \
        > result/otus.sintax

    # 方法2. 觉得筛选不合理可以不筛选
    # cp result/raw/otu* result/

    #可选统计方法：OTU表简单统计 Summary OTUs table
    usearch -otutab_stats result/otutab.txt \
      -output result/otutab.stat
    cat result/otutab.stat
    #注意最小值、分位数，或查看result/raw/otutab_nonBac.stat中样本详细数据量，用于重采样

### 5.3 等量抽样标准化

    # Normlize by subsample

    #使用vegan包进行等量重抽样，输入reads count格式Feature表result/otutab.txt
    #可指定输入文件、抽样量和随机数，输出抽平表result/otutab_rare.txt和多样性alpha/vegan.txt
    mkdir -p result/alpha
    Rscript ${db}/script/otutab_rare.R --input result/otutab.txt \
      --depth 327 --seed 1 \
      --normalize result/otutab_rare.txt \
      --output result/alpha/vegan.txt
    usearch -otutab_stats result/otutab_rare.txt \
      -output result/otutab_rare.stat
    cat result/otutab_rare.stat


## 6. α多样性 alpha diversity

### 6.1. 计算α多样性 calculate alpha diversity

    # 使用USEARCH计算14种alpha多样性指数(Chao1有错勿用)
    #details in http://www.drive5.com/usearch/manual/alpha_metrics.html
    usearch -alpha_div result/otutab_rare.txt \
      -output result/alpha/alpha.txt

### 6.2. 计算稀释丰富度 calculate rarefaction richness

    #稀释曲线：取1%-100%的序列中OTUs数量，每次无放回抽样
    #Rarefaction from 1%, 2% .. 100% in richness (observed OTUs)-method without_replacement https://drive5.com/usearch/manual/cmd_otutab_subsample.html
    usearch -alpha_div_rare result/otutab_rare.txt \
      -output result/alpha/alpha_rare.txt \
      -method without_replacement
    #预览结果
    head -n2 result/alpha/alpha_rare.txt
    #样本测序量低出现非数值"-"的处理，详见常见问题8
    sed -i "s/-/\t0.0/g" result/alpha/alpha_rare.txt

### 6.3. 筛选高丰度菌 Filter by abundance

    #计算各特征的均值，有组再求分组均值，需根据实验设计metadata.txt修改组列名
    #输入文件为feautre表result/otutab.txt，实验设计metadata.txt
    #输出为特征表按组的均值-一个实验可能有多种分组方式
    #-h显示脚本帮助(参数说明)
    Rscript ${db}/script/otu_mean.R -h
    #scale是否标准化，zoom标准化总和，all输出全部样本均值，type计算类型mean或sum
    Rscript ${db}/script/otu_mean.R --input result/otutab.txt \
      --metadata result/metadata.txt \
      --group Group --thre 0 \
      --scale TRUE --zoom 100 --all TRUE --type mean \
      --output result/otutab_mean.txt
    # 结果为全部和各组均值
    head -n3 result/otutab_mean.txt

    #如以平均丰度>0.1%筛选，可选0.5或0.05，得到每个组的OTU组合，i=3为从第三列开始，要保留all列的话改为i=2,这里为了进行四组比较的韦恩图绘制保留了all列
    awk 'BEGIN{OFS=FS="\t"}{if(FNR==1) {for(i=2;i<=NF;i++) a[i]=$i; print "OTU","Group";} \
        else {for(i=2;i<=NF;i++) if($i>0.1) print $1, a[i];}}' \
        result/otutab_mean.txt > result/alpha/otu_group_exist.txt
    head result/alpha/otu_group_exist.txt
    cut -f 2 result/alpha/otu_group_exist.txt | sort | uniq -c
    # 试一试：不同丰度下各组有多少OTU/ASV
    # 可在 http://ehbio.com/test/venn/ 中绘图并显示各组共有和特有维恩或网络图
    # 也可在 http://www.ehbio.com/ImageGP 绘制Venn、upSetView和Sanky

## 7. β多样性 Beta diversity

    #结果有多个文件，需要目录
    mkdir -p result/beta/
    #基于OTU构建进化树 Make OTU tree, 6.7s
    usearch -cluster_agg result/otus.fa -treeout result/otus.tree
    #1s生成5种距离矩阵：bray_curtis, euclidean, jaccard, manhatten, unifrac
    usearch -beta_div result/otutab_rare.txt -tree result/otus.tree \
    -filename_prefix result/beta/


## 8. 物种注释分类汇总

    #OTU对应物种注释2列格式：去除sintax中置信值，只保留物种注释，替换:为_，删除引号
    cut -f 1,4 result/otus.sintax \
      |sed 's/\td/\tk/;s/:/__/g;s/,/;/g;s/"//g' \
      > result/taxonomy2.txt
    head -n3 result/taxonomy2.txt

    #OTU对应物种8列格式：注意注释是非整齐
    #生成物种表格OTU/ASV中空白补齐为Unassigned
    awk 'BEGIN{OFS=FS="\t"}{delete a; a["k"]="Unassigned";a["p"]="Unassigned";a["c"]="Unassigned";a["o"]="Unassigned";a["f"]="Unassigned";a["g"]="Unassigned";a["s"]="Unassigned";\
      split($2,x,";");for(i in x){split(x[i],b,"__");a[b[1]]=b[2];} \
      print $1,a["k"],a["p"],a["c"],a["o"],a["f"],a["g"],a["s"];}' \
      result/taxonomy2.txt > temp/otus.tax
    sed 's/;/\t/g;s/.__//g;' temp/otus.tax|cut -f 1-8 | \
      sed '1 s/^/OTUID\tKingdom\tPhylum\tClass\tOrder\tFamily\tGenus\tSpecies\n/' \
      > result/taxonomy.txt
    head -n3 result/taxonomy.txt

    #统计门纲目科属，使用 rank参数 p c o f g s，为phylum, class, order, family, genus, species缩写
    mkdir -p result/tax
    for i in p c o f g s;do
      usearch -sintax_summary result/otus.sintax \
      -otutabin result/otutab_rare.txt -rank ${i} \
      -output result/tax/sum_${i}.txt
    done
    sed -i 's/(//g;s/)//g;s/\"//g;s/\#//g;s/\/Chloroplast//g' result/tax/sum_*.txt
    # 列出所有文件
    wc -l result/tax/sum_*.txt
    head -n3 result/tax/sum_g.txt


## 9. 有参定量特征表

    # 比对Greengenes97% OTUs比对，用于PICRUSt/Bugbase功能预测
    mkdir -p result/gg/

    #方法1. usearch比对更快，但文件超限报错选方法2
    # 默认10核以下使用1核，10核以上使用10核
    usearch -otutab temp/filtered.fa -otus ${db}/gg/97_otus.fa \
    	-otutabout result/gg/otutab.txt -threads 4
    # 比对率80.8%, 4核1m, 内存使用731Mb
    head -n3 result/gg/otutab.txt

    # #方法2. vsearch比对，更准更慢，但并行24-96线程更强
    # vsearch --usearch_global temp/filtered.fa --db ${db}/gg/97_otus.fa \
       --otutabout result/gg/otutab.txt --id 0.97 --threads 12
    # 比对率83.83%, 12核4m24s

    #统计
    usearch -otutab_stats result/gg/otutab.txt -output result/gg/otutab.stat
    cat result/gg/otutab.stat


## 10. 空间清理及数据提交

    #删除中间大文件
    rm -rf temp/*.fq

    # 分双端统计md5值，用于数据提交
    cd seq
    md5sum *.fastq > ../result/md5sum.txt
    cat ../result/md5sum.txt

# EasyAmplicon2 base DADA2
## 加载必要包

    library(dada2)
    library(tidyverse)

## 设置路径

    path <- "fastq"  # 修改为你的fastq所在路径
    metadata <- read.delim("metadata.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)
    
    # 获取样本ID
    samples <- metadata$SampleID
    fastq_files <- file.path(path, paste0(samples, ".fastq.gz"))
    names(fastq_files) <- samples
    
    # 创建过滤后文件保存目录
    filt_path <- file.path(path, "filtered")
    dir.create(filt_path, showWarnings = FALSE)
    filt_files <- file.path(filt_path, paste0(samples, "_filt.fastq.gz"))
    
## Step 1: 过滤和修剪（可调整 maxEE 和 minLen）

    out <- filterAndTrim(
      fwd = fastq_files,
      filt = filt_files,
      maxN = 0,
      maxEE = 2,
      minLen = 1000,
      matchIDs = TRUE,
      multithread = TRUE,
      compress = TRUE
    )

## Step 2: 学习错误模型
    err <- learnErrors(filt_files, multithread = TRUE, randomize = TRUE)

## Step 3: 去噪（每个样本单独）
    dada_list <- vector("list", length(samples))
    names(dada_list) <- samples
    
    for (i in seq_along(samples)) {
      cat("正在处理样本:", samples[i], "\n")
      derep <- derepFastq(filt_files[i])
      dada_list[[i]] <- dada(derep, err = err, multithread = TRUE)
    }

## Step 4: 构建ASV表
    seqtab <- makeSequenceTable(dada_list)
    dim(seqtab)

## Step 5: 去嵌合体
    seqtab.nochim <- removeBimeraDenovo(seqtab, method = "consensus", multithread = TRUE)
    dim(seqtab.nochim)

## Step 6: 保存ASV表和序列
    asv_seqs <- colnames(seqtab.nochim)
    asv_headers <- paste0("ASV", seq_along(asv_seqs))
    names(asv_seqs) <- asv_headers

## 保存FASTA
    asv_fasta <- "ASVs.fasta"
    writeLines(paste0(">", asv_headers, "\n", asv_seqs), asv_fasta)

## 保存ASV表
    asv_table <- as.data.frame(t(seqtab.nochim))
    colnames(asv_table) <- samples
    rownames(asv_table) <- asv_headers
    write.csv(asv_table, "ASV_table.csv")

## 保存metadata
    write.csv(metadata, "sample_metadata.csv", row.names = FALSE)

## Step 7: 物种注释（可选，使用silva或GTDB数据库）
    # taxa <- assignTaxonomy(seqtab.nochim, "silva_nr99_v138.1_train_set.fa.gz", multithread = TRUE)

## Step 8: 构建phyloseq对象（可选）
    # library(phyloseq)
    # ps <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE), sample_data(metadata), tax_table(taxa))
    # saveRDS(ps, file = "phyloseq_object.rds")


# R语言多样性和物种组成分析

## 1. Alpha多样性

### 1.1 Alpha多样性箱线图

    # 查看帮助,后续所有R代码均可用此方法查看帮助信息
    Rscript ${db}/script/alpha_div_box.R -h
    # 完整参数，多样性指数可选richness chao1 ACE shannon simpson invsimpson
    Rscript ${db}/script/alpha_div_box.R --input result/alpha/vegan.txt   --metadata result/metadata.txt   --alpha_index richness,chao1,ACE,shannon,simpson,invsimpson   --group Group   --out_prefix result/alpha/vegan
    #--alpha_index控制绘制多样性指数类型，可绘制单个

### 1.2 稀释曲线

    Rscript ${db}/script/alpha_rare_curve.R \
      --input result/alpha/alpha_rare.txt --design result/metadata.txt \
      --group Group --output result/alpha/ \
      --width 120 --height 78

### 1.3 多样性维恩图

    # 支持2-4组比较
    Rscript ${db}/script/venn.R --input result/alpha/otu_group_exist.txt --groups feces,plaque,saliva --output result/alpha/venn.pdf

## 2. Beta多样性

### 2.1 距离矩阵热图pheatmap

    # 以bray_curtis为例，-f输入文件,-h是否聚类TRUE/FALSE,-u/v为宽高英寸
    bash ${db}/script/sp_pheatmap.sh \
      -f result/beta/bray_curtis.txt \
      -H 'TRUE' -u 6 -v 5
    # 添加分组注释，如2，4列的基因型和地点
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

    # 输入文件，选择分组，输出文件
    Rscript ${db}/script/beta_PCoA.R  --input result/otutab.txt   --metadata result/metadata.txt   --group Group   --output result/beta/PCoa.pdf
      
### 2.3 限制性主坐标分析CPCoA

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

### 3.2 弦/圈图circlize

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

    # 以属为例（genus，g），绘制丰度前15的属水平丰度气泡图；输入物种丰度表和样本metadata文件，输出分组物种丰度气泡图
    i=g
    Rscript ${db}/script/tax_bubble.R \
    -i result/tax/sum_${i}.txt \
    -g result/metadata.txt \
    -c Group -w 7 -e 4 -s 15 -n 15 \
    -o result/tax/sum_g.bubble.pdf



# 4、差异比较

## 1. R语言差异分析

### 1.1 差异比较

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

    详见volcano.R

### 1.3 热图

    # 输入compare.R的结果，筛选列数，指定元数据和分组、物种注释，图大小英寸和字号
    bash ${db}/script/compare_heatmap.sh -i result/compare/${compare}.txt -l 7 \
       -d result/metadata.txt -A Group \
       -t result/taxonomy.txt \
       -w 12 -h 20 -s 14 \
       -o result/compare/${compare}


### 1.4 三元图
    Rscript ${db}/script/ternary_plot.R   \
    --input result/tax/sum_p.txt   \
    --metadata result/metadata.txt  \
    --group Group   \
    --taxlevel Phylum \
    --output result/compare/ternary_p.pdf   \
    --topn 10

## 2. STAMP差异分析图

    详见stamp.R

### 2.1 生成输入文件(备选)

    Rscript ${db}/script/format2stamp.R -h
    mkdir -p result/stamp
    Rscript ${db}/script/format2stamp.R --input result/otutab.txt \
      --taxonomy result/taxonomy.txt --threshold 0.01 \
      --output result/stamp/tax
    # 可选Rmd文档见result/format2stamp.Rmd



# 3. 功能预测

## 3.1 PICRUSt2环境导出和导入
    
    # 方法1. 直接安装
    n=picrust2
    conda create -n ${n} -c bioconda -c conda-forge ${n}=2.3.0_b
    # 加载环境
    conda activate ${n}

    # 方法2. 导出安装环境
    cd ~/db/conda/
    # 设置环境名
    n=picrust2
    conda activate ${n}
    # 打包环境为压缩包
    conda pack -n ${n} -o ${n}.tar.gz
    # 导出软件安装列表
    conda env export > ${n}.yml
    # 添加权限，方便下载和别人使用
    chmod 755 ${n}.*
    
    # 方法3. 导入安装环境，如qiime2 humann2 meta(包括picurst)
    n=picrust2
    # 复制安装包，或下载我的环境打包
    wget -c ftp://download.nmdc.cn/tools/conda/${n}.tar.gz
    # 指定安装目录并解压
    condapath=~/miniconda3
    mkdir -p ${condapath}/envs/${n}
    tar -xvzf ${n}.tar.gz -C ${condapath}/envs/${n}
    # 激活环境并初始化
    source ${condapath}/envs/${n}/bin/activate
    conda unpack

## 3.2 PICRUSt2功能预测

    # (可选)PICRUSt2(Linux/Windows下Linux子系统，要求>16GB内存)
    # 安装参考附录6的方式直接下载安装包并解压即可使用
    
    # 加载环境
    conda activate picrust2
    # 进入工作目录，服务器要修改工作目录
    wd=/mnt/d/amplicon/result/picrust2
    mkdir -p ${wd} && cd ${wd}
    # 运行流程，内存15.7GB，耗时12m
    picrust2_pipeline.py -s ../otus.fa -i ../otutab.txt -o ./out -p 8
    # 添加EC/KO/Pathway注释
    cd out
    add_descriptions.py -i pathways_out/path_abun_unstrat.tsv.gz -m METACYC \
    -o pathways_out/path_abun_unstrat_descrip.tsv.gz
    add_descriptions.py -i EC_metagenome_out/pred_metagenome_unstrat.tsv.gz -m EC \
      -o EC_metagenome_out/pred_metagenome_unstrat_descrip.tsv.gz
    add_descriptions.py -i KO_metagenome_out/pred_metagenome_unstrat.tsv.gz -m KO \
      -o KO_metagenome_out/pred_metagenome_unstrat_descrip.tsv.gz 
    # KEGG按层级合并
    db=/mnt/c/EasyMicrobiome/
    zcat KO_metagenome_out/pred_metagenome_unstrat.tsv.gz > KEGG.KO.txt
    python3 ${db}/script/summarizeAbundance.py \
      -i KEGG.KO.txt \
	    -m ${db}/kegg/KO1-4.txt \
	    -c 2,3,4 -s ',+,+,' -n raw \
	    -o KEGG
    # 统计各层级特征数量
    wc -l KEGG*
    #功能预测结果可视化
    # 设置比较组的名称变量，这里比较的是"feces"（粪便）和"plaque"（牙菌斑）
    compare="feces-plaque"
    
    # 第一个R脚本：执行两组间的差异分析
    # 调用R脚本compare.R
    Rscript ${db}/script/compare.R \
    --input KEGG.PathwayL2.raw.txt \
    --design /d/EasyAmplicon_paper_materials/PacBio/result/metadata.txt \
    --group Group \
    --compare ${compare} \
    --threshold 0 \
    --method wilcox \
    --pvalue 0.05 \
    --fdr 0.2 \
    --output ./
    
    # 第二个R脚本：生成层级结构的条形图
    # 调用R脚本compare_hierarchy_facet.R
    Rscript ${db}/script/compare_hierarchy_facet.R \
      --input ${compare}.txt \
      --data MeanA \
      --annotation ../${l}.anno.txt \
      --output ../${compare}.MeanA.bar.pdf






# 4. Evolution进化树

    cd ${wd}
    mkdir -p result/tree
    cd ${wd}/result/tree

## 4.1 筛选高丰度/指定的特征

    #方法1. 按丰度筛选特征，一般选0.001或0.005，且OTU数量在30-150个范围内
    #统计特征表中ASV数量，如总计1234个
    tail -n+2 ../otutab_rare.txt | wc -l
    #按相对丰度0.2%筛选高丰度OTU
    usearch -otutab_trim ../otutab_rare.txt \
        -min_otu_freq 0.002 \
        -output otutab.txt
    #统计筛选OTU表特征数量，总计~105个
    tail -n+2 otutab.txt | wc -l

    #方法2. 按数量筛选
    # #按丰度排序，默认由大到小
    # usearch -otutab_sortotus ../otutab_rare.txt  \
    #     -output otutab_sort.txt
    # #提取高丰度中指定Top数量的OTU ID，如Top100,
    # sed '1 s/#OTU ID/OTUID/' otutab_sort.txt \
    #     | head -n101 > otutab.txt

    #修改特征ID列名
    sed -i '1 s/#OTU ID/OTUID/' otutab.txt
    #提取ID用于提取序列
    cut -f 1 otutab.txt > otutab_high.id

    # 筛选高丰度菌/指定差异菌对应OTU序列
    usearch -fastx_getseqs ../otus.fa -labels otutab_high.id \
        -fastaout otus.fa
    head -n 2 otus.fa

    ## 筛选OTU对物种注释
    awk 'NR==FNR{a[$1]=$0} NR>FNR{print a[$1]}' ../taxonomy.txt \
        otutab_high.id > otutab_high.tax

    #获得OTU对应组均值，用于样本热图
    #依赖之前otu_mean.R计算过按Group分组的均值
    awk 'NR==FNR{a[$1]=$0} NR>FNR{print a[$1]}' ../otutab_mean.txt otutab_high.id \
        | sed 's/#OTU ID/OTUID/' > otutab_high.mean
    head -n3 otutab_high.mean

    #合并物种注释和丰度为注释文件
    cut -f 2- otutab_high.mean > temp
    paste otutab_high.tax temp > annotation.txt
    head -n 3 annotation.txt

## 4.2 构建进化树

    # 起始文件为 result/tree目录中 otus.fa(序列)、annotation.txt(物种和相对丰度)文件
    # Muscle软件进行序列对齐，3s
    muscle -in otus.fa -out otus_aligned.fas

    ### 方法1. 利用IQ-TREE快速构建ML进化树，2m
    # rm -rf iqtree
    # mkdir -p iqtree
    # iqtree -s otus_aligned.fas \
        -bb 1000 -redo -alrt 1000 -nt AUTO \
        -pre iqtree/otus

    ### 方法2. FastTree快速建树(Linux)
    # 注意FastTree软件输入文件为fasta格式的文件，而不是通常用的Phylip格式。输出文件是Newick格式。
    # 该方法适合于大数据，例如几百个OTUs的系统发育树！
    # Ubuntu上安装fasttree可以使用`apt install fasttree`
    fasttree -gtr -nt otus_aligned.fas > otus.nwk

## 4.3 进化树可视化





# 常见问题


## 1. 数据过大无法使用usearch聚类或去噪,替换vsearch

    # 仅限usearch免费版受限时，可通过提高minuniquesize参数减少非冗余数据量。OTU/ASV过万下游分析等待时间过长，确保OTU/ASV数据小于5000，一般不会受限，而且也有利于下游开展快速分析。

    # 备选vsearch聚类生成OTU，但无自动de novo去嵌合功能。输入2155条序列，聚类后输出661。

    cd /c/amplicon/FAQ/03feature
    # 重命名relabel、按相似id=97%聚类，不屏蔽qmask
    # 记录输入sizein和输出频率sizeout
    vsearch --cluster_size uniques.fa  \
     --relabel OTU_ --id 0.97 \
     --qmask none --sizein --sizeout \
     --centroids otus_raw.fa 


    # 再de novo去嵌合。55个嵌合，606个非嵌合。把OTU_1都去除了，没有Usearch内置去嵌合的方法合理。

    # 自身比对去嵌合
    vsearch --uchime_denovo otus_raw.fa \
        --nonchimeras otus.fa
    # 删除序列频率
    sed -i 's/;.*//' otus.fa

## 2. 读长计数(Read counts)标准化为相对丰度

    cd /c/amplicon/FAQ/04norm
    # 求取各个OTU在样品中的丰度频率(标准化为总和1)
    usearch -otutab_counts2freqs otutab.txt \
        -output otutab_freq.txt

## 3. 运行R提示Permission denied
 
    # 例如write.table保存表时，报错信息示例如下：意思是写入文件无权限，一般为目标文件正在被打开，请关闭相关文件后重试

    Error in file(file, ifelse(append, "a", "w")) :
    Calls: write.table -> file
    : Warning message:
    In file(file, ifelse(append, "a", "w")) :
      'result/raw/otutab_nonBac.txt': Permission denied

## 4. 文件批量命名

    # 如我们有文件A1和A2，编写一个样本名对应目标名的表格metadata.txt，检查样本名是否唯一，使用awk进行批量改名

    cd /c/amplicon/FAQ/06rename
    # (可选)快速生成文件列表，用于编辑metadata.txt，如A1.fq修改为WT1.fastq，以此类推，参考metadata.bak.txt
    ls *R1_001.fastq > metadata.txt
    # 编辑列表，第二名为最终命名，确定名称唯一
    # 转换行尾换行符
    sed -i 's/\r//' metadata.txt
    # 检查手动命名列2是否唯一
    cut -f 2 metadata.txt|wc -l
    cut -f 2 metadata.txt|sort|uniq|wc -l
    # 如果两次结果一致，则命名非冗余
    # 可选移动mv，复制cp，硬链ln，或软链ln -s
    # 此处使用复制cp
    awk '{system("cp "$1" "$2"_R1.fastq")}' metadata.txt
    awk '{system("cp "$1" "$2"_R2.fastq")}' metadata.txt
    
    mv /raw/DC*fastq > /seq

## 5. Rstudio中Terminal找不到Linux命令

    # 需要把 C:\Program Files\Git\usr\bin 目录添加到系统环境变量
    # 文件资源管理器——此电脑——属性——高级系统设置——环境变量——系统变量——Path——编辑——新建——填写“C:\Program Files\Git\usr\bin”——确定——确定——确定
    # 注意win10系统是一个目录一行；win7中多个目录用分号分隔，注意向后添加目录

## 6. usearch -alpha_div_rare结果前两行出现“-”

    #问题：抽样0时补“-”，且缺失制表符

    #处理：替换“-”为"制作符\t+0"即可恢复

    cd /c/amplicon/FAQ/08rare
    sed "s/-/\t0.0/g" alpha_rare_wrong.txt\
        > alpha_rare.txt

## 7. 物种注释otus.sintax方向全为“-”，需要序列取反向互补

    #是原始序列方向错误，将filtered.fa序列需要取反向互补。再从头开始分析

    cd /c/amplicon/FAQ/09revcom
    vsearch --fastx_revcomp filtered_RC.fa \
      --fastaout filtered.fa

## 8. windows换行符查看和删除

    #Windows换行符为换行($)+^M，等于Linux换行+mac换行。分析数据中以linux格式为通用标准，因此windows中如excel编写并存为文本文件(制表符分隔)(*.txt)的表格，行尾有不可见的^M符号，导致分析出错。可使用cat -A命令查看此符号，可用sed删除。

    cd /c/amplicon/FAQ/10^M
    # 查看行尾是否有^M
  	cat -A metadata.txt
  	# 删除^M，并写入新文件
  	sed 's/\r//' metadata.txt > metadata.mod.txt
  	# 检查是否成功
  	cat -A metadata.mod.txt
  	
  	# 直接原文件删除
  	sed -i 's/\r//' metadata.txt

## 9. UNITE数据库分析报错

    #USEARCH使用UNITE下载的utax数据库，提示各种错误

    cd /c/amplicon/FAQ/11unite
    # 解压Unite的useach使用物种注释库
    gunzip -c utax_reference_dataset_all_04.02.2020.fasta.gz > unite.fa
    # 对ITS序列进行注释，默认阈值0.8
    usearch --sintax  otus.fa \
      --db unite.fa \
      --tabbedout otus.sintax --strand plus
       --sintax_cutoff 0.6

    #报错信息如下：
    ---Fatal error---
    Missing x: in name >JN874928|SH1144646.08FU;tax=d:Metazoa,p:Cnidaria,c:Hydrozoa,o:Trachylina,f:,g:Craspedacusta,s:Craspedacusta_sowerbii_SH1144646.08FU;
    “Unprintable ASCII character no 195 on or right before line 236492”
    
    # 分析原因为分类级存在空缺。可用sed补全即可解决
    # 分类级存在空缺，sed补全
    sed -i 's/,;/,Unnamed;/;s/:,/:Unnamed,/g' unite.fa
    # 再运行前面usearch --sintax命令
    #注：vsearch有问题，推荐用usearch，结尾添加--strand plus才能成功运行

## 10. Windows的Linux子系统本地安装qiime2

    # 详见 qiime2/pipeline_qiime2.sh
    n=qiime2-2023.2
    # 安装包下载链接 
    wget -c ftp://download.nmdc.cn/tools/conda/${n}.tar.gz
    # 新环境安装
    mkdir -p ~/miniconda3/envs/${n}
    tar -xzf ${n}.tar.gz -C ~/miniconda3/envs/${n}
    # 激活并初始化环境
    conda activate ${n}
    conda unpack
