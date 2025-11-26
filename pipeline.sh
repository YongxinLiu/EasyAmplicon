[TOC]

# 易扩增子EasyAmplicon 2

    # 作者 Authors: 刘永鑫(Yong-Xin Liu), 陈同(Tong Chen), 罗豪(Hao Luo)等
    # 版本 Version: v2.00
    # 更新 Update: 2026-04-11
    # 系统要求 System requirement: Windows 10+ / Mac OS 10.12+ / Ubuntu 20.04+
    # 引文 Reference: Liu, et al. 2023. EasyAmplicon: An easy-to-use, open-source, reproducible, and community-based pipeline for amplicon data analysis in microbiome research. iMeta 2: e83. https://doi.org/10.1002/imt2.83
    # Yousuf, et al. 2024. Unveiling microbial communities with EasyAmplicon: A user-centric guide to perform amplicon sequencing data analysis. iMetaOmics 1: e42. https://doi.org/10.1002/imo2.42
    # Luo, et al. 2025. EasyAmplicon 2: Expanding PacBio and Nanopore Long Amplicon Sequencing Analysis Pipeline for Microbiome. Advanced Science 12: https://doi.org/10.1002/advs.202512447

    # 设置工作(work directory, wd)和软件数据库(database, db)目录
    # 添加环境变量，并进入工作目录 Add environmental variables and enter work directory
    # **每次打开Rstudio必须运行下面4行 Run it**，可选替换${db}为EasyMicrobiome安装位置
    wd=/d/amplicon
    db=/d/EasyMicrobiome
    PATH=$PATH:${db}/win
    mkdir -p ${wd}
    cd ${wd}


## 1. 起始文件 start files

    # 1. 分析流程pipeline.sh
    # 2. 样本元信息metadata.txt，保存于result目录
    # 3. 测序数据fastq文件保存于seq目录，通常以`.fq.gz`结尾，每个样品一对文件
    # 4. 创建临时文件存储目录，分析结束可删除
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
    cat -A result/metadata.txt | head -n3

### 1.2. 测序数据 sequencing data

    # # 本段代码可在RStudio中Ctrl + Shift + C 取消注释“#”后运行
    # # (可选)下载测序数据，按GSA的CRA(批次)和CRR(样品)编号下载数据
    # # 示例下载单个文件并改名
    # mkdir -p seq
    # wget -c ftp://download.big.ac.cn/gsa/CRA002352/CRR117575/CRR117575_f1.fq.gz -O seq/KO1_1.fq.gz
    # # 按实验设计编号批量下载并改名
    # awk '{system("wget -c ftp://download.big.ac.cn/gsa/"$5"/"$6"/"$6"_f1.fq.gz -O seq/"$1"_1.fq.gz")}' \
    #     <(tail -n+2 result/metadata.txt)
    # awk '{system("wget -c ftp://download.big.ac.cn/gsa/"$5"/"$6"/"$6"_r2.fq.gz -O seq/"$1"_2.fq.gz")}' \
    #     <(tail -n+2 result/metadata.txt)

    # 公司返回的测序结果，通常为一个样品一对fq/fastq.gz格式压缩文件
    # 文件名与样品名务必对应：不一致时手工修改，批量改名见"常见问题6"
    # 如果测序数据是.gz的压缩文件，有时需要使用gunzip解压后使用，vsearch通常可直接读取压缩文件
    # gunzip seq/*.gz
    # zless按页查看压缩文件，空格翻页、q退出；head默认查看前10行，-n指定行
    ls -sh seq/
    zless seq/KO1_1.fq.gz | head -n4 
    # 每行太长，指定查看每行的1-60个字符
    zless seq/KO1_1.fq | head | cut -c 1-60
    # 统计测序数据，依赖seqkit程序
    seqkit stat seq/KO1_1.fq.gz
    # 批量统计测序数据并汇总表
    seqkit stat seq/*.fq.gz > result/seqkit.txt
    head result/seqkit.txt

### 1.3. 流程和数据库 pipeline & database

    # 数据库第一次使用必须解压，以后可跳过此段

    # usearchs可用16S/18S/ITS数据库：RDP, SILVA和UNITE，本地文件位置 ${db}/usearch/
    # usearch数据库database下载页: http://www.drive5.com/usearch/manual/sintax_downloads.html
    # 解压16S RDP数据库，gunzip解压缩，seqkit stat统计
    # 保留原始压缩文件
    gunzip -c ${db}/usearch/rdp_16s_v18.fa.gz > ${db}/usearch/rdp_16s_v18.fa
    seqkit stat ${db}/usearch/rdp_16s_v18.fa # 2.1万条序列
    # 解压ITS UNITE数据库，需自行从官网https://unite.ut.ee 或网盘db/amplicon/usearch中下载
    # gunzip -c ${db}/usearch/utax_reference_dataset_all_19.02.2025.fasta.gz >${db}/usearch/unite.fa
    # seqkit stat ${db}/usearch/unite.fa # 26.6万
    # Greengene数据库用于功能注释: ftp://greengenes.microbio.me/greengenes_release/gg_13_5/gg_13_8_otus.tar.gz
    # 默认解压会删除原文件，-c指定输出至屏幕，> 写入新文件(可改名)
    gunzip -c ${db}/gg/97_otus.fasta.gz > ${db}/gg/97_otus.fa
    seqkit stat ${db}/gg/97_otus.fa
    # 根据序列和物种注释制作usearch格式数据库，参考 https://mp.weixin.qq.com/s/3FCId5KuE6bigH4G2YCJVg

## 2. 序列合并和重命名 reads merge and rename

### 2.1 合并双端序列并按样品重命名 Merge pair-end reads and rename

    #测试.以WT1单样品合并为例
    #time统计计算时间，real是物理时间，user是计算时间，sys是硬件等待时间
    time vsearch --fastq_mergepairs seq/WT1_1.fq.gz \
      --reverse seq/WT1_2.fq.gz \
      --fastqout temp/WT1.merged.fq \
      --relabel WT1.
    head temp/WT1.merged.fq

    #依照实验设计批处理并合并
    #tail -n+2去表头，cut -f1取第一列，获得样本列表；18个样本x1.5万对序列合并8s
    #Win下复制Ctrl+C为Linux下中止，为防止异常中断，结尾添加&转后台，无显示后按回车继续
    
    # 一部分电脑 rush 不支持，运行时调度失败，请使用 for 循环部分
    # for 循环部分是放入后台运行的，点完 run 之后，看上去程序已运行完，实际没运行完，而是正在运行中。
    # 不要急于运行后面的程序。
    # 之前课程，有发现每次运行结果都不一样，就是因为 for 循环部分没运行完，只生成了部分数据，导致后面
    # 每个样品 reads 数每次运行都会不一致。
    #方法1.for循环顺序处理
    time for i in `tail -n+2 result/metadata.txt|cut -f1`;do
      vsearch --fastq_mergepairs seq/${i}_1.fq.gz --reverse seq/${i}_2.fq.gz \
      --fastqout temp/${i}.merged.fq --relabel ${i}.
    done &

    # 一部分电脑 rush 不支持，运行时调度失败，请使用 for 循环部分
    #方法2.rush并行处理，任务数jobs(j),2可加速1倍4s；建议设置2-4
    time tail -n+2 result/metadata.txt | cut -f 1 | \
     rush -j 2 "vsearch --fastq_mergepairs seq/{}_1.fq.gz --reverse seq/{}_2.fq.gz \
      --fastqout temp/{}.merged.fq --relabel {}."
    # 检查最后一个文件前10行中样本名
    head temp/`tail -n+2 result/metadata.txt | cut -f 1 | tail -n1`.merged.fq | grep ^@
    
    ##方法3.不支持压缩文件时解压再双端合并
     # time tail -n+2 result/metadata.txt | cut -f 1 | \
     #   rush -j 1 "vsearch --fastq_mergepairs <(zcat seq/{}_1.fq.gz) --reverse <(zcat seq/{}_2.fq.gz) \
     #    --fastqout temp/{}.merged.fq --relabel {}."
     #  time for i in `tail -n+2 result/metadata.txt|cut -f1`;do
     #     vsearch --fastq_mergepairs <(zcat seq/${i}_1.fq.gz) --reverse <(zcat seq/${i}_2.fq.gz) \
     #     --fastqout temp/${i}.merged.fq --relabel ${i}.
     #   done &
      
### 2.2 (可选)单端文件改名 Single-end reads rename

    # # 单个序列改名示例
    # i=WT1
    # gunzip -c seq/${i}_1.fq.gz > seq/${i}.fq
    # usearch -fastx_relabel seq/${i}.fq -fastqout temp/${i}.merged.fq -prefix ${i}.
    # 
    # # 批量改名，需要有单端fastq文件，且解压(usearch不支持压缩格式)
    # gunzip seq/*.gz
    # time for i in `tail -n+2 result/metadata.txt|cut -f1`;do
    #   usearch -fastx_relabel seq/${i}.fq -fastqout temp/${i}.merged.fq -prefix ${i}.
    # done &
    # # vsearch大数据方法参考“常见问题2”

### 2.3 改名后序列整合 integrate renamed reads

    #合并所有样品至同一文件
    cat temp/*.merged.fq > temp/all.fq
    #查看文件大小223M，软件不同版本结果略有差异
    ls -lsh temp/all.fq
    # 查看序列名，“.”之前是否为样本名，样本名绝不允许有点 (".")
    # 样本名有点 (.) 的一个显著特征是生成的特征表会很大，特征表里面列很多，导致后面分析出现内存不足。
    # 后面分析获得特征表后要看一眼有没有问题，遇到内存不足问题，也要回头来排查。
    head -n 6 temp/all.fq|cut -c1-60


## 3. 切除引物与质控 Cut primers and quality filter

    # 左端10bp标签+19bp上游引物V5共为29，右端V7为18bp下游引物
    # Cut barcode 10bp + V5 19bp in left and V7 18bp in right
    # 务必清楚实验设计和引物长度，引物已经去除可填0，27万条序列14s
    time vsearch --fastx_filter temp/all.fq \
      --fastq_stripleft 29 --fastq_stripright 18 \
      --fastq_maxee_rate 0.01 \
      --fastaout temp/filtered.fa
    # 查看文件了解fa文件格式
    head temp/filtered.fa


## 4. 去冗余挑选OTU/ASV Dereplicate and cluster/denoise

### 4.1 序列去冗余 Dereplicate

    # 并添加miniuniqusize最小为8或1/1M，去除低丰度噪音并增加计算速度
    # -sizeout输出丰度, --relabel必须加序列前缀更规范, 1s
    vsearch --derep_fulllength temp/filtered.fa \
      --minuniquesize 10 --sizeout --relabel Uni_ \
      --output temp/uniques.fa 
    #高丰度非冗余序列非常小(500K~5M较适合)，名称后有size和频率
    ls -lsh temp/uniques.fa
    # Uni_1;size=6423  - 去冗余后序列的名字 Uni_1；该序列在所有样品测序数据中出现 6423 次
    # 为出现最多的序列。
    head -n 2 temp/uniques.fa

### 4.2 聚类OTU/去噪ASV Cluster or denoise

    #有两种方法：推荐unoise3去噪获得单碱基精度ASV，备选传统的97%聚类OTU(属水平精度)
    #usearch两种特征挑选方法均自带de novo去嵌合体
    #-minsize二次过滤，控制OTU/ASV数量至1-5千，方便下游统计分析

    #方法1. 97%聚类OTU，适合大数据/ASV规律不明显/reviewer要求
    #结果耗时1s, 产生508 OTUs, 去除126 chimeras
    # usearch -cluster_otus temp/uniques.fa -minsize 10 \
    #  -otus temp/otus.fa \
    #  -relabel OTU_

    #方法2. ASV去噪 Denoise: predict biological sequences and filter chimeras
    #6s, 1530 good, 41 chimeras, 序列百万条可能需要几天/几周
    usearch -unoise3 temp/uniques.fa -minsize 10 \
      -zotus temp/zotus.fa
    #修改序列名：Zotu为改为ASV方便识别
    sed 's/Zotu/ASV_/g' temp/zotus.fa > temp/otus.fa
    head -n 2 temp/otus.fa

    #方法3. 数据过大无法使用usearch时，备选vsearch方法见"常见问题3"

### 4.3 基于参考去嵌合 Reference-based chimera detect

    # 不推荐，容易引起假阴性，因为参考数据库无丰度信息
    # 而de novo时要求亲本丰度为嵌合体16倍以上防止假阴性
    # 因为已知序列不会被去除，数据库选择越大越合理，假阴性率最低
    mkdir -p result/raw

    # 方法1. vsearch+rdp去嵌合(快但容易假阴性)
    vsearch --uchime_ref temp/otus.fa \
      -db ${db}/usearch/rdp_16s_v18.fa \
      --nonchimeras result/raw/otus.fa
    # RDP: 7s, 143 (9.3%) chimeras; SILVA：9m, 151 (8.4%) chimeras
    # Win vsearch结果添加了windows换行符^M需删除，mac不要执行此命令
    sed -i 's/\r//g' result/raw/otus.fa

    # 方法2. 不去嵌合
    # cp -f temp/otus.fa result/raw/otus.fa


## 5. 特征表构建和筛选 Feature table create and filter

    # OTU和ASV统称为特征(Feature)，它们的区别是：
    # OTU通常按97%聚类后挑选最高丰度或中心的代表性序列；
    # ASV是基于序列进行去噪(排除或校正错误序列，并挑选丰度较高的可信序列)作为代表性序列

### 5.1 生成特征表

    # id(1)：100%相似度比对49.45%序列，1m50s
    # id(0.97)：97%相似度比对83.66%序列，1m10s(更高数据使用率，更快)
    time vsearch --usearch_global temp/filtered.fa \
      --db result/raw/otus.fa \
      --id 0.97 --threads 4 \
    	--otutabout result/raw/otutab.txt 
    #212862 of 268019 (79.42%)可比对
    # vsearch结果windows用户删除换行符^M校正为标准Linux格式
    sed -i 's/\r//' result/raw/otutab.txt
    head -n6 result/raw/otutab.txt | cut -f 1-6 |cat -A
    # csvtk统计表行列
    # 这里一定看好列数，是不是等于你的样品数；如果不等，一般是样品命名存在问题，具体看上面解释
    csvtk -t stat result/raw/otutab.txt


### 5.2 物种注释和去除质体和非细菌 Remove plastid and non-Bacteria

    # 物种注释-去除质体和非细菌/古菌并统计比例(可选)
    # RDP物种注释(rdp_16s_v18)更快，但缺少完整真核来源数据,可能不完整，耗时15s;
    # SILVA数据库(silva_16s_v123.fa)更好注释真核、质体序列，极慢耗时3h起
    # 置信阈值通常0.6/0.8，vserch最低0.1/usearch可选0输出最相似物种注释用于观察潜在分类
    vsearch --sintax result/raw/otus.fa \
      --db ${db}/usearch/rdp_16s_v18.fa \
      --sintax_cutoff 0.1 \
      --tabbedout result/raw/otus.sintax 
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
      --depth 10000 --seed 1 \
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

    #如以平均丰度>0.1%筛选，可选0.5或0.05，得到每个组的OTU组合
    awk 'BEGIN{OFS=FS="\t"}{if(FNR==1) {for(i=2;i<=NF;i++) a[i]=$i; print "OTU","Group";} \
        else {for(i=2;i<=NF;i++) if($i>0.05) print $1, a[i];}}' \
        result/otutab_mean.txt > result/alpha/otu_group_exist.txt
    head result/alpha/otu_group_exist.txt
    cut -f 2 result/alpha/otu_group_exist.txt | sort | uniq -c
    # 试一试：不同丰度下各组有多少OTU/ASV
    # 可在 http://ehbio.com/test/venn/ 中绘图并显示各组共有和特有维恩或网络图
    # 也可在 https://www.bic.ac.cn/ImageGP 或 https://www.bic.ac.cn/BIC/#/ 绘制Venn、upSetView和Sanky

## 7. β多样性 Beta diversity

    #结果有多个文件，需要目录
    mkdir -p result/beta/
    #基于OTU构建进化树 Make OTU tree, 4s
    usearch -cluster_agg result/otus.fa -treeout result/otus.tree
    #生成5种距离矩阵：bray_curtis, euclidean, jaccard, manhatten, unifrac
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

    #统计门纲目科属，使用 rank参数 p c o f g，为phylum, class, order, family, genus缩写
    mkdir -p result/tax
    for i in p c o f g;do
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

    # usearch比对更快，但文件超限报错选附录14 vsearch比对
    # 默认10核以下使用1核，10核以上使用10核
    usearch -otutab temp/filtered.fa -otus ${db}/gg/97_otus.fa \
    	-otutabout result/gg/otutab.txt -threads 4
    # 比对率80.0%, 1核11m，4核3m，10核2m，内存使用743Mb
    head -n3 result/gg/otutab.txt

    #统计
    usearch -otutab_stats result/gg/otutab.txt -output result/gg/otutab.stat
    cat result/gg/otutab.stat


## 10. 空间清理及数据提交

    #删除中间大文件
    rm -rf temp/*.fq

    # 分双端统计md5值，用于数据提交
    cd seq
    md5sum *_1.fq.gz > md5sum1.txt
    md5sum *_2.fq.gz > md5sum2.txt
    paste md5sum1.txt md5sum2.txt | awk '{print $2"\t"$1"\t"$4"\t"$3}' | sed 's/*//g' > ../result/md5sum.txt
    rm md5sum*
    cd ..
    cat result/md5sum.txt
    

# R语言多样性和物种组成分析

## 1. Alpha多样性

### 1.1 Alpha多样性箱线图

    # 在线绘图平台 https://www.bic.ac.cn/BIC/#/ 提供更多定制参数和绘制灵活性
    # 查看帮助
    Rscript ${db}/script/alpha_boxplot.R -h
    # 完整参数，多样性指数可选richness chao1 ACE shannon simpson invsimpson
    Rscript ${db}/script/alpha_boxplot.R --alpha_index richness \
      --input result/alpha/vegan.txt --design result/metadata.txt \
      --group Group --output result/alpha/ \
      --width 89 --height 59
    # 使用循环绘制6种常用指数
    for i in `head -n1 result/alpha/vegan.txt|cut -f 2-`;do
      Rscript ${db}/script/alpha_boxplot.R --alpha_index ${i} \
        --input result/alpha/vegan.txt --design result/metadata.txt \
        --group Group --output result/alpha/ \
        --width 89 --height 59
    done
    mv alpha_boxplot_TukeyHSD.txt result/alpha/
    # 每个指数下两两比较的P值
    cat result/alpha/alpha_boxplot_TukeyHSD.txt

    # Alpha多样性柱状图+标准差
    Rscript ${db}/script/alpha_barplot.R --alpha_index richness \
      --input result/alpha/vegan.txt --design result/metadata.txt \
      --group Group --output result/alpha/ \
      --width 89 --height 59

### 1.2 稀释曲线

    Rscript ${db}/script/alpha_rare_curve.R \
      --input result/alpha/alpha_rare.txt --design result/metadata.txt \
      --group Group --output result/alpha/ \
      --width 120 --height 59

### 1.3 多样性维恩图

    # 三组比较:-f输入文件,-a/b/c/d/g分组名,-w/u为宽高英寸,-p输出文件名后缀
    bash ${db}/script/sp_vennDiagram.sh \
      -f result/alpha/otu_group_exist.txt \
      -a WT -b KO -c OE \
      -w 3 -u 3 \
      -p WT_KO_OE
    # 四组比较，图和代码见输入文件目录，运行目录为当前项目根目录
    bash ${db}/script/sp_vennDiagram.sh \
      -f result/alpha/otu_group_exist.txt \
      -a WT -b KO -c OE -d All \
      -w 3 -u 3 \
      -p WT_KO_OE_All
    # EVenn在线绘制维恩图 https://www.ehbio.com/test/venn 

## 2. Beta多样性

### 2.1 距离矩阵热图pheatmap

    # 在线绘图平台 https://www.bic.ac.cn/BIC/#/ 提供更多定制参数和绘制灵活性
    # 添加分组注释，如2，4列的基因型和地点
    cut -f 1-2 result/metadata.txt > temp/group.txt
    # 以bray_curtis为例，-f输入文件,-h是否聚类TRUE/FALSE,-u/v为宽高英寸
    # -P添加行注释文件，-Q添加列注释
    bash ${db}/script/sp_pheatmap.sh \
      -f result/beta/bray_curtis.txt \
      -H 'TRUE' -u 6.9 -v 5.6 \
      -P temp/group.txt -Q temp/group.txt

### 2.2 主坐标分析PCoA

    # 在线绘图平台 https://www.bic.ac.cn/BIC/#/ 提供更多定制参数和绘制灵活性
    # 输入文件，选择分组，输出文件，图片尺寸mm，统计见beta_pcoa_stat.txt
    Rscript ${db}/script/beta_pcoa.R \
      --input result/beta/bray_curtis.txt --design result/metadata.txt \
      --group Group --label FALSE --width 89 --height 59 \
      --output result/beta/bray_curtis.pcoa.pdf
    # 添加样本标签 --label TRUE
    Rscript ${db}/script/beta_pcoa.R \
      --input result/beta/bray_curtis.txt --design result/metadata.txt \
      --group Group --label TRUE --width 89 --height 59 \
      --output result/beta/bray_curtis.pcoa.label.pdf
    mv beta_pcoa_stat.txt result/beta/
    # 查看组间两两比较P值，概率事件每次略有不同
    cat result/beta/beta_pcoa_stat.txt
      
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

    # 在线绘图平台 https://www.bic.ac.cn/BIC/#/ 提供更多定制参数和绘制灵活性
    # 以门(p)水平为例，结果包括output.sample/group.pdf两个文件
    Rscript ${db}/script/tax_stackplot.R \
      --input result/tax/sum_p.txt --design result/metadata.txt \
      --group Group --color ggplot --legend 7 --width 89 --height 59 \
      --output result/tax/sum_p.stackplot
    # 修改颜色--color ggplot, manual1(22), Paired(12) or Set3(12)
    Rscript ${db}/script/tax_stackplot.R \
      --input result/tax/sum_p.txt --design result/metadata.txt \
      --group Group --color Set3 --legend 10 --width 181 --height 119 \
      --output result/tax/sum_p.stackplotSet3
      
    # 批量绘制输入包括p/c/o/f/g共5级
    for i in p c o f g; do
    Rscript ${db}/script/tax_stackplot.R \
      --input result/tax/sum_${i}.txt --design result/metadata.txt \
      --group Group --output result/tax/sum_${i}.stackplot \
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

### 3.3 树图treemap(参考)

    # 多层级包含物种关系，输入特征表和物种注释，输出树图
    # 指定包含特征数量和图片宽高，100个ASV耗时12s
    Rscript ${db}/script/tax_maptree.R \
      --input result/otutab.txt --taxonomy result/taxonomy.txt \
      --output result/tax/tax_maptree.pdf \
      --topN 100 --width 183 --height 118



# 24、差异比较

## 1. R语言差异分析

### 1.1 差异比较

    # Error in file(file, ifelse(append, "a", "w")),输出目录不存在，创建目录即可
    mkdir -p result/compare/
    # 输入特征表、元数据；指定分组列名、比较组和丰度
    # 选择方法 wilcox/t.test/edgeR、pvalue和fdr和输出目录
    compare="KO-WT"
    Rscript ${db}/script/compare.R \
      --input result/otutab.txt --design result/metadata.txt \
      --group Group --compare ${compare} --threshold 0.1 \
      --method edgeR --pvalue 0.05 --fdr 0.2 \
      --output result/compare/

### 1.2 火山图

    # 输入compare.R的结果，输出火山图带数据标签，可指定图片大小
    Rscript ${db}/script/compare_volcano.R \
      --input result/compare/${compare}.txt \
      --output result/compare/${compare}.volcano.pdf \
      --width 89 --height 59

### 1.3 热图

    # 输入compare.R的结果，筛选列数，指定元数据和分组、物种注释，图大小英寸和字号
    bash ${db}/script/compare_heatmap.sh -i result/compare/${compare}.txt -l 7 \
       -d result/metadata.txt -A Group \
       -t result/taxonomy.txt \
       -w 8 -h 5 -s 7 \
       -o result/compare/${compare}

### 1.4 曼哈顿图

    # i差异比较结果,t物种注释,p图例,w宽,v高,s字号,l图例最大值
    # 图例显示不图，可增加高度v为119+即可，后期用AI拼图为KO-WT.heatmap.emf
    bash ${db}/script/compare_manhattan.sh -i result/compare/${compare}.txt \
       -t result/taxonomy.txt \
       -p result/tax/sum_p.txt \
       -w 183 -v 59 -s 7 -l 10 \
       -o result/compare/${compare}.manhattan.p.pdf
    # 上图只有6个门，切换为纲c和-L Class展示细节
    bash ${db}/script/compare_manhattan.sh -i result/compare/${compare}.txt \
       -t result/taxonomy.txt \
       -p result/tax/sum_c.txt \
       -w 183 -v 59 -s 7 -l 10 -L Class \
       -o result/compare/${compare}.manhattan.c.pdf
    # 显示完整图例，再用AI拼图
    bash ${db}/script/compare_manhattan.sh -i result/compare/${compare}.txt \
       -t result/taxonomy.txt \
       -p result/tax/sum_c.txt \
       -w 183 -v 149 -s 7 -l 10 -L Class \
       -o result/compare/${compare}.manhattan.c.legend.pdf

### 1.5 单个特征的绘制

    # 筛选显示差异ASV，按KO组丰度降序列，取ID展示前10
    awk '$4<0.05' result/compare/KO-WT.txt | sort -k7,7nr | cut -f1 | head
    # 差异OTU细节展示
    Rscript ${db}/script/alpha_boxplot.R --alpha_index ASV_2 \
      --input result/otutab.txt --design result/metadata.txt \
      --transpose TRUE --scale TRUE \
      --width 89 --height 59 \
      --group Group --output result/compare/feature_ 
    # ID不存在会报错： Error in data.frame(..., check.names = FALSE) : 参数值意味着不同的行数: 0, 18  Calls: alpha_boxplot -> cbind -> cbind -> data.frame
    
    # 指定某列排序：按属丰度均值All降序
    csvtk -t sort -k All:nr result/tax/sum_g.txt | head
    # 差属细节展示
    Rscript ${db}/script/alpha_boxplot.R --alpha_index Lysobacter \
      --input result/tax/sum_g.txt --design result/metadata.txt \
      --transpose TRUE \
      --width 89 --height 59 \
      --group Group --output result/compare/feature_

### 1.5 三元图

  #参考示例见：result\compare\ternary\ternary.Rmd 文档
  #备选教程[246.三元图的应用与绘图实战](https://mp.weixin.qq.com/s/3w3ncpwjQaMRtmIOtr2Jvw)


## 2. STAMP输入文件准备

### 2.1 生成输入文件

    Rscript ${db}/script/format2stamp.R -h
    mkdir -p result/stamp
    Rscript ${db}/script/format2stamp.R --input result/otutab.txt \
      --taxonomy result/taxonomy.txt --threshold 0.01 \
      --output result/stamp/tax
    # 可选Rmd文档见result/format2stamp.Rmd

### 2.2 绘制扩展柱状图和表

    compare="KO-WT"
    # 替换ASV(result/otutab.txt)为属(result/tax/sum_g.txt)
    # 代码优化：补充图例、实验vs对照，柱从上到下顺序，目前为WT上KO下
    Rscript ${db}/script/compare_stamp.R \
      --input result/stamp/tax_5Family.txt --metadata result/metadata.txt \
      --group Group --compare ${compare} --threshold 0.1 \
      --method "t.test" --pvalue 0.05 --fdr "none" \
      --width 280 --height 159 \
      --output result/stamp/${compare}
    # 可选Rmd文档见result/CompareStamp.Rmd

## 3. LEfSe输入文件准备

    ### 3.1. 命令行生成文件
    # 可选命令行生成输入文件
    Rscript ${db}/script/format2lefse.R -h
    mkdir -p result/lefse
    # threshold控制丰度筛选以控制作图中的枝数量
    Rscript ${db}/script/format2lefse.R --input result/otutab.txt \
      --taxonomy result/taxonomy.txt --design result/metadata.txt \
      --group Group --threshold 0.4 \
      --output result/lefse/LEfSe

    ### 3.2 Rmd生成输入文件(可选)
    #1. result目录中存在otutab.txt, metadata.txt, taxonomy.txt三个文件；
    #2. Rstudio打开EasyAmplicon中format2lefse.Rmd，另存至result目录并Knit生成输入文件和可重复计算网页；

    ### 3.3 LEfSe分析
    #方法1. 打开LEfSe.txt并在线提交 https://www.bic.ac.cn/BIC/#/analysis?tool_type=tool&page=b%27MzY%3D%27
    #方法2. LEfSe本地分析(限Linux系统、选学)，参考代码见附录
    #方法3. LEfSe官网在线使用


# 25、QIIME 2分析流程

    # 代码详见 qiime2/pipeline_qiime2.sh


# 31、功能预测

## 1. PICRUSt功能预测

### PICRUSt 1.0

    # 方法1. 使用 https://www.bic.ac.cn/ImageGP 在线分析 gg/otutab.txt
    # 方法2. Linux服务器用户可参考"附录2. PICRUSt功能预测"实现软件安装和分析
    # 然后结果使用STAMP/R进行差异比较
    # R语言绘图
    # 输入文件格式调整
    l=L2
    sed '/# Const/d;s/OTU //' result/picrust/all_level.ko.${l}.txt > result/picrust/${l}.txt
    num=`head -n1 result/picrust/${l}.txt|wc -w`
    paste <(cut -f $num result/picrust/${l}.txt) <(cut -f 1-$[num-1] result/picrust/${l}.txt) \
      > result/picrust/${l}.spf
    cut -f 2- result/picrust/${l}.spf > result/picrust/${l}.mat.txt
    awk 'BEGIN{FS=OFS="\t"} {print $2,$1}' result/picrust/${l}.spf | sed 's/;/\t/' | sed '1 s/ID/Pathway\tCategory/' \
      > result/picrust/${l}.anno.txt
    # 差异比较
    compare="KO-WT"
    Rscript ${db}/script/compare.R \
      --input result/picrust/${l}.mat.txt --design result/metadata.txt \
      --group Group --compare ${compare} --threshold 0 \
      --method wilcox --pvalue 0.05 --fdr 0.2 \
      --output result/picrust/
    # 可对结果${compare}.txt筛选
    # 绘制指定组(A/B)的柱状图，按高分类级着色和分面
    Rscript ${db}/script/compare_hierarchy_facet.R \
      --input result/picrust/${compare}.txt \
      --data MeanA \
      --annotation result/picrust/${l}.anno.txt \
      --output result/picrust/${compare}.MeanA.bar.pdf
    # 绘制两组显著差异柱状图，按高分类级分面
    Rscript ${db}/script/compare_hierarchy_facet2.R \
      --input result/picrust/${compare}.txt \
      --pvalue 0.05 --fdr 0.1 \
      --annotation result/picrust/${l}.anno.txt \
      --output result/picrust/${compare}.bar.pdf
      
### PICRUSt 2.0
    
    # 软件安装见附录6. PICRUSt环境导出和导入

    # (可选)PICRUSt2(Linux/Windows下Linux子系统，要求>16GB内存)
    # 安装参考附录5的方式直接下载安装包并解压即可使用
    
    # Linux中加载conda环境
    conda activate picrust2
    # 进入工作目录，服务器要修改工作目录
    wd=/mnt/d/amplicon/result/picrust2
    mkdir -p ${wd} && cd ${wd}
    # 运行流程，内存15.7GB，耗时12m
    picrust2_pipeline.py -s ../otus.fa -i ../otutab.txt -o ./out -p 8
    # 添加EC/KO/Pathway注释
    cd out
    add_descriptions.py -i pathways_out/path_abun_unstrat.tsv.gz -m METACYC \
      -o METACYC.tsv
    add_descriptions.py -i EC_metagenome_out/pred_metagenome_unstrat.tsv.gz -m EC \
      -o EC.tsv
    add_descriptions.py -i KO_metagenome_out/pred_metagenome_unstrat.tsv.gz -m KO \
      -o KO.tsv
    # KEGG按层级合并
    db=/mnt/d/EasyMicrobiome/
    python3 ${db}/script/summarizeAbundance.py \
      -i KO.tsv \
	    -m ${db}/kegg/KO1-4.txt \
	    -c 2,3,4 -s ',+,+,' -n raw \
	    -o KEGG
    # 统计各层级特征数量
    wc -l KEGG*
    # 可视化见picrust2文件夹中ggpicrust2.Rmd

## 2. 元素循环FAPROTAX

    ## 方法1. 在线分析，推荐使用 http://www.bic.ac.cn/ImageGP/index.php/Home/Index/FAPROTAX.html 在线分析
    
    ## 方法2. Linux下分析、如QIIME 2环境下

    # 设置工作目录
    wd=/mnt/d/amplicon/result/faprotax/
    mkdir -p ${wd} && cd ${wd}
    # 设置脚本目录
    sd=/mnt/d/EasyMicrobiome/script/FAPROTAX_1.2.7

    ### 1. 软件安装
    # 注：软件已经下载至 EasyMicrobiome/script目录，在qiime2环境下运行可满足依赖关系
    #(可选)下载软件新版本，以1.2.7版为例， 2023/7/14更新数据库
    #wget -c https://pages.uoregon.edu/slouca/LoucaLab/archive/FAPROTAX/SECTION_Download/MODULE_Downloads/CLASS_Latest%20release/UNIT_FAPROTAX_1.2.7/FAPROTAX_1.2.7.zip
    #解压
    #unzip FAPROTAX_1.2.7.zip
    #新建一个python3环境并配置依赖关系，或进入qiime2 python3环境
    conda activate qiime2-amplicon-2024.2
    # source /home/silico_biotech/miniconda3/envs/qiime2/bin/activate
    #测试是否可运行，弹出帮助即正常工作
    python $sd/collapse_table.py

    ### 2. 制作输入OTU表
    #txt转换为biom json格式
    biom convert -i ../otutab_rare.txt -o otutab_rare.biom --table-type="OTU table" --to-json
    #添加物种注释
    biom add-metadata -i otutab_rare.biom --observation-metadata-fp ../taxonomy2.txt \
      -o otutab_rare_tax.biom --sc-separated taxonomy \
      --observation-header OTUID,taxonomy
    #指定输入文件、物种注释、输出文件、注释列名、属性列名

    ### 3. FAPROTAX功能预测
    #python运行collapse_table.py脚本、输入带有物种注释OTU表tax.biom、
    #-g指定数据库位置，物种注释列名，输出过程信息，强制覆盖结果，结果文件和细节
    #下载faprotax.txt，配合实验设计可进行统计分析
    #faprotax_report.txt查看每个类别中具体来源哪些OTUs
    python ${sd}/collapse_table.py -i otutab_rare_tax.biom \
      -g ${sd}/FAPROTAX.txt \
      --collapse_by_metadata 'taxonomy' -v --force \
      -o faprotax.txt -r faprotax_report.txt

    ### 4. 制作OTU对应功能注释有无矩阵
    # 对ASV(OTU)注释行，及前一行标题进行筛选
    grep 'ASV_' -B 1 faprotax_report.txt | grep -v -P '^--$' > faprotax_report.clean
    # faprotax_report_sum.pl脚本将数据整理为表格，位于public/scrit中
    perl ${sd}/../faprotax_report_sum.pl -i faprotax_report.clean -o faprotax_report
    # 查看功能有无矩阵，-S不换行
    less -S faprotax_report.mat

## 3. Bugbase细菌表型预测

    ### 1. Bugbase命令行分析
    cd ${wd}/result
    bugbase=${db}/script/BugBase
    rm -rf bugbase/
    # 脚本已经优化适合R4.0，biom包更新为biomformat
    Rscript ${bugbase}/bin/run.bugbase.r -L ${bugbase} \
      -i gg/otutab.txt -m metadata.txt -c Group -o bugbase/

    ### 2. 其它可用分析
    # 使用 http://www.bic.ac.cn/ImageGP/index.php/Home/Index/BugBase.html
    # 官网，https://bugbase.cs.umn.edu/ ，有报错，不推荐
    # Bugbase细菌表型预测Linux，详见附录3. Bugbase细菌表型预测


# 32、MachineLearning机器学习

    # RandomForest包使用的R代码见advanced/RandomForestClassification和RandomForestRegression
    ## Silme2随机森林/Adaboost使用代码见EasyMicrobiome/script/slime2目录中的slime2.py，详见附录4
    # 使用实战(使用QIIME 2的Python3环境，以在Windows中为例)
    conda activate qiime2-amplicon-2024.2
    cd /mnt/d/EasyMicrobiome/script/slime2
    #使用adaboost计算10000次(16.7s)，推荐千万次
    ./slime2.py otutab.txt design.txt --normalize --tag ab_e4 ab -n 10000
    #使用RandomForest计算10000次(14.5s)，推荐百万次，支持多线程
    ./slime2.py otutab.txt design.txt --normalize --tag rf_e4 rf -n 10000


# 33、Evolution进化树

    cd ${wd}
    mkdir -p result/tree
    cd ${wd}/result/tree

## 1. 筛选高丰度/指定的特征

    #方法1. 按丰度筛选特征，一般选0.001或0.005，且OTU数量在30-150个范围内
    #统计特征表中ASV数量，如总计1609个
    tail -n+2 ../otutab_rare.txt | wc -l
    #按相对丰度0.2%筛选高丰度OTU
    usearch -otutab_trim ../otutab_rare.txt \
        -min_otu_freq 0.002 \
        -output otutab.txt
    #统计筛选OTU表特征数量，总计~81个
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

## 2. 构建进化树

    # 起始文件为 result/tree目录中 otus.fa(序列)、annotation.txt(物种和相对丰度)文件
    # Muscle软件进行序列对齐，3s
    muscle -in otus.fa -out otus_aligned.fas

    ### 方法1. 利用IQ-TREE快速构建ML进化树，2m
    rm -rf iqtree
    mkdir -p iqtree
    iqtree -s otus_aligned.fas \
        -bb 1000 -redo -alrt 1000 -nt AUTO \
        -pre iqtree/otus

    ### 方法2. FastTree快速建树(Linux)
    # 注意FastTree软件输入文件为fasta格式的文件，而不是通常用的Phylip格式。输出文件是Newick格式。
    # 该方法适合于大数据，例如几百个OTUs的系统发育树！
    # Ubuntu上安装fasttree可以使用`apt install fasttree`
    # fasttree -gtr -nt otus_aligned.fas > otus.nwk

## 3. 进化树美化

    # 访问http://itol.embl.de/，上传otus.nwk，再拖拽下方生成的注释方案于树上即美化

    ## 方案1. 外圈颜色、形状分类和丰度方案
    # annotation.txt OTU对应物种注释和丰度，
    # -a 找不到输入列将终止运行（默认不执行）-c 将整数列转换为factor或具有小数点的数字，-t 偏离提示标签时转换ID列，-w 颜色带，区域宽度等， -D输出目录，-i OTU列名，-l OTU显示名称如种/属/科名，
    # cd ${wd}/result/tree
    Rscript ${db}/script/table2itol.R -a -c double -D plan1 -i OTUID -l Genus -t %s -w 0.5 annotation.txt
    # 生成注释文件中每列为单独一个文件

    ## 方案2. 生成丰度柱形图注释文件
    Rscript ${db}/script/table2itol.R -a -d -c none -D plan2 -b Phylum -i OTUID -l Genus -t %s -w 0.5 annotation.txt

    ## 方案3. 生成热图注释文件
    Rscript ${db}/script/table2itol.R -c keep -D plan3 -i OTUID -t %s otutab.txt

    ## 方案4. 将整数转化成因子生成注释文件
    Rscript ${db}/script/table2itol.R -a -c factor -D plan4 -i OTUID -l Genus -t %s -w 0 annotation.txt

    # 树iqtree/otus.contree在 http://itol.embl.de/ 上展示，拖拽不同Plan中的文件添加树注释

    # 返回工作目录
    cd ${wd}

## 4. 进化树可视化

   https://www.bic.ac.cn/BIC/#/ 提供了更简易的可视化方式

# 附加视频

    # 目录 Supp，网课有对应视频(可能编号不同，找关键字)


## S1. 网络分析R/CytoGephi

    # 目录 Supp/S1NetWork

## S2. 溯源和马尔可夫链

    # 目录 Supp/S2SourcetrackerFeastMarkov


## S11、网络分析ggClusterNet

    # 代码：advanced/ggClusterNet/Practice.Rmd

## S12、Microeco包数据可视化

    # 代码：advanced/microeco/Practice.Rmd


# 附录：Linux服务器下分析(选学)

    #注：Windows下可能无法运行以下代码，推荐在Linux，或Windows下Linux子系统下conda安装相关程序

## 1. LEfSe分析

    mkdir -p ~/amplicon/lefse
    cd ~/amplicon/lefse
    # format2lefse.Rmd代码制作或上传输入文件LEfSe.txt
    # 安装lefse
    # conda install lefse

    #格式转换为lefse内部格式
    lefse-format_input.py LEfSe.txt input.in -c 1 -o 1000000
    #运行lefse
    run_lefse.py input.in input.res
    #绘制物种树注释差异
    lefse-plot_cladogram.py input.res cladogram.pdf --format pdf
    #绘制所有差异features柱状图
    lefse-plot_res.py input.res res.pdf --format pdf
    #绘制单个features柱状图(同STAMP中barplot)
    head input.res #查看差异features列表
    lefse-plot_features.py -f one --feature_name "Bacteria.Firmicutes.Bacilli.Bacillales.Planococcaceae.Paenisporosarcina" \
       --format pdf input.in input.res Bacilli.pdf
    #批量绘制所有差异features柱状图，慎用(几百张差异结果柱状图阅读也很困难)
    mkdir -p features
    lefse-plot_features.py -f diff --archive none --format pdf \
      input.in input.res features/


## 2. PICRUSt功能预测

    #推荐使用 http://www.bic.ac.cn/BIC/#/analysis?tool_type=tool&page=b%27Mzk%3D%27 在线分析
    #有Linux服务器用户可参考以下代码搭建本地流程

    n=picrust
    conda create -n ${n} ${n} -c bioconda -y

    wd=/mnt/d/amplicon
    cd $wd/result/gg
    # 启动环境
    conda activate picrust
    #上传gg/otutab.txt至当前目录
    #转换为OTU表通用格式，方便下游分析和统计
    biom convert -i otutab.txt \
        -o otutab.biom \
        --table-type="OTU table" --to-json

    # 设置数据库目录，如 /mnt/d/db
    db=~/db
    #校正拷贝数，30s, 102M
    normalize_by_copy_number.py -i otutab.biom \
        -o otutab_norm.biom \
        -c ${db}/picrust/16S_13_5_precalculated.tab.gz
    #预测宏基因组KO表，3m,1.5G，biom方便下游归类，txt方便查看分析
    predict_metagenomes.py -i otutab_norm.biom \
        -o ko.biom \
        -c ${db}/picrust/ko_13_5_precalculated.tab.gz
    predict_metagenomes.py -f -i otutab_norm.biom \
        -o ko.txt \
        -c ${db}/picrust/ko_13_5_precalculated.tab.gz

    #按功能级别分类汇总, -c输出KEGG_Pathways，分1-3级
    sed  -i '/# Constru/d;s/#OTU //' ko.txt
    num=`head -n1 ko.txt|wc -w`
    paste <(cut -f $num ko.txt) <(cut -f 1-$[num-1] ko.txt) > ko.spf
    for i in 1 2 3;do
      categorize_by_function.py -f -i ko.biom -c KEGG_Pathways -l ${i} -o pathway${i}.txt
      sed  -i '/# Const/d;s/#OTU //' pathway${i}.txt
      paste <(cut -f $num pathway${i}.txt) <(cut -f 1-$[num-1] pathway${i}.txt) > pathway${i}.spf
    done
    wc -l *.spf


## 3. Bugbase细菌表型预测

    ### 1. 软件安装(己整合到EasyMicrobiome中，原代码需要更新才能在当前运行)
    #有两种方法可选，推荐第一种，可选第二种，仅需运行一次
    # #方法1. git下载，需要有git
    # git clone https://github.com/knights-lab/BugBase
    # #方法2. 下载并解压
    # wget -c https://github.com/knights-lab/BugBase/archive/master.zip
    # mv master.zip BugBase.zip
    # unzip BugBase.zip
    # mv BugBase-master/ BugBase

    cd BugBase
    #安装依赖包
    export BUGBASE_PATH=`pwd`
    export PATH=$PATH:`pwd`/bin
    #安装了所有依赖包
    run.bugbase.r -h
    #测试数据
    run.bugbase.r -i doc/data/HMP_s15.txt -m doc/data/HMP_map.txt -c HMPBODYSUBSITE -o output


    ### 2. 准备输入文件
    cd ~/amplicon/result
    #输入文件：基于greengene OTU表的biom格式(本地分析支持txt格式无需转换)和mapping file(metadata.txt首行添加#)
    #上传实验设计+刚才生成的otutab_gg.txt
    #生成在线分析使用的biom1.0格式
    biom convert -i gg/otutab.txt -o otutab_gg.biom --table-type="OTU table" --to-json
    sed '1 s/^/#/' metadata.txt > MappingFile.txt
    #下载otutab_gg.biom 和 MappingFile.txt用于在线分析

    ### 3. 本地分析
    export BUGBASE_PATH=`pwd`
    export PATH=$PATH:`pwd`/bin
    run.bugbase.r -i otutab_gg.txt -m MappingFile.txt -c Group -o phenotype/

## 4. Silme2随机森林/Adaboost

    #下载安装
    # cd ~/software/
    # wget https://github.com/swo/slime2/archive/master.zip
    # mv master.zip slime2.zip
    # unzip slime2.zip
    # mv slime2-master/ slime2
    # cp slime2/slime2.py ~/bin/
    # chmod +x ~/bin/slime2.py
    #安装依赖包
    # sudo pip3 install --upgrade pip
    # sudo pip3 install pandas
    # sudo pip3 install sklearn


## 5. PICRUSt2环境导出和导入
    
    # 方法1. 下载安装包并解压 
    # 下载安装包，备用链接见百度云：https://pan.baidu.com/s/1Ikd_47HHODOqC3Rcx6eJ6Q?pwd=0315
    wget -c ftp://download.nmdc.cn/tools/conda/picrust2.tar.gz
    # 指定安装目录并解压
    mkdir -p ~/miniconda3/envs/picrust2
    tar -xvzf picrust2.tar.gz -C ~/miniconda3/envs/picrust2
    # 激活环境并初始化
    conda activate picrust2
    conda unpack

    # 方法2. 直接安装或打包安装环境
    n=picrust2
    conda create -n ${n} -c bioconda -c conda-forge ${n}=2.3.0_b
    # 加载环境
    conda activate ${n}
    # 打包环境(可选)
    conda pack -n ${n} -o ${n}.tar.gz


# 常见问题

## 1. 文件phred质量错误——Fastq质量值64转33

    # 使用head查看fastq文件，phred64质量值多为小写字母，需要使用vsearch的--fastq_convert命令转换为通用的phred33格式。

    cd /c/amplicon/FAQ/01Q64Q33
    # 预览phred64格式，注意看第4行质量值多为小写字母
    head -n4 test_64.fq
    # 转换质量值64编码格式为33
    vsearch --fastq_convert test_64.fq \
        --fastq_ascii 64 --fastq_asciiout 33 \
        --fastqout test.fq 
    # 查看转换后33编码格式，质量值多为大写字母
    head -n4 test.fq

    # 如果是Ion torrent测序结果，由于是非主流测序平台，需要公司转换帮助转换为标准的Phred33格式文件才可以使用。

## 2. 序列双端已经合并——单端序列添加样本名

    # 扩增子分析要求序列名为样品名+序列编号，双端序列在合并同时可直接添加样本名。单端序列，或双端合并的序列需单独添加。这里使用vsearch的--fastq_convert命令中的--relabel参加添加样本名

    cd /c/amplicon/FAQ/02relabel
    # 查看文件序列名
    head -n1 test.fq
    # 序列按样本重命名
    vsearch --fastq_convert test.fq \
        --relabel WT1. \
        --fastqout WT1.fq
    # 查看重命名结果
    head -n1 WT1.fq

## 3. 数据过大无法使用usearch聚类或去噪,替换vsearch

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

## 4. 读长计数(Read counts)标准化为相对丰度

    cd /c/amplicon/FAQ/04norm
    # 求取各个OTU在样品中的丰度频率(标准化为总和1)
    usearch -otutab_counts2freqs otutab.txt \
        -output otutab_freq.txt

## 5. 运行R提示Permission denied
 
    # 例如write.table保存表时，报错信息示例如下：意思是写入文件无权限，一般为目标文件正在被打开，请关闭相关文件后重试

    Error in file(file, ifelse(append, "a", "w")) :
    Calls: write.table -> file
    : Warning message:
    In file(file, ifelse(append, "a", "w")) :
      'result/raw/otutab_nonBac.txt': Permission denied

## 6. 文件批量命名

    # 如我们有文件A1和A2，编写一个样本名对应目标名的表格metadata.txt，检查样本名是否唯一，使用awk进行批量改名

    cd /c/amplicon/FAQ/06rename
    # (可选)快速生成文件列表，用于编辑metadata.txt，如A1.fq修改为WT1.fastq，以此类推，参考metadata.bak.txt
    ls *.fq > metadata.txt
    # 编辑列表，第二名为最终命名，确定名称唯一
    # 转换行尾换行符
    sed -i 's/\r//' metadata.txt
    # 检查手动命名列2是否唯一
    cut -f 2 metadata.txt|wc -l
    cut -f 2 metadata.txt|sort|uniq|wc -l
    # 如果两次结果一致，则命名非冗余
    # 可选移动mv，复制cp，硬链ln，或软链ln -s
    # 此处使用复制cp
    awk '{system("cp "$1" "$2)}' metadata.txt

## 7. Rstudio中Terminal找不到Linux命令

    # 需要把 C:\Program Files\Git\usr\bin 目录添加到系统环境变量
    # 文件资源管理器——此电脑——属性——高级系统设置——环境变量——系统变量——Path——编辑——新建——填写“C:\Program Files\Git\usr\bin”——确定——确定——确定
    # 注意win10系统是一个目录一行；win7中多个目录用分号分隔，注意向后添加目录

## 8. usearch -alpha_div_rare结果前两行出现“-”

    #问题：抽样0时补“-”，且缺失制表符

    #处理：替换“-”为"制作符\t+0"即可恢复

    cd /c/amplicon/FAQ/08rare
    sed "s/-/\t0.0/g" alpha_rare_wrong.txt\
        > alpha_rare.txt

## 9. 物种注释otus.sintax方向全为“-”，需要序列取反向互补

    #是原始序列方向错误，将filtered.fa序列需要取反向互补。再从头开始分析

    cd /c/amplicon/FAQ/09revcom
    vsearch --fastx_revcomp filtered_RC.fa \
      --fastaout filtered.fa

## 10. windows换行符查看和删除

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

## 11. UNITE数据库分析报错

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

## 12. Windows的Linux子系统本地安装qiime2

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

## 13. RDP 16-18注释结果比较

    # 统计序列中门的数量，从60降为39
    grep '>' ${db}/usearch/rdp_16s_v16_sp.fa|cut -f2 -d ';'|cut -f1-2 -d ','|sort|uniq|wc -l
    grep '>' ${db}/usearch/rdp_16s_v18.fa|cut -f2 -d ';'|cut -f1-2 -d ','|sort|uniq|wc -l
    # 统计序列中属的数量，从2517增长为3061
    grep '>' ${db}/usearch/rdp_16s_v16_sp.fa|cut -f2 -d ';'|cut -f1-6 -d ','|sort|uniq|wc -l
    grep '>' ${db}/usearch/rdp_16s_v18.fa|cut -f2 -d ';'|cut -f1-6 -d ','|sort|uniq|wc -l

    cd /c/amplicon/FAQ/13rdp16_18
    # 门由15个降为13个
    tail -n+2 rdp16_sintax.txt|cut -f3|sort|uniq -c|wc -l
    tail -n+2 rdp18_sintax.txt|cut -f3|sort|uniq -c|wc -l
    # 属由176个降为144个
    tail -n+2 rdp16_sintax.txt|cut -f7|sort|uniq -c|wc -l
    tail -n+2 rdp18_sintax.txt|cut -f7|sort|uniq -c|wc -l  


## 14. usearch生成OTU表小样本比vsearch更快

    # usearch生成特征表，小样本(<30)快；但大样本受限且多线程效率低，83.2%, 4核17s
    time usearch -otutab temp/filtered.fa \
      -otus result/raw/otus.fa \
      -threads 4 \
      -otutabout result/raw/otutab.txt

    # vsearch比对，更准更慢，但并行24-96线程更强
    vsearch --usearch_global temp/filtered.fa --db ${db}/gg/97_otus.fa \
      --otutabout result/gg/otutab.txt --id 0.97 --threads 12
    # 比对率81.04%, 1核30m, 12核7m

## 14. vsearch不可用/无结果
    # Windows系统更新导致vserach缺失dll依赖关系
    # 下载新版vsearch，https://github.com/torognes/vsearch 右侧release找最新下载链接，如 vsearch-2.30.0-win-x86_64.zip
    # 软件vsearch.exe配套带有libbz2.dll和zlib1.dll，替换EasyMicrobiome中win目录即可

# 版本更新记录

- 2021/4/3 EasyAmplicon 1.11:
    - R包amplicon升级为 1.11.0，解决metadata两列报错的问题。
    - 调整课程顺序，每天上午9点-12点2节，下午1点半-6点3节。
    - 提供附加课程Supp目录。
- 2021/7/23 EasyAmplicon 1.12:
    - R运行环境升级为4.1.0，配套有4.1.zip的全套包
    - R包amplicon升级为 1.12.0，alpha_boxplot去掉Y轴的index
    - alpha_boxplot.R增加标准化、转置等参数，可用于绘制任何特征箱线图
    - beta_pcoa/cpcoa.R增加控制椭圆、标签显示等参数
    - tax_stackplot.R增加多种配色方案
    - picurst流程更新，并提供打包的conda下载
    - picurst2新增KO合并为KEGG通路1-3级代码，并提供打包的conda下载
    - 随机森林：提供分类级筛选、随机数筛选、可视化代码
- 2021/10/15 EasyAmplicon 1.13:
    - R运行环境升级为4.1.1，配套有4.1.zip的最新全套包
    - 元数据方差分解PERMANOVA：在Diversity-tutorial.Rmd中Beta多样性分析中新增adonis计算变量对群落的解析率和显著性分析
    - 树图treemap升级后无颜色，改为代码供参考，并在Diversity_tutrial.Rmd中删除此部分
    - alpha_boxplot输出无默认目录，可指定文件名头，添加无ID报错注释
- 2022/1/7 EasyAmplicon 1.14:
    - R运行环境升级为4.1.2，配套有4.1.zip的最新全套包
    - RStudio更新为2021.09.1
    - 文涛重写amplicon包中tax_maptree函数，不依赖其他包，解决无法着色问题
    - EasyMicrobiome中添加compare_stamp.R脚本，直接差异比较绘制STAMP扩展柱状图；代码详见result/CompareStamp.Rmd
    - EasyMicrobiome中添加compare_hierarchy_facet.R和compare_hierarchy_facet2.R，展示KEGG的1，2级总览和差异
    - 更新高级分析目录advanced：包括环境因子、马尔可无链、网络模块、网络比较、随机森林分类、随机森林回归、微生态等
- 2023/2/3 EasyAmplicon 1.18:
    - R运行环境升级为4.2.2，配套有4.2.zip的最新全套包
    - RStudio更新为2022.12.0
    - amplicon、EasyAmplicon和EasyMicrobiome更新为1.18
    - QIIME 2更新为v2023.2
    - vsearch更新为v2.22.1
    - 新增ggClusterNet课程-文涛
- 2023/10/13 EasyAmplicon 1.20:
    - R运行环境升级为4.3.1，配套有4.3.zip的最新全套包
    - RStudio更新为2023.12.0
    - amplicon、EasyAmplicon和EasyMicrobiome更新为1.20
    - QIIME 2更新为v2023.7，数据库更新为greengene2 2022.10
    - 新增ggpicrust2分析picrust2结果可视化
    - 更新FAPROTAX为1.2.7
- 2025/4/11 EasyAmplicon 1.23:
    - R运行环境升级为4.4.3，配套有4.4.zip的最新全套包
    - RStudio更新为2024.12.1
    - amplicon、EasyAmplicon和EasyMicrobiome更新为1.23
    - QIIME 2更新为v2024.10，数据库更新为python 3.10版学习的greengene2 2022.10
    - 更新vsearch至2.30，新增2个dll文件，解决windows中运行无结果问题

- 2025/10/16 EasyAmplicon 2.00:
    - EasyAmplicon 2.0文章在Advanced Science发表
    - PacBio三代测序扩增子流程；
    - Nanopore三代测序扩增子流程；
    - Nanopore三代测序扩增子流程；
    - vsearch v2.15.2更新为v2.30.1，新增libbz2.dll和zlib1.dll，支持压缩文件处理
    - vsearch v2.15.2 update to v2.30.1, add libbz2.dll和zlib1.dll and support zip and gz file


- ITS测序使用UNITE分析流程；
- 18S测序使用SILVA分析流程，含SILVA格式；
- ggClusterNet 2网络和多界分析；
- QIIME 2中greengene2的最新数据格式化；
- STAMP中脚本缺少图例，还有KO-WT顺序调整；

每季度视频课程安排：http://www.ehbio.com/trainLongTerm/TrainLongTerm/amplicongenomeLearnGuide.html

使用此脚本，请引用下文：

If used this script, please cited:

Hao Luo, et al. 2025. EasyAmplicon 2: Expanding PacBio and Nanopore Long Amplicon Sequencing Analysis Pipeline for Microbiome. Advanced Science 12: https://doi.org/https://doi.org/10.1002/advs.202512447

Yousuf, et al. 2024. Unveiling microbial communities with EasyAmplicon: 
A user-centric guide to perform amplicon sequencing data analysis. 
iMetaOmics 1: e42. https://doi.org/10.1002/imo2.42

Yong-Xin Liu, Lei Chen, Tengfei Ma, et al. 2023. 
EasyAmplicon: An easy-to-use, open-source, reproducible, and community-based pipeline for amplicon data analysis in microbiome research. 
iMeta 2: e83. https://doi.org/10.1002/imt2.83

Copyright 2016-2025 Yong-Xin Liu <liuyongxin@caas.cn>, Tao Wen <taowen@njau.edu.cn>, Tong Chen <chent@nrc.ac.cn>