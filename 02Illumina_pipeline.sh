[TOC]

# 易扩增子EasyAmplicon

    # 作者 Authors: 刘永鑫(Yong-Xin Liu), 陈同(Tong Chen)等
    # 版本 Version: v1.20
    # 更新 Update: 2023-10-13
    # 系统要求 System requirement: Windows 10+ / Mac OS 10.12+ / Ubuntu 20.04+
    # 引文 Reference: Liu, et al. 2023. EasyAmplicon: An easy-to-use, open-source, reproducible, and community-based
    # pipeline for amplicon data analysis in microbiome research. iMeta 2: e83. https://doi.org/10.1002/imt2.83


    # 设置工作(work directory, wd)和软件数据库(database, db)目录
    # 添加环境变量，并进入工作目录 Add environmental variables and enter work directory
    # **每次打开Rstudio必须运行下面4行 Run it**，可选替换${db}为EasyMicrobiome安装位置
    wd=/c/amplicon
    db=/c/EasyMicrobiome
    PATH=$PATH:${db}/win
    cd ${wd}


## 1. 起始文件 start files

    # 1. 分析流程pipeline.sh
    # 2. 样本元信息metadata.txt，保存于result目录
    # 3. 测序数据fastq文件保存于seq目录，通常以`.fastq.gz`结尾，每个样品一对文件
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

    # 公司返回的测序结果，通常为一个样品一对fq/fastq.gz格式压缩文件
    # 文件名与样品名务必对应：不一致时手工修改，批量改名见"常见问题6"
    # 如果测序数据是.gz的压缩文件，有时需要使用gunzip解压后使用，vsearch通常可直接读取压缩文件
    # gunzip seq/*.gz
    # 批量统计测序数据并汇总表
    seqkit stat seq/*.fastq.gz > result/seqkit.txt
    head result/seqkit.txt

### 1.3. 流程和数据库 pipeline & database

    # 数据库第一次使用必须解压，以后可跳过此段

    # usearchs可用16S/18S/ITS数据库：RDP, SILVA和UNITE，本地文件位置 ${db}/usearch/
    # usearch数据库database下载页: http://www.drive5.com/usearch/manual/sintax_downloads.html
    # 解压16S RDP数据库，gunzip解压缩，seqkit stat统计
    # 保留原始压缩文件
    gunzip -c ${db}/usearch/rdp_16s_v18.fa.gz > ${db}/usearch/rdp_16s_v18.fa
    seqkit stat ${db}/usearch/rdp_16s_v18.fa # 2.1万条序列
    # 解压ITS UNITE数据库，需自行从官网或网盘db/amplicon/usearch中下载
    # gunzip -c ${db}/usearch/utax_reference_dataset_all_25.07.2023.fasta.gz >${db}/usearch/unite.fa
    seqkit stat ${db}/usearch/unite.fa # 32.6万
    # Greengene数据库用于功能注释: ftp://greengenes.microbio.me/greengenes_release/gg_13_5/gg_13_8_otus.tar.gz
    # 默认解压会删除原文件，-c指定输出至屏幕，> 写入新文件(可改名)
    gunzip -c ${db}/gg/97_otus.fasta.gz > ${db}/gg/97_otus.fa
    seqkit stat ${db}/gg/97_otus.fa


## 2. 序列合并和重命名 reads merge and rename

### 2.1 合并双端序列并按样品重命名 Merge pair-end reads and rename

    #依照实验设计批处理并合并
    #tail -n+2去表头，cut -f1取第一列，获得样本列表；18个样本x1.5万对序列合并8s
    #Win下复制Ctrl+C为Linux下中止，为防止异常中断，结尾添加&转后台，无显示后按回车继续
    
    # 一部分电脑 rush 不支持，运行时调度失败，请使用 for 循环部分
    # for 循环部分是放入后台运行的，点完 run 之后，看上去程序已运行完，实际没运行完，而是正在运行中。
    # 不要急于运行后面的程序。
    # 之前课程，有发现每次运行结果都不一样，就是因为 for 循环部分没运行完，只生成了部分数据，导致后面
    # 每个样品 reads 数每次运行都会不一致。
    #方法1.for循环顺序处理
    # time for i in `tail -n+2 result/metadata.txt|cut -f1`;do
    #   vsearch --fastq_mergepairs seq/${i}_1.fastq.gz --reverse seq/${i}_2.fastq.gz \
    #   --fastqout temp/${i}.merged.fq --relabel ${i}.
    # done &

    # 一部分电脑 rush 不支持，运行时调度失败，请使用 for 循环部分
    #方法2.rush并行处理，任务数jobs(j),2可加速1倍4s；建议设置2-4
    time tail -n+2 result/metadata.txt | cut -f 1 | \
     rush -j 2 "vsearch --fastq_mergepairs seq/{}_1.fastq.gz --reverse seq/{}_2.fastq.gz \
      --fastqout temp/{}.merged.fq --relabel {}."
    # 检查最后一个文件前10行中样本名
    head temp/`tail -n+2 result/metadata.txt | cut -f 1 | tail -n1`.merged.fq | grep ^@
    
    ##方法3.不支持压缩文件时解压再双端合并
    #  time tail -n+2 result/metadata.txt | cut -f 1 | \
    #    rush -j 1 "vsearch --fastq_mergepairs <(zcat seq/{}_1.fastq.gz) --reverse <(zcat seq/{}_2.fastq.gz) \
    #     --fastqout temp/{}.merged.fq --relabel {}."
    # 
    #   time for i in `tail -n+2 result/metadata.txt|cut -f1`;do
    #      vsearch --fastq_mergepairs <(zcat seq/${i}_1.fastq.gz) --reverse <(zcat seq/${i}_2.fastq.gz) \
    #      --fastqout temp/${i}.merged.fq --relabel ${i}.
    #    done &
      
### 2.2 (可选)单端文件改名 Single-end reads rename

    # # 单个序列改名示例
    # i=WT1
    # gunzip -c seq/${i}_1.fastq.gz > seq/${i}.fq
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

    # 前向引物18bp，反向引物20bp
    # 务必清楚实验设计和引物长度，引物已经去除可填0
    time vsearch --fastx_filter temp/all.fq \
      --fastq_stripleft 18 --fastq_stripright 20 \
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
    # Uni_1;size=6050  - 去冗余后序列的名字 Uni_1；该序列在所有样品测序数据中出现 6050 次
    # 为出现最多的序列。
    head -n 2 temp/uniques.fa

### 4.2 聚类OTU/去噪ASV Cluster or denoise

    #有两种方法：推荐unoise3去噪获得单碱基精度ASV，备选传统的97%聚类OTU(属水平精度)
    #usearch两种特征挑选方法均自带de novo去嵌合体
    #-minsize二次过滤，控制OTU/ASV数量至1-5千，方便下游统计分析

    #方法1. 97%聚类OTU，适合大数据/ASV规律不明显/reviewer要求
    # usearch -cluster_otus temp/uniques.fa -minsize 10 \
    #  -otus temp/otus.fa \
    #  -relabel OTU_

    #方法2. ASV去噪 Denoise: predict biological sequences and filter chimeras
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
    # 可自行下载silva并解压(替换rdp_16s_v18.fa为silva_16s_v123.fa)，极慢但理论上更好
    vsearch --uchime_ref temp/otus.fa \
      -db ${db}/usearch/rdp_16s_v18.fa \
      --nonchimeras result/raw/otus.fa
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
    # id(0.97)：97%相似度比对(更高数据使用率，更快)
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
    awk 'BEGIN{OFS=FS="\t"}{if(FNR==1) {for(i=3;i<=NF;i++) a[i]=$i; print "OTU","Group";} \
        else {for(i=3;i<=NF;i++) if($i>0.1) print $1, a[i];}}' \
        result/otutab_mean.txt > result/alpha/otu_group_exist.txt
    head result/alpha/otu_group_exist.txt
    cut -f 2 result/alpha/otu_group_exist.txt | sort | uniq -c
    # 试一试：不同丰度下各组有多少OTU/ASV
    # 可在 http://ehbio.com/test/venn/ 中绘图并显示各组共有和特有维恩或网络图
    # 也可在 http://www.ehbio.com/ImageGP 绘制Venn、upSetView和Sanky

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
    vsearch --usearch_global temp/filtered.fa --db ${db}/gg/97_otus.fa \
      --otutabout result/gg/otutab.txt --id 0.97 --threads 8
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
    md5sum *_1.fastq.gz > md5sum1.txt
    md5sum *_2.fastq.gz > md5sum2.txt
    paste md5sum1.txt md5sum2.txt | awk '{print $2"\t"$1"\t"$4"\t"$3}' | sed 's/*//g' > ../result/md5sum.txt
    rm md5sum*
    cd ..
    cat result/md5sum.txt
