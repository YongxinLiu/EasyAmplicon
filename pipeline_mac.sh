[TOC]

#!/bin/bash

    # 作者 Authors: Yong-Xin Liu, Tong Chen, Liang Chen, Xin Zhou
    # 版本 version: v1.80
    # 更新 Update: 2020-5-28

    # 设置软件(software/binary, bin)、数据库(database, db)和工作目录(work directory, wd)
    # 务必根据自己的情况修改以下路径
    bin=~/amplicon/02Software
    db=~/amplicon/03Database
    # Linux服务器用户仅在家目录下有权限，目录如：~/amplicon，名称随意但尽量同样本名一样规范
    wd=~/amplicon
    # 用户切换至工作目录
    cd ${wd}
    # bin=02Software
    # db=03Database
    # wd="."

# 22、扩增子分析流程 16S Amplicon pipeline (De novo + Reference)

    # 系统要求 System: Windows 10 / Linux Ubuntu 18.04 / Mac OS 10.12+
    # 依赖软件 Sofware，添加 C:\amplicon\02Software\win 至环境变量
    # https://github.com/torognes/vsearch/releases
    # vsearch v2.14.1
    # https://www.drive5.com/usearch/download.html
    # usearch v10.0.240
    # R 3.6.x
    # Rstudio 1.2.5001
    # gitforwidnows 2.23.0 (Windows only)

    # 运行前准备
    # 1. 将整个amplicon目录复制到Windows C盘根目录(C:/) 或 Mac/Linux服务器家目录(~/)
    # 2. 学员按U盘`01PPT/11扩增子软件安装和测试手册.pdf`课件说明安装软件并添加环境变量
    # 3. 本节至少包括测序数据seq/*.fq.gz、样本元数据result/metadata.txt和流程脚本pipeline.sh
    # 3. Rstudio打开pipeline.sh文件，设置默认目录为C:/amplicon 或 Terminal切换至工作目录
    # Linux服务器用户Chrome访问IP地址:8787登陆Rstudio网页版，选中代码按Ctrl+Shift+C切换注释


## 1. 了解工作目录和起始文件

    #1. 原始测序数据保存于seq目录，通常以.fq.gz结尾，每个样品一对文件
    mkdir -p seq
    ls -lsh seq
    #2. 样本信息，即实验设计 metadata.txt，保存在最终结果result目录
    mkdir -p result
    #3. 分析流程pipeline.sh，每个项目复制一份，再进行个性化修改
    #4. 创建临时文件存储目录，分析结束可删除
    mkdir -p temp

### 1.1. metadata.txt 实验设计文件

    #查看前3行
    cat -tev result/metadata.txt | head -n3
    #windows用户如果结尾有^M，运行sed命令去除，并cat -A检查结果
    # sed -i 's/\r/\n/' metadata.txt
    # cat -A metadata.txt | head -n3

### 1.2. seq/*.fq.gz 原始测序数据

    # 公司返回的测序结果，通常为一个样品一对fq/fastq.gz格式压缩文件
    # 文件名格式不一致，参见常见问题6，进行简单批量改名；大量复杂批量改名推荐 https://github.com/shenwei356/brename
    # zless按页查看压缩文件，鼠标切换至下方代码区，空格翻页、q退出，注意输入法保持英文状态
    # zless seq/KO1_1.fq.gz
    # 如果测序数据是.gz结尾的压缩文件，使用gunzip解压，用完后再压缩节省空间
    # 5S
    time gunzip seq/*.gz
    # less按页查看，空格翻页、q退出；head查看前10行，-n指定行
    ls seq/
    head seq/KO1_1.fq
    # cut -c 1-60 seq/KO1_1.fq | head

### 1.3. pipeline.sh 流程依赖数据库

    # 数据库第一次使用必须解压，解压过可跳过此段

    # usearch使用16S数据库： RDP, SILVA and UNITE，本地文件位置 ${db}/usearch/
    # usearch数据库database下载页: http://www.drive5.com/sintax
    # 解压rdp用于物种注释和silva用于去嵌合体，*代表任意字符，选择以fa.gz结尾的文件
    time gunzip -f -k ${db}/usearch/*.fa.gz # 14s
    # 查看注释数据库的格式
    head -n2 ${db}/usearch/rdp_16s_v16_sp.fa

    # QIIME greengene 13_8有参数据库用于功能注释: ftp://greengenes.microbio.me/greengenes_release/gg_13_5/gg_13_8_otus.tar.gz
    gunzip -f -k ${db}/gg/*.fasta.gz


## 2. 合并双端序列并按样品重命名 Merge paired reads and label samples

    # 以WT1单样品合并为例，只为测试，以后路过
    # time vsearch --fastq_mergepairs seq/WT1_1.fq \
    #   --reverse seq/WT1_2.fq \
    #   --fastqout temp/WT1.merged.fq \
    #   --relabel WT1.

    #依照实验设计批处理并合并
    #RStudio中运行长时间异常中断，因为系统复制Ctrl+C为Linux下中止命令，由其它软件如词典取词引起
    #tail -n+2去表头，cut -f 1取第一列，即获得样本列表；18个5万对序列样本合并27s
    
    # tail -n+2 result/metadata.txt
    # tail -n+2 result/metadata.txt | cut -f 1
    # for i in `tail -n+2 result/metadata.txt | cut -f 1`;do echo ${i}; done
    
    time for i in `tail -n+2 result/metadata.txt | cut -f 1`;do
      vsearch --fastq_mergepairs seq/${i}_1.fq --reverse seq/${i}_2.fq \
      --fastqout temp/${i}.merged.fq --relabel ${i}.
    done &
    
    # & 代表转后台，防止程序被中断，Terminal右上角有小红点，代表正在运行
    # 运行结束后，注意按一下会车键
    
    # 己经合并样本直接改名，接入分析流程；替换上面vsearch --fastq_mergepairs命令为
    # usearch -fastx_relabel seq/${i}.fq -fastqout temp/${i}.merged.fq -prefix ${i}.
    # 另一种方法参考“常见问题”——“2. 序列双端已经合并——单端序列重命名/添加样本名”

    #合并所有样品至同一文件
    cat temp/*.merged.fq > temp/all.fq
    #查看文件大小 663 Mb，不同软件版本结果略有差异。
    ls -lsh temp/all.fq
    # 评估是否合理：合理继续，不合理回头检查原因
    head -n 4 temp/all.fq

## 3. 切除引物与质控 Cut primers and quality filter

    # Cut barcode 10bp + V5 19bp in left and V7 18bp in right
    # 左端10bp Barcode+19bp上游引物为29，右端为18bp下游引物填18
    # 务必清楚实验设计和引物长度，已经去除可填0。
    # 76万条序列37s
    time vsearch --fastx_filter temp/all.fq \
      --fastq_stripleft 29 --fastq_stripright 18 \
      --fastq_maxee_rate 0.01 \
      --fastaout temp/filtered.fa

    # 查看文件前2行，了解fa文件格式
    head -n 2 temp/filtered.fa
    # 筛选序列名，并显示前4行，检查序列名是否合格
    # 序列名字需符合 样品名.编号 格式，如sampl1.1 sample2.1
    # 如果名字格式不符合，获得的OTU表列数会特别多
    grep '>' temp/filtered.fa | head -n4


## 4. 去冗余和生成ASV/OTU Dereplication and cluster/denoise

### 4.1 序列去冗余 Dereplication和去低丰度序列

    # 去冗余 Find unique read sequences and abundances
    # 并添加miniuniqusize最小为10或1/1M，去除低丰度，增加计算速度
    # -sizeout输出丰度
    # --relabel必须加序列前缀，否则文件头不正常
    time vsearch --derep_fulllength temp/filtered.fa \
      --output temp/uniques.fa --relabel Uni --minuniquesize 10 --sizeout
    #2s，高丰度非冗余序列非常小(<2Mb)
    ls -lsh temp/uniques.fa
    head -n 2 temp/uniques.fa
    
    # >Uni1;size=17963
    # GTAGTCCACGCCGTAAACGGTGGGCGCTAGATGTGGGGACCTTCCACGGTTTCTGCGTCGCAGCTAACGCATTAAGCGCC

### 4.2 生成OTU/ASV Cluster OTUs / denoise ASV

    #有两种方法：推荐unoise3去噪(种/株水平精度ASV)，传统的97%聚类OTU (属水平精度)，供备选。

    #usearch自带de novo去嵌合体

    #方法1. 97%聚类OTU，不推荐，除非reviewer要求/数据量太大/ASV规律不明显
    # usearch -cluster_otus temp/uniques.fa \
    #  -otus temp/otus.fa \
    #  -relabel OTU_
    #2s, 14Mb, 878 OTUs, 320 chimeras

    #方法2. ASV去噪 Denoise: predict biological sequences and filter chimeras
    usearch -unoise3 temp/uniques.fa \
      -zotus temp/zotus.fa
    #19s, 47Mb, 2920 good, 227 chimeras
    #修改序列名：格式调整Zotu为通用的ASV，方便识别；文件名叫otus与流程兼容
    sed 's/Zotu/ASV_/g' temp/zotus.fa > temp/otus.fa
    head -n 2 temp/otus.fa
    
    # >ASV_1
    # GTAGTCCACGCCGTAAACGGTGGGCGCTAGATGTGGGGACCTTCCACGGTTTCTGCGTCGCAGCTAACGCATTAAGCGCC
    
    #方法3. 数据过大无法使用usearch时，备选vsearch见"常见问题3"


### 4.3 基于参考去嵌合 Reference-based chimera detect

    # 不推荐，容易引起假阴性，因为参考数据库无丰度信息，
    # 而de novo时要求亲本丰度为嵌合体16倍以上防止假阴性
    # 数据库选择越大越合理，减少假阴性，
    # usearch免费版受限只能用rdp trainning set, vsearch可用silva更合理
    
    mkdir -p result/raw

    # 方法1. vsearch+rdp去嵌合(快但容易假阴性)，或 
    # silva去嵌合(silva_16s_v123.fa)，推荐(慢，耗时15m ~ 3h，但更好)
    time vsearch --uchime_ref temp/otus.fa \
      -db ${db}/usearch/rdp_16s_v16_sp.fa \
      --nonchimeras result/raw/otus.fa
    # RDP: 9s, 250 (8.6%) chimeras; SILVA：10m, 255 (8.7%) chimeras

    # Win用户注意：vsearch去嵌合后每行添加了windows换行符^M，需删除
    # mac不要执行此命令
    # sed -i 's/\r//g' result/raw/otus.fa

    # 方法2. 不去嵌合
    # cp -f temp/otus.fa result/raw/otus.fa



## 5. 特征表生成和常用操作 Feature table

    # 注OTU和ASV统称为特征(Feature)，它们的区别是：
    # OTU通常按97%聚类后挑选最高丰度或中心的代表性序列；
    # ASV是基于序列进行去噪(排除或校正错误序列，并挑选丰度较高的可信序列)作为代表性序列

### 5.1 生成特征表 Creat Feature table

    # 方法1. vsearch生成特征表
    time vsearch --usearch_global temp/filtered.fa --db result/raw/otus.fa \
    	--otutabout result/raw/otutab.txt --id 0.97 --threads 4
    #652036 of 761432 (85.63%)可比对，耗时3m4s
    
    # windows用户删除换行符^M
    # sed -i 's/\r//' result/raw/otutab.txt
    # head result/raw/otutab.txt |cat -A
    # 如果是mac，运行下面的
    head result/raw/otutab.txt | cat -etv

    # 方法2. usearch生成特征表，小样本(<30)快；但大样本(>100)多线程效率低于vsearch
    # usearch -otutab temp/filtered.fa -otus result/raw/otus.fa \
    #   -otutabout result/raw/otutab.txt -threads 4
    # 84.1%, 4核1m46s；12核36s

### 5.2 物种注释-去除质体和非细菌/古菌并统计比例(可选) Remove plastid and non-Bacteria

    # RDP物种注释(rdp_16s_v16_sp)更快，但缺少完整真核来源数据,可能不完整，耗时3s;
    # SILVA数据库(silva_16s_v123.fa)更好注释真核、线粒体和叶绿体等序列，数据库是RDP的100多倍，
    # 12线程仍需73m，普通电脑3-5小时
    time vsearch --sintax result/raw/otus.fa --db ${db}/usearch/rdp_16s_v16_sp.fa \
      --tabbedout result/raw/otus.sintax --sintax_cutoff 0.6

    # 原始特征表行数
    wc -l result/raw/otutab.txt
    #R脚本选择细菌古菌(真核)、去除叶绿体、线粒体并统计比例；输出筛选并排序的OTU表
    #输入为OTU表result/raw/otutab.txt和物种注释result/raw/otus.sintax
    #输出筛选并排序的特征表result/otutab.txt和统计污染比例文件result/raw/otutab_nonBac.txt
    #真菌ITS数据，请改用otutab_filter_nonFungi.R脚本
    Rscript ${bin}/script/otutab_filter_nonBac.R
    # 筛选后特征表行数
    wc -l result/otutab.txt

    #过滤特征表对应序列
    cut -f 1 result/otutab.txt | tail -n+2 > temp/otutab.id
    usearch -fastx_getseqs result/raw/otus.fa \
        -labels temp/otutab.id \
        -fastaout result/otus.fa
    #过滤特征表对应序列注释
    awk 'NR==FNR{a[$1]=$0}NR>FNR{print a[$1]}'\
        result/raw/otus.sintax temp/otutab.id \
        > result/otus.sintax
    #补齐末尾列
    sed -i 's/\t$/\td:Unassigned/' result/otus.sintax
    head -n2 result/otus.sintax

    # 方法2. 觉得筛选不对可以不筛选
    # cp result/raw/otu* result/

    #可选统计方法：OTU表简单统计 Summary OTUs table
    usearch -otutab_stats result/otutab.txt \
      -output result/otutab.stat
    cat result/otutab.stat
    #注意最小值、1/4分位数，用于标准化

### 5.3 等量抽样标准化 normlize by subsample

    #使用vegan包进行等量重抽样，输入reads count格式Feature表result/otutab.txt
    #可指定输入文件、抽样量和随机数，输出抽平表result/otutab_rare.txt和多样性alpha/vegan.txt
    mkdir -p result/alpha
    Rscript ${bin}/script/otutab_rare.R -h
    Rscript ${bin}/script/otutab_rare.R --input result/otutab.txt \
      --depth 30000 --seed 1 \
      --normalize result/otutab_rare.txt \
      --output result/alpha/vegan.txt


## 6. Alpha多样性 Alpha diversity

### 6.1. 计算多样性指数 Calculate alpha diversity index
    #Calculate all alpha diversity index,
    #details in http://www.drive5.com/usearch/manual/alpha_metrics.html
    usearch -alpha_div result/otutab_rare.txt \
      -output result/alpha/alpha.txt

### 6.2. 计算稀释过程的丰富度变化 Rarefaction
    #稀释曲线：取1%-100%的序列中OTUs数量 Rarefaction from 1%, 2% .. 100% in richness (observed OTUs)-method fast / with_replacement / without_replacement https://drive5.com/usearch/manual/cmd_otutab_subsample.html
    time usearch -alpha_div_rare result/otutab_rare.txt \
      -output result/alpha/alpha_rare.txt -method without_replacement

### 6.3. 筛选各组高丰度菌用于比较

    #按组求均值，需根据实验设计metadata.txt修改组列名
    #输入文件为feautre表result/otutab.txt，实验设计metadata.txt
    #输出为特征表按组的均值-一个实验可能有多种分组方式
    Rscript ${bin}/script/otu_mean.R --input result/otutab.txt \
      --design result/metadata.txt \
      --group Group \
      --output result/otutab_mean.txt
      
    head result/otutab_mean.txt

    #如以平均丰度频率高于千分之一(0.1%)为筛选标准，得到每个组的OTU组合
    awk 'BEGIN{OFS=FS="\t"}{if(FNR==1) {for(i=2;i<=NF;i++) a[i]=$i;} \
        else {for(i=2;i<=NF;i++) if($i>0.1) print $1, a[i];}}' \
        result/otutab_mean.txt > result/alpha/otu_group_exist.txt
    head result/alpha/otu_group_exist.txt
    # 结果可以直接在http://www.ehbio.com/ImageGP绘制Venn、upSetView和Sanky


## 7. Beta多样性 Beta diversity

    #结果有多个文件，需要目录
    mkdir -p result/beta/
    #基于OTU构建进化树 Make OTU tree
    time usearch -cluster_agg result/otus.fa -treeout result/otus.tree # 8s
    #生成5种距离矩阵：bray_curtis, euclidean, jaccard, manhatten, unifrac
    time usearch -beta_div result/otutab_rare.txt -tree result/otus.tree \
      -filename_prefix result/beta/ # 1s


## 8. 物种注释结果分类汇总 Taxonomy summary

    #统计门纲目科属，使用 rank参数 p c o f g，为phylum, class, order, family, genus缩写
    mkdir -p result/tax
    for i in p c o f g;do
      time usearch -sintax_summary result/otus.sintax \
      -otutabin result/otutab_rare.txt \
      -rank ${i} \
      -output result/tax/sum_${i}.txt
    done
    #Taxonomy中异常字符
    sed -i 's/(//g;s/)//g;s/\"//g;s/\#//g;s/\/Chloroplast//g' result/tax/sum_*.txt
    
    # 列出所有文件
    ls result/tax/sum_*.txt
    
    # 查看文件内容
    head result/tax/sum_g.txt

    #OTU对应物种注释2列格式：去除sintax中置信值，只保留物种注释，替换:为_，删除引号
    cut -f 1,4 result/otus.sintax \
      |sed 's/\td/\tk/;s/:/__/g;s/,/;/g;s/"//g;s/\/Chloroplast//' \
      > result/taxonomy2.txt
    head -n3 result/taxonomy2.txt

    #OTU对应物种8列格式：
    #注意注释是非整齐的，由于新物种只是相近而不同
    #生成物种表格：注意OTU/ASV中会有末知为空白，补齐分类未知新物种为Unassigned
    awk 'BEGIN{OFS=FS="\t"}{delete a; a["k"]="Unassigned";a["p"]="Unassigned";a["c"]="Unassigned";a["o"]="Unassigned";a["f"]="Unassigned";a["g"]="Unassigned";a["s"]="Unassigned";\
      split($2,x,";");for(i in x){split(x[i],b,"__");a[b[1]]=b[2];} \
      print $1,a["k"],a["p"],a["c"],a["o"],a["f"],a["g"],a["s"];}' \
      result/taxonomy2.txt > temp/otus.tax
    # mac自带的sed不识别tab \t
    sed 's/;/\t/g;s/.__//g;' temp/otus.tax|cut -f 1-8 | \
      sed '1 s/^/OTUID\tKingdom\tPhylum\tClass\tOrder\tFamily\tGenus\tSpecies\n/' \
      > result/taxonomy.txt
    head -n3 result/taxonomy.txt


## 9. 有参比对——功能预测，如Greengenes，可用于picurst, bugbase分析

    mkdir -p result/gg/
    #与GG所有97% OTUs比对，用于功能预测

    #方法1. usearch比对更快，但文件超限报错选方法2
    usearch -otutab temp/filtered.fa -otus ${db}/gg/97_otus.fasta \
    	-otutabout result/gg/otutab.txt -threads 4
    	
    head -n3 result/gg/otutab.txt
    #79.9%, 4核时6m

    # #方法2. vsearch比对，更准但更慢，但并行更强
    # time vsearch --usearch_global temp/filtered.fa --db ${db}/gg/97_otus.fasta \
    #   --otutabout result/gg/otutab.txt --id 0.97 --threads 4
    # #80.87%, 49m, 594Mb

    #统计
    usearch -otutab_stats result/gg/otutab.txt -output result/gg/otutab.stat
    cat result/gg/otutab.stat


## 10. 空间清理及数据提交

    #删除中间大文件
    rm -rf temp/*.fq
    #原始数据及时压缩节省空间并上传数据中心备份, 54s
    gzip seq/*

    # 分双端统计md5值，用于数据提交
    cd seq
    md5sum *_1.fq.gz > md5sum1.txt
    md5sum *_2.fq.gz > md5sum2.txt
    paste md5sum1.txt md5sum2.txt | awk '{print $2"\t"$1"\t"$4"\t"$3}' | sed 's/*//g' > ../result/md5sum.txt
    rm md5sum*
    cd ..
    cat result/md5sum.txt



# 23、LEfSe、STAMP统计分析和可视化

## 1. LEfSe windows + 网页分析

    #1. 23LEfSe目录中准备otutab.txt, metadata.txt, taxonomy.txt三个文件；
    #2. Rstudio打开format2lefse.Rmd并Knit生成输入文件和可重复计算网页；
    # 可选命令行生成输入文件
    Rscript ${bin}/script/format2lefse.R -h
    Rscript ${bin}/script/format2lefse.R --input result/otutab.txt \
      --taxonomy result/taxonomy.txt --design result/metadata.txt \
      --group Group --threshold 0.1 \
      --output result/LEfSe.txt
    # MAC下可能会出现这个错误
    # Error in `$<-.data.frame`(`*tmp*`, "Phylum", value = character(0)) :
    # replacement has 0 rows, data has 163
    # Calls: format2lefse -> $<- -> $<-.data.frame
    # 停止执行
    # 这是文件result/taxonomy.txt格式不对导致的
    # MAC的sed不支持\t，建议安装gnu-sed，或打开文件自己修改
    #3. 打开LEfSe.txt并在线提交 http://www.ehbio.com/ImageGP/index.php/Home/Index/LEFSe.html

## 2. STAMP文件生成

    #1. 23STAMP目录中准备otutab.txt和taxonomy.txt文件；
    #2. Rstudio打开format2stamp.Rmd并Knit生成输入文件和可重复计算网页；
    # 可选命令行生成输入文件
    Rscript ${bin}/script/format2stamp.R -h
    mkdir -p result/stamp
    Rscript ${bin}/script/format2stamp.R --input result/otutab.txt \
      --taxonomy result/taxonomy.txt --threshold 0.1 \
      --output result/stamp/tax
    #3. 打开LEfSe.txt并在线提交 http://www.ehbio.com/ImageGP/index.php/Home/Index/LEFSe.html

## 3. LEfSe Linux服务器下分析(选学)

    #注：Windows下无法运行
    mkdir -p ~/amplicon/23LEfSe
    cd ~/amplicon/23LEfSe
    #上传 LEfSe.txt 文件

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


# 24 R语言差异分析

    # 准备门水平曼哈顿图例，n10代表前10门作为图例
    tail -n+2 result/tax/sum_p.txt | grep -v 'Unassigned' \
      | cut -f 1 | head -n10 > result/tax/tax_phylum.top
    # 差异比较，指定列名和两组名，用短线连接两组，如A-B
    Rscript ${bin}/script/compare.R \
      --input result/otutab.txt --design result/metadata.txt \
      --taxonomy result/taxonomy.txt --topNtax result/tax/tax_phylum.top \
      --group Group --compare KO-WT --output result/compare/ \
      --pvalue 0.05 --fdr 0.1 --width 8 --height 5

# 24 R语言多样性分析

## 1. Alpha多样性箱线图
    # 查看帮助
    Rscript ${bin}/script/alpha_boxplot.R -h
    # 完整参数
    Rscript ${bin}/script/alpha_boxplot.R --alpha_index richness \
      --input result/alpha/alpha.txt --design result/metadata.txt \
      --group Group --output result/alpha/ \
      --width 89 --height 59



# 31. PICRUSt功能预测 (Linux下分析、选学)

    #推荐使用 http://www.ehbio.com/ImageGP 在线分析
    #有Linux服务器用户可参考以下代码搭建本地流程
    cd ~/amplicon/31PICRUSt

    #上传gg/otutab.txt至当前目录
    #转换为OTU表通用格式，方便下游分析和统计
    biom convert -i otutab.txt \
        -o otutab.biom \
        --table-type="OTU table" --to-json
    #校正拷贝数
    normalize_by_copy_number.py -i otutab.biom \
        -o otutab_norm.biom \
        -c /db/picrust/16S_13_5_precalculated.tab.gz
    #预测宏基因组KO表，biom方便下游归类，txt方便查看分析
    predict_metagenomes.py -i otutab_norm.biom \
        -o ko.biom \
        -c /db/picrust/ko_13_5_precalculated.tab.gz
    predict_metagenomes.py -f -i otutab_rare.biom \
        -o ko.txt \
        -c /db/picrust/ko_13_5_precalculated.tab.gz

    #按功能级别分类汇总, -c输出KEGG_Pathways，分1-3级
    sed  -i '/#Constru/d;s/#OTU //' ko.txt
    num=`tail -n1 ko.txt|wc -w`
    for i in 1 2 3;do
      categorize_by_function.py -f -i ko.biom -c KEGG_Pathways -l ${i} -o ko${i}.txt
      sed  -i '/#Constru/d;s/#OTU //' ko${i}.txt
      paste <(cut -f $num ko${i}.txt) <(cut -f 1-$[num-1] ko${i}.txt) > ko${i}.spf
    done
    wc -l ko*.spf


# 32. Bugbase细菌表型预测


## 1. Bugbase在Windows下分析
    bugbase=${bin}/script/BugBase
    rm -rf result/bugbase/
    Rscript ${bugbase}/bin/run.bugbase.r -L ${bugbase} \
      -i result/gg/otutab.txt -m result/metadata.txt -c Group -o result/bugbase/

## 2. 使用 http://www.ehbio.com/ImageGP 在线分析

## 3. Linux下安装和命令行

### 1. 软件安装(仅一次)

    #有两种方法可选，推荐第一种，可选第二种，仅需运行一次

    #方法1. git下载，需要有git
    git clone https://github.com/knights-lab/BugBase

    #方法2. 下载并解压
    wget https://github.com/knights-lab/BugBase/archive/master.zip
    mv master.zip BugBase.zip
    unzip BugBase.zip
    mv BugBase-master/ BugBase

    #安装依赖包
    cd BugBase
    export BUGBASE_PATH=`pwd`
    export PATH=$PATH:`pwd`/bin
    #安装了所有依赖包
    run.bugbase.r -h
    #测试数据
    run.bugbase.r -i doc/data/HMP_s15.txt -m doc/data/HMP_map.txt -c HMPBODYSUBSITE -o output


### 2. 准备输入文件

    cd ~/amplicon/32BugBase
    #输入文件：基于greengene OTU表的biom格式(本地分析支持txt格式无需转换)和mapping file(metadata.txt首行添加#)
    #上传实验设计+刚才生成的otutab_gg.txt
    #生成在线分析使用的biom1.0格式
    biom convert -i ../result/gg/otutab.txt -o otutab_gg.biom --table-type="OTU table" --to-json
    sed '1 s/^/#/' ../result/metadata.txt > MappingFile.txt
    #下载otutab_gg.biom 和 MappingFile.txt用于在线分析

### 3. 本地分析

    export BUGBASE_PATH=`pwd`
    export PATH=$PATH:`pwd`/bin
    run.bugbase.r -i otutab_gg.txt -m MappingFile.txt -c Group -o phenotype/


# 32. FAPROTAXS元素循环和Bugbase表型预测

## 方法1. 在线分析

    #推荐使用 http://www.ehbio.com/ImageGP 在线分析

## 方法2. Linux下分析、选学

    cd amplicon/32FAPROTAX

### 1. 软件安装

    #下载软件1.1版， June 10, 2017更新数据库
    wget -c https://pages.uoregon.edu/slouca/LoucaLab/archive/FAPROTAX/SECTION_Download/MODULE_Downloads/CLASS_Latest%20release/UNIT_FAPROTAX_1.2/FAPROTAX_1.2.zip
    #解压
    unzip FAPROTAX_1.2.zip

    #测试是否可运行，弹出帮助即正常工作
    python FAPROTAX_1.2/collapse_table.py

    #如果报错，一般提示缺少numpy，可使用conda安装依赖包
    conda install numpy
    conda install biom

### 2. 制作输入OTU表

    #txt转换为biom json格式
    biom convert -i otutab_rare.txt -o otutab_rare.biom --table-type="OTU table" --to-json
    #添加物种注释
    biom add-metadata -i otutab_rare.biom --observation-metadata-fp taxonomy2.txt \
      -o otutab_rare_tax.biom --sc-separated taxonomy \
      --observation-header OTUID,taxonomy
    #指定输入文件、物种注释、输出文件、注释列名、属性列名

### 3. FAPROTAX功能预测

    #python运行collapse_table.py脚本、输入带有物种注释OTU表tax.biom、
    #-g指定数据库位置，物种注释列名，输出过程信息，强制覆盖结果，结果文件和细节
    #下载faprotax.txt，配合实验设计可进行统计分析
    #faprotax_report.txt查看每个类别中具体来源哪些OTUs
    python FAPROTAX_1.2/collapse_table.py -i otutab_rare_tax.biom \
      -g FAPROTAX_1.2/FAPROTAX.txt \
      --collapse_by_metadata 'taxonomy' -v --force \
      -o faprotax.txt -r faprotax_report.txt

### 4. 制作OTU对应功能注释有无矩阵

    # 对OTU注释行，及前一行标题进行筛选
    grep 'ASV_' -B 1 faprotax_report.txt | grep -v -P '^--$' > faprotax_report.clean
    # 筛选Perl脚本将数据整理为表格，搜索我的github(YongxinLiu)或32FAPROTAX目录
    ./faprotax_report_sum.pl -i faprotax_report.clean -o faprotax_report
    # 查看功能有无矩阵，-S不换行
    less -S faprotax_report.mat



# 33、MachineLearning机器学习

## Silme2随机森林/Adaboost

    # 使用实战
    cd 33MachineLearning/slime2
    #使用adaboost计算10000次(16.7s)，推荐千万次
    ./slime2.py otutab.txt design.txt --normalize --tag ab_e4 ab -n 10000
    #使用RandomForest计算10000次(14.5s)，推荐百万次，支持多线程
    ./slime2.py otutab.txt design.txt --normalize --tag rf_e4 rf -n 10000
    cd ../../
  
    #下载安装
    cd ~/software/
    wget https://github.com/swo/slime2/archive/master.zip
    mv master.zip slime2.zip
    unzip slime2.zip
    mv slime2-master/ slime2
    cp slime2/slime2.py ~/bin/
    chmod +x ~/bin/slime2.py

    #安装依赖包
    sudo pip3 install --upgrade pip
    sudo pip3 install pandas
    sudo pip3 install sklearn



# 34. Evolution进化树

    cd ${wd}
    mkdir -p result/tree
    cd ${wd}/result/tree

## 1. 筛选高丰度、指定ASV序列

    #方法1. 按丰度筛选：筛选树高丰度OTU，一般选0.001或0.005，且OTU数量在30-150个范围内
    #统计OTU表中OTU数量，如总计2631个
    tail -n+2 ../otutab_rare.txt | wc -l
    #按相对丰度0.2%筛选高丰度OTU
    usearch -otutab_trim ../otutab_rare.txt \
        -min_otu_freq 0.002 \
        -output otutab.txt
    #统计筛选OTU表特征数量，总计79或80个
    tail -n+2 otutab.txt | wc -l
    #提取ID用于提取序列
    cut -f 1 otutab.txt | sed '1 s/#OTU ID/OTUID/' > otutab_high.id

    #方法2. 按数量筛选
    #按丰度排序，默认由大到小
    usearch -otutab_sortotus ../otutab_rare.txt  \
        -output otutab_sort.txt
    #提取高丰度中指定Top数量的OTU ID，如Top100,
    sed '1 s/#OTU ID/OTUID/' otutab_sort.txt \
        | head -n101 > otutab.txt
    cut -f 1 otutab.txt > otutab_high.id

    #筛选高丰度菌/指定差异菌对应OTU序列
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
	  time muscle -in otus.fa -out otus_aligned.fas

    ### 方法1. 利用IQ-TREE快速构建ML进化树，2m
    mkdir -p iqtree
    time iqtree -s otus_aligned.fas \
        -bb 1000 -redo -alrt 1000 -nt AUTO \
        -pre iqtree/otus

    ### 方法2. FastTree快速建树(Linux)
    # 注意FastTree软件输入文件为fasta格式的文件，而不是通常用的Phylip格式。输出文件是Newick格式。
    # 该方法适合于大数据，例如几百个OTUs的系统发育树！
    # Ubuntu上安装fasttree可以使用`apt install fasttree`
    # fasttree -gtr -nt otus_aligned.fa > otus.nwk


## 3. 进化树美化

    # 访问http://itol.embl.de/，上传otus.nwk，再拖拽下方生成的注释方案于树上即美化

    ## 方案1. 外圈颜色、形状分类和丰度方案
    # annotation.txt OTU对应物种注释和丰度，
    # -a 找不到输入列将终止运行（默认不执行）-c 将整数列转换为factor或具有小数点的数字，-t 偏离提示标签时转换ID列，-w 颜色带，区域宽度等， -D输出目录，-i OTU列名，-l OTU显示名称如种/属/科名，
    # cd ${wd}/result/tree
    Rscript ${bin}/script/table2itol.R -a -c double -D plan1 -i OTUID -l Genus -t %s -w 0.5 annotation.txt
    # 生成注释文件中每列为单独一个文件

    ## 方案2. 生成丰度柱形图注释文件
    Rscript ${bin}/script/table2itol.R -a -d -c none -D plan2 -b Phylum -i OTUID -l Genus -t %s -w 0.5 annotation.txt

    ## 方案3. 生成热图注释文件
    Rscript ${bin}/script/table2itol.R -c keep -D plan3 -i OTUID -t %s otutab.txt

    ## 方案4. 将整数转化成因子生成注释文件
    Rscript ${bin}/script/table2itol.R -a -c factor -D plan4 -i OTUID -l Genus -t %s -w 0 annotation.txt

    # 返回工作目录
    cd ${wd}



# 常见问题

## 1. 文件phred质量错误——Fastq质量值64转33

    #查看64位格式文件，质量值多为小写字母
    head -n4 FAQ/Q64Q33/test_64.fq
    #转换质量值64编码格式为33
    vsearch --fastq_convert FAQ/Q64Q33/test_64.fq \
        --fastqout FAQ/test.fq \
        --fastq_ascii 64 --fastq_asciiout 33
    #查看转换后33编码格式，质量值多为大写字母
    head -n4 FAQ/test.fq

## 2. 序列双端已经合并——单端序列重命名/添加样本名

    #查看文件序列名
    head -n1 FAQ/test.fq
    #序列按样本命名，并输出到新文件夹
    mkdir -p FAQ/relabel
    vsearch --fastq_convert FAQ/test.fq \
        --fastqout FAQ/relabel/WT1.fq --relabel WT1.
    #查看转换后33编码格式，质量值多为大写字母
    head -n1 FAQ/relabel/WT1.fq

## 3. 数据过大无法使用usearch聚类或去噪-vsearch

    #备选vsearch生成OTU，但无自动de novo去嵌合功能
    #仅限usearch免费版受限时(可通过提高minuniquesize参数减少数据量)使用，不推荐
    #重命名、相似97%，不屏蔽，输入和输出count
    vsearch --cluster_size temp/uniques.fa  \
     --centroids temp/otus.fa \
     --relabel OTU_ --id 0.97 --qmask none --sizein --sizeout
    #5s Clusters: 1062
    #vsearch还需连用--uchime3_denovo


## 4. Reads count整数值如何标准化为相对丰度

    #求取各个OTU在对应样品的丰度频率
    usearch -otutab_counts2freqs result/otutab_rare.txt \
        -output result/otutab_rare_freq.txt

## 5. 运行R提示write.table Permission denied

    #例如报错信息示例如下：
    Error in file(file, ifelse(append, "a", "w")) :
    Calls: write.table -> file
    : Warning message:
    In file(file, ifelse(append, "a", "w")) :
      'result/raw/otutab_nonBac.txt': Permission denied
    #翻译为写入文件无权限，一般为目标文件正在被打开，请关闭重试

## 6. 文件批量命名

    # 注意修改路径
    cd /c/project/seq
    ls > ../filelist.txt
    # 编辑列表，第二名为最终命名，确定名称唯一
    # 检查手动命名是否唯一
    cut -f 2 ../filelist.txt |wc -l
    cut -f 2 ../filelist.txt | sort | uniq |wc -l
    # 如果两次结果一致，则命名非冗余
    awk '{system("mv "$1" "$2)}' ../filelist.txt

## 7. Rstudio中Terminal找不到Linux命令

    # 需要把 C:\Program Files\Git\usr\bin 目录添加到系统环境变量
    # 注意win10系统是一个目录一行；win7中多个目录用分号分隔，注意向后添加目录


## 8. 测试均值丢失组
    cd /c/amplicon/FAQ/merge
    #按组求均值，需根据实验设计metadata.txt修改组列名
    #输入文件为feautre表result/otutab.txt，实验设计metadata.txt
    #输出为特征表按组的均值-一个实验可能有多种分组方式
    Rscript /c/amplicon/22Pipeline/script/otu_mean.R

    #如以平均丰度频率高于0.05%为筛选标准，得到每个组的OTU组合
    awk 'BEGIN{OFS=FS="\t"}{if(FNR==1) {for(i=2;i<=NF;i++) a[i]=$i;} \
        else {for(i=2;i<=NF;i++) if($i>0.05) print $1, a[i];}}' \
        result/otutab_mean_Genotype.txt > alpha/otu_group_exist.txt
    # 结果可以直接在http://www.ehbio.com/ImageGP绘制Venn、upSetView和Sanky

## 9. usearch/vsearch 生成OTU表时无法匹配
    #是原始序列方向错误，将序列需要取反向互补
    vsearch --fastx_revcomp ../FAQ/filtered_test.fa \
      --fastaout ../FAQ/filtered_test_RC.fa
    # 再分析
    usearch -otutab ../FAQ/filtered_test_RC.fa -otus db/gg/97_otus.fasta \
    	-otutabout gg/otutab.txt -threads 6

## 10. 检查文件windows换行符和删除

    cd /c/amplicon/FAQ/190614_ITS_taxSum_0field
  	i=g
    usearch -sintax_summary sintax.txt \
    -otutabin otutab_rare.txt \
    -rank ${i} \
    -output sum_${i}.txt
