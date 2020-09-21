[TOC]

    # 名称 Name: 易扩增子(EasyAmplicon)
    # 作者 Authors: 刘永鑫(Yong-Xin Liu), 陈同(Tong Chen)，周欣(Xin Zhou)等
    # 版本 Version: v1.10
    # 更新 Update: 2020-12-04

    # 设置流程EasyAmplicon(ea)和工作目录(work directory, wd)，添加环境变量并进入wd
    # **每次打开Rstudio必须运行下面4行**
    ea=/d/home/github/EasyAmplicon
    wd=/c/test
    PATH=$PATH:${ea}/win
    cd ${wd}


# 22、扩增子分析流程 Analysis pipeline for 16S amplicon 

    # 系统要求 System requirement: Windows 10 / Mac OS 10.12+ / Ubuntu 18.04+
    # 依赖软件和数据 Sofware and database dependencies: 复制public目录到C盘
    # gitforwidnows 2.28.0 http://gitforwindows.org/(Windows only)
    # R 4.0.2 https://www.r-project.org/
    # Rstudio 1.3.1056 https://www.rstudio.com/products/rstudio/download/#download
    # vsearch v2.15.0 https://github.com/torognes/vsearch/releases
    # usearch v10.0.240 https://www.drive5.com/usearch/download.html

    # 运行前准备
    # 1. 将EasyAmplicon目录下发 /复制到指定位置，如Widnows的C盘(C:/) 或 Mac/Linux服务器家目录(~/)
    # 2. 按EasyAmplicon中Readme.md中说明安装依赖软件和数据库
    # 3. 准备流程脚本pipeline.sh、样本元数据metadata.txt和测序数据seq/*.fq.gz
 
## 1. 准备输入数据(测试数据是流程正对照)

    #1. 分析流程pipeline.sh，每个项目复制一份，可修改实现分析过程全记录(可重复分析)
    #2. 样本元数据/实验设计metadata.txt
    wget -c http://210.75.224.110/github/EasyAmplicon/data/metadata.txt
    head -n2 metadata.txt
    #3. 原始测序数据保存于seq目录，通常以`.fq.gz`结尾，每个样品一对文件
    # 下载测试数据，按元数据中GSA的CRA(批次)和CRR(样品)编号从公共数据库下载
    # 示例下载单个文件并改名
    # wget -c ftp://download.big.ac.cn/gsa/CRA002352/CRR117575/CRR117575_f1.fq.gz -O seq/KO1_1.fq.gz
    # 按实验设计编号批量下载并改名，速度单位10M/s，家500K/s
    mkdir -p seq
    awk '{system("wget -c ftp://download.big.ac.cn/gsa/"$5"/"$6"/"$6"_f1.fq.gz -O seq/"$1"_1.fq.gz")}' \
        <(tail -n+2 metadata.txt)
    awk '{system("wget -c ftp://download.big.ac.cn/gsa/"$5"/"$6"/"$6"_r2.fq.gz -O seq/"$1"_2.fq.gz")}' \
        <(tail -n+2 metadata.txt)
    ls -sh seq
    #4. 创建临时和结果目录，临时目录分析结束可删除
    mkdir -p temp result

### 1.1. metadata.txt 实验设计文件

    #cat查看前3行，-A显示符号
    cat -A metadata.txt | head -n3
    #windows用户如果结尾有^M，运行sed命令去除，并cat -A检查结果
    # sed -i 's/\r/\n/' metadata.txt
    # cat -A metadata.txt | head -n3

### 1.2. seq/*.fq.gz 原始测序数据

    # 公司返回的测序结果，通常为一个样品一对fq/fastq.gz格式压缩文件
    # 文件名与样品名务必对应：不一致时手工修改，批量改名见"常见问题6"
    # 测序数据是.gz的压缩文件，可以使用zcat/zless查看，连用head查看几行
    # zless/less按页查看，空格翻页、q退出；head查看前10行，-n指定行
    zcat seq/KO1_1.fq.gz|head -n4
    # gz压缩文件可以使用gunzip解压作为USEARCH输入，vsearch可直接读取gz文件作为输入
    # gunzip seq/*.gz
    ls -sh seq/
    # 每行太长，指定查看每行的1-60个字符
    # cut -c 1-60 seq/KO1_1.fq | head -n4

### 1.3. pipeline.sh 流程依赖数据库

    # 数据库第一次使用需要解压，解压过可跳过此段

    # usearchs可用16S/18S/ITS数据库：RDP, SILVA和UNITE，本地文件位置 ${ea}/usearch/
    # usearch数据库database下载页: http://www.drive5.com/sintax
    # 解压rdp用于物种注释和silva用于去嵌合体，*代表任意字符，选择以fa.gz结尾的文件
    # time gunzip ${ea}/usearch/*.fa.gz  # 51s
    # 解压缩数据库
    gunzip ${ea}/usearch/rdp_16s_v16_sp.fa.gz
    # 查看注释数据库的格式
    head -n2 ${ea}/usearch/rdp_16s_v16_sp.fa

    # QIIME greengene 13_8有参数据库用于功能注释: ftp://greengenes.microbio.me/greengenes_release/gg_13_5/gg_13_8_otus.tar.gz
    gunzip ${ea}/gg/97_otus.fasta.gz
    # 不用时压缩结省空间
    # gzip ${ea}/gg/97_otus.fasta


## 2. 合并双端序列并按样品重命名 Merge paired reads and label samples

    # 以WT1单样品合并为例
    vsearch --fastq_mergepairs seq/WT1_1.fq.gz \
      --reverse seq/WT1_2.fq.gz \
      --fastqout temp/WT1.merged.fq \
      --relabel WT1.

    #依照实验设计批处理并合并。
    #tail -n+2去表头，cut -f 1取第一列，即获得样本列表。5万对序列的18个样本合并耗时约2分钟。
    #Windows下复制命令Ctrl+C在Linux下终止，为防止异常中断，结尾添加&使命令转后台。
    for i in `tail -n+2 metadata.txt | cut -f 1`;do
      vsearch --fastq_mergepairs seq/${i}_1.fq.gz --reverse seq/${i}_2.fq.gz \
      --fastqout temp/${i}.merged.fq --relabel ${i}.
    done &

    # 己经合并样本直接改名，接入分析流程；替换上面vsearch --fastq_mergepairs命令为
    # usearch -fastx_relabel seq/${i}.fq -fastqout temp/${i}.merged.fq -prefix ${i}.
    # 另一种方法参考“常见问题2”

    #合并所有样品至同一文件
    cat temp/*.merged.fq > temp/all.fq
    #查看文件大小634M，软件不同版本结果略有差异
    ls -lsh temp/all.fq
    #查看序列名.之前是否为样本名，样本名绝不允许有.
    head -n 6 temp/all.fq|cut -c1-60


## 3. 引物切除和质量控制 Cut primers and quality filter

    # 采用等长的方式切除引物，引物外侧如果有标签，标签的长度需要计算在内。
    # 如本例的结果为左端标签10 bp + 正向引物V5 19 bp共为29 bp，反向引物V7 18 bp。
    # Cut barcode 10bp + V5 19bp in left and V7 18bp in right
    # 注意：务必清楚实验设计中引物和标签长度，如果引物已经去除，可在下方参数处填0表示无需去除。76万条序列37s
    vsearch --fastx_filter temp/all.fq \
      --fastq_stripleft 29 --fastq_stripright 18 \
      --fastq_maxee_rate 0.01 \
      --fastaout temp/filtered.fa

    # 查看文件了解fa文件格式
    head temp/filtered.fa


## 4. 序列去冗余并挑选代表序列(OTU/ASV) Dereplicate and cluster/denoise

### 4.1 序列去冗余 Dereplication

    # 去冗余数据量起码降低1个数量级，减小下游分析工作量，也更适合基于丰度鉴定真实特征序列，如可操作分类单元(OTUs)或扩增序列变化(ASVs)。
    # 参数--miniuniqusize设置使用序列的最小出现次数，默认为8，此处设置为10，推荐最小为总数据量的百万分之一，可实现去除低丰度噪音并增加计算速度。
    # -sizeout输出丰度, --relabel添加序列前缀。
    vsearch --derep_fulllength temp/filtered.fa \
      --output temp/uniques.fa \
      --relabel Uni --minuniquesize 10 --sizeout
    #高丰度非冗余序列非常小(<2Mb),名称后有size和频率
    ls -lsh temp/uniques.fa
    head -n 2 temp/uniques.fa

### 4.2 聚类OTU/去噪ASV Cluster OTUs / denoise ASV

    # OTU和ASV统称为特征(Feature)，它们的区别是：
    # OTU通常按97%聚类后挑选最高丰度或中心的代表性序列；
    # ASV是基于序列进行去噪(排除或校正错误序列，并挑选丰度较高的可信序列)作为代表性序列

    #有两种方法：推荐unoise3去噪获得单碱基精度ASV，传统的97%聚类OTU (属水平精度)供备选
    #usearch两种特征挑选方法均自带de novo去嵌合体

    #方法1. 97%聚类OTU，适合大数据/ASV规律不明显/reviewer要求
    #结果耗时6s, 产生878 OTUs, 去除320 chimeras
    # usearch -cluster_otus temp/uniques.fa \
    #  -otus temp/otus.fa \
    #  -relabel OTU_

    #方法2. ASV去噪 Denoise: predict biological sequences and filter chimeras
    #59s, 2920 good, 227 chimeras
    usearch -unoise3 temp/uniques.fa \
      -zotus temp/zotus.fa
    #修改序列名：Zotu为改为ASV方便识别
    sed 's/Zotu/ASV_/g' temp/zotus.fa > temp/otus.fa
    head -n 2 temp/otus.fa

    #方法3. 数据过大无法使用usearch时，备选vsearch方法见"常见问题3"

### 4.3 基于参考去嵌合 Reference-based chimera detect

    # 不推荐，容易引起假阴性，因为参考数据库无丰度信息，
    # 而de novo去嵌合时要求亲本丰度为嵌合体16倍以上防止假阴性
    # 因为已知序列不会被去除，数据库选择越大越合理，假阴性率最低
    mkdir -p result/raw

    # 方法1. vsearch+rdp去嵌合(快但容易假阴性)，或
    # silva去嵌合(silva_16s_v123.fa)，推荐(慢，耗时15m ~ 3h，但更好)
    vsearch --uchime_ref temp/otus.fa \
      -db ${ea}/usearch/rdp_16s_v16_sp.fa \
      --nonchimeras result/raw/otus.fa
    # RDP: 51s, 250 (8.6%) chimeras; SILVA：10m, 255 (8.7%) chimeras
    # Win vsearch结果添加了windows换行符^M需删除，mac不要执行此命令
    sed -i 's/\r//g' result/raw/otus.fa

    # 方法2. 不去嵌合
    cp -f temp/otus.fa result/raw/otus.fa


## 5. 特征表生成和筛选 Feature table create & filter

### 5.1 生成特征表 Creat Feature table

    # 方法1. usearch生成特征表，小样本(<30)快；但大样本受限且多线程效率低，84.1%, 4核1m
    # usearch -otutab temp/filtered.fa -otus result/raw/otus.fa \
    #   -otutabout result/raw/otutab.txt -threads 4

    # 方法2. vsearch生成特征表
    vsearch --usearch_global temp/filtered.fa --db result/raw/otus.fa \
    	--otutabout result/raw/otutab.txt --id 0.97 --threads 4
    #660953 of 761432 (86.80%)可比对，耗时8m
    # windows用户删除换行符^M
    sed -i 's/\r//' result/raw/otutab.txt
    head -n3 result/raw/otutab.txt |cat -A


### 5.2 物种注释-去除质体和非细菌/古菌并统计比例(可选) Remove plastid and non-Bacteria

    # RDP物种注释(rdp_16s_v16_sp)数据库小，分类速度极快，但缺少完整真核来源数据，耗时15s;
    # SILVA数据库(silva_16s_v123.fa)更好注释真核、质体序列，3h
    # --sintax_cutoff设置分类的可信度阈值，常用0.6/0.8，越大注释比例越低。
    vsearch --sintax result/raw/otus.fa --db ${ea}/usearch/rdp_16s_v16_sp.fa \
      --tabbedout result/raw/otus.sintax --sintax_cutoff 0.6

    # 原始特征表行数
    wc -l result/raw/otutab.txt
    # 为去除16S rDNA测序中的非特异性扩增和按需分配污染，
    #细菌和古菌(原核生物)、以及去除叶绿体和线粒体并统计比例，输出筛选并按丰度排序的OTU表。
    #输入为特征表(result/raw/otutab.txt)和物种注释(result/raw/otus.sintax)，
    #输出筛选并排序的特征表(result/otutab.txt)、统计污染比例文件(result/raw/otutab_nonBac.txt)和过滤细节(otus.sintax.discard)。
    #注：真菌ITS数据，请改用otutab_filter_nonFungi.R脚本，只筛选注释为真菌的序列。
    # 查看脚本帮助，请运行Rscript ${ea}/script/otutab_filter_nonBac.R -h
    Rscript ${ea}/script/otutab_filter_nonBac.R \
      --input result/raw/otutab.txt \
      --taxonomy result/raw/otus.sintax \
      --output result/otutab.txt\
      --stat result/raw/otutab_nonBac.stat \
      --discard result/raw/otus.sintax.discard
    # 筛选后特征表行数，2921-2904，删除17个ASV
    wc -l result/otutab.txt

    #过滤特征表对应序列
    cut -f 1 result/otutab.txt | tail -n+2 > result/otutab.id
    usearch -fastx_getseqs result/raw/otus.fa \
        -labels result/otutab.id -fastaout result/otus.fa
    #过滤特征表对应序列注释
    awk 'NR==FNR{a[$1]=$0}NR>FNR{print a[$1]}'\
        result/raw/otus.sintax result/otutab.id \
        > result/otus.sintax
    #补齐末尾列
    sed -i 's/\t$/\td:Unassigned/' result/otus.sintax
    # head -n2 result/otus.sintax


    # 方法2. 觉得筛选不合理可以不筛选
    # cp result/raw/otu* result/

    #可选统计方法：OTU表简单统计 Summary OTUs table
    usearch -otutab_stats result/otutab.txt \
      -output result/otutab.stat
    cat result/otutab.stat
    #注意最小值、分位数，或查看result/raw/otutab_nonBac.stat中样本详细数据量，用于重采样

### 5.3 等量抽样标准化 normlize by subsample

    #使用vegan包进行等量重抽样，输入reads count格式Feature表result/otutab.txt
    #可指定输入文件、抽样量和随机数，输出抽平表result/otutab_rare.txt和多样性alpha/vegan.txt
    mkdir -p result/alpha
    Rscript ${ea}/script/otutab_rare.R --input result/otutab.txt \
      --depth 32086 --seed 1 \
      --normalize result/otutab_rare.txt \
      --output result/alpha/vegan.txt
    usearch -otutab_stats result/otutab_rare.txt \
      -output result/otutab_rare.stat
    cat result/otutab_rare.stat


## 6. Alpha多样性 Alpha diversity

### 6.1. 计算多样性指数 Calculate alpha diversity index
    #Calculate all alpha diversity index(Chao1有错误勿用)
    #details in http://www.drive5.com/usearch/manual/alpha_metrics.html
    usearch -alpha_div result/otutab_rare.txt \
      -output result/alpha/alpha.txt

### 6.2. 计算稀释过程的丰富度变化 Rarefaction
    #稀释曲线：取1%-100%的序列中OTUs数量，20s
    #Rarefaction from 1%, 2% .. 100% in richness (observed OTUs)-method fast / with_replacement / without_replacement https://drive5.com/usearch/manual/cmd_otutab_subsample.html
    usearch -alpha_div_rare result/otutab_rare.txt \
      -output result/alpha/alpha_rare.txt \
      -method without_replacement

### 6.3. 筛选各组高丰度菌用于比较

    #计算各特征的均值，有组再求分组均值，需根据实验设计metadata.txt修改组列名
    #输入文件为feautre表result/otutab.txt，实验设计metadata.txt
    #输出为特征表按组的均值-一个实验可能有多种分组方式
    Rscript ${ea}/script/otu_mean.R --input result/otutab.txt \
      --design metadata.txt \
      --group Group --thre 0 \
      --output result/otutab_mean.txt
    head -n3 result/otutab_mean.txt

    #如以平均丰度频率高于千一(0.1%)为筛选标准，可选千五或万五，得到每个组的OTU组合
    awk 'BEGIN{OFS=FS="\t"}{if(FNR==1) {for(i=2;i<=NF;i++) a[i]=$i;} \
        else {for(i=2;i<=NF;i++) if($i>0.1) print $1, a[i];}}' \
        result/otutab_mean.txt > result/alpha/otu_group_exist.txt
    head result/alpha/otu_group_exist.txt
    # 结果可以直接在http://www.ehbio.com/ImageGP绘制Venn、upSetView和Sanky
    # 可在http://ehbio.com/test/venn/中绘图并显示各组共有和特有元素

## 7. Beta多样性 Beta diversity

    #结果有多个文件，需要目录
    mkdir -p result/beta/
    #基于OTU构建进化树 Make OTU tree, 30s
    usearch -cluster_agg result/otus.fa -treeout result/otus.tree
    #生成5种距离矩阵：bray_curtis, euclidean, jaccard, manhatten, unifrac, 3s
    usearch -beta_div result/otutab_rare.txt -tree result/otus.tree \
      -filename_prefix result/beta/ # 1s


## 8. 物种注释及分类汇总 Taxonomy summary

    #OTU对应物种注释2列格式：去除sintax中置信值，只保留物种注释，替换:为_，删除引号
    cut -f 1,4 result/otus.sintax \
      |sed 's/\td/\tk/;s/:/__/g;s/,/;/g;s/"//g;s/\/Chloroplast//' \
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
    ls -sh result/tax/sum_*.txt
    head -n3 result/tax/sum_g.txt


## 9. 有参比对——功能预测，如Greengene，可用于picurst, bugbase分析

    mkdir -p result/gg/
    #与GG所有97% OTUs比对，用于功能预测

    #方法1. usearch比对更快，但文件超限报错选方法2
    usearch -otutab temp/filtered.fa -otus ${ea}/gg/97_otus.fasta \
    	-otutabout result/gg/otutab.txt -threads 12
    #79.9%, 4核时8m；12核5m
    head -n3 result/gg/otutab.txt

    # #方法2. vsearch比对，更准更慢，但并行更强
    # vsearch --usearch_global temp/filtered.fa --db ${ea}/gg/97_otus.fasta \
    #   --otutabout result/gg/otutab.txt --id 0.97 --threads 12
    #80.9%, 12cores 20m, 1core 1h, 594Mb

    #统计
    usearch -otutab_stats result/gg/otutab.txt  \
      -output result/gg/otutab.stat
    cat result/gg/otutab.stat


## 10. 空间清理及数据提交

    #删除中间大文件
    rm -rf temp/*.fq
    
    #短期不用数据库压缩节省空间
    gzip ${ea}/usearch/rdp_16s_v16_sp.fa
    gzip ${ea}/gg/97_otus.fasta

    # 原始序列统计md5值，用于数据提交
    cd seq
    md5sum *_1.fq.gz > md5sum1.txt
    md5sum *_2.fq.gz > md5sum2.txt
    paste md5sum1.txt md5sum2.txt | awk '{print $2"\t"$1"\t"$4"\t"$3}' | \
      sed 's/*//g' > ../result/md5sum.txt
    rm md5sum*
    cd ..
    cat result/md5sum.txt



# 23、R语言多样性和物种分析


## 1. Alpha多样性

### 1.1 Alpha多样性箱线图
    # 查看帮助
    Rscript ${ea}/script/alpha_boxplot.R -h
    # 完整参数，多样性指数可选richness chao1 ACE shannon simpson invsimpson
    Rscript ${ea}/script/alpha_boxplot.R --alpha_index richness \
      --input result/alpha/vegan.txt --design metadata.txt \
      --group Group --output result/alpha/ \
      --width 89 --height 59
    # 使用循环绘制6种常用指数
    for i in `head -n1 result/alpha/vegan.txt|cut -f 2-`;do
      Rscript ${ea}/script/alpha_boxplot.R --alpha_index ${i} \
        --input result/alpha/vegan.txt --design metadata.txt \
        --group Group --output result/alpha/ \
        --width 89 --height 59
    done

    # Alpha多样性柱状图+标准差
    Rscript ${ea}/script/alpha_barplot.R --alpha_index richness \
      --input result/alpha/vegan.txt --design metadata.txt \
      --group Group --output result/alpha/ \
      --width 89 --height 59

### 1.2 稀释曲线
    Rscript ${ea}/script/alpha_rare_curve.R \
      --input result/alpha/alpha_rare.txt --design metadata.txt \
      --group Group --output result/alpha/ \
      --width 89 --height 59

### 1.3 多样性维恩图
    # 三组比较:-f输入文件,-a/b/c/d/g分组名,-w/u为宽高英寸,-p输出文件名后缀
    bash ${ea}/script/sp_vennDiagram.sh \
      -f result/alpha/otu_group_exist.txt \
      -a WT -b KO -c OE \
      -w 3 -u 3 \
      -p WT_KO_OE
    # 四组比较，图和代码见输入文件目录，运行目录为当前项目根目录
    bash ${ea}/script/sp_vennDiagram.sh \
      -f result/alpha/otu_group_exist.txt \
      -a WT -b KO -c OE -d All \
      -w 3 -u 3 \
      -p WT_KO_OE_All

## 2. Beta多样性

### 2.1 距离矩阵热图pheatmap
    # 以bray_curtis为例，-f输入文件,-h是否聚类TRUE/FALSE,-u/v为宽高英寸
    bash ${ea}/script/sp_pheatmap.sh \
      -f result/beta/bray_curtis.txt \
      -H 'TRUE' -u 5 -v 5
    # 添加分组注释，如2，4列的基因型和地点
    cut -f 1-2 metadata.txt > temp/group.txt
    # -P添加行注释文件，-Q添加列注释
    bash ${ea}/script/sp_pheatmap.sh \
      -f result/beta/bray_curtis.txt \
      -H 'TRUE' -u 8.9 -v 5.6 \
      -P temp/group.txt -Q temp/group.txt
    # 距离矩阵与相关类似，可尝试corrplot或ggcorrplot绘制更多样式
    # - [绘图相关系数矩阵corrplot](http://mp.weixin.qq.com/s/H4_2_vb2w_njxPziDzV4HQ)
    # - [相关矩阵可视化ggcorrplot](http://mp.weixin.qq.com/s/AEfPqWO3S0mRnDZ_Ws9fnw)

### 2.2 主坐标分析PCoA
    # 输入文件，选择分组，输出文件，图片尺寸mm，统计见beta_pcoa_stat.txt
    Rscript ${ea}/script/beta_pcoa.R \
      --input result/beta/bray_curtis.txt --design metadata.txt \
      --group Group --output result/beta/bray_curtis.pcoa.pdf \
      --width 89 --height 59

### 2.3 限制性主坐标分析CPCoA
    Rscript ${ea}/script/beta_cpcoa.R \
      --input result/beta/bray_curtis.txt --design metadata.txt \
      --group Group --output result/beta/bray_curtis.cpcoa.pdf \
      --width 89 --height 59

## 3. 物种组成Taxonomy

### 3.1 堆叠柱状图Stackplot
    # 以门(p)水平为例，结果包括output.sample/group.pdf两个文件
    Rscript ${ea}/script/tax_stackplot.R \
      --input result/tax/sum_p.txt --design metadata.txt \
      --group Group --output result/tax/sum_p.stackplot \
      --legend 5 --width 89 --height 59
    # 批量绘制输入包括p/c/o/f/g共5级
    for i in p c o f g; do
    Rscript ${ea}/script/tax_stackplot.R \
      --input result/tax/sum_${i}.txt --design metadata.txt \
      --group Group --output result/tax/sum_${i}.stackplot \
      --legend 8 --width 89 --height 59; done

### 3.2 弦图(圈图)circlize
    # 以纲(class,c)为例，绘制前5组
    i=c
    Rscript ${ea}/script/tax_circlize.R \
      --input result/tax/sum_${i}.txt --design metadata.txt \
      --group Group --legend 5
    # 结果位于当前目录circlize.pdf(随机颜色)，circlize_legend.pdf(指定颜色+图例)
    # 移动并改名与分类级一致
    mv circlize.pdf result/tax/sum_${i}.circlize.pdf
    mv circlize_legend.pdf result/tax/sum_${i}.circlize_legend.pdf

### 3.3 树图treemap/maptree
    # 多层级包含物种关系，输入特征表和物种注释，输出树图，指定包含特征数量和图片宽高
    Rscript ${ea}/script/tax_maptree.R \
      --input result/otutab.txt --taxonomy result/taxonomy.txt \
      --output result/tax/tax_maptree.pdf \
      --topN 200 --width 183 --height 118


# 24、差异比较

## 24-1. R语言差异分析

### 1.1 差异比较Difference comparison
    mkdir -p result/compare/
    # 输入特征表、元数据；指定分组列名、比较组和丰度
    # 选择方法wilcox/t.test/edgeR、pvalue和fdr和输出目录
    compare="KO-WT"
    Rscript ${ea}/script/compare.R \
      --input result/otutab.txt --design metadata.txt \
      --group Group --compare ${compare} --threshold 0.1 \
      --method edgeR --pvalue 0.05 --fdr 0.2 \
      --output result/compare/

### 1.2 火山图
    # 输入compare.R的结果，输出火山图带数据标签，可指定图片大小
    Rscript ${ea}/script/compare_volcano.R \
      --input result/compare/${compare}.txt \
      --output result/compare/${compare}.txt.volcano.pdf \
      --width 89 --height 59

### 1.3 热图
    # 输入compare.R的结果，筛选列数，指定元数据和分组、物种注释，图大小英寸和字号
    bash ${ea}/script/compare_heatmap.sh -i result/compare/${compare}.txt -l 7 \
       -d metadata.txt -A Group \
       -t result/taxonomy.txt \
       -w 8 -h 5 -s 7 \
       -o result/compare/${compare}.txt

### 1.4 曼哈顿图
    # i差异比较结果,t物种注释,p图例,w宽,v高,s字号,l图例最大值
    bash ${ea}/script/compare_manhattan.sh -i result/compare/${compare}.txt \
       -t result/taxonomy.txt \
       -p result/tax/sum_p.txt \
       -w 183 -v 59 -s 7 -l 10 \
       -o result/compare/${compare}.manhattan.p.pdf
    # 上图只有6个门，切换为纲c和-L Class展示细节
    bash ${ea}/script/compare_manhattan.sh -i result/compare/${compare}.txt \
       -t result/taxonomy.txt \
       -p result/tax/sum_c.txt \
       -w 130 -v 59 -s 7 -l 10 -L Class \
       -o result/compare/${compare}.manhattan.c.pdf

## 24-2. STAMP输入文件准备

### 2.1 命令行生成输入文件
    Rscript ${ea}/script/format2stamp.R -h
    mkdir -p result/stamp
    Rscript ${ea}/script/format2stamp.R --input result/otutab.txt \
      --taxonomy result/taxonomy.txt --threshold 0.1 \
      --output result/stamp/tax

### 2.2 Rmd生成输入文件
    #1. 24compare/stamp目录中准备otutab.txt和taxonomy.txt文件；
    #2. Rstudio打开format2stamp.Rmd，设置参数；
    #3. 点击Knit在当前目录生成stamp输入文件和可重复计算网页。


## 24-3. LEfSe输入文件准备

### 3.1. 命令行生成文件

    # 可选命令行生成输入文件
    Rscript ${ea}/script/format2lefse.R -h
    mkdir -p result/lefse
    Rscript ${ea}/script/format2lefse.R --input result/otutab.txt \
      --taxonomy result/taxonomy.txt --design metadata.txt \
      --group Group --threshold 0.1 \
      --output result/lefse/LEfSe

### 3.2 Rmd生成输入文件
    #1. 24Compare/LEfSe目录中准备otutab.txt, metadata.txt, taxonomy.txt三个文件；
    #2. Rstudio打开format2lefse.Rmd并Knit生成输入文件和可重复计算网页；

### 3.3 LEfSe分析
    #方法1. 打开LEfSe.txt并在线提交 http://www.ehbio.com/ImageGP/index.php/Home/Index/LEFSe.html
    #方法2. LEfSe本地分析(限Linux服务器、选学)，参考代码见附录
    #方法3. LEfSe官网在线使用



# 25、QIIME 2分析流程
    # 代码详见 25QIIME2/pipeline_qiime2.sh



# 31、功能预测


## 1. PICRUSt功能预测

    # 方法1. 使用 http://www.ehbio.com/ImageGP 在线分析 gg/otutab.txt
    # 方法2. Linux服务器用户可参考附录2
    # 然后结果使用STAMP/R进行差异比较

    # (可选)PICRUSt2(Linux/Windows下Linux子系统，要求>16GB内存)
    # 安装
    conda create -n picrust2 -c bioconda -c conda-forge picrust2=2.3.0_b
    # 加载环境
    conda activate picrust2
    # 进入工作目录
    wd=/mnt/c/amplicon/result/
    cd ${wd}
    # 运行流程
    picrust2_pipeline.py -s otus.fa -i otutab.txt \
      -o picrust2 -p 8
    # 添加EC/KO/Pathway注释
    add_descriptions.py -i EC_metagenome_out/pred_metagenome_unstrat.tsv.gz -m EC \
      -o EC_metagenome_out/pred_metagenome_unstrat_descrip.tsv.gz
    add_descriptions.py -i KO_metagenome_out/pred_metagenome_unstrat.tsv.gz -m KO \
      -o KO_metagenome_out/pred_metagenome_unstrat_descrip.tsv.gz 
    add_descriptions.py -i pathways_out/path_abun_unstrat.tsv.gz -m METACYC \
      -o pathways_out/path_abun_unstrat_descrip.tsv.gz  

## 2. 元素循环FAPROTAX

    ## 方法1. 在线分析，推荐使用 http://www.ehbio.com/ImageGP 在线分析
    ## 方法2. Linux下分析、选学，详见附录3

## 3. Bugbase细菌表型预测


### 1. Bugbase命令行分析
    bugbase=${ea}/script/BugBase
    rm -rf result/bugbase/
    # 脚本已经优化适合R4.0，biom包更新为biomformat
    Rscript ${bugbase}/bin/run.bugbase.r -L ${bugbase} \
      -i result/gg/otutab.txt -m metadata.txt -c Group -o result/bugbase/

### 2. 其它可用分析
    # 使用 http://www.ehbio.com/ImageGP
    # 官网，https://bugbase.cs.umn.edu/ ，有报错，不推荐
    # Bugbase细菌表型预测Linux，详见附录4. Bugbase细菌表型预测



# 33、MachineLearning机器学习

    # RandomForest包使用的R代码见33MachineLearning目录中的RF_classification和RF_regression
    ## Silme2随机森林/Adaboost使用代码见33MachineLearning目录中的slime2，或附录5


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
    #统计筛选OTU表特征数量，总计80个
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
	  muscle -in otus.fa -out otus_aligned.fas

    ### 方法1. 利用IQ-TREE快速构建ML进化树，2m
    mkdir -p iqtree
    iqtree -s otus_aligned.fas \
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
    Rscript ${ea}/script/table2itol.R -a -c double -D plan1 -i OTUID -l Genus -t %s -w 0.5 annotation.txt
    # 生成注释文件中每列为单独一个文件

    ## 方案2. 生成丰度柱形图注释文件
    Rscript ${ea}/script/table2itol.R -a -d -c none -D plan2 -b Phylum -i OTUID -l Genus -t %s -w 0.5 annotation.txt

    ## 方案3. 生成热图注释文件
    Rscript ${ea}/script/table2itol.R -c keep -D plan3 -i OTUID -t %s otutab.txt

    ## 方案4. 将整数转化成因子生成注释文件
    Rscript ${ea}/script/table2itol.R -a -c factor -D plan4 -i OTUID -l Genus -t %s -w 0 annotation.txt

    # 返回工作目录
    cd ${wd}


# 附录：Linux服务器下分析(选学)

    #注：Windows下可能无法运行以下代码，推荐在Linux下conda安装相关程序

## 1. LEfSe分析

    mkdir -p ~/amplicon/24Compare/LEfSe
    cd ~/amplicon/24Compare/LEfSe
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

    #推荐使用 http://www.ehbio.com/ImageGP 在线分析
    #有Linux服务器用户可参考以下代码搭建本地流程
    cd ~/amplicon/result
    mkdir -p picrust

    # 安装picurst

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


## 3. FAPROTAXS元素循环

    # 设置工作目录
    wd=/mnt/c/amplicon/result/
    cd ${wd}
    # 设置脚本目录
    script=/mnt/c/public/script/

### 1. 软件安装

    # 注：软件已经下载至public/script目录，在qiime2环境下运行可满足依赖关系
    #(可选)下载软件1.2.3版， May 1, 2020更新数据库
    #wget -c https://pages.uoregon.edu/slouca/LoucaLab/archive/FAPROTAX/SECTION_Download/MODULE_Downloads/CLASS_Latest%20release/UNIT_FAPROTAX_1.2.3/FAPROTAX_1.2.3.zip
    #解压
    #unzip FAPROTAX_1.2.3.zip
    #(可选)依赖关系，可使用conda安装依赖包
    #conda install numpy
    #conda install biom

    #新建一个python3环境，或进入python3环境，如qiime2
    conda activate qiime2-2020.6
    #测试是否可运行，弹出帮助即正常工作
    python ${script}/FAPROTAX_1.2.3/collapse_table.py


### 2. 制作输入OTU表

    conda activate qiime2-2020.6
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
    python ${script}/FAPROTAX_1.2.3/collapse_table.py -i otutab_rare_tax.biom \
      -g ${script}/FAPROTAX_1.2.3/FAPROTAX.txt \
      --collapse_by_metadata 'taxonomy' -v --force \
      -o faprotax.txt -r faprotax_report.txt

### 4. 制作OTU对应功能注释有无矩阵

    # 对ASV(OTU)注释行，及前一行标题进行筛选
    grep 'ASV_' -B 1 faprotax_report.txt | grep -v -P '^--$' > faprotax_report.clean
    # faprotax_report_sum.pl脚本将数据整理为表格，位于public/scrit中
    perl ${script}/faprotax_report_sum.pl -i faprotax_report.clean -o faprotax_report
    # 查看功能有无矩阵，-S不换行
    less -S faprotax_report.mat


## 4. Bugbase细菌表型预测

### 1. 软件安装(仅一次)

    #有两种方法可选，推荐第一种，可选第二种，仅需运行一次

    #方法1. git下载，需要有git
    git clone https://github.com/knights-lab/BugBase

    #方法2. 下载并解压
    wget -c https://github.com/knights-lab/BugBase/archive/master.zip
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

## 5. Silme2随机森林/Adaboost

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

    # 使用实战
    cd 33MachineLearning/slime2
    #使用adaboost计算10000次(16.7s)，推荐千万次
    ./slime2.py otutab.txt design.txt --normalize --tag ab_e4 ab -n 10000
    #使用RandomForest计算10000次(14.5s)，推荐百万次，支持多线程
    ./slime2.py otutab.txt design.txt --normalize --tag rf_e4 rf -n 10000
    cd ../../




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

    head -n3 metadata.txt|cat -A
    sed -i 's/\r//' metadata.txt
    head -n3 metadata.txt|cat -A


## 11. UNITE数据库分析报错问题

    # Unprintable ASCII character no 195 on or right before line 236492
    cd /c/amplicon/FAQ/ITS_unite
    # 分类级存在空缺，sed补全
    sed -i 's/,;/,Unnamed;/;s/:,/:Unnamed,/g' ${ea}/usearch/utax_reference_dataset_all_04.02.2020.fasta
    # vsearch有问题，改为usearch，结尾添加--strand plus即可
    usearch --sintax  otus.fa \
      --db ${ea}/usearch/utax_reference_dataset_all_04.02.2020.fasta \
      --tabbedout otus.sintax --sintax_cutoff 0.6 --strand plus


# 版本更新记录

> v1.01 2018/03/10 建立以USEARCH为核心的扩增子分析流程，并采用RStudio 1.1新增Terimal结合GitForWindows为使用环境

> v1.02 2018/09/14 添加VSEARCH解决USEARCH免费版文件大小限制问题，新增Windows下R包移植方法

> v1.03 2019/01/11 绘图代码整理为函数，建立 amplicon包 https://github.com/microbiota/amplicon

> v1.04 2019/04/12 新增beta_cpcoa_dis函数，支持基于距离矩阵计算限制性PCoA

> v1.05 2019/07/19 更新format2lefse函数，修改默认分组Genotype为Group，支持设置输出文件名

> v1.06 2019/10/25 全面测试流程代码和R包，优化稳定性；陈同测试Mac系统下兼容性

> v1.07 2020/01/03 更新格式转换函数format2stamp/lefse/lefse2

> v1.08 2020/05/29 发表综述可供本流程引用Yong-Xin Liu, Yuan Qin, Tong Chen, Meiping Lu, Xubo Qian, Xiaoxuan Guo & Yang Bai. (2020). A practical guide to amplicon and metagenomic analysis of microbiome data. Protein & Cell 11, doi: https://doi.org/10.1007/s13238-020-00724-8，优化amplicon包中的多个统计绘图函数alpha_barplot/MicroTest

> v1.09 2020/08/14 优化amplicon包中的多个统计绘图函数alpha_boxplot/alpha_rare_all/alpha_sample_rare/alpha_rare_curve/BetaDiv/pairMicroTest/beta_pcoa/beta_cpcoa/beta_rda-cca/tax_circlize/tax_maptree/tax_stackplot

> v1.10 2020/12/04 Github公共发布1.10版，并发表中文使用手册于Bio-Protocol中Microbiome eBook，优化amplicon中多个函数beta_pcoa_stat/tax_stack_cluster/tax_wordcloud



