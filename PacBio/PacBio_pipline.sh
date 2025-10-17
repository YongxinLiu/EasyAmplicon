[TOC]

# EasyAmplicon 2 PacBio Pipeline分析流程 (Usearch/Vsearch)

    # Author作者: Yong-xin Liu(刘永鑫),,Tong Chen(陈同), Hao Luo(罗豪), Defeng Bai(白德凤), et al.
    # Update更新时间: 2025-10-17
    # Version版本: 2.0

    # Set the working directory (wd) and the software/database directory (db)
    # Most of the pipeline runs in linux bash
    # **每次打开Rstudio必须运行下面4行 Run it**，可选替换${db}为EasyMicrobiome安装位置
    # **The following 4 lines must be run every time you open RStudio**, you can optionally replace ${db} with your EasyMicrobiome installation location
    wd=/mnt/d/amplicon2/PacBio
    db=/mnt/d/EasyMicrobiome
    PATH=$PATH:${db}/script:${db}/linux
    mkdir -p $wd && cd ${wd}
    mkdir -p seq result temp 

## 1. Start files起始文件

    # 1. Analysis pipeline 分析流程: pipeline.sh
    # 2. Sample metadata 样本元信息: result/metadata.txt
    # 3. Sequencing data测序数据: `seq/*.fq.gz`结尾
    
### 1.1. Metadata/Design元数据/实验设计

    # Prepare the sample metadata file准备样本元数据: result/metadata.txt
    # csvtk统计表行(样本数，不含表头)列数，-t设置列分隔为制表符，默认为;
    # Use csvtk to count the rows (number of samples, excluding header) and columns of the table. -t sets the column separator to tab, default is comma.
    csvtk -t stat result/metadata.txt
    # 元数据至少3列，首列为样本ID(SampleID)，结尾列为描述(Description)
    # The metadata must have at least 3 columns: the first column is the sample ID (SampleID), and the last column is the description (Description).
    # cat查看文件，-A显示符号，"|"为管道符实现命令连用，head显示文件头，-n3控制范围前3行
    # Use 'cat -A' to view the file with non-printing characters. '|' is a pipe to chain commands. 'head -n3' displays the first 3 lines.
    cat -A result/metadata.txt | head -n3
    # windows用户结尾有^M，运行sed命令去除，再用cat -A检查
    # For Windows users, if there is a ^M at the end of the line, run the sed command to remove it, and then check with 'cat -A'.
    sed -i 's/\r//g' result/metadata.txt
    cat -A result/metadata.txt | head -n3

### 1.2. Sequencing data测序数据

    # 测序数据下载百度网盘链接：https://pan.baidu.com/s/1Ikd_47HHODOqC3Rcx6eJ6Q?pwd=0315；文件路径：db/pacbio/seq
    # Baidu Net Disk link for sequencing data download: https://pan.baidu.com/s/1Ikd_47HHODOqC3Rcx6eJ6Q?pwd=0315; File path: db/pacbio/seq
    # 公司返回的三代（Pacbio或者Nanopore）全长扩增子测序结果，通常为一个样品只有一个fq/fastq.gz格式压缩文件
    # The full-length amplicon sequencing results from third-generation sequencing (PacBio or Nanopore) usually consist of one compressed fq/fastq.gz file per sample.
    # 如果测序数据是.gz的压缩文件，有时需要使用gunzip解压后使用，vsearch通常可直接读取压缩文件
    # If the sequencing data is in a .gz compressed file, it sometimes needs to be decompressed with gunzip before use. Vsearch can usually read compressed files directly.
    gunzip seq/*.gz
    # zless按页查看压缩文件，空格翻页、q退出；head默认查看前10行，-n指定行
    # Use 'zless' to view compressed files page by page (space to scroll, q to quit); 'head' displays the first 10 lines by default, -n specifies the number of lines.
    ls -sh seq/
    head seq/HY1_1.fq
    # 每行太长，指定查看每行的1-60个字符
    # The lines are too long, so we view only the first 60 characters of each line.
    head seq/HY1_1.fq | cut -c 1-60
    # 统计测序数据，依赖seqkit程序
    # Use seqkit to get statistics of the sequencing data.
    seqkit stat seq/HY1_1.fq
    # 批量统计测序数据并汇总表
    # Batch process all sequencing files and summarize the statistics in a table.
    seqkit stat seq/*.fq > result/seqkit.txt
    head result/seqkit.txt

### 1.3. Pipeline & Database 流程和数据库

    # 数据库第一次使用必须解压，以后可跳过此段
    # The database must be decompressed before the first use. This step can be skipped later.

    # usearch可用的16S/18S/ITS数据库：
    # Available 16S/18S/ITS databases for usearch:

    # 解压Silva数据库，需自行从官网或百度网盘 (https://pan.baidu.com/s/1Ikd_47HHODOqC3Rcx6eJ6Q?pwd=0315; db/amplicon/silva) 下载
    # Decompress the Silva database. You need to download it from the official website or the Baidu Net Disk (https://pan.baidu.com/s/1Ikd_47HHODOqC3Rcx6eJ6Q?pwd=0315; db/amplicon/silva).

    # 将SILVA数据库SILVA_modified.fasta保存在${db}/usearch/
    # Save SILVA_modified.fasta to ${db}/usearch/

    # 将SILVA数据库silva_nr99_v138.1_train_DADA2.fa.gz保存在${db}/DADA2/
    # Save silva_nr99_v138.1_train_DADA2.fa.gz to ${db}/DADA2/

    # 将其他可选数据库从百度网盘 (https://pan.baidu.com/s/1Ikd_47HHODOqC3Rcx6eJ6Q?pwd=0315; db/amplicon/usearch)保存到${db}/usearch/目录，例如：
    # Save other optional databases to ${db}/usearch/, for example:
    #   sintax_defalut_emu_database.fasta.gz
    #   sintax_ncbi_database.fasta.gz
    #   gtdb_sintax_database.fasta.gz
    gunzip ${db}/usearch/sintax_defalut_emu_database.fasta.gz
    seqkit stat ${db}/usearch/sintax_defalut_emu_database.fasta    
    # Greengene数据库用于功能注释: ftp://greengenes.microbio.me/greengenes_release/gg_13_5/gg_13_8_otus.tar.gz
    # Greengenes database for functional annotation: ftp://greengenes.microbio.me/greengenes_release/gg_13_5/gg_13_8_otus.tar.gz
    # 默认解压会删除原文件，-c指定输出至屏幕，> 写入新文件(可改名)
    # Decompression deletes the original file by default. -c specifies output to the screen, > writes to a new file (can be renamed).
    gunzip -c ${db}/gg/97_otus.fasta.gz > ${db}/gg/97_otus.fa
    seqkit stat ${db}/gg/97_otus.fa

## 2. Reads rename and merge序列重命名和合并

    # Example of renaming a single sequence file 单个序列改名测试
    i=HY1_1
    usearch -fastx_relabel seq/${i}.fq -fastqout temp/${i}.rename.fq -prefix ${i}.
    head temp/${i}.rename.fq | cut -c1-60
    
    # Batch reads rename 批量序列重命名
    time for i in `tail -n+2 result/metadata.txt|cut -f1`;do
     usearch -fastx_relabel seq/${i}.fq -fastqout temp/${i}.rename.fq -prefix ${i}.
    done &
    # For vsearch method refer to "FAQ 2"

### 2.2 Merge renamed reads改名后序列整合

    cat temp/*.rename.fq > temp/all.fastq
    # File size查看文件大小353 M
    ls -lsh temp/all.fastq

## 3. Cut primers and quality control切除引物与质控

    conda activate easyamplicon2
    # Cutadapt去除接头："前向引物...反向引物", 12s
    # Cutadapt to remove adapters: "forward_primer...reverse_primer", ~12s
    # 引物错误率控制<10%，job任务数0为自动选择，去掉末匹配引物，剩余~243M
    time cutadapt -g "AGAGTTTGATCCTGGCTCAG...AAGTCSTAACAAGGTADCCSTA" \
      --error-rate=0.1 -j 0 \
      --discard-untrimmed \
      -o temp/allfilter.fastq temp/all.fastq
    ls -lsh temp/allfilter.fastq

    # 细菌16S片段长度通常约为1500bp，为避免过长或过短序列干扰，使用vsearch进行片段长度筛选。筛选长度在1200~1800 bp之间的序列，输出为fasta格式
    # The length of bacterial 16S fragments is usually around 1500 bp. To avoid interference from excessively long or short sequences, use vsearch for fragment length filtering. Filter sequences between 1200 and 1800 bp and output in FASTA format.
    vsearch --fastx_filter temp/allfilter.fastq --fastq_minlen 1200 --fastq_maxlen 1800 \
      --fastaout temp/filtered.fa --fastq_qmax 93
    head temp/filtered.fa

## 4. Dereplicate and cluster/denoise去冗余挑选OTU/ASV

### 4.1 Dereplicate sequences序列去冗余

    #去冗余：完全一致的序列合并，统计size信息，输出唯一序列，-minsize表示保留序列的最小条数,-sizeout输出丰度,--fasta_width 0输出FASTA时不换行（单行序列）减少文件体积 --relabel必须加序列前缀更规范, 1s
    # Dereplicate: Merge identical sequences, count size information, output unique sequences. -minsize specifies the minimum number of sequences to keep. -sizeout outputs abundance. --fasta_width 0 outputs FASTA without line breaks (single-line sequence) to reduce file size. --relabel must be added with a sequence prefix for better standardization. Takes 1s.
    time vsearch --derep_fulllength temp/filtered.fa --fasta_width 0 --sizeout --relabel Uni_  \
      --output temp/uniques.fa --minuniquesize 2 --threads 8

    # View查看uniques.fa, 4.1M
    ls -lsh temp/uniques.fa
    head temp/uniques.fa | cut -c 1-60
    # Uni_1;size=2148  - Unique sequence Uni_1 appears 165 times

### 4.2 Cluster聚类OTU/Denoise去噪ASV

    #有两种方法：推荐unoise3去噪获得单碱基精度ASV，备选传统的97%聚类OTU(属水平精度)
    # There are two methods: unoise3 denoising to obtain single-base precision ASVs is recommended. The alternative is traditional 97% clustering for OTUs (genus-level precision).
    # -minsize二次过滤，控制OTU/ASV数量至1-5千，方便下游统计分析
    # -minsize for secondary filtering, to control the number of OTUs/ASVs to 1-5 thousand, which is convenient for downstream statistical analysis.

    # Method方法1. ASV 去噪 Denoise: predict biological sequences and filter chimeras
    # 1s, 331 good, 18 chimeras
    time usearch -unoise3 temp/uniques.fa -minsize 2 \
      -zotus temp/zotus.fa
      
    # Method2. 97% OTU clustering 聚类
    # 3s, produced 264 OTUs, and removed 4 chimeras.
    # time usearch -cluster_otus temp/uniques.fa -minsize 2 \
    #   -otus temp/otus.fa \
    #   -relabel OTU_
 
    #修改序列名：Zotu改为ASV方便识别
    # Rename sequences: Change Zotu to ASV for easy identification.
    sed 's/Zotu/ASV_/g' temp/zotus.fa > temp/otus.fa
    head -n 2 temp/otus.fa

    #方法3. 数据过大无法使用usearch时，备选vsearch方法见"常见问题3"
    # Method 3. When the data is too large to use usearch, see "FAQ 3" for the alternative vsearch method.

### 4.3 Reference-based chimera detection 基于参考去嵌合

    mkdir -p result/raw

    # Option 1. Chimera removal with vsearch + SILVA
    vsearch --uchime_ref temp/otus.fa \
      -db ${db}/usearch/SILVA_modified.fasta \
      --nonchimeras result/raw/otus.fa
    # 1m, 19 chimeras

    # Option 2. Do not remove chimeras不去嵌合
    cp -f temp/otus.fa result/raw/otus.fa

## 5. Feature table 特征表构建和筛选

    # OTU和ASV统称为特征(Feature)，它们的区别是：
    # OTUs and ASVs are collectively referred to as Features. Their differences are:
    # OTU通常按97%聚类后挑选最高丰度或中心的代表性序列；
    # OTUs are usually representative sequences selected with the highest abundance or from the center after 97% clustering.
    # ASV是基于序列进行去噪(排除或校正错误序列，并挑选丰度较高的可信序列)作为代表性序列
    # ASVs are representative sequences obtained by denoising based on sequences (excluding or correcting erroneous sequences and selecting credible sequences with higher abundance).

### 5.1 Generate a feature table 生成特征表

    # vsearch generates a feature table
    # id(1): 100% similarity alignment: 22148 of 86454 (25.62%), real time 8m16s
    # time vsearch --usearch_global temp/filtered.fa \
    #   --db result/raw/otus.fa \
    #   --id 1 --threads 4 \
    # 	--otutabout result/raw/otutab.txt
    	
    # id(0.97)推荐：(更高数据使用率，更快)82408 of 86454 (95.32%), time 2m53s
    time vsearch --usearch_global temp/filtered.fa \
      --db result/raw/otus.fa \
      --id 0.97 --threads 4 \
    	--otutabout result/raw/otutab.txt
    head -n6 result/raw/otutab.txt | cut -f 1-6
    # csvtk count the rows and columns 统计表行列
    # 这里一定看好列数，是不是等于你的样品数；如果不等，一般是样品命名存在问题，具体看上面解释
    # Be sure to check the number of columns here. Is it equal to your number of samples? If not, there is generally a problem with the sample naming. See the explanation above for details.
    csvtk -t stat result/raw/otutab.txt

### 5.2 Taxonomic annotation, and/or removal of plastids and non-Bacteria 物种注释，且/或去除质体和非细菌

    # 方法1. vsearch sintax注释
    # Method1. sintax annotation
    # Taxonomic annotation - remove plastids and non-bacteria/archaea and calculate the proportion (optional).
    # SILVA数据库(SILVA_138.2_SSURef_NR99_tax_silva.fasta)更好注释真核、质体序列
    # The SILVA database (SILVA_138.2_SSURef_NR99_tax_silva.fasta) is better for annotating eukaryotic and plastid sequences.
    # 修改SILVA数据库(SILVA_138.2_SSURef_NR99_tax_silva.fasta)格式以适应注释代码（百度网盘下载已转换好的数据库文件SILVA_modified.fasta）
    # Modify the format of the SILVA database (SILVA_138.2_SSURef_NR99_tax_silva.fasta) to adapt to the annotation code (download the converted database file SILVA_modified.fasta from Baidu Net Disk).
    # python database_silva.py 输出数据库文件SILVA_modified.fasta，包含种和株系的信息
    # python database_silva.py outputs the database file SILVA_modified.fasta, which contains species and strain information.
    # sintax_defalut_emu_database.fasta is from emu pipeline
    # 置信阈值通常0.6/0.8，vserch最低0.1/usearch可选0输出最相似物种注释用于观察潜在分类,1m58s
    # The confidence threshold is usually 0.6/0.8. The minimum for vsearch is 0.1. usearch can optionally be set to 0 to output the most similar species annotation for observing potential classifications. Takes 1m 58s.
    vsearch --sintax result/raw/otus.fa \
      --db ${db}/usearch/sintax_defalut_emu_database.fasta \
      --sintax_cutoff 0.1 \
      --tabbedout result/raw/otus.sintax 
    head result/raw/otus.sintax 

    # Number of rows in the original feature table, 272
    wc -l result/raw/otutab.txt
    # R脚本选择细菌古菌(真核)、去除叶绿体、线粒体并统计比例；输出筛选并排序的OTU表
    # R script to select bacteria/archaea (eukaryotes), remove chloroplasts, mitochondria, and calculate proportions; output a filtered and sorted OTU table.
    # 输入为OTU表result/raw/otutab.txt和物种注释result/raw/otus.sintax
    # The input is the OTU table result/raw/otutab.txt and the taxonomic annotation result/raw/otus.sintax.
    # 输出筛选并排序的特征表result/otutab.txt和
    # The output is the filtered and sorted feature table result/otutab.txt and
    # 统计污染比例文件result/raw/otutab_nonBac.txt和过滤细节otus.sintax.discard
    # the contamination proportion file result/raw/otutab_nonBac.txt and the filtering details otus.sintax.discard.
    # 真菌ITS数据，请改用otutab_filter_nonFungi.R脚本，只筛选真菌
    # For fungal ITS data, please use the otutab_filter_nonFungi.R script to filter only fungi.
    # linux下运行不成功，可切换到windows git bash中运行 db=/d/EasyMicrobiome
    Rscript ${db}/script/otutab_filter_nonBac.R -h # 显示参数说明 / Display parameter description
    Rscript ${db}/script/otutab_filter_nonBac.R \
      --input result/raw/otutab.txt \
      --taxonomy result/raw/otus.sintax \
      --output result/otutab.txt\
      --stat result/raw/otutab_nonBac.stat \
      --discard result/raw/otus.sintax.discard
    # Number of rows in the feature table after filtering 筛选后特征表行数, 272
    wc -l result/otutab.txt
    #过滤特征表对应序列
    # Filter the corresponding sequences in the feature table
    cut -f 1 result/otutab.txt | tail -n+2 > result/otutab.id
    usearch -fastx_getseqs result/raw/otus.fa \
        -labels result/otutab.id -fastaout result/otus.fa
    #过滤特征表对应序列注释
    # Filter the corresponding sequence annotations in the feature table
    awk 'NR==FNR{a[$1]=$0}NR>FNR{print a[$1]}'\
        result/raw/otus.sintax result/otutab.id \
        > result/otus.sintax

    # 可选. 怀疑筛选不合理跳过
    # 可选. If you think the filtering is unreasonable, you can skip it.
    # cp result/raw/otu* result/

    # Summary OTU table
    usearch -otutab_stats result/otutab.txt \
      -output result/otutab.stat
    cat result/otutab.stat
    #注意最小值、分位数，或查看result/raw/otutab_nonBac.stat中样本详细数据量，用于重采样
    # Pay attention to the minimum value, quantiles, or view the detailed sample data volume in result/raw/otutab_nonBac.stat for resampling.

    # Method方法2. emu annotation注释
    # cutadapt by each sample
    for fq in temp/*.rename.fq; do
      sample=$(basename "$fq" .rename.fq)
      cutadapt -g "AGAGTTTGATCCTGGCTCAG...AAGTCSTAACAAGGTADCCSTA" \
        --error-rate=0.1 --action=trim --rc \
        -m 1000 -M 1800 -j 8 \
        -o "temp/${sample}.filtered.fastq" "$fq";
    done
    
    # emu quantify (可选)
    # for fastq in temp/*.filtered.fastq; do
    # sample=$(basename "$fastq" .filtered.fastq)
    # # Create a directory for each sample 为每个样本创建输出目录
    # mkdir -p "result/Emu/silva/$sample"
    # # emu abundance丰度估计
    # emu abundance "$fastq" \
    # --type map-ont \
    # --db /mnt/c/EasyMicrobiome/Silva_Emu \
    # --output-dir "result/Emu/silva/$sample" \
    # --threads 4
    # done

### 5.3 等量抽样标准化
### 5.3 Normalization by subsampling

    # Normlize by subsample
    # 通过子抽样进行标准化

    #使用vegan包进行等量重抽样，输入reads count格式Feature表result/otutab.txt
    # Use the vegan package for equal resampling. The input is the feature table in reads count format: result/otutab.txt.
    #可指定输入文件、抽样量和随机数，输出抽平表result/otutab_rare.txt和多样性alpha/vegan.txt
    # You can specify the input file, sampling depth, and random seed. The output is the rarefied table result/otutab_rare.txt and the diversity file alpha/vegan.txt.
    mkdir -p result/alpha
    Rscript ${db}/script/otutab_rare.R --input result/otutab.txt \
      --depth 7333 --seed 1 \
      --normalize result/otutab_rare.txt \
      --output result/alpha/vegan.txt
    usearch -otutab_stats result/otutab_rare.txt \
      -output result/otutab_rare.stat
    cat result/otutab_rare.stat

## 6. α多样性
## 6. Alpha diversity

### 6.1. 计算α多样性
### 6.1. Calculate alpha diversity

    # 使用USEARCH计算14种alpha多样性指数(Chao1有错勿用)
    # Use USEARCH to calculate 14 alpha diversity indices (Chao1 has errors, do not use).
    #details in http://www.drive5.com/usearch/manual/alpha_metrics.html
    usearch -alpha_div result/otutab_rare.txt \
      -output result/alpha/alpha.txt

### 6.2. 计算稀释丰富度
### 6.2. Calculate rarefaction richness

    #稀释曲线：取1%-100%的序列中OTUs数量，每次无放回抽样
    # Rarefaction curve: Take the number of OTUs in 1%-100% of the sequences, sampling without replacement each time.
    #Rarefaction from 1%, 2% .. 100% in richness (observed OTUs)-method without_replacement https://drive5.com/usearch/manual/cmd_otutab_subsample.html
    usearch -alpha_div_rare result/otutab_rare.txt \
      -output result/alpha/alpha_rare.txt \
      -method without_replacement
    #预览结果
    # Preview the results
    head -n2 result/alpha/alpha_rare.txt
    #样本测序量低出现非数值"-"的处理，详见常见问题8
    # For samples with low sequencing depth, non-numeric values "-" may appear. See FAQ 8 for details on how to handle this.
    sed -i "s/-/\t0.0/g" result/alpha/alpha_rare.txt

### 6.3. 筛选高丰度菌
### 6.3. Filter by abundance

    #计算各特征的均值，有组再求分组均值，需根据实验设计metadata.txt修改组列名
    # Calculate the mean of each feature. If there are groups, calculate the group means. You need to modify the group column name according to the experimental design metadata.txt.
    #输入文件为feautre表result/otutab.txt，实验设计metadata.txt
    # The input files are the feature table result/otutab.txt and the experimental design metadata.txt.
    #输出为特征表按组的均值-一个实验可能有多种分组方式
    # The output is the mean of the feature table by group. An experiment may have multiple grouping methods.
    #-h显示脚本帮助(参数说明)
    # -h displays the script help (parameter description).
    Rscript ${db}/script/otu_mean.R -h
    #scale是否标准化，zoom标准化总和，all输出全部样本均值，type计算类型mean或sum
    # scale: whether to standardize; zoom: standardize the sum; all: output the mean of all samples; type: calculation type, mean or sum.
    Rscript ${db}/script/otu_mean.R --input result/otutab.txt \
      --metadata result/metadata.txt \
      --group Group --thre 0 \
      --scale TRUE --zoom 100 --all TRUE --type mean \
      --output result/otutab_mean.txt
    # 结果为全部和各组均值
    # The result is the mean of all samples and each group.
    head -n3 result/otutab_mean.txt

    #如以平均丰度>0.1%筛选，可选0.5或0.05，得到每个组的OTU组合，i=3为从第三列开始，要保留all列的话改为i=2,这里为了进行四组比较的韦恩图绘制保留了all列
    # For example, filter by an average abundance of >0.1% (you can choose 0.5 or 0.05) to get the OTU combination for each group. i=3 means starting from the third column. If you want to keep the 'all' column, change it to i=2. Here, the 'all' column is kept for drawing a Venn diagram for four-group comparison.
    awk 'BEGIN{OFS=FS="\t"}{if(FNR==1) {for(i=2;i<=NF;i++) a[i]=$i; print "OTU","Group";} \
        else {for(i=2;i<=NF;i++) if($i>0.1) print $1, a[i];}}' \
        result/otutab_mean.txt > result/alpha/otu_group_exist.txt
    head result/alpha/otu_group_exist.txt
    cut -f 2 result/alpha/otu_group_exist.txt | sort | uniq -c
    # 试一试：不同丰度下各组有多少OTU/ASV
    # Try it: How many OTUs/ASVs are there in each group at different abundances?
    # 可在 http://ehbio.com/test/venn/ 中绘图并显示各组共有和特有维恩或网络图
    # You can draw and display the common and unique Venn or network diagrams of each group at http://ehbio.com/test/venn/.
    # 也可在 http://www.ehbio.com/ImageGP 绘制Venn、upSetView和Sanky
    # You can also draw Venn, upSetView, and Sankey diagrams at http://www.ehbio.com/ImageGP.

## 7. β多样性
## 7. Beta diversity

    #结果有多个文件，需要目录
    # The results have multiple files, so a directory is needed.
    mkdir -p result/beta/
    #基于OTU构建进化树 Make OTU tree, 6.7s
    # Build a phylogenetic tree based on OTUs. Takes 6.7s.
    usearch -cluster_agg result/otus.fa -treeout result/otus.tree
    #1s生成5种距离矩阵：bray_curtis, euclidean, jaccard, manhatten, unifrac
    # Generate 5 distance matrices in 1s: bray_curtis, euclidean, jaccard, manhattan, unifrac.
    usearch -beta_div result/otutab_rare.txt -tree result/otus.tree \
    -filename_prefix result/beta/

## 8. 物种注释分类汇总
## 8. Summary of taxonomic annotation

    #OTU对应物种注释2列格式：去除sintax中置信值，只保留物种注释，替换:为_，删除引号
    # 2-column format for OTU corresponding species annotation: remove the confidence value in sintax, keep only the species annotation, replace ":" with "_", and delete quotation marks.
    cut -f 1,4 result/otus.sintax \
      |sed 's/\td/\tk/;s/:/__/g;s/,/;/g;s/"//g' \
      > result/taxonomy2.txt
    head -n3 result/taxonomy2.txt

    #OTU对应物种8列格式：注意注释是非整齐
    # 8-column format for OTU corresponding species: Note that the annotation is not uniform.
    #生成物种表格OTU/ASV中空白补齐为Unassigned
    # In the generated species table, fill in the blanks in OTU/ASV with "Unassigned".
    awk 'BEGIN{OFS=FS="\t"}{delete a; a["k"]="Unassigned";a["p"]="Unassigned";a["c"]="Unassigned";a["o"]="Unassigned";a["f"]="Unassigned";a["g"]="Unassigned";a["s"]="Unassigned";\
      split($2,x,";");for(i in x){split(x[i],b,"__");a[b[1]]=b[2];} \
      print $1,a["k"],a["p"],a["c"],a["o"],a["f"],a["g"],a["s"];}' \
      result/taxonomy2.txt > temp/otus.tax
    sed 's/;/\t/g;s/.__//g;' temp/otus.tax|cut -f 1-8 | \
      sed '1 s/^/OTUID\tKingdom\tPhylum\tClass\tOrder\tFamily\tGenus\tSpecies\n/' \
      > result/taxonomy.txt
    head -n3 result/taxonomy.txt

    #统计门纲目科属，使用 rank参数 p c o f g s，为phylum, class, order, family, genus, species缩写
    # Count phylum, class, order, family, genus, species. Use the rank parameter p, c, o, f, g, s, which are abbreviations for phylum, class, order, family, genus, species.
    mkdir -p result/tax
    for i in p c o f g s;do
      usearch -sintax_summary result/otus.sintax \
      -otutabin result/otutab_rare.txt -rank ${i} \
      -output result/tax/sum_${i}.txt
    done
    sed -i 's/(//g;s/)//g;s/"//g;s/\#//g;s/\/Chloroplast//g' result/tax/sum_*.txt
    # 列出所有文件
    # List all files
    wc -l result/tax/sum_*.txt
    head -n3 result/tax/sum_s.txt

## 9. 有参定量特征表
## 9. Reference-based quantitative feature table

    # 比对Greengenes97% OTUs比对，用于PICRUSt/Bugbase功能预测
    # Align to Greengenes 97% OTUs for PICRUSt/Bugbase functional prediction.
    mkdir -p result/gg/

    #方法1. usearch比对更快，但文件超限报错选方法2
    # Method 1. usearch alignment is faster, but if the file size exceeds the limit and an error is reported, choose method 2.
    # 默认10核以下使用1核，10核以上使用10核
    # By default, 1 core is used for less than 10 cores, and 10 cores are used for more than 10 cores.
    usearch -otutab temp/filtered.fa -otus ${db}/gg/97_otus.fa \
    	-otutabout result/gg/otutab.txt -threads 4
    # 比对率80.8%, 4核1m, 内存使用731Mb
    # Alignment rate 80.8%, 4 cores, 1m, memory usage 731Mb.
    head -n3 result/gg/otutab.txt

    # #方法2. vsearch比对，更准更慢，但并行24-96线程更强
    # # Method 2. vsearch alignment is more accurate but slower, but it is more powerful with 24-96 parallel threads.
    # vsearch --usearch_global temp/filtered.fa --db ${db}/gg/97_otus.fa \
    #    --otutabout result/gg/otutab.txt --id 0.97 --threads 12
    # 比对率83.83%, 12核4m24s
    # Alignment rate 83.83%, 12 cores, 4m 24s.

    #统计
    # Statistics
    usearch -otutab_stats result/gg/otutab.txt -output result/gg/otutab.stat
    cat result/gg/otutab.stat

## 10. 空间清理及数据提交
## 10. Workspace cleanup and data submission

    #删除中间大文件
    # Delete large intermediate files
    rm -rf temp/*.fq

    # 分双端统计md5值，用于数据提交
    # Calculate md5 values for both ends for data submission
    cd seq
    md5sum *.fastq > ../result/md5sum.txt
    cat ../result/md5sum.txt

```
# 以上是基于Vsearh、Usearch的流程，接下来我们还提供利用DADA2 R软件包处理PacBio扩增子reads获得精准的ASVs，可以得到物种及其多样性信息
# The above is the pipeline based on Vsearh and Uusearch. Next, we also provide a method to process PacBio amplicon reads using the DADA2 R package to obtain accurate ASVs, which can provide species and diversity information.
# DADA2使用朴素贝叶斯分类算法结合数据库注释ASVs，可能会出现假阳性
# DADA2 uses a naive Bayesian classification algorithm combined with a database to annotate ASVs, which may produce false positives.
# 用户可根据自己的需求选择
# Users can choose according to their own needs.
```
# EasyAmplicon 2 PacBio Pipeline (DADA2)
# EasyAmplicon 2 PacBio 分析流程 (DADA2)
```
  # 在您第一次运行此流程之前，必须确保已经安装了DADA2及其他脚本依赖的R包。这个步骤只需要做一次。
  # Before you run this pipeline for the first time, you must ensure that DADA2 and other R packages that the script depends on have been installed. This step only needs to be done once.
  # DADA2包托管在Bioconductor上，需要通过BiocManager来安装。请打开您的R或Rstudio，在控制台(Console)中输入并执行以下命令：
  # The DADA2 package is hosted on Bioconductor and needs to be installed via BiocManager. Please open your R or Rstudio and enter and execute the following commands in the Console:
  # 首先，安装Bioconductor的核心管理工具 BiocManager
  # First, install Bioconductor's core management tool, BiocManager
  if (!requireNamespace("BiocManager", quietly = TRUE))
      install.packages("BiocManager")
  # 然后，通过 BiocManager 安装 DADA2
  # Then, install DADA2 via BiocManager
  BiocManager::install("dada2")
  # 此外，该脚本还需要 argparse 包来解析命令行参数，通过常规方式安装
  # In addition, this script also requires the argparse package to parse command-line arguments, which can be installed in the conventional way.
  install.packages("argparse")
  # 注意： 安装DADA2时，控制台可能会出现大量的编译和安装信息，请耐心等待其完成。如果遇到问题，请根据提示信息检查您的R语言环境或网络连接。
  # Note: When installing DADA2, a large amount of compilation and installation information may appear in the console. Please wait patiently for it to complete. If you encounter problems, please check your R language environment or network connection according to the prompt.
```
     # Define paths for input, output, metadata, and the R script
       # INPUT_DIR: Directory containing the raw FASTQ sequence files.
       # OUTPUT_DIR: Directory where the DADA2 analysis results will be saved.
       # METADATA: Path to the metadata file, which should contain sample information.
       # SCRIPT_PATH: Path to the R script that will be executed.
       # TAXONOMY_DB: Path to the taxonomy database file.
     # 定义输入、输出、元数据和R脚本的路径
       # INPUT_DIR: 包含原始FASTQ序列文件的目录。
       # OUTPUT_DIR: 用于保存DADA2分析结果的目录。
       # METADATA: 元数据文件的路径，应包含样本信息。
       # SCRIPT_PATH: 将要执行的R脚本的路径。
       # TAXONOMY_DB: 物种注释数据库文件的路径。
       INPUT_DIR="/mnt/d/EasyAmplicon2/PacBio/seq"
       OUTPUT_DIR="/mnt/d/EasyAmplicon2/PacBio/result/pacbio_dada2_output"
       METADATA="/mnt/d//EasyAmplicon2/PacBio/result/metadata.txt"
       SCRIPT_PATH="${db}/script/DADA2_PB.R"
       TAXONOMY_DB="${db}/DADA2/silva_nr99_v138.1_train_DADA2.fa.gz"
       # Run the R script using Rscript
       # This command executes the DADA2_PB.R script with several command-line arguments.
       # --input_dir: Specifies the directory with the input FASTQ files.
       # --output_dir: Specifies the directory to store the output results.
       # --metadata_file: Provides the path to the metadata file.
       # --taxonomy_db: Provides the path to the taxonomy database.
       # --threads: Sets the number of threads for parallel processing to improve performance.
       # 使用 Rscript 运行 R 脚本
       # 该命令使用多个命令行参数执行 DADA2_PB.R 脚本。
       # --input_dir: 指定输入FASTQ文件所在的目录。
       # --output_dir: 指定存放输出结果的目录。
       # --metadata_file: 提供元数据文件的路径。
       # --taxonomy_db: 提供物种注释数据库的路径。
       # --threads: 设置用于并行处理的线程数，以提高性能。
       Rscript ${db}/script/DADA2_PB.R \
           --input_dir ${INPUT_DIR} \
           --output_dir ${OUTPUT_DIR} \
           --metadata_file ${METADATA} \
           --taxonomy_db ${TAXONOMY_DB} \
           --threads 8

       # --- DADA2_PB.R --help information ---
       # --- DADA2_PB.R 的 --help 帮助信息 ---
       # usage: DADA2_PB.R [-h] --input_dir INPUT_DIR --output_dir OUTPUT_DIR
       #                   --metadata_file METADATA_FILE [--fwd_primer FWD_PRIMER]
       #                   [--rev_primer REV_PRIMER] [--min_len MIN_LEN]
       #                   [--max_len MAX_LEN] [--max_ee MAX_EE] [--threads THREADS]
       #
       # DADA2 pipeline for PacBio data. (用于PacBio数据的DADA2流程)
       #
       # options:
       #   -h, --help            show this help message and exit
       #   --input_dir INPUT_DIR
       #                         Path to the directory containing input FASTQ files.
          (输入FASTQ文件所在目录的路径)
       #   --output_dir OUTPUT_DIR
       #                         Path to the directory where results will be saved.
          (用于保存结果的目录路径)
       #   --metadata_file METADATA_FILE
       #                         Path to the metadata file. Must contain a 'SampleID' column.
          (元数据文件路径，必须包含'SampleID'列)
       #   --fwd_primer FWD_PRIMER
       #                         Forward primer sequence (e.g., 27F). (正向引物序列，例如27F)
       #   --rev_primer REV_PRIMER
       #                         Reverse primer sequence (e.g., 1492R). (反向引物序列，例如1492R)
       #   --min_len MIN_LEN     Minimum length of reads to keep after filtering.
          (过滤后保留序列的最小长度)
       #   --max_len MAX_LEN     Maximum length of reads to keep after filtering.
          (过滤后保留序列的最大长度)
       #   --max_ee MAX_EE       Maximum expected errors allowed. (允许的最大期望误差)
       #   --threads THREADS     Number of threads for multithreading. (用于多线程的线程数)
       #   --taxonomy_db TAXONOMY_DB    Path to the taxonomy database file (e.g., SILVA, GTDB). (物种注释数据库文件路径，例如SILVA、GTDB)
# PICRUSt2功能预测
# PICRUSt2 Functional Prediction

    # (可选)PICRUSt2(Linux/Windows下Linux子系统，要求>16GB内存)
    # (Optional) PICRUSt2 (Linux/Linux subsystem on Windows, requires >16GB memory)
    # 方法1. 直接安装
    # Method 1. Direct installation
    n=picrust2
    conda create -n ${n} -c bioconda -c conda-forge ${n}=2.3.0_b
    # 方法2. 导入安装环境(推荐)
    # Method 2. Import installation environment (recommended)
    n=picrust2
    # 复制安装包，或下载我的环境打包
    # Copy the installation package, or download my packaged environment
    wget -c ftp://download.nmdc.cn/tools/conda/${n}.tar.gz
    # 指定安装目录并解压
    # Specify the installation directory and decompress
    condapath=~/miniconda3
    mkdir -p ${condapath}/envs/${n}
    tar -xvzf ${n}.tar.gz -C ${condapath}/envs/${n}
    # 激活环境并初始化
    # Activate the environment and initialize
    conda activate picrust2
    conda unpack
    # 加载环境
    # Load the environment
    conda activate picrust2
    # 进入工作目录，服务器要修改工作目录
    # Enter the working directory. The working directory needs to be modified on the server.
    mkdir -p /mnt/c/EasyAmplicon2/PacBio/result/picrust2
    cd /mnt/c/EasyAmplicon2/PacBio/result/picrust2
    # 运行流程，内存15.7GB，耗时12m；内存不足会导致程序中断，out目录里只有一个intermediate/文件夹，需使用大内存设备运行
    # Run the pipeline. Memory: 15.7GB, time: 12m. Insufficient memory will cause the program to be interrupted. The out directory will only have an intermediate/ folder. A device with large memory is required.
    picrust2_pipeline.py -s ../otus.fa -i ../otutab.txt -o ./out -p 8
    # 添加EC/KO/Pathway注释
    # Add EC/KO/Pathway annotations
    cd out
    add_descriptions.py -i pathways_out/path_abun_unstrat.tsv.gz -m METACYC \
    -o pathways_out/path_abun_unstrat_descrip.tsv.gz
    add_descriptions.py -i EC_metagenome_out/pred_metagenome_unstrat.tsv.gz -m EC \
      -o EC_metagenome_out/pred_metagenome_unstrat_descrip.tsv.gz
    add_descriptions.py -i KO_metagenome_out/pred_metagenome_unstrat.tsv.gz -m KO \
      -o KO_metagenome_out/pred_metagenome_unstrat_descrip.tsv.gz 
    # KEGG按层级合并
    # Merge KEGG by hierarchy
    db=/mnt/c/EasyMicrobiome/
    zcat KO_metagenome_out/pred_metagenome_unstrat.tsv.gz > KEGG.KO.txt
    python3 ${db}/script/summarizeAbundance.py \
      -i KEGG.KO.txt \
        -m ${db}/kegg/KO1-4.txt \
        -c 2,3,4 -s ',+,+,' -n raw \
        -o KEGG
    # 统计各层级特征数量
    # Count the number of features at each level
    wc -l KEGG*

# Evolution进化树
# Evolution Tree

    #切换gitbash
    # Switch to gitbash
    cd ${wd}/PacBio
    mkdir -p result/tree
    cd ${wd}/PacBio/result/tree

## 筛选高丰度/指定的特征
## Filter high-abundance/specified features

    #方法1. 按丰度筛选特征，一般选0.001或0.005，且OTU数量在30-150个范围内
    # Method 1. Filter features by abundance, usually 0.001 or 0.005, and the number of OTUs is in the range of 30-150.
    #统计特征表中ASV数量，如总计1234个
    # Count the number of ASVs in the feature table, e.g., a total of 1234.
    tail -n+2 ../otutab_rare.txt | wc -l
    #按相对丰度0.2%筛选高丰度OTU
    # Filter high-abundance OTUs by a relative abundance of 0.2%.
    usearch -otutab_trim ../otutab_rare.txt \
        -min_otu_freq 0.002 \
        -output otutab.txt
    #统计筛选OTU表特征数量，总计~105个
    # Count the number of features in the filtered OTU table, a total of ~105.
    tail -n+2 otutab.txt | wc -l

    #方法2. 按数量筛选
    # Method 2. Filter by quantity
    # #按丰度排序，默认由大到小
    # # Sort by abundance, from largest to smallest by default.
    # usearch -otutab_sortotus ../otutab_rare.txt  \
    #     -output otutab_sort.txt
    # #提取高丰度中指定Top数量的OTU ID，如Top100,
    # # Extract the specified number of top OTU IDs from the high-abundance OTUs, e.g., Top100.
    # sed '1 s/#OTU ID/OTUID/' otutab_sort.txt \
    #     | head -n101 > otutab.txt

    #修改特征ID列名
    # Modify the feature ID column name
    sed -i '1 s/#OTU ID/OTUID/' otutab.txt
    #提取ID用于提取序列
    # Extract IDs for sequence extraction
    cut -f 1 otutab.txt > otutab_high.id

    # 筛选高丰度菌/指定差异菌对应OTU序列
    # Filter the OTU sequences corresponding to high-abundance bacteria/specified differential bacteria.
    usearch -fastx_getseqs ../otus.fa -labels otutab_high.id \
        -fastaout otus.fa
    head -n 2 otus.fa

    ## 筛选OTU对物种注释
    ## Filter OTU for species annotation
    awk 'NR==FNR{a[$1]=$0} NR>FNR{print a[$1]}' ../taxonomy.txt \
        otutab_high.id > otutab_high.tax

    #获得OTU对应组均值，用于样本热图
    # Get the mean of the OTU corresponding group for the sample heatmap.
    #依赖之前otu_mean.R计算过按Group分组的均值
    # Depends on the mean calculated by group in the previous otu_mean.R.
    awk 'NR==FNR{a[$1]=$0} NR>FNR{print a[$1]}' ../otutab_mean.txt otutab_high.id \
        | sed 's/#OTU ID/OTUID/' > otutab_high.mean
    head -n3 otutab_high.mean

    #合并物种注释和丰度为注释文件
    # Merge species annotation and abundance into an annotation file.
    cut -f 2- otutab_high.mean > temp
    paste otutab_high.tax temp > annotation.txt
    head -n 3 annotation.txt

## 构建进化树
## Build a phylogenetic tree

    # 起始文件为 result/tree目录中 otus.fa(序列)、annotation.txt(物种和相对丰度)文件
    # The starting files are otus.fa (sequences) and annotation.txt (species and relative abundance) in the result/tree directory.
    # Muscle软件进行序列对齐，3s
    # Use Muscle for sequence alignment, takes 3s.
    muscle -in otus.fa -out otus_aligned.fas

    ### 方法1. 利用IQ-TREE快速构建ML进化树，2m
    ### Method 1. Use IQ-TREE to quickly build an ML phylogenetic tree, takes 2m.
    # rm -rf iqtree
    # mkdir -p iqtree
    # iqtree -s otus_aligned.fas \
    # -bb 1000 -redo -alrt 1000 -nt AUTO \
    # -pre iqtree/otus

    ### 方法2. FastTree快速建树(Linux)
    ### Method 2. FastTree for fast tree building (Linux)
    # 注意FastTree软件输入文件为fasta格式的文件，而不是通常用的Phylip格式。输出文件是Newick格式。
    # Note that the input file for FastTree is in FASTA format, not the commonly used Phylip format. The output file is in Newick format.
    # 该方法适合于大数据，例如几百个OTUs的系统发育树！
    # This method is suitable for big data, such as phylogenetic trees with hundreds of OTUs!
    # Ubuntu上安装fasttree可以使用`sudo apt install fasttree`
    # On Ubuntu, you can use `sudo apt install fasttree` to install fasttree.
    cd result/tree
    fasttree -gtr -nt otus_aligned.fas > otus.nwk