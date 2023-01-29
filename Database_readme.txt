Instructions for use:

1. For the first time, please decompress the data files you need to use;
2. If you need to use it frequently, it is recommended not to compress. Each compression and decompression will waste a lot of time; compression is only to save time during transmission;
3. Before analysis, please visit the official database of each database to check whether there is an update of the database version;

1. GreenGene database (gg)

FTP address: ftp://greengenes.microbio.me/greengenes_release

The latest version of the QIIME format database download link: 13_8 represents the update in August 2013

ftp://greengenes.microbio.me/greengenes_release/gg_13_5/gg_13_8_otus.tar.gz # 305 MB, updated on 2013/8/15, the sequence is the same as 13_5, only some species annotations are corrected

After decompression, there are only three files in the 97 directory:

97_otu_taxonomy.txt.gz # Species annotation of 97% cluster representative sequences, 1.2 MB
97_otus.fasta.gz # 97% cluster representative sequence, 29 MB
97_otus.tree.gz # Multiple sequence tree file of 97% cluster representative sequences, 929 KB


2. USEARCH format database (usearch)

Download page: http://www.drive5.com/sintax

16S

rdp_16s_v18.fa.gz # 20,000 sequences for genus annotation, newly version

rdp_16s_v16_sp.fa.gz # 13,000 sequences for species annotation,

silva_16s_v123.fa.gz # 1.6 million for removing chimeras to prevent false negatives. The more complete the chimera, the better


3. UNITE database - ITS data

Download page: https://unite.ut.ee/

utax_reference_dataset_all_29.11.2022.fasta.gz # UNITE newly version fungal sequences and annotations

The USEARCH data is a secondary arrangement and may not be updated in time. You can download the latest and full versions from the UNITE official website.

If there is a problem, please refer to Appendix ## 11. UNITE database analysis error report in Pipeline.sh

使用说明：

1. 第一次使用请先解压缩需要用到的数据文件；
2. 如果需要经常使用，建议不要压缩，每次压缩和解压会浪费很多时间；压缩只为在传输时节约时间；
3. 分析前请访问各数据库官方，查询是否存在数据库版本的更新；

一、GreenGene数据库(gg)

FTP地址：ftp://greengenes.microbio.me/greengenes_release

最新版 QIIME 格式数据库下载链接：13_8代表13年8月更新

ftp://greengenes.microbio.me/greengenes_release/gg_13_5/gg_13_8_otus.tar.gz # 305 MB, 2013/8/15更新，序列与13_5一致，只是修正部分物种注释

解压后，通用只有97目录下的三个文件：

97_otu_taxonomy.txt.gz  # 97%聚类代表序列的物种注释，1.2 MB
97_otus.fasta.gz  # 97%聚类代表序列，29 MB
97_otus.tree.gz # 97%聚类代表序列的多序列树文件，929 KB


二、USEARCH格式数据库(usearch)

下载页面：http://www.drive5.com/sintax

16S

	rdp_16s_v18.fa.gz # 20,000万条序列，用于物种注释

	rdp_16s_v16_sp.fa.gz # 1.3万条序列，用于物种注释，物种注释是越准越好

	silva_16s_v123.fa.gz # 160万条，用于去嵌合体，防止假阴性去除嵌合体是越全越好


三、UNITE数据库 - ITS 数据

下载页面：https://unite.ut.ee/

	utax_reference_dataset_all_29.11.2022.fasta.gz # UNITE 最新版真菌序列及注释

USEARCH数据是二次整理，不一定及时更新，可至UNITE官网下载最新版、完整版

有问题，详见Pipeline.sh中附录 ## 11. UNITE数据库分析报错问题



