# EasyAmplicon (易扩增子)

English Version in README.md

**An easy using, open-resource, reproducible, and community-based pipeline for amplicon data analysis in microbiome** (简单易学易用、开源、可重复且社区支持的扩增子数据分析流程)

Version (版本)：EasyAmplicon v1.13

Update (更新时间)：2021/9/6

## Pipeline Manual and file description (流程使用和文件介绍)

Using RStudio open the pipeline.sh

Files description:

- Readme.md # Introduction and install
- pipeline.sh # Command-line analysis for Windows and Linux
- pipeline_mac.sh # Command-line analysis for MacOS
- result/ # Example data
- result/Diversity-tutorial.Rmd # Interactive analysis in R and output reproducible report in HTML format

## What can we do? (结果展示)

- Analysis and visualization of microbiome data, especially for 16S rDNA amplicon;
- From raw data into feature tables;
- Support 20+ analysis methods and publish-read visualization;
- Finish your project at your laptop in 3 hours;
- Chinese manual and video supported.

![image](http://210.75.224.110/Note/R/amplicon/fig1.png)

**Figure 1. Pipeline of EasyAmplicon for analysis pair-end amplicon data.**


![image](http://210.75.224.110/Note/R/amplicon/fig2.png)

**Figure 2. Examples of visualizations.**

## 安装

系统要求 System requirement: Windows 10 / Mac OS 10.12+ / Ubuntu 16.04+

### 软件安装

- R语言环境R 4.x.x，支持R代码运行，下载适合自己系统的最新版并安装：https://www.r-project.org/ 
- R语言开发环境RStudio 1.x.xxxx，用于执行流程，下载适合自己系统的最新版并安装 ：https://www.rstudio.com/products/rstudio/download/#download
- (可选，仅Windows用户)Git Bash命令行环境的Git for Windows 2.xx.x，支持在Windows系统中运行Shell语言，下载并安装最新版：http://gitforwindows.org/
- 易扩增子流程EasyAmplicon v1.13，包括分析流程代码、测序数据和示例结果(分析的正对照)，可用 git clone https://github.com/YongxinLiu/EasyAmplicon 下载项目，或网页中下载打包的最新版并解压。
- 易微生物组EasyMicribome v1.13，提供易扩增子流程依赖的常用软件、脚本和数据库，可用 git clone https://github.com/YongxinLiu/EasyMicrobiome 下载项目，或网页中下载打包的最新版并解压，如Windows用户建议下载到C盘根目录。已经整合的重点软件如下：
    - usearch v10.0.240 https://www.drive5.com/usearch/download.html
    - vsearch v2.15.2 https://github.com/torognes/vsearch/releases

### (可选)扩展软件和数据库

- R包批量安装：R包通常会在使用过程中自动安装最新版，可选批量下载编绎好的R包合集。[Windows R 4.1](http://210.75.224.110/db/R/4.1.zip)，[Mac R 4.1](http://210.75.224.110/db/R/4.1_mac.zip)
- 16S数据库：16S常用RDP/SILVA/GreenGene/EzBioCloud数据库进行物种注释，默认下载了RDP和EzBioCloud。可以从上述数据库官网下载并整理为USEARCH使用的格式，此处推荐从[USEARCH官网](http://www.drive5.com/sintax)下载USEARCH兼容格式的数据库。默认流程使用体积小巧的RDP v18训练集数据库 (rdp_16s_v18.fa.gz)，并已保存于EasyMicrobiome/usearch目录中。可选GreenGenes 13.5 (gg_16s_13.5.fa.gz)和SILVA (silva_16s_v123.fa.gz) 数据库，从[USEARCH官网](http://www.drive5.com/sintax)根据需要下载并保存于usearch目录中。此外，如果要开展PICRUSt和Bugbase功能预测分析，还需要使用GreenGenes数据库13.5中按97%聚类的OTU序列 (己保存于流程gg目录中97_otus.fasta.gz)。该数据源于[GreenGenes官方](ftp://greengenes.microbio.me/greengenes_release)，解压后选择其中的97_otus.fasta保存于gg目录下
- ITS数据库：研究真菌或真核生物采用转录间隔区 (Intergenic Transcribed Spacer) 测序，需要使用UNITE数据库，目前最新版已经保存于EasyMicrobiome\usearch目录(unite_all_2021.5.10.gz)。如流程中数据库没有及时更新，可在UNITE官网 (https://unite.ut.ee/) 下载适合USEARCH的最新版注释数据库。官方数据库存在格式问题，详细常见pipeline.sh中附录常用问题

## 运行

1. 准确输入数据：典型的扩增子测序起始文件包括测序数据和元数据两类。

测序数据(*.fq.gz)为seq/目录中的成对fastq/fq文件，通常采用.gz的压缩格式保存节省空间。元数据(metadata.txt)为按样本编号对应的分组、时间、地点等描述信息。EasyAmplicon项目中有准备好的demo数据用于测试分析流程是否可以正常工作(正对照)，同时提供标准格式的参考模板，指导用户准备标准的输入数据。

然后用RStudio打开pipeline.sh即可开始分析。新项目在准备好测序数据和元数据后，复制EasyAmplicon中的pipeline.sh至新项目文件夹并用RStudio打开即可开始分析之旅。

2. 开始分析流程

参考Pipeline.sh中的代码，按说明设置工作目录、脚本和数据库(EasyMicrobiome)位置等，在RStudio中逐行或逐段选择代码并运行(Run)即可完成整套分析流程。

主要数据分析步骤如下：
- 合并双端序列并按样品重命名
- 引物切除和质量控制
- 序列去冗余并挑选代表序列
- 特征表生成和筛选
- Alpha多样性计算
- Beta多样性计算
- 物种注释分类汇总
- 有参分析和功能预测
- 空间清理及数据提交
- STAMP和LEfSe软件输入文件准备

每步骤参数和结果的详细解读，详见 《易扩增子：易用、可重复和跨平台的扩增子分析流程》https://bio-protocol.org/bio101/e2003641




## FAQ (常见问题)

