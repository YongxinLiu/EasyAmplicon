# 易扩增子 (EasyAmplicon)

英文版见(English Version in) README.md

版本：1.20

更新日期：2023/10/13

简单易学易用、开源、可重复且社区支持的扩增子数据分析流程

第一次使用参考`安装`段落下载并安装流程及依赖关系。使用RStudio打开`EasyAmplicon/pipeline.sh`即可逐行完成扩增子分析。

## 简介 (Introduction)

文件描述:

-   Readme.md \# 英文版帮助
-   Readme_cn.md \# 中文版帮助
-   pipeline.sh \# Windows或Linux版命令行分析流程
-   pipeline_mac.sh \# MacOS版命令行分析流程
-   result/ \# 示例结果(正对照)
-   result/Diversity.Rmd \# 交互式多样性可重复分析代码，可编译为HTML网页或Word文档报告

主要功能：

-   分析和可视化微生物组数据，尤其是16S rDNA扩增子测序
-   从原始数据到特征表的端对端
-   支持20余种分析方法，并生成出版级图表
-   在普通个人电脑上3小完成示例项目
-   中、英文双语帮助文档，中文视频教程支持确保可重复

![image](https://github.com/YongxinLiu/EasyAmplicon/blob/master/result/Figure1.jpg)

**图1. 易扩增子分析双端扩增子数据的流程**

![image](https://github.com/YongxinLiu/EasyAmplicon/blob/master/result/Figure2.jpg)

**图2. 结果的部分可视化示例**

## 安装(Install)

系统要求 System requirement: Windows 10+ / Mac OS 10.12+ / Ubuntu 20.04+

安装视频教程：https://www.bilibili.com/video/BV1Cb411f7La/

### 依赖软件环境(Install Dependency)

请安装与你操作系统一致的软件

-   R语言环境R 4.x.x，支持R代码运行：<https://www.r-project.org/>，推荐下载Rtools实现包的源码安装
-   R语言开发环境RStudio 2023.xx.x，用于执行流程：<https://posit.co/download/rstudio-desktop/>
-   STAMP v2.1.3 特征表统计和分析图型界面软件 <http://kiwi.cs.dal.ca/Software/STAMP>
-   (可选，仅Windows用户)Git Bash命令行环境的Git for Windows
    2.xx.x，支持在Windows系统中运行Shell语言，下载并安装最新版：<http://gitforwindows.org/>

以最常用的Windows系统(87.5%)为例，你可以快速下载上面的软件安装包：[Git for Windows](https://gitforwindows.org/)、[R](https://mirrors.tuna.tsinghua.edu.cn/CRAN/bin/windows/base/)、[RStudio](https://posit.co/download/rstudio-desktop/)、[STAMP](https://github.com/dparks1134/STAMP/releases/download/v2.1.3/STAMP_2_1_3.exe)，[合集见百度网盘db/win目录](https://pan.baidu.com/s/1Ikd_47HHODOqC3Rcx6eJ6Q?pwd=0315)。

- (可选，推荐)R包的快速安装

在R语言的统计和可视化中会使用超过500个R包，安装不仅费时费力，而且经常出错或依赖其他编绎工具。为方便大家使用，我们提供了编绎好的R包合集下载，如 [Windows版](https://pan.baidu.com/s/1Ikd_47HHODOqC3Rcx6eJ6Q?pwd=0315)、[Mac版](https://pan.baidu.com/s/1Ikd_47HHODOqC3Rcx6eJ6Q?pwd=0315 )。你可以下载解压后，将 `4.3` 目录移动至 `C:\Users\[$UserName]\AppData\Local\R\win-library\`中即完成安装。

### 安装易扩增子 (Install EasyAmplicon)

-   易扩增子流程EasyAmplicon，包括分析流程代码、测序数据和示例结果(分析的正对照)， <https://github.com/YongxinLiu/EasyAmplicon>
-   易微生物组EasyMicribome，提供易扩增子流程依赖的常用软件、脚本和数据库，  <https://github.com/YongxinLiu/EasyMicrobiome>

下载以上项目至C或D盘，并解压。以下提供三种下载方式可选(让你永远留有后手)

- 方法1. 网页下载。访问项目主页，点击 `Code` -- `Download`，选择下载位置，开始下载
- 方法2. git命令行下载。直接生成目录，无需解压。`git clone https://github.com/YongxinLiu/EasyAmplicon`和`git clone https://github.com/YongxinLiu/EasyMicrobiome`。 注：提示`fatal: unable to access`可能只是网络问题，重试或改天重试一般可解决，或找代理或朋友帮忙下载。
- 方法3. 直接从国内百度网盘链接中db/soft目录下载: [EasyAmplicon](https://pan.baidu.com/s/1Ikd_47HHODOqC3Rcx6eJ6Q?pwd=0315 )、[EasyMicrobiome](https://pan.baidu.com/s/1Ikd_47HHODOqC3Rcx6eJ6Q?pwd=0315 )

### (可选)扩展软件和数据库

-   Rtools：用于从源码进行R包的安装，windows版本见：https://cran.rstudio.com/bin/windows/Rtools/
-   16S数据库：16S常用RDP/SILVA/GreenGene/EzBioCloud数据库进行物种注释，默认下载了RDP和EzBioCloud。可以从上述数据库官网下载并整理为USEARCH使用的格式，此处推荐从[USEARCH官网](http://www.drive5.com/sintax)下载USEARCH兼容格式的数据库。默认流程使用体积小巧的RDP
    v18训练集数据库
    (rdp_16s_v18.fa.gz)，并已保存于EasyMicrobiome/usearch目录中。可选GreenGenes
    13.5 (gg_16s_13.5.fa.gz)和SILVA (silva_16s_v123.fa.gz)
    数据库，从[USEARCH官网](http://www.drive5.com/sintax)根据需要下载并保存于usearch目录中。此外，如果要开展PICRUSt和Bugbase功能预测分析，还需要使用GreenGenes数据库13.5中按97%聚类的OTU序列
    (己保存于流程gg目录中97_otus.fasta.gz)。该数据源于[GreenGenes官方](ftp://greengenes.microbio.me/greengenes_release)，解压后选择其中的97_otus.fasta保存于gg目录下
-   ITS数据库：研究真菌或真核生物采用转录间隔区 (Intergenic Transcribed
    Spacer)
    测序，需要使用UNITE数据库，目前最新版已经保存于EasyMicrobiome\usearch目录。如流程中数据库没有及时更新，可在UNITE官网
    (<https://unite.ut.ee/>)
    下载适合USEARCH的最新版注释数据库。官方数据库存在格式问题，详细常见pipeline.sh中附录常用问题

## 快速运行(Quick Start)

使用视频教程：https://www.bilibili.com/video/BV1is4y157Ms/

1.  准确输入数据：典型的扩增子测序起始文件包括测序数据和元数据两类。

测序数据(\*.fq.gz)为seq/目录中的成对fastq/fq文件，通常采用.gz的压缩格式保存节省空间。元数据(metadata.txt)为按样本编号对应的分组、时间、地点等描述信息。EasyAmplicon项目中有准备好的demo数据用于测试分析流程是否可以正常工作(正对照)，同时提供标准格式的参考模板，指导用户准备标准的输入数据。

新项目在准备好测序数据和元数据后，复制EasyAmplicon中的pipeline.sh至新项目文件夹。然后用RStudio打开pipeline.sh即可开始分析之旅。

2.  开始分析流程

参考Pipeline.sh中的代码，按说明设置工作目录(work directory)、脚本和数据库(EasyMicrobiome)位置等，在RStudio中逐行或逐段选择代码并运行(Run)即可完成整套分析流程。

主要数据分析步骤如下： - 合并双端序列并按样品重命名 - 引物切除和质量控制
- 序列去冗余并挑选代表序列 - 特征表生成和筛选 - Alpha多样性计算 -
Beta多样性计算 - 物种注释分类汇总 - 有参分析和功能预测 -
空间清理及数据提交 - STAMP和LEfSe软件输入文件准备

每步骤参数和结果的详细解读，详见
《易扩增子：易用、可重复和跨平台的扩增子分析流程》<https://bio-protocol.org/bio101/e2003641>

## 常见问题 (FAQ)

pipeline.sh 中的常见问题
注：.sh脚本全部为markdown格式，使用有道Note或VSCode，阅读体验更佳。

## 更新日志 (Change log)

2023/3/11 1.18.1
解决备用下载链接失效问题？视频转移到B站链接，下载文件提供百度云链接。

2023/6/4 1.19
R和Rtools更新为4.3，RStudio更新为2023.03.1

## 引文 (Citation)

使用此脚本，请引用下文：

**Yong-Xin Liu**, Lei Chen, Tengfei Ma, Xiaofang Li, Maosheng Zheng, Xin Zhou, Liang Chen, Xubo Qian, Jiao Xi, Hongye Lu, Huiluo Cao, Xiaoya Ma, Bian Bian, Pengfan Zhang, Jiqiu Wu, Ren-You Gan, Baolei Jia, Linyang Sun, Zhicheng Ju, Yunyun Gao, **Tao Wen**, **Tong Chen**. 2023. EasyAmplicon: An easy-to-use, open-source, reproducible, and community-based pipeline for amplicon data analysis in microbiome research. **iMeta** 2: e83. https://doi.org/10.1002/imt2.83

版本所有 2016-2023 刘永鑫(Yong-Xin Liu) <liuyongxin@caas.cn>, 文涛(Tao Wen) <taowen@njau.edu.cn>, 陈同(Tong Chen) <chent@nrc.ac.cn>
