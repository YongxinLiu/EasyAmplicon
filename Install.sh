
[TOC]

# EasyAmplicon 2 Software and Database Installation Tutorial(软件与数据库安装教程)

    # Author作者: Yong-xin Liu(刘永鑫),,Tong Chen(陈同), Hao Luo(罗豪), Defeng Bai(白德凤), et al.
    # Update更新时间: 2025-10-17
    # Version版本: 2.0

Software download软件下载:
* GitHub: [https://github.com/YongxinLiu/EasyAmplicon](https://github.com/YongxinLiu/EasyAmplicon)*

---
## 0. Install安装 Git、R、RStudio

Please install the dependency software according to your system (Win/Mac/Linux).

# R 4.x.x is recommended for running R scripts
https://www.r-project.org/. It is also recommended that Rtools be installed for source code packages.
RStudio 2025.xx.x is a integrated development environment for R https://posit.co/download/rstudio-desktop/

# Git（仅Windows用户需单独安装）
# Git (Windows users need to install it separately)
Git for Windows 2.xx.x (Windows only) http://gitforwindows.org/

# R packages quick install
The statistics and visualization may require > 500 R packages. Installation is time-consuming and may also rely on other compilation tools. 
You can download all needed R packages in https://pan.baidu.com/s/1Ikd_47HHODOqC3Rcx6eJ6Q?pwd=0315 db/win/4.x.zip or db/mac/R4.2_mac_libraryX86_64.zip, then unzip and take the `4.x` folder in C:\Users\[$UserName]\AppData\Local\R\win-library\
    # 注意：如显示缺少某个R包，可以通过以下方法单独安装
    # Note: If an R package is missing, you can install it separately using the following methods
    # 例如DADA2包托管在Bioconductor上，需要通过BiocManager来安装。
    # For example, the DADA2 package is hosted on Bioconductor and needs to be installed via BiocManager.
    # 请打开您的R或Rstudio，在控制台(Console)中输入并执行以下命令：
    # Please open your R or Rstudio, and enter and execute the following commands in the Console:
    # 首先，安装Bioconductor的核心管理工具 BiocManager
     # First, install BiocManager's core management tool, BiocManager
    if (!requireNamespace("BiocManager", quietly = TRUE))
        install.packages("BiocManager")
    # 然后，通过 BiocManager 安装 DADA2
    # Then, install DADA2 via BiocManager
    BiocManager::install("dada2")
    # 此外，一些R包可以通过常规方式安装，比如需要 argparse 包来解析命令行参数，就可以使用：
    # In addition, some R packages can be installed in the conventional way. For example, if you need the argparse package to parse command line arguments, you can use:
    install.packages("argparse")
    # 安装各类R包时，控制台可能会出现大量的编译和安装信息，请耐心等待其完成。如果遇到问题，请根据提示信息检查您的R语言环境或网络连接。
    # When installing various R packages, a large amount of compilation and installation information may appear in the console. Please wait patiently for it to complete. If you encounter problems, please check your R language environment or network connection according to the prompt.



## 1. EasyAmplicon 2 conda enviroment

需要在Linux环境下和Windows中Ubuntu子系统中下运行

```bash
# 软件安装目录（推荐conda）
# Software installation directory (conda is recommended)
soft=~/miniconda3
# 数据库保存位置，Linux服务器建议~/db，子系统推荐 /mnt/d/EasAmplicon2/db，无d盘改为c
# Database storage location
wd=/mnt/d/amplicon2
db=/mnt/d/amplicon2/db

mkdir -p ${db}
export PATH=${soft}/bin:${soft}/condabin:${PATH}
echo $PATH
cd $wd
```

安装 Conda 和常用工具包 Conda Installation

```bash
# 下载最新版miniconda3 v24.9.2 , 安装日期2024/11/12, 141.47 Mb
# Download the latest version of miniconda3 v24.9.2, installation date 2024/11/12, 141.47 Mb
    wget -c https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
    
    # 安装，-b批量，-f无提示，-p目录，许可协议打yes
    # Install, -b for batch, -f for no prompt, -p for directory, type yes for license agreement
    bash Miniconda3-latest-Linux-x86_64.sh -b -f 
    # 激活，然后关闭终端重开，提示符前出现(base)即成功
    # Activate, then close and reopen the terminal, success if (base) appears before the prompt
    ~/miniconda3/condabin/conda init
    source ~/.bashrc
    # 查看版本，conda 25.5.1, python 3.13.5
    # Check version, conda 25.5.1, python 3.13.5
    conda -V  # 25.5.1
    python --version  # 3.13.5
    # Add frequently used channels添加常用频道
    conda config --add channels bioconda # 生物软件 / Bioconda for bioinformatics software
    conda config --add channels conda-forge # Highest priority / Conda-forge has the highest priority
    # Show channels 查看软件安装渠道
    conda config --show channels
  
    # conda默认配置文件为 ~/.condarc 查看配置文件位置
    # The default conda configuration file is ~/.condarc. Check the configuration file location
    # 你使用的是新版本 Conda（≥24），它要求用户必须手动接受各个源的服务条款（ToS），否则不能使用对应的频道。
    # You are using a new version of Conda (≥24), which requires users to manually accept the Terms of Service (ToS) for each source, otherwise the corresponding channel cannot be used.
    # conda tos accept --override-channels --channel https://repo.anaconda.com/pkgs/main
    # conda tos accept --override-channels --channel https://repo.anaconda.com/pkgs/r
    # mamba是用于管理环境的 CLI 工具。相比于 conda，mamba 是用 c++重写了 conda 的部分功能，运行效率显著提高，可以进行并行的下载，使用 rpm 包管理工具中的 libsolv，可以更快的解决环境依赖问题。
    # mamba is a CLI tool for managing environments. Compared to conda, mamba rewrites some of conda's functions in C++, significantly improving operational efficiency. It can perform parallel downloads and uses libsolv from the rpm package management tool to resolve environment dependencies faster.
    conda tos accept --override-channels --channel https://repo.anaconda.com/pkgs/main
    conda tos accept --override-channels --channel https://repo.anaconda.com/pkgs/r
    conda install mamba -y
    mamba install pandas -y
    mamba install conda-pack -y
    
    #conda config --set channel_priority strict #设置严格的仓库优先级（最好不要使用）/ Set strict channel priority (better not to use)
    #conda config --set channel_priority flexible #禁用仓库优先级 / Disable channel priority
    
    conda config --show-sources
    # 查看虚拟环境列表 
    # List conda environments
    conda env list
```
更多conda中文安装使用教程参考：[Nature Method：Bioconda解决生物软件安装的烦恼](https://mp.weixin.qq.com/s/SzJswztVB9rHVh3Ak7jpfA)
For more Chinese tutorials on conda installation and use, please refer to: [Nature Method: Bioconda solves the trouble of biological software installation](https://mp.weixin.qq.com/s/SzJswztVB9rHVh3Ak7jpfA)
[一文掌握Conda软件安装：虚拟环境、软件通道、加速solving、跨服务器迁移](https://mp.weixin.qq.com/s/tKAU09_w7Cu7khA9M2EGEQ)
[Mastering Conda Software Installation in One Article: Virtual Environments, Software Channels, Accelerating Solving, and Cross-Server Migration](https://mp.weixin.qq.com/s/tKAU09_w7Cu7khA9M2EGEQ)

---

## 2. 软件安装 Software Installation

### 2.1 Install安装 EasyAmplicon 2

下载EasyAmplicon最新版 https://github.com/YongxinLiu/EasyAmplicon 作为正对照模板

```bash
git clone https://github.com/YongxinLiu/EasyAmplicon
```

### 3. Install dependency 安装依赖 EasyMicrobiome

EasyAmplicon 2 依赖EasyMicrobiome，其包括众多脚本、软件和数据库的集合，网址：https://github.com/YongxinLiu/EasyMicrobiome
EasyAmplicon 2 depends on EasyMicrobiome, which includes a collection of many scripts, software and databases. Website: https://github.com/YongxinLiu/EasyMicrobiome

    # 方法1. 网页中下载
    # Method 1. Download from the webpage
    # https://github.com/YongxinLiu/EasyMicrobiome 中Code Download ZIP下载压缩包，上传至服务器，并解压
    # In https://github.com/YongxinLiu/EasyMicrobiome, click Code -> Download ZIP to download the compressed package, upload it to the server, and unzip it.
    unzip EasyMicrobiome-master.zip
    mv EasyMicrobiome-master EasyMicrobiome
    
    # 方法2. 备用链接下载
    # Method 2. Download from an alternative link
    wget -c ftp://download.nmdc.cn/tools/soft/EasyMicrobiome.tar.gz
    tar -xvzf EasyMicrobiome.tar.gz
    
    # 方法3. git下载，需安装git，注意网络问题
    # Method 3. Download with git, git needs to be installed, pay attention to network issues
    git clone https://github.com/YongxinLiu/EasyMicrobiome
    # 旧版更新
    # Update old version
    cd EasyMicrobiome && git pull && cd ../
    
    # 软件安装
    # Software installation
    # 添加linux命令可执行权限
    # Add executable permission to linux commands
    # cd /mnt/d/
    chmod +x `pwd`/EasyMicrobiome/linux/* `pwd`/EasyMicrobiome/script/*
    # 添加环境变量
    # Add environment variables
    echo "export PATH=\"$PATH:`pwd`/EasyMicrobiome/linux:`pwd`/EasyMicrobiome/script\"" >> ~/.bashrc
    source ~/.bashrc
    echo $PATH
```
### 4. Create environment and install core tools创建环境并安装核心工具

```bash
# Create, install and activate the EasyAmplicon2 environment
# **注：直接安装、下载解压安装，二选一。一种方法不成功，尝试另一种。**
# **Note: Choose one of the two options: direct installation or download and unzip for installation. If one method fails, try the other.**
cd $wd

# Method方法1.Download & install下载安装(recommended推荐)
### Specify conda file name指定conda文件名
s=easyamplicon2
soft=~/miniconda3
### Download and install 百度网盘下载链接：Baidu Net Disk：https://pan.baidu.com/s/1Ikd_47HHODOqC3Rcx6eJ6Q?pwd=0315
# File path文件路径：db/amplicon/amplicon2.tar.gz
### Specify installation directory 指定安装目录，~21s(SSD)
mkdir -p ${soft}/envs/${s}
time tar -xvzf ${s}.tar.gz -C ${soft}/envs/${s}
### 启动环境
### Activate the environment
conda activate ${s}
### 初始化环境 Initialize the environment
conda unpack

# Method方法2.Direct installation直接安装
conda env create -f amplicon2.yaml
conda activate amplicon2

## 方法3.下载singularity
## Method 3. Download and run with Singularity (recommended)
### Download Singularity image下载Singularity镜像
# 百度网盘下载链接 Baidu Net Disk: https://pan.baidu.com/s/1Ikd_47HHODOqC3Rcx6eJ6Q?pwd=0315  
# 文件路径 File path: db/amplicon/amplicon2.sif  
### Specify storage directory指定存放目录
mkdir -p ~/singularity/amplicon2
mv amplicon2.sif ~/singularity/amplicon2/
### Run the environment运行环境
singularity exec ~/singularity/amplicon2/amplicon2.sif bash
```

## 4. 数据库下载 Database Installation

**注：直接下载转换好格式的数据库、下载原始数据库，二选一。一种方法不成功，尝试另一种。**
**Note: Choose one of two options: direct installation or download and unzip for installation. If one method fails, try the other.**

### 方法1. 下载数据库 Download the pre-formatted database

数据库百度网盘下载链接：Baidu Net Disk：https://pan.baidu.com/s/1Ikd_47HHODOqC3Rcx6eJ6Q?pwd=0315 

文件路径 File path:  
- db/amplicon/silva/SILVA_modified.fasta  
- db/amplicon/silva/silva_nr99_v138.1_train_DADA2.fa.gz  
- db/amplicon/usearch/gtdb_sintax_database.fasta.gz  
- db/amplicon/usearch/sintax_defalut_emu_database.fasta.gz  
- db/amplicon/usearch/sintax_ncbi_database.fasta.gz  
- db/amplicon/GTDB
- db/amplicon/silva/Silva_Emu

全部存放在EasyMicrobiome/usearch文件夹中

### 方法2. 下载原始参考数据库及格式化Download the original reference databases (manual conversion required)

Silva: https://www.arb-silva.de/current-release/Exports
  - SILVA_138.2_SSURef_NR99_tax_silva.fasta.gz  
 
GTDB: https://data.gtdb.aau.ecogenomic.org/releases/release226/226.0/  
  - bac120_taxonomy_r226.tsv; bac120_ssu_reps_r226.fna.gz ; ssu_all_r226.fna.gz

NCBI: https://ftp.ncbi.nlm.nih.gov/blast/db/  
  - 16S_ribosomal_RNA.tar.gz  

Emu: https://osf.io/56uf7/files/osfstorage  
  - emu_default.tar.gz  

## 5. 添加环境变量Add environment variables (optional)

```bash
echo "export PATH=\"$PATH:${db}/amplicon2\"" >> ~/.bashrc
echo "export PATH=\"$PATH:${soft}/bin\"" >> ~/.bashrc
source ~/.bashrc
```

---


## 6. 软件打包Software packaging (optional)

```bash
mkdir -p ~/project/amplicon2/package
cd ~/project/amplicon2/package
conda activate amplicon2
conda pack -f -n amplicon2 -o amplicon2.tar.gz
```

---

## 7. 参考链接 Reference

* [amplicon2 GitHub](https://github\.com/YongxinLiu/EasyAmplicon)
* [QIIME2 官方教程 / Official Tutorial](https://docs.qiime2.org/)
* [PICRUSt2 文档 / Documentation](https://github.com/picrust/picrust2/wiki)
* [SILVA 数据库 / Database](https://www.arb-silva.de/)
* [UNITE 真菌数据库 / Fungal Database](https://unite.ut.ee/)
* [GTDB 官网 / Official Website](https://gtdb.ecogenomic.org/)
* [RDP Classifier](https://sourceforge.net/projects/rdp-classifier/)

---