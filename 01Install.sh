[TOC]

# EasyAmplicon2 软件与数据库安装教程

*Version: 2.0*
*OS: Ubuntu 22.04+ / CentOS 8+*

*Homepage: [https://github\.com/YongxinLiu/EasyAmplicon](https://github.com/YongxinLiu/EasyAmplicon)*

软件下载：

* GitHub：[https://github\.com/YongxinLiu/EasyAmplicon](https://github.com/YongxinLiu/EasyAmplicon)

---

## 一、初始化环境 Initialization

```bash
# 软件安装目录（推荐conda）
soft=~/miniconda3
# 数据库保存位置
db=~/db

mkdir -p ${soft} ${db}
export PATH=${soft}/bin:${soft}/condabin:${PATH}
echo $PATH
```

---

## 二、安装 Conda 和常用工具包 Conda Installation

```bash
# 下载最新版miniconda3 v24.9.2 , 安装日期2024/11/12, 141.47 Mb  
    wget -c https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
    
    # 安装，-b批量，-f无提示，-p目录，许可协议打yes
    bash Miniconda3-latest-Linux-x86_64.sh -b -f 
    # 激活，然后关闭终端重开，提示符前出现(base)即成功
    ~/miniconda3/condabin/conda init
    source ~/.bashrc
    # 查看版本，conda 24.9.2, python 3.12.2
    conda -V  # 24.9.2
    python --version  # 3.12.2
    # 添加常用频道
    conda config --add channels bioconda # 生物软件
    conda config --add channels conda-forge # Highest priority
    
    # conda默认配置文件为 ~/.condarc 查看配置文件位置
    conda install mamba -c conda-forge -c bioconda -y
    mamba install pandas -c conda-forge -c bioconda -y
    mamba install conda-pack -c conda-forge -c bioconda -y
    
    #conda config --set channel_priority strict #设置严格的仓库优先级（最好不要使用）
    #conda config --set channel_priority flexible #禁用仓库优先级
    
    conda config --show-sources
    # 查看虚拟环境列表 
    conda env list

更多conda中文安装使用教程参考：[Nature Method：Bioconda解决生物软件安装的烦恼](https://mp.weixin.qq.com/s/SzJswztVB9rHVh3Ak7jpfA)
[一文掌握Conda软件安装：虚拟环境、软件通道、加速solving、跨服务器迁移](https://mp.weixin.qq.com/s/tKAU09_w7Cu7khA9M2EGEQ)
```

---

## 三、软件工具安装 Software Installation

### 1. 安装 Git、R、RStudio
```bash
Please install the dependency software according to your system (Win/Mac/Linux).

# R 4.x.x is recommended for running R scripts
https://www.r-project.org/. It is also recommended that Rtools be installed for source code packages.
RStudio 2025.xx.x is a integrated development environment for R https://posit.co/download/rstudio-desktop/

# Git（仅Windows用户需单独安装）
Git for Windows 2.xx.x (Windows only) http://gitforwindows.org/

# R packages quick install
The statistics and visualization may require > 500 R packages. Installation is time-consuming and may also rely on other compilation tools. 
You can download all needed R packages in https://pan.baidu.com/s/1Ikd_47HHODOqC3Rcx6eJ6Q?pwd=0315 db/win/4.x.zip or db/mac/R4.2_mac_libraryX86_64.zip, then unzip and take the `4.x` folder in C:\Users\[$UserName]\AppData\Local\R\win-library\

```


### 2. 安装 EasyAmplicon2 

```bash
cd ${db}
git clone https://github\.com/YongxinLiu/EasyAmplicon
cd EasyAmplicon2
chmod +x *.sh
```
### 3. 安装 EasyMicrobiome

```bash
EasyMetagenome依赖EasyMicrobiome，其包括众多脚本、软件和数据库的集合，网址：https://github.com/YongxinLiu/EasyMicrobiome
    
    # 方法1. 网页中下载
    # https://github.com/YongxinLiu/EasyMicrobiome 中Code Download ZIP下载压缩包，上传至服务器，并解压
    unzip EasyMicrobiome-master.zip
    mv EasyMicrobiome-master EasyMicrobiome
    
    # 方法2. 备用链接下载
    wget -c ftp://download.nmdc.cn/tools/soft/EasyMicrobiome.tar.gz
    tar -xvzf EasyMicrobiome.tar.gz
    
    # 方法3. git下载，需安装git
    git clone https://github.com/YongxinLiu/EasyMicrobiome
    # 旧版更新
    cd EasyMicrobiome && git pull && cd ../
    
    # 软件安装
    # 添加linux命令可执行权限
    chmod +x `pwd`/EasyMicrobiome/linux/* `pwd`/EasyMicrobiome/script/*
    # 添加环境变量
    echo "export PATH=\"\$PATH:`pwd`/EasyMicrobiome/linux:`pwd`/EasyMicrobiome/script\"" >> ~/.bashrc
    source ~/.bashrc
    echo $PATH
```
### 4. 创建环境并安装核心工具

```bash
# 创建并激活环境
conda create -n easyamplicon2 -y
conda activate easyamplicon2

# 安装核心工具
mamba install cutadapt -y
mamba install -c bioconda picrust2 -y
mamba install -c bioconda emu -y
```


## 四、数据库下载 Database Installation

### 1. SILVA 数据库

```bash
mkdir -p ${db}/silva
cd ${db}/silva
wget -c https://www.arb-silva.de/fileadmin/silva_databases/current/Exports/SILVA_138.1_SSURef_NR99_tax_silva.fasta.gz
gunzip SILVA_138.1_SSURef_NR99_tax_silva.fasta.gz
```

### 2. GTDB 数据库

```bash
mkdir -p ${db}/gtdb
cd ${db}/gtdb
wget https://data.gtdb.ecogenomic.org/releases/latest/auxillary_files/gtdb_taxonomy.tsv
```

### 3. UNITE 数据库（真菌）

```bash
mkdir -p ${db}/unite
cd ${db}/unite
wget https://unite.ut.ee/sh_files/UNITE_public_mothur_01.12.2023.fasta.gz
gunzip UNITE_public_mothur_01.12.2023.fasta.gz
```

### 4. GreenGene 数据库

```bash
mkdir -p ${db}/greengene
cd ${db}/greengene
wget https://greengenes.secondgenome.com/downloads/13_5/gg_13_5_otus.tar.gz
tar xvzf gg_13_5_otus.tar.gz
```

### 5. RDP 分类器数据库

```bash
mkdir -p ${db}/rdp
cd ${db}/rdp
wget https://sourceforge.net/projects/rdp-classifier/files/RDP_Classifier_2.13.zip
unzip RDP_Classifier_2.13.zip
```

---

## 五、添加环境变量（可选）

```bash
echo "export PATH=\"\$PATH:${db}/EasyAmplicon2\"" >> ~/.bashrc
echo "export PATH=\"\$PATH:${soft}/bin\"" >> ~/.bashrc
source ~/.bashrc
```

---


## 六、软件打包（可选）

```bash
mkdir -p ~/project/EasyAmplicon2/package
cd ~/project/EasyAmplicon2/package
conda activate easyamplicon2
conda pack -f -n easyamplicon2 -o easyamplicon2.tar.gz
```

---

## 七、参考链接 Reference

* [EasyAmplicon2 GitHub](https://github\.com/YongxinLiu/EasyAmplicon)
* [QIIME2 官方教程](https://docs.qiime2.org/)
* [PICRUSt2 文档](https://github.com/picrust/picrust2/wiki)
* [SILVA 数据库](https://www.arb-silva.de/)
* [UNITE 真菌数据库](https://unite.ut.ee/)
* [GTDB 官网](https://gtdb.ecogenomic.org/)
* [RDP Classifier](https://sourceforge.net/projects/rdp-classifier/)

---

