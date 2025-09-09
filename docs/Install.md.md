# Install Dependency (安装依赖软件)
## All the software backup can be found in

FTP: [Filezilla](https://filezilla-project.org/index.php) visiting FTP download.nmdc.cn in anonymous. In tools directory, you can find all the software and packages in amplicon and different system supporting such mac, win
Baidu Net Disk：https://pan.baidu.com/s/1Ikd_47HHODOqC3Rcx6eJ6Q?pwd=0315

## Please install the dependency software according with you system (Win/Mac/Linux).

-- R 4.x.x for run R scripts https://www.r-project.org/, also recommended install Rtools for install source code packages.

-- RStudio 2025.xx.x is a integrated development environment for R https://posit.co/download/rstudio-desktop/

-- STAMP v2.1.3 http://kiwi.cs.dal.ca/Software/STAMP

-- Git for Windows 2.xx.x (Windows only) http://gitforwindows.org/

-- R packages quick install

-- The statistics and visualization may require > 500 R packages. Installation is time-consuming and may also rely on other compilation tools. You can download all needed R packages in https://pan.baidu.com/s/1Ikd_47HHODOqC3Rcx6eJ6Q?pwd=0315 db/win/4.x.zip or db/mac/R4.2_mac_libraryX86_64.zip, then unzip and take the 4.x folder in C:\Users[$UserName]\AppData\Local\R\win-library\

Note: If an R package is missing, you can install it separately using the following methods
For example, the DADA2 package is hosted on Bioconductor and needs to be installed via BiocManager.
Please open your R or Rstudio, and enter and execute the following commands in the Console:

```r
# First, install BiocManager's core management tool, BiocManager
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
# Then, install DADA2 via BiocManager
BiocManager::install("dada2")
# In addition, some R packages can be installed in the conventional way. For example, if you need the argparse package to parse command line arguments, you can use:
install.packages("argparse")
# When installing various R packages, a large amount of compilation and installation information may appear in the console. Please wait patiently for it to complete. If you encounter problems, please check your R language environment or network connection according to the prompt.
```


## Install EasyAmplicon 2(安装易扩增子2)

### Method 1. Visit the GitHub homepage, Code -- Download

EasyAmplicon pipeline (Positive control) https://github.com/YongxinLiu/EasyAmplicon

EasyMicrobiome include scripts and databases https://github.com/YongxinLiu/EasyMicrobiome

Download the the project in C: or D:, then unzip (keep the directoray name exact the software name)

### Method 2. Download by the mirror site in BaiduNetDisk: https://pan.baidu.com/s/1Ikd_47HHODOqC3Rcx6eJ6Q?pwd=0315 db/soft/EasyAmplicon.tar.gz or EasyMicrobiome.tar.gz

### Method 3. git clone https://github.com/YongxinLiu/EasyAmplicon and git clone https://github.com/YongxinLiu/EasyMicrobiome. Note: fatal: unable to access can retry.

## Install Conda (安装Conda)

```bash
# Download the latest version of miniconda3 v24.9.2, installation date 2024/11/12, 141.47 Mb
wget -c https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh

# Install, -b for batch, -f for no prompt, -p for directory, type yes for license agreement
bash Miniconda3-latest-Linux-x86_64.sh -b -f 

# Activate, then close and reopen the terminal, success if (base) appears before the prompt
~/miniconda3/condabin/conda init
source ~/.bashrc

# Check version, conda 25.5.1, python 3.13.5
conda -V  # 25.5.1
python --version  # 3.13.5

# Add frequently used channels
conda config --add channels bioconda # 生物软件 / Bioconda for bioinformatics software
conda config --add channels conda-forge # Highest priority / Conda-forge has the highest priority

# The default conda configuration file is ~/.condarc. Check the configuration file location

# You are using a new version of Conda (≥24), which requires users to manually accept the Terms of Service (ToS) for each source, otherwise the corresponding channel cannot be used.
conda tos accept --override-channels --channel https://repo.anaconda.com/pkgs/main
conda tos accept --override-channels --channel https://repo.anaconda.com/pkgs/r

# mamba is a CLI tool for managing environments. Compared to conda, mamba rewrites some of conda's functions in C++, significantly improving operational efficiency. It can perform parallel downloads and uses libsolv from the rpm package management tool to resolve environment dependencies faster.
conda install mamba -y
mamba install pandas -y
mamba install conda-pack -y

conda config --set channel_priority strict #设置严格的仓库优先级（最好不要使用）/ Set strict channel priority (better not to use)
conda config --set channel_priority flexible #禁用仓库优先级 / Disable channel priority

conda config --show-sources

# List conda environments
conda env list
# For more Chinese tutorials on conda installation and use, please refer to: [Nature Method: Bioconda solves the trouble of biological software installation](https://mp.weixin.qq.com/s/SzJswztVB9rHVh3Ak7jpfA)

[Mastering Conda Software Installation in One Article: Virtual Environments, Software Channels, Accelerating Solving, and Cross-Server Migration](https://mp.weixin.qq.com/s/tKAU09_w7Cu7khA9M2EGEQ)

```

## Create, install and activate the easyamplicon2 environment

**Note: Choose one of the two options: direct installation or download and unzip for installation. If one method fails, try the other.**
cd EasyAmplicon2

```bash

## Method 1. Direct installation
conda env create -f EasyAmplicon2.yaml
conda activate easyamplicon2

## Method 2. Download and install (recommended)
### Specify conda file name
s=easyamplicon2
soft=~/miniconda3

### Download and install

# 百度网盘下载链接：Baidu Net Disk：https://pan.baidu.com/s/1Ikd_47HHODOqC3Rcx6eJ6Q?pwd=0315
# 文件路径：db/amplicon/easyamplicon2.tar.gz
# File path: db/amplicon/easyamplicon2.tar.gz

### Specify installation directory
mkdir -p ${soft}/envs/${s}
tar -xvzf ${s}.tar.gz -C ${soft}/envs/${s}

### Activate the environment
conda activate ${s}

### Initialize the environment

### The easyamplicon2 environment contains most of the analysis software
conda unpack

```