

# EasyAmplicon (易扩增子)

The Chinese version in (中文版见) [README_cn.md](./README_cn.md)

**An easy using, open-resource, reproducible, and community-based pipeline for amplicon data analysis in microbiome**

Version：v1.14

Update：2022/1/7

## Pipeline Manual and file description (流程使用和文件介绍)

Using RStudio open the pipeline.sh

Files description:

- Readme.md # Introduction and install
- pipeline.sh # Command-line analysis for Windows and Linux
- pipeline_mac.sh # Command-line analysis for MacOS
- result/ # Example data
- result/Diversity.Rmd # Interactive analysis in R and output reproducible report in HTML format

## What can we do? (结果展示)

- Analysis and visualization of microbiome data, especially for 16S rDNA amplicon;
- From raw data into feature tables;
- Support 20+ analysis methods and publish-ready visualization;
- Finish your project at your laptop in 3 hours;
- Chinese/English manual and video supported.

![image](http://210.75.224.110/Note/R/amplicon/fig1.png)

**Figure 1. Pipeline of EasyAmplicon for analysis pair-end amplicon data.**


![image](http://210.75.224.110/Note/R/amplicon/fig2.png)

**Figure 2. Examples of visualizations.**

## Install (安装)

### Install Dependency (安装依赖软件)

Please install the dependency software according with you system (Win/Mac/Linux).

- R 4.x.x for run R scripts https://www.r-project.org/
- RStudio 2021.xx.x is a integrated development environment for R https://www.rstudio.com/products/rstudio/download/#download
- STAMP v2.1.3 http://kiwi.cs.dal.ca/Software/STAMP 
- Git for Windows 2.xx.x (Windows only) http://gitforwindows.org/

As the windows users (87.5%) as example, you can quickly batch download above [R 4.1.2](https://mirrors.tuna.tsinghua.edu.cn/CRAN/bin/windows/base/R-4.1.2-win.exe), [RStudio 2021.09.1](https://download1.rstudio.org/desktop/windows/RStudio-2021.09.1-372.exe), [STAMP v2.1.3](https://github.com/dparks1134/STAMP/releases/download/v2.1.3/STAMP_2_1_3.exe), [Git 2.34.1](https://github.com/git-for-windows/git/releases/download/v2.34.1.windows.1/Git-2.34.1-64-bit.exe) in [zip](http://210.75.224.110/db/win.tar.gz). 

- R packages quick install

The statistics and visualization may require > 500 R packages. Installation is time-consuming and may also rely on other compilation tools. You can download all needed R packages in [zip for windows](http://210.75.224.110/R/4.1.zip) or [Mac](http://210.75.224.110/db/R/4.1_mac.zip), then unzip and take the `4.1` folder in C:\Users\[$UserName]\Documents\R\win-library\.

### Install EasyAmplicon (安装易扩增子)

- EasyAmplicon pipeline (Positive control) https://github.com/YongxinLiu/EasyAmplicon
- EasyMicrobiome include scripts and databases https://github.com/YongxinLiu/EasyMicrobiome

Download the the project in C: or D:, then unzip (three methods keep always have backup)

- Method 1. Visit the homepage, Code -- Download
- Method 2. Download by the mirror site in China: [EasyAmplicon](http://210.75.224.110/db/EasyAmplicon.tar.gz) [EasyMicrobiome](http://210.75.224.110/db/EasyMicrobiome.tar.gz)
- Method 3. `git clone https://github.com/YongxinLiu/EasyAmplicon` and `git clone https://github.com/YongxinLiu/EasyMicrobiome`. Note: `fatal: unable to access` can retry.

## Quick Start (快速运行)

Using Windows 10 as example:

1. Open RStudio, set termianl as Git Bash (Tools -- Global Options -- Terminal -- New termianls -- Git Bash -- OK)
2. File -- Open File -- `EasyAmplicon` folder -- pipeline.sh (windows/linux) or pipeline_mac.sh (mac)
3. Setup the `work directory`(wd), and `EasyMicrobiome directory`(db), then run each line by click run in top right corner

## Example dataset (示例数据)

- seq/ # raw sequencing in zipped fastq format, backup can download by metadata in GSA https://ngdc.cncb.ac.cn/gsa/
- result/ # Example data and figures for standard pipeline, such as alpha, beta, tax
- advanced/ # Example of advanced analysis, included data, scripts and output figures

## FAQ (常见问题)

Frequenty Asked Questions in pipeline.sh

Note: All the .sh script is writting in markdown format, using Youdao Note or VSCode for better reading experience.

## Citation (引文)

The paper is in submittion.

Please cite us review:

Yong-Xin Liu, Yuan Qin, Tong Chen, Meiping Lu, Xubo Qian, Xiaoxuan Guo & Yang Bai. (2021). A practical guide to amplicon and metagenomic analysis of microbiome data. Protein & Cell 12, 315-330, doi: https://doi.org/10.1007/s13238-020-00724-8


