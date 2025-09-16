# EasyAmplicon2 Snakemake Pipeline

这是一个基于Snakemake的三代扩增子分析流程。

## 特征 Features

- **全自动化流程**: 从原始FASTQ文件到最终分析结果
- **模块化设计**: 每个分析步骤都是独立的Snakemake规则
- **可配置参数**: 通过config.yaml文件轻松调整参数
- **环境管理**: 使用conda环境确保软件版本一致性
- **质量控制**: 内置多个质量控制步骤
- **多样性分析**: α和β多样性分析
- **物种注释**: 基于SILVA/自定义数据库的物种分类
- **可选模块**: PICRUSt2功能预测和系统发育树构建

## 文件结构 Directory Structure

```
EasyAmplicon2/
├── Snakefile                   # 主工作流文件
├── config/
│   └── config.yaml             # 配置文件
├── environment.yml             # conda环境文件
├── scripts/                    # R脚本
│   ├── otutab_filter_nonBac.R
│   └── otutab_rare.R
├── schemas/                    # 数据验证schema
│   └── samples.schema.yaml
├── seq/                        # 原始测序数据目录
│   ├── HY1_1.fq.gz
│   └── ...
├── metadata.txt               # 样本元数据文件
└── results/                   # 分析结果目录
```

## 安装 Installation

```bash
    **注：直接安装、下载解压安装，二选一。一种方法不成功，尝试另一种。**

    cd EasyAmplicon2

    ## 方法1.直接安装
    conda env create -f EasyAmplicon2.yaml
    conda activate easyamplicon2

    ## 方法2.下载安装(推荐)
    ### 指定conda文件名
    s=easyamplicon2
    soft=~/miniconda3

    ### 下载安装
    百度网盘下载链接：Baidu Net Disk：https://pan.baidu.com/s/1Ikd_47HHODOqC3Rcx6eJ6Q?pwd=0315
    文件路径：db/amplicon/easyamplicon2.tar.gz

    ### 指定安装目录
    mkdir -p ${soft}/envs/${s}
    tar -xvzf ${s}.tar.gz -C ${soft}/envs/${s}

    ### 启动环境
    conda activate ${s}

    ### 初始化环境
    ### easyamplicon2环境包含了大部分分析软件
    conda unpack
```


## 使用方法 Usage

### 1. 准备输入文件

- 将测序数据放在`seq/`目录下，文件命名格式为`{SampleID}.fq.gz`
- 编辑`metadata.txt`文件，包含样本信息
- 根据需要修改`config/config.yaml`配置文件

### 2. 运行流程

```bash
# 激活easyamplicon2环境（包含snakemake）
conda activate easyamplicon2


# 方法1: 直接运行snakemake
snakemake -s Snakefile --cores 8

# 方法2: 先检查，再运行
snakemake -s Snakefile -n  # 干跑检查
snakemake -s Snakefile --cores 8 --printshellcmds  # 正式运行

```

### 3. 主要输出文件

- `results/feature_table/otutab_raw.txt`: 特征表(OTU/ASV表)
- `results/taxonomy/taxonomy.txt`: 物种注释表
- `results/alpha_diversity/alpha.txt`: α多样性指数
- `results/beta_diversity/`: β多样性距离矩阵
- `results/taxonomic_summary/`: 各分类水平汇总表
- `results/summary_report.txt`: 分析总结报告

## 配置说明 Configuration

主要配置项在`config/config.yaml`中：

- `threads`: 并行线程数
- `primers`: 引物序列
- `databases`: 数据库路径
- `rarefaction.depth`: 稀释深度
- `run_picrust2`: 是否运行功能预测
- `build_tree`: 是否构建系统发育树

## 故障排除 Troubleshooting

### 常见问题

1. **内存不足**: 调整配置文件中的`resources`设置
2. **样本名包含点号**: 确保样本ID不包含点号"."
3. **数据库路径错误**: 检查config.yaml中的数据库路径
4. **USEARCH许可**: 确保USEARCH安装正确且有使用权限

### 获取帮助

- 检查Snakemake日志文件
- 使用`snakemake --detailed-summary`查看详细状态
- 查看各个conda环境是否正确安装


