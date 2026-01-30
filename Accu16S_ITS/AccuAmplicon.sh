[TOC]

# 扩增子绝对定量分析流程 For the absolute quantification of amplicon sequencing using Accu16S/ITS

    # 版本Version: v1.00, 2025/10/29
    # 操作系统：Operation System: Linux Ubuntu 22.04+/ CentOS 7.7+
    # 以下在windows 11下的Linux子系统运行
    
# 0.软件安装Software installation

    # blat, fastx_toolkit, seqtk, parallel
    conda create -n meta_env -y
    conda activate meta_env
    conda install -c bioconda blat fastx_toolkit seqtk parallel -y
    # 另外安装基础工具（通常已有）
    conda install -c conda-forge coreutils gzip awk gawk -y
    
    # 安装seqtk（如果未安装）
    conda install -c bioconda seqtk
    
    # 安装qiime2
    # 安装参考https://mp.weixin.qq.com/s/FZRBp3xFzN6aZvu5rBDoDg
    wget -c http://www.imeta.science/db/qiime2/qiime2-amplicon-ubuntu-latest-conda.yml
    conda env create \
    --name qiime2-amplicon-2025.4 \
    --file qiime2-amplicon-ubuntu-latest-conda.yml
    # 检查是否安装成功
    conda activate qiime2-amplicon-2025.4
    qiime --help
    
    
    ## 准备工作 Preparing
    # 设置工作路径
    cd ./AccuAmplicon
    mkdir -p seq temp result
    
    # 元数据细节优化
    # 转换Windows回车为Linux换行，去除空格
    sed -i 's/\r//;s/ //g' result/metadata.txt
    
    
# 1.拆分spike Split spike

    ## 01-原始序列处理
    # 启动工作环境
    conda activate meta_env

    # fastq.gz文件批量转为fasta文件
    mkdir -p temp/fasta_out
    for fq in seq/*.fastq.gz; do
        sample=$(basename "$fq" .fastq.gz)
        seqtk seq -A "$fq" > "temp/fasta_out/${sample}.fasta"
    done


    ## 02-并行blat比对
    fasta_dir=temp/fasta_out  # 输入 fasta 文件夹   
    spikein_db=./db/16S_SPIKEIN_with_10.fasta # spike-in 数据库
    blat_output=temp/blat_output # 输出文件夹
    BLAT=/usr/local/bin/blat
    mkdir -p "$blat_output"
    for fasta_file in "$fasta_dir"/*.fasta; do
        sample=$(basename "$fasta_file" .fasta)
        echo "Running BLAT for sample: $sample"
        # 执行比对并过滤比对长度 > 100 的结果
        $BLAT "$spikein_db" "$fasta_file" -out=blast8 -fastMap stdout \
            | awk '$4 > 100' > "$blat_output/${sample}_IR.map.xls"
    done

    # spike-in 数据库
    spiken_db=./db/16S_SPIKEIN_with_10.fasta
    blat_exact=temp/blat_exact

    mkdir -p "$blat_exact"

    # 提取 spike-in ID 并排序
    check_spike=$(grep ">" "$spiken_db" | awk '{print $1}' | sed 's/>//' | sort -V)

    # 为 R1 和 R2 创建统计文件并写入表头
    for Rtype in R1 R2; do
        out_file="$blat_exact/spikein_${Rtype}_stat.txt"
        echo -e "sample_id\tdirection\tspikein_id\tspikein_num" > "$out_file"
    done

    # 测试输出所有 spike-in ID
    echo "Spike-in IDs:"
    echo "$check_spike"
    
    
    ## 03-根据数据拆分策略loose/strict提取符合条件的序列标记为内参序列
    # 输入变量
    metadata="./result/metadata.txt"
    blat_output="temp/blat_output"                 # BLAT 输出目录
    blat_exact="temp/blat_exact"                   # 筛选后输出目录
    identity=96
    length=200

    mkdir -p "$blat_exact"
    
    # 读取 SampleID 列（跳过表头）
    samples=($(tail -n +2 "$metadata" | cut -f1))
    for sample in "${samples[@]}"; do
        for read_type in R1 R2; do
            input_blats="$blat_output/${sample}_${read_type}_IR.map.xls"
            select_blats="$blat_exact/select_${sample}_${read_type}_IR.map.xls"
            echo "Processing $sample $read_type ..."
            # 筛选符合条件的 reads
            awk -v id="$identity" -v len="$length" '
            {
                if ($3 >= id && $4 >= len) print $0
            }' "$input_blats" > "$select_blats"
            # 统计每个 spike-in 出现次数
            echo "$sample $read_type spike-in counts:"
            awk '{counts[$2]++} END{for(s in counts) print s, counts[s]}' "$select_blats"
        done
    done

    # 只提取 spike-in ID 列，循环处理每个样本
    for sample in "${samples[@]}"; do
        echo "Processing sample: $sample"
        awk '{print $1}' "$blat_output/${sample}_R1_IR.map.xls" \
            | sort \
            | uniq \
            > "$blat_exact/${sample}_spike_ids.txt"
    done


    ## 04-根据标记的内参序列对原数据进行拆分
    # 输入输出文件夹
    fastq_dir="./seq"
    client_final="temp/client_rawdata"
    spike_final="temp/spike_rawdata"
    blat_exact="temp/blat_exact"
    mkdir -p "$client_final" "$spike_final" "$blat_exact"

    # 只提取 spike-in ID 列，循环处理每个样本
    for sample in "${samples[@]}"; do
        spike_ids_file="$blat_exact/${sample}_spike_ids.txt"
        for read_type in R1 R2; do
            input_fastq="$fastq_dir/${sample}_${read_type}.fastq.gz"
            client_fastq="$client_final/${sample}_${read_type}.fastq"
            spike_fastq="$spike_final/${sample}_${read_type}.fastq"
            echo "Processing $sample $read_type ..."
            zcat "$input_fastq" | awk -v spike="$spike_ids_file" '
            BEGIN {
                while (getline < spike) { spikes[$1]=1 }
            }
            {
                if (NR%4==1) {
                    header=$0
                    read_id=$1
                    sub(/^@/, "", read_id)
                    sub(/\/[12]$/, "", read_id)
                    if (read_id in spikes) { out=2 } else { out=1 }
                } else if (NR%4==2) { seq=$0 }
                else if (NR%4==3) { plus=$0 }
                else if (NR%4==0) {
                    qual=$0
                    if (out==2) print header"\n"seq"\n"plus"\n"qual >> "'"$spike_fastq"'"
                    else print header"\n"seq"\n"plus"\n"qual >> "'"$client_fastq"'"
                }
            }'
            # 压缩输出文件
            gzip -f "$client_fastq" "$spike_fastq"
        done
    done
    
    
    ## 05-统计输出各样本R1/R2端 各内参reads 并检测有无添加量为0的内参
    # 输入文件
    blat_exact="temp/blat_exact"
    check_spike="$blat_exact/check_spike.txt" # 所有 spike-in ID列表，一行一个
    
    # 输出文件
    mkdir -p "$blat_exact"
    stat_R1="$blat_exact/spikein_R1_stat.txt"
    stat_R2="$blat_exact/spikein_R2_stat.txt"
    
    # 写表头
    echo -e "sample_id\tdirection\tspikein_id\tspikein_num" > "$stat_R1"
    echo -e "sample_id\tdirection\tspikein_id\tspikein_num" > "$stat_R2"
    
    # 主循环
    for sample in "${samples[@]}"; do
        echo "Processing sample: $sample"
        for Rtype in R1 R2; do
            select_file="$blat_exact/select_${sample}_${Rtype}_IR.map.xls"
            stat_file="$blat_exact/spikein_${Rtype}_stat.txt"
            if [ ! -s "$select_file" ]; then
                echo "文件不存在或为空: $select_file"
                continue
            fi
            echo "Counting spike-ins for $Rtype"
            # 统计每个 spike-in 的 reads 数量，并输出
            awk -v sample="$sample" -v Rtype="$Rtype" -v spikes="$check_spike" '
            BEGIN {
                # 读取 spike-in ID 列表
                while ((getline < spikes) > 0) {
                    gsub(/\r/, "", $1)
                    spike_ids[$1] = 0
                }
            }
            {
                # 第二列假设为 spike-in ID
                spike_ids[$2]++
            }
            END {
                for (id in spike_ids) {
                    n = spike_ids[id]
                    print sample "\t" Rtype "\t" id "\t" n
                    if (n == 0)
                        print "[ERROR] " sample "\t" Rtype "\t" id " spikein数量为 0" > "/dev/stderr"
                }
            }' "$select_file" >> "$stat_file"
        done
    done
    
    

# 2.获取ASV Get ASV

    # 启动运行环境
    conda activate qiime2-amplicon-2025.4
    qiime --help
    
    ## 2.1 Spike数据
    # 01-创建manifest文件
    # 输入与输出路径
    input_dir="temp/spike_rawdata"
    output_dir="temp/new_fq"
    manifest_file="temp/spike_rawdata/manifest"
    
    # 创建输出文件夹
    mkdir -p "${output_dir}"
    
    # 遍历输入文件
    for item in "${input_dir}"/*.fastq.gz; do
      base=$(basename "$item")
      echo "Processing: $base"
      # 提取样本ID与读向编号
      num=$(echo "$base" | sed -E 's/.*_([12])\.fastq\.gz/\1/')
      id=$(echo "$base" | sed -E 's/(SRR[0-9]*)_.*/\1/')
      # 根据编号复制并重命名
      if [[ "$num" == "1" ]]; then
        cp "$item" "${output_dir}/${id}"
      else
        cp "$item" "${output_dir}/${id}"
      fi
    done
    
    # 生成 manifest 文件
    > "$manifest_file"
    echo -e "sample-id\tforward-absolute-filepath\treverse-absolute-filepath" >> "$manifest_file"
    
    # 遍历样本生成 manifest 内容
    for id in $(ls "${output_dir}" | sed -E 's/_R[12]\.fastq\.gz//' | sort | uniq); do
      f_path="$(pwd)/${output_dir}/${id}_R1.fastq.gz"
      r_path="$(pwd)/${output_dir}/${id}_R2.fastq.gz"
      echo -e "${id}\t${f_path}\t${r_path}" >> "$manifest_file"
    done
    echo "Manifest file generated successfully: $(pwd)/${manifest_file}"
    
    
    # 02-qiime2导入和质控
    mkdir -p temp/qiime2_output/spike
    qiime tools import \
        --type 'SampleData[PairedEndSequencesWithQuality]' \
        --input-path temp/spike_rawdata/manifest \
        --output-path temp/qiime2_output/spike/demux.qza \
        --input-format PairedEndFastqManifestPhred33V2
    
    # 数据质量查看
    qiime demux summarize \
        --i-data temp/qiime2_output/spike/demux.qza \
        --o-visualization temp/qiime2_output/spike/demux-summary.qzv
    
    
    # 03-qiime2去除接头和引物
    # 正向
    qiime cutadapt trim-paired \
        --i-demultiplexed-sequences temp/qiime2_output/spike/demux.qza \
        --p-error-rate 0.15 \
        --p-cores 4 \
        --p-no-indels \
        --p-adapter-f CTGTCTCTTATACACATCTGACGCTGCCGACGA \
        --p-front-f GACTACHVGGGTATCTAATCC \
        --p-adapter-r CTGTCTCTTATACACATCTCCGAGCCCACGAGAC \
        --p-front-r GGACTACHVGGGTWTCTAAT  \
        --o-trimmed-sequences temp/qiime2_output/spike/primer-trimmed-demux.qza
    
    # 反向
    qiime cutadapt trim-paired \
        --i-demultiplexed-sequences temp/qiime2_output/spike/primer-trimmed-demux.qza \
        --p-error-rate 0.15 \
        --p-cores 4 \
        --p-no-indels \
        --p-adapter-f CTGTCTCTTATACACATCTCCGAGCCCACGAGAC \
        --p-front-f GACTACHVGGGTATCTAATCC \
        --p-adapter-r CTGTCTCTTATACACATCTGACGCTGCCGACGA \
        --p-front-r GGACTACHVGGGTWTCTAAT  \
        --o-trimmed-sequences temp/qiime2_output/spike/primer-trimmed-demux2.qza
    
    # 查看结果
    qiime demux summarize \
        --i-data temp/qiime2_output/spike/primer-trimmed-demux2.qza \
        --o-visualization temp/qiime2_output/spike/primer-trimmed-demux.qzv
    
    
    # 04-qiime2中dada2去噪、去嵌合及特征表构建
    # DADA2 参数（根据实验调整）
    # trim_left_f=0   # 正向前端截断长度
    # trim_left_r=0   # 反向前端截断长度
    # trunc_len_f=223   # 正向 trunc 长度
    # trunc_len_r=223   # 反向 trunc 长度
    qiime dada2 denoise-paired \
        --i-demultiplexed-seqs temp/qiime2_output/spike/primer-trimmed-demux2.qza \
        --p-n-threads 2 \
        --p-trim-left-f 0 --p-trim-left-r 0 \
        --p-trunc-len-f 223 --p-trunc-len-r 223 \
        --o-table temp/qiime2_output/spike/dada2-table2.qza \
        --o-representative-sequences temp/qiime2_output/spike/dada2-rep-seqs2.qza \
        --o-denoising-stats temp/qiime2_output/spike/denoising-stats2.qza
        
    # 导出 DADA2 结果
    # DADA2 输出路径
    output="result/dada2"
    dada2_path="$output/5.dada2_spike"
    mkdir -p $dada2_path
    rep_seqs="$dada2_path/dada2-rep-seqs2.qza"
    feature_table="$dada2_path/dada2-table2.qza"
    denoise_stats="$dada2_path/denoising-stats2.qza"
    
    cp ./temp/qiime2_output/spike/dada2-rep-seqs2.qza $dada2_path/dada2-rep-seqs2.qza
    cp ./temp/qiime2_output/spike/dada2-table2.qza $dada2_path/dada2-table2.qza
    cp ./temp/qiime2_output/spike/denoising-stats2.qza $dada2_path/denoising-stats2.qza
    
    qiime tools export --input-path "$denoise_stats" --output-path "$dada2_path/denoise_stats"
    qiime tools export --input-path "$rep_seqs" --output-path "$dada2_path/rep_seqs"
    qiime tools export --input-path "$feature_table" --output-path "$dada2_path/feature_table"
    
    # BIOM 转 TSV
    biom convert \
      -i "$dada2_path/feature_table/feature-table.biom" \
      -o "$dada2_path/feature_table/feature-table.txt" \
      --to-tsv
    
    echo "DADA2 denoise & export 完成"
    echo "输出目录: $dada2_path"
    echo "feature-table.txt: $dada2_path/feature-table.txt"
    
    
    # 05-qiime2结果整理
    # 定义输出路径
    output="./result"   #请改成实际路径
    final_path="${output}/qiime2/spike/6.final"
    dada2_path="${output}/dada2"  # 假设 DADA2 结果在这里，请按实际情况修改
    
    # 创建 final 目录
    mkdir -p "${final_path}"
    
    # 拷贝结果文件
    cp -rf "${dada2_path}/5.dada2_spike/rep_seqs/dna-sequences.fasta" "${final_path}/"
    cp -rf "${dada2_path}/5.dada2_spike/feature_table/feature-table.txt" "${final_path}/"
    
    # 提取统计信息（去掉包含 'numeric' 的行）
    grep -v "numeric" "${dada2_path}/5.dada2_spike/denoise_stats/stats.tsv" > "${final_path}/seqStat.txt"
    
    
    
    ## 2.2 client数据
    # 01-创建manifest文件
    # 输入与输出路径
    input_dir="temp/client_rawdata"
    output_dir="temp/new_fq_client"
    manifest_file="temp/client_rawdata/manifest_client"
    
    # 创建输出文件夹
    mkdir -p "${output_dir}"
    
    # 遍历输入文件
    for item in "${input_dir}"/*.fastq.gz; do
      base=$(basename "$item")
      echo "Processing: $base"
      # 提取样本ID与读向编号
      num=$(echo "$base" | sed -E 's/.*_([12])\.fastq\.gz/\1/')
      id=$(echo "$base" | sed -E 's/(SRR[0-9]*)_.*/\1/')
      # 根据编号复制并重命名
      if [[ "$num" == "1" ]]; then
        cp "$item" "${output_dir}/${id}"
      else
        cp "$item" "${output_dir}/${id}"
      fi
    done
    
    # 生成 manifest 文件
    > "$manifest_file"
    echo -e "sample-id\tforward-absolute-filepath\treverse-absolute-filepath" >> "$manifest_file"
    
    # 遍历样本生成 manifest 内容
    for id in $(ls "${output_dir}" | sed -E 's/_R[12]\.fastq\.gz//' | sort | uniq); do
      f_path="$(pwd)/${output_dir}/${id}_R1.fastq.gz"
      r_path="$(pwd)/${output_dir}/${id}_R2.fastq.gz"
      echo -e "${id}\t${f_path}\t${r_path}" >> "$manifest_file"
    done
    echo "Manifest file generated successfully: $(pwd)/${manifest_file}"
    
    
    # 02-qiime2导入和质控
    mkdir -p temp/qiime2_output/client
    qiime tools import \
        --type 'SampleData[PairedEndSequencesWithQuality]' \
        --input-path temp/client_rawdata/manifest_client \
        --output-path temp/qiime2_output/client/demux_client.qza \
        --input-format PairedEndFastqManifestPhred33V2
    
    # 数据质量查看
    qiime demux summarize \
        --i-data temp/qiime2_output/client/demux_client.qza \
        --o-visualization temp/qiime2_output/client/demux-summary-client.qzv
    
    
    # 03-qiime2去除接头和引物
    # 正向
    qiime cutadapt trim-paired \
        --i-demultiplexed-sequences temp/qiime2_output/client/demux_client.qza \
        --p-error-rate 0.15 \
        --p-cores 4 \
        --p-no-indels \
        --p-adapter-f CTGTCTCTTATACACATCTGACGCTGCCGACGA \
        --p-front-f GACTACHVGGGTATCTAATCC \
        --p-adapter-r CTGTCTCTTATACACATCTCCGAGCCCACGAGAC \
        --p-front-r GGACTACHVGGGTWTCTAAT  \
        --o-trimmed-sequences temp/qiime2_output/client/primer-trimmed-demux-client.qza
    
    # 反向
    qiime cutadapt trim-paired \
        --i-demultiplexed-sequences temp/qiime2_output/client/primer-trimmed-demux-client.qza \
        --p-error-rate 0.15 \
        --p-cores 4 \
        --p-no-indels \
        --p-adapter-f CTGTCTCTTATACACATCTCCGAGCCCACGAGAC \
        --p-front-f GACTACHVGGGTATCTAATCC \
        --p-adapter-r CTGTCTCTTATACACATCTGACGCTGCCGACGA \
        --p-front-r GGACTACHVGGGTWTCTAAT  \
        --o-trimmed-sequences temp/qiime2_output/client/primer-trimmed-demux2-client.qza
    
    # 查看结果
    qiime demux summarize \
        --i-data temp/qiime2_output/client/primer-trimmed-demux2-client.qza \
        --o-visualization temp/qiime2_output/client/primer-trimmed-demux-client.qzv
    
    
    # 04-qiime2中dada2去噪、去嵌合及特征表构建
    # DADA2 参数（根据实验调整）
    # trim_left_f=0   # 正向前端截断长度
    # trim_left_r=0   # 反向前端截断长度
    # trunc_len_f=223   # 正向 trunc 长度
    # trunc_len_r=223   # 反向 trunc 长度
    qiime dada2 denoise-paired \
        --i-demultiplexed-seqs temp/qiime2_output/client/primer-trimmed-demux2-client.qza \
        --p-n-threads 2 \
        --p-trim-left-f 0 --p-trim-left-r 0 \
        --p-trunc-len-f 223 --p-trunc-len-r 223 \
        --o-table temp/qiime2_output/client/dada2-table2-client.qza \
        --o-representative-sequences temp/qiime2_output/client/dada2-rep-seqs2-client.qza \
        --o-denoising-stats temp/qiime2_output/client/denoising-stats2-client.qza
        
    # 导出 DADA2 结果
    # DADA2 输出路径
    #dada2_path="./output/5.dada2_client"
    #rep_seqs="$dada2_path/dada2-rep-seqs2-client.qza"
    #feature_table="$dada2_path/dada2-table2-client.qza"
    #denoise_stats="$dada2_path/denoising-stats2-client.qza"
    
    output="result/dada2"
    dada2_path="$output/5.dada2_client"
    mkdir -p $dada2_path
    rep_seqs="$dada2_path/dada2-rep-seqs2-client.qza"
    feature_table="$dada2_path/dada2-table2-client.qza"
    denoise_stats="$dada2_path/denoising-stats2-client.qza"
    
    cp ./temp/qiime2_output/client/dada2-rep-seqs2-client.qza $dada2_path/dada2-rep-seqs2-client.qza
    cp ./temp/qiime2_output/client/dada2-table2-client.qza $dada2_path/dada2-table2-client.qza
    cp ./temp/qiime2_output/client/denoising-stats2-client.qza $dada2_path/denoising-stats2-client.qza
    
    qiime tools export --input-path "$denoise_stats" --output-path "$dada2_path/denoise_stats"
    qiime tools export --input-path "$rep_seqs" --output-path "$dada2_path/rep_seqs"
    qiime tools export --input-path "$feature_table" --output-path "$dada2_path/feature_table"
    
    # BIOM 转 TSV
    biom convert \
      -i "$dada2_path/feature_table/feature-table.biom" \
      -o "$dada2_path/feature_table/feature-table.txt" \
      --to-tsv
    
    echo "DADA2 denoise & export 完成"
    echo "输出目录: $dada2_path"
    echo "feature-table.txt: $dada2_path/feature-table.txt"
    
    
    # 05-qiime2结果整理
    # 定义输出路径
    output="./result"   #请改成实际路径
    final_path="${output}/qiime2/client/6.final"
    dada2_path="${output}/dada2"  #假设 DADA2 结果在这里，请按实际情况修改
    
    # 创建 final 目录
    mkdir -p "${final_path}"
    
    # 拷贝结果文件
    cp -rf "${dada2_path}/5.dada2_client/rep_seqs/dna-sequences.fasta" "${final_path}/"
    cp -rf "${dada2_path}/5.dada2_client/feature_table/feature-table.txt" "${final_path}/"
    
    # 提取统计信息（去掉包含 'numeric' 的行）
    grep -v "numeric" "${dada2_path}/5.dada2_client/denoise_stats/stats.tsv" > "${final_path}/seqStat.txt"



# 3.绝对定量

    # 定义输出路径（请替换为你的实际路径）
    output="./result/qiime2"  # 输出主目录
    final="${output}/6.final" # 最终结果目录
    SPIKEIN_INFO="./db/spikein_filter.txt"  # spike-in 信息文件（如果后面要用）
    amplicon_area="16S" # 举例，可按实际修改
    f_primer="ACTCCTACGGGAGGCAGCAG" # 正向引物
    r_primer="GGACTACHVGGGTWTCTAAT" # 反向引物
    input_client="./result/qiime2/client" # 客户输入目录
    input_spike="./result/qiime2/spike" # spike-in 输入目录
    
    # 创建目录 
    mkdir -p "$final"
    
    # 定义文件路径
    client_otu_txt="${input_client}/6.final/feature-table.txt"
    client_rep_fasta="${input_client}/6.final/dna-sequences.fasta"
    spike_fasta="${input_spike}/6.final/dna-sequences.fasta"
    spike_otu_txt="${input_spike}/6.final/feature-table.txt"
    
    ## 01-根据样本名提取 ir_copies及dna_quantity
    # 提取 ir_copies
    # 在去之前需要检查feature_table.txt文件，确保文件正确，这里是去掉了feature_table.txt文件的第一行
    otu_file="${input_client}/6.final/feature-table2.txt"
    raw_ir_copies="./result/IR_copies2.txt"
    target_ir_copies="${output}/6.final/0.spike_IR_copies.txt"
    
    # 提取 OTU 表样本名（按 tab 安全读取，保留原始顺序）
    IFS=$'\t' read -r -a _header <<< "$(head -1 "$otu_file")"
    # _header[0] 是第一列（OTU ID），其后为样本名
    samples=("${_header[@]:1}")
    
    # 用 awk 提取 IR_Copies（显式 tab 分割，去除 \r）
    awk -v samples="${samples[*]}" -F $'\t' -v OFS=$'\t' '
    BEGIN{
        # 把 samples 字符串按空格拆成数组（样本名中含空格的情形已被上面处理）
        n = split(samples, s_arr, " ");
        for(i=1;i<=n;i++){
            gsub(/\r/,"",s_arr[i])
        }
    }
    NR==1{
        # 去掉 header 中的 CR 字符并建立映射
        for(i=1;i<=NF;i++){
            key=$i
            gsub(/\r/,"",key)
            header[key]=i
        }
        next
    }
    NR>1{
        # 处理每一行：第一列作为 otu
        gsub(/\r/,"",$1)
        otu=$1
        printf "%s", otu
        for(i=1;i<=n;i++){
            sample=s_arr[i]
            col = header[sample]
            if(!col){
                # 尝试带前缀的名字
                col = header["GS_YD_" sample]
            }
            if(!col){
                # 如果仍然找不到，输出 NA 并在 stderr 打印警告（便于调试）
                printf "\tNA"
                # 把哪个样本没找到输出到 stderr（只在第一次遇到时输出）
                if(!warned[sample]++){
                    print "WARNING: sample not found in raw_ir_copies header -> " sample > "/dev/stderr"
                }
            } else {
                printf "\t%s", $col
            }
        }
        print ""
    }' "$raw_ir_copies" > "$target_ir_copies"
    
    # 添加表头（保持 tab）
    {
        echo -e "#OTU ID\t$(head -1 "$otu_file" | cut -f2-)"
        cat "$target_ir_copies"
    } > tmp && mv tmp "$target_ir_copies"
    echo "[INFO] IR copies 提取完成：$target_ir_copies"
    
    
    # 提取dna_quantity
    otu_file="${input_client}/6.final/feature-table2.txt"
    raw_dna_quantity="./result/DNA_quantity2.txt"
    target_dna_quantity="${output}/6.final/0.DNA_quantity.txt"
    
    # 获取样本名（严格按制表符读取）
    IFS=$'\t' read -r -a header <<< "$(head -1 "$otu_file")"
    samples=("${header[@]:1}")   # 跳过第一列 OTU_ID
    
    # 生成 DNA_quantity 新文件
    {
        # 打印原始表头
        head -1 "$raw_dna_quantity"
        
        # 用 awk 按制表符匹配样本
        awk -F '\t' -v OFS='\t' -v samples="${samples[*]}" '
        BEGIN {
            n = split(samples, s_arr, " ")
            for (i = 1; i <= n; i++) {
                sample = s_arr[i]
                gsub(/\r/, "", sample)
                sample_set[sample] = 1
            }
        }
        NR > 1 {
            gsub(/\r/, "", $1)
            otu = $1
            if (otu in sample_set) {
                print $0
            } else {
                # 尝试带前缀 GS_YD_
                gsyd = "GS_YD_" otu
                if (gsyd in sample_set) {
                    $1 = gsyd
                    print $0
                }
            }
        }' "$raw_dna_quantity"
    
    } > "$target_dna_quantity"
    echo "[INFO] DNA quantity 提取完成：$target_dna_quantity"
    
    
    ## 02-blastn 两两比对，找到数据库中与输入序列最相似的序列
    # blastn安装
    # 从网页下载然后解压后安装成功
    # wget ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/ncbi-blast-2.17.0+-x64-linux.tar.gz
    # tar -zxvf ncbi-blast-2.17.0+-x64-linux.tar.gz
    export PATH=$PATH:/mnt/f/001_absolute_quantification/Accu16S_ITS/ncbi-blast-2.17.0+/bin
    source ~/.bashrc
    blastn -version
    
    # 定义变量
    output="./result/qiime2/6.final"
    BLASTN="/mnt/f/001_absolute_quantification/Accu16S_ITS/ncbi-blast-2.17.0+/bin/blastn"    # blastn 完整路径
    MAKEBLASTDB="/mnt/f/001_absolute_quantification/Accu16S_ITS/ncbi-blast-2.17.0+/bin/makeblastdb" # makeblastdb 完整路径
    spike_fasta="./result/qiime2/spike/6.final/dna-sequences.fasta"
    SPIKEIN_DB="./db/16S_SPIKEIN.fasta"
    THREAD=8   # 根据CPU线程数设置
    
    # 创建输出目录
    mkdir -p "$output"
    
    # 判断数据库是否已建立
    if [ ! -f "${SPIKEIN_DB}.nhr" ]; then
        echo "[INFO]数据库索引不存在，开始建立数据库..."
        $MAKEBLASTDB -in "${SPIKEIN_DB}.fasta" -dbtype nucl -out "$SPIKEIN_DB"
        if [ $? -ne 0 ]; then
            echo "[ERROR]数据库建立失败，请检查 $SPIKEIN_DB.fasta 是否存在"
            exit 1
        fi
    fi
    
    # 运行 blastn
    echo "[INFO] 开始运行 blastn..."
    $BLASTN \
      -query "$spike_fasta" \
      -db "$SPIKEIN_DB" \
      -evalue 1e-50 \
      -outfmt 6 \
      -max_target_seqs 100000000 \
      -perc_identity 90 \
      -num_threads "$THREAD" \
      > "$output/1.spike.blast"
    
    if [ $? -eq 0 ]; then
        echo "[INFO] blastn 运行完成，结果保存在 $output/1.spike.blast"
    else
        echo "[ERROR] blastn 运行失败"
        exit 1
    fi
    
    # 生成pair_id.txt
    blastn_result="$output/1.spike.blast"         # 输入：BLASTN 比对结果文件
    spikein_info="./db/spikein_filter.txt"       # 输入：内参对应信息
    pair_id="./result/qiime2/6.final/2.pair_id.txt"  # 输出：OTU 与 GS 对应结果
    # 输出文件表头
    echo -e "OTU_ID\tGs_BSI_ID" > "$pair_id"
    # awk 匹配并增加筛选条件
    awk -F'\t' '
    FNR==NR {spike[$1]=1; next}                  # 先读取 spikein_info 的第一列
    ($2 in spike) && ($3 >= 96) && ($4 >= 400) { # 条件：第二列匹配 & 第3列>=96 & 第4列>=200
        print $1"\t"$2
    }' "$spikein_info" "$blastn_result" >> "$pair_id"


    # 整理统计blastn结果
    # 进入R语言界面
    R
    getwd()
    
    # 配置文件路径
    output <- "./result/qiime2/6.final"
    blastn_result <- file.path(output, "1.spike.blast")
    spike_otu_txt <- "./result/qiime2/spike/6.final/feature-table2.txt"
    
    pair_id_file <- file.path(output, "2.pair_id.txt")
    match_spike_otu <- file.path(output, "2.spike_otu.txt")
    mismatch_spike_otu <- file.path(output, "2.mismatch_spike_otu.txt")
    sample_sum_file <- file.path(output, "2.sample_sum.txt")
    
    # 加载软件包
    library(data.table)
    
    # 读取数据
    pair_id <- fread(pair_id_file, header = TRUE, colClasses = "character")
    spike_otu <- fread(spike_otu_txt, header = TRUE, colClasses = "character")
    
    # 去掉 OTU_ID 前后空格
    pair_id$OTU_ID <- trimws(pair_id$OTU_ID)
    spike_otu[[1]] <- trimws(spike_otu[[1]])
    
    # 设置行名并删除第一列
    rownames(spike_otu) <- spike_otu[[1]]
    
    # 匹配 OTU_ID
    matched_ids <- pair_id$OTU_ID
    cat("[INFO] 总 OTU 数量:", nrow(spike_otu), "\n")
    cat("[INFO] 匹配 ID 数量:", length(matched_ids), "\n")
    cat("[INFO] 匹配 ID 示例:", head(matched_ids), "\n")
    
    # 实际匹配
    match_spike <- spike_otu[rownames(spike_otu) %in% matched_ids, , drop = FALSE]
    cat("[INFO] 实际匹配行数:", nrow(match_spike), "\n")
    if (nrow(match_spike) == 0) {
      unmatched <- setdiff(matched_ids, rownames(spike_otu))
      cat("[WARN] 以下 OTU_ID 没匹配上:\n")
      print(unmatched)
      stop("[ERROR] 没有找到匹配的 OTU，请检查 OTU_ID\n")
    }
    
    # 添加 OTU_ID 列
    colnames(match_spike)[1] = "OTU"
    pair_id = as.data.frame(pair_id)
    colnames(pair_id)[1] <- "OTU"
    
    # 合并 pair_id 信息
    match_spike <- merge(pair_id, match_spike, by = "OTU", all.x = FALSE)
    match_spike <- as.data.frame(match_spike)
    
    # 聚合：合并相同 Gs_BSI_ID 的行并求和
    samples <- colnames(spike_otu)[2:4]
    data2 <- match_spike[, c("Gs_BSI_ID", samples)]
    
    # 假设 data2 第一列是 Gs_BSI_ID，其余列需要聚合
    numeric_cols <- setdiff(names(data2), "Gs_BSI_ID")
    data2[numeric_cols] <- lapply(data2[numeric_cols], as.numeric)
    match_spike_agg <- aggregate(. ~ Gs_BSI_ID, data = data2, FUN = sum)
    write.table(match_spike_agg, file = match_spike_otu, sep = "\t", quote = FALSE, row.names = FALSE)
    
    # mismatch 部分
    all_ids <- rownames(spike_otu)
    mismatch_ids <- setdiff(all_ids, matched_ids)
    mismatch_spike <- spike_otu[rownames(spike_otu) %in% mismatch_ids, , drop = FALSE]
    mismatch_spike <- cbind(OTU_ID = rownames(mismatch_spike), mismatch_spike)
    write.table(mismatch_spike, file = mismatch_spike_otu, sep = "\t", quote = FALSE, row.names = FALSE)
    
    # 样本统计
    spike_otu2 <- spike_otu[, -1]
    #spike_otu2 <- as.numeric(spike_otu2)
    spike_otu2 <- as.data.frame(lapply(spike_otu2, as.numeric))
    total_reads <- colSums(spike_otu2)
    matched_reads <- colSums(match_spike_agg[, samples, drop = FALSE])
    percent_match <- matched_reads * 100 / total_reads
    
    sample_sum_df <- data.frame(
      Sample = samples,
      Totle_spikein_reads = total_reads,
      Match_spikein_reads = matched_reads,
      Match_spikein_percent = percent_match
    )
    write.table(sample_sum_df, file = sample_sum_file, sep = "\t", quote = FALSE, row.names = FALSE)
    cat("[INFO] ir_detection 分析完成:", sample_sum_file, "\n")
    
    # Ctrl + A + D 退出R语言环境
    
    
    ## 03-根据内参添加量，计算拟合曲线  从而计算绝对拷贝数  DNA水平绝对拷贝数  SAMPLE 水平绝对拷贝数
    # 根据内参添加量计算绝对拷贝数
    # 输入与输出路径定义
    output="./result/qiime2/6.final"
    otu_txt="./result/qiime2/spike/6.final/feature-table2.txt"  # OTU 表（可根据实际路径修改）
    match_spike_otu="$output/2.spike_otu.txt" # 匹配的 spike OTU 表
    spike_IR_copies="$output/0.spike_IR_copies.txt" # 内参拷贝数
    DNA_quantity="$output/0.DNA_quantity.txt" # DNA 浓度信息
    
    # 输出文件路径（供下游步骤使用）
    spike_reads_percent="$output/3.spike_reads_percent.txt"
    standard_curve_formula="$output/3.standard_curve_formula.txt"
    absolute_otu_txt="$output/3.absolute_otu.txt"
    unit_dna_otu_copies="$output/3.unit_dna_otu_copies.txt"
    unit_sample_otu_copies="$output/3.unit_sample_otu_copies.txt"
    
    # 计算脚本路径（请修改为你真实路径）
    CALCULATE_ABSOLUTE_COPIES="./script/quantitation_calculate_absolute_copies.R"
    
    # 执行 R 脚本计算绝对拷贝数
    Rscript "$CALCULATE_ABSOLUTE_COPIES" \
        -o "$otu_txt" \
        -s "$match_spike_otu" \
        -i "$spike_IR_copies" \
        -d "$DNA_quantity" \
        -p 3 \
        --output "$output"
    
    # 检查标准曲线文件是否生成
    if [[ ! -s "$standard_curve_formula" ]]; then
        echo "[ERROR] 未生成标准曲线文件：$standard_curve_formula" >&2
        echo "请检查 R 脚本 calculate_absolute_copies.R 是否执行成功。" >&2
        exit 1
    else
        echo "[INFO] 标准曲线文件检查通过：$standard_curve_formula"
    fi
    
    
    ##04- 整理最终结果
    # 定义路径变量（根据你的项目路径修改）
    output="./result/qiime2/6.final"
    final="result/6.final"
    otu_txt="./result/qiime2/spike/6.final/feature-table2.txt"
    rep_fasta="./result/qiime2/spike/6.final/dna-sequences.fasta"
    match_spike_otu="${output}/2.spike_otu.txt"
    absolute_otu_txt="${output}/3.absolute_otu.txt"
    standard_curve_formula="${output}/3.standard_curve_formula.txt"
    spike_reads_percent="${output}/3.spike_reads_percent.txt"
    unit_dna_otu_copies="${output}/3.unit_dna_otu_copies.txt"
    unit_sample_otu_copies="${output}/3.unit_sample_otu_copies.txt"
    
    # 创建输出目录
    mkdir -p "$final"
    
    # 定义复制函数（带错误检查）
    copy_file() {
        local src=$1
        local dest=$2
        if [[ -e "$src" ]]; then
            cp -rf "$src" "$dest" || {
                echo "[ERROR] 复制失败：$src → $dest" >&2
                exit 1
            }
        else
            echo "[WARNING] 未找到文件：$src" >&2
        fi
    }
    
    # 开始复制文件
    copy_file "$otu_txt" "$final/"
    copy_file "$rep_fasta" "$final/"
    copy_file "$match_spike_otu" "$final/spike_otu.txt"
    copy_file "$absolute_otu_txt" "$final/absolute_otu.txt"
    copy_file "$standard_curve_formula" "$final/standard_curve_formula.txt"
    copy_file "$spike_reads_percent" "$final/spike_reads_percent.txt"
    copy_file "$unit_dna_otu_copies" "$final/unit_dna_otu_copies.txt"
    
    # 仅当存在时才复制 unit_sample_otu_copies
    if [[ -e "$unit_sample_otu_copies" ]]; then
        copy_file "$unit_sample_otu_copies" "$final/unit_sample_otu_copies.txt"
    fi
    
    # 打印结果提示
    echo ""
    echo "[OK] 所有结果已整理至目录：$final"



# 4. 最终结果
    
    # 最终结果路径
    otu_table=./result/6.final/feature-table2.txt
    otu_fasta=./result/6.final/dna-sequences.fasta
    formula=./result/6.final/standard_curve_formula.txt
    absolute_otu=./result/6.final/absolute_otu.txt
    absolute_otu_dna=./result/6.final/unit_dna_otu_copies.txt
    absolute_otu_sample=./result/6.final/unit_sample_otu_copies.txt
    spike_table=./result/6.final/spike_otu.txt
    spike_percent=./result/6.final/spike_reads_percent.txt
    
    echo "ASV 丰度表： $otu_table"
    echo "ASV 序列： $otu_fasta"
    echo "拟合曲线公式： $formula"
    echo "绝对拷贝数： $absolute_otu"
    echo "DNA水平单位绝对拷贝数： $absolute_otu_dna"
    echo "样本水平单位绝对拷贝数： $absolute_otu_sample"
    echo "spike 丰度表： $spike_table"
    echo "spike 在下机数据中的含量： $spike_percent"

