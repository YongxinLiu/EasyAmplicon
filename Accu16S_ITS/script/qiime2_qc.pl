$| = 1;
# Perl 系统包
use warnings;
use strict;
use File::Spec;
use Getopt::Long;
use HTML::TableExtract;
use Parallel::ForkManager;
use List::Util qw/sum/;
use Sort::Key::Natural qw(natsort);
use Bio::SeqIO;
use Cwd 'abs_path';

# 常量
use constant SCRIPT_DIR => (File::Spec->splitpath(File::Spec->rel2abs($0)))[1];
use constant PWD => $ENV{"PWD"};

# 变量
my $id = `id`;
my ($uid, $gid) = $id =~ /uid=(\d+).+?gid=(\d+)/;

my $DEFAULT_THREAD = 40;
my $DEFAULT_MAXEE  = 2;
my $DEFAULT_TRUNCQ = 2;
my $DEFAULT_QIIME2               = "export TZ='Pacific/Auckland' &&  source /home/genesky/software/conda/4.8.3/bin/activate /home/genesky/software/qiime2/2023.2";


# 参数输入
my $ARGV_INFO = join " ", @ARGV;
my ( $config, $input, $output, $samples, $THREAD, $QIIME2, $if_help);
GetOptions(
    "config|c=s"         => \$config,
    "input|i=s"          => \$input,
    "output|o=s"         => \$output,
    "samples|s=s"        => \$samples,
    "thread=i"           => \$THREAD,

    "qiime2=s"           => \$QIIME2,
    "help|h"             => \$if_help,
);
die "
Options: 必填

        --config/-c            参数配置文件
                               F_Primer:正向测序引物,必须填写
                               R_Primer:反向测序引物,必须填写
                               p-trunc-len-f: 自定义 QIIME2 质控时R1序列截取保留的长度
                               p-trunc-len-r: 自定义 QIIME2 质控时R2序列截取保留的长度
                               p-trim-left-f: 自定义 QIIME2 质控时R1序列左侧截取开始位置
                               p-trim-left-r: 自定义 QIIME2 质控时R2序列左侧截取开始位置
        --input/-i             fastq数据路径，包含 sample_R1.fastq.gz sample_R1.fastq.gz
        --output/-o            结果输出路径
        --samples/-s           选定要分析的样本名称，多个样本用逗号分隔，例如 'sample1,sample2'

Options: 可选

        --thread               使用线程数 (default: $DEFAULT_THREAD)

        --qiime2               更改软件 qiime2 版本 (default: '$DEFAULT_QIIME2')
        --help/-h              查看帮助文档
\n" if (defined $if_help or not defined $config or not defined $input or not defined $output or not defined $samples);
system qq{mkdir -p $output} if not -d $output;

$input  = abs_path($input);
$output = abs_path($output);
$THREAD               = $DEFAULT_THREAD if (not defined $THREAD);
$QIIME2               = $DEFAULT_QIIME2        if (not defined $QIIME2);

my @samples = split /,/, $samples;
my @parameter     = ("F_Primer","R_Primer");
my %config        = read_config($config, \@parameter);
my $f_primer      = $config{"F_Primer"}{"value"};
my $r_primer      = $config{"R_Primer"}{"value"};
my $trunc_len_f   = exists $config{"p-trunc-len-f"} ? $config{"p-trunc-len-f"}{"value"} : "";
my $trunc_len_r   = exists $config{"p-trunc-len-r"} ? $config{"p-trunc-len-r"}{"value"} : "";
my $trim_left_f   = exists $config{"p-trim-left-f"} ? $config{"p-trim-left-f"}{"value"} : 0;
my $trim_left_r   = exists $config{"p-trim-left-r"} ? $config{"p-trim-left-r"}{"value"} : 0;
$f_primer .= "N"; # 去除引物后的第一个碱基（测序错误率较高）
$r_primer .= "N";
###################################################################### 初始化


#######################################################################主程序

my $qiime2_tmp = "$output/tmp_qiime2";
system qq{rm -rf $qiime2_tmp} if -d $qiime2_tmp;
system qq{mkdir -p $qiime2_tmp};

###1.准备原始数据文件路径
my $manifest_raw = "$output/1.rawdata-manifest.txt";
prepare_rawdata_manifest($input, $manifest_raw, \@samples);

###1.原始数据转格式
my $raw_qiime2 = "$output/1.raw_data";
system qq{$QIIME2 && export TMPDIR=$qiime2_tmp && qiime tools import --type 'SampleData[PairedEndSequencesWithQuality]' --input-path $manifest_raw --input-format PairedEndFastqManifestPhred33 --output-path $raw_qiime2};

###2.去接头和引物---正向
my $trim_forward        = "$output/2.trim_adapter_and_primer_forward";
my $trim_forward_export = "$output/2.trim_adapter_and_primer_forward_export";
system qq{$QIIME2 && export TMPDIR=$qiime2_tmp && qiime cutadapt trim-paired --p-error-rate 0.15 --p-cores $THREAD --p-discard-untrimmed --i-demultiplexed-sequences $raw_qiime2.qza --p-adapter-f CTGTCTCTTATACACATCTGACGCTGCCGACGA --p-front-f $f_primer --p-adapter-r CTGTCTCTTATACACATCTCCGAGCCCACGAGAC --p-front-r $r_primer --o-trimmed-sequences $trim_forward};
system qq{$QIIME2 && export TMPDIR=$qiime2_tmp && qiime tools export --input-path $trim_forward.qza --output-path $trim_forward_export};

###3.去接头和引物---反向
my $trim_reverse        = "$output/3.trim_adapter_and_primer_reverse";
my $trim_reverse_export = "$output/3.trim_adapter_and_primer_reverse_export";
system qq{$QIIME2 && export TMPDIR=$qiime2_tmp && qiime cutadapt trim-paired --p-error-rate 0.15 --p-cores $THREAD --p-discard-untrimmed --i-demultiplexed-sequences $raw_qiime2.qza --p-adapter-f CTGTCTCTTATACACATCTCCGAGCCCACGAGAC --p-front-f $r_primer --p-adapter-r CTGTCTCTTATACACATCTGACGCTGCCGACGA --p-front-r $f_primer --o-trimmed-sequences $trim_reverse};
system qq{$QIIME2 && export TMPDIR=$qiime2_tmp && qiime tools export --input-path $trim_reverse.qza --output-path $trim_reverse_export};

###4.重新导入去引物和接头后序列
my $trimmed_fastq = "$output/4.trimmed_adapter_primer_fastq";
my $trimmed       = "$output/4.trimmed_adapter_primer";
my $trimmed_summarize         = "$output/4.trimmed_adapter_primer_summarize";
my $trimmed_summarize_export  = "$output/4.trimmed_adapter_primer_summarize_export";
system qq{mkdir -p $trimmed_fastq} if not -d $trimmed_fastq;

combine_raw_data($trim_forward_export, $trim_reverse_export, $trimmed_fastq, $input, $THREAD);  # 合并正向切除引物和反向切除引物fastq, 并计算各样本引物检测成功序列比例


system qq{$QIIME2 && export TMPDIR=$qiime2_tmp && qiime tools import --type 'SampleData[PairedEndSequencesWithQuality]' --input-path $trimmed_fastq/NO-AP-MANIFEST --input-format PairedEndFastqManifestPhred33 --output-path $trimmed};
system qq{$QIIME2 && export TMPDIR=$qiime2_tmp && qiime demux summarize --i-data $trimmed.qza --o-visualization $trimmed_summarize};
system qq{$QIIME2 && export TMPDIR=$qiime2_tmp && qiime tools export --input-path $trimmed_summarize.qzv --output-path $trimmed_summarize_export};


###5.DADA2 fastq质量过滤、merge、去嵌合体、聚ASV、ASV生成丰度表
my $dada2_path    = "$output/5.dada2";
my $rep_seqs      = "$dada2_path/rep_seqs";
my $feature_table = "$dada2_path/feature_table";
my $denoise_stats = "$dada2_path/denoise_stats";

system qq{mkdir -p $dada2_path} if not -d $dada2_path;
system qq{rm -rf $dada2_path/feature-table.biom} if -e "$dada2_path/feature-table.biom";  # 对于已存在的文件,biom convert无法重定向

system qq{$QIIME2 && export TMPDIR=$qiime2_tmp && qiime dada2 denoise-paired --p-n-threads $THREAD  --i-demultiplexed-seqs $trimmed.qza --p-trim-left-f $trim_left_f --p-trim-left-r $trim_left_r --p-trunc-len-f $trunc_len_f --p-trunc-len-r $trunc_len_r  --o-representative-sequences $rep_seqs --o-table $feature_table --o-denoising-stats $denoise_stats};

system qq{$QIIME2 && export TMPDIR=$qiime2_tmp && qiime tools export --input-path $denoise_stats.qza --output-path $dada2_path};
system qq{$QIIME2 && export TMPDIR=$qiime2_tmp && qiime tools export --input-path $rep_seqs.qza      --output-path $dada2_path};
system qq{$QIIME2 && export TMPDIR=$qiime2_tmp && qiime tools export --input-path $feature_table.qza --output-path $dada2_path};
system qq{$QIIME2 && export TMPDIR=$qiime2_tmp && biom convert -i $dada2_path/feature-table.biom -o $dada2_path/feature-table.txt --to-tsv};

change_feature_name("$dada2_path/dna-sequences.fasta", "$dada2_path/feature-table.txt", $dada2_path);

##6.final  可以直接拷贝给客户
my $final_path = "$output/6.final";
system qq{mkdir -p $final_path} if not -d $final_path;

my @commands;
push @commands, ['cp -rf', "$dada2_path/otu.fasta", "$final_path"];
push @commands, ['cp -rf', "$dada2_path/otu_table.txt", "$final_path"];
push @commands, ['grep -v numeric', "$dada2_path/stats.tsv", "> $final_path/seqStat.txt"];

check_copy_and_run(\@commands);


print "[OK] 运行正常\n";

#######################################################################子程序
sub check_empty_compressed_file{
    my @files = @_;
    my $err = "";
    foreach my $file(@files){
        my $flag = `less $file | head -n 1 | wc -l`;
        if($flag == 0){
            $err .= "[ERROR] $file 文件为空文件\n";
        }
    }
    return $err;
}

# -e 文件或目录存在 -s 文件或目录存在且不为0(返回字节数) # 空文件夹字节数4096 -f 文件为普通文件
sub check_copy_and_run{
    my $command_arr = shift;
    foreach my $command(@$command_arr){
            system qq{$$command[0] $$command[1] $$command[2]};
    }

}


sub combine_raw_data {
    my ($forward_trim, $reverse_trim, $combine_trim, $raw_input, $thread) = @_;
    my %forward_hash = read_manifest($forward_trim);
    my %reverse_hash = read_manifest($reverse_trim);
    
    open COM, ">$combine_trim/NO-AP-MANIFEST";
    print COM "sample-id,absolute-filepath,direction\n";
    open XLS, ">$combine_trim/trim_succeed.xls";
    print XLS "sample\traw_reads\tsuccess_trim_reads\tsuccess_trim_perc\n";
    
    my $ForkManager = Parallel::ForkManager->new($thread);
    foreach my $sample (keys %forward_hash) {
        my $pid = $ForkManager->start and next;

        # 合并正向切除引物和反向切除引物fastq, 同时记录匹配成功的序列，从原始序列中筛选输出引物匹配失败的序列
        my $failed_output = "$combine_trim/failed"; system qq{mkdir -p $failed_output} if not -d $failed_output;
        my($succeed_reads, $raw_reads) = read_and_combine($forward_hash{$sample}{"forward"}, $reverse_hash{$sample}{"reverse"}, "$combine_trim/$sample\_R1.fastq.gz", "$raw_input/$sample\_R1.fastq.gz", "$failed_output/$sample\_R1.fastq.gz");
        read_and_combine($forward_hash{$sample}{"reverse"}, $reverse_hash{$sample}{"forward"}, "$combine_trim/$sample\_R2.fastq.gz", "$raw_input/$sample\_R2.fastq.gz", "$failed_output/$sample\_R2.fastq.gz");

        # 计算引物匹配成功的序列比例
        my $succeed_perc  = 100*$succeed_reads/$raw_reads;
        print XLS "$sample\t$raw_reads\t$succeed_reads\t$succeed_perc\n";
        print COM "$sample,$combine_trim/$sample\_R1.fastq.gz,forward\n";
        print COM "$sample,$combine_trim/$sample\_R2.fastq.gz,reverse\n";
        
        $ForkManager->finish;
    }
    $ForkManager->wait_all_children;
    close COM;
    close XLS;
}

sub read_manifest{
    my $export_dir = shift;
    my %hash;
    open MANI,"$export_dir/MANIFEST" or die "[ERROR]: can not open $export_dir/MANIFEST !!";
    while (<MANI>){
        next if /sample-id,filename,direction/;
        $_ =~ s/[\r\n]//g;
        my @tmp = split /,/,$_;
        my $abspath = join "/", "$export_dir", $tmp[1];
        $hash{$tmp[0]}{$tmp[2]} = $abspath;
    }
    close MANI;
    return %hash;
}

sub read_and_combine{
    my($gz_1, $gz_2, $combine_gz, $raw_gz, $fail_gz) = @_;
    
    # 读入并合并trim_forward trim_reverse
    open GZ, "zcat $gz_1 $gz_2 |";
    open TRIM_COM, "|gzip > $combine_gz";
    
    my %succeed;
    while(my $id = <GZ>){
        my $seq = <GZ>;
        my $plus = <GZ>;
        my $score = <GZ>;
        map{$_ =~ s/[\r\n]//g}($id, $seq, $plus, $score);
        my $real_id = (split/\s+/, $id)[0];
        $real_id =~ s/\/1$//;
        $real_id =~ s/\/2$//;
        next if (exists $succeed{$real_id});
        $succeed{$real_id} = 1;
        print TRIM_COM "$id\n$seq\n$plus\n$score\n";
    }
    close GZ;
    close TRIM_COM;
    
    # 筛选输出未检测出引物的序列
    open RAW, "zcat $raw_gz |";
    open FAIL, "|gzip > $fail_gz";
    my $raw_reads;
    while(my $id = <RAW>){
        my $seq = <RAW>;
        my $plus = <RAW>;
        my $score = <RAW>;
        $raw_reads++;
        map{$_ =~ s/[\r\n]//g}($id, $seq, $plus, $score);
        my $real_id = (split/\s+/, $id)[0];
        $real_id =~ s/\/1$//;
        $real_id =~ s/\/2$//;
        next if exists $succeed{$real_id};
        print FAIL "$id\n$seq\n$plus\n$score\n";
    }
    close RAW;
    close FAIL;

    my $err = check_empty_compressed_file($combine_gz);
    die $err."检查下机数据是否异常以及引物是否写错\n" if $err ne "";

    return((scalar keys %succeed), $raw_reads);
}

sub prepare_rawdata_manifest {
    ###生成manifest配置文件
    my ($input, $manifest, $samples) = @_;
    open PATH, ">$manifest";
    print PATH "sample-id,absolute-filepath,direction\n";

    foreach my $sample(@$samples){
        $sample =~ s/\s//g;
        $sample =~ s/[\r\n]//g;
        print PATH "$sample,$input/$sample\_R1.fastq.gz,forward\n";
        print PATH "$sample,$input/$sample\_R2.fastq.gz,reverse\n";
    }
    close PATH;
}

sub change_feature_name {
    my ($seq_fasta, $feature_table, $output) = @_;
    my $prefix = "ASV";
    my $cnt = 1;

    my %feature_name_hash;
    open TABLE, $feature_table or die "无法打开文件 $feature_table\n";
    open NEW_T, ">$output/otu_table.txt" or die "无法写入文件 $output/otu_table.txt\n";
    open MAP,   ">$output/feature_HASHid_vs_ASVid.txt" or die "无法写入文件 $output/feature_HASHid_vs_ASVid.txt\n";
    while (<TABLE>) {
        $_ =~ s/[\r\n]//g;
        next if /^# Constructed from biom file/;  
        if(/^#OTU ID/){
            print NEW_T "$_\n" ;
        }else{
            my @data_t = split /\t/, $_;
            my $new_name = $prefix.$cnt;
            $feature_name_hash{$data_t[0]} = $new_name;
            my $new_data = join "\t", $new_name, map{sprintf("%.0f", $data_t[$_])}1..$#data_t;
            print NEW_T "$new_data\n";
            print MAP "$data_t[0]\t$new_name\n";
            $cnt++;            
        }
    }
    close TABLE;
    close NEW_T;
    close MAP;


    open NEW_F, ">$output/otu.fasta" or die "无法写入 $output/otu.fasta\n";
    my $FASTA_SEQ = Bio::SeqIO->new(-file => $seq_fasta, -format=>'Fasta') or die "不存在文件 $seq_fasta\n";
    while(my $inSeq = $FASTA_SEQ->next_seq){
        my $id  = $inSeq->display_id;
        my $seq = $inSeq->seq;  
        print NEW_F ">$feature_name_hash{$id}\n$seq\n" if(exists $feature_name_hash{$id});
    }
    close NEW_F;

}


##读取config并检查必须参数
sub read_config{
    my $file           = shift;
    my $parameter_list = shift;
    my %hash;
    open FILE, $file or die "[ERROR] 无法读取配置文件, $file\n";
    while(my $line = <FILE>){
        $line =~ s/#.+//g;
        $line =~ s/^\s+//;
        $line =~ s/\s+$//;
        my @dat = split /\s*=\s*/, $line, 2;
        my ($a, $b) = @dat;
        next if(not defined $b);
        next if($a =~ /ENV/);
        $hash{$a}{"value"} = $b;
    }
    close FILE;
    
    ##check_parameter
    my $err = "";
    foreach my $parameter(@$parameter_list){
        $err .= "[ERROR] 没有定义 $parameter\n" if(not exists $hash{$parameter});
    }
    die $err if($err ne "");
    return %hash;
}

##读取table
sub read_table{
    my $file = shift;
    my %hash;
    if($file eq ""){
        return %hash;
    }else{
        open IN,$file or die "[ERROR] Cannot open $file\n";
        my $title = <IN>;
        $title =~ s/[\r\n]//g;
        my @head = split/\t/,$title;
        while(<IN>){
            $_ =~ s/[\r\n]//g;
            my @tmp = split/\t/,$_;
            foreach my $i(1..$#head){
                my $value = exists $tmp[$i] ? $tmp[$i] : "";
                $hash{$tmp[0]}{$head[$i]} = $tmp[$i];
            }
        }
        close IN;
        return %hash;
    }
}

sub seq_distribution{
    ###统计 feature reads
    my $feature_path = shift;
    my %hash_seq_length_sum;
    open FEAT, "$feature_path/feature-table.txt", or die "cannot open $feature_path/feature-table.txt, please check it !\n";
    while (<FEAT>) {
        next if /^#\sConstructed\sfrom\sbiom\sfile/;
        chomp $_;
        next if (/^#OTU/) ;
        my @split_feat = split /\t/, $_;
        my $featuresum = 0;
        for my $j (1..$#split_feat){
            $featuresum += $split_feat[$j];
        }    
        $hash_seq_length_sum{$split_feat[0]}{"sum"} = $featuresum;
    }
    close FEAT;

    ###统计per feature seq length
    my $max = "NA";
    my $min = "NA";
    open SEQ, "$feature_path/dna-sequences.fasta" or die "cannot open $feature_path/dna-sequences.fasta, please check it !\n";
    while (my $id = <SEQ>) {
        my $seq = <SEQ>;
        chomp $id;
        chomp $seq;
        $id =~ s/^>//g;
        $hash_seq_length_sum{$id}{"length"} = length($seq);
        $max = length($seq) if $max eq "NA";
        $min = length($seq) if $min eq "NA";
        $max = ($max <= length($seq))? length($seq) : $max;
        $min = ($min >= length($seq))? length($seq) : $min;
    }
    close SEQ;

    ###划分区段
    my @range_arr;
    my %range_hash;
    my $cnt = 0;
    if (($max-$min) > 80) {
        my $step = sprintf "%.0f", (($max-$min)/40)+0.6;
        for my $i (0..39){
            push @range_arr, ($min+($step*$i))."-".($min+($step*($i+1))-1);
            my $key = ($min+($step*$i))."-".($min+($step*($i+1))-1);
            $range_hash{$cnt}{$key} = 0;
            $cnt++;
        }
    }else{
        my $step = sprintf "%.0f", (($max-$min)/20)+0.6;
        for my $i (0..19){
            push @range_arr, ($min+($step*$i))."-".($min+($step*($i+1))-1);
            my $key = ($min+($step*$i))."-".($min+($step*($i+1))-1);
            $range_hash{$cnt}{$key} = 0;
            $cnt++;
        }    
    }

    ####区段统计
    foreach my $feat_id (keys %hash_seq_length_sum) {
        for my $i (0..$#range_arr){
            my ($s, $e) = split "-", $range_arr[$i];
            if( $hash_seq_length_sum{$feat_id}{"length"} >= $s and $hash_seq_length_sum{$feat_id}{"length"} <= $e){
                $range_hash{$i}{$range_arr[$i]} += $hash_seq_length_sum{$feat_id}{"sum"};            
            }
        }
    }

    ###输出统计结果
    open OUT, qq{>$feature_path/reads.length.distribution.xls} or die "Can not write into $feature_path, please check it!\n";
    for my $i (0..$#range_arr){
        foreach my $ragne_id (keys %{$range_hash{$i}}){
            print OUT $ragne_id."\t".$range_hash{$i}{$range_arr[$i]}."\n";
        }
    }
    close OUT;
}

sub get_final_trunc_len{
    my($ref_file, $final_file, $trunc_len_f, $trunc_len_r, $trim_left_f, $trim_left_r) = @_;
    
    my %hash;
    open IN, $ref_file or die "[ERROR] Cannot open $ref_file\n";
    while(<IN>){
        $_ =~ s/[\r\n]//g;
        my($key, $value) = split/\t/, $_;
        $hash{$key} = $value;
    }
    close IN;
    
    my $final_trunc_len_f = $trunc_len_f eq "" ? $hash{"--p-trunc-len-f"} : $trunc_len_f;
    my $final_trunc_len_r = $trunc_len_r eq "" ? $hash{"--p-trunc-len-r"} : $trunc_len_r;
    
    open OUT, ">$final_file" or die "[ERROR] Cannot write into $final_file\n";
    print OUT "--p-trim-left-f $trim_left_f\n--p-trim-left-r $trim_left_r\n--p-trunc-len-f $final_trunc_len_f\n--p-trunc-len-r $final_trunc_len_r\n";
    close OUT;
    
    return($final_trunc_len_f, $final_trunc_len_r);
}