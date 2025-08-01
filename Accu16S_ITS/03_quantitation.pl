$| = 1;
# Perl 系统包
use warnings;
use strict;
use File::Spec;
use Getopt::Long;
use List::Util qw/sum/;
use Sort::Key::Natural qw(natsort);
use FindBin qw($RealBin $RealScript);

# 常量
use constant SCRIPT_DIR => (File::Spec->splitpath(File::Spec->rel2abs($0)))[1];
use constant PWD => $ENV{"PWD"};


my $DEFAULT_THREAD = 40;
my $DEFAULT_MEMORY = 50;
my $DEFAULT_ALLOW_SPIKEIN_ZERO   = "FALSE";
my $DEFAULT_SPIKEIN_ZERO_WARNING = "Gs_BSI7,Gs_BSI8,Gs_BSI9,Gs_ITS-SPK5,Gs_ITS-SPK6,Gs_ITS-SPK7";
my $DEFAULT_SPIKEIN_DB_16S       = "/home/genesky/database_new/self_build_database/MGS/SpikeIN/16S_SPIKEIN.fasta";
my $DEFAULT_SPIKEIN_DB_ITS       = "/home/genesky/database_new/self_build_database/MGS/SpikeIN_ITS/new/ITS_SPIKEIN.fasta";
my $DEFAULT_BLASTN               = "/home/genesky/software/blast+/2.9.0/bin/blastn";
my $DEFAULT_SPIKEIN_INFO         = "/home/genesky/database_new/self_build_database/MGS/spikein_set.txt";
my $DEFAULT_SOFT_R_SCRIPT        = "/home/genesky/software/r/3.5.1/bin/Rscript";
my $DEFAULT_SOFT_R_LIB           = "/home/genesky/database/r/3.5.1";
my $DEFAULT_SOFT_R_SCRIPT4       = "/home/genesky/software/r/4.0.3/bin/Rscript";
my $DEFAULT_SOFT_R_LIB4          = "/home/genesky/software/r/4.0.3/lib64/R/library";
my $CALCULATE_ABSOLUTE_COPIES    = SCRIPT_DIR."/quantitation_calculate_absolute_copies.r";

# 参数输入
my $ARGV_INFO = join " ", @ARGV;
my ($config, $input, $copy_prediction, $output, $SPIKEIN_DB, $thread, $no_check, 
    $SPIKEIN_INFO, $BLASTN, $SOFT_R_SCRIPT, $PRIMERDB, $SOFT_R_LIB, $SOFT_R_SCRIPT4, $SOFT_R_LIB4, $if_help);
GetOptions(

    "config|c=s"           => \$config,
    "input|i=s"            => \$input,
    "copy-prediction=s"    => \$copy_prediction,
    "spikein-db=s"         => \$SPIKEIN_DB,
    "output|o=s"           => \$output,
    "no-check"             => \$no_check,
    "thread=i"             => \$thread,
    "spikein-info|si=s"    => \$SPIKEIN_INFO,
    "blastn=s"             => \$BLASTN,
    "Rscript=s"            => \$SOFT_R_SCRIPT,
    "R_Lib=s"              => \$SOFT_R_LIB,
    "Rscript4=s"           => \$SOFT_R_SCRIPT4,
    "R_Lib4=s"             => \$SOFT_R_LIB4,
    "help|h"               => \$if_help,
);
die "
Options: 必填

        --config/-c                参数配置文件，用到的参数有：
                                   F_Primer  : 正向测序引物, 用于提取spikein-info信息, 必须填写
                                   R_Primer  : 反向测序引物, 用于提取spikein-info信息, 必须填写
                                   Amplicon_Area : 扩增区域, 可选 V4 / V3V4 / V4V5 / ITS1 / ITS2, 用于提取spikein-info信息, 必须填写
                                   Amplicon_Type : 扩增子类型，可选 16S / 18S / ITS1 / ITS2, 用于提取spikein-info信息
                                   IR_Copies    :IR_copies文件,里面包含样本spikein的添加量信息,实验提供, 必须填写
                                   DNA_Quantity :DNA_quantity文件,里面包含各样本DNA上机量, DNA抽提量和抽提DNA所用样本量, 实验提供, 必须填写
                                   Allow_Spikein_Zero : 是否允许样本中某条内参检出值为0, TRUE/FALSE, (default: '$DEFAULT_ALLOW_SPIKEIN_ZERO')
                                   Spikein_Zero_Warning : 当这些内参在样本中检出值为0时，只给出Warning警告信息，不报错中断。(default: '$DEFAULT_ALLOW_SPIKEIN_ZERO')
                                   Drop_Spikein: 在绝对定量计算过程中，线性拟合不纳入指定的spike
                                   Gene_Copy_Corr_Method: 基因拷贝数矫正方法，可选：RDP/PICRUST2/NONE, NONE表示不进行拷贝数矫正，必须填写
        --input/-i                 输入数据路径, 该目录下必须同时包含: *_otu.txt, *_rep.fasta, *_spike_otu.txt, *_spike_rep.fasata
        --spikein-db               引物信息数据库,必须填写 ( ITS:$DEFAULT_SPIKEIN_DB_ITS   16S : $DEFAULT_SPIKEIN_DB_16S)
        --spikein-info/-si         spikein的位置信息, 必须填写( eg.$DEFAULT_SPIKEIN_INFO )
        --output/-o                数据结果路径, intm_result/quantitation

Options: 可选

        --copy-prediction          拷贝数预测数据路径，当 Gene_Copy_Corr_Method 不为 NONE 时必须定义，该目录下必须包好：*_copynumber.txt
        --thread                   使用线程数 (default: $DEFAULT_THREAD)
        --no-check                 不计算空闲CPU数量
        --blastn                   blastn软件，(default: $DEFAULT_BLASTN)
        --Rscript                  更改软件 Rscript 版本 (default: '$DEFAULT_SOFT_R_SCRIPT')
        --R_Lib                    更改软件 R Lib路径 (default: $DEFAULT_SOFT_R_LIB)
        --Rscript4                 更改软件 Rscript 版本 (default: '$DEFAULT_SOFT_R_SCRIPT4')
        --R_Lib4                   更改软件 R Lib路径 (default: $DEFAULT_SOFT_R_LIB4)
        --help/-h                  查看帮助文档
\n" if (defined $if_help or not defined $config or not defined $input or not defined $SPIKEIN_DB  or not defined $SPIKEIN_INFO or not defined $output);
$thread        = $DEFAULT_THREAD if (not defined $thread);
$SOFT_R_SCRIPT = $DEFAULT_SOFT_R_SCRIPT   if (not defined $SOFT_R_SCRIPT);
$SOFT_R_LIB    = $DEFAULT_SOFT_R_LIB      if (not defined $SOFT_R_LIB);
$SOFT_R_SCRIPT4 = $DEFAULT_SOFT_R_SCRIPT4 if (not defined $SOFT_R_SCRIPT4);
$SOFT_R_LIB4    = $DEFAULT_SOFT_R_LIB4    if (not defined $SOFT_R_LIB4);
$SPIKEIN_INFO   = $DEFAULT_SPIKEIN_INFO    if (not defined $SPIKEIN_INFO);
$BLASTN         = $DEFAULT_BLASTN          if (not defined $BLASTN);

my @parameter          = ("Amplicon_Type", "F_Primer", "R_Primer", "IR_Copies", "DNA_Quantity", "Gene_Copy_Corr_Method");
my %config             = read_config_with_check($config, \@parameter);
my $amplicon_area      = $config{"Amplicon_Area"}{"value"};
my $amplicon_type      = $config{"Amplicon_Type"}{"value"};
my $f_primer           = $config{"F_Primer"}{"value"};
my $r_primer           = $config{"R_Primer"}{"value"};
my $copy_corr_method   = $config{"Gene_Copy_Corr_Method"}{"value"};
my $ir_copies            = exists $config{"IR_Copies"}    ? $config{"IR_Copies"}{"value"}    : "";
my $dna_quantity         = exists $config{"DNA_Quantity"} ? $config{"DNA_Quantity"}{"value"} : "";
my $allow_spikein_zero   = exists $config{"Allow_Spikein_Zero"} ? $config{"Allow_Spikein_Zero"}{"value"} : $DEFAULT_ALLOW_SPIKEIN_ZERO;
my $spikein_zero_warning = exists $config{"Spikein_Zero_Warning"} ? $config{"Spikein_Zero_Warning"}{"value"} : $DEFAULT_SPIKEIN_ZERO_WARNING;
my $drop_spikein         = exists $config{"Drop_Spikein"} ? $config{"Drop_Spikein"}{"value"} : "null";

die "[ERROR] Gene_Copy_Corr_Method = $copy_corr_method 时必须定义拷贝数预测数据路径 --copy-prediction\n" if($copy_corr_method ne "NONE" and not defined $copy_prediction);
###################################################################### 初始化


# 运行参数

my $RUN_INFO = "
---
Command: perl $RealBin/$RealScript $ARGV_INFO
---
[SET] 软件 BLASTN        : $BLASTN
[SET] 软件 Rscript       : $SOFT_R_SCRIPT
[SET] 软件 R 数据库      : $SOFT_R_LIB
[SET] 软件 Rscript       : $SOFT_R_SCRIPT4
[SET] 软件 R 数据库      : $SOFT_R_LIB4
[SET] 参数配置文件       : $config
[SET] 输入数据路径       : $input
[SET] 正向测序引物       : $f_primer
[SET] 反向测序引物       : $r_primer
[SET] 扩增区域           : $amplicon_area
[SET] Gene_Copy_Corr_Method : $copy_corr_method
[SET] IR_copies文件      : $ir_copies
[SET] DNA_Quantity文件   : $dna_quantity
[SET] spikein fasta      : $SPIKEIN_DB
[SET] spikein的位置信息  : $SPIKEIN_INFO
[SET] 结果输出目录       : $output
";

system "mkdir -p $output" if not -d $output;
open SAVE, ">>$output/quatitation_run.info"; print SAVE $RUN_INFO; close SAVE;
print $RUN_INFO;

###################################################################### 主程序
my $final = "$output/final";system "mkdir -p $final" if not -d $final;
my %spikein_info    = select_spikein_info($SPIKEIN_INFO, $amplicon_area, $amplicon_type, $f_primer, $r_primer);

foreach my $prefix('client'){
    my $otu_txt         = "$input/$prefix\_otu.txt";
    my $rep_fasta       = "$input/$prefix\_rep.fasta";
    my $spike_fasta     = "$input/$prefix\_spike_rep.fasta";
    my $spike_otu_txt   = "$input/$prefix\_spike_otu.txt";
    next if not -e $spike_otu_txt;
    die "[ERROR] 不存在文件 $spike_fasta\n"     if not "$spike_fasta";
    
    # 0. 根据样本名提取 ir_copies 及 dna_quantity (client, y_sample)
    my $spike_IR_copies = "$output/$prefix\.0.spike_IR_copies.txt";
    my $DNA_quantity    = "$output/$prefix\.0.DNA_quantity.txt";
    extract_ir_copiese($otu_txt, $ir_copies, $spike_IR_copies);
    extract_dna_quantity($otu_txt, $dna_quantity, $DNA_quantity);

    # 1.blastn 两两比对，找到数据库中与输入序列最相似的序列
    my $blastn_result          = "$output/$prefix\.1.spike.blast";
    system qq{$BLASTN -query $spike_fasta -db $SPIKEIN_DB -evalue 1e-50 -outfmt 6 -max_target_seqs 100000000 -perc_identity 90 -num_threads 8 > $blastn_result};
    
    # 2.整理统计blastn结果
    my $pair_id                = "$output/$prefix\.2.pair_id.txt";  # OTU_ID 和 GS_ID 的对应关系
    my $match_spike_otu        = "$output/$prefix\.2.spike_otu.txt";  # 合并匹配到同一内参的OTU序列
    my $mismatch_spike_otu     = "$output/$prefix\.2.mismatch_spike_otu.txt";  # 未匹配到内参的OTU序列
    my $sample_sum             = "$output/$prefix\.2.sample_sum.txt";  # 统计各样本所有内参序列size数
    my $ir_detection_rscript   = "$output/$prefix\.2.ir_detection.r";
    ir_detection($blastn_result, $spike_otu_txt, \%spikein_info, $pair_id, $match_spike_otu, $mismatch_spike_otu, $sample_sum, $ir_detection_rscript);
    check_spike($match_spike_otu, $SPIKEIN_INFO, \%spikein_info, $allow_spikein_zero, $spikein_zero_warning);  # 检测各样本内参检出值是否非0
    
    # 3.根据内参添加量，计算拟合曲线  从而计算绝对拷贝数  DNA水平绝对拷贝数  SAMPLE 水平绝对拷贝数
    my $spike_reads_persent    = "$output/$prefix\.3.spike_reads_persent.txt";  # 成功匹配内参otu序列数和所有otu序列数统计
    my $standard_curve_formula = "$output/$prefix\.3.standard_curve_formula.txt";  # 拟合曲线方程
    my $absolute_otu_txt       = "$output/$prefix\.3.absolute_otu.txt";  # 计算绝对拷贝数后otu
    my $unit_dna_otu_copies    = "$output/$prefix\.3.unit_dna_otu_copies.txt";  # DNA水平绝对拷贝数
    my $unit_sample_otu_copies = "$output/$prefix\.3.unit_sample_otu_copies.txt";  # SAMPLE水平绝对拷贝数
    system qq{$SOFT_R_SCRIPT4 $CALCULATE_ABSOLUTE_COPIES -o $otu_txt -s $match_spike_otu -i $spike_IR_copies -d $DNA_quantity -p $prefix.3 --output $output --drop_spikein $drop_spikein --lib $SOFT_R_LIB4 &> $output/$prefix\_quatitation_run.log};
    formula_check($standard_curve_formula);
    
    # 4.基因拷贝数矫正 RDP/PICRUST2/NONE
    my $corrected_unit_dna_otu_copies    = "$output/$prefix\.4.corrected_unit_dna_otu_copies.txt";  # DNA水平基因拷贝数矫正后拷贝数
    my $corrected_unit_sample_otu_copies = "$output/$prefix\.4.corrected_unit_sample_otu_copies.txt"; # SAMPLE水平基因拷贝数矫正后拷贝数
    my $gene_copy_correct_rscript        = "$output/$prefix\.4.gene_copy_correct.r";

    my @commands;
    if (uc $copy_corr_method ne "NONE") {
        my $estimated_copynumber = "$copy_prediction/$prefix\_copynumber.txt";
        gene_copy_correct($unit_dna_otu_copies, $estimated_copynumber, $corrected_unit_dna_otu_copies, $gene_copy_correct_rscript);
        gene_copy_correct($unit_sample_otu_copies, $estimated_copynumber, $corrected_unit_sample_otu_copies, $gene_copy_correct_rscript) if -e $unit_sample_otu_copies;
        push @commands, ['cp -rf', "$estimated_copynumber",             "$final/$prefix\_copynumber.txt"];
        push @commands, ['cp -rf', "$corrected_unit_dna_otu_copies",    "$final/$prefix\_corrected_unit_dna_otu_copies.txt"];
        push @commands, ['cp -rf', "$corrected_unit_sample_otu_copies", "$final/$prefix\_corrected_unit_sample_otu_copies.txt"] if -e $corrected_unit_sample_otu_copies;
    }
    
    ##整理结果
    push @commands, ['cp -rf', "$otu_txt",                "$final"];
    push @commands, ['cp -rf', "$rep_fasta",              "$final"];
    push @commands, ['cp -rf', "$match_spike_otu",        "$final/$prefix\_spike_otu.txt"];
    push @commands, ['cp -rf', "$absolute_otu_txt",       "$final/$prefix\_absolute_otu.txt"];
    push @commands, ['cp -rf', "$standard_curve_formula", "$final/$prefix\_standard_curve_formula.txt"];
    push @commands, ['cp -rf', "$spike_reads_persent",    "$final/$prefix\_spike_reads_persent.txt"];
    push @commands, ['cp -rf', "$unit_dna_otu_copies",    "$final/$prefix\_unit_dna_otu_copies.txt"];
    push @commands, ['cp -rf', "$unit_sample_otu_copies", "$final/$prefix\_unit_sample_otu_copies.txt"] if -e $unit_sample_otu_copies;
    
    my $err = check_copy_and_run(\@commands);
    die $err if $err ne "";
}
print "\n[OK] 运行完成\n";

###################################################################### 子程序
sub check_copy_and_run{
    my $command_arr = shift;
    my $err = "";
    foreach my $command(@$command_arr){
        my $err_info = "";
        if(-f $$command[1]){ #文件存在，只检测为空文件的情况
            if($$command[1] =~ /.pdf$/){
                $err_info .= check_result_pdf($$command[1]);
            }else{
                $err_info .= check_result_txt($$command[1]);
            }
        }else{
            $err_info .= check_result_dir($$command[1]);
        }

        if($err_info eq ""){
            system qq{$$command[0] $$command[1] $$command[2]};
        }else{
            $err .= $err_info;
        }

    }
    return $err;
}

sub gene_copy_correct {
    my ($otu_file, $gene_copy_file, $new_otu_file, $gene_copy_correct_rscript) = @_;
    my $rscript = qq{
.libPaths('$SOFT_R_LIB')
options(scipen = 200)
otu      = read.table('$otu_file', header = TRUE, check.name = F, quote = "", comment.char = "", sep = "\\t", fill = T)
genecopy = read.table('$gene_copy_file', header = TRUE, check.name = F, quote = "", comment.char = "", sep = "\\t", fill = T)
rownames(otu)  = as.character(otu[, 1])
rownames(genecopy) = as.character(genecopy\$OTU_id)

lost_ids = rownames(otu)[! rownames(otu) %in% rownames(genecopy)]
if(length(lost_ids) > 0){
    message("[ERROR] 拷贝数信息缺失：", paste(lost_ids, collapse = ","))
    q()
}
genecopy = genecopy[rownames(otu), , drop = F]

result = list()
for (i in 1:nrow(otu)) {
    count = round(otu[i, 2:(ncol(otu)-9)]/genecopy[i, 2], 0)
    # result = rbind(result, count)
    result[[i]] = count
}

result = data.frame(matrix(unlist(result), ncol=length(result[[1]]), byrow=T), stringsAsFactors=FALSE)

result_data = data.frame(OTUid = otu[,1], result, OTUsize = rowSums(result), otu[,(ncol(otu)-7):ncol(otu)])
colnames(result_data) = c("#OTU ID", colnames(otu)[2:(ncol(otu)-9)], "OTUsize", colnames(otu)[(ncol(otu)-7):ncol(otu)])
write.table(result_data, file = '$new_otu_file', quote = FALSE, sep = "\\t", row.names = F)
    };
    open SAVE, ">$gene_copy_correct_rscript"; print SAVE $rscript; close SAVE;
    system "$SOFT_R_SCRIPT $gene_copy_correct_rscript";
}

sub formula_check {
    my $file = shift;
    my @Warns;
    open FOA, "$file" or die "[ERROR]Cannot open $file\n";
    my $head = <FOA>;
    while(<FOA>){
        $_ =~ s/[\r\n]//g;
        my @dat = split /\t/, $_;
        push @Warns, $dat[0] if ($dat[2] < 0.95);
    }
    close FOA;
    print "[Warnings] 样本的R2值小于0.95, @Warns!\n" if (scalar @Warns != 0)
}

sub ir_detection{
    my ($blastn_result, $spike_otu_txt, $spikein_info, $pair_id, $match_spike_otu, $mismatch_spike_otu, $sample_sum, $ir_detection_rscript) = @_;
    open BLASTN, $blastn_result or die "[ERROR] Cannot open $blastn_result\n";
    open ID,">$pair_id" or die "[ERROR]Cannot write into $pair_id\n";
    print ID "OTU_ID\tGs_BSI_ID\n";
    my $pass_spikein = 0;
    while (<BLASTN>) {
        $_ =~ s/[\r\n]//g;
        my @arr = split /\t/, $_;
        if($arr[4] ~~ @{$$spikein_info{$arr[1]}{"no_map_num"}} and $arr[5] ~~ @{$$spikein_info{$arr[1]}{"gap_num"}} and $arr[6] ~~ @{$$spikein_info{$arr[1]}{"seq_start_index"}} and $arr[7] ~~ @{$$spikein_info{$arr[1]}{"seq_end_index"}}){
            next if($amplicon_area =~ /_FL$/ and ($arr[2] < 96 or $arr[3] < 1300));
            print ID "$arr[0]\t$arr[1]\n";
            $pass_spikein++;
        }
    }
    close BLASTN;
    close ID;

    die "[ERROR] 没有筛选到符合条件的内参比对结果，检查扩增子区域及spikein-info信息是否正确\n" if $pass_spikein == 0;

    my $rscript = qq{
.libPaths('$SOFT_R_LIB')
pair_id_data       = read.table('$pair_id', header = TRUE,row.names = 1, check.name = F, quote = "", comment.char = "", sep = "\\t", fill = T)
spikein_data       = read.table('$spike_otu_txt', header = TRUE, row.names = 1, check.name = F, quote = "", comment.char = "", sep = "\\t", fill = T)
match_spikein_data = spikein_data[rownames(pair_id_data),,drop = F]

##没有匹配到内参序列的OTU
mismatch_spike    = spikein_data[setdiff(rownames(spikein_data),rownames(pair_id_data)),,drop = F]
mismatch_spike    = cbind(rownames(mismatch_spike), mismatch_spike)
colnames(mismatch_spike)[1] = "OTU_ID"
mismatch_spike    = mismatch_spike[order(mismatch_spike\$OTU_ID),,drop = F]
write.table(mismatch_spike, file = '$mismatch_spike_otu', quote = FALSE, sep = "\\t", row.names = F)


##合并匹配到同一内参序列的otu
match_spike_otu = cbind(pair_id_data,match_spikein_data)
match_spike_otu = aggregate(match_spike_otu[,2:ncol(match_spike_otu)],by = match_spike_otu["Gs_BSI_ID"],FUN=sum)
colnames(match_spike_otu) = c("Gs_BSI_ID",colnames(match_spikein_data))##只有一个样本时
match_spike_otu = match_spike_otu[order(match_spike_otu\$Gs_BSI_ID),,drop = F] ##排序
write.table(match_spike_otu, file = '$match_spike_otu', quote = FALSE, sep = "\\t", row.names = F)


##以样本为单位统计spike_otu
samplename    = colnames(match_spikein_data)
sum_sample    = rbind(sapply(samplename, function(t){sum(spikein_data[,t])}),
                      sapply(samplename, function(t){sum(match_spike_otu[,t])}),
                      sapply(samplename, function(t){sum(match_spike_otu[,t])*100/sum(spikein_data[,t])}))
sum_sample<-rbind(colnames(sum_sample),sum_sample)
rownames(sum_sample) = c("Sample","Totle_spikein_reads", "Match_spikein_reads", "Match_spikein_percent(%)")
write.table(t(sum_sample), file = '$sample_sum', sep = "\\t",row.names = F, col.names = T, quote = F)
    };
    open SAVE, ">$ir_detection_rscript"; print SAVE $rscript; close SAVE;
    system "$SOFT_R_SCRIPT $ir_detection_rscript";
}

sub select_spikein_info {
    my ($SPIKEIN_INFO, $amplicon_area, $amplicon_type, $f_primer, $r_primer) = @_;
    my $spikein_pos_name = $amplicon_type eq '18S' ? "$amplicon_type\_$amplicon_area" : $amplicon_area;
    my %spikein_info;
    open IN, $SPIKEIN_INFO or die "[ERROR] Cannot open $SPIKEIN_INFO\n";
    my $title = <IN>;
    $title =~ s/[\r\n]//g;
    my @split_title = split/\t/, $title;
    my %index;
    map{$index{$split_title[$_]} = $_}(0..$#split_title);
    while (<IN>){
        $_ =~ s/[\r\n]//g;
        my @data = split /\t/, $_;
        if ($data[$index{"region"}] eq $spikein_pos_name) {
            @{$spikein_info{$data[$index{"Gs_BSI_ID"}]}{"no_map_num"}}      = 0..$data[$index{"no_map_num"}];##blast 结果允许的错配数
            @{$spikein_info{$data[$index{"Gs_BSI_ID"}]}{"gap_num"}}         = 0..$data[$index{"gap_num"}]; ##blast 结果允许的gap数
            @{$spikein_info{$data[$index{"Gs_BSI_ID"}]}{"seq_start_index"}} = ($data[$index{"seq_F_start_index"}] - 0)..($data[$index{"seq_F_start_index"}] + 4); ##blast结果中， reads的起始位点 
            my $combine_length_no_p = $data[$index{"length"}] - length($f_primer) - length($r_primer);
            @{$spikein_info{$data[$index{"Gs_BSI_ID"}]}{"seq_end_index"}}  = ($combine_length_no_p - 3)..($combine_length_no_p + 3); ##blast结果中 reads的终止位点 
        }
    }
    close IN;
    return %spikein_info;
}

sub check_spike{
    my($match_spike_otu, $SPIKEIN_INFO, $spikein_info, $allow_spikein_zero, $spikein_zero_warning) = @_;

    my %zero_warn;
    map{$zero_warn{$_} = 1}(split/\s*,\s*/, $spikein_zero_warning);

    # 应检测出的内参 1.去掉允许未检测出的内参 2.去掉多余的内参序列 Gs_BSI10 或 Gs_ITS-SPK-Long
    my $ref_spike_num = scalar(keys %$spikein_info);
    my @ref_rm = ("Gs_BSI10", "Gs_ITS-SPK-Long", (keys %zero_warn));
    map{$ref_spike_num-- if($_ ~~ @ref_rm)}(keys %$spikein_info);
    my $real_spike_num = `cat $match_spike_otu|wc -l ` - 1;
    die "[ERROR] 应检测出 $ref_spike_num 条内参，实际检测出 $real_spike_num 条\n具体信息请查看文件 $match_spike_otu\n" if ( $real_spike_num < $ref_spike_num);
    ##检测各样本内参检出值是否非0
    open IN, $match_spike_otu or die"[ERROR] Cannot open $match_spike_otu\n";
    my $title = <IN>;
    $title =~ s/[\r\n]//g;
    my @split_title = split/\t/, $title;
    my @zero_spikein;
    while(<IN>){
        $_ =~ s/[\r\n]//g;
        my @data = split/\t/,$_;
        foreach my $i(1..$#split_title){
            push @zero_spikein, "$split_title[$i],$data[0]" if $data[$i] == 0;
        }
    }
    close IN;
    
    my $warn_info = "";
    my $err_info = "";
    if(scalar @zero_spikein > 0){
        foreach my $info(@zero_spikein){
            my($sample, $spikein_name) = split/,/, $info;
            if(exists $zero_warn{$spikein_name} or $allow_spikein_zero eq "TRUE"){
                $warn_info .= "样本 $sample 中内参 $spikein_name 的检出值为0\n";
            }else{
                $err_info .= "样本 $sample 中内参 $spikein_name 的检出值为0\n";
            }
        }
    }
    
    if($warn_info ne ""){
        print $warn_info."[Warnings] 具体信息请查看文件 $match_spike_otu\n";
    }
    
    if($err_info ne ""){
        die $err_info."[ERROR] 具体信息请查看文件 $match_spike_otu\n";
    }
}

sub extract_ir_copiese{
    my ($otu_file, $raw_ir_copies, $target_ir_copies) =@_;
    ##获取需要样本的下标
    open OTU ,$otu_file or die "[ERROR] Cannot open $otu_file\n";
    my $title = <OTU>;
    $title =~ s/[\r\n]//g;
    my @split_title = split/\t/, $title;
    my $index = ($title =~ /Taxonomy/)? @split_title - 10 : @split_title - 1;
    close OTU;
    # 读取IR_Copies.txt
    my %hash;
    open IR, $raw_ir_copies or die "[ERROR] Cannot open $raw_ir_copies\n";
    my $line1 = <IR>;
    $line1 =~ s/[\r\n]//g;
    my @head = split/\t/, $line1;
    my @GBS;
    while(my $line = <IR>){
        $line =~ s/[\r\n]//g;
        next if $line eq "";
        my @data = split/\t/, $line;
        map{$hash{$head[$_]}{$data[0]} = $data[$_]}(1..$#head);
        push @GBS, $data[0];  # 顺序输出
    }
    close IR;

    # 对目标样本进行提取
    open NEW_IR, ">$target_ir_copies" or die "[ERROR] Cannot write into $target_ir_copies\n";
    print NEW_IR "#OTU ID\t".join("\t", @split_title[1..$index])."\n";
    foreach my $GS(natsort @GBS){
        print NEW_IR $GS;
        foreach my $sample(@split_title[1..$index]){
            # 判断是不是 GS_YD，如果是，则用其配对客户样本数据填充
            if($sample =~/^GS_YD/){
                my $sample_client = $sample;
                $sample_client =~s/^GS_YD_//;
                $hash{$sample}{$GS} = $hash{$sample_client}{$GS};
            }
            die "[ERROR] IR_copies.txt 中缺失样本 $sample 的 $GS 值\n" if not exists $hash{$sample}{$GS};
            print NEW_IR "\t".$hash{$sample}{$GS};
        }
        print NEW_IR "\n";
    }
    close NEW_IR;
}

sub extract_dna_quantity{
    my ($otu_file, $raw_dna_quantity, $target_dna_quantity) = @_;
    ##获取需要样本的下标
    open OTU ,$otu_file or die "[ERROR] Cannot open $otu_file\n";
    my $title = <OTU>;
    $title =~ s/[\r\n]//g;
    my @split_title = split/\t/, $title;
    my $index = ($title =~ /Taxonomy/)? @split_title - 10 : @split_title - 1;
    my @samples = @split_title[1..$index];
    close OTU;
    ##读取DNA_quantity.txt并输出需要的样本
    open IN, $raw_dna_quantity or die "[ERROR] Cannot open $raw_dna_quantity\n";
    open OUT, ">$target_dna_quantity" or die "[ERROR] Cannot write into $target_dna_quantity\n";
    my $line1 = <IN>;
    print OUT $line1;
    while(<IN>){
        $_ =~ s/[\r\n]//g;
        my @data = split/\t/, $_;
        print OUT join("\t", @data)."\n" if ($data[0] ~~ @samples);
        # 判定样本是否有配对的 GS_YD ，有的话，再输出一份 GS_YD 的版本
        my $sample_gs_yd = "GS_YD_$data[0]";
        if($sample_gs_yd ~~ @samples){
            $data[0] = $sample_gs_yd;
            print OUT join("\t", @data)."\n";
        }
    }
    close IN;
    close OUT;
}

sub read_config_with_check{
    my $file = shift;
    my $parameter_list = shift;
    my $err = "";
    my $env_id;
    my %hash;
    open FILE, $file or die "[ERROR] 无法读取配置文件, $file\n";
    while(my $line = <FILE>){
        $line =~ s/#.+//g;
        $line =~ s/^\s+//;
        $line =~ s/\s+$//;
        my @dat = split /\s*=\s*/, $line, 2;
        my ($a, $b) = @dat;
        next if(not defined $b or $b eq "");

        if($a eq "Sample_Group"){

            my($group_id, $value) = split /\s*=\s*/, $b, 2;  # 多个分组文件，不同的rename
            next if (not defined $value or $value eq "");

            # Sample_Group = G2 = /path/to/group.txt  G2只能出现一次
            if(exists $hash{$a}{$group_id}){
                $err .= "[ERROR] Sample_Group 对应多个重复的 \'$group_id\'\n";
                next;
            }else{
                $hash{$a}{$group_id} = $value;
            }

        }elsif($a eq "ENV"){

            my($group_id, $value) = split /\s*=\s*/, $b, 2;  # 多个分组文件，不同的rename
            next if (not defined $value or $value eq "");
            if(exists $hash{$a}{$group_id}){
                $env_id++;
                $hash{$a}{$group_id}{$env_id} = $value;
            }else{
                $env_id = 1;
                $hash{$a}{$group_id}{$env_id} = $value;
            }
        }elsif($a eq "Quality_levels"){
            $hash{$a}{$b} = 1;
        }else{
            $hash{$a}{"value"} = $b;
        }
    }
    close FILE;

    # check_parameter
    @$parameter_list = () if(exists $hash{"Platform"}{"value"} and $hash{"Platform"}{"value"} eq "gsonline");
    foreach my $parameter(@$parameter_list){
         $err .= "[ERROR] 没有定义 $parameter\n" if(not exists $hash{$parameter});
    }

    die $err if($err ne "");
    return %hash;
}
