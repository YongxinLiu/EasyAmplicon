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


my $DEFAULT_THREAD = 5;
my $DEFAULT_SPIKEIN_DB_16S       = "/home/genesky/database_new/self_build_database/MGS/SpikeIN/16S_SPIKEIN.fasta";
my $DEFAULT_BLASTN               = "/home/genesky/software/blast+/2.9.0/bin/blastn";
my $DEFAULT_SPIKEIN_INFO         = "/home/genesky/database_new/self_build_database/MGS/spikein_set.txt";
my $DEFAULT_SOFT_R_SCRIPT        = "/home/genesky/software/r/4.0.3/bin/Rscript";
my $CALCULATE_ABSOLUTE_COPIES    = SCRIPT_DIR."/quantitation_calculate_absolute_copies.r";

# 参数输入
my $ARGV_INFO = join " ", @ARGV;
my ($config, $input_client, $input_spike, $copy_prediction, $output, $SPIKEIN_DB, $thread, 
    $SPIKEIN_INFO, $BLASTN, $SOFT_R_SCRIPT, $PRIMERDB, $if_help);
GetOptions(

    "config|c=s"           => \$config,
    "input_client=s"            => \$input_client,
    "input_spike=s"            => \$input_spike,
    "copy-prediction=s"    => \$copy_prediction,
    "spikein-db=s"         => \$SPIKEIN_DB,
    "output|o=s"           => \$output,
    "thread=i"             => \$thread,
    "spikein-info|si=s"    => \$SPIKEIN_INFO,
    "blastn=s"             => \$BLASTN,
    "Rscript=s"            => \$SOFT_R_SCRIPT,
    "help|h"               => \$if_help,
);
die "
Options: 必填

        --config/-c                参数配置文件，用到的参数有：
                                   F_Primer  : 正向测序引物, 用于提取spikein-info信息, 必须填写
                                   R_Primer  : 反向测序引物, 用于提取spikein-info信息, 必须填写
                                   Amplicon_Area : 扩增区域, 可选 V4 / V3V4 / V4V5 / ITS1 / ITS2, 用于提取spikein-info信息, 必须填写
                                   IR_Copies    :IR_copies文件,里面包含样本spikein的添加量信息,实验提供, 必须填写
                                   DNA_Quantity :DNA_quantity文件,里面包含各样本DNA上机量, DNA抽提量和抽提DNA所用样本量, 实验提供, 必须填写
        --input_client             菌群ASV输入数据路径, 该目录下必须同时包含: otu.fasta otu_table.txt
        --input_spike              spike ASV输入数据路径, 该目录下必须同时包含: otu.fasta otu_table.txt
        --spikein-db               spike数据库 (  16S : $DEFAULT_SPIKEIN_DB_16S)
        --spikein-info/-si         spikein的位置信息, 必须填写( eg.$DEFAULT_SPIKEIN_INFO )
        --output/-o                数据结果路径, intm_result/quantitation

Options: 可选

        --thread                   使用线程数 (default: $DEFAULT_THREAD)
        --blastn                   blastn软件，(default: $DEFAULT_BLASTN)
        --Rscript                  更改软件 Rscript 版本 (default: '$DEFAULT_SOFT_R_SCRIPT')
        --help/-h                  查看帮助文档
\n" if (defined $if_help or not defined $config or not defined $input_spike or not defined $input_client or not defined $SPIKEIN_DB  or not defined $SPIKEIN_INFO or not defined $output);
$thread        = $DEFAULT_THREAD if (not defined $thread);
$SOFT_R_SCRIPT = $DEFAULT_SOFT_R_SCRIPT   if (not defined $SOFT_R_SCRIPT);
$BLASTN         = $DEFAULT_BLASTN          if (not defined $BLASTN);

my @parameter          = ("F_Primer", "R_Primer", "IR_Copies", "DNA_Quantity", "Amplicon_Area");
my %config             = read_config_with_check($config, \@parameter);
my $amplicon_area      = $config{"Amplicon_Area"}{"value"};
my $f_primer           = $config{"F_Primer"}{"value"};
my $r_primer           = $config{"R_Primer"}{"value"};
my $ir_copies            = exists $config{"IR_Copies"}    ? $config{"IR_Copies"}{"value"}    : "";
my $dna_quantity         = exists $config{"DNA_Quantity"} ? $config{"DNA_Quantity"}{"value"} : "";


###################################################################### 初始化
system "mkdir -p $output" if not -d $output;

###################################################################### 主程序
my $final = "$output/final";system "mkdir -p $final" if not -d $final;
my %spikein_info    = select_spikein_info($SPIKEIN_INFO, $amplicon_area, $f_primer, $r_primer);

my $otu_txt         = "$input_client/otu_table.txt";
my $rep_fasta       = "$input_client/otu.fasta";
my $spike_fasta     = "$input_spike/otu.fasta";
my $spike_otu_txt   = "$input_spike/otu_table.txt";
next if not -e $spike_otu_txt;
die "[ERROR] 不存在文件 $spike_fasta\n"     if not "$spike_fasta";

# 0. 根据样本名提取 ir_copies 及 dna_quantity 
my $spike_IR_copies = "$output/0.spike_IR_copies.txt";
my $DNA_quantity    = "$output/0.DNA_quantity.txt";
extract_ir_copiese($otu_txt, $ir_copies, $spike_IR_copies);
extract_dna_quantity($otu_txt, $dna_quantity, $DNA_quantity);
# 1.blastn 两两比对，找到数据库中与输入序列最相似的序列
my $blastn_result          = "$output/1.spike.blast";
system qq{$BLASTN -query $spike_fasta -db $SPIKEIN_DB -evalue 1e-50 -outfmt 6 -max_target_seqs 100000000 -perc_identity 90 -num_threads $thread > $blastn_result};

# 2.整理统计blastn结果
my $pair_id                = "$output/2.pair_id.txt";  # OTU_ID 和 GS_ID 的对应关系
my $match_spike_otu        = "$output/2.spike_otu.txt";  # 合并匹配到同一内参的OTU序列
my $mismatch_spike_otu     = "$output/2.mismatch_spike_otu.txt";  # 未匹配到内参的OTU序列
my $sample_sum             = "$output/2.sample_sum.txt";  # 统计各样本所有内参序列size数
my $ir_detection_rscript   = "$output/2.ir_detection.r";
ir_detection($blastn_result, $spike_otu_txt, \%spikein_info, $pair_id, $match_spike_otu, $mismatch_spike_otu, $sample_sum, $ir_detection_rscript);

# 3.根据内参添加量，计算拟合曲线  从而计算绝对拷贝数  DNA水平绝对拷贝数  SAMPLE 水平绝对拷贝数
my $spike_reads_percent    = "$output/3.spike_reads_percent.txt";  # 成功匹配内参otu序列数和所有otu序列数统计
my $standard_curve_formula = "$output/3.standard_curve_formula.txt";  # 拟合曲线方程
my $absolute_otu_txt       = "$output/3.absolute_otu.txt";  # 计算绝对拷贝数后otu
my $unit_dna_otu_copies    = "$output/3.unit_dna_otu_copies.txt";  # DNA水平绝对拷贝数
my $unit_sample_otu_copies = "$output/3.unit_sample_otu_copies.txt";  # SAMPLE水平绝对拷贝数
system qq{$SOFT_R_SCRIPT $CALCULATE_ABSOLUTE_COPIES -o $otu_txt -s $match_spike_otu -i $spike_IR_copies -d $DNA_quantity -p 3 --output $output  };
formula_check($standard_curve_formula);

my @commands;
##整理结果
push @commands, ['cp -rf', "$otu_txt",                "$final"];
push @commands, ['cp -rf', "$rep_fasta",              "$final"];
push @commands, ['cp -rf', "$match_spike_otu",        "$final/spike_otu.txt"];
push @commands, ['cp -rf', "$absolute_otu_txt",       "$final/absolute_otu.txt"];
push @commands, ['cp -rf', "$standard_curve_formula", "$final/standard_curve_formula.txt"];
push @commands, ['cp -rf', "$spike_reads_percent",    "$final/spike_reads_percent.txt"];
push @commands, ['cp -rf', "$unit_dna_otu_copies",    "$final/unit_dna_otu_copies.txt"];
push @commands, ['cp -rf', "$unit_sample_otu_copies", "$final/unit_sample_otu_copies.txt"] if -e $unit_sample_otu_copies;

my $err = check_copy_and_run(\@commands);

print "\n[OK] 运行完成\n";

###################################################################### 子程序
sub check_copy_and_run{
    my $command_arr = shift;
    foreach my $command(@$command_arr){
        system qq{$$command[0] $$command[1] $$command[2]};
    }

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
    my ($SPIKEIN_INFO, $amplicon_area, $f_primer, $r_primer) = @_;
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
        if ($data[$index{"region"}] eq $amplicon_area) {
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


sub extract_ir_copiese{
    my ($otu_file, $raw_ir_copies, $target_ir_copies) =@_;
    ##获取需要样本的下标
    open OTU ,$otu_file or die "[ERROR] Cannot open $otu_file\n";
    my $title = <OTU>;
    $title =~ s/[\r\n]//g;
    my @split_title = split/\t/, $title;
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
    print NEW_IR "#OTU ID\t".join("\t", @split_title[1..$#split_title])."\n";
    foreach my $GS(natsort @GBS){
        print NEW_IR $GS;
        foreach my $sample(@split_title[1..$#split_title]){
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
    my @samples = @split_title[1..$#split_title];
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
    foreach my $parameter(@$parameter_list){
         $err .= "[ERROR] 没有定义 $parameter\n" if(not exists $hash{$parameter});
    }

    die $err if($err ne "");
    return %hash;
}
