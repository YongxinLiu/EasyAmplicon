$| = 1; 
# Perl 系统包
use warnings;
use strict;
use File::Spec;
use Getopt::Long;
use Parallel::ForkManager;
use List::Util qw/sum/;

# 常量
use constant SCRIPT_DIR => (File::Spec->splitpath(File::Spec->rel2abs($0)))[1];
use constant PWD => $ENV{"PWD"};

# 变量
my $DEFAULT_THREAD = 40;
my $DEFAULT_MEMORY = 50;
my $DEFAULT_SEPARATE_STRATEGY = "L";
my $DEFAULT_FASTQ_TO_FASTA = "/home/genesky/software/fastx_toolkit/0.0.13/fastq_to_fasta";
my $DEFAULT_SPIKEIN_16S_DB = "/home/genesky/ganb/database_pipeline/mgs/spikein/20250115/16S_SPIKEIN_with_10.fasta";
my $DEFAULT_BLAT           = "/home/genesky/software/blat/37/blat";

# 参数输入
my $ARGV_INFO = join " ", @ARGV;
my ($fastq_dir, $output, $samples, $spiken_db, $self_primer, $no_check, $thread, $SPIKEIN_INFO, $FASTQ_TO_FASTA, $BLAT, $if_help);
GetOptions(
    "fastq_dir|f=s"         => \$fastq_dir,
    "output|o=s"         => \$output,
    "spikein-db=s"       => \$spiken_db,
    "samples|s=s"        => \$samples,
    "self-primer:s"      => \$self_primer,
    "no-check:s"         => \$no_check,
    "thread=i"           => \$thread,
    "spikein-info|si=s"  => \$SPIKEIN_INFO,
    "blat=s"             => \$BLAT,
    "fastq-to-fasta=s"   => \$FASTQ_TO_FASTA,
    "help|h"             => \$if_help,
);
die "
Options: 必填

        --fastq_dir/-f         Fastq_Dir: 原数据路径,fastq.gz
        --output/-o            结果输出路径, intm_result/qc
        --samples/-s           选择指定样本质控, 如'sample1,sample2'
        --spikein-db           spikein数据库,16S:$DEFAULT_SPIKEIN_16S_DB


Options: 可选

        --no-check                 不计算空闲CPU数量
        --thread                   使用线程数 (default: $DEFAULT_THREAD)
        --fastq-to-fasta           更改软件 fastq-to-fasta  版本 (default: '$DEFAULT_FASTQ_TO_FASTA')
        --blat                     更改软件 blat 版本(default: $DEFAULT_BLAT)  
        --help/-h                  查看帮助文档
\n" if (defined $if_help or not defined $fastq_dir or not defined $output or not defined $samples or not defined $spiken_db);
$thread = $DEFAULT_THREAD if (not defined $thread);
$DEFAULT_MEMORY *= $thread/$DEFAULT_THREAD;
$FASTQ_TO_FASTA  = $DEFAULT_FASTQ_TO_FASTA  if (not defined $FASTQ_TO_FASTA);
$BLAT           = $DEFAULT_BLAT             if (not defined $BLAT);

my @samples = split /,/, $samples;


###################################################################### 主程序
my $fasta_output  = "$output/blat_input";  system qq{mkdir -p $fasta_output}  if not -d $fasta_output;
my $blat_output   = "$output/blat_result"; system qq{mkdir -p $blat_output} if not -d $blat_output;
my $blat_exact    = "$output/blat_exact";   system qq{mkdir -p $blat_exact}   if not -d $blat_exact;
my $client_final  = "$output/client_rawdata";                  system qq{mkdir -p $client_final} if not -d $client_final;
my $spike_final   = "$output/spike_rawdata";                   system qq{mkdir -p $spike_final}  if not -d $spike_final;

# 并行序列，进行blat比对
my $pm1 = Parallel::ForkManager -> new($thread);
foreach my $sample(@samples){
    foreach my $direction("R1", "R2"){
        my $pid = $pm1 -> start and next;
        
        ## 1. 数据预处理 fastq -> fasta 并统计原始reads数
        my $fasta_file = "$fasta_output/$sample\_$direction\.fasta";
        system qq{gunzip -c "$fastq_dir/$sample\_$direction\.fastq.gz"|$FASTQ_TO_FASTA -n -v -o $fasta_file -Q 33 &> $fasta_output/$sample\_$direction\_run.log};

        ## 2. blat 比对
        my $blat_result = "$blat_output/$sample\_$direction\.IR.map.xls";
        system qq{$BLAT $spiken_db $fasta_file -out=blast8 -fastMap stdout | awk '\$4 \> 100' > $blat_result};        
        $pm1 -> finish;
    }
}
$pm1 -> wait_all_children;

my @check_spike = `grep ">" $spiken_db | awk -F " " '{print \$1}'| sed 's/>//'|sort -V`; map{$_ =~ s/[\r\n]//g}@check_spike;
my %hashHandle;
for("R1", "R2"){
    open $hashHandle{$_}, ">$blat_exact/spikein_$_\_stat.txt";
    print {$hashHandle{$_}} "sample_id\tdirection\tspikein_id\tspikein_num\n";
}

# 并行样本，对blat结果进行读取和筛选（只要有一端符合筛选条件，则认为该条序列为spikein序列, 须同时考虑R1, R2端）
my $pm2 = Parallel::ForkManager -> new($thread);

foreach my $sample(@samples){
    my $pid = $pm2 -> start and next;
    
    my $blat_R1        = "$blat_output/$sample\_R1.IR.map.xls";
    my $blat_R2        = "$blat_output/$sample\_R2.IR.map.xls";
    my $select_IR_map_R1 = "$blat_exact/select_$sample\_R1.IR.map.xls";
    my $select_IR_map_R2 = "$blat_exact/select_$sample\_R2.IR.map.xls";

    ## 3.根据数据拆分策略loose/strict 提取符合条件的序列标记为内参序列
    my (%hashExact, %count_spike);
    
    loose_extrac_id($blat_R1, $select_IR_map_R1, 96, 200, \%hashExact, \%count_spike, "R1");
    loose_extrac_id($blat_R2, $select_IR_map_R2, 96, 200, \%hashExact, \%count_spike, "R2");
    
    ## 4. 根据标记的内参序列对原数据进行拆分
    split_data("$fastq_dir/$sample\_R1.fastq.gz", \%hashExact, "$client_final/$sample\_R1.fastq.gz", "$spike_final/$sample\_R1.fastq.gz");
    split_data("$fastq_dir/$sample\_R2.fastq.gz", \%hashExact, "$client_final/$sample\_R2.fastq.gz", "$spike_final/$sample\_R2.fastq.gz");
    
    ## 5. 统计输出各样本R1/R2端 各内参reads 并检测有无添加量为0的内参
    my $err = "";
    foreach my $Rtype("R1", "R2"){
        foreach my $spikein_id(@check_spike){
            if(exists $count_spike{$Rtype}{$spikein_id}){
                print {$hashHandle{$Rtype}} "$sample\t$Rtype\t$spikein_id\t$count_spike{$Rtype}{$spikein_id}\n";
            }else{
                print {$hashHandle{$Rtype}} "$sample\t$Rtype\t$spikein_id\t0\n";
                $err .= "[ERROR] $sample\t$Rtype\t$spikein_id spikein数量为 0\n";
            }
        }
    }
    print $err;
    
    $pm2 -> finish;
}
$pm2 -> wait_all_children;

for("R1","R2"){
    close $hashHandle{$_}
}

print "[OK] 运行正常\n";

################################################################################################


sub split_data{
    my($raw_fastq, $select_hash, $client_fastq, $spike_fastq) = @_;
    
    open GZ, qq{zcat $raw_fastq |};
    open CLIENT, "|gzip > $client_fastq" or die "[ERROR] Cannot write into $client_fastq\n";
    open SPIKE,  "|gzip > $spike_fastq"  or die "[ERROR] Cannot write into $spike_fastq\n";
    while(my $id = <GZ>){
        my $seq = <GZ>;
        my $plus = <GZ>;
        my $quality = <GZ>;
        my $ir_id  = (split/\s/, $id)[0];
        $ir_id =~ s/\/1$//;
        $ir_id =~ s/\/2$//;
        if(exists $select_hash -> {$ir_id}){
            next if($select_hash -> {$ir_id} eq "Gs_BSI10" or $select_hash -> {$ir_id} eq "Gs_ITS-SPK-Long");  # 剔除内参Gs_BSI10, Gs_ITS-SPK-Long
            print SPIKE "$id$seq$plus$quality";
        }else{
            print CLIENT "$id$seq$plus$quality";
        }
    }
    close GZ;
    close CLIENT;
    close SPIKE;
}


sub loose_extrac_id{
    my($total_IR_map, $select_IR_map, $identity, $length, $hash_exact, $hash_count, $Rtype) = @_;
    open IR, $total_IR_map or die "[ERROR] Cannot open $total_IR_map\n";
    open SELECT, ">$select_IR_map" or die "[ERROR] Cannot write into $select_IR_map\n";
    
    while(<IR>){
        $_ =~ s/[\r\n]//g;
        my @tmp = split /\t/,$_;
        $tmp[0] =~ s/\/1$//;
        $tmp[0] =~ s/\/2$//;
        my $id  = "@".$tmp[0];
        if($tmp[1] eq "Gs_BSI10" and $tmp[2] >= $identity and $tmp[3] >= $length and $tmp[6] < 30 and $tmp[8] > $tmp[9]){
            $hash_exact -> {$id} = $tmp[1] if not exists $hash_exact -> {$id};  # 只存入一次，以R1端为最终判断内参
            $hash_count -> {$Rtype}{$tmp[1]}++;
            print SELECT "$_\n";
        }
        
        if($tmp[1] ne "Gs_BSI10" and $tmp[2] >= $identity and $tmp[3] >= $length){
            $hash_exact -> {$id} = $tmp[1] if not exists $hash_exact -> {$id};
            $hash_count -> {$Rtype}{$tmp[1]}++;
            print SELECT "$_\n";
        }
    }
    
    close IR;
    close SELECT;
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
