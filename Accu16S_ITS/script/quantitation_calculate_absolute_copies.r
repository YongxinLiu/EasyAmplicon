#!/home/genesky/software/r/4.0.3/bin/Rscript
options(scipen=200)
library("argparse")
# 参数解析
parse_args <- function(){

    parser <- ArgumentParser(description='绝对拷贝数计算', formatter_class="argparse.RawTextHelpFormatter", epilog=" 结果文件包括： \\n prefix.spike_reads_persent.txt 内参占比统计 \\n prefix.standard_curve_formula.txt 标准曲线方程 \\n prefix.absolute_otu.txt 绝对拷贝数 \\n prefix.unit_dna_otu_copies.txt 每单位DNA绝对拷贝数 \\n prefix.unit_sample_otu_copies.txt 每单位样本绝对拷贝数 \\n")

    parser$add_argument('-o', '--otufile', type='character', metavar='otufile', required=TRUE,
        help = "otu文件")

    parser$add_argument('-s', '--spikein_otufile', type='character', metavar='spikein_otufile', required=TRUE,
        help = "内参otu文件 spikein_otu.txt")

    parser$add_argument('-i', '--ir_copiesfile', type='character', metavar='ir_copiesfile', required=TRUE,
        help = "IR_copies文件,里面包含样本spikein的添加量信息,实验提供")

    parser$add_argument('-d', '--dna_quantityfile', type='character', metavar='dna_quantityfile', required=TRUE,
        help = "DNA_quantity文件,里面包含各样本DNA上机量, DNA抽提量和抽提DNA所用样本量, 实验提供")

    parser$add_argument('--output', type='character', metavar='output dir', required=TRUE,
        help = " 结果输出目录")

    parser$add_argument('-p', '--prefix', type='character', metavar='prefix', required=TRUE,
        help = "结果文件前缀")

    parser$add_argument('--drop_spikein', type='character', metavar='eg. Gs_ITS-SPK2,Gs_ITS-SPK4,Gs_ITS-SPK6', default="null",
        help = "在绝对定量计算过程中，线性拟合不纳入的内参，多条内参用英文逗号 , 隔开 [default: %(default)s] ")

    args <- parser$parse_args()

    return(args)
}

# 保存命令行信息，方便debug
OUT_COMMAND <- function(log_file){
    # 1. 时间
    time_now = format(Sys.time(), "%Y-%m-%d %H:%M:%OS")
    # 2. R路径
    rscript = gsub('/lib64/R/bin/exec/R$', '/bin/Rscript', commandArgs(trailingOnly = FALSE)[1], perl=TRUE)
    # 3. R脚本路径
    r_script = ""
    for(arg in commandArgs(trailingOnly = FALSE)){
        if(grepl('^--file=', arg, perl=TRUE)){
            r_script = gsub('^--file=', '', arg, perl=TRUE)
            break
        }
    }
    # 4. 参数组合
    parameter = c()
    for(arg in commandArgs(trailingOnly = TRUE)){
        if(grepl('^-', arg, perl=TRUE)){
            parameter = c(parameter, arg)
        }else{
            # 所有的值添加 单引号，防止有特殊字符导致传参异常
            parameter = c(parameter, paste0('\'', arg, '\''))
        }
    }

    # 5. 输出
    # 命令行
    command = paste(c(rscript, r_script, parameter), collapse = " ")
    # 测试用 argparse
    command_argparse = paste(commandArgs(trailingOnly = TRUE), collapse='\',\'')
    command_argparse = paste0('args <- parser$parse_args(c(\'', command_argparse, '\'))')
    
    write(paste0("[Command]: ", time_now, '\n', command, '\n', command_argparse, '\n\n'), file = log_file, append = TRUE)
}

args = parse_args()
if(!file.exists(args$output)){
    dir.create(args$output)
}

OUT_COMMAND(paste0(args$output, '/command.info'))

otufile          = normalizePath(args$otufile)
spikein_otufile  = normalizePath(args$spikein_otufile)
ir_copiesfile    = normalizePath(args$ir_copiesfile)
dna_quantityfile = normalizePath(args$dna_quantityfile)
output           = normalizePath(args$output)

setwd(output)


###导入数据
otu_data          = read.table(otufile, header = TRUE, row.names = 1, check.name = F, quote = "", comment.char = "", sep = "\t", fill = T)
spikein_data      = read.table(spikein_otufile, header = TRUE, row.names = 1, check.name = F, quote = "", comment.char = "", sep = "\t", fill = T)
ir_copies_data    = read.table(ir_copiesfile, header = TRUE, row.names = 1, check.name = F, quote = "", comment.char = "", sep = "\t", fill = T)
dna_quantity_data = read.table(dna_quantityfile, header = TRUE, row.names = 1, check.name = F, quote = "", comment.char = "", sep = "\t", fill = T, colClasses = "character")

# 去除线性拟合指定不纳入指定的spike
if(args$drop_spikein != "null"){
    drop_spikein = gsub(" ", "", args$drop_spikein)  # 去除空格
    drop_spikein = unlist(strsplit(drop_spikein, ","))
}else{
    drop_spikein = c()
}

####提取相应数据并统一顺序
samples               = colnames(spikein_data)
sub_spikein_data      = spikein_data[, samples, drop = F]
sub_otu_data          = otu_data[, samples, drop = F]
sub_ir_copies_data    = ir_copies_data[rownames(sub_spikein_data), samples, drop = F]
sub_dna_quantity_data = t(dna_quantity_data[samples, ,drop = F])

###统计占比
stat_data = rbind(colSums(sub_spikein_data) + colSums(sub_otu_data),
                  colSums(sub_spikein_data), 
                  round(100*colSums(sub_spikein_data)/(colSums(sub_spikein_data) + colSums(sub_otu_data)), 2))
stat_data = rbind(colnames(stat_data),stat_data)
rownames(stat_data) = c("Sample","Totle_reads","Spikein_reads", "Spikein_percent(%)")
write.table(t(stat_data), file = paste(args$prefix,".spike_reads_percent.txt", sep = ""), sep = "\t",row.names = F, col.names = T, quote = F)

########计算
STAT <- function(ir_copies, spike_otu, otu, type, drop_spikein){

    result <- sapply(1:ncol(ir_copies), function(t){

        keep = (spike_otu[,t] > 0) & (!rownames(spike_otu) %in% drop_spikein)
        if(sum(keep) == 0){
            message("存在样本所有内参含量为0，无法计算绝对拷贝数", colnames(ir_copies)[t])
            q()
        }
        x = log10(ir_copies[keep, t])  # spikein 添加量
        y = log10(spike_otu[keep, t])  # spikein 检出量
        lm = summary(lm(y~x))

        ##方程  y = a*x + b，r为r2拟合值
        if(length(x) == 1){
            a = y/x;
            b = 0;
        }else{
            b = round(lm$coefficients[[1]], 4)
            a = round(lm$coefficients[[2]], 4)
        }
        r = lm$r.squared
        if (b<0){
            combine = paste("y=", a, "x", b, sep = "")
        }else{
            combine = paste("y=", a, "x+", b, sep = "")
        }
        formula = c(colnames(ir_copies)[t], combine, r)

        ##计算绝对拷贝数
        otulog = log10(otu[, t])  
        mi = (otulog-b)/a
        absolute_otu = round(10^mi)

        ##选择性返回结果
        if(type == "formula"){
            return(formula)
        }else{
            return(absolute_otu)
        }    
    })
    return(result)
}
###=============算方程=====================
formula      = t(STAT(sub_ir_copies_data, sub_spikein_data, sub_otu_data, "formula", drop_spikein))  
colnames(formula) = c("sample", "linear equation", "R2")
write.table(formula, file = paste(args$prefix,".standard_curve_formula.txt", sep = ""), quote = FALSE, sep = "\t", row.names = F)

####===========算绝对拷贝数================
absolute_otu = data.frame(STAT(sub_ir_copies_data, sub_spikein_data, sub_otu_data, "copies", drop_spikein))
colnames(absolute_otu) = samples


absolute_otu_result = data.frame(taxo = rownames(otu_data), absolute_otu)
colnames(absolute_otu_result) = as.character(c("#OTU ID",colnames(absolute_otu)))
rownames(absolute_otu_result) = rownames(otu_data)

write.table(absolute_otu_result, file = paste(args$prefix,".absolute_otu.txt", sep = ""), quote = FALSE, sep = "\t", row.names = F)

####=========== 绝对量的单位计算================
sub_absolute_otu = absolute_otu_result[ , samples, drop = F]

### DNA水平绝对拷贝数：copies/DNA_ng
result_dna <- sapply(1:ncol(sub_dna_quantity_data), function(j) {
    template_DNA = as.numeric(sub_dna_quantity_data[1,j])   ##模板DNA量
    template_copies = as.numeric(sub_absolute_otu[,j])    ##绝对拷贝数
    dna_copies = round(template_copies/template_DNA)
    return(dna_copies)
})
colnames(result_dna) = samples

result_dna_data = data.frame(taxo = rownames(sub_absolute_otu), result_dna)
colnames(result_dna_data) = as.character(c("#OTU ID", colnames(result_dna)))
write.table(result_dna_data , file = paste(args$prefix,".unit_dna_otu_copies.txt", sep = ""), sep = "\t",row.names = F,  quote = FALSE)

### SAMPLE 水平绝对拷贝数：copies/quantity_g
if(ncol(dna_quantity_data) == 3){
    result_sample <- sapply(1:ncol(sub_dna_quantity_data), function(j) {
        template_DNA = as.numeric(sub_dna_quantity_data[1,j])   ##模板DNA量
        template_copies = as.numeric(sub_absolute_otu[,j])    ##绝对拷贝数
        sample_DNA = as.numeric(sub_dna_quantity_data[2,j])           ##DNA总量
        sample_quantity = as.numeric(sub_dna_quantity_data[3,j])      ##所用样本量
        sample_copies = round(template_copies / (template_DNA * (sample_quantity / sample_DNA)))
        #sample_copies = round(((template_copies*sample_DNA)/template_DNA)/sample_quantity)
        return(sample_copies)
    })
    colnames(result_sample) = samples

    result_sample_data = data.frame(taxo = rownames(sub_absolute_otu), result_sample)
    colnames(result_sample_data) = as.character(c("#OTU ID", colnames(result_sample)))
    write.table(result_sample_data , file = paste(args$prefix,".unit_sample_otu_copies.txt",sep = ""), sep = "\t",row.names = F,  quote = FALSE)
}
