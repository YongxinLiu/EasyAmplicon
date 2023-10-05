[TOC]

# For amplicon sequencing EasyAmplicon

    # Authors: Yong-Xin Liu(刘永鑫), Tong Chen(陈同)et al.
    # Version: v1.20
    # Update: 2023-10-13
    # System requirement: Windows 10+ / Mac OS 10.12+ / Ubuntu 20.04+
    # Reference: Liu, et al. 2023. EasyAmplicon: An easy-to-use, open-source, reproducible, and community-based
    # pipeline for amplicon data analysis in microbiome research. iMeta 2: e83. https://doi.org/10.1002/imt2.83


    # Setting working directory (work directory, wd) and software database directory (database, db)
    # Add environmental variables and enter the work directory 
    # Every time you open RStudio, you must run the following 4 lines. Run it. You can optionally replace ${db} with the installation location of EasyMicrobiome.
    wd=/c/EasyAmplicon
    db=/c/EasyMicrobiome
    PATH=$PATH:${db}/win
    cd ${wd}


## 1. start files

    # 1. Analysis Process pipeline.sh
    # 2. Sample metadata "metadata.txt" saved in the "result" directory
    # 3. Sequencing data fastq files saved in the "seq" directory, typically ending with .fq.gz, with one pair of files per sample
    # 4. Create a temporary file storage directory that can be deleted after the analysis is completed
    mkdir -p seq result temp 

### 1.1. metadata/Experimental design

    # Prepare sample metadata in "result/metadata.txt
    # Use "csvtk" to count the number of rows (samples, excluding the header) and columns in a table. Use the -t option to set the column delimiter as a tab (default is comma);
    csvtk -t stat result/metadata_raw.txt
    # The metadata should have at least 3 columns. The first column should be the sample ID (SampleID), and the last column should be the description (Description)
    # You can use the cat command to view a file, and the -A option to display symbols. The '|' symbol is used as a pipe to connect commands. The head command is used to display the file header, and the -n3 option limits the output to the first 3 lines.
    cat -A result/metadata_raw.txt | head -n3
    # Windows users end with ^M, run the sed command to remove, and then check with cat -A
    sed 's/\r//' result/metadata_raw.txt > result/metadata.txt
    cat -A result/metadata.txt | head -n3

### 1.2. sequencing data

    # # This code can be run after Ctrl + Shift + C uncomment "#" in RStudio
    # # Optionally, download sequencing data, by GSA's CRA (batch) and CRR (sample) number
    # # The sample downloads a single file and renames it
    # mkdir -p seq
    # wget -c ftp://download.big.ac.cn/gsa/CRA002352/CRR117575/CRR117575_f1.fq.gz -O seq/KO1_1.fq.gz
    # # Download and rename in bulk by experimental design number
    # awk '{system("wget -c ftp://download.big.ac.cn/gsa/"$5"/"$6"/"$6"_f1.fq.gz -O seq/"$1"_1.fq.gz")}' \
    #     <(tail -n+2 result/metadata.txt)
    # awk '{system("wget -c ftp://download.big.ac.cn/gsa/"$5"/"$6"/"$6"_r2.fq.gz -O seq/"$1"_2.fq.gz")}' \
    #     <(tail -n+2 result/metadata.txt)

    # The sequencing results returned by the company are usually a pair of FQ/FastQ compressed files in .gz format for one sample
    # The file name and the sample name must be corresponding: manually modify it when it is inconsistent, and see "FAQ 6" for batch name change
    # If the sequencing data is a compressed file of .gz, sometimes you need to use gunzip to decompress and use, vsearch can usually read the compressed file directly
    # gunzip seq/*.gz
    # zless view compressed files by page, space page turn, q exit; head looks at the first 10 rows by default, and -n specifies the rows
    ls -sh seq/
    zless seq/KO1_1.fq.gz | head -n4 
    # Each line is too long, specify to view 1-60 characters per line
    zless seq/KO1_1.fq | head | cut -c 1-60
    #  sequencing data, relying on the seqkit program
    seqkit stat seq/KO1_1.fq.gz
    # Batch count sequencing data and summarize tables
    seqkit stat seq/*.fq.gz > result/seqkit.txt
    head result/seqkit.txt

### 1.3.Pipeline & databases

    # The database must be decompressed the first time it is used, and you can skip this section later

    # usearches available 16S/18S/ITS databases: RDP, SILVA and UNITE, local file location ${db}/usearch/
    # Usearch database database download page: http://www.drive5.com/usearch/manual/sintax_downloads.html
    # Decompress 16S RDP database, gunzip decompress, seqkit stat statistics
    #  Keep the original compressed file
    gunzip -c ${db}/usearch/rdp_16s_v18.fa.gz > ${db}/usearch/rdp_16s_v18.fa
    seqkit stat ${db}/usearch/rdp_16s_v18.fa # 2.1万sequences
    # To decompress the ITS UNITE database, you need to download it from the official website or the network disk db/amplicon/usearch
    # gunzip -c ${db}/usearch/utax_reference_dataset_all_29.11.2022.fasta.gz >${db}/usearch/unite.fa
    seqkit stat ${db}/usearch/unite.fa # 32.6万
    # The Greengene database is used for feature annotations: ftp://greengenes.microbio.me/greengenes_release/gg_13_5/gg_13_8_otus.tar.gz
    # The default decompression will delete the original file, -c specifies the output to the screen, > Write a new file (can be renamed):
    gunzip -c ${db}/gg/97_otus.fasta.gz > ${db}/gg/97_otus.fa
    seqkit stat ${db}/gg/97_otus.fa


## 2. Reads merge and rename

### 2.1 Merge pair-end reads and rename

    #Test. Take the combination of WT1 single samples as an example
    #Time counts the calculation time, real is the physical time, user is the calculation time, sys is the hardware wait time
    # time vsearch --fastq_mergepairs seq/WT1_1.fq.gz \
    #   --reverse seq/WT1_2.fq.gz \
    #   --fastqout temp/WT1.merged.fq \
    #   --relabel WT1.
    # head temp/WT1.merged.fq

    #Batch batches and merging according to the experimental design
    #tail -n+2 to the table header, cut -f1 to take the first column, get the sample list; 18 samples x 15,000 pairs of sequences combined for 8 s
    #Copy Ctrl+C under Win to terminate under Linux, in order to prevent abnormal interruption, add & turn background at the end, press Enter to continue after no display
    
    # Some computers are not supported by rush, and the runtime scheduling fails, please use the for loop part
    # The looping section is placed in the background for execution. After clicking 'run,' the program appears to have completed running, but it hasn't actually finished; it is still running。
    # Take your time before running the following program。
    # In the previous lessons, it was noticed that the results were different each time the program was run. This was because the 'for' loop section didn't complete its execution, only generating partial data, which affected the subsequent steps
    # The number of reads for each sample will be inconsistent with each run。
    #Method 1: Sequential Processing with a 'for' Loop
    # time for i in `tail -n+2 result/metadata.txt|cut -f1`;do
    #   vsearch --fastq_mergepairs seq/${i}_1.fq.gz --reverse seq/${i}_2.fq.gz \
    #   --fastqout temp/${i}.merged.fq --relabel ${i}.
    # done &

    # Some computers do not support rush, and scheduling fails during execution. Please use the 'for' loop section instead
    #Method 2: Parallel Processing with 'rush,' utilizing job count (j). With 2 jobs, speedup by 1x, 4 seconds; it is recommended to set it between 2-4
    time tail -n+2 result/metadata.txt | cut -f 1 | \
     rush -j 3 "vsearch --fastq_mergepairs seq/{}_1.fq.gz --reverse seq/{}_2.fq.gz \
      --fastqout temp/{}.merged.fq --relabel {}."
    # Check the sample names in the first 10 lines of the last file
    head temp/`tail -n+2 result/metadata.txt | cut -f 1 | tail -n1`.merged.fq | grep ^@
    
    ##Method 3: If compression is not supported, decompress and then merge in paired-end mode
    #  # gunzip seq/*.fq.gz
    #  time tail -n+2 result/metadata.txt | cut -f 1 | \
    #    rush -j 1 "vsearch --fastq_mergepairs <(zcat seq/{}_1.fq.gz) --reverse <(zcat seq/{}_2.fq.gz) \
    #     --fastqout temp/{}.merged.fq --relabel {}."
    # 
    #   time for i in `tail -n+2 result/metadata.txt|cut -f1`;do
    #      vsearch --fastq_mergepairs <(zcat seq/${i}_1.fq.gz) --reverse <(zcat seq/${i}_2.fq.gz) \
    #      --fastqout temp/${i}.merged.fq --relabel ${i}.
    #    done &
      
### 2.2 (Optional) Rename single-end files

    # # Example of renaming a single sequence
    # i=WT1
    # gunzip -c seq/${i}_1.fq.gz > seq/${i}.fq
    # usearch -fastx_relabel seq/${i}.fq -fastqout temp/${i}.merged.fq -prefix ${i}.
    # 
    # # To perform batch renaming, you need uncompressed single-end FASTQ files (usearch does not support compressed formats)
    # gunzip seq/*.gz
    # time for i in `tail -n+2 result/metadata.txt|cut -f1`;do
    #   usearch -fastx_relabel seq/${i}.fq -fastqout temp/${i}.merged.fq -prefix ${i}.
    # done &
    # # Refer to 'FAQ 2' for the large data processing method using VSEARCH”

### 2.3 integrate sequences after renaming 

    #Merge all samples into a single file
    cat temp/*.merged.fq > temp/all.fq
    #View file size: 223MB. There are slight differences in results due to different software versions.
    ls -lsh temp/all.fq
    # Check if the sequence name before the '.' is the sample name. Sample names must not contain dots ('.')
    # One notable characteristic of sample names with dots (.) is that the generated feature table becomes significantly large, with numerous columns, leading to insufficient memory during subsequent analyses。
    # After obtaining the feature table from the subsequent analysis, it's important to take a look to check for any issues. If you encounter memory-related problems, you'll need to go back and troubleshoot。
    head -n 6 temp/all.fq|cut -c1-60


## 3. Trim primers and perform quality control 

    # The left end has a 10bp tag + 19bp upstream primer V5, totaling 29bp. The right end has the V7 downstream primer of 18bp
    # Cut barcode 10bp + V5 19bp in left and V7 18bp in right
    # Be sure to understand the experimental design and primer lengths. If primers have been removed, you can enter 0. Processing 270,000 sequences will take around 14 seconds
    time vsearch --fastx_filter temp/all.fq \
      --fastq_stripleft 29 --fastq_stripright 18 \
      --fastq_maxee_rate 0.01 \
      --fastaout temp/filtered.fa
    # View the file to understand the FASTA file format
    head temp/filtered.fa


## 4. Remove redundancy and select OTU/ASV 

### 4.1 Sequence de-redundancy 

    # And add a miniuniqusize of at least 8 or 1/1M to remove low-abundance noise and increase calculation speed
    # -sizeout outputs abundances, and --relabel must be used with sequence prefixes for better organization. Takes 1 second
    vsearch --derep_fulllength temp/filtered.fa \
      --minuniquesize 10 --sizeout --relabel Uni_ \
      --output temp/uniques.fa 
    #High-abundance non-redundant sequences are very small (500K~5M is suitable), with size and frequency appended to the names
    ls -lsh temp/uniques.fa
    # Uni_1;size=6423" - The name of the sequence after dereplication. This sequence appeared 6423 times in the sequencing data across all samples
    # is the sequence that appears the most。
    head -n 2 temp/uniques.fa

### 4.2 Cluster OTUs or denoise ASVs

    #There are two methods: It is recommended to use unoise3 for denoising to obtain ASVs with single-base accuracy. Alternatively, you can use the traditional 97% clustering for OTUs (genus-level accuracy) as an alternative
    #Both feature selection methods in usearch come with built-in de novo chimera removal.
    #-minsize secondary filtering, control the number of OTU/ASV to 1-5,000, convenient for downstream statistical analysis

    #Method 1: 97% clustering for OTUs, suitable for large datasets/when ASV patterns are not clear/reviewer's request
    #The process took 1 second, resulting in 508 OTUs with 126 chimeras removed.
    # usearch -cluster_otus temp/uniques.fa -minsize 10 \
    #  -otus temp/otus.fa \
    #  -relabel OTU_

    #Method 2. ASV Denoise: predict biological sequences and filter chimeras
    #In 6 seconds, 1530 sequences were identified as good, with 41 chimeras removed. For datasets with millions of sequences, the process might take several days or weeks.
    usearch -unoise3 temp/uniques.fa -minsize 10 \
      -zotus temp/zotus.fa
    #Modify sequence names: Change "Zotu" to "ASV" for easier identification
    sed 's/Zotu/ASV_/g' temp/zotus.fa > temp/otus.fa
    head -n 2 temp/otus.fa

    #
Method 3: When the data is too large to use usearch, consider the alternative vsearch method as described in "Common Issue 3"

### 4.3 Perform reference-based chimera detection.

    # Not recommended, as it can lead to false negatives due to the absence of abundance information in the reference database
    # In contrast, de novo chimera detection requires parent abundance to be at least 16 times that of the chimera to prevent false negatives
    # Since known sequences will not be removed, choosing a larger reference database is more appropriate to minimize the false negative rate
    mkdir -p result/raw

    # Method 1: Using vsearch in combination with RDP for chimera removal (fast but prone to false negatives)
    # You can download and extract the SILVA database (replacing rdp_16s_v18.fa with silva_16s_v123.fa). This process is very slow, but theoretically provides better results.
    vsearch --uchime_ref temp/otus.fa \
      -db ${db}/usearch/rdp_16s_v18.fa \
      --nonchimeras result/raw/otus.fa
    # RDP: 7s, 143 (9.3%) chimeras; SILVA：9m, 151 (8.4%) chimeras
    # If you're using Windows and vsearch results have added Windows line endings (^M), you need to remove them. However, if you're using a Mac, you don't need to execute this command
    sed -i 's/\r//g' result/raw/otus.fa

    # Method 2: Do not perform chimera removal
    # cp -f temp/otus.fa result/raw/otus.fa


## 5. Feature table construction and filtering

    # OTUs and ASVs are collectively referred to as "features". Their differences are as follows：
    # OTUs are typically representative sequences selected from the highest abundance or centroid after 97% clustering；
    # ASVs are representative sequences based on sequence denoising (exclusion or correction of erroneous sequences, selecting credible sequences with higher abundance)

### 5.1 Generate the feature table

    # Method 1: Generate feature table using usearch. It's fast for small samples (<30), but for large samples, it's limited and multi-threading efficiency is low (83.2%, 4 cores, 17 seconds)
    # time usearch -otutab temp/filtered.fa \
    #   -otus result/raw/otus.fa \
    #   -threads 4 \
    #   -otutabout result/raw/otutab.txt

    # Method 2: Generate feature table using vsearch
    # id(1): 49.45% of sequences aligned with 100% similarity, taking 1 minute and 50 seconds.
    # id(0.97): 83.66% of sequences aligned with 97% similarity, taking 1 minute and 10 seconds. (Higher data utilization leads to faster results.)
    time vsearch --usearch_global temp/filtered.fa \
      --db result/raw/otus.fa \
      --id 0.97 --threads 4 \
    	--otutabout result/raw/otutab.txt 
    #212,862 out of 268,019 (79.42%) are alignable
    # For Windows users, you should remove the newline characters (^M) from vsearch results and correct them to the standard Linux format.
    sed -i 's/\r//' result/raw/otutab.txt
    head -n6 result/raw/otutab.txt | cut -f 1-6 |cat -A
    # To count rows and columns in a table using "csvtk"
    # Here, it's important to carefully check the number of columns and whether it matches the number of your samples. If they're not equal, there might be issues with sample naming. Refer to the explanations above for more details.
    csvtk -t stat result/raw/otutab.txt


### 5.2  Species annotation, and/or removal of plastids and non-bacteria 
    # Species Annotation - Remove plastids and non-bacterial/archaea and count the proportions (optional)
    # RDP species annotation (rdp_16s_v18) is faster, but lacks complete eukaryotic source data, which may be incomplete and take 15s;
    # SILVA database (silva_16s_v123.fa) is better at annotating eukaryotic and plastid sequences, which is extremely slow and takes up to 3h
    # The confidence threshold is typically 0.6/0.8, with a minimum of 0.1/usearch for vserch, optional, and 0 outputs most similar species annotations for observing potential taxonomy
    vsearch --sintax result/raw/otus.fa \
      --db ${db}/usearch/rdp_16s_v18.fa \
      --sintax_cutoff 0.1 \
      --tabbedout result/raw/otus.sintax 
    head result/raw/otus.sintax | cat -A
    sed -i 's/\r//' result/raw/otus.sintax

    # Method1. The number of original feature table rows
    wc -l result/raw/otutab.txt
    #R script selects bacterial archaea (eukaryotes), removes chloroplasts, mitochondria and counts the proportions; Output filtered and sorted OTU tables
    #The input is the OTU table result/raw/otutab.txt and the species annotation result/raw/otus.sintax
    #Output filtered and sorted feature tables result/otutab.txt sum
    #Statistical contamination ratio file result/raw/otutab_nonBac.txt and filter details otus.sintax.discard
    #For fungal ITS data, use the otutab_filter_nonFungi.R script instead, only screening fungi
    # Rscript ${db}/script/otutab_filter_nonBac.R -h # Displays a description of the parameter
    Rscript ${db}/script/otutab_filter_nonBac.R \
      --input result/raw/otutab.txt \
      --taxonomy result/raw/otus.sintax \
      --output result/otutab.txt\
      --stat result/raw/otutab_nonBac.stat \
      --discard result/raw/otus.sintax.discard
    # The number of filtered feature table rows
    wc -l result/otutab.txt
    #Filter the sequence corresponding to the feature table
    cut -f 1 result/otutab.txt | tail -n+2 > result/otutab.id
    usearch -fastx_getseqs result/raw/otus.fa \
        -labels result/otutab.id -fastaout result/otus.fa
    #Filter feature tables corresponding to sequence annotations
    awk 'NR==FNR{a[$1]=$0}NR>FNR{print a[$1]}'\
        result/raw/otus.sintax result/otutab.id \
        > result/otus.sintax

    # Method2. If you feel that the screening is unreasonable, you can not screen
    # cp result/raw/otu* result/

    #Optional statistical method: Summary OTUs table
    usearch -otutab_stats result/otutab.txt \
      -output result/otutab.stat
    cat result/otutab.stat
    #Note the minimum, quantile, or look at the amount of sample detail in result/raw/otutab_nonBac.stat for resampling

### 5.3 Equalization of sampling

    # Normlize by subsample

    #Use the vegan package for equal resampling and enter the reads count format Feature table result/otutab.txt
    #You can specify the input file, sample size, and random number, and output the equalizer result/otutab_rare.txt and diversity alpha/vegan .txt
    mkdir -p result/alpha
    Rscript ${db}/script/otutab_rare.R --input result/otutab.txt \
      --depth 10000 --seed 1 \
      --normalize result/otutab_rare.txt \
      --output result/alpha/vegan.txt
    usearch -otutab_stats result/otutab_rare.txt \
      -output result/otutab_rare.stat
    cat result/otutab_rare.stat


## 6. alpha (α) diversity

### 6.1. calculate alpha (α) diversity

    # Use USEARCH to calculate 14 alpha diversity indices (don't use Chao1 if you make a mistake
    #details in http://www.drive5.com/usearch/manual/alpha_metrics.html
    usearch -alpha_div result/otutab_rare.txt \
      -output result/alpha/alpha.txt

### 6.2. calculate rarefaction richness

    #Dilution Curve: Take the number of OTUs (Operational Taxonomic Units) from sequences ranging from 1% to 100%, with non-replacement sampling each time.
    #Rarefaction from 1%, 2% .. 100% in richness (observed OTUs)-method without_replacement https://drive5.com/usearch/manual/cmd_otutab_subsample.html
    usearch -alpha_div_rare result/otutab_rare.txt \
      -output result/alpha/alpha_rare.txt \
      -method without_replacement
    #Preview of Results
    head -n2 result/alpha/alpha_rare.txt
    #Handling of non-numeric "-" values due to low sample sequencing, see FAQ #8 for details
    sed -i "s/-/\t0.0/g" result/alpha/alpha_rare.txt

### 6.3. Filtering for high-abundance bacteria

    #Calculate the mean of each feature, and if there are groups, calculate the group means. Group column names should be modified based on the experimental design in metadata.txt.
    #The input files are the feature table (result/otutab.txt) and the experimental design metadata (metadata.txt)
    #The output is the feature table with group-wise means. A single experiment may have multiple grouping methods
    #The "-h" option displays the script's help, providing explanations for the parameters
    Rscript ${db}/script/otu_mean.R -h
    #The parameter "scale" determines whether to perform scaling, "zoom" normalizes the total sum, "all" outputs mean for all samples, and "type" specifies the calculation type as "mean" or "sum"
    Rscript ${db}/script/otu_mean.R --input result/otutab.txt \
      --metadata result/metadata.txt \
      --group Group --thre 0 \
      --scale TRUE --zoom 100 --all TRUE --type mean \
      --output result/otutab_mean.txt
    # The output includes both the overall mean and the mean values for each group
    head -n3 result/otutab_mean.txt

    #By filtering with an average abundance threshold of >0.1%, which can be selected as 0.5% or 0.05%, you will obtain the OTU combinations for each group
    awk 'BEGIN{OFS=FS="\t"}{if(FNR==1) {for(i=3;i<=NF;i++) a[i]=$i; print "OTU","Group";} \
        else {for(i=3;i<=NF;i++) if($i>0.1) print $1, a[i];}}' \
        result/otutab_mean.txt > result/alpha/otu_group_exist.txt
    head result/alpha/otu_group_exist.txt
    cut -f 2 result/alpha/otu_group_exist.txt | sort | uniq -c
    # Give it a try: How many OTUs/ASVs does each group have at different abundances?
    # Using the following link http://ehbio.com/test/venn/ You can plot and visualize Venn or network diagrams for shared and unique components of each group at different abundances 
    # Using the following link http://www.ehbio.com/ImageGP Create Venn, UpSetView, and Sankey diagrams

## 7. Beta(β) diversity

    #The results generate multiple files and require a directory
    mkdir -p result/beta/
    #Building an evolutionary tree based on OTUs
    usearch -cluster_agg result/otus.fa -treeout result/otus.tree
    #Generate five distance matrices：bray_curtis, euclidean, jaccard, manhatten, unifrac
    usearch -beta_div result/otutab_rare.txt -tree result/otus.tree \
      -filename_prefix result/beta/


## 8. Compilation of species annotation classifications

    #OTU corresponding species annotation in a 2-column format: remove confidence values from sintax, retain species annotations, replace ":" with "_", and remove quotation marks
    cut -f 1,4 result/otus.sintax \
      |sed 's/\td/\tk/;s/:/__/g;s/,/;/g;s/"//g' \
      > result/taxonomy2.txt
    head -n3 result/taxonomy2.txt

    #OTU corresponding species in an 8-column format: please note that annotations are not uniform
    #Generate a species table where blanks in the OTU/ASV are filled with "Unassigned"
    awk 'BEGIN{OFS=FS="\t"}{delete a; a["k"]="Unassigned";a["p"]="Unassigned";a["c"]="Unassigned";a["o"]="Unassigned";a["f"]="Unassigned";a["g"]="Unassigned";a["s"]="Unassigned";\
      split($2,x,";");for(i in x){split(x[i],b,"__");a[b[1]]=b[2];} \
      print $1,a["k"],a["p"],a["c"],a["o"],a["f"],a["g"],a["s"];}' \
      result/taxonomy2.txt > temp/otus.tax
    sed 's/;/\t/g;s/.__//g;' temp/otus.tax|cut -f 1-8 | \
      sed '1 s/^/OTUID\tKingdom\tPhylum\tClass\tOrder\tFamily\tGenus\tSpecies\n/' \
      > result/taxonomy.txt
    head -n3 result/taxonomy.txt

    #Count phylum, class, order, family, and genus using the rank parameters p, c, o, f, g respectively.
    mkdir -p result/tax
    for i in p c o f g;do
      usearch -sintax_summary result/otus.sintax \
      -otutabin result/otutab_rare.txt -rank ${i} \
      -output result/tax/sum_${i}.txt
    done
    sed -i 's/(//g;s/)//g;s/\"//g;s/\#//g;s/\/Chloroplast//g' result/tax/sum_*.txt
    #List all files
    wc -l result/tax/sum_*.txt
    head -n3 result/tax/sum_g.txt


## 9. There is a reference-based quantitative feature table

    # Perform alignment against the Greengenes 97% OTUs database for use in PICRUSt/Bugbase functional prediction
    mkdir -p result/gg/

    #Method 1. Using USEARCH for alignment is faster, but encountering file size limit errors. Choose method 2
    # By default, utilize 1 core for systems with less than 10 cores, and allocate 10 cores for systems with 10 or more cores
    usearch -otutab temp/filtered.fa -otus ${db}/gg/97_otus.fa \
    	-otutabout result/gg/otutab.txt -threads 4
    # Alignment rate: 80.0%, 1 core took 11 minutes, 4 cores took 3 minutes, 10 cores took 2 minutes, and memory usage was 743MB
    head -n3 result/gg/otutab.txt

    # #Method 2. vsearch alignment is more accurate but slower. However, it exhibits stronger performance with parallelization using 24-96 threads
    # vsearch --usearch_global temp/filtered.fa --db ${db}/gg/97_otus.fa \
    #   --otutabout result/gg/otutab.txt --id 0.97 --threads 12
    # Alignment rate: 81.04%, 1 core took 30 minutes, and 12 cores took 7 minutes

    #Statistics
    usearch -otutab_stats result/gg/otutab.txt -output result/gg/otutab.stat
    cat result/gg/otutab.stat


## 10. Space cleanup and data submission.

    #Delete intermediate large files
    rm -rf temp/*.fq

    # Calculate MD5 values for paired-end data to use for data submission
    cd seq
    md5sum *_1.fq.gz > md5sum1.txt
    md5sum *_2.fq.gz > md5sum2.txt
    paste md5sum1.txt md5sum2.txt | awk '{print $2"\t"$1"\t"$4"\t"$3}' | sed 's/*//g' > ../result/md5sum.txt
    rm md5sum*
    cd ..
    cat result/md5sum.txt
    

# RAnalysis of linguistic diversity and species composition


## 1. Alpha diversity

### 1.1 Box plot of alpha diversity.

    # View Help
    Rscript ${db}/script/alpha_boxplot.R -h
    # Full parameters, diversity index optional richness chao1 ACE shannon simpson invsimpson
    Rscript ${db}/script/alpha_boxplot.R --alpha_index richness \
      --input result/alpha/vegan.txt --design result/metadata.txt \
      --group Group --output result/alpha/ \
      --width 89 --height 59
    # Use loops to plot 6 commonly used indices
    for i in `head -n1 result/alpha/vegan.txt|cut -f 2-`;do
      Rscript ${db}/script/alpha_boxplot.R --alpha_index ${i} \
        --input result/alpha/vegan.txt --design result/metadata.txt \
        --group Group --output result/alpha/ \
        --width 89 --height 59
    done
    mv alpha_boxplot_TukeyHSD.txt result/alpha/

    # Alpha diversity histogram + standard deviation
    Rscript ${db}/script/alpha_barplot.R --alpha_index richness \
      --input result/alpha/vegan.txt --design result/metadata.txt \
      --group Group --output result/alpha/ \
      --width 89 --height 59

### 1.2  Dilution curve

    Rscript ${db}/script/alpha_rare_curve.R \
      --input result/alpha/alpha_rare.txt --design result/metadata.txt \
      --group Group --output result/alpha/ \
      --width 120 --height 59

### 1.3 Diversity Venn diagram

    # Three groups of comparisons: -f input file, -a/b/c/d/g group name, -w/u for width and height inches, -p output file name suffix
    bash ${db}/script/sp_vennDiagram.sh \
      -f result/alpha/otu_group_exist.txt \
      -a WT -b KO -c OE \
      -w 3 -u 3 \
      -p WT_KO_OE
    # Four sets of comparisons, the diagram and code are shown in the input file directory, and the running directory is the root directory of the current project
    bash ${db}/script/sp_vennDiagram.sh \
      -f result/alpha/otu_group_exist.txt \
      -a WT -b KO -c OE -d All \
      -w 3 -u 3 \
      -p WT_KO_OE_All
    #  plots Venn diagrams online https://www.ehbio.com/test/venn

## 2. Beta diversity

### 2.1 Distance matrix heatmap 

    # Taking bray_curtis as an example, -f input file, -h whether the clustering is TRUE/FALSE, -u/v is width and height inches
    bash ${db}/script/sp_pheatmap.sh \
      -f result/beta/bray_curtis.txt \
      -H 'TRUE' -u 6 -v 5
    # Add grouping notes such as genotype and location for columns 2,4
    cut -f 1-2 result/metadata.txt > temp/group.txt
    # -P to add row comment files, -Q to add column comments
    bash ${db}/script/sp_pheatmap.sh \
      -f result/beta/bray_curtis.txt \
      -H 'TRUE' -u 6.9 -v 5.6 \
      -P temp/group.txt -Q temp/group.txt
    # The distance matrix is similar to correlation, try corrplot or ggcorrplot to plot more styles
    # - [Plot the correlation coefficient matrix corrplot](http://mp.weixin.qq.com/s/H4_2_vb2w_njxPziDzV4HQ)
    # - [Correlation matrix visualization ggcorrplot](http://mp.weixin.qq.com/s/AEfPqWO3S0mRnDZ_Ws9fnw)

### 2.2 Primary coordinate analysis PCoA

    # Input file, choose grouping, output file, image size in mm, and statistical results in beta_pcoa_stat.txt
    Rscript ${db}/script/beta_pcoa.R \
      --input result/beta/bray_curtis.txt --design result/metadata.txt \
      --group Group --label FALSE --width 89 --height 59 \
      --output result/beta/bray_curtis.pcoa.pdf
    # Add sample labels with the parameter --label TRUE
    Rscript ${db}/script/beta_pcoa.R \
      --input result/beta/bray_curtis.txt --design result/metadata.txt \
      --group Group --label TRUE --width 89 --height 59 \
      --output result/beta/bray_curtis.pcoa.label.pdf
    mv beta_pcoa_stat.txt result/beta/
      
### 2.3 Constrained Principal Coordinates Analysis CPCoA

    Rscript ${db}/script/beta_cpcoa.R \
      --input result/beta/bray_curtis.txt --design result/metadata.txt \
      --group Group --output result/beta/bray_curtis.cpcoa.pdf \
      --width 89 --height 59
    # Add sample labels with the parameter --label TRUE
    Rscript ${db}/script/beta_cpcoa.R \
      --input result/beta/bray_curtis.txt --design result/metadata.txt \
      --group Group --label TRUE --width 89 --height 59 \
      --output result/beta/bray_curtis.cpcoa.label.pdf
      
## 3. Species composition Taxonomy

### 3.1 Stacked bar chart Stackplot

    # Taking the phylum (p) level as an example, the results will include two files: "output.sample.pdf" and "output.group.pdf"
    Rscript ${db}/script/tax_stackplot.R \
      --input result/tax/sum_p.txt --design result/metadata.txt \
      --group Group --color ggplot --legend 7 --width 89 --height 59 \
      --output result/tax/sum_p.stackplot
    # Modify colors using ggplot manual1 (22), Paired (12), or Set3 (12) color palettes
    Rscript ${db}/script/tax_stackplot.R \
      --input result/tax/sum_p.txt --design result/metadata.txt \
      --group Group --color Paired --legend 12 --width 181 --height 119 \
      --output result/tax/sum_p.stackplotPaired
      
    # Batch draw using input that includes phylum (p), class (c), order (o), family (f), and genus (g) for a total of 5 levels
    for i in p c o f g; do
    Rscript ${db}/script/tax_stackplot.R \
      --input result/tax/sum_${i}.txt --design result/metadata.txt \
      --group Group --output result/tax/sum_${i}.stackplot \
      --legend 8 --width 89 --height 59; done

### 3.2 Chord/circular diagram using Circlize

    # Taking class (c) as an example, draw the top 5 groups.
    i=c
    Rscript ${db}/script/tax_circlize.R \
      --input result/tax/sum_${i}.txt --design result/metadata.txt \
      --group Group --legend 5
    # The results are located in the current directory: "circlize.pdf" (with random colors) and "circlize_legend.pdf" (with specified colors and legend)
    # Move and rename to match the classification level
    mv circlize.pdf result/tax/sum_${i}.circlize.pdf
    mv circlize_legend.pdf result/tax/sum_${i}.circlize_legend.pdf

### 3.3 Treemap visualization (for reference)

    # Hierarchical treemap visualization depicting species relationships. Input feature table and species annotation, output treemap
    # Specify the number of features and image dimensions. Generating the treemap for 100 ASVs took 12 seconds
    Rscript ${db}/script/tax_maptree.R \
      --input result/otutab.txt --taxonomy result/taxonomy.txt \
      --output result/tax/tax_maptree.pdf \
      --topN 100 --width 183 --height 118

# 24、Comparison of differences

## 1. R language difference analysis

### 1.1 Comparison of differences 

    mkdir -p result/compare/
    # Input feature table, metadata; Specify the grouping column name, comparison group, and abundance
    # Select the methods wilcox/t.test/edgeR, pvalue, and fdr and output directories
    compare="KO-WT"
    Rscript ${db}/script/compare.R \
      --input result/otutab.txt --design result/metadata.txt \
      --group Group --compare ${compare} --threshold 0.1 \
      --method edgeR --pvalue 0.05 --fdr 0.2 \
      --output result/compare/
    # Common error: Error in file(file, ifelse(append, "a", "w")): Unable to open link Calls: write.table -> file
    # Solution: The output directory does not exist, just create a directory

### 1.2 Volcano map

    #Enter compare. As a result of R, output volcano map with data labels, you can specify the image size
    Rscript ${db}/script/compare_volcano.R \
      --input result/compare/${compare}.txt \
      --output result/compare/${compare}.volcano.pdf \
      --width 89 --height 59

### 1.3 Heat map

    # Enter compare. R results, filter the number of columns, specify metadata and grouping, species annotations, figure size inches and font size
    bash ${db}/script/compare_heatmap.sh -i result/compare/${compare}.txt -l 7 \
       -d result/metadata.txt -A Group \
       -t result/taxonomy.txt \
       -w 8 -h 5 -s 7 \
       -o result/compare/${compare}

### 1.4 Map of manhattan

    # idifference comparison results, t species annotation, p legend, w width, v height, s size, l legend maximum
    # The legend shows no figure, you can increase the height v to 119+, and the AI puzzle is KO-WT.heatmap.emf in the later stage
    bash ${db}/script/compare_manhattan.sh -i result/compare/${compare}.txt \
       -t result/taxonomy.txt \
       -p result/tax/sum_p.txt \
       -w 183 -v 59 -s 7 -l 10 \
       -o result/compare/${compare}.manhattan.p.pdf
       
          # There are only 6 gates in the picture above, switching to class c and -L class to show the details
    bash ${db}/script/compare_manhattan.sh -i result/compare/${compare}.txt \
       -t result/taxonomy.txt \
       -p result/tax/sum_c.txt \
       -w 183 -v 59 -s 7 -l 10 -L Class \
       -o result/compare/${compare}.manhattan.c.pdf
    # Show the full legend and use AI puzzles
    bash ${db}/script/compare_manhattan.sh -i result/compare/${compare}.txt \
       -t result/taxonomy.txt \
       -p result/tax/sum_c.txt \
       -w 183 -v 149 -s 7 -l 10 -L Class \
       -o result/compare/${compare}.manhattan.c.legend.pdf

### 1.5Drawing of individual features

    # Differential ASV was screened and displayed in descending order by abundance of KO group, and the top 10 were displayed by ID
    awk '$4<0.05' result/compare/KO-WT.txt | sort -k7,7nr | cut -f1 | head
    # Differential OTU details display
    Rscript ${db}/script/alpha_boxplot.R --alpha_index ASV_2 \
      --input result/otutab.txt --design result/metadata.txt \
      --transpose TRUE --scale TRUE \
      --width 89 --height 59 \
      --group Group --output result/compare/feature_ 
    # If the ID does not exist, an error will be reported： Error in data.frame(..., check.names = FALSE) : Parameter values mean different number of rows: 0, 18  Calls: alpha_boxplot -> cbind -> cbind -> data.frame
    # Specify a column sort: All descending by mean genus abundance
    csvtk -t sort -k All:nr result/tax/sum_g.txt | head
    # Poor details are shown
    Rscript ${db}/script/alpha_boxplot.R --alpha_index Lysobacter \
      --input result/tax/sum_g.txt --design result/metadata.txt \
      --transpose TRUE \
      --width 89 --height 59 \
      --group Group --output result/compare/feature_

### 1.5 Ternary Diagram

  #For reference examples, see: result\compare\ternary\ternary. RMD documentation
  #Alternative Tutorial[246.Application and Drawing Practice of Ternary Diagrams](https://mp.weixin.qq.com/s/3w3ncpwjQaMRtmIOtr2Jvw)

## 2. STAMP input file preparation

### 2.1 Generate an input file

    Rscript ${db}/script/format2stamp.R -h
    mkdir -p result/stamp
    Rscript ${db}/script/format2stamp.R --input result/otutab.txt \
      --taxonomy result/taxonomy.txt --threshold 0.01 \
      --output result/stamp/tax
    # Optional RMD documentation see result/format2stamp. Rmd

### 2.2 Plot extended column charts and tables

    compare="KO-WT"
    # Replace ASV (result/otutab.txt) with genus (result/tax/sum_g.txt)
    Rscript ${db}/script/compare_stamp.R \
      --input result/stamp/tax_5Family.txt --metadata result/metadata.txt \
      --group Group --compare ${compare} --threshold 0.1 \
      --method "t.test" --pvalue 0.05 --fdr "none" \
      --width 189 --height 159 \
      --output result/stamp/${compare}
    # Optional RMD documentation can be found in result/compareStamp.Rmd

## 3. LEfSe input file preparation

    ### 3.1. Command line make-file
    # Optional command line to generate input files
    Rscript ${db}/script/format2lefse.R -h
    mkdir -p result/lefse
    # Threshold controls abundance filtering to control the number of branches in the map
    Rscript ${db}/script/format2lefse.R --input result/otutab.txt \
      --taxonomy result/taxonomy.txt --design result/metadata.txt \
      --group Group --threshold 0.4 \
      --output result/lefse/LEfSe

    ### 3.2  RMD generation input file (optional)
    #1. The result directory contains three files: otutab.txt, metadata.txt, taxonomy .txt；
    #2. Rstudio opens format2lefse.rmd in EasyAmplicon, saves to the result directory and Knit generates input files and repeatable web pages；
    
	### 3.3 LEfSe Analysis 
    #Method 1. Open the LEfSe .txt and submit online https://www.bic.ac.cn/BIC/#/analysis?page=b%27MzY%3D%27
    #Method 2. LEfSe local analysis (Linux system only, optional learning), see the appendix for the reference code
    #Method 3. The LEfSe website is available online


# 25、QIIME 2 Analyze the process

    # See the code for details qiime2/pipeline_qiime2.sh


# 31、Feature prediction

## 1. PICRUSt Feature prediction

    # PICRUSt 1.0
    # Method 1. Utilize http://www.ehbio.com/ImageGP for online analysis of gg/otutab.txt
    # Method 2. Linux server users can refer to "Appendix 2. PICRUSt Capability Prediction "enables software installation and analysis
    # The results were then compared using STAMP/R

    # R language drawing
    # Input file format adjustments
    l=L2
    sed '/# Const/d;s/OTU //' result/picrust/all_level.ko.${l}.txt > result/picrust/${l}.txt
    num=`head -n1 result/picrust/${l}.txt|wc -w`
    paste <(cut -f $num result/picrust/${l}.txt) <(cut -f 1-$[num-1] result/picrust/${l}.txt) \
      > result/picrust/${l}.spf
    cut -f 2- result/picrust/${l}.spf > result/picrust/${l}.mat.txt
    awk 'BEGIN{FS=OFS="\t"} {print $2,$1}' result/picrust/${l}.spf | sed 's/;/\t/' | sed '1 s/ID/Pathway\tCategory/' \
      > result/picrust/${l}.anno.txt
    # Comparison of differences
    compare="KO-WT"
    Rscript ${db}/script/compare.R \
      --input result/picrust/${l}.mat.txt --design result/metadata.txt \
      --group Group --compare ${compare} --threshold 0 \
      --method wilcox --pvalue 0.05 --fdr 0.2 \
      --output result/picrust/
    # The result ${compare} can be filtered .txt
    # Generate a bar plot(histogram) for the specified groups (A/B), color-coded by higher taxonomic levels and faceted.
    Rscript ${db}/script/compare_hierarchy_facet.R \
      --input result/picrust/${compare}.txt \
      --data MeanA \
      --annotation result/picrust/${l}.anno.txt \
      --output result/picrust/${compare}.MeanA.bar.pdf
    # Generate bar plots depicting significant differences between two groups, with faceting based on higher taxonomic levels
    Rscript ${db}/script/compare_hierarchy_facet2.R \
      --input result/picrust/${compare}.txt \
      --pvalue 0.05 --fdr 0.1 \
      --annotation result/picrust/${l}.anno.txt \
      --output result/picrust/${compare}.bar.pdf
      
    # PICRUSt 2.0
    # Software installation, Appendix 6: Exporting and Importing the PICRUSt Environment
    # Utilization, Appendix 7: Functional Prediction with PICRUSt2

## 2. Looping through elements of FAPROTAX

    ## Method1. Online analysis, recommended for use http://www.bic.ac.cn/ImageGP/index.php/Home/Index/FAPROTAX.html online analysis 
    ## Analysis on Linux, such as within the QIIME 2 environment, is explained in detail in Appendix 3

## 3. Bugbase Bacterial phenotype prediction

    ### 1. Bugbase command-line analysis
    cd ${wd}/result
    bugbase=${db}/script/BugBase
    rm -rf bugbase/
    # The script has been optimized for R 4.0, and the "biom" package has been updated to "biomformat"
    Rscript ${bugbase}/bin/run.bugbase.r -L ${bugbase} \
      -i gg/otutab.txt -m metadata.txt -c Group -o bugbase/

    ### 2. Other available analyses
    # Using http://www.bic.ac.cn/ImageGP/index.php/Home/Index/BugBase.html
    # Official website，https://bugbase.cs.umn.edu/ ，There are errors, not recommended
    # Bacterial phenotype prediction using Bugbase on Linux, please refer to Appendix 4: Bugbase Bacterial Phenotype Prediction


# 32、MachineLearning

    # The R code for using the RandomForest package can be found in the "advanced/RandomForestClassification and RandomForestRegression section
    ## The code for using Random Forest and Adaboost with Silme2 can be found in the "EasyMicrobiome/script/slime2" directory, specifically in the "slime2.py" file. Please refer to Appendix 5 for more details.

# 33、Evolutionary tree

    cd ${wd}
    mkdir -p result/tree
    cd ${wd}/result/tree

## 1. Filtering for high abundance/specified features

    #Method1. Filter features based on abundance, typically selecting 0.001 or 0.005, with the OTU count falling within the range of 30 to 150
    #Count the number of ASVs in the feature table, for example, a total of 1609
    tail -n+2 ../otutab_rare.txt | wc -l
    #Filter high-abundance OTUs based on a relative abundance threshold of 0.2%
    usearch -otutab_trim ../otutab_rare.txt \
        -min_otu_freq 0.002 \
        -output otutab.txt
    #Count the number of features in the filtered OTU table; there are a total of approximately 81 features
    tail -n+2 otutab.txt | wc -l

    #Method2. Filter based on quantity/count
    # #Sort by abundance, defaulting to descending order (from largest to smallest)
    # usearch -otutab_sortotus ../otutab_rare.txt  \
    #     -output otutab_sort.txt
    # #Extract the OTU IDs from the top specified high-abundance OTUs, such as the top 100,
    # sed '1 s/#OTU ID/OTUID/' otutab_sort.txt \
    #     | head -n101 > otutab.txt

    #Modify the column name of the feature ID
    sed -i '1 s/#OTU ID/OTUID/' otutab.txt
    #Extract IDs for sequence retrieval
    cut -f 1 otutab.txt > otutab_high.id

    # Filter high-abundance bacteria/specify differential bacteria-associated OTU sequences
    usearch -fastx_getseqs ../otus.fa -labels otutab_high.id \
        -fastaout otus.fa
    head -n 2 otus.fa

    ## Filter OTUs based on species annotation
    awk 'NR==FNR{a[$1]=$0} NR>FNR{print a[$1]}' ../taxonomy.txt \
        otutab_high.id > otutab_high.tax

    #Obtain group mean values corresponding to OTUs for use in generating a sample heatmap
    #依赖Dependent on the previously calculated group mean values using the "otu_mean.R" script
    awk 'NR==FNR{a[$1]=$0} NR>FNR{print a[$1]}' ../otutab_mean.txt otutab_high.id \
        | sed 's/#OTU ID/OTUID/' > otutab_high.mean
    head -n3 otutab_high.mean

    #Merge species annotation and abundance into an annotation file
    cut -f 2- otutab_high.mean > temp
    paste otutab_high.tax temp > annotation.txt
    head -n 3 annotation.txt

## 2. Build an evolutionary tree

    # The starting files are the otus.fa (sequence), annotation.txt (species and relative abundance) files in the result/tree directory
    # Muscle software for sequence alignment, 3s
    muscle -in otus.fa -out otus_aligned.fas

    ### Method1. Use IQ-TREE to quickly build ML evolutionary tree, 2m
    rm -rf iqtree
    mkdir -p iqtree
    iqtree -s otus_aligned.fas \
        -bb 1000 -redo -alrt 1000 -nt AUTO \
        -pre iqtree/otus

    ### Method2. FastTree(Linux)
    # Note that the FastTree software input file is a fasta format file, not the usual Phylip format. The output file is in Newick format。
    # This method is suitable for big data, such as phylogenetic trees of hundreds of OTUs！
    # To install fasttree on Ubuntu, you can use `apt install fasttree`
    # fasttree -gtr -nt otus_aligned.fas > otus.nwk

## 3. Evolutionary Tree Beautification

    # Visit http://itol.embl.de/，upload otus.nwk，and then drag and drop the annotation scheme generated below to be beautiful on the tree
    ## Scheme 1. Outer circle color, shape classification and abundance scheme
    # annotation.txt OTU corresponds to species annotation and abundance，
    # -a If the input column is not found, the operation will be terminated (default is not executed) -c Convert the integer column to a factor or a number with a decimal point, -t Convert the ID column when it deviates from the prompt label, -w Color band, area width, etc., -D output directory, -i OTU column name, -l OTU display name such as species/genus/family name,
    # cd ${wd}/result/tree
    Rscript ${db}/script/table2itol.R -a -c double -D plan1 -i OTUID -l Genus -t %s -w 0.5 annotation.txt
    # Each column in the generated comment file is a separate file

    ## Scheme 2. Generate abundance histogram annotation file
    Rscript ${db}/script/table2itol.R -a -d -c none -D plan2 -b Phylum -i OTUID -l Genus -t %s -w 0.5 annotation.txt

    ## Scheme 3. Generate heatmap annotation files
    Rscript ${db}/script/table2itol.R -c keep -D plan3 -i OTUID -t %s otutab.txt

    ## Scheme 4. Convert integers into factors to generate annotation files
    Rscript ${db}/script/table2itol.R -a -c factor -D plan4 -i OTUID -l Genus -t %s -w 0 annotation.txt

    # The tree iqtree/otus.contree is displayed on http://itol.embl.de/, drag and drop files in different Plans to add tree comments

    # return to working directory
    cd ${wd}

## 4. Evolutionary tree visualization

   https://www.bic.ac.cn/BIC/#/ Provides an easier way to visualize

# Additional video

    # Catalog Supp, online courses have corresponding videos (the numbers may be different, look for keywords)


## S1. Network Analysis R/CytoGephi

    # Directory Supp/S1NetWork

## S2. Traceability and Markov chains

    # Directory Supp/S2SourcetrackerFeastMarkov


## S11. Network analysis ggClusterNet

    # Code: advanced/ggClusterNet/Practice.Rmd

## S12, Microeco package data visualization

    # Code: advanced/microeco/Practice.Rmd


# Appendix: Analysis under Linux Server (Optional)

    #Note: The following code may not be able to run under Windows. It is recommended to install related programs on Linux or conda under the Linux subsystem under Windows.

## 1. LEfSe analysis

    mkdir -p ~/amplicon/lefse
    cd ~/amplicon/lefse
    # format2lefse.Rmd code production or upload input file LEfSe.txt
    # install lefse
    # conda install lefse

    #Format conversion to lefse internal format
    lefse-format_input.py LEfSe.txt input.in -c 1 -o 1000000
    #run lefse
    run_lefse.py input.in input.res
    #Draw a species tree to annotate differences
    lefse-plot_cladogram.py input.res cladogram.pdf --format pdf
    #Draw a histogram of all differential features
    lefse-plot_res.py input.res res.pdf --format pdf
    #Draw a single features histogram (same as barplot in STAMP)
    head input.res #View list of difference features
    lefse-plot_features.py -f one --feature_name "Bacteria.Firmicutes.Bacilli.Bacillales.Planococcaceae.Paenisporosarcina" \
       --format pdf input.in input.res Bacilli.pdf
    #Draw all difference feature histograms in batches, use with caution (it is also difficult to read hundreds of difference result histograms)
    mkdir -p features
    lefse-plot_features.py -f diff --archive none --format pdf \
      input.in input.res features/


## 2. PICRUSt function prediction

    #It is recommended to use http://www.bic.ac.cn/BIC/#/analysis?tool_type=tool&page=b%27Mzk%3D%27 online analysis  
    #Users with Linux servers can refer to the following code to build a local process

    n=picrust
    conda create -n ${n} ${n} -c bioconda -y

    wd=/mnt/d/amplicon
    cd $wd/result/gg
    # startup environment
    conda activate picrust
    #Upload gg/otutab.txt to the current directory
    #Convert to the common format of OTU table, which is convenient for downstream analysis and statistics
    biom convert -i otutab.txt \
        -o otutab.biom \
        --table-type="OTU table" --to-json

    # Set the database directory, such as /mnt/d/db
    db=~/db
    #Corrected copy number, 30s, 102M
    normalize_by_copy_number.py -i otutab.biom \
        -o otutab_norm.biom \
        -c ${db}/picrust/16S_13_5_precalculated.tab.gz
    #Predict metagenomic KO table, 3m, 1.5G, biom for downstream classification, txt for viewing and analysis
    predict_metagenomes.py -i otutab_norm.biom \
        -o ko.biom \
        -c ${db}/picrust/ko_13_5_precalculated.tab.gz
    predict_metagenomes.py -f -i otutab_norm.biom \
        -o ko.txt \
        -c ${db}/picrust/ko_13_5_precalculated.tab.gz

    #Classify and summarize by function level, -c output KEGG_Pathways, divided into 1-3 levels
    sed  -i '/# Constru/d;s/#OTU //' ko.txt
    num=`head -n1 ko.txt|wc -w`
    paste <(cut -f $num ko.txt) <(cut -f 1-$[num-1] ko.txt) > ko.spf
    for i in 1 2 3;do
      categorize_by_function.py -f -i ko.biom -c KEGG_Pathways -l ${i} -o pathway${i}.txt
      sed  -i '/# Const/d;s/#OTU //' pathway${i}.txt
      paste <(cut -f $num pathway${i}.txt) <(cut -f 1-$[num-1] pathway${i}.txt) > pathway${i}.spf
    done
    wc -l *.spf


## 3. FAPROTAXS element cycle

    # set working directory
    wd=/mnt/d/amplicon/result/faprotax/
    mkdir -p ${wd} && cd ${wd}
    # Set script directory
    sd=/mnt/d/EasyMicrobiome/script/FAPROTAX_1.2.6

    ### 1. Software Installation
    # Note: The software has been downloaded to the EasyAmplicon/script directory, and running in the qiime2 environment can satisfy the dependencies
    #(Optional) Download the new version of the software, take version 1.2.6 as an example, update the database on July 14, 2022
    #wget -c https://pages.uoregon.edu/slouca/LoucaLab/archive/FAPROTAX/SECTION_Download/MODULE_Downloads/CLASS_Latest%20release/UNIT_FAPROTAX_1.2.6/FAPROTAX_1.2.6.zip
    #unzip
    #unzip FAPROTAX_1.2.6.zip
    #(Optional) dependencies, you can use conda to install dependency packages
    #conda install numpy
    #conda install biom
    # View conda environment name and location
    # conda env list
    #Create a new python3 environment and configure dependencies, or enter the qiime2 python3 environment
    conda activate qiime2-2023.2
    # source /home/silico_biotech/miniconda3/envs/qiime2/bin/activate
    #Test whether it can be run, and the help will pop up and it will work normally
    python $sd/collapse_table.py

    ### 2. Make input OTU table
    #txt to biom json format
    biom convert -i ../otutab_rare.txt -o otutab_rare.biom --table-type="OTU table" --to-json
    #Add species annotation
    biom add-metadata -i otutab_rare.biom --observation-metadata-fp ../taxonomy2.txt \
      -o otutab_rare_tax.biom --sc-separated taxonomy \
      --observation-header OTUID,taxonomy
    #Specify input file, species annotation, output file, annotation column names, attribute column names

    ### 3. FAPROTAX Functional Prediction
    #python runs the collapse_table.py script, enters the OTU table tax.biom with species annotations,
    #-g specify database location, species annotation column names, output process information, force overwrite results, result files and details
    #Download faprotax.txt for statistical analysis with experimental design
    #faprotax_report.txt View which OTUs are from specific sources in each category
    python ${sd}/collapse_table.py -i otutab_rare_tax.biom \
      -g ${sd}/FAPROTAX.txt \
      --collapse_by_metadata 'taxonomy' -v --force \
      -o faprotax.txt -r faprotax_report.txt

    ### 4. Create matrix with or without functional annotations corresponding to OTUs
    # Filter the ASV (OTU) comment line and the title of the previous line
    grep 'ASV_' -B 1 faprotax_report.txt | grep -v -P '^--$' > faprotax_report.clean
    # The faprotax_report_sum.pl script organizes the data into tables and is located in public/scrit
    perl ${sd}/../faprotax_report_sum.pl -i faprotax_report.clean -o faprotax_report
    # 查看功能有无矩阵，-S不换行
    less -S faprotax_report.mat


## 4. Bugbase bacterial phenotype prediction

    ### 1.Software installation (it has been integrated into EasyMicrobiome, the original code needs to be updated to run now)
    #There are two options, the first is recommended, the second is optional, and only needs to be run once
    # #Method 1. git download, need to have git
    # git clone https://github.com/knights-lab/BugBase
    # #Method 2. Download and extract
    # wget -c https://github.com/knights-lab/BugBase/archive/master.zip
    # mv master.zip BugBase.zip
    # unzip BugBase.zip
    # mv BugBase-master/ BugBase

    cd BugBase
    #Install dependencies
    export BUGBASE_PATH=`pwd`
    export PATH=$PATH:`pwd`/bin
    #All dependencies installed
    run.bugbase.r -h
    #Test Data
    run.bugbase.r -i doc/data/HMP_s15.txt -m doc/data/HMP_map.txt -c HMPBODYSUBSITE -o output


    ### 2. Prepare input file
    cd ~/amplicon/result
    #Input file: biom format based on greengene OTU table (local analysis supports txt format without conversion) and mapping file (add # to the first line of metadata.txt)
    #Upload experimental design + otutab_gg.txt just generated
    #Generate biom1.0 format for online analysis
    biom convert -i gg/otutab.txt -o otutab_gg.biom --table-type="OTU table" --to-json
    sed '1 s/^/#/' metadata.txt > MappingFile.txt
    #Download otutab_gg.biom and MappingFile.txt for online analysis

    ### 3. Local Analysis 
    export BUGBASE_PATH=`pwd`
    export PATH=$PATH:`pwd`/bin
    run.bugbase.r -i otutab_gg.txt -m MappingFile.txt -c Group -o phenotype/

## 5. Silme2 Random Forest/Adaboost

    #Download and install
    # cd ~/software/
    # wget https://github.com/swo/slime2/archive/master.zip
    # mv master.zip slime2.zip
    # unzip slime2.zip
    # mv slime2-master/ slime2
    # cp slime2/slime2.py ~/bin/
    # chmod +x ~/bin/slime2.py
    #Install dependencies
    # sudo pip3 install --upgrade pip
    # sudo pip3 install pandas
    # sudo pip3 install sklearn

    # Use actual combat (use the Python3 environment of QIIME 2, take Windows as an example)
    conda activate qiime2-2023.2
    cd /mnt/d/EasyMicrobiome/script/slime2
    #Use adaboost to calculate 10,000 times (16.7s), recommend tens of millions of times
    ./slime2.py otutab.txt design.txt --normalize --tag ab_e4 ab -n 10000
    #Use RandomForest to calculate 10,000 times (14.5s), recommend millions of times, support multi-threading
    ./slime2.py otutab.txt design.txt --normalize --tag rf_e4 rf -n 10000


## 6. PICRUSt2 environment export and import
    
    # Method 1. Direct installation
    n=picrust2
    conda create -n ${n} -c bioconda -c conda-forge ${n}=2.3.0_b
    # load environment
    conda activate ${n}

    # Method 2. Export the installation environment
    cd ~/db/conda/
    # Set the environment name
    n=picrust2
    conda activate ${n}
    # The packaging environment is a compressed package
    conda pack -n ${n} -o ${n}.tar.gz
    # Export a list of software installations
    conda env export > ${n}.yml
    # Add permissions for easy download and use by others
    chmod 755 ${n}.*
    
    # Method 3. Import the installation environment, such as qiime2 humann2 meta (including picurst)
    n=picrust2
    # Copy the installation package, or download my environment package
    wget -c ftp://download.nmdc.cn/tools/conda/${n}.tar.gz
    # Specify the installation directory and unzip it
    condapath=~/miniconda3
    mkdir -p ${condapath}/envs/${n}
    tar -xvzf ${n}.tar.gz -C ${condapath}/envs/${n}
    # Activate the environment and initialize
    source ${condapath}/envs/${n}/bin/activate
    conda unpack

## 7. PICRUSt2 function prediction

    # (Optional) PICRUSt2 (Linux subsystem under Linux/Windows, requires >16GB memory)
    # For installation, refer to Appendix 6 to directly download the installation package and decompress it to use.
    
    # load environment
    conda activate picrust2
    # Enter the working directory, the server needs to modify the working directory
    wd=/mnt/d/amplicon/result/picrust2
    mkdir -p ${wd} && cd ${wd}
    # Run the process, the memory is 15.7GB, and it takes 12m
    picrust2_pipeline.py -s ../otus.fa -i ../otutab.txt -o ./out -p 8
    # Add EC/KO/Pathway notes
    cd out
    add_descriptions.py -i pathways_out/path_abun_unstrat.tsv.gz -m METACYC \
      -o pathways_out/path_abun_unstrat_descrip.tsv.gz
    add_descriptions.py -i EC_metagenome_out/pred_metagenome_unstrat.tsv.gz -m EC \
      -o EC_metagenome_out/pred_metagenome_unstrat_descrip.tsv.gz
    add_descriptions.py -i KO_metagenome_out/pred_metagenome_unstrat.tsv.gz -m KO \
      -o KO_metagenome_out/pred_metagenome_unstrat_descrip.tsv.gz 
    # KEGG Merged by Hierarchy
    db=/mnt/d/EasyMicrobiome/
    zcat KO_metagenome_out/pred_metagenome_unstrat.tsv.gz > KEGG.KO.txt
    python3 ${db}/script/summarizeAbundance.py \
      -i KEGG.KO.txt \
	    -m ${db}/kegg/KO1-4.txt \
	    -c 2,3,4 -s ',+,+,' -n raw \
	    -o KEGG
    # Count the number of features at each level
    wc -l KEGG*



# Common problem

## 1. File phred quality error - Fastq quality value 64 to 33

    # Use head to view fastq files. Most of the phred64 quality values ​​​​are lowercase letters. You need to use the --fastq_convert command of vsearch to convert them to the general phred33 format.

    cd /c/amplicon/FAQ/01Q64Q33
    # Preview the phred64 format, pay attention to see that the quality value in line 4 is mostly lowercase letters
    head -n4 test_64.fq
    # Convert quality value 64 to encoding format 33
    vsearch --fastq_convert test_64.fq \
        --fastq_ascii 64 --fastq_asciiout 33 \
        --fastqout test.fq 
    # View the converted 33 encoding format, the quality value is mostly uppercase letters
    head -n4 test.fq

    # If it is an Ion torrent sequencing result, since it is a non-mainstream sequencing platform, it needs the help of the company to convert it into a standard Phred33 format file before it can be used.

## 2. The paired-end sequences have been merged, and now you are adding sample names to the single-end sequences

    # For amplicon analysis, the sequence names need to be in the format of sample name + sequence number. For paired-end sequences, you can directly add the sample name during merging. However, for single-end sequences or merged paired-end sequences, you need to add the sample name separately.Here, you can use the --relabel parameter in vsearch's --fastq_convert command to add the sample name. 
    cd /c/amplicon/FAQ/02relabel
    # To view the sequence names in a file
    head -n1 test.fq
    # Rename sequences according to samples
    vsearch --fastq_convert test.fq \
        --relabel WT1. \
        --fastqout WT1.fq
    # View the results of renaming
    head -n1 WT1.fq

## 3. The data is too large to use USEARCH for clustering or denoising. Replace it with vsearch

    # When restricted to the limitations of the free version of USEARCH, you can reduce the non-redundant data size by increasing the "minuniquesize" parameter. If you're dealing with over 10,000 OTUs/ASVs and the downstream analysis is taking too long, ensure that the OTU/ASV data is under 5000. Typically, this won't be restricted, and it's also advantageous for faster downstream analyses
    # Using vsearch for clustering to generate OTUs is an option, but it lacks automatic de novo chimera removal. With an input of 2155 sequences, the clustering process resulted in 661 output clusters.

    cd /c/amplicon/FAQ/03feature
    # Rename using "relabel" and cluster at a similarity threshold of 97% without masking qmask
    # Record the input size ("sizein") and the output frequency ("sizeout")
    vsearch --cluster_size uniques.fa  \
     --relabel OTU_ --id 0.97 \
     --qmask none --sizein --sizeout \
     --centroids otus_raw.fa 

    # Then, perform de novo chimera removal. Out of the 55 chimeras, 606 are non-chimeric. Removing all OTU_1 sequences seems reasonable since Usearch does not have built-in chimera removal methods。
    # Self-alignment to remove chimeras
    vsearch --uchime_denovo otus_raw.fa \
        --nonchimeras otus.fa
    # Removing sequence frequencies
    sed -i 's/;.*//' otus.fa

## 4. Normalization of read counts to relative abundance

    cd /c/amplicon/FAQ/04norm
    # Calculating the abundance frequency of each OTU in the samples (normalized to a total sum of 1).
    usearch -otutab_counts2freqs otutab.txt \
        -output otutab_freq.txt

## 5. Running R prompts Permission denied
 
    #  For example, when using write.table to save a table, the error message might be as follows: the meaning is that there is no permission to write to the file, usually because the target file is currently open. Please close the relevant files and try again.

    Error in file(file, ifelse(append, "a", "w")) :
    Calls: write.table -> file
    : Warning message:
    In file(file, ifelse(append, "a", "w")) :
      'result/raw/otutab_nonBac.txt': Permission denied

## 6. Batch naming of files

    # If you have files A1 and A2, you can create a metadata table "metadata.txt" that maps sample names to target names. To check if sample names are unique and rename them in batch using awk.

    cd /c/amplicon/FAQ/06rename
    # (Optional) Quickly generate a file list for editing "metadata.txt" by renaming files, for example, changing "A1.fq" to "WT1.fastq" and so on. You can refer to "metadata.bak.txt" for guidance.
    ls *.fq > metadata.txt
    # Edit the list, where the second name will be the final naming. Ensure that the names are unique
    # Convert line endings
    sed -i 's/\r//' metadata.txt
    # Check if the manually named column 2 is unique
    cut -f 2 metadata.txt|wc -l
    cut -f 2 metadata.txt|sort|uniq|wc -l
    # If the two results are consistent, then the naming is non-redundant
    # You can choose to use mv (move), cp (copy), ln (hard link), or ln -s (symbolic link)
    # Here, we are using cp (copy)
    awk '{system("cp "$1" "$2)}' metadata.txt

## 7. In RStudio, the Terminal cannot find Linux commands

    # Need to add the directory "C:\Program Files\Git\usr\bin" to the system environment variables.
    # File Explorer - This PC - Properties - Advanced System Settings - Environment Variables - System Variables - Path - Edit - Add New - Enter "C:\Program Files\Git\usr\bin" - OK - OK - OK
    # Note that in Windows 10, each directory should be on a separate line. In Windows 7, multiple directories should be separated by semicolons. Make sure to append the directory towards the end.

## 8. In the "usearch -alpha_div_rare" results, "-" appears in the first two lines.

    #Issue: Fill with "-" when sampling is 0, and missing tab

    #Solution: Replace "-" with "tab character\t+0" to restore.

    cd /c/amplicon/FAQ/08rare
    sed "s/-/\t0.0/g" alpha_rare_wrong.txt\
        > alpha_rare.txt

## 9. The species annotation in otus.sintax is all "-" and requires reverse complementation of the sequences

    #The original sequence direction is incorrect, and the sequences in filtered.fa need to be reverse complemented. Restart the analysis from the scratch

    cd /c/amplicon/FAQ/09revcom
    vsearch --fastx_revcomp filtered_RC.fa \
      --fastaout filtered.fa

## 10. View and delete windows line breaks

    #In Windows, line breaks are represented as a combination of newline ($)+^M, which is equivalent to a newline in Linux plus a carriage return in macOS. When analyzing data, the Linux format is often considered the standard. Therefore, if you create a table using Excel in Windows and save it as a text file (*.txt) with tab separation, you might encounter invisible ^M symbols at the end of each line. These symbols can lead to analysis errors.You can use the cat -A command to view these symbols, and you can utilize sed to remove them.

    cd /c/amplicon/FAQ/10^M
    # Check if there is ^M at the end of the line
  	cat -A metadata.txt
  	# remove ^M, and write to new file
  	sed 's/\r//' metadata.txt > metadata.mod.txt
  	# Check for success
  	cat -A metadata.mod.txt
  	
  	# Delete the original file directly
  	sed -i 's/\r//' metadata.txt

## 11. Error in UNITE database analysis

    #USEARCH uses the utax database downloaded by UNITE, prompting various errors
	
    cd /c/amplicon/FAQ/11unite
    # Unzip Unite's useeach to use the species annotation library
    gunzip -c utax_reference_dataset_all_04.02.2020.fasta.gz > unite.fa
    # Annotate the ITS sequence, the default threshold is 0.8
    usearch --sintax  otus.fa \
      --db unite.fa \
      --tabbedout otus.sintax --strand plus
       --sintax_cutoff 0.6

    #The error message is as follows:
    ---Fatal error---
    Missing x: in name >JN874928|SH1144646.08FU;tax=d:Metazoa,p:Cnidaria,c:Hydrozoa,o:Trachylina,f:,g:Craspedacusta,s:Craspedacusta_sowerbii_SH1144646.08FU;
    “Unprintable ASCII character no 195 on or right before line 236492”
    
    # The reason for the analysis is that there are gap at the classification level. It can be solved by sed completion
    # There is a gap in the category level, sed completes it
    sed -i 's/,;/,Unnamed;/;s/:,/:Unnamed,/g' unite.fa
    # Then run the previous use search --sintax command
    #Note: There is a problem with vsearch, it is recommended to use usearch, add --strand plus at the end to run successfully

## 12. Install qiime2 locally on the Linux subsystem of Windows

    # See qiime2/pipeline_qiime2.sh for details
    n=qiime2-2023.2
    # Installation package download link 
    wget -c ftp://download.nmdc.cn/tools/conda/${n}.tar.gz
    # New environment installation
    mkdir -p ~/miniconda3/envs/${n}
    tar -xzf ${n}.tar.gz -C ~/miniconda3/envs/${n}
    # Activate and initialize the environment
    conda activate ${n}
    conda unpack

## 13. RDP 16-18 Annotated Results Comparison

    # Number of gates in a stat sequence, reduced from 60 to 39
    grep '>' ${db}/usearch/rdp_16s_v16_sp.fa|cut -f2 -d ';'|cut -f1-2 -d ','|sort|uniq|wc -l
    grep '>' ${db}/usearch/rdp_16s_v18.fa|cut -f2 -d ';'|cut -f1-2 -d ','|sort|uniq|wc -l
    # Count the number of genera in the sequence, increased from 2517 to 3061
    grep '>' ${db}/usearch/rdp_16s_v16_sp.fa|cut -f2 -d ';'|cut -f1-6 -d ','|sort|uniq|wc -l
    grep '>' ${db}/usearch/rdp_16s_v18.fa|cut -f2 -d ';'|cut -f1-6 -d ','|sort|uniq|wc -l

    cd /c/amplicon/FAQ/13rdp16_18
    # The number of doors was reduced from 15 to 13
    tail -n+2 rdp16_sintax.txt|cut -f3|sort|uniq -c|wc -l
    tail -n+2 rdp18_sintax.txt|cut -f3|sort|uniq -c|wc -l
    # Genus reduced from 176 to 144
    tail -n+2 rdp16_sintax.txt|cut -f7|sort|uniq -c|wc -l
    tail -n+2 rdp18_sintax.txt|cut -f7|sort|uniq -c|wc -l  

# Version Update History

- 2021/4/3 EasyAmplicon 1.11:
    - Upgraded R package 'amplicon' to version 1.11.0, addressing the issue with two-column metadata errors。
    - Adjusted course sequence: 2 sessions from 9:00 AM to 12:00 PM and 3 sessions from 1:30 PM to 6:00 PM every day。
    - Please provide the supplementary course materials in the 'Supp' directory。
- 2021/7/23 EasyAmplicon 1.12:
    - Upgraded R runtime environment to 4.1.0, accompanied by the complete set of packages in 4.1.zip
    - Upgraded R package 'amplicon' to version 1.12.0, removed Y-axis index from alpha_boxplot
    - alpha_boxplot.R has been updated to include parameters for normalization, transposition, and more. It can now be used to create boxplots for any feature
    - beta_pcoa/cpcoa.R has been updated to include parameters for controlling ellipses, label display, and more
    - tax_stackplot.R has been updated to include multiple color scheme options
    - The PICURST workflow has been updated, and a packaged Conda download is provided
    - Picrust2 has been updated to include the merging of KOs into KEGG pathways at levels 1-3, and a packaged Conda download is provided
    - Random Forest: Providing code for taxonomic level filtering, random number filtering, and visualization
- 2021/10/15 EasyAmplicon 1.13:
    - Upgraded R runtime environment to 4.1.1, accompanied by the latest complete set of packages in 4.1.zip
    - Metadata Variance Decomposition PERMANOVA: Added adonis calculation of variable contribution to community composition and significance analysis in Beta diversity analysis in Diversity-tutorial.Rmd
    - After upgrading the treemap, colors are no longer present. The code has been modified for reference, and the corresponding section has been removed from Diversity_tutorial.Rmd
    - In alpha_boxplot, output no longer defaults to a specific directory. You can now specify a filename prefix. Additionally, an error annotation has been added for cases with no ID
- 2022/1/7 EasyAmplicon 1.14:
    - Upgraded R runtime environment to 4.1.2, accompanied by the latest complete set of packages in 4.1.zip
    - RStudio has been updated to version 2021.09.1
    - Wentao has rewritten the 'tax_maptree' function in the amplicon package, making it independent of other packages and addressing the coloring issue
    - Added the 'compare_stamp.R' script to EasyMicrobiome, allowing direct differential comparison and visualization of extended bar plots using STAMP. For detailed code, please refer to 'result/CompareStamp.Rmd
    - Added 'compare_hierarchy_facet.R' and 'compare_hierarchy_facet2.R' scripts to EasyMicrobiome, showcasing an overview and differential analysis of KEGG pathways at levels 1 and 2
    - Updated the 'advanced' analysis directory to include environmental factors, Markov chain Monte Carlo, network modules, network comparison, random forest classification, random forest regression, and microbiota analysis
- 2023/2/3 EasyAmplicon 1.18:
    - Upgraded R runtime environment to 4.2.2, accompanied by the latest complete set of packages in 4.2.zip
    - RStudio has been updated to version 2022.12.0
    - amplicon," "EasyAmplicon," and "EasyMicrobiome" have been updated to version 1.18
    - QIIME 2 has been updated to version v2023.2
    - vsearch has been updated to version v2.22.1
    - Added 'ggClusterNet' course - by Wentao
    

Quarterly Video Course Schedule：http://www.ehbio.com/trainLongTerm/TrainLongTerm/amplicongenomeLearnGuide.html

When using this script, please reference the following passage：

If used this script, please cited:

**Yong-Xin Liu**, Yuan Qin, **Tong Chen**, et. al. A practical guide to amplicon and metagenomic analysis of microbiome data. **Protein Cell**, 2021(12) 5:315-330, doi: [10.1007/s13238-020-00724-8](https://doi.org/10.1007/s13238-020-00724-8)

Copyright 2016-2023 Yong-Xin Liu <liuyongxin@caas.cn>, Tao Wen <taowen@njau.edu.cn>, Tong Chen <chent@nrc.ac.cn>