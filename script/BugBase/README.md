## BugBase Usage and Installation

### Standard Users
Standard users can analyze their microbiome samples with the default phenotypes, as well as the KEGG pathways. 

### Dependencies

BugBase relies on R packages. When you are ready to try BugBase, the following R packages will be looked for on your computer, and installed for you if they are missing.

* optparse
* beeswarm
* RColorBrewer
* reshape2
* plyr
* grid
* gridExtra
* ggplot2
* RJSONIO
* biom
* Matrix
* labeling
* digest

### Installation
#### Mac OS
You can download BugBase here, and then add the paths to BugBase and the BugBase bin to your `~/.bash_profile` file. This is what is in an example `~/.bash_profile` modification looks like:

```
export BUGBASE_PATH=/Path/to/my/BugBase
export PATH=$PATH:$BUGBASE_PATH/bin
```

Note: you will need change the paths to match your system and which directory you have downloaded BugBase into. You might need to put it in ~/.bashrc instead of ~/.bash_profile depending on your system. After adding these paths to the .bash_profile or ~/.bashrc, reopen the terminal or login again.

To check your install, type the following in the command line.  It will first install any missing packages, and then it will print the options available to run BugBase.

```
run.bugbase.r -h 
```

### Demo
BugBase comes with a test dataset in the bugbase/doc/data directory. To analyze this data set you would type the following:

```
run.bugbase.r -i $BUGBASE_PATH/doc/data/HMP_s15.txt -m $BUGBASE_PATH/doc/data/HMP_map.txt -c HMPBODYSUBSITE -o output
```

You can view other options with `run.bugbase.r -h`.

### Using BugBase 

BugBase has one main command, `run.bugbase.r`, that will:
-	Normalize your OTU table according to 16S copy number (WGS data will not be normalized)
-	Plot the variance in phenotype possession for thresholds 0-1
-	Determine which threshold to set for each microbiome phenotype
-	Determine the proportion of each microbiome with a given phenotype
-	Plot the proportions of the microbiome with a given phenotype
-	Statistically analyze the microbiome phenotype proportions according the treatment groups specified, or by using regression for continuous data
-	Plot OTU contributions for each phenotype

BugBase's `run.bugbase.r` parameters:
```
Required
	-i	Input OTU table, picked against the Greengenes database (16S) or IMG (Shotgun) (.txt or .biom (json))
	-o	Output directory name
	
Optional
	-m	Mapping file (tab-delimited text file)
	-c	Map column header to plot by (which column denotes treatment groups)
	-w	Data is shotgun metagenomic data (picked against IMG database)
	-a 	Plot all samples (no stats will be run)
	-x	Output prediction files only, no plots will be made
	-g 	Specify subset of groups in map column to plot (list, comma separated)
	-z 	Data is of type continuous 
	-C 	Use coefficient of variance instead of variance to determine thresholds
	-l 	Centered log-ratio transform the data instead of using relative abundance
	-t	Taxa level to plot OTU contributions by (number 1-7)
	-T 	Specify a threshold to use for all traits (number 0-1)
	-k 	Use the KEGG modules instead of default traits (Note: you must specify which modules!)
	-p 	List modules or traits to predict (comma separated list, no spaces)
	-u	Use a user-define trait table. Absolute file path must be specified
	
```

### BugBase compatible OTU tables

BugBase takes in QIIME compatible OTU tables in the classic (.txt) or json version (.biom).  Below are some instructions regarding OTU table generation for downstream BugBase use.

16S data:
- Closed reference OTU picking with the Greengenes 97% representative set
- Do not include taxonomy if in the classic form (.txt)
- Counts, not relative abundance
- You can use OTU tables generated from many programs, such as [NINJA-OPS.](https://github.com/GabeAl/NINJA-OPS "NINJA") or [QIIME](http://qiime.org/ "qiime")

Shotgun data:
- Closed reference OTU picking with the IMG reference sequences

To generate a BugBase compatible OTU table from shotgun data, please follow the steps below:

1. Download the UTree release specific to your operating system by downloading and unzipping version 0.0.1 of UTree [here.](https://github.com/knights-lab/UTree/releases/tag/v0.0.1 "UTree")
2. Install miniconda, if you don't already have it. Miniconda can be found [here.](https://conda.io/miniconda.html "miniconda")
3. Download and unzip directory for the SHOGUN-BugBase database (IMG reference sequences and maps) needed for OTU picking [here.](https://drive.google.com/open?id=0ByVmiknyDGaiM3M0dDBJMkZuZDg "shogun-bugbase-db") 

4. Finish installing SHOGUN for BugBase with the following commands:
```
# create a shogun conda environment
conda create -n shogun python=3   

# activate the shogun environment
source activate shogun 

# install shogun for BugBase
pip install git+https://github.com/knights-lab/SHOGUN.git@v0.0.1 --no-cache-dir --upgrade

pip install git+https://github.com/knights-lab/NINJA-utils

pip install git+https://github.com/knights-lab/NINJA-dojo

# deactivate the environment
source deactivate
```

5. Run OTU picking with the following commands.  Update the `shogun_bugbase` command to be specific to the filepaths for your input sequences and the SHOGUN-BugBase database you downloaded.  Your input sequences should be in one directory, with one .fna file per sample. The name of each .fna file should be the name of the sample it corresponds to. Once OTU picking is complete, you will have an OTU table in classic format (.txt) called 'taxa_counts.txt' within the output directory specified.
```
# activate the shogun environment
source activate shogun 

# run OTU picking with shogun
shogun_bugbase -i path_to_sequences -o output_path -u path_to_shogun_bugbase_db

# deactivate the shogun environment
source deactivate
```

### Creating user-defined trait tables for predictions

You can create your own traits of interest for BugBase predictions using KEGG orthologies.  To do so, you will need the following:

- A file for each trait that consists of the KO IDs involved in the trait/pathway, one KO ID per line
- A directory that houses the trait files mentioned about (one per trait), each trait file name will be used as the trait name (exactly) in the BugBase table created

You can create your user-defined custom BugBase input table using `make.user.table.r` that will:
- create intermediate files for each trait specified
- merge all intermediate tables into one table that has each trait of interest as a column
- create a final BugBase input file will be called "Custom_BugBase_Traits.txt" and it will be located in directory you specified as the input

BugBase parameters for `make.user.table.r`:
```
Required
	-i     path to directory housing the files that list the KO IDs per trait
	
Optional
	-w	traits are for shotgun metagenomic sequencing, default is 16S
```

To create a custom BugBase input and run the BugBase predictions:

```
#Example:

make.user.table.r -i directory_with_trait_files

run.bugbase.r -i path/to/OTU_Table.biom -m path/to/map.txt -c metadata_column -u directory_with_trait_files/Custom_BugBase_Traits.txt -o output_name

```

=======
