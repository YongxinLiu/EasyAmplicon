# table2itol

## About

Interactive Tree of Life ([iTOL](http://itol.embl.de/)) is a popular tool for
displaying phylogenetic trees and associated information. The `table2itol.R`
script makes it easy to generate iTOL annotations from spreadsheet files.

## Features

* Works with [CSV](https://en.wikipedia.org/wiki/Delimiter-separated_values),
  OpenOffice, LibreOffice and Microsoft Excel files.
* Supports iTOL domains, colour strips, simple bars, gradients, binary data,
  heat maps, and texts.
* Partially supports iTOL branch annotation (currently work in progress).
* By default selects the appropriate visualisation from the data type of each
  input column but this can be modified by the user.
* Provides carefully chosen colour vectors for up to 40 levels and optionally
  combines them with symbols for maximizing contrast.
* The default colour vectors can be replaced by user-defined colour vectors.
* Can be used either interactively on any operating system on which R is
  running, or non-interactively using the command line of a UNIX-like system.

## Prerequisites

* A recent (>= 3.2.0) version of [R](https://cran.r-project.org/).
* The [optparse](https://CRAN.R-project.org/package=optparse) package for R if 
  you want to run the script in non-interactive mode or if you want to read the
  help message (which is placed in `tests/table2itol_help.txt`, too).
* The [plotrix](https://CRAN.R-project.org/package=plotrix) package for R if
  you want to generate branch annotations from continuous numeric data.
* The [readxl](https://CRAN.R-project.org/package=readxl) package for R if
  you want to apply the script to Microsoft Excel files.
* The [readODS](https://CRAN.R-project.org/package=readODS) package for R if
  you want to apply the script to Libreoffice or Openoffice
  [ods](https://en.wikipedia.org/wiki/OpenDocument) files.
* The [yaml](https://CRAN.R-project.org/package=yaml) package for R if
  you want to define colour vectors yourself.

*Please note that explaining how to correctly install R is beyond the scope of
this manual, and please do not contact the `table2itol.R` authors about this
issue. There is plenty of online material available elsewhere. As for the
installation of R packages see the FAQ below.*

## Installation

First, obtain the script as indicated on its GitHub page.

### Command-line use

The following explanations are for *non-experts*; there is nothing special with 
running this script in command-line mode on UNIX-like systems. First, if
necessary make the script executable:

`chmod +x table2itol.R`

Then call:

`./table2itol.R`

to obtain the help message. If this yields an error, see the **troubleshooting**
chapter.

Optionally place the script in a folder that is contained in the `$PATH`
variable, e.g.

`install table2itol.R ~/bin`

or even

`sudo install table2itol.R /usr/local/bin`

if you have sudo permissions. Then you can call the script by just entering

`table2itol.R`

### Interactive use

Open R or [RStudio](https://www.rstudio.com/) or whatever interface to R you
are using, then enter at the console:

```R
source("table2itol.R")
```

provided the script is located in the current working directory as given by
`getwd()`. Alternatively, first use `setwd()` to move to the directory in which
`table2itol.R` resides or enter the full path to the location of the script.

When loading the script it shows the usual help message and an indication that 
you are running it in interactive mode. When doing so, you might need to modify 
the arguments of the function much like command-line users might need to apply
certain command-line options. For instance, in analogy to entering:

`./table2itol.R --na-strings X --identifier Tip --label Name ann1.tsv ann2.tsv`

on the command line of a UNIX-like system, you would enter within R the
following:

```R
source("table2itol.R")
create_itol_files(infiles = c("ann1.tsv", "ann2.tsv"),
  identifier = "Tip", label = "Name", na.strings = "X")
```

The analogy should be obvious, hence for details on the arguments of 
`create_itol_files` see the help message. The arguments of the function are 
identical to the long version of the arguments of the script, subjected to the 
replacement of dashes by dots to yield syntactic names. The sole mandatory 
argument of the function is `infiles`, whose value is identical to the
positional arguments of the script. With some basic knowledge of R it is thus
easy to set up customized scripts that set the arguments for your input files
and generate the intended output.

## Examples

Exemplars for input table files are found within the `tests/INPUT` folder. A
list of examples for calling `table2itol.R` is found in `tests/examples.txt`.

*Experts only*: On a UNIX-like system you can run these examples by calling 
`tests/run_tests.sh` provided a modern
[Bash](https://www.gnu.org/software/bash/) is installed. The versions of R and
the R packages used for testing by the maintainer are found in the file 
`tests/R_settings.txt`.

## Troubleshooting

Some commonly encountered error messages are mentioned in the following. Note
that you might actually get these error messages in a language other than
English (e.g., your own language) or with other minor modifications.

### Command-line use

#### Bad interpreter

`/usr/local/bin/Rscript: bad interpreter: No such file or directory`

Solution: Enter

`locate Rscript`

and watch the output. If it is empty, you must install
[R](https://cran.r-project.org/) first. If you instead obtained a location such
as `/usr/bin/Rscript` you could do the following:

`sudo ln -s /usr/bin/Rscript /usr/local/bin/Rscript`

if you had sudo permissions. Alternatively, within the first line of the script
replace `/usr/local/bin/Rscript` by `/usr/bin/Rscript` or wherever your
`Rscript` executable is located. A third option is to leave the script as-is and
enter `Rscript table2itol.R` instead of `./table2itol.R` or whatever location of
the script you are using. But this is less convenient in the long run.

*Please note that explaining how to correctly install R is beyond the scope of
this manual, and please do not contact the `table2itol.R` authors about this
issue. There is plenty of online material available elsewhere.*

### Command-line or interactive use

#### Missing R package

`there is no package called 'optparse'`

Solution: Install the [optparse](https://CRAN.R-project.org/package=optparse) 
package for R. (It is not an absolute requirement in interactive mode but 
without it you would not see the help message. However, this message is placed
in `tests/table2itol_help.txt` anyway.)

`there is no package called 'plotrix'`

Solution: Install the [plotrix](https://CRAN.R-project.org/package=plotrix)
package for R. (It is only needed if you want to create branch annotations
from continuous numeric data.)

`there is no package called 'readODS'`

Solution: Install the [readODS](https://CRAN.R-project.org/package=readODS)
package for R. (It is only needed if you want to apply the script to
[ods](https://en.wikipedia.org/wiki/OpenDocument) files.)

`there is no package called 'readxl'`

Solution: Install the [readxl](https://CRAN.R-project.org/package=readxl)
package for R. (It is only needed if you want to apply the script to Microsoft
Excel files.)

`there is no package called 'yaml'`

Solution: Install the [yaml](https://CRAN.R-project.org/package=yaml)
package for R. (It is only needed if you want to use the script in conjunction
with your own colour vectors.)

*Please note that explaining how to correctly install R is beyond the scope of
this manual, and please do not contact the `table2itol.R` authors about this
issue. There is plenty of online material available elsewhere. As for the
installation of or R packages see the FAQ  below.*

#### Outdated R version

`need a newer version of R, 3.2.0 or higher`

Solution: Install a newer version of [R](https://cran.r-project.org/).

*Please note that explaining how to correctly install R is beyond the scope of
this manual, and please do not contact the `table2itol.R` authors about this
issue. There is plenty of online material available elsewhere.*

#### The script generates not enough output files

Solution: Watch the warnings and error messages generated by the script. Without
any input files, the script *should* not generate any output. The script would
also skip input files or single tables if they failed to contain columns you
have requested. You can use the `--abort` option to let the script immediately
stop in such cases, then look up the last error message in this manual. But even
without `--abort` the script generates warnings when data sets get skipped.

#### The script generates too many output files

Solution: Accept as a design decision that the scripts generates one file for 
each input column (except for the tip identifier column and when a heatmap is
created). Since you can still decide to not upload (some of) the generated files
to iTOL and also deselect data sets within iTOL, we believe it would not make
much sense to also include a selection mechanism within the `table2itol.R`
script. As last resort you could also reduce the number of input columns.
However, if you are mainly concerned about the script cluttering up your working
directory with files, simply consider using the `--directory` option to place
all output files in a dedicated directory. An empty argument to this option
causes the script to place every output file in the directory in which the
respective inout file resides.

#### A column is requested but missing

`selected column 'ID' does not exist`

Solution: Use the `--identifier` option to set the name of the tip identifier
column.

`selected column 'Label' does not exist`

Solution: Use the `--label` option to set the name of the tip label column.

#### A name clash of output file names occurs

`name clash: file [...] has already been generated`

Solution: Name the columns distinctly in distinct tables within the same file
and/or call the `table2itol.R` script individually for each input file, maybe
best with distinct values of `--directory`.

Since the name of each output file is generated from the name of the respective
input column and the resulting kind of visualisation, columns from distinct
input tables but with the same name and the same resulting kind of visualisation
would yield only a single output file. Instead of silently overwriting the
earlier ones, the script stops with an informative error message.

## Frequently asked questions

#### How can I install the required R packages?

There are several ways to obtain, install and update R packages. From within R,
we found the following approach to be convenient:

```R
source("http://bioconductor.org/biocLite.R")
biocLite(c("optparse", "plotrix", "readODS", "readxl", "yaml"))
```

This should under normal circumstances install or update all R packages
recommended for `table2itol.R`.

*Please note that explaining in greater detail how to correctly install R
packages is beyond the scope of this manual, and please do not contact the
`table2itol.R` authors about this issue. There is plenty of online material
available elsewhere. The `table2itol.R` authors cannot guarantee that the script
provided by [BioConductor](http://bioconductor.org/) works as expected.*

#### How can I generate other kinds of visualisation from integer columns?

For generating distinct kinds of visualisation from distinct integer columns,
run the script several times with distinct values of `--conversion`. Since the
name of each output file is generated from the name of the respective input
column *and* the resulting kind of visualisation, nothing of importance will be
overwritten (but see the section on name clashes between distinct spreadsheets).

#### How can I generate other kinds of visualisation from non-integer numbers?

You can create bar charts instead of gradients from numbers with decimal points
by using the `--double-to-bars` option. Note that the number of decimal points
of the range as shown in the legend can be modified using the `--precision`
option, as usual.

#### How can I define my own colour vectors?

For replacing the default colour vectors by other colour vectors, use the 
`--colour-file` option. Its argument must be the name of a file in 
[YAML](http://yaml.org/) format. The `tests/INPUT` folder contains an example
file with user-defined colours.

Names for the colour vectors are optional in such files but increase
readability. Not all colour vectors need to be defined, only those that should
replace certain default colour vectors. Assignment is solely by vector length.

An attempt is made to standardize the input colours, yielding hexadecimal codes 
understood by iTOL. Thus all kinds of colour specifications can be used that are
accepted by the `col2rgb` function. Non-interpretable and duplicate colours
yield an error. Call `colors()` to obtain the list of human-readable colour
names accepted by R. You might also want to visit the 
[colorbrewer](http://colorbrewer2.org) web site for generating useful colour 
vectors.

#### How can I distinctly assign symbols?

For assigning symbols (to be used in combinations of symbols and colours)
according to certain groups, use the `--emblems` option. This does not
individually assign symbols but the symbols will then be consistently assigned
according to some column. For triggering the earlier use of symbols (instead of
using more colours), play around with `--favour` and `--max-size`.

#### How can I assign specific colours to binary data?

Binary data are interpreted by the script either as factors with two levels or 
as logical vectors. Logical vectors yield distinct kinds of output files, and 
the script already picks distinct colours and symbols for distinct logical 
vectors in turn. To get columns understood as logical vectors their fields must 
(in addition to NA values, if any) only contain values from one of the following
pairs: `0`/`1`, `true`/`false`, `t`/`f`, `yes`/`no`, `y`/`n` or `on`/`off`. This
considerably extends the values recognized by base R as logical vectors on
input. Moreover, case differences do not matter, but you cannot mix any of these
variants. These conversions get turned off when using `--conversion keep` or
`--conversion double`.

Colours defined using `--colour-file` only play a role for columns treated as 
factor, not for those yielding a logical vector. To set the colours used for end
points as gradients as well as for binary data, use the `--gradient-file`
argument. In contrast, one cannot modify the symbols used; this would not make
much sense though, since iTOL understands only certain symbols anyway.

#### How can I generate a heat map?

A heat map is automatically generated from a data frame when all columns (except
for the special columns addressed using `--identifier`, `--label` and 
`--background`) are numeric (uniformly of mode double or uniformly of mode 
integer). These columns get merged into a single file whose name is derived from
the name of the identifier column. You might need the `--conversion` argument to
obtain columns that are uniformly numeric.

#### The colours appear too strong. What can I do?

The argument `--opacity` can be used to defined the so-called [alpha 
channel](https://en.wikipedia.org/wiki/Alpha_compositing) of each colour. Set 
the opacity to a value between 0 and 1 to obtain more transparent colours (1 
means fully opaque and 0 means fully transparent, which is not normally very 
useful). Also note that some colours can be modified with few clicks in iTOL 
itself. But instead of attempting to modify the files generated by
`table2itol.R` with tools such as `sed` or by hand we suggest to contact the
`table2itol.R` authors with an according feature request.

