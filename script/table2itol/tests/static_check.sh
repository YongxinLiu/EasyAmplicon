#!/bin/bash


################################################################################
#
# static_check.sh -- Bash script for checking the table2itol.R script.
#
# (C) since 2017 by Markus Goeker (markus [DOT] goeker [AT] dsmz [DOT] de)
#
# This program is distributed under the terms of the Gnu Public License. See
# http://www.gnu.org/licenses/gpl.html for further information.
#
################################################################################


set -eu


################################################################################


function static_check
{
  local infile outfile tmpfile
  declare -i errs=0
  tmpfile=$(mktemp --suffix ."$FUNCNAME")
  outfile=$1
  shift
  for infile; do
    if R --interactive --slave > /dev/null <<-______EOF

# the file to inspect
source("$infile")

# collect all "::" calls and output their package versions
relevant_package_versions <- function(env = globalenv()) {
  packages_called <- function(expr, result) {
    rec_collect <- function(x, result, wanted) {
      if (!is.recursive(x))
        return(result)
      if (is.call(x) && identical(x[[1L]], wanted))
        result[[deparse(x[[2L]])]] <- TRUE
      lapply(x, rec_collect, result, wanted)
      result
    }
    rec_collect(expr, result, as.symbol("::"))
  }
  result <- new.env(TRUE, emptyenv())
  for (thing in lapply(ls(env), get))
    if (is.function(thing))
      packages_called(body(thing), result)
  result <- sort.int(names(result))
  result <- sapply(X = result, FUN = packageVersion, simplify = FALSE)
  c(R = as.character(getRversion()), vapply(result, as.character, ""))
}
cat(formatDL(relevant_package_versions()), sep = "\n", file = "$outfile")

# use codetools for a static check but remove irrelevant issues
customised_check <- function(env = globalenv()) {
  problems <- character()
  codetools::checkUsageEnv(env = env, all = TRUE, report = function(s)
    problems <<- c(problems, s), suppressParamAssigns = TRUE)
  # filter out known, irrelevant issues detected in table2itol.R
  ok <- paste0("\\\\b(emit_\\\\w+|GENERATED_FILES)\\\\W+(assigned but|",
    "parameter\\\\W+(ids|outdir|x))\\\\W+may not be used\\\\b")
  problems[!grepl(ok, problems, FALSE, TRUE)]
}
cat(paste0(customised_check(), collapse = "\n"), file = "$tmpfile")

quit("no", 0L)

______EOF
    then
      true
    else
      let errs+=1
      echo 'ERROR: call of R did not result' >&2
    fi
    if [ -s "$tmpfile" ]; then
      let errs+=1
      echo "PROBLEMS in file $infile detected by static check:"
      cat "$tmpfile"
    else
      echo "File $infile seems alright."
    fi
    echo
  done
  rm -f "$tmpfile"
  [ $errs -gt 0 ] && return 1 || return 0
}


################################################################################


[ -z "${0%/*}" ] || cd "${0%/*}"


[ $# -gt 0 ] || set -- ../table2itol.R


for infile; do
  helpfile=${infile##*/}
  helpfile=${helpfile%.*}_help.txt
  "$infile" -h > "$helpfile" || true
done


static_check R_settings.txt "$@"


