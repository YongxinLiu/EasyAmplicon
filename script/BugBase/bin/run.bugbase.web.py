#!/usr/bin/env python

# Runs BugBase analysis:
# - Uses a single-cell approach to predict microbiome phenotypes
# - plots phenotype relative abundances
# - prints statisical analyses
# - plots thresholds
# - plots otu contributions

# Options
# -t    Specify which taxa level to plot otu contributions by (number 1-7) 
# -p    Specify list of traits (phenotypes) to test (list, comma separated)
# -T    Specify a threshold to use for all traits (number 0-1)
# -g    Specify subset of groups in map column to plot (list, comma separated)
# -u  Use a user-define trait table. Absolute file path must be specified
# -k  Use the kegg pathways instead of default traits
# -z    Data is of type continuous 
# -C    Use covariance instead of variance to determine thresholds
# -a    Plot all samples (no stats will be run)
import site
import sys
import os
import random
import operator
import csv
from subprocess import Popen, PIPE, STDOUT
from optparse import OptionParser

# These are the environment paths added by the module:

# Python
os.environ['PATH'] = os.environ['PATH'] + '/soft/python-2.7/bin'

# PyCogent
os.environ['PYTHONPATH'] = ':/web/research/bugbase.cs.umn.edu/site-packages/PyCogent-1.5.3'

# R
#os.environ['PATH'] = os.environ['PATH'] + ':/soft/r/3.1.2/linux_x86_64/bin'
#os.environ['PATH'] = os.environ['PATH'] + ':/soft/r/3.1.2/linux_x86_64/rstudio-0.98.1087/bin'
os.environ['PATH'] = os.environ['PATH'] + ':/soft/r/3.2.2/linux_x86_64/bin'
os.environ['PATH'] = os.environ['PATH'] + ':/soft/r/3.2.2/linux_x86_64/rstudio-0.99.486/bin'

# BugBase
os.environ['PATH'] = os.environ['PATH'] + ':/web/research/bugbase.cs.umn.edu/BugBase/bin'
os.environ['BUGBASE_PATH'] = '/web/research/bugbase.cs.umn.edu/BugBase/'

def make_option_parser():
	parser = OptionParser(usage="usage: %prog [options] filename",
		version="%prog 1.0")
	parser.add_option("-v","--verbose",
		action="store_true",
		default=True,
		help="Verbose output (default %default)",)
	parser.add_option("-p","--print_only",
		action="store_true",
		default=False,
		help="Only print commands (default %default)",)
	parser.add_option("-i", "--input_OTU",
		default=None,
		type='string',
		help="OTU table (required)")
	parser.add_option("-m", "--mapping_file",
		default=None,
		type='string',
		help="mapping file (required)")
	parser.add_option("-o", "--output",
		default=".",
		type='string',
		help="Output directory (default %default)")
	parser.add_option("-c","--map_column",
		default=None,
		type='string',
		help="name of column that lists treatment groups")
	parser.add_option("-P","--phenotypes",
		default=None,
		type='string',
		help="specific phenotypes to predict, comma separated list, no spaces")
	parser.add_option("-x","--predict",
		default=False,
		action="store_true",
		help="only output the prediction table, do not make plots")
	parser.add_option("-T","--threshold",
		default=None,
		type='float',
		help="threshold (0 to 1) you would like to set for the trait listed")
	parser.add_option("-g","--groups",
		default=None,
		type='string',
		help="treatment groups you would like to plot, separated by commas with no spaces")
	parser.add_option("-t","--taxalevel",
		default=None,
		type='float',
		help="taxa level to plot otu contributions by, default is 2 (phylum)")
	parser.add_option("-a","--plot_all",
		action="store_true",
		default=False,
		help="Plot all samples without treatment-group seperation and no statistics")
	parser.add_option("-C","--cov",
		action="store_true",
		default=False,
		help="use coefficient of variance instead of variance")
	parser.add_option("-l","--clr",
		action="store_true",
		default=False,
		help="use centered log-ratio transformation instead of relative abundance")
	parser.add_option("-k","--kegg",
		action="store_true",
		default=False,
		help="use kegg pathway table")
	parser.add_option("-z","--continuous",
		action="store_true",
		default=False,
		help="Plot data according to a continuous variable (picked against IMG)")
	parser.add_option("-w","--wgs",
		action="store_true",
		default=False,
		help="Data is whole genome shotgun data")
	return parser

def run_commands(commands, print_only=False, verbose=True, error_on_fail=True):
	return_vals = []

	# Run all commands
	for cmd in commands:
		print cmd
		if not print_only:
			proc = Popen(cmd,shell=True,universal_newlines=True,stdout=PIPE,stderr=PIPE)
			stdout, stderr = proc.communicate()

		# If requested, prints all output from the program
			if verbose:
				print stdout
				print stderr
			if error_on_fail == True and proc.returncode != 0:
				print stdout
				print stderr
				raise ValueError('Command failed: ' + cmd)

			return_vals.append(proc.returncode)
	return(return_vals)

		
if __name__ == '__main__':

	parser = make_option_parser()
	(options, args) = parser.parse_args()

	# Name user inputs
	output = " -o /web/research/bugbase.cs.umn.edu/results/" + options.output

	otu_table = " -i /web/research/bugbase.cs.umn.edu/uploads/" + options.input_OTU 
	
	if options.mapping_file:
		mapping = " -m /web/research/bugbase.cs.umn.edu/uploads/" + options.mapping_file
		map_file = "/web/research/bugbase.cs.umn.edu/uploads/" + options.mapping_file
	else:
		mapping = ""
		
	if options.map_column:
		column = " -c " + options.map_column
		col_name = options.map_column
	else:
		column = ""

	if options.phenotypes:
		phenos = " -p " + options.phenotypes
	else:
		phenos = ""

	if options.predict:
		predict_only = " -x"
	else:
		predict_only = ""

	if options.threshold:
		thres = " -T %s" %(options.threshold)
	else:
		thres = ""

	if options.groups:
		groups = " -g " + options.groups
		groups_list = options.groups.split(",")
	else:
		groups = ""

	if options.taxalevel:
		taxa = " -t %s" %(options.taxalevel)
	else:
		taxa = ""

	if options.plot_all is True:
		plot_all = " -a"
	else:
		plot_all = ""

	if options.cov is True:
		cov = " -C"
	else:
		cov = ""

	if options.clr is True:
		clr = " -l"
	else:
		clr = ""

	if options.kegg is True:
		kegg = " -k"
	else:
		kegg = ""

	if options.continuous is True:
		continuous = " -z"
	else:
		continuous = ""

	if options.wgs is True:
		wgs = " -w"
	else:
		wgs = ""

	if options.plot_all is False:
		if options.predict is False:
			if options.mapping_file is None:
				print "[ERROR_MESSAGE]Mapping file must be specified"
				sys.exit()
			if options.map_column is None:
				print "[ERROR_MESSAGE]Column header must be specified"
				sys.exit()

			# Make sure map column is valid
			with open(map_file, 'rU') as input_map:
				reader = csv.reader(input_map, delimiter='\t')
				headers = reader.next()
			if col_name in headers:
				print col_name + " was specified as map column header\n"
			else:
				print "[ERROR_MESSAGE]Column header specified does not exist in mapping file\n"
				print "[ERROR_MESSAGE]These are the available column headers: "+ ', '.join(headers)
				sys.exit()
		
			# If groups are specified, check they are valid
			if options.groups is not None:
				groups_avail = []
				with open(map_file, 'rU') as input_map:
					reader = csv.reader(input_map, delimiter='\t')
					headers = reader.next()
					column_index = headers.index(col_name)
					for row in reader:
						name = str(row[column_index])
						groups_avail.append(name)   
				for group_defined in groups_list:
					if group_defined in groups_avail:
						if len(groups_list) <= 1:
							print "[ERROR_MESSAGE]A minimum of two groups must be tested"
							sys.exit()
					else:
						groups_avail = list(set(groups_avail))
						print "[ERROR_MESSAGE]Groups specified do not exist in mapping file"
						print "[ERROR_MESSAGE]These are the groups available under " + col_name + " header: " + ', '.join(groups_avail)
						sys.exit()
	
	# If threshold is user-specified, state what will be used
	#if options.threshold is not None:
	#	print "A user-specified threshold of %s will be used for all traits" %(options.threshold)

	# If kegg is specified, check that modules are valid
	if options.kegg is True:
		if options.predict is False:
			if options.phenotypes is None:
				print "[ERROR_MESSAGE]When using KEGG, modules to plot must be specified with a comma separated list"

	commands = []
						
	# Run the run.bugase.web.r script
	
	#cmd = "/soft/r/3.2.2/linux_x86_64_old/bin/Rscript /web/research/bugbase.cs.umn.edu/BugBase/bin/run.bugbase.web.r -h" 
	cmd = "Rscript /web/research/bugbase.cs.umn.edu/BugBase/bin/run.bugbase.web.r" + otu_table + mapping + column + phenos + predict_only + taxa + groups + thres + plot_all + cov + clr + kegg + wgs + continuous + output 
	#cmd = "/soft/r/3.2.2/linux_x86_64_old/bin/Rscript /web/research/bugbase.cs.umn.edu/BugBase/bin/test.r" + otu_table + mapping + column + phenos + predict_only + taxa + groups + thres + plot_all + cov + kegg + wgs+ continuous + output
	commands.append(cmd)
	
	# Run commands
	return_vals = run_commands(commands, print_only=options.print_only, verbose=options.verbose)
		
	print "[SUCCESSFUL] Data is analyzed"
