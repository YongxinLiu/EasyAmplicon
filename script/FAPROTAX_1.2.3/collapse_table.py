#!/usr/bin/env python3
#
# Please read the copyright notice and license agreement below before using or sharing this script.
#
# This python script collapses a data table's rows or columns, down to a custom set of row-groups.
# This is an improved version of the previous script "collapse_table_rows.py" (not published).
# The input table can either be in BIOM format or in classical tabular format (e.g. tab-separated).
# The output table will be in classical tabular format.
#
# Use as in the following examples:
#	./collapse_table.py -i taxonomic_otu_table.biom -g FAPROTAX_database.txt -o functional_otu_table.tsv -r report.txt --group_leftovers_as 'other' --normalize_collapsed 'columns_after_collapsing' -v
#	./collapse_table.py -i taxonomic_otu_table.tsv --groups_file FAPROTAX_database.txt -f -o functional_otu_table.tsv -r report.txt --column_names_are_in last_comment_line  --keep_header_comments --non_numeric consolidate -v --row_names_are_in_column "taxonomy" --omit_columns 0 --normalize_collapsed columns_before_collapsing --group_leftovers_as 'other'
# See the help menu (-h) of the script for more information.
#
# Tested on python 3.7.6, Mac OS 10.13.6
#
# Script version: 1.3.2
# May 01, 2020
#
# Copyright (c) 2020, Stilianos Louca
#
# LICENSE AGREEMENT
# All rights reserved.
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
#
#    * Redistributions of source code must retain the above copyright
#      notice, this list of conditions and the following disclaimer.
#    * Redistributions in binary form must reproduce the above copyright
#      notice, this list of conditions and the following disclaimer in
#      the documentation and/or other materials provided with the distribution
#    * Neither the name of the original author (Stilianos Louca), nor the names
#      of other contributors may be used to endorse or promote products derived
#      from this software without specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
# ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
# LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
# CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
# SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
# INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
# CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
# ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
# POSSIBILITY OF SUCH DAMAGE.


# DEPENDENCIES
import string
import time
import sys, os
import argparse
import numpy
import shutil
import re
import fnmatch
import gzip

from numpy import NaN


# OPTIONAL DEPENDENCIES
try:
	import biom
	HAS_BIOME_MODULE = True;
except ImportError as e:
	HAS_BIOME_MODULE = False;
	name = 'biom'
	if e.message != 'No module named %s'%(name):
		print("WARNING: An error occurred while importing the module '%s': %s\n         This script will only work with classical (e.g. tab-separated) tables"%(name,e.message))
	else:
		print("WARNING: The module '%s' was not found.\n         This script will only accept classical (e.g. tab-separated) tables"%(name))


try:
	import json
except ImportError as e:
	name = 'json'
	if e.message != 'No module named %s'%(name):
		print("WARNING: An error occurred while importing the module '%s': %s\n         This script may not work with tables in BIOM JSON format"%(name,e.message))
	else:
		print("WARNING: The module '%s' was not found.\n         This script may not work with tables in BIOM JSON format"%(name))


try:
	import h5py
except ImportError as e:
	name = 'h5py'
	if e.message != 'No module named %s'%(name):
		print("WARNING: An error occurred while importing the module '%s': %s\n         This script may not work with tables in BIOM HDF5 format"%(name,e.message))
	else:
		print("WARNING: The module '%s' was not found.\n         This script may not work with tables in BIOM HDF5 format"%(name))


# Python 2 & 3 universality issues
try:
  str
except NameError:
  str = str
  

try:
    xrange
except NameError:
    xrange = range



############################
######## AUXILIARY #########



# split a list of indices corresponding to given sortingKeys[] and scores[] into multiple partitions[][]
# partitions[p][] contains the indices of all items with scoreThresholds[p]<=score<scoreThresholds[p+1] (the last+1 scoreThreshold being formally INFTY)
# within each partitions[p][], items are sorted alphabetically according to their sortingKeys[]
# Note: If scoreThresholds[] is not ascending, then partitions may overlap.
def partitionIndexListByScores(	indices,
								sortingKeys,
								scores,
								scoreThresholds):
	NP = len(scoreThresholds);
	NI = len(indices);
	# create partitions:
	partitions = [[i for i in indices if (scoreThresholds[p]<=scores[i]) and (p>=NP-1 or scores[i]<scoreThresholds[p+1])] for p in range(NP)];
	# sort alphabetically within partitions:
	partitions = [sorted(partitions[p], key=lambda i: sortingKeys[i].lower()) for p in range(NP)]
	return partitions;

								

def whichPrefix(haystack, candidate_prefixes):
	for i,p in enumerate(candidate_prefixes):
		if(haystack.startswith(p)): return i;
	return -1;


set_operations_keywords = ['add_group:', 'subtract_group:', 'intersect_group:']


def object2str(x):
	if(isinstance(x, str)):
		return x
	elif(isinstance(x, float)):
		return "%.10g"%(x)
	else:
		return str(x)


#####################################
######## DATA INPUT #########


	
def is_biom_file(path):
	path = path.lower()
	return (path.endswith(".biom") or path.endswith(".biom.gz") or path.endswith(".hbiom") or path.endswith(".hbiom.gz") or path.endswith(".jbiom") or path.endswith(".jbiom.gz"))


def is_hbiom_file(path):
	path = path.lower()
	return (path.endswith(".hbiom") or path.endswith(".hbiom.gz"))


# save a BIOM table to file
# format can be "auto", "BIOM" or "HBIOM"
def save_biom_table(table, filepath, format):
	generator_name = os.path.basename(__file__)
	if((format=="HBIOM") or ((format=="auto") and is_hbiom_file(filepath))):
		# write as HDF5-BIOM
		if(filepath.lower().endswith(".gz")):
			# first write to a temporary file, without outside gzip compression, then compress the temporary file
			# this is needed because the h5py module writes directly to the file and cannot work wrapped around a gzip-file handle
			# since we are responsible for deleting the temp-file after we're done, enclose any actions on it in a try/finally statement, in case we run into an exception along the way
			[temp_id, temp_path] = tempfile.mkstemp(text=False)
			os.close(temp_id)
			try:
				with h5py.File(temp_path, 'w') as fout:
					table.to_hdf5(h5grp=fout, generated_by=generator_name, compress=True)
				# gzip-compress the temporary file to its final destination
				with open(temp_path, 'rb') as fin, gzip.open(filepath, 'wb') as fout:
					fout.writelines(fin)
			finally:
				# delete the temporary file, even if an exception occurs
				os.remove(temp_path)
		else:
			with h5py.File(filepath, 'w') as fout:
				table.to_hdf5(h5grp=fout, generated_by=generator_name, compress=True)
	else:
		# write as JSON-BIOM
		with (gzip.open(filepath,'wt') if filepath.lower().endswith(".gz") else open(filepath,'wt')) as fout:
			table.to_json(generated_by=generator_name, direct_io=fout)



	
def is_number(s):
    try:
        float(s)
        return True
    except ValueError:
        return False
        
        
	
def is_non_nan_number(s):
	return (s.lower() is not 'nan') and is_number(s);
	
	
def is_number_or_nan(s):
	s = s.lower();
	if((s=='nan') or (s=='na') or (s=='null')): return True;
	return is_number(s);
	
	
def float_or_nan(s):
    try:
        f = float(s)
        return f
    except ValueError:
        return NaN
      
        
def float_or_zero_if_nan(s):
    try:
        f = float(s)
        return f
    except ValueError:
        return 0	

	
numpy_float_or_nan = numpy.vectorize(float_or_nan)



def filter_index_list(	N, 
						only_items, 		# list of integers, can be None
						omit_items):		# list of integers, can be None
	items_to_keep=list(range(N))
	if(only_items is not None):
		items_to_keep = [i for i in only_items if i<N]
		items_to_keep.sort();	
	if(omit_items is not None):
		items_to_keep = [i for i in items_to_keep if not (i in omit_items)]
	return items_to_keep




def split_comments(line, comment_prefix):
	if(comment_prefix is ""): return [line];
	pos = line.find(comment_prefix);
	if(pos<0): return [line];
	return [line[0:pos], line[pos+len(comment_prefix):]];



def split_at_first_whitespace(line):
	parts = line.split(None, 1)
	if(len(parts)==1): return parts[0],""
	else: return parts[0],parts[1]



#####################################################
# Filtering samples & observations by name & metadata

# return a list of indices corresponding to items to keep
def filter_name_list(	names,
						items_to_keep,				# starting pool to further filter. This will typically be range(len(names)), but can also be a subset thereof in case of iterative filtering. Can also be None (equivalent to range(len(names)))
						only_items_with_name, 		# list or set of strings, can be None. Name wildcards to keep.
						omit_items_with_name,		# list or set of strings, can be None. Name wildcards to exclude.
						case_sensitive):
	if(items_to_keep is None): items_to_keep=list(range(len(names)))
	if(not case_sensitive): names = [name.lower() for name in names];
	if(only_items_with_name is not None):
		if(not case_sensitive): only_items_with_name = [name.lower() for name in only_items_with_name];
		items_to_keep = [i for i in items_to_keep if (next((w for w in only_items_with_name if fnmatch.fnmatchcase(names[i],w)),-1)>=0)]
	if(omit_items_with_name is not None):
		if(not case_sensitive): omit_items_with_name = [name.lower() for name in omit_items_with_name];
		items_to_keep = [i for i in items_to_keep if (next((w for w in omit_items_with_name if fnmatch.fnmatchcase(names[i],w)),-1)<0)]
	return items_to_keep


def get_metadata_predicates(metadata_predicate_list_str):
	if(metadata_predicate_list_str==""): return None,None;
	predicates = list(filter(len, metadata_predicate_list_str.split(';')))
	predicate_parts = [predicate.split(':',1) for predicate in predicates]
	metadata_names = [pp[0] for pp in predicate_parts if (len(pp)>1)]
	metadata_values = [pp[1].split(',') for pp in predicate_parts if (len(pp)>1)]
	return metadata_names,metadata_values


def filter_by_name_and_metadata(names,
								metadata,
								keep,			# a starting pool of indices to keep. Can be None (equivalent to range(len(names))).
								only_names,
								omit_names,
								only_by_metadata,
								omit_by_metadata,
								case_sensitive):
		# filter by name
		only_names = (None if (only_names=="") else only_names.split(','))
		omit_names = (None if (omit_names=="") else omit_names.split(','))
		keep = filter_name_list(names, keep, only_names, omit_names, case_sensitive)
		
		# filter by metadata
		if((metadata is not None) and (only_by_metadata!="" or omit_by_metadata!="")):		
			only_by_metadata_names,only_by_metadata_values = get_metadata_predicates(only_by_metadata)
			omit_by_metadata_names,omit_by_metadata_values = get_metadata_predicates(omit_by_metadata)
			if((only_by_metadata_names is not None) and (only_by_metadata_values is not None)):
				for mi,metadata_name in enumerate(only_by_metadata_names):
					filter_metadata = [(str(metadata[c][metadata_name]) if (metadata_name in metadata[c]) else "") for c in range(len(names))]
					keep = filter_name_list(filter_metadata, keep, only_by_metadata_values[mi], None, case_sensitive)
			if((omit_by_metadata_names is not None) and (omit_by_metadata_values is not None)):
				for mi,metadata_name in enumerate(omit_by_metadata_names):
					filter_metadata = [(str(metadata[c][metadata_name]) if (metadata_name in metadata[c]) else "") for c in range(len(names))]
					keep = filter_name_list(filter_metadata, keep, None, omit_by_metadata_values[mi], case_sensitive)
		
		return keep


#####################################################





	

def read_table(	file, 
				delimiter, 
				comment_prefix, 
				row_names_are_in_column, # can either be a column name, or an index if column_names_are_in=='none', or empty (no row names available). Only applies to classical tables.
				column_names_are_in, 
				only_columns, 	# list of integers (indices). Can be None (include all columns)
				omit_columns,	# list of integers (indices). Can be None (don't explicitly exclude any columns)
				only_rows, 		# list of integers (indices). Can be None (include all rows)
				omit_rows,		# list of integers (indices). Can be None (don't explicitly exclude any rows)
				only_samples,
				omit_samples,
				only_observations,
				omit_observations,
				only_samples_by_metadata,
				omit_samples_by_metadata,
				only_observations_by_metadata,
				omit_observations_by_metadata,
				case_insensitive_names,
				verbose,
				verbose_prefix):
	if(is_biom_file(file)):
		# read .biom table
		if(not HAS_BIOME_MODULE):
			print("%sERROR: The input table '%s' appears to be in BIOM format, but the BIOM module is not available"%(verbose_prefix,file))
			sys.exit(1)
		try:
			biom_table = biom.load_table(file);
			#biom_table = biom.table.Table.from_json(json.load(open(file,'rU')))
			#biom_table = biom.table.Table.from_hdf5(h5py.File(file,'r'))
		except Exception as e:
			if e.message.endswith('does not appear to be a BIOM file!'):
				print("%sERROR loading BIOM table '%s'.\n%sPerhaps not really in BIOM format, or perhaps the HDF5 or JSON module is not installed?"%(verbose_prefix,file,verbose_prefix))
				sys.exit(1);
			else:
				print("%sERROR loading biom table: %s. Perhaps not in BIOM format?"%(verbose_prefix,e.message))
				sys.exit(1);

		column_names = list(biom_table.ids('sample'))
		row_names = list(biom_table.ids('observation'))
		sample_metadata_for_filtering = (None if (only_samples_by_metadata=="" and omit_samples_by_metadata=="") else [biom_table.metadata(id=id, axis='sample') for id in column_names])
		observation_metadata_for_filtering = (None if (only_observations_by_metadata=="" and omit_observations_by_metadata=="") else [biom_table.metadata(id=id, axis='observation') for id in row_names])
		found_row_names = True
		found_column_names = True
		summary = "Original BIOM table contained %d rows & %d columns"%(len(row_names),len(column_names))
		
		# filter columns
		keep_columns = filter_index_list(len(column_names), only_columns, omit_columns)
		keep_columns = filter_by_name_and_metadata(	column_names,
													sample_metadata_for_filtering,
													keep_columns,
													only_samples,
													omit_samples,
													only_samples_by_metadata,
													omit_samples_by_metadata,
													(not case_insensitive_names))
		column_names = [column_names[c] for c in keep_columns]
		
		# filter rows
		keep_rows = filter_index_list(len(row_names), only_rows, omit_rows)
		keep_rows = filter_by_name_and_metadata(	row_names,
													observation_metadata_for_filtering,
													keep_rows,
													only_observations,
													omit_observations,
													only_observations_by_metadata,
													omit_observations_by_metadata,
													(not case_insensitive_names))
		row_names = [row_names[r] for r in keep_rows]

		# extract filtered BIOM table
		table = [[x[c] for c in keep_columns] for x in biom_table.iter_data(axis='observation')]
		delimiter = (delimiter if delimiter!="" else "\t")
		header_lines = ("# .biom table '%s'\n#Name" % file)+delimiter + delimiter.join(biom_table.ids('sample')) + "\n"		
		table = [table[r] for r in keep_rows]
		summary += "\nAfter filtering rows & columns based on names and metadata, obtained a table comprising %d rows & %d columns"%(len(keep_rows),len(keep_columns))

		
		# extract column & row metadata in original format BIOM (i.e. as list of dictionaries)
		row_metadata 	= [biom_table.metadata(id=row_name, axis='observation') for row_name in row_names]
		column_metadata = [biom_table.metadata(id=column_name, axis='sample') for column_name in column_names]
		
		# convert table to strings for uniform handling later on
		table = [["%.10g"%table[r][c] for c in range(len(table[r]))] for r in range(len(table))]
		
		if(verbose): print(verbose_prefix+"Read %d rows and %d columns from .biom file '%s'" % (len(row_names),len(column_names),file))
			
	else:
		# read classical table file
		row_metadata = None;
		column_metadata = None;
		linecount=0
		rowcount=0
		header_lines=[]
		last_comment_line=None
		all_column_names = None # not just the focal ones
		max_column_needed = (0 if (only_columns is None) else max(only_columns));
		table=[]
		if(not (only_rows is None)): only_rows = set(only_rows);
		row_names = []
		row_names_column = None
		column_names = None
		keep_columns = None
		found_row_names = False
		found_column_names = False
		with (gzip.open(file,'rt') if file.endswith(".gz") else open(file,'rt')) as fin:
			# read table line by line
			for line in fin:
				linecount+=1
				if(line[-1]=="\n"): line = line[0:-1] # remove trailing newline char
				parts = split_comments(line, comment_prefix)
				if((rowcount==0) and (len(parts)>1) and (parts[0].strip()=='')): 
					# pure comment line in header
					last_comment_line = parts[1];
					last_comment_line_count = linecount;
					header_lines.append(parts[1]);
				line = parts[0] # ignore comments
				if(line == ""): continue
			
				# this line seems to be a data line
				# see if we should interpret it (or the last comment line) as a header
				if((all_column_names is None) and (rowcount==0) and (not (last_comment_line is None)) and (column_names_are_in=='last_comment_line')):
					# interpret last comment line as header
					all_column_names = (last_comment_line.split() if (delimiter=="") else last_comment_line.split(delimiter))
					all_column_names = [name.strip() for name in all_column_names] # remove flanking whitespace
					header_lines = header_lines[0:-1]; # omit last comment line from file header lines
					if(len(all_column_names)<=max_column_needed):
						print(verbose_prefix+"ERROR: Need at least %d columns, but only found %d in header line %d" % (max_column_needed+1,len(all_column_names),last_comment_line_count))
						sys.exit(1)
					else:
						found_column_names = True;

				elif((all_column_names is None) and (rowcount==0) and (column_names_are_in=='first_data_line')):
					# interpret data line as header
					all_column_names = (line.split() if (delimiter=="") else line.split(delimiter))
					all_column_names = [name.strip() for name in all_column_names] # remove flanking whitespace
					if(len(all_column_names)<=max_column_needed):
						print(verbose_prefix+"ERROR: Need at least %d columns, but only found %d in line %d" % (max_column_needed+1,len(all_column_names),linecount))
						sys.exit(1)
					else:
						found_column_names = True
						continue;
				
				elif(all_column_names is None):
					# set column names to column indices
					all_column_names = [str(c) for c in range(len(line.split() if (delimiter=="") else line.split(delimiter)))]
			
				rowcount+=1
			
				parts = (line.split(delimiter) if (delimiter!="") else line.split())
				if(len(parts)<=max_column_needed):
					print(verbose_prefix+"ERROR: Need at least %d columns, but only found %d in line %d" % (max_column_needed+1,len(parts),linecount))
					sys.exit(1)
						
				# figure out where row names are
				if(row_names_column is None):
					if(row_names_are_in_column==""): 
						# row names are nowhere
						row_names_column = -1;
					elif(column_names_are_in!='none'):
						# columns have actual names, so search for the right column by name
						if(case_insensitive_names): row_names_are_in_column = row_names_are_in_column.lower()
						row_names_column = next((c for c in range(len(all_column_names)) if (all_column_names[c].lower()==row_names_are_in_column if case_insensitive_names else all_column_names[c]==row_names_are_in_column)), -1);
						if(row_names_column<0): 
							print(verbose_prefix+"ERROR: Unknown column '%s' specified for row names" % row_names_are_in_column)
							sys.exit(1)
						else: found_row_names = True;
					else:
						# column names are just column indices
						row_names_column = int(row_names_are_in_column);
						if(row_names_column>=len(parts)):
							print(verbose_prefix+"ERROR: Column index for row names exceeds number of available columns (%d)" % (len(parts)))
							sys.exit(1)
						else: found_row_names = True
						
						
				# figure out which columns to include
				if(keep_columns is None): 
					keep_columns = filter_index_list(len(all_column_names), only_columns, (omit_columns if (row_names_column<0) else [row_names_column]+([] if (omit_columns is None) else omit_columns)));
					column_names = [all_column_names[c] for c in keep_columns]
				
				# check if we should include this row
				if((only_rows is not None) and (not ((rowcount-1) in only_rows))): continue;
				if((omit_rows is not None) and ((rowcount-1) in omit_rows)): continue;
				
				# extract data, and possibly the row name, from this line
				if(row_names_column<0):
					# none of the columns contains row-names
					table.append([parts[c] for c in keep_columns])
					row_names.append(str(rowcount-1));
				else: 
					# one of the columns contains row names, so omit
					table.append([parts[c] for c in keep_columns])
					row_names.append(parts[row_names_column])
				if((len(table)>1) and (len(table[-1])!=len(table[-2]))): 
					print(verbose_prefix+"ERROR: Number of columns (%d) in line %d is inconsistent with previous data line" % (len(parts),linecount))
					sys.exit(1)

		if(len(table)==0):
			summary = "Original classical table is empty"
		else:
			summary = "Original classical table contained %d rows and %d columns\nAfter filtering rows & columns based on indices, obtained a table containing %d rows & %d columns"%(rowcount,len(all_column_names),len(table),len(column_names))
			if(verbose): print(verbose_prefix+"Loaded %d out of %d rows amongst %d lines, and %d columns, from file '%s'" % (len(table), rowcount, linecount, len(column_names), file))
			
	return 	table, row_names, column_names, found_row_names, found_column_names, header_lines, row_metadata, column_metadata, summary
	




def parse_group_name_and_metadata(line):
	name, metadata_unparsed = split_at_first_whitespace(line)
	if(metadata_unparsed==""): return name, None;
	metadata_entries = metadata_unparsed.split(';')
	NM = len(metadata_entries)
	metadata_keys 	= [None,]*NM
	metadata_values = [None,]*NM
	for m,metadata_entry in enumerate(metadata_entries):
		parts = metadata_entry.split(':',1)
		metadata_keys[m] = parts[0].strip()
		if(len(parts)==1):
			metadata_values[m] = "";
		else:
			metadata_values[m] = parts[1].split(',')
			if(len(metadata_values[m])==1): metadata_values[m] = metadata_values[m][0].strip();
			else: metadata_values[m] = [value.strip() for value in metadata_values[m]]
	metadata = {metadata_keys[m]:metadata_values[m] for m in range(NM)}
	return name, metadata




def read_groups_from_list(groups_list,allow_duplicate_groups,verbose_prefix):
	group_members 		= [];
	group_names 		= [];
	all_members			= [];
	group_name2index	= {};
	if(groups_list!=""):
		groups_str = re.split(r'(?<!\\),', groups_list)				# split at commas but exclude escaped commas
		groups_str = [s.replace('\,',',') for s in groups_str]		# un-escape escaped commas
		for group_str in groups_str:
			parts = re.split(r'(?<!\\):', group_str)				# split at colons but exclude escaped colons
			parts = [s.replace('\:',':') for s in parts]			# un-escape escaped colons
			if(parts[0] is ""):
				print(verbose_prefix+"ERROR: Empty group name found in groups list%s" % (" (with members '%s', ..)"%parts[1] if (len(parts)>0) else ""))
				sys.exit(1);
			next_member_index = len(all_members)
			all_members.extend((p,p) for p in parts[1:])
			group_name = parts[0]
			if(group_name in group_name2index):
				if(allow_duplicate_groups):
					# merge with existing group
					group_index = group_name2index[group_name]
				else:
					print("%sERROR: Duplicate group names: The group '%s' was defined multiple times. Consider allowing duplicate groups."%(verbose_prefix,group_name))
					sys.exit(1)
			else:
				# this is a new group
				group_names.append(group_name)
				group_index = len(group_names)-1
				group_name2index[group_name] = group_index
				group_members.append([])
			# add members to group
			group_members[group_index].extend(list(range(next_member_index,next_member_index+len(parts)-1)))
	group_metadata = [None,] * len(group_names)
	return group_members, group_names, group_metadata, all_members;
	
	


# returns group_names[], group_members[][], all_members[], all_unique_members[]
# 	Each group_members[g][] contains member indices (int) refering to all_members[], as well as group set operations (lists, each containing two ints)
# 	Each set operation is a list [operation_type, referenced_group], where operation_type (int) is interpreted according to set_operations_keywords[]
#	all_members[] is a redundant list of all defined members. Each entry is a 2-tuple (name/wildcard/regex, original_line)
#	all_unique_members[] is a set of strings representing unique member_names/wildcards/regex
def read_groups(groups_file,
				additional_list_front, 	# additional comma-separated list of groups, each of which is a colon-separated list of [group name and] member names
				additional_list_back, 	# similar to additional_list_front, but these groups are appended to the list (instead of prepended)
				delim, 
				comment_prefix, 
				no_group_names_in_file,
				single_line_groups,
				allow_group_set_operations,
				allow_duplicate_groups,# if true, duplicate definitions of groups are simply merged; otherwise, an error is thrown if a group is defined multiple times
				verbose,
				verbose_prefix):
		
	all_members			= [];	# all_members[m] is a tuple (member_name/wildcard/regex, original_line)
	group_members		= [];	# group_members[g][] is a list containing indices (pointing to all_members[]) or pairs [operation, referenced_group] corresponding to group operations
	group_names			= [];
	group_metadata		= [];
	group_name2index	= {};
	all_unique_members 	= set();
	# read groups from file
	if(groups_file!=""):
		with (gzip.open(groups_file,'rt') if (os.path.splitext(groups_file)[1]==".gz") else open(groups_file,'rt')) as fin:
			linecount=0
			group_name=None;
			if(not single_line_groups):
				for original_line in fin.readlines():
					linecount+=1
					if(original_line.strip() is ""):
						# expect new group (or end of file)
						group_name = None;
						continue;
					original_line = original_line.rstrip("\n")		
					line = (original_line if (comment_prefix is "") else original_line.split(comment_prefix))[0].strip() # ignore comments and strip preceding & trailing white space
					if((len(line)>0) and (line[0] in ['"','\'']) and (line[-1] in ['"','\''])): line = line[1:-1]; # if enclosed in quotes, extract quoted part
					if(line==""): continue;
					if(group_name is None):
						# create new group or find existing group
						if(no_group_names_in_file):
							# use a number for group name. Hence, groups are never duplicates
							group_members.append([]);
							group_name = str(len(group_names));
							group_names.append(group_name);
							group_metadata.append(None);
							group_index = len(group_names)-1;
							group_name2index[group_name] = group_index;
						else:
							# the first line in the group is its name
							group_name, this_group_metadata = parse_group_name_and_metadata(line);
							if(group_name in group_name2index):
								# group has been redefined before
								if(allow_duplicate_groups):
									# merge with existing group
									group_index = group_name2index[group_name]
								else:
									print("%sERROR: Duplicate group names: The following group was defined multiple times: '%s'. Consider allowing duplicate groups."%(verbose_prefix,group_name))
									sys.exit(1)
							else:
								# this is a new group
								group_members.append([]);
								group_names.append(group_name);
								group_metadata.append(this_group_metadata);
								group_index = len(group_names)-1;
								group_name2index[group_name] = group_index;
							continue;
					# line is either a member (name/wildcard/regex) to be added, or a group set operation
					if(allow_group_set_operations and (whichPrefix(line,set_operations_keywords)>=0)):
						# member is a set operation, convert into a number and extract target of set operation
						operation = whichPrefix(line,set_operations_keywords)
						line = line[len(set_operations_keywords[operation]):]
						if(line in group_name2index):
							referenced_group = group_name2index[line]
							if(referenced_group==group_index):
								print(verbose_prefix+"ERROR: Self-operation encountered in group '"+group_name+"'")
								sys.exit(1)
							else:
								group_members[group_index].append([operation,referenced_group])
						else:
							print(verbose_prefix+"ERROR: Unknown group '"+line+"', referenced by set operation in group '"+group_name+"'. Referenced groups need to be defined beforehand.")
							sys.exit(1)
					else:
						# add new member to current group
						all_members.append((line, original_line))
						group_members[group_index].append(len(all_members)-1);
						all_unique_members.add(line);
					
			else:
				for line in fin.readlines():
					linecount+=1
					line = (line if (comment_prefix is "") else line.split(comment_prefix))[0].rstrip() # ignore comments and strip right white space
					if(line == ""): continue
					parts = list(filter(len, (line.split(delim) if (delim!="") else line.split())))
					if(no_group_names_in_file):
						group_name = str(len(group_names))
						group_names.append(group_name)
						group_name2index[group_name] = len(group_names)-1
						members = parts
					else:
						group_name = parts[0]
						members = parts[1:]
						if(group_name in group_name2index):
							# group has been redefined before
							if(allow_duplicate_groups):
								# merge with existing group
								group_index = group_name2index[group_name]
							else:
								print("%sERROR: Duplicate group names: The group '%s' was defined multiple times. Consider allowing duplicate groups."%(verbose_prefix,group_name))
								sys.exit(1)
						else:
							# this is a new group
							group_names.append(group_name)
							group_index = len(group_names)-1
							group_name2index[group_name] = group_index
							group_members.append([]);
					for member in members:
						# member is either a name/wildcard/regex to be added, or a group set operation
						if(allow_group_set_operations and (whichPrefix(member,set_operations_keywords)>=0)):
							# member is a set operation
							operation = whichPrefix(member,set_operations_keywords)
							member = member[len(set_operations_keywords[operation]):]
							if(member in group_name2index):
								referenced_group = group_name2index[member]
								if(referenced_group==len(group_names)-1):
									print(verbose_prefix+"ERROR: Self-operation encountered in group '"+group_name+"'")
									sys.exit(1)
								else:
									group_members[group_index].append([operation,referenced_group])
							else:
								print(verbose_prefix+"ERROR: Unknown group '"+member+"', referenced by set operation in group '"+group_name+"'. Referenced groups need to be defined beforehand.")
								sys.exit(1)
						else:
							# add new member to current group
							all_members.append((member, member))
							group_members[group_index].append(len(all_members)-1);
							all_unique_members.add(member);
				group_metadata = [None,]*len(group_names); # no group metadata in single-line group format
			if(verbose): print(verbose_prefix+"Read %d lines from groups file, found %d groups with %d members (%d unique members)" % (linecount, len(group_members), len(all_members), len(all_unique_members)))
		
	# prepend/append additional groups defined in lists
	# Note: Currently groups that are defined multiple times across various sources (i.e. in at least two of additional_list_front, additional_list_back and groups_file) will not be merged (i.e. the option allow_duplicate_groups is only implemented for within any single source of groups)
	additional_group_members_front, additional_group_names_front, additional_group_metadata_front, additional_all_members_front = read_groups_from_list(additional_list_front, allow_duplicate_groups, verbose_prefix);
	additional_group_members_back, additional_group_names_back, additional_group_metadata_back, additional_all_members_back = read_groups_from_list(additional_list_back, allow_duplicate_groups, verbose_prefix);
	NAGF = len(additional_group_names_front);
	temp_NAM = len(all_members)
	group_members 	= [[(temp_NAM + m) for m in members] for members in additional_group_members_front] + group_members + [[(temp_NAM + len(additional_all_members_front) + m) for m in members] for members in additional_group_members_back]
	all_members.extend(additional_all_members_front)
	all_members.extend(additional_all_members_back)
	group_members 	= [[(member if isinstance(member, int) else [member[0],member[1]+NAGF]) for member in g] for g in group_members] # shift indices of referenced groups in set operations, to account for prefixed additional groups
	group_names   	= additional_group_names_front + group_names + additional_group_names_back
	group_metadata 	= additional_group_metadata_front + group_metadata + additional_group_metadata_back
	all_unique_members.update(member[0] for member in additional_all_members_front)
	all_unique_members.update(member[0] for member in additional_all_members_back)
	group_metadata = [(gm if (gm is not None) else {}) for gm in group_metadata]
		
	return group_members, group_names, group_metadata, all_members, all_unique_members
	
	
	
# Calculate effective number of members per group, taking into account set operations
# This is only approximate, as it is hard to define the number of member for a group that includes set operations.	
def calculate_effective_number_of_members_per_group(group_members,N_all_members):
	NG = len(group_members)
	group_contains_member = numpy.zeros([NG,N_all_members], dtype=bool)
	for g,members in enumerate(group_members):
		for member in members:
			if(isinstance(member,int)):
				group_contains_member[g,member] = True
			else:
				operation 	 = member[0]
				target_group = member[1]
				if(operation==0): 	group_contains_member[g,group_contains_member[target_group,:]] = True 	# add (set union)
				elif(operation==1): group_contains_member[g,group_contains_member[target_group,:]] = False 	# subtract (set difference)
				else: 				group_contains_member[g,numpy.logical_or(numpy.logical_not(group_contains_member[g,:]),numpy.logical_not(group_contains_member[target_group,:]))] = False		# set intersection
	effective_number_of_members_per_group = list(numpy.sum(group_contains_member,axis=1,dtype=float))
	return effective_number_of_members_per_group
	
	
	
	
	
	
def find_matches_to_words_expression(expression, candidates, valid_word_symbols):
	words = list(filter(len,expression.split('*')));
	matches = [];
	for ci in range(len(candidates)):
		candidate = candidates[ci];
		LC = len(candidate);
		next_start = 0;
		found_word = True;
		for word in words:
			found_word = False;
			while(True):
				pos = candidate.find(word, next_start);
				if(pos<0):
					found_word = False;
					break; # word not found, so drop this candidate
				next_start = pos+len(word);
				if((pos>0) and (candidate[pos-1].isalnum() or (candidate[pos-1] in valid_word_symbols))):
					continue;	# instance does not seem to be a complete word
				elif((next_start<LC) and (candidate[next_start].isalnum() or (candidate[next_start] in valid_word_symbols))):
					continue;	# instance does not seem to be a complete word
				else:
					found_word = True;
					break;
			if(found_word): continue;
			else: break;
		if(found_word): matches.append(ci);
	return matches;
			
			
	
# record_labels is a list of record labels (e.g. taxonomy)
# group_members is a list of lists, each of which containing the members of a group (or group set operations)
# all_members[] is a list of all members somehow referred to in various groups. Each 
#
# group members can be:
#	a. an integer index pointing to an entry in all_members[]
#	b. a group set operation defined as a list [operation_type, referenced_group]
#
# Each entry in all_members[] is a 2-tuple (member, original_line), where member can be:
#	a. A name, i.e. simple string matches modulo flanking whitespace, e.g. "norBC" (iff group_members_defined_as=='match')
#      For example, 'norBC' will match '  norBC' but not 'abba norBC' or 'abbanorBC'
#	b. A shell wildcard expression, e.g. "*cox*"  (iff group_members_defined_as=='wildcards')
#	c. A regular expressions, e.g. ".*nor[BC]"  (iff group_members_defined_as=='regex')
#	d. A sequence of complete words separated by *, where * can represent additional intermediate words wildcard (iff group_members_defined_as=='words'). 
#      For example, '*Proteobacteria*europaea*' is matched by 'Proteobacteria;europaea' and 'Bacteria:Proteobacteria:Nitrosomonas:europaea', 
#      but not by 'Proteobacteria:Nitroeuropaea' or 'Proteobacteria:1europaea' or 'Proteobacteriaeuropaea'.
def assign_records_to_groups(	group_names, 
								group_members,
								all_members,
								record_labels, 
								case_insensitive, 
								valid_word_symbols, # only relevant if group_members_defined_as=='words'
								group_members_defined_as): # one of 'match', 'regex', 'words', 'wildcards', 
	NG = len(group_members);
	NR = len(record_labels);
	if(case_insensitive): record_labels = [name.lower() for name in record_labels]
	record_labels = [name.strip() for name in record_labels] # remove flanking whitespace
	
	if(group_members_defined_as=='match'):
		# create records labels hash for faster matching
		# label2records[label] will be a list of record indices with the given label
		label2records = {}
		for r in range(NR):
			if(record_labels[r] in label2records):
				label2records[record_labels[r]].append(r)
			else:
				label2records[record_labels[r]] = [r]
	
	group_to_records 		= [set() for g in range(NG)];
	group_members_used 		= [[] for g in range(NG)]	# group_members_used[g] has a similar structure to group_members[g], but only includes members that were actually used (or representing group set operations)
	group_members_unused 	= [[] for g in range(NG)]	# group_members_unused[g] has a similar structure to group_members[g], but only includes members that were not used
	for g in range(NG):
		for m in group_members[g]:
			if(isinstance(m, int)):
				member = all_members[m][0]
				# member is name/wildcard/regex
				temp_previous_NR = len(group_to_records[g])
				if(group_members_defined_as=='match'):
					member = member.strip();
					if(case_insensitive): member = member.lower()
					if(member in label2records):
						added = label2records[member]
					else:
						added = []
				if(group_members_defined_as=='regex'):
					regex = re.compile(member,(re.IGNORECASE if case_insensitive else 0))
					added = [r for r in range(NR) if (regex.search(record_labels[r]) is not None)]
				elif(group_members_defined_as=='wildcards'):
					if(case_insensitive): member = member.lower()
					added = [r for r in range(NR) if fnmatch.fnmatchcase(record_labels[r],member)]
				elif(group_members_defined_as=='words'):
					if(case_insensitive): member = member.lower()
					added = find_matches_to_words_expression(member, record_labels, valid_word_symbols)
				group_to_records[g].update(added)
				if(len(added)>0):
					# this member caused the inclusion of more records in the group, so add it to group_members_used[g][]
					group_members_used[g].append(m)
				else:
					# this member was not used, so add it to group_members_unused[g][]
					group_members_unused[g].append(m)
			else:
				# member is a group set operation, i.e. a list [operation_type, referenced_group]
				operation = m[0];
				group_members_used[g].append(m);
				if(operation==0): group_to_records[g].update(group_to_records[m[1]]) 					# add (set union)
				elif(operation==1): group_to_records[g].difference_update(group_to_records[m[1]]) 		# subtract (set difference)
				else: group_to_records[g].intersection_update(group_to_records[m[1]]) 					# set intersection
				
	# figure out leftover records
	is_leftover = [True]*NR;
	for g in range(NG):
		for r in group_to_records[g]: 
			is_leftover[r] = False;
	leftover_records = set(r for r in range(NR) if is_leftover[r])
	
	return group_to_records, leftover_records, group_members_used, group_members_unused
	
	


def consolidate_categorial(values, non_consensus_value):
	if(len(set(values))==1): return values[0];
	else: return non_consensus_value;
	
	
def XOR(a,b):
	return (a and (not b)) or (b and (not a));
	
	
def get_date_time():
	return time.strftime("%Y.%m.%d") + " " + time.strftime("%H:%M:%S")



def get_shell_command():
	arguments = sys.argv;
	arguments = [("'"+a+"'" if (a=="" or re.search(r"[\s*#:;]", a)) else ("'"+'\\\\t'+"'" if (a=='\\t') else a)) for a in arguments];
	return ' '.join(arguments);



def arbitrary_metadata_values_to_record_name(values, group_members_defined_as):
	if(group_members_defined_as=='words'):
		return ["" if (value is None) else ('*'.join(value) if isinstance(value,list) else str(value)) for value in values];
	else:
		return ["" if (value is None) else (''.join(value) if isinstance(value,list) else str(value)) for value in values];



def find_duplicates_in_list(L):
	L_unique = set()
	duplicates = set()
	for l in L:
		if(l in L_unique): duplicates.add(l)
		else: L_unique.add(l)
	return list(duplicates)



def normalize_table(table,normalization):
	if(normalization=="columns"):
		sums = numpy.nansum(table, axis=0, keepdims=True);
		sums[sums==0] = 1.0
		table = table/sums
	elif(normalization=="rows"):
		sums = numpy.nansum(table, axis=1, keepdims=True);
		sums[sums==0] = 1.0
		table = table/sums;
	elif(normalization=="none"):
		pass;
	return table;




def  save_subtable(	file,
					table_id,					# e.g. 'focal_group_01' or 'anammox_subtable'
					asBIOM,						# if false, classical (tabular) format will be used instead of BIOM
					table,						# numpy array of size NR * NC
					row_names,					# if asBIOM==True, then these will be observation IDs
					column_names,				# if asBIOM==true, then these will be sample IDs
					header_lines,				# additional header lines to include from original table file. May be None. Only relevant if asBIOM==False
					row_metadata, 				# may be None, only relevant for BIOM output files
					column_metadata,			# may be None, only relevant for BIOM output files
					normalize,					# one of 'none', 'rows', 'columns'
					include_summary_comments,
					column_names_in,			# can be 'none', 'last_comment_line' or 'first_data_line'
					include_numbers_in_column_header,
					delimiter,
					comment_prefix,
					verbose,
					verbose_prefix):
	NR = table.shape[0]
	NC = table.shape[1]
	
	# normalize output table if needed
	table = normalize_table(table, normalize);

	# write table
	if(asBIOM):
		# BIOM format
		biom_table = biom.table.Table(	table, 
										observation_ids=row_names, 
										sample_ids=column_names, 
										observation_metadata=row_metadata, 
										sample_metadata=column_metadata,
										table_id=table_id);
		save_biom_table(biom_table, file, "BIOM")

	else:
		# Classical (tabular) format
		with (gzip.open(file,'wt') if file.lower().endswith(".gz") else open(file,'wt')) as fout:
			if(include_summary_comments): fout.write("%s %s: subtable of '%s'\n%s Generated on: %s\n%s Used command:\n%s   %s\n%s\n%s Summary: %d rows, %d columns\n" % (comment_prefix,table_id, args.input_table,comment_prefix,get_date_time(),comment_prefix,comment_prefix,get_shell_command(),comment_prefix,comment_prefix,len(row_names),len(column_names)))
			if(normalize!='none'): fout.write("%s %s are normalized to unit sum\n"%(comment_prefix,normalize));
			fout.write("%s\n"%(comment_prefix));
			if(header_lines is not None): fout.write("%s Original (uncollapsed) table header:\n%s%s\n"%(comment_prefix,comment_prefix,("\n"+comment_prefix).join(header_lines)))
			delimiter = (delimiter if delimiter!="" else "\t")
			fout.write((comment_prefix if column_names_in=='last_comment_line' else '') + 'record' + delimiter + delimiter.join([(str(c+1)+":" if include_numbers_in_column_header else '') + column_names[c] for c in range(NC)])+"\n")
			for r in range(NR):
				fout.write(row_names[r]+delimiter+delimiter.join([object2str(table[r,c]) for c in range(NC)])+"\n")



def is_cultured_taxon(name):
	name = name.lower();
	if('uncultured' in name): return False;
	if('metagenome' in name): return False;
	if('unknown' in name): return False;
	if('unidentified' in name): return False;
	if('other' in name): return False;
	return True;



def get_jaccard_index(set1, set2):
	if((len(set1)==0) and (len(set2)==0)): return 1.0;
	N = len(set1.intersection(set2));
	return N/float(len(set1) + len(set2) - N); 




##########################################################
# multi-layer sorting



def get_order_from_tree(IDs, file_path):
	if(not os.path.exists(file_path)):
		raise IOError("Missing tree file, expected at: %s" % (file_path))
	tree 		= Phylo.read(file_path, 'newick')
	ID2index 	= {ID:i for i,ID in enumerate(IDs)}
	order 		= [None]*len(IDs)
	# for details on the various tree functions, see: http://biopython.org/DIST/docs/api/Bio.Phylo.BaseTree-pysrc.html
	leafs		= tree.find_clades(target=None,terminal=True,order='preorder')
	next_item	= 0
	for leaf in leafs:
		if(leaf.name in ID2index):
			index = ID2index[leaf.name]
			if(index<0):
				raise RuntimeError("Leaf ID '%s' found multiple times in tree file"%(leaf.name))
			order[next_item] 	= index
			ID2index[leaf.name] = -1
			next_item += 1
	
	if(next_item<len(IDs)):
		# some IDs were not found in the tree file. Stack all leftover IDs at the end
		for ID,index in ID2index.items():
			if(index>=0):
				order[next_item] = index
				next_item += 1
	return order



def atoi(text):
	try:
		return (int(text) if text.isdigit() else text)
	except Exception:
		return NaN
		
		
def natural_sorting_for_value_index_pair(value_index_pair):
	if(isinstance(value_index_pair[0],str)):
		return [atoi(c) for c in re.split('(\d+)', value_index_pair[0])]
	else:
		return value_index_pair[0]


def natural_sorting(value):
	if(isinstance(value,str)):
		return [atoi(c) for c in re.split('(\d+)', value)]
	else:
		return value


# turns a list of strings into a list of numbers, provided that all entries are valid numbers
# also returns a boolean, indicating whether a conversion took place
def list2num_if_sensible(L):
	if(next((x for x in L if (not is_number(x))), -1)>=0): return L,False;
	else: return [float(x) for x in L], True;



# sort indices[block_start],indices[block_start+1],...  according to values[indices[block_start]], values[indices[block_start+1]], ...
# indices[] points to entries in values[], so that roughly values[indices[]] is sorted
def sort_block_by_indices_in_situ(	indices, 
									values,
									block_start, 
									block_end, 
									sorting_direction, 		# can be 'f' (forward),'r' (reverse), 'n' (natural), 'nr' (natural-reverse) or 'c' (custom)
									custom_direction_list):	# (list of strings) only relevant if sorting_direction=='c'
	if(block_end<=block_start): return indices; # nothing to do
	local_values = [values[indices[i]] for i in range(block_start,block_end+1)]		
	if(sorting_direction=='c'):
		# customize order relationship
		local_values = [str(v).lower() for v in local_values]
		anchors = {anchor.lower():a for a,anchor in enumerate(custom_direction_list)}
		local_values = [(anchors[value] if (value in anchors) else (len(anchors)+v)) for v,value in enumerate(local_values)]
	else:
		local_values = [(v.lower() if isinstance(v,str) else v) for v in local_values]
	if(sorting_direction in ['n','nr']):
		local_values, numeric = list2num_if_sensible(local_values)
		if(numeric):
			temp = sorted([(v,indices[i+block_start]) for i,v in enumerate(local_values)], reverse=(sorting_direction=='nr'))
		else:
			temp = sorted([(v,indices[i+block_start]) for i,v in enumerate(local_values)], key=natural_sorting_for_value_index_pair, reverse=(sorting_direction=='nr'))
	else:
		temp = sorted([(v,indices[i+block_start]) for i,v in enumerate(local_values)], reverse=(sorting_direction=='r'))
	for i,t in enumerate(temp):
		indices[i+block_start] = t[1]
	return indices




def define_subgroups_by_attribute(order,old_group_assignments,attributes):
	N = len(old_group_assignments)
	new_group_asignments = [None,]*N
	new_group_asignments[order[0]] = 0;
	for i in range(1,N):
		new_group_asignments[order[i]] = new_group_asignments[order[i-1]]
		if((old_group_assignments[order[i-1]]!=old_group_assignments[order[i]]) or (attributes[order[i-1]]!=attributes[order[i]])):
			new_group_asignments[order[i]] += 1;
	return new_group_asignments


def get_blocks_from_group_assignments(order,group_assignments):
	block_starts = []
	block_ends = []
	block_start = 0
	N = len(group_assignments)
	for n in range(1,N):
		if(group_assignments[order[n-1]]!=group_assignments[order[n]]):
			# encountered a new group
			block_starts.append(block_start)
			block_ends.append(n-1)
			block_start = n
	# add last group
	block_starts.append(block_start)
	block_ends.append(N-1)
	return block_starts, block_ends
	
	


def sort_multilayered(	order_specification,	# (string) semicolon-separated list of <data_type>:<direction> pairs, where <data_type> can be 'label','given','average','non_zeros', 'data':<data_name>, or 'metadata':<metadata_name>  and <direction> can be 'f','r' or 'c':<comma-separated_list of values>
						IDs,
						labels,
						table,
						data_names,				# list of data names (column names when sorting rows, row names when sorting columns)
						all_metadata,
						table_axis,				# 0 if we're sorting rows, 1 if we're sorting columns
						verbose_prefix):
	N = len(labels)
	if(order_specification==""): return None;
	
	# parse order specification
	layers = list(filter(len,order_specification.split(';')))
	NL = len(layers)
	layer_types = [None,]*NL
	values = [None,]*NL
	directions = [None,]*NL
	custom_direction_lists = [None,]*NL
	for l in range(NL):
		parts = layers[l].split(':')
		next_p = 0;
		layer_types[l] = parts[next_p]
		if(len(parts)<2):
			print("%sERROR in order specification '%s': Missing ':' in layer '%s'"%(verbose_prefix,order_specification,layers[l]))
			sys.exit(1)
		if(layer_types[l]=='metadata'):
			next_p += 1;
			metadata_name = parts[next_p]
			values[l] = [("" if ((all_metadata is None) or (all_metadata[n] is None) or (metadata_name not in all_metadata[n])) else str(all_metadata[n][metadata_name])).lower() for n in range(N)]
			values[l], dummy = list2num_if_sensible(values[l])
		elif(layer_types[l]=='data'):
			next_p += 1;
			data_name = parts[next_p]
			data_index = next((d for d in range(len(data_names)) if (data_names[d]==data_name)), -1);
			if(data_index<0): 
				print("%sERROR in order specification '%s': Unknown data '%s' specified in layer '%s'"%(verbose_prefix,order_specification,data_name,layers[l]))
				sys.exit(1)
			values[l] = ([table[n,data_index] for n in range(N)] if (table_axis==0) else [table[data_index,n] for n in range(N)])
			values[l], dummy = list2num_if_sensible(values[l])
		elif(layer_types[l]=='given'):
			values[l] = list(range(N))
		elif(layer_types[l]=='label'):
			values[l] = labels
		elif(layer_types[l]=='average_value'):
			values[l] = numpy.mean(table,axis=(1-table_axis))
		elif(layer_types[l]=='average_absolute_value'):
			values[l] = numpy.nanmean(numpy.abs(table),axis=(1-table_axis))
		elif(layer_types[l]=='non_zeros'):
			values[l] = numpy.sum(table!=0,axis=(1-table_axis))
		else:
			print("%sERROR: Unknown layer type '%s' for sorting"%(verbose_prefix,layer_types[l]))
			sys.exit(1)
		next_p += 1
		if(next_p>=len(parts)):
			print("%sERROR in order specification '%s': Missing direction in layer '%s'"%(verbose_prefix,order_specification,layers[l]))
			sys.exit(1)			
		directions[l] = parts[next_p]
		if(directions[l] not in ['f','r','c','n','nr']):
			print("%sERROR in order specification '%s': Invalid direction '%s' in layer '%s'"%(verbose_prefix,order_specification,directions[l],layers[l]))
			sys.exit(1)
		if(directions[l]=='c'):
			next_p += 1
			if(len(parts)<=next_p):
				print("%sERROR in order specification '%s': Missing custom direction list in layer '%s'"%(verbose_prefix,order_specification,layers[l]))
				sys.exit(1)
			custom_direction_lists[l] = parts[next_p].split(',')
		
		
	# sort values according one layer at a time
	group_assignments = [0,]*N
	order = list(range(N))
	for l in range(NL):
		block_starts, block_ends = get_blocks_from_group_assignments(order,group_assignments)
		for b in range(len(block_starts)):
			order = sort_block_by_indices_in_situ(	order, 
													values[l],
													block_starts[b], 
													block_ends[b], 
													directions[l],
													custom_direction_lists[l])
		# renew (finer) group assignments based on new (smaller) groups
		group_assignments = define_subgroups_by_attribute(order,group_assignments,values[l])
	return order




	
#####################################################
# Checking output files


def check_output_file(file,force,verbose,verbose_prefix):
	if(os.path.exists(file)):
		if(force): 
			if(os.path.isdir(file)): shutil.rmtree(file)
			else: os.remove(file)
			if(verbose): print("%sNote: Replacing output file '%s'" % (verbose_prefix,file))
		else:
			print("%sERROR: Output file '%s' already exists\n%s       Cowardly refusing to continue.\n%s       Use --force to ignore this message" % (verbose_prefix,file,verbose_prefix,verbose_prefix))
			sys.exit(1)
	else:
		parent = os.path.dirname(file);
		if((not parent is '') and (not os.path.exists(parent))): os.makedirs(parent)
	
	
def check_output_file_append(file,force,verbose,verbose_prefix):
	if(os.path.exists(file)):
		if(os.path.isdir(file)):
			if(force): 
				shutil.rmtree(file)
				if(verbose): print("%sNote: Replacing directory '%s' with an output file" % (verbose_prefix,file))
			else:
				print("%sERROR: Output file '%s' already exists as a directory\n%s       Cowardly refusing to continue.\n%s       Use --force to ignore this message" % (verbose_prefix,file,verbose_prefix,verbose_prefix))
				sys.exit(1)
	else:
		parent = os.path.dirname(file);
		if((not parent is '') and (not os.path.exists(parent))): os.makedirs(parent)
			
			
			
def check_output_dir(dir,force,verbose,verbose_prefix):
	if(os.path.isfile(dir)):
		if(args.force): 
			if(verbose): print("%sNote: Replacing file '%s' with output directory" % (verbose_prefix,dir))
			shutil.rmtree(dir)
			os.makedirs(dir)
		else:
			print("%sERROR: Output directory '%s' is a file.\n%s       Cowardly refusing to continue.\n%s       Use --force to avoid this message" % (verbose_prefix,dir,verbose_prefix,verbose_prefix))
			sys.exit(1)
	elif(not os.path.exists(dir)): os.makedirs(dir)

#####################################################



if __name__ == '__main__':
	# parse command line arguments
	parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter, description="Collapse table rows (or columns if --collapse_columns_instead_of_rows is set) down to a defined set of novel groups. Collapse is either by averaging or summing. Row names don't need to be unique. For example rows may be 'named' using OTU taxonomy and multiple OTUs may have the same taxonomy.\n\nA group's member list can include set operations (addition/subtraction) with previously defined groups. Addition (set union) is indicated by the prefix 'add_group:', subtraction (set difference) by the prefix 'subtract_group:', intersection by the prefix 'intersect_group:'. For example, a group 'mammal_pathogens' might include another group via 'add_group:human_pathogens'. If --no_group_names_in_file is set, groups can be referred to by number (starting at 0). Group set operations are only allowed for groups defined in the groups file (see option --groups_file).", epilog="Copyright (c) 2016, Stilianos Louca. Please read the license agreement in the header of this script.")
	parser.add_argument('-i','--input_table', action='store', required=True, help="Path to a classical table file (TSV, CSV or similar) or a BIOM file. Row names are stored in a column specified by --row_names_are_in_column, the rest of the columns contain numerical or string data. Non-numerical entries will be interpreted according to --non_numeric. Alternatively, a .biom observation table can be given, in which case observation IDs will be taken as 'row names'.");
	parser.add_argument('-g','--input_groups_file', action='store', default='', help="Path to a file defining the groups by which to collapse the table. The file should list one group member per line, with each group spanning several lines (if the flag --single_line_groups is not set) or one group per row and spanning multiple columns (if the flag --single_line_groups is set). If neither --single_line_groups nor --no_group_names_in_file are set, group names can be followed by whitespace and an optional semicolon-separated list of metadata entries of the format <key>:<value> or <key>:<value1>,<value2>,<value3> and so on. If --single_line_groups is set, then the first column contains the group name (unless --no_group_names_in_file is set), the rest of the columns contain the names of the original records (rows, or columns if --collapse_columns_instead_of_rows is set) that are to be grouped together. Groups can be overlapping, i.e. share members. Group member names can be wildcard expressions (shell style) if --group_members_defined_as is 'wildcards', or regular expressions (python style) if --group_members_defined_as is 'regex' or a list of complete words if --group_members_defined_as is 'words'. Leave blank (default) to not load any groups from file (also see options --groups_list_front and --groups_list_back).");
	parser.add_argument('-o','--out_collapsed', action='store', default="", type=str, help='Path to optional output collapsed table.');
	parser.add_argument('-r','--out_report', default='', help='Path to optional output report file. If provided, this file will list the rows (or columns) representing each group.');
	parser.add_argument('-l','--out_log', action='store', default='', help="Path to optional log file, to which summary comments are to be appended.");
	parser.add_argument('-s','--out_sub_tables_dir', action='store', default='', help="Path to optional directory, to which sub-tables shall be saved. Each subtable contains a subset of the original input table, with rows (or columns, if --collapse_columns_instead_of_rows is set) corresponding to a particular group. Leave this empty (default) to not save sub-tables.");
	parser.add_argument('--out_groups2records_table', default='', help='Path to optional output table listing group-record associations for the particular data set. By default, this table will list original records as rows, and groups as columns, with each cell being 1 (association) or 0 (no association), but see option --normalize_groups2records_table.');
	parser.add_argument('--out_groups2records_table_dense', default='', help='Path to optional classical (tabular) output table listing groups associated with each record. This table will list original records as rows, and groups assigned to each record as comma-separated lists. If row metadata is used for collapsing (instead of row names), that metadata is also listed in a separate column.');
	parser.add_argument('--out_collapsed_deconvoluted_table_records_vs_data', default='', help='Path to optional output BIOM table listing record-data table after all transformations. This BIOM table will have a similar structure to the original input_table (numerical part only), but entries will be adjusted such that the matrix product out_collapsed_deconvoluted_table_groups_vs_records * out_collapsed_deconvoluted_table_records_vs_sample  yields the collapsed output table.');
	parser.add_argument('--out_collapsed_deconvoluted_table_groups_vs_records', default='', help='Path to optional output BIOM table listing group-record association table. This BIOM table will such that the matrix product out_collapsed_deconvoluted_table_groups_vs_records * out_collapsed_deconvoluted_table_records_vs_data  yields the collapsed output table (group vs data).');
	parser.add_argument('--out_group_overlaps', default='', help='Path to optional output table listing group overlaps (Jaccard similarity index) in terms of shared records from the input table. An overlap of 1.0 betwee two groups means that the exact same records were assigned to both groups. An overlap of 0.0 means the two groups have no records in common and at least one of them is non-empty.');
	parser.add_argument('--out_group_definitions_used', default='', help='Path to optional output file listing group definitions relevant to the input data, i.e. members of each group that were actually matched to at least one record.');
	parser.add_argument('--out_group_definitions_unused', default='', help='Path to optional output file listing group definitions irrelevant to the input data, i.e. members of each group that were not matched to any record.');
	parser.add_argument('--out_all_members', default='', help='Path to optional output text file, for saving all unique members encountered in the loaded group definitions. This file will list one member per row.');
	parser.add_argument('--groups_list_front', action='store', default='', help="Comma-separated list of additional groups (in addition to groups in --groups_file). Each group begins with a non-empty group name, followed by a colon and a colon-separated list of member names. For example, 'nitrifiers:amo:nxr,denitrifiers:narGHIJ:norBC:nosZ' defines two groups, 'nitrifiers' and 'denitrifiers'. Leave blank (default) for no additional groups. To include a comma or colon in group or member names, escape them with a back-slash. Additional groups specified using --groups_list_front are prepended to the list of groups loaded from file.");
	parser.add_argument('--groups_list_back', action='store', default='', help="Similar to groups_list_front, with the differebce that groups are appended to the list of groups loaded from file.");
	parser.add_argument('-d','--row_names_are_in_column', default="", help="Column containing row names (as referred to in the groups unless --collapse_columns_instead_of_rows is set). If column names are available, this specifies a column by name, otherwise it specifies a column by index (starting at 0). If this is empty ('', default), then row names are set to row indices (starting at 0) and all data in the table is taken as is. Only relevant for classical input tables.");
	parser.add_argument('--collapse_by_metadata', default="", help="Metadata by which to collapse rows or columns, in case of a BIOM input table. If empty (default), then row (observation) or column (sample) names are use instead, depending on the flag --collapse_columns_instead_of_rows.");
	parser.add_argument('--no_group_names_in_file', action='store_true', dest="no_group_names_in_file", default=False, help="If set, then a group's name is merely its order in the 'groups' file (starting at 0). Otherwise, the first member of a group is taken as the group name and not as a member (default).");
	parser.add_argument('--single_line_groups', action='store_true', dest="single_line_groups", default=False, help='If not set, then each group in the groups file is defined across multiple lines (one member per line), with groups separated by at least one blank or purely whitespace line. Each group starts with a group name, unless --no_group_names_in_file is set. Prefixed or suffixed whitespace in each line is ignored, unless enclosed within quotes. If --single_line_groups, then each group is defined on a single line, with members delimited by --group_delimiter.');
	parser.add_argument('--disable_group_set_operations', action='store_true', dest="disable_group_set_operations", default=False, help="Don't allow group set operations in the groups file.");
	parser.add_argument('--dont_parse_group_metadata', action='store_true', dest="dont_parse_group_metadata", default=False, help="Don't parse metadata specifications in the groups file, and thus consider the entire line as a group name. By default, group names can be followed by a list of metadata entries, which will be included in the various output tables (if in BIOM format). Not relevant if --single_line_groups or --no_group_names_in_file is set.");
	parser.add_argument('--allow_duplicate_groups', action='store_true', dest="allow_duplicate_groups", default=False, help="Allow duplicate group definitions, i.e. merge all member lists of a group. By default, an error is thrown if a group is defined multiple times.");
	parser.add_argument('--group_members_defined_as', action='store', dest="group_members_defined_as", default='words', choices=['match','wildcards','regex','words'], help="How group members are defined. 'match' means members must be matched exactly. 'wildcards' means group members are interpreted as shell wildcard expressions (e.g. '*europaea*' matches everything containing 'europaea'). 'regex' means group members are regular expressions (python style). 'words' means members are defined as a list of complete words (separated by '*') that are to be matched (e.g. 'E*coli' will match 'E coli' and 'E;coli' but not 'Escherichia coli'). What constitutes a valid word can be modified via the option --valid_word_symbols. In all cases, flanking whitespace is ignored. (default: %(default)s)");
	parser.add_argument('--valid_word_symbols', action='store', default='-', help="Symbols allowed in words (e.g. in species names), in addition to alphanumeric characters. Only relevant if --group_members_defined_as is 'words'. (default: %(default)s)");

	parser.add_argument('--column_names_are_in', default='first_data_line', choices=['none','last_comment_line','first_data_line'], help="Where the table header (listing column names) is located. The option 'last_comment_line' means that column names are in the last comment line prior to any data. Only relevant for classical tables. If set to 'none', then column names are column indices (starting at 0).");
	parser.add_argument('--table_delimiter', default='\\t', help="Column delimiter in input table (if classical) and output tables. For tabs use '\\t' (default). Set to empty ('') for elastic whitespace.");
	parser.add_argument('--group_delimiter', default='\\t', help="Column delimiter in groups file, if the flag --single_line_groups is set. For tabs use '\\t' (default). Set to empty ('') for elastic whitespace. Irrelevant if --single_line_groups is not set.");
	parser.add_argument('-c','--comment_prefix', default='#', help='Comment prefix assumed in input files (default: %(default)s).');

	parser.add_argument('--include_leftovers_in_groups2records_table', action='store_true', dest="include_leftovers_in_groups2records_table", default=False, help='Include leftover records (i.e. not assigned to any group) in the group2records table (associated with the leftovers group).');
	parser.add_argument('--include_leftovers_in_group_overlaps_table', action='store_true', dest="include_leftovers_in_group_overlaps_table", default=False, help='Include leftover records (i.e. not assigned to any group) in the group_overlaps table (associated with the leftovers group).');

	# options for filtering input table
	parser.add_argument('--only_columns', default='', help="Comma-separated list (no spaces) of column indices (0-based) in the input table to include. Set to empty to not explicitly omit any column (default).");
	parser.add_argument('--omit_columns', default='', help='Comma-separated list (no spaces) of column indices (0-based) in the input table to ignore. Set to empty to not explicitly omit any column (default).');
	parser.add_argument('--only_rows', default='', help="Comma-separated list (no spaces) of row indices (0-based) in the input table to include. Set to empty to not explicitly omit any row (default). If --column_names_are_in is 'first_data_line' and the input table is in classical format, then the header line is not counted.");
	parser.add_argument('--omit_rows', default='', help="Comma-separated list (no spaces) of row indices (0-based) in the input table to ignore. Set to empty to not explicitly omit any row (default). If --column_names_are_in is 'first_data_line' and the input table is in classical format, then the header line is not counted.");
	parser.add_argument('--only_records', default="", help='Comma-separated list (no spaces) of record names to include when collapsing. Set to "" (empty) to include all records (default). Usage of wildcards (*) is allowed.');
	parser.add_argument('--omit_records', default="", help='Comma-separated list (no spaces) of record names to omit when collapsing. Set to "" (empty) to not explicitly omit any records (default), but also see option --only_samples. Usage of wildcards (*) is allowed.');
	parser.add_argument('--only_samples', default="", help='Comma-separated list (no spaces) of sample names to include when collapsing a BIOM table. Set to "" (empty) to include all samples (default). Usage of wildcards (*) is allowed. Samples in BIOM tables correspond to columns.');
	parser.add_argument('--omit_samples', default="", help='Comma-separated list (no spaces) of sample names to omit when collapsing a BIOM table. Set to "" (empty) to not explicitly omit any samples (default), but also see option --only_samples. Usage of wildcards (*) is allowed. Samples in BIOM tables correspond to columns.');
	parser.add_argument('--only_observations', default="", help='Comma-separated list (no spaces) of observation names to include when collapsing a BIOM table. Set to "" (empty) to include all observations (default). Usage of wildcards (*) is allowed. Observations in BIOM tables correspond to rows.');
	parser.add_argument('--omit_observations', default="", help='Comma-separated list (no spaces) of observation names to omit when collapsing a BIOM table. Set to "" (empty) to not explicitly omit any observations (default), but also see option --only_observations. Usage of wildcards (*) is allowed. Observations in BIOM tables correspond to rows.');
	parser.add_argument('--only_samples_by_metadata', default='', help="Metadata predicates for samples to include when collapsing a BIOM table (in addition to any filtering specified by --only_samples and --omit_samples). Specified as semicolon-separated list of predicates. Each predicate specifies a metadata category, followed by a colon and a comma-separated list of permissible values for that metadata. For example, 'continent:africa,europe;type:soil,mud' will only allow soil and mud samples from either europe or africa (i.e. ';' means AND). Leave empty to not explicitly omit any samples (default).");
	parser.add_argument('--omit_samples_by_metadata', default='', help="Metadata predicates for samples to omit when collapsing a BIOM table (in addition to any filtering specified by --only_samples and --omit_samples). Specified as semicolon-separated list of predicates. Each predicate specifies a metadata category, followed by a colon and a comma-separated list of non-permissible values for that metadata. For example, 'continent:asia;type:lake,pond' will exclude any lake or pond samples, as well as any samples from asia (i.e. ';' means OR). Leave empty to not explicitly omit any samples (default).");
	parser.add_argument('--only_observations_by_metadata', default='', help="Metadata predicates for observations to include when collapsing a BIOM table (in addition to any filtering specified by --only_observations and --omit_observations). Specified as semicolon-separated list of predicates. Each predicate specifies a metadata category, followed by a colon and a comma-separated list of permissible values (or wildcards) for that metadata. For example, 'genus:Nitrosomonas;shape:rod,disc' will only allow rod- or disc-shaped observations whose genus is Nitrosomonas. Leave empty to not explicitly omit any observations (default).");
	parser.add_argument('--omit_observations_by_metadata', default='', help="Metadata predicates for observations to omit when collapsing a BIOM table (in addition to any filtering specified by --only_observations and --omit_observations). Specified as semicolon-separated list of predicates. Each predicate specifies a metadata category, followed by a colon and a comma-separated list of non-permissible values (or wildcards) for that metadata. For example, 'genus:Nitrosomonas;shape:rod,disc' will exclude any rod- or disc-shaped observations, as well as any observations whose genus is Nitrosomonas. Leave empty to not explicitly omit any observations (default).");

	parser.add_argument('--record_order', default="given:f", help="How to sort records in the input table prior to collapsing, specified as a semicolon-separated list of ordering layers. Each layer specifies one particular algorithm for sorting, applied on top of previous sorting algorithms to resolve remaining ambiguities. Each layer is specified as a colon-separated list of the format <layer_type>:<direction>, where <layer_type> can be one of 'given','label','average_value','average_absolute_value','non_zeros','data:<data_name>','metadata':<metadata_name>, and <direction> can be one of 'f' for forward,'r' for reverse, 'n' for natural, 'nr' for natural-reverse or 'c':<custom_direction_list> for a custom order. <custom_direction_list> is a comma-separated list of values specifying a custom order to be applied on the particular layer. For example, the --observation_order specification 'label:f;metadata:depth:r;metadata:country:c:USA,Canada,Germany' specifies 3 successive sorting layers. If your table contains the columns 'A' and 'B', then 'data:A:f;data:B;f' will sort rows first by values in A, then by values in B. (default: %(default)s)");

	parser.add_argument('--non_numeric', default='consolidate', choices=['ignore', 'consolidate', 'first'], help="How to treat data with non-numeric contributions when collapsing a group. 'consolidate' (default) means the value assigned to a group will be NA unless all contributions are the same (e.g. if all member of a group have color 'red', the group will be 'red', but if some members are 'green', then the group will be 'NA'). 'ignore' means non-numeric data are not counted. 'first' means the first value listed is kept. If the output table is in BIOM format, consolidated data will be included as group metadata.");
	parser.add_argument('--average', default='none', choices=['none','across_records','across_group_members','across_used_group_members','maximum','minimum','minimum_across_records'], help="How to consolidate (e.g. average) group scores in collapsed table. 'none' (default) means that a group's score is the sum of scores across all records matched by any of the group's members. 'across_records' means a group's score is divided by the number of records assigned to the group. 'across_group_members' means a group's score is divided by the number of members of the group. 'across_used_group_members' means a group's score is divided by the number of members of the group that matched at least one record. 'maximum' and 'minimum' use the max & min value across members of a group (unmatched members are assumed to have value 0). 'maximum_across_records' and 'minimum_across_records' uses the max & min across all records assigned to a group.");
	parser.add_argument('--transpose_collapsed', action='store_true', dest="transpose_collapsed", default=False, help="Transpose collapsed output table, compared to the input table. By default, if the input table is classical and the output table is BIOM, then columns (or rows) in the input table correspond to samples (or observations) in the output table. This relation is reversed if --transpose_output is set.");
	parser.add_argument('--transpose_sub_tables', action='store_true', dest="transpose_sub_tables", default=False, help="Transpose output table, compared to the input table. By default, if the input table is classical and the output table is BIOM, then columns (or rows) in the input table correspond to samples (or observations) in the output table. This relation is reversed if --transpose_output is set.");
	parser.add_argument('--collapse_columns_instead_of_rows', action='store_true', dest="collapse_columns_instead_of_rows", default=False, help='Transpose input table, i.e. collapse columns instead of rows. ');
	parser.add_argument('--missing_entry', default='NA', help='Value to be used for undefined categorial values in output table (e.g. categorial data for unrepresented groups). (default: %(default)s)');
	parser.add_argument('--group_title', default='group', help="Column title to be used for groups (e.g. 'functional_group' or 'COG'). Irrelevant if --collapse_columns_instead_of_rows XOR --transpose_output is set. (default: %(default)s)");
	parser.add_argument('--data_title', default='', help="Column title to be used for listing data names (e.g. 'sample'). Irrelevant if not(--collapse_columns_instead_of_rows XOR --transpose_output). (default: %(default)s)");
	parser.add_argument('--keep_header_comments', action='store_true', dest="keep_header_comments", default=False, help="Keep lines in classical table header that start with the comment prefix and precede any data. The comment prefix must be the first sequence in the line. Only relevant for classical (tabular) output tables (see option --output_format_collapsed).");
	parser.add_argument('--include_summary_comments', action='store_true', dest="include_summary_comments", default=False, help="Include summary comments in the output table. This is independent of the optional log file (see option --log) or report file (see option --report). Only relevant for classical (tabular) output tables (see option --output_format_collapsed).");
	parser.add_argument('--case_sensitive', action='store_true', dest="case_sensitive", default=False, help="Column or row name matching should be case sensitive, instead of case insensitive.");
	parser.add_argument('--group_leftovers_as', default='', help="Records not included in any group shall be summarized as a separate group with this name. Leave empty (default) to ignore rows not included in any group.");
	parser.add_argument('--include_numbers_in_column_header', action='store_true', dest="include_numbers_in_column_header", default=False, help="Include column numbers in output table header (e.g. '1:depth 2:year')");
	parser.add_argument('--omit_unrepresented_groups', action='store_true', dest="omit_unrepresented_groups", default=False, help="Omit groups not represented in the input table.");
	parser.add_argument('--avoid_creating_empty_BIOM_tables', action='store_true', dest="avoid_creating_empty_BIOM_tables", default=False, help="Dont create certain BIOM output tables if empty (even if technically correct). This may be useful if some downstream programs can't handle empty BIOM tables as input.");
	parser.add_argument('--compress_sub_tables', action='store_true', dest="compress_sub_tables", default=False, help="Compress output sub-tables using gzip. Only relevant if --out_sub_tables_dir was specified.");

	# options for report format
	parser.add_argument('--identify_cultured_taxa', action='store_true', dest="identify_cultured_taxa", default=False, help="In the report, identify records corresponding to cultured taxa (e.g. taxon names not including 'uncultured' or 'metagenome' or 'unknown') in the report. This only makes sense if the table is collapsed based on taxonomic annotations.");
	parser.add_argument('--partition_each_group_by_scores', default='', help="Comma-separated list (no spaces) of score-thresholds (numbers between 0 and 1) by which to partition each group in the report file. The score of an entry is its total number of hits across all data columns, divided by the total number of hits across all entries and data columns. E.g. if this is '0.01,0.1,0.5' then each group is split into 3 partitions, containing entries with a score 0.01-0.1, a score 0.1-0.5 and a score >0.5, respectively. Make sure to include '0' in the list if you really want all entries listed. Within each partition, observations are listed alphabetically.");
	parser.add_argument('--report_list_full_records', action='store_true', dest="report_list_full_records", default=False, help="List the full record in the report, instead of just the record's row name.");
	parser.add_argument('--report_only_groups', default='', help="Optional comma-separated list of group names to include in the report. If empty (default), all groups are listed in the report.");
	parser.add_argument('--report_omit_groups', default='', help="Optional comma-separated list of group names to omit from the report. If empty (default), all groups are listed in the report.");

	parser.add_argument('-n','--normalize_collapsed', default="none", choices=["none", "columns_before_collapsing", "rows_before_collapsing","columns_after_collapsing","rows_after_collapsing","columns_before_collapsing_excluding_unassigned","rows_before_collapsing_excluding_unassigned"], help="How to normalize the table (before_collapsing normalizes the input table prior to processing, after_collapsing normalizes the output table). If the output table is in BIOM format, then columns correspond to 'samples' and rows to 'observations'. (default: '%(default)s)'");
	parser.add_argument('--normalize_sub_tables', default="none", choices=["none", "columns", "rows"], help="How to normalize sub-tables corresponding to individual groups (see option --sub_tables_dir). If the output tables are in BIOM format, then columns correspond to 'samples' and rows to 'observations'. (default: '%(default)s)'");
	parser.add_argument('--normalize_groups2records_table', default="none", choices=["none", "columns", "rows"], help="How to normalize groups2records association table (see option --groups2records_table). For BIOM output, columns correspond to 'samples' and rows to 'observations'. (default: '%(default)s)'");

	parser.add_argument('--output_format_collapsed', default="auto", choices=["auto", "BIOM", "HBIOM", "classical"], help="File format for the collapsed table. If 'auto' (default), then the output format will be BIOM iff the extension is .biom. For BIOM tables, the BIOM JSON type is always used.");
	parser.add_argument('--output_format_sub_tables', default="auto", choices=["auto", "BIOM", "HBIOM", "classical"], help="File format for the sub-tables corresponding to individual groups (see option --sub_tables_dir). If 'auto' (default), then the sub-table format will be the same as the input format (classical or BIOM). For BIOM tables, the BIOM JSON type is always used.");
	parser.add_argument('--output_format_groups2records_table', default="auto", choices=["auto", "BIOM", "HBIOM", "classical"], help="File format for the groups2records_table (see option --out_groups2records_table). If 'auto' (default), then the groups2records_table format will be BIOM iff the extension is .biom. For BIOM tables, the BIOM JSON type is always used.");
	parser.add_argument('--output_format_group_overlaps', default="auto", choices=["auto", "BIOM", "HBIOM", "classical"], help="File format for the group_overlaps (see option --out_group_overlaps). If 'auto' (default), then the group_overlaps format will be BIOM iff the extension is .biom. For BIOM tables, the BIOM JSON type is always used.");

	parser.add_argument('--verbose_prefix', default="  ", help="Line prefix to be used for standard output messages. This may be useful if the script is part of another pipeline. (default: '%(default)s)'");
	parser.add_argument('-f','--force', action='store_true', dest="force", default=False, help='Replace existing output files without warning.');
	parser.add_argument('-v','--verbose', action='store_true', dest="verbose", default=False, help='Show lots of information.');
	args = parser.parse_args()
	if(args.table_delimiter=='\\t'): args.table_delimiter='\t';
	if(args.group_delimiter=='\\t'): args.group_delimiter='\t';
				
	# basic file checking
	if(not os.path.isfile(args.input_table)):
		print(args.verbose_prefix+"ERROR: Input table '%s' does not exist or is not a file" % (args.input_table))
		sys.exit(1)
	if((args.input_groups_file!='') and (not os.path.isfile(args.input_groups_file))):
		print(args.verbose_prefix+"ERROR: Groups file '%s' does not exist or is not a file" % (args.input_groups_file))
		sys.exit(1)
	if(args.input_groups_file=='' and args.groups_list_front=='' and args.groups_list_back==''):
		print(args.verbose_prefix+"ERROR: Missing non-empty --groups_file or --groups_list_front or --groups_list_back")
		sys.exit(1)
	if((args.average in ['maximum','minimum','minimum_across_records','maximum_across_records']) and args.out_collapsed_deconvoluted_table_groups_vs_records!=""):
		print(args.verbose_prefix+"ERROR: A deconvoluted groups_vs_records table makes no sense when choosing --average='%s'"%(args.average))
		sys.exit(1)
	if(args.out_sub_tables_dir!=''): check_output_dir(args.out_sub_tables_dir,args.force,args.verbose,args.verbose_prefix)
	if(args.out_collapsed!=""): check_output_file(args.out_collapsed,args.force,args.verbose,args.verbose_prefix)
	if(args.out_report!=""): check_output_file(args.out_report,args.force,args.verbose,args.verbose_prefix)
	if(args.out_log!=""): check_output_file_append(args.out_log,args.force,args.verbose,args.verbose_prefix)
	if(args.out_groups2records_table!=""): check_output_file(args.out_groups2records_table,args.force,args.verbose,args.verbose_prefix)
	if(args.out_groups2records_table_dense!=""): check_output_file(args.out_groups2records_table_dense,args.force,args.verbose,args.verbose_prefix)
	if(args.out_group_definitions_used!=""): check_output_file(args.out_group_definitions_used,args.force,args.verbose,args.verbose_prefix)
	if(args.out_group_definitions_unused!=""): check_output_file(args.out_group_definitions_unused,args.force,args.verbose,args.verbose_prefix)
	if(args.out_group_overlaps!=""): check_output_file(args.out_group_overlaps,args.force,args.verbose,args.verbose_prefix)
	if(args.out_collapsed_deconvoluted_table_records_vs_data!=""): check_output_file(args.out_collapsed_deconvoluted_table_records_vs_data,args.force,args.verbose,args.verbose_prefix)
	if(args.out_collapsed_deconvoluted_table_groups_vs_records!=""): check_output_file(args.out_collapsed_deconvoluted_table_groups_vs_records,args.force,args.verbose,args.verbose_prefix)
	if(args.out_all_members!=""): check_output_file(args.out_all_members,args.force,args.verbose,args.verbose_prefix)
	#if(args.collapse_columns_instead_of_rows): args.row_names_are_in_column=""; # disable row-name extraction if collapsing columns
	saveCollapsedAsBIOM 		= ((args.output_format_collapsed.endswith('BIOM')) or (args.output_format_collapsed=='auto' and is_biom_file(args.out_collapsed)));
	saveSubtablesAsBIOM 		= ((args.output_format_collapsed.endswith('BIOM')) or (args.output_format_collapsed=='auto' and is_biom_file(args.input_table)));
	saveGroups2recordsAsBIOM 	= ((args.output_format_groups2records_table.endswith('BIOM')) or (args.output_format_groups2records_table=='auto' and is_biom_file(args.out_groups2records_table)));
	saveGroupOverlapsAsBIOM 	= ((args.output_format_group_overlaps.endswith('BIOM')) or (args.output_format_group_overlaps=='auto' and is_biom_file(args.out_group_overlaps)));
	output_delimiter 			= ("\t" if args.table_delimiter=="" else args.table_delimiter)
	
	
		
	# figure out which rows & columns to include (in case of a tabular table)
	only_columns = [int(c) for c in args.only_columns.split(',') if is_number(c)]
	if(len(only_columns)==0): only_columns = None;
	else: only_columns.sort()
	
	omit_columns = [int(c) for c in args.omit_columns.split(',') if is_number(c)]
	if(len(omit_columns)==0): omit_columns = None;
	
	only_rows = [int(r) for r in args.only_rows.split(',') if is_number(r)]
	if(len(only_rows)==0): only_rows = None;
	else: only_rows.sort()	

	omit_rows = [int(r) for r in args.omit_rows.split(',') if is_number(r)]
	if(len(omit_rows)==0): omit_rows = None;
	
	
	# read records table
	# row_metadata and column_metadata will only be relevant in case of BIOM files
	if(args.verbose): print(args.verbose_prefix+"Reading input table..")
	table,\
	row_names,\
	column_names,\
	found_row_names,\
	found_column_names,\
	header_lines,\
	row_metadata,\
	column_metadata,\
	input_table_summary = read_table(	args.input_table, 
										args.table_delimiter, 
										args.comment_prefix, 
										args.row_names_are_in_column, 
										args.column_names_are_in,
										only_columns,
										omit_columns,
										only_rows,
										omit_rows,
										args.only_samples,
										args.omit_samples,
										args.only_observations,
										args.omit_observations,
										args.only_samples_by_metadata,
										args.omit_samples_by_metadata,
										args.only_observations_by_metadata,
										args.omit_observations_by_metadata,
										(not args.case_sensitive),
										args.verbose,
										args.verbose_prefix)
	if(len(table)==0): 
		print(args.verbose_prefix+"WARNING: Input table is empty. Exiting")
		sys.exit(0)
	NR = len(table)
	NC = (0 if len(table)==0 else len(table[0]));
	
		
	# normalize input table prior to collapsing if requested
	if(args.normalize_collapsed=="columns_before_collapsing"):
		sums = [sum([(float(table[r][c]) if is_non_nan_number(table[r][c]) else 0) for r in range(NR)]) for c in range(NC)];
		sums = [(1 if s==0 else s) for s in sums];
		table = [[("%.10g"%(float(table[r][c])/sums[c]) if is_non_nan_number(table[r][c]) else table[r][c]) for c in range(NC)] for r in range(NR)];
	elif(args.normalize_collapsed=="rows_before_collapsing"):
		sums = [sum([(float(table[r][c]) if is_non_nan_number(table[r][c]) else 0) for c in range(NC)]) for r in range(NR)];
		sums = [(1 if s==0 else s) for s in sums];
		table = [[("%.10g"%(float(table[r][c])/sums[r]) if is_non_nan_number(table[r][c]) else table[r][c]) for c in range(NC)] for r in range(NR)];
	
			
	# transpose input table if needed, so that in the remaining analysis rows=records and columns=data_fields(e.g. samples)
	if(args.collapse_columns_instead_of_rows):
		table 				= [[table[r][c] for r in range(NR)] for c in range(NC)];
		full_records		= [output_delimiter.join(table[c]) for c in range(NC)]
		data_names 			= row_names;
		found_data_names	= found_row_names;
		record_metadata		= column_metadata;	# metadata of the collapsed axis. may be None
		data_metadata		= row_metadata;		# metadata of the non-collapsed axis. may be None
		record_names 		= column_names;
		
	else:	
		if(not found_column_names): data_names=["source_column_"+cn for cn in column_names]
		else: data_names 	= column_names
		full_records		= [output_delimiter.join(table[r]) for r in range(NR)]
		found_data_names	= found_column_names;
		record_metadata		= row_metadata;	# metadata of the collapsed axis. may be None
		data_metadata 		= column_metadata;	# metadata of the non-collapsed axis. may be None
		record_names 		= row_names;

	# sort input table if needed
	if(args.record_order!="given:f"):
		if(args.verbose): print("%sSorting records.."%(args.verbose_prefix))
		table_numerical = numpy_float_or_nan(numpy.array(table,dtype=object))
		record_order = sort_multilayered(	args.record_order,
											IDs = record_names,
											labels = record_names,
											table = table_numerical,
											data_names = data_names,
											all_metadata = record_metadata,
											table_axis = 0,
											verbose_prefix = args.verbose_prefix)
		if(record_order is not None):
			record_names 	= [record_names[i] for i in record_order]
			record_metadata = (None if (record_metadata is None) else [record_metadata[i] for i in record_order])
			table 			= [table[r] for r in record_order]
						
						
	# prepare record labels by which to collapse
	if((args.collapse_by_metadata!="") and (record_metadata is not None)):
		# use specific metadata for collapsing, instead of record names
		NR = len(table);
		record_labels = [(record_metadata[r][args.collapse_by_metadata] if ((record_metadata[r] is not None) and (args.collapse_by_metadata in record_metadata[r])) else None) for r in range(NR)]
		if(next((n for n,name in enumerate(record_labels) if (name is not None)),-1)<0):
			if(args.verbose): print("%sWARNING: All record names appear to be missing.\n%s         Maybe metadata name '%s' is wrong?"%(args.verbose_prefix,args.verbose_prefix,args.collapse_by_metadata))
		record_labels = arbitrary_metadata_values_to_record_name(record_labels, args.group_members_defined_as)
		collapsed_records_by_name = False
	else:
		record_labels = record_names;
		collapsed_records_by_name = True
		
	
	# filter records if needed
	records_to_keep = filter_by_name_and_metadata(	record_names,
													None,
													None,
													args.only_records,
													args.omit_records,
													"",
													"",
													args.case_sensitive)
	if(len(records_to_keep)<len(record_names)):
		if(args.verbose): print("%sFiltering out %d records.."%(args.verbose_prefix,len(record_names)-len(records_to_keep)))
		record_names 	= (None if (record_names is None) else [record_names[r] for r in records_to_keep])
		record_metadata = (None if (record_metadata is None) else [record_metadata[r] for r in records_to_keep])
		record_labels 	= [record_labels[r] for r in records_to_keep]
		full_records	= [full_records[r] for r in records_to_keep]
		table 			= [table[r] for r in records_to_keep]
	input_table_summary += "\nAfter (potentially) filtering out records based on name, obtained a table comprising %d records & %d data entries per record"%(len(records_to_keep),len(data_names))
		
	
	NR = len(record_names)
	ND = len(data_names)
	output_columns_are_groups = XOR(args.transpose_collapsed, args.collapse_columns_instead_of_rows);	
	
	# output table with have size NG*ND (or ND*NE if args.transpose_collapsed XOR args.collapse_columns_instead_of_rows)	
	

	
	# read groups
	if(args.verbose): print(args.verbose_prefix+"Reading groups..")
	group_members,\
	group_names,\
	group_metadata,\
	all_members,\
	all_unique_members = read_groups(	args.input_groups_file, 
										args.groups_list_front, 
										args.groups_list_back,
										args.group_delimiter, 
										args.comment_prefix, 
										args.no_group_names_in_file,
										args.single_line_groups,
										(not args.disable_group_set_operations),
										args.allow_duplicate_groups,
										args.verbose,
										args.verbose_prefix)
	N_groups_loaded = len(group_names)
	NG = N_groups_loaded
	N_all_members = len(all_members)
	N_all_unique_members = len(all_unique_members)
	if(next((len(g) for g in group_members if len(g)>0),-1)<0): print(args.verbose_prefix+"WARNING: All groups have empty member lists")
	duplicate_group_names = find_duplicates_in_list(group_names)
	if(len(duplicate_group_names)>0):
		print("%sERROR: Duplicate group names. The following %d group names were used multiple times:\n%s  %s\n"%(args.verbose_prefix,len(duplicate_group_names),args.verbose_prefix,("\n"+args.verbose_prefix+"  ").join(duplicate_group_names)))
		sys.exit(1)
		
		
	# Calculate effective number of members per group, taking into account set operations
	# This is only approximate, as it is hard to define the number of members for a group that includes set operations
	effective_number_of_members_per_group = calculate_effective_number_of_members_per_group(group_members,N_all_members)
	
	
	# figure out which records to count to which group
	if(args.verbose): print("%sAssigning %s to groups.."%(args.verbose_prefix,"columns" if args.collapse_columns_instead_of_rows else "rows"))
	group_to_records,\
	leftover_records,\
	group_members_used,\
	group_members_unused = assign_records_to_groups(group_names, 
													group_members,
													all_members,
													record_labels, 
													(not args.case_sensitive), 
													args.valid_word_symbols,
													args.group_members_defined_as);

	# Calculate effective number of used members per group, taking into account set operations
	effective_number_of_used_members_per_group = calculate_effective_number_of_members_per_group(group_members_used,N_all_members)
	
	# check if some groups are not represented in input table, and filter out if needed
	N_leftovers 			= len(leftover_records);
	N_annotations 			= sum(len(group_to_records[g]) for g in range(NG))
	is_represented 			= [len(group_to_records[g])>0 for g in range(NG)]
	N_groups_represented 	= sum(is_represented[g] for g in range(NG))
	if((NG>0) and args.omit_unrepresented_groups and (N_groups_represented<NG)):
		# some groups are not represented in input table, so filter out
		# make sure to also correct group references in group set operations
		if(args.verbose): 
			print(args.verbose_prefix+"Note: No entries found for the following groups:")
			for g in range(NG): 
				if(not is_represented[g]): print("        '"+group_names[g]+"'")
		# correct group references before deleting
		deleted_groups_until = [0,]*NG;
		deleted_groups_until[0] = (0 if is_represented[0] else 1);
		for g in range(1,NG):
			deleted_groups_until[g] = deleted_groups_until[g-1]
			if(not is_represented[g]): deleted_groups_until[g] += 1;
		for g in range(NG):
			group_members[g] = [m if isinstance(m,int) else ([m[0],-1] if (not is_represented[m[1]]) else [m[0],m[1]-deleted_groups_until[m[1]]]) for m in group_members[g]]
			group_members_used[g] = [m if isinstance(m,int) else ([m[0],-1] if (not is_represented[m[1]]) else [m[0],m[1]-deleted_groups_until[m[1]]]) for m in group_members_used[g]]
			group_members_unused[g] = [m if isinstance(m,int) else ([m[0],-1] if (not is_represented[m[1]]) else [m[0],m[1]-deleted_groups_until[m[1]]]) for m in group_members_unused[g]]
		# delete groups
		effective_number_of_members_per_group = [effective_number_of_members_per_group[g] for g in range(NG) if is_represented[g]];
		effective_number_of_used_members_per_group = [effective_number_of_used_members_per_group[g] for g in range(NG) if is_represented[g]];
		group_to_records 	= [group_to_records[g] for g in range(NG) if is_represented[g]];
		group_members 		= [group_members[g] for g in range(NG) if is_represented[g]];
		group_members_used	= [group_members_used[g] for g in range(NG) if is_represented[g]];
		group_members_unused= [group_members_unused[g] for g in range(NG) if is_represented[g]];
		group_names 		= [group_names[g] for g in range(NG) if is_represented[g]];
		group_metadata		= [group_metadata[g] for g in range(NG) if is_represented[g]];
		NG 					= len(group_names)
			
	# include leftover records as an additional group if needed
	if(not(args.group_leftovers_as=="") and (len(leftover_records)>0)):
		group_names.append(args.group_leftovers_as);
		effective_number_of_members_per_group.append(1.0);
		effective_number_of_used_members_per_group.append(1.0);
		group_metadata.append({});
		group_to_records.append(leftover_records);
		group_members.append([])
		group_members_used.append([]);
		group_members_unused.append([]);
		NG = len(group_names)
		NG_excluding_leftovers = NG-1	# for routines that only want to access actual groups, keep track of which ones these were
	else:
		NG_excluding_leftovers = NG
		
	generic_summary_info = "# Generated on: %s\n# Used command:\n#   %s\n#\n# Summary:\n#   Loaded %d groups comprising %d members (%d unique members)\n#   Established %d assignments of records to groups\n#   %d out of %d records could not be assigned to any group (%s)\n#   Number of groups represented: %d\n#"%(get_date_time(),get_shell_command(),N_groups_loaded,N_all_members,N_all_unique_members,N_annotations,N_leftovers,NR,('leftovers' if (args.group_leftovers_as=="") else args.group_leftovers_as),N_groups_represented)

			
	# write group-record associations table if requested
	if(args.out_groups2records_table!=""):
		if(args.verbose): print(args.verbose_prefix+"Saving group-record associations as table..")
		NG_for_groups2records = (NG if args.include_leftovers_in_groups2records_table else NG_excluding_leftovers)
		group_record_associations = numpy.zeros((NR, NG_for_groups2records), dtype=float)
		for g in range(NG_for_groups2records):
			for r in group_to_records[g]:
				group_record_associations[r,g] = 1.0;
		group_record_associations = normalize_table(group_record_associations,args.normalize_groups2records_table);
		if(saveGroups2recordsAsBIOM):
			# BIOM format
			biom_table = biom.table.Table(	group_record_associations, 
											observation_ids=record_names, 
											sample_ids=group_names[0:NG_for_groups2records], 
											observation_metadata=record_metadata, 
											sample_metadata=group_metadata[0:NG_for_groups2records]);
			save_biom_table(biom_table, args.out_groups2records_table, args.output_format_groups2records_table)
		else:
			# classical format
			with open(args.out_groups2records_table,'w') as fout:
				fout.write("# Group-record associations within table '%s'\n# Based on groups table '%s'\n%s\n" % (args.input_table,args.input_groups_file,generic_summary_info))
				fout.write("record%s%s\n"%(output_delimiter,output_delimiter.join(group_names[0:NG_for_groups2records])));
				for r in range(NR):
					fout.write("%s%s%s\n"%(record_names[r],output_delimiter,output_delimiter.join('%.10g'%group_record_associations[r,g] for g in range(NG_for_groups2records))))
			del group_record_associations # clear
			

	# write dense lists of groups associated with each record, as classical table
	if(args.out_groups2records_table_dense!=""):
		if(args.verbose): print(args.verbose_prefix+"Saving group-record associations as dense table (lists of groups per record)..")
		record2groups = [[] for r in range(NR)]
		for g in range(NG_excluding_leftovers):
			for r in group_to_records[g]:
				record2groups[r] += [g]
		with open(args.out_groups2records_table_dense,'w') as fout:
			fout.write("# Groups associated with records, within table '%s'\n# Based on groups table '%s'\n%s\n" % (args.input_table,args.input_groups_file,generic_summary_info))
			fout.write("record%sgroup%s\n"%(output_delimiter,("" if collapsed_records_by_name else output_delimiter+"label")));
			for r in range(NR):
				fout.write("%s%s%s%s\n"%(record_names[r], output_delimiter, ",".join(group_names[g] for g in record2groups[r]), ("" if collapsed_records_by_name else output_delimiter+record_labels[r])))
		del record2groups # clear
				
				
	# write group overlaps (Jaccard similarity index) as table if requested
	if(args.out_group_overlaps!=""):
		if(args.verbose): print(args.verbose_prefix+"Saving group overlaps as table..")
		NG_for_group_overlaps = (NG if args.include_leftovers_in_group_overlaps_table else NG_excluding_leftovers)
		group_overlaps = numpy.zeros((NG_for_group_overlaps, NG_for_group_overlaps), dtype=float)
		for g1 in range(NG_for_group_overlaps):
			for g2 in range(g1,NG_for_group_overlaps):
				group_overlaps[g1,g2] = get_jaccard_index(group_to_records[g1], group_to_records[g2])
				group_overlaps[g2,g1] = group_overlaps[g1,g2];
		if(saveGroupOverlapsAsBIOM):
			# BIOM format
			biom_table = biom.table.Table(	group_overlaps, 
											observation_ids=group_names[0:NG_for_group_overlaps], 
											sample_ids=group_names[0:NG_for_group_overlaps], 
											observation_metadata=group_metadata[0:NG_for_group_overlaps], 
											sample_metadata=group_metadata[0:NG_for_group_overlaps]);
			save_biom_table(biom_table, args.out_group_overlaps, args.output_format_group_overlaps)
		else:
			# classical table
			with open(args.out_group_overlaps,'w') as fout:
				fout.write("# Group overlaps (Jarrard similarity index) in terms of records shared from the table '%s'\n# Based on groups table '%s'\n%s\n" % (args.input_table,args.input_groups_file,generic_summary_info))
				fout.write("%s%s%s\n"%(args.group_title,output_delimiter,output_delimiter.join(group_names[0:NG_for_group_overlaps])));
				for g1 in range(NG_for_group_overlaps):
					fout.write("%s%s%s\n"%(group_names[g1],output_delimiter,output_delimiter.join('%.10g'%group_overlaps[g1,g2] for g2 in range(NG_for_group_overlaps))))
		group_overlaps = None; # clear


	# split input table into a numerical and a categorial part
	interpret_data_as_numbers = [True]*ND
	for d in range(ND):
		for g in range(NG):
			all_numeric_for_this_group = (next((r for r in group_to_records[g] if not is_number_or_nan(table[r][d])), -1)<0);
			if(not (all_numeric_for_this_group or args.non_numeric=='ignore')):
				interpret_data_as_numbers[d] = False;
	table_numerical_part 		= numpy.array([[float_or_nan(table[r][d]) for d in range(ND) if interpret_data_as_numbers[d]] for r in range(NR)])
	numerical_data_names 		= [data_names[d] for d in range(ND) if interpret_data_as_numbers[d]]		
	numerical_data_metadata 	= (None if (data_metadata is None) else [data_metadata[d] for d in range(ND) if interpret_data_as_numbers[d]])
	ND_numerical				= len(numerical_data_names)
	table_categorial_part		= [[table[r][d] for d in range(ND) if (not interpret_data_as_numbers[d])] for r in range(NR)]
	categorial_data_names 		= [data_names[d] for d in range(ND) if (not interpret_data_as_numbers[d])]
	categorial_data_metadata 	= (None if (data_metadata is None) else [data_metadata[d] for d in range(ND) if (not interpret_data_as_numbers[d])])
	ND_categorial				= len(categorial_data_names)
	
	# normalize input table (numerical part) after removing unassigned, if requested
	if(args.normalize_collapsed in ["columns_before_collapsing_excluding_unassigned", "rows_before_collapsing_excluding_unassigned"]):
		assigned_records = [r for r in range(NR) if (not r in leftover_records)]
		if(args.normalize_collapsed=="columns_before_collapsing_excluding_unassigned"):
			sums = numpy.nansum(table_numerical_part[assigned_records,:], axis=0);
			sums[sums==0] = 1.0
			table_numerical_part /= sums
		elif(args.normalize_collapsed=="rows_before_collapsing_excluding_unassigned"):
			sums = numpy.nansum(table_numerical_part[assigned_records,:], axis=1);
			sums[sums==0] = 1.0
			table_numerical_part /= sums[:,numpy.newaxis];
		
	# define right-deconvoluted part of the collapsed table
	collapsed_deconvoluted_table_records_vs_data = numpy.copy(table_numerical_part)
	
	# define left-deconvoluted part of the collapsed table
	if(args.out_collapsed_deconvoluted_table_groups_vs_records!=""):
		collapsed_deconvoluted_table_groups_vs_records = numpy.zeros((NG, NR), dtype=float)
		for g in range(NG):
			collapsed_deconvoluted_table_groups_vs_records[g,list(group_to_records[g])] = 1.0;
	else:
		collapsed_deconvoluted_table_groups_vs_records = None
				
	# collapse table
	if(args.verbose): print(args.verbose_prefix+"Collapsing table..")
	collapsed_table_numerical_part 	= numpy.zeros([NG, ND_numerical])
	collapsed_table_categorial_part = [None,]*NG
	for g in range(NG):
		records_for_this_group = list(group_to_records[g])
		if(len(group_to_records[g])==0):
			# this group is not represented in the table
			if(args.verbose): print(args.verbose_prefix+"Note: No entry found for group '%s'" % (group_names[g]))
			collapsed_table_categorial_part[g] 	= [args.missing_entry]*ND_categorial
		else:
			# numerical part
			collapsed_table_numerical_part[g,:] = numpy.nansum(table_numerical_part[records_for_this_group,:],axis=0)
			if(args.average=='across_records'):
				collapsed_table_numerical_part[g,:] /= float(len(records_for_this_group))
				collapsed_deconvoluted_table_groups_vs_records[g,:] /= float(len(records_for_this_group))
			elif(args.average=='across_group_members'):
				collapsed_table_numerical_part[g,:] /= max(1.0,effective_number_of_members_per_group[g])
				collapsed_deconvoluted_table_groups_vs_records[g,:] /= max(1.0,effective_number_of_members_per_group[g])
			elif(args.average=='across_used_group_members'):
				collapsed_table_numerical_part[g,:] /= max(1.0,effective_number_of_used_members_per_group[g])
				collapsed_deconvoluted_table_groups_vs_records[g,:] /= max(1.0,effective_number_of_used_members_per_group[g])
			elif(args.average=='maximum_across_records'):
				collapsed_table_numerical_part[g,:] = numpy.nanmax(table_numerical_part[records_for_this_group,:],axis=0)
			elif(args.average=='minimum_across_records'):
				collapsed_table_numerical_part[g,:] = numpy.nanmin(table_numerical_part[records_for_this_group,:],axis=0)
			elif(args.average=='maximum'):
				collapsed_table_numerical_part[g,:] = numpy.nanmax(table_numerical_part[records_for_this_group,:],axis=0)
				if(len(group_members[g])>len(group_members_used[g])):				
					collapsed_table_numerical_part[g,:] = numpy.maximum(collapsed_table_numerical_part[g,:], 0)
			elif(args.average=='minimum'):
				collapsed_table_numerical_part[g,:] = numpy.nanmin(table_numerical_part[records_for_this_group,:],axis=0)
				if(len(group_members[g])>len(group_members_used[g])):				
					collapsed_table_numerical_part[g,:] = numpy.minimum(collapsed_table_numerical_part[g,:], 0)
				
			# categorial part
			collapsed_table_categorial_part[g] = [None,]*ND_categorial
			if(args.non_numeric=='consolidate'):
				for dc in range(ND_categorial):
					collapsed_table_categorial_part[g][dc] = consolidate_categorial([table_categorial_part[r][dc] for r in records_for_this_group], args.missing_entry);
			elif(args.non_numeric=='first'):
				for dc in range(ND_categorial):
					collapsed_table_categorial_part[g][dc] = table_categorial_part[records_for_this_group[0]][dc]
				
		
	# At this point collapsed_table_numerical_part[g,dn] (float) is value of numerical data dn for group with name group_names[g]. Can be NaN.
	# Similarly, collapsed_table_categorial_part[g][dc] (string) is value of categorial data dc for group with name group_names[g].
	if(args.verbose): print(args.verbose_prefix+"Assigned %d records to groups, %d records were leftovers"%(NR-N_leftovers,N_leftovers))
	
			
	# normalize collapsed numerical table, if requested
	# also apply changes to deconvoluted parts (if available)
	if(args.normalize_collapsed in ["columns_after_collapsing", "rows_after_collapsing"]):
		if(XOR(args.normalize_collapsed=="columns_after_collapsing",output_columns_are_groups)):
			sums = numpy.nansum(collapsed_table_numerical_part, axis=0);
			sums[sums==0] = 1.0
			collapsed_table_numerical_part /= sums
			if(collapsed_deconvoluted_table_records_vs_data is not None): collapsed_deconvoluted_table_records_vs_data /= sums;
		elif(XOR(args.normalize_collapsed=="rows_after_collapsing",output_columns_are_groups)):
			sums = numpy.nansum(collapsed_table_numerical_part, axis=1);
			sums[sums==0] = 1.0
			collapsed_table_numerical_part /= sums[:,numpy.newaxis];
			if(collapsed_deconvoluted_table_groups_vs_records is not None): collapsed_deconvoluted_table_groups_vs_records /= sums[:,numpy.newaxis]


	# write right deconvoluted collapsed table (numerical part)
	if((args.out_collapsed_deconvoluted_table_records_vs_data!="") and (collapsed_deconvoluted_table_records_vs_data is not None)):
		if(args.verbose): print("%sWriting right deconvoluted collapsed table (records vs numerical data).."%(args.verbose_prefix))
		biom_table = biom.table.Table(	collapsed_deconvoluted_table_records_vs_data, 
										observation_ids=record_names, 
										sample_ids=numerical_data_names, 
										observation_metadata=record_metadata, 
										sample_metadata=numerical_data_metadata);
		save_biom_table(biom_table, args.out_collapsed_deconvoluted_table_records_vs_data, "BIOM")
	
	# write left deconvoluted collapsed table (numerical part)
	if((args.out_collapsed_deconvoluted_table_groups_vs_records!="") and (collapsed_deconvoluted_table_groups_vs_records is not None)):
		if(args.verbose): print("%sWriting left deconvoluted collapsed table (groups vs records).."%(args.verbose_prefix))
		biom_table = biom.table.Table(	collapsed_deconvoluted_table_groups_vs_records, 
										observation_ids=group_names, 
										sample_ids=record_names, 
										observation_metadata=group_metadata, 
										sample_metadata=record_metadata);
		save_biom_table(biom_table, args.out_collapsed_deconvoluted_table_groups_vs_records, "BIOM")
		
	
	# write collapsed table
	if(args.out_collapsed!=""):
		if(args.verbose): print("%sWriting collapsed table '%s'.."%(args.verbose_prefix,args.out_collapsed))
		if(NG==0 and args.verbose): print("%sWARNING: Collapsed table is empty"%(args.verbose_prefix))
		if(saveCollapsedAsBIOM):
			# BIOM format
			if((NG>0) or (not args.avoid_creating_empty_BIOM_tables)):
				collapsed_table_categorial_part_as_group_metadata = [{categorial_data_names[dc]:collapsed_table_categorial_part[g][dc] for dc in range(ND_categorial)} for g in range(NG)]
				for g in range(NG): collapsed_table_categorial_part_as_group_metadata[g].update(group_metadata[g]); # add any intrinsic group metadata
				if(output_columns_are_groups):
					# groups are 'samples' in BIOM table
					biom_table = biom.table.Table(	numpy.transpose(collapsed_table_numerical_part), 
													observation_ids=numerical_data_names,
													sample_ids=group_names,
													observation_metadata=numerical_data_metadata, # metadata of the non-collapsed axis. may be None
													sample_metadata=collapsed_table_categorial_part_as_group_metadata);
				else:
					# groups are 'observations' in BIOM table
					biom_table = biom.table.Table(	collapsed_table_numerical_part, 
													observation_ids=group_names, 
													sample_ids=numerical_data_names, 
													observation_metadata=collapsed_table_categorial_part_as_group_metadata, 
													sample_metadata=numerical_data_metadata);		# metadata of the non-collapsed axis. may be None
				save_biom_table(biom_table, args.out_collapsed, args.output_format_collapsed)

		else:
			# Classical (tabular) format
			with (gzip.open(args.out_collapsed,'wt') if args.out_collapsed.lower().endswith(".gz") else open(args.out_collapsed,'wt')) as fout:
				if(args.include_summary_comments): fout.write("# Collapsed table '%s'\n%s\n" % (args.input_table,generic_summary_info))
				if((header_lines is not None) and args.keep_header_comments): fout.write("# Original (uncollapsed) table header:\n#%s\n"%("\n#".join(header_lines)))
				if(output_columns_are_groups):
					omit_data_names = ((data_names is None) or (not found_data_names))
					fout.write(('#' if args.column_names_are_in=='last_comment_line' else '') + ("" if omit_data_names else args.data_title+output_delimiter) + output_delimiter.join([(str(g+(0 if omit_data_names else 1))+":" if args.include_numbers_in_column_header else '') + group_names[g] for g in range(NG)])+"\n")
					dc = -1; dn = -1;
					for d in range(ND):
						if(interpret_data_as_numbers[d]):
							dn += 1;
							fout.write(("" if omit_data_names else numerical_data_names[dn]+output_delimiter)+output_delimiter.join("%.10g"%collapsed_table_numerical_part[g,dn] for g in range(NG))+"\n")
						else:
							dc += 1;
							fout.write(("" if omit_data_names else categorial_data_names[dc]+output_delimiter)+output_delimiter.join(collapsed_table_categorial_part[g][dc] for g in range(NG))+"\n")
	
				else:
					fout.write(('#' if args.column_names_are_in=='last_comment_line' else '') + ("0:" if (args.include_numbers_in_column_header and (args.group_title!="")) else '') + args.group_title + output_delimiter + output_delimiter.join([(str(d+1)+":" if args.include_numbers_in_column_header else '') + data_names[d] for d in range(ND)])+"\n")
					for g in range(NG):
						fout.write(group_names[g]);
						dc = -1; dn = -1;
						for d in range(ND):
							if(interpret_data_as_numbers[d]):
								dn += 1;
								fout.write(output_delimiter+"%.10g"%collapsed_table_numerical_part[g,dn])
							else:
								dc += 1;
								fout.write(output_delimiter+collapsed_table_categorial_part[g][dc])
						fout.write("\n");


	
	
		
	# save sub-tables of input table corresponding to individual groups
	if(args.out_sub_tables_dir!=""):
		if(args.verbose): print(args.verbose_prefix+"Saving sub-tables corresponding to groups..")
		transpose_sub_tables_wrt_table = XOR(args.transpose_sub_tables, args.collapse_columns_instead_of_rows);	
		for g in range(NG):
			if(len(group_to_records[g])!=0):
				if(saveSubtablesAsBIOM):
					sub_table = table_numerical_part[list(group_to_records[g]),:]
					if(args.avoid_creating_empty_BIOM_tables and saveSubtablesAsBIOM and numpy.all(sub_table==0)):
						if(args.verbose): print("%sNote: Omitting empty output sub-table '%s', as requested"%(args.verbose_prefix,group_names[g]))
						continue;
					sub_table_data_names	= numerical_data_names
					sub_table_data_metadata	= numerical_data_metadata
				else:
					sub_table = numpy.array([table[r] for r in group_to_records[g]],dtype="object")
					sub_table_data_names = data_names
					sub_table_data_metadata = data_metadata
				sub_table_record_names 		= [record_names[r] for r in group_to_records[g]]
				sub_table_record_metadata	= (None if (record_metadata is None) else [record_metadata[r] for r in group_to_records[g]])
				sub_table_file = os.path.join(args.out_sub_tables_dir,group_names[g]+(".biom" if saveSubtablesAsBIOM else ".txt")+(".gz" if args.compress_sub_tables else ""));
				check_output_file(sub_table_file,args.force,args.verbose,args.verbose_prefix)
				save_subtable(	sub_table_file,
								group_names[g],
								saveSubtablesAsBIOM,
								sub_table,
								(sub_table_data_names if transpose_sub_tables_wrt_table else sub_table_record_names),
								(sub_table_record_names if transpose_sub_tables_wrt_table else sub_table_data_names),
								(header_lines if args.keep_header_comments else None),
								(sub_table_data_metadata if transpose_sub_tables_wrt_table else sub_table_record_metadata),
								(sub_table_record_metadata if transpose_sub_tables_wrt_table else sub_table_data_metadata),
								args.normalize_sub_tables,
								args.include_summary_comments,
								args.column_names_are_in,
								args.include_numbers_in_column_header,
								args.table_delimiter,
								"#",
								args.verbose,
								args.verbose_prefix)
	
	
	
	# write log if needed
	if(args.out_log!=""):
		with (gzip.open(args.out_log,'wt') if args.out_log.lower().endswith(".gz") else open(args.out_log,'wt')) as fout:
			fout.write("%s:\n  Table '%s' collapsed according to groups defined in '%s'\n  Generated on: %s\n  Used command: %s\n" % (args.out_collapsed,args.input_table,args.input_groups_file,get_date_time(),get_shell_command()))
			if(args.out_sub_tables_dir!=""): fout.write("Sub-tables corresponding to groups saved to '%s'\n"%(args.out_sub_tables_dir));
			fout.write("\n");
	

	# save list of all members found in groups definitions
	if(args.out_all_members):
		if(args.verbose): print("%sSaving list of all %d group members found in groups definitions, to file.."%(args.verbose_prefix,N_all_unique_members))
		with (gzip.open(args.out_all_members,'wt') if args.out_all_members.lower().endswith(".gz") else open(args.out_all_members,'wt')) as fout:
			fout.write("# List of all members found in group definitions, loaded from file: '%s'\n%s\n"%(args.input_groups_file,generic_summary_info));
			for member in all_unique_members:
				fout.write("%s\n"%member)
	
	
	# create report file
	if(args.out_report!=""):
		# calculate record scores if needed
		if(args.partition_each_group_by_scores!=""): 
			if(args.verbose): print(args.verbose_prefix+"Calculating scores (fractions)..")
			scoreType = "fraction"
			record_scores = numpy.nansum(table_numerical_part,axis=1)
			total_score = sum(record_scores);
			if(total_score>0): record_scores /= float(total_score);
			# prepare score thresholds list
			scoreThresholds = sorted([float(s) for s in args.partition_each_group_by_scores.split(',') if is_non_nan_number(s)]);
			NP = len(scoreThresholds);
			partition_names = [("%g<=%s<%g"%(scoreThresholds[p],scoreType,scoreThresholds[p+1]) if (p+1<NP) else "%s>=%g"%(scoreType,scoreThresholds[p])) for p in range(NP)]
		else:
			NP = 0;
			
		# write groups to report
		if(args.verbose): print(args.verbose_prefix+"Writing report '%s'.." % args.out_report)
		with (gzip.open(args.out_report,'wt') if args.out_report.lower().endswith(".gz") else open(args.out_report,'wt')) as fout:
			fout.write("# Report for collapsed table: %s\n# Collapsed table generated on: %s\n# Original table: %s\n# Details:\n#  %s\n# Used command:\n#   %s\n" % (args.out_collapsed,get_date_time(),args.input_table,input_table_summary.replace("\n","\n#  "),get_shell_command()))
			if(args.out_sub_tables_dir!=""): fout.write("# Sub-tables corresponding to groups were saved to '%s'\n"%(args.out_sub_tables_dir));
			fout.write("#\n# Summary of group assignments:\n");
			for g in range(NG):
				fout.write("#   %s: %d records\n" % (group_names[g], len(group_to_records[g])))
			fout.write("#\n# Loaded %d groups comprising %d members (%d unique members)\n# Established %d assignments of records to groups\n# %d out of %d records (%g %%) were assigned to at least one group\n# %d out of %d records (%g %%) could not be assigned to any group (%s)\n# %d groups were represented (i.e. associated with at least one record)\n# Detailed group assignments are listed below\n"%(N_groups_loaded,N_all_members,N_all_unique_members,N_annotations,(NR-N_leftovers),NR,(NR-N_leftovers)/(0.01*NR),N_leftovers,NR,N_leftovers/(0.01*NR),('leftovers' if (args.group_leftovers_as=="") else "'"+args.group_leftovers_as+"'"),N_groups_represented));
			
			# filter groups if needed
			include_groups = filter_by_name_and_metadata(	group_names,
															None,
															None,
															args.report_only_groups,
															args.report_omit_groups,
															"",
															"",
															args.case_sensitive)
			if(len(include_groups)<NG):
				fout.write("# Note: Only a subset of %d groups is listed below, as requested\n"%(len(include_groups)))
			fout.write("\n")
			
			# list records assigned to each group
			for g in include_groups:		
				if(NP>0):
					# partition group and sort within:
					partitions = partitionIndexListByScores(group_to_records[g], record_labels, record_scores, scoreThresholds);
				
					# write group by partitions:
					fout.write("# %s (%d records):" % (group_names[g], len(group_to_records[g])));
					for p in range(NP-1,-1,-1): # reverse-traverse partition list (highest to lowest score)
						if(len(partitions[p])>0):
							fout.write("\n#  %s (%d records):\n" % (partition_names[p],len(partitions[p])))
							for r in partitions[p]:
								if(args.report_list_full_records):
									fout.write("      %s%s%s" % (record_names[r],output_delimiter,full_records[r]))
								else:
									fout.write("      %s%s"%(record_names[r],("" if collapsed_records_by_name else " (%s)"%record_labels[r])))
								fout.write("  # %s%s: %g\n"%((("(cultured) " if is_cultured_taxon(record_labels[r]) else "(uncultured) ") if args.identify_cultured_taxa else ""),scoreType,record_scores[r]))
					fout.write("\n\n");
				else:
					# don't partition group, but still sort:
					this_group_sorted_records = sorted(group_to_records[g], key=lambda r: record_labels[r].lower())
					fout.write("# %s (%d records):\n" % (group_names[g], len(group_to_records[g])))
					for r in this_group_sorted_records:
						if(args.report_list_full_records):
							fout.write("    %s%s%s" % (record_names[r],output_delimiter,full_records[r]))						
						else:
							fout.write("    %s" % (record_names[r]))
							if(not collapsed_records_by_name): fout.write(" (%s)"%(record_labels[r]));
						if(args.identify_cultured_taxa): fout.write("  # "+("(cultured)" if is_cultured_taxon(record_labels[r]) else "(uncultured)"))
						fout.write("\n");
					fout.write("\n\n")
	
	
	# make a list of group definitions (group member lists) actually used for the input data
	if(args.out_group_definitions_used!=""):
		if(args.verbose): print(args.verbose_prefix+"Saving subset of group definitions used in this analysis..")
		with (gzip.open(args.out_group_definitions_used,'wt') if args.out_group_definitions_used.lower().endswith(".gz") else open(args.out_group_definitions_used,'wt')) as fout:
			fout.write("# Group definitions relevant to input table '%s'\n%s\n"%(args.input_table,generic_summary_info));
			for g in range(NG):
				fout.write("%s\n# %s\n"%(group_names[g], '- '*(1+len(group_names[g])//2)))
				for m in group_members_used[g]:
					if(isinstance(m, int)):
						fout.write("%s\n"%(all_members[m][1]))
					else:
						# member is a group set operation, i.e. a list [operation_type, referenced_group]
						if(m[1]<0):
							# for some reason the referenced group does not exist anymore (may have been removed from the list)
							pass;
						else:
							fout.write("%s%s\n"%(set_operations_keywords[m[0]],group_names[m[1]]))					
				fout.write("\n\n")

	
	# make a list of group definitions (group member lists) not used for the input data, i.e. not matched to any record
	if(args.out_group_definitions_unused!=""):
		if(args.verbose): print(args.verbose_prefix+"Saving subset of group definitions not used in this analysis..")
		with (gzip.open(args.out_group_definitions_unused,'wt') if args.out_group_definitions_unused.lower().endswith(".gz") else open(args.out_group_definitions_unused,'wt')) as fout:
			fout.write("# Group definitions irrelevant to input table: '%s'\n%s\n"%(args.input_table,generic_summary_info));
			for g in range(NG):
				fout.write("%s\n# %s\n"%(group_names[g], '- '*(1+len(group_names[g])//2)))
				for m in group_members_unused[g]:
					if(isinstance(m, int)):
						fout.write("%s\n"%(all_members[m][1]))
					else:
						# member is a group set operation, i.e. a list [operation_type, referenced_group]
						if(m[1]<0):
							# for some reason the referenced group does not exist anymore (may have been removed from the list)
							pass;
						else:
							fout.write("%s%s\n"%(set_operations_keywords[m[0]],group_names[m[1]]))					
				fout.write("\n\n")
	
	sys.exit(0)
	
	
	
	
	
	
	
	
	