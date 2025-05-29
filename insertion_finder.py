#!/usr/bin/env python
# -*- coding: utf-8 -*-
import time
import traceback
import warnings
import sys
def warn_with_traceback(message, category, filename, lineno, file=None, line=None):
    log = file if hasattr(file,'write') else sys.stderr
    traceback.print_stack(file=log)
    log.write(warnings.formatwarning(message, category, filename, lineno, line))
    log.flush()
warnings.showwarning = warn_with_traceback
from Bio import SeqIO
from Bio.Blast.Applications import NcbiblastnCommandline
from datetime import datetime
import argparse
import math
import multiprocessing
import os
import random
import re
import shutil
import subprocess
from concurrent.futures import ProcessPoolExecutor
from concurrent.futures import as_completed
start_time = datetime.now()
call=os.path.abspath(os.getcwd())
param=dict()
version="2.6.0"
help = 'insertion_finder v{} - element insertion finder in a genome through a BLAST search\n'.format(version)
help = help + '(c) 2021. Arthur Gruber & Giuliana Pola\n'
help = help + 'For the latest version acess: https://github.com/GiulianaPola/insertion_finder\n'
help = help + 'Usage: insertion_finder.py -q <query file> -d <database file> -run local\n'
help = help + '\tinsertion_finder.py -q <query file> -d <database file> -run local -tab <BLASTn table file>\n'
help = help + '\tinsertion_finder.py -q <query file> -run web -tab <BLASTn table file>\n'
help = help + '\tinsertion_finder.py -q <query file> -run web\n'
help = help + '\nMandatory parameters:\n'
help = help + '-q <fasta or multifasta file>\tSequence to search with\n'
help = help + '-d <multifasta file>\tDatabase to BLAST against\n'
help = help + '-run <local|web>\tchoice of running local or web BLAST search\n'
help = help + '\nOptional parameters:\n'
help = help + '-tab <table file>\tBLASTn search result table (fields: qseqid,sseqid,qcovs,qlen,slen,qstart,qend)\n'
help = help + '-org <int>\tTaxid(s) to restrict the database of the BLASTn search\n'
help = help + "-out <path>\tOutput directory (default: output_dir1)\n"
help = help + "-enddist <int>\tMaximum distance between block tip and query tip in base pairs(bp) (default: 50)\n"
help = help + "-minlen <int>\tMinimum element's length in base pairs(bp) (default: 5000)\n"
help = help + "-maxlen <int>\tMaximum element's length in base pairs(bp) (default: 50000)\n"
help = help + "-mincov <int>\tMinimum % query coverage per subject (default: 30)\n"
help = help + "-maxcov <int>\tMaximum % query coverage per subject (default: 90)\n"
help = help + "-cpu <int>\tNumber of threads to execute the blastn search (default: 18)\n"
help = help + "-color <int>\tThe RGB color of the element that is shown by the feature table, three integers between 0 and 255 separated by commas (default: 255,0,0)"
parser = argparse.ArgumentParser(add_help=False)
parser.add_argument('-q')
parser.add_argument('-d')
parser.add_argument('-out')
parser.add_argument('-org')
parser.add_argument('-color')
parser.add_argument('-cpu')
parser.add_argument('-tab')
parser.add_argument('-enddist')
parser.add_argument('-minlen')
parser.add_argument('-maxlen')
parser.add_argument('-mincov')
parser.add_argument('-maxcov')
parser.add_argument('-run')
parser.add_argument('-version', action='store_true')
parser.add_argument('-h', '--help', action='store_true')
args = parser.parse_args()
def is_fasta(filename):
    with open(filename, "r") as handle:
        fasta = SeqIO.parse(handle, "fasta")
        return any(fasta)
def rename(i,name,typ):
  path=''
  if '/' in name:
    path=os.path.split(name)[0]
    name=os.path.split(name)[1]
  newname=os.path.join(path, name)
  if typ=='dir':
    while os.path.isdir(newname):
      i+=1
      newname=os.path.join(path, str(name+"_"+str(i)))
  elif typ=='file':
    while os.path.isfile(newname):
      i+=1
      newname=os.path.join(path, str(name+"_"+str(i)))
  return newname
def blast_parse(tab):
  hits=[]
  qid=[]
  columns=''
  file = open(tab, "r")
  for line in file.readlines():
    line=line.replace("  "," ")
    line=line.replace("\n","")
    if not re.search('^#',line): 
      hits.append(line.split('\t'))
    elif re.search('^# Fields:',line): 
      columns=line.split(': ')[1] 
      columns=columns.replace("\n","").split(', ')
    if re.search('^# Query:',line): 
      qid.append(line.split(': ')[1]) 
  return qid,columns,hits
def validate_args(args):
  valid=True
  if args.q==None:
    print("Missing query file (-q)!")
    valid=False
  elif not args.q==None:
   if not os.path.isfile(args.q):
     print("Query file (-q) not exist!")
     valid=False
   else:
     if is_fasta(args.q):
       param['q']=os.path.realpath(args.q)
     else:
       print("Invalid query file (-q)!")
       valid=False
  if args.run==None and args.tab==None:
    print("Missing the choice (-run) of BLAST search: 'local' or 'web'!")
    valid=False
  elif not args.run==None:
    try:
      str(args.run).lower()
    except Exception:
      traceback.print_exc(file=sys.stdout)
      print("ERROR: BLASTn run choice (-run) must be string: 'local' or 'web'!\n")
      valid=False
    else:
      if str(args.run).lower()=='web':
        param['run']='web'
      elif str(args.run).lower()=='local':
        param['run']='local'
      else:
        print("Invalid BLASTn run choice (-run), must be 'local' or 'web'!")
        valid=False
  if str(args.run).lower()=='local': 
    if args.d==None:
      print("Missing database file (-d)!")
      valid=False
    else:
      if not os.path.isfile(args.d):
        print("Database file (-d) not exist!")
        valid=False
      else:
        if is_fasta(args.d):
          param['d']=os.path.realpath(args.d)
        else:
          print("Invalid database file (-d), invalid formatting!")
          valid=False
  if valid==True and not args.run==None:
    if str(args.run).lower()=='local':
      if args.cpu==None:
        param['cpu']=18
      else:
        try:
          int(args.cpu)
        except Exception:
          traceback.print_exc(file=sys.stdout)
          print("ERROR: Number of threads (-cpu) is not integer!\n")
          valid=False
        else:
          if int(args.cpu)>=1:
            param['cpu']=int(args.cpu)
          else:
            print("Number of threads (-cpu) must be greater than or equal to 1!")
            valid=False
  if valid==True:
    if not args.tab==None:
      if not os.path.isfile(args.tab):
        print("BLASTn table file (-tab) not exist!")
        valid=False
      else:
        missing=[]
        qid,columns,hits=blast_parse(args.tab)
        if qid==[] or columns==[]:
          print("Invalid BLASTn table (-tab)!")
          valid=False
        else:
          for col in ['query id', 'subject id', '% query coverage per subject', 'query length','q. start', 'q. end']:
            if not col in columns:
              valid=False
              missing.append(col)                
        if valid==False: 
          if not missing: 
            pass 
          elif len(missing)==1:
            print("Column {} is not in table (-tab) {}!".format(missing[0],args.tab))
          else:
            print("Columns {} is not in table (-tab) {}!".format(','.join(missing),args.tab))
        else:
          param['columns']=columns
          param['qid']=qid
          param['hits']=hits
          param['tab']=os.path.realpath(args.tab)
  if valid==True and not args.org==None:
    try:
      taxids=[]
      for id_val in args.org.split(','): 
        int(id_val)
        taxids.append("txid{}[ORGN]".format(id_val))
    except Exception:
      traceback.print_exc(file=sys.stdout)
      print("ERROR: Taxid(s) must be integers separated by commas!\n")
      valid=False
    else:
      param['org']=' OR '.join(taxids)
  if valid==True:
    if args.enddist==None:
      param['enddist']=50
    else:
      try:
        int(args.enddist)
      except Exception:
        traceback.print_exc(file=sys.stdout)
        print("ERROR: The maximum distance (-enddist) is not integer!\n")
        valid=False
      else:
        param['enddist']=int(args.enddist)
  if valid==True:
    if args.mincov==None and args.maxcov==None:
      param['mincov']=30
      param['maxcov']=90
    elif not args.mincov==None and not args.maxcov==None:
      try:
        int(args.mincov)
      except Exception:
        traceback.print_exc(file=sys.stdout)
        print("Minimum query coverage (-mincov) is not integer!\n")
        valid=False
      try:
        int(args.maxcov)
      except Exception:
        traceback.print_exc(file=sys.stdout)
        print("ERROR: Maximum query coverage (-maxcov) is not integer!\n")
        valid=False
      if valid==True:
        if int(args.mincov)>100 or int(args.mincov)<0:
          print("ERROR: Minimum query coverage (-mincov) must be between 0 and 100!")
          valid=False
        if int(args.maxcov)>100 or int(args.maxcov)<0: 
          print("ERROR: Maximum query coverage (-maxcov) must be between 0 and 100!")
          valid=False
        if valid==True and int(args.mincov)<int(args.maxcov):
          param['mincov']=int(args.mincov)
          param['maxcov']=int(args.maxcov)
        elif valid==True and int(args.mincov)>int(args.maxcov): 
          param['maxcov']=int(args.mincov)
          param['mincov']=int(args.maxcov)
        elif valid==True and int(args.mincov)==int(args.maxcov): 
          print("Minimum query coverage (-mincov) parameter must be smaller than the maximum query coverage (-maxcov)!")
          valid=False
    elif args.mincov==None and not args.maxcov==None:
      try:
        int(args.maxcov)
      except Exception:
        traceback.print_exc(file=sys.stdout)
        print("ERROR: Maximum query coverage (-maxcov) is not integer!\n")
        valid=False
      else: 
        if int(args.maxcov)<=100 and int(args.maxcov)>=0:
          param['maxcov']=int(args.maxcov)
        else:
          print("ERROR: Maximum query coverage (-maxcov) must be between 0 and 100!\n")
          valid=False
    elif args.maxcov==None and not args.mincov==None:
      try:
        int(args.mincov)
      except Exception:
        traceback.print_exc(file=sys.stdout)
        print("ERROR: Minimum query coverage (-mincov) is not integer!\n")
        valid=False
      else: 
        if int(args.mincov)<=100 and int(args.mincov)>=0:
          param['mincov']=int(args.mincov)
        else:
          print("Minimum query coverage (-mincov) must be between 0 and 100!")
          valid=False
  if valid==True:    
    if args.minlen==None and args.maxlen==None:
      param['minlen']=4000
      param['maxlen']=50000
    elif not args.minlen==None and not args.maxlen==None:
      try:
        int(args.minlen)
      except Exception:
        traceback.print_exc(file=sys.stdout)
        print("ERROR: Minimum element length (-minlen) is not integer!\n")
        valid=False
      try:
        int(args.maxlen)
      except Exception:
        traceback.print_exc(file=sys.stdout)
        print("ERROR: Maximum element length (-maxlen) is not integer!\n")
        valid=False
      if valid==True:
        if int(args.minlen)<0:
          print("ERROR: Minimum element length (-minlen) must be greater than or equal to 0!")
          valid=False
        if int(args.maxlen)<0: 
          print("ERROR: Maximum element length (-maxlen) must be greater than or equal to 0!")
          valid=False
        if valid==True and int(args.minlen)<int(args.maxlen):
          param['minlen']=int(args.minlen)
          param['maxlen']=int(args.maxlen)
        elif valid==True and int(args.minlen)>int(args.maxlen): 
          param['maxlen']=int(args.minlen)
          param['minlen']=int(args.maxlen)
        elif valid==True and int(args.minlen)==int(args.maxlen): 
          print("ERROR: Minimum element length (-minlen) must be smaller than the maximum element length (-maxlen)!")
          valid=False
    elif args.minlen==None and not args.maxlen==None:
      try:
        int(args.maxlen)
      except Exception:
        traceback.print_exc(file=sys.stdout)
        print("ERROR: Maximum element length (-maxlen) is not integer!\n")
        valid=False
      else:
        if int(args.maxlen)>=0:
          param['maxlen']=int(args.maxlen)
        else:
          print("Maximum element length (-maxlen) must be greater than or equal to 0!")
          valid=False
    elif args.maxlen==None and not args.minlen==None:
      try:
        int(args.minlen)
      except Exception:
        traceback.print_exc(file=sys.stdout)
        print("ERROR: Minimum element length (-minlen) is not integer!\n")
        valid=False
      else:
        if int(args.minlen)>=0:
          param['minlen']=int(args.minlen)
        else:
          print("Minimum element length (-minlen) must be greater than or equal to 0!")
          valid=False
  if valid==True:
    if args.color==None:
      param['color']=[255,0,0]
    else:
        try:
          if args.color.count(",")==2:
            param['color']=[]
            for num in args.color.split(","):
              n=int(num)
              if n>=0 and n<=255:
                param['color'].append(int(num))
              else:
                print("The RGB color of the element (-color) must be three integers between 0 and 255!")
                valid=False
                break
          else:
            print("The RGB color of the element (-color) must be three integers between 0 and 255!")
            valid=False
        except Exception:
          traceback.print_exc(file=sys.stdout)
          print("ERROR: The RGB color of the element (-color) must be three integers between 0 and 255!\n")
          valid=False
  if valid==True:
    try:
      os.chdir(call) 
      if args.out is None:
          param["out"] = os.path.join(call, 'output_dir')
      else:
          param["out"] = os.path.join(call, args.out) if os.path.basename(args.out) == args.out else os.path.abspath(args.out)
      original_path = param["out"]
      counter = 1 
      
      if os.path.exists(param["out"]):
          counter = 2
          param["out"] = original_path + "_{}".format(counter)
          while os.path.exists(param["out"]):
              counter += 1
              param["out"] = original_path + "_{}".format(counter)
      
      os.makedirs(param["out"])
      
    except Exception:
      traceback.print_exc(file=sys.stdout)
      print("ERROR: Output directory (-out) not valid or could not be created!\n")
      valid=False
    else:
      print("Creating output directory (-out) {}...".format(param["out"]))
      
  return valid,param
def split_fasta(input_file):
    MAX_ALLOWED_BATCHES = 3
    MAX_BP_DONT_SPLIT_THRESHOLD = 200000  
    THRESHOLD_DONT_SPLIT_MAX_SEQS = 50    
    output_dir_parts = None 
    try:
        if not os.path.isfile(input_file):
            if 'log' in globals() or 'log' in locals() and hasattr(log, 'write'):
                log.write("\nERROR: Input file '{}' not found.".format(input_file))
                log.flush()
            raise IOError("Input file '{}' not found.".format(input_file))
        
        all_records = list(SeqIO.parse(open(input_file, "r"), "fasta")) 
        total_num_sequences = len(all_records)
        if total_num_sequences == 0:
            if 'log' in globals() or 'log' in locals() and hasattr(log, 'write'):
                log.write("\nInput file '{}' contains no sequences or is not a valid FASTA file.".format(input_file))
                log.flush()
            return []
        total_bp_all_sequences = sum(len(r.seq) for r in all_records)
        if total_bp_all_sequences == 0 and total_num_sequences > 0:
            if 'log' in globals() or 'log' in locals() and hasattr(log, 'write'):
                log.write("\nInput file '{}' contains sequences with total 0 base pairs. Processing as a single batch.".format(input_file))
                log.flush()
            return [input_file] 
        is_small_enough_by_seq_count = total_num_sequences <= THRESHOLD_DONT_SPLIT_MAX_SEQS
        is_small_enough_by_bp = total_bp_all_sequences <= MAX_BP_DONT_SPLIT_THRESHOLD
        if is_small_enough_by_seq_count and is_small_enough_by_bp:
            if 'log' in globals() or 'log' in locals() and hasattr(log, 'write'):
                log.write("\nInput file '{}' is small enough ({} sequences, {} bp), not splitting.".format(input_file, total_num_sequences, total_bp_all_sequences))
                log.flush()
            return [input_file]
        if not ('param' in globals() or 'param' in locals()) or not param.get("out"):
            err_msg = "\nERROR: param['out'] is not defined. Cannot create query_parts directory."
            if 'log' in globals() or 'log' in locals() and hasattr(log, 'write'): log.write(err_msg); log.flush()
            raise ValueError(err_msg.strip())
        output_dir_parts = os.path.join(param["out"], "query_parts")
        if not os.path.exists(output_dir_parts):
            try:
                os.makedirs(output_dir_parts)
            except OSError as e_mkdir:
                if 'log' in globals() or 'log' in locals() and hasattr(log, 'write'):
                    log.write("\nERROR: Could not create query_parts directory '{}': {}. Splitting aborted.".format(output_dir_parts, str(e_mkdir)))
                    log.flush()
                raise 
        output_files_created = []
        record_idx = 0 
        for batch_num_idx in range(MAX_ALLOWED_BATCHES): 
            batch_num_for_file = batch_num_idx + 1
            if record_idx >= total_num_sequences: 
                break
            num_batches_left_to_create = MAX_ALLOWED_BATCHES - batch_num_idx
            bp_still_to_distribute = sum(len(r.seq) for r in all_records[record_idx:])
            if bp_still_to_distribute == 0 and record_idx < total_num_sequences: 
                 if batch_num_idx == MAX_ALLOWED_BATCHES -1 or all(len(r.seq) == 0 for r in all_records[record_idx:]):
                     target_bp_for_this_batch = 0 
                 else: 
                     target_bp_for_this_batch = 1 
            elif num_batches_left_to_create > 0 :
                 target_bp_for_this_batch = bp_still_to_distribute / float(num_batches_left_to_create) 
            else: 
                 target_bp_for_this_batch = bp_still_to_distribute 
            current_batch_records = []
            current_batch_bp = 0
            while record_idx < total_num_sequences:
                record = all_records[record_idx]
                record_len = len(record.seq)
                if not current_batch_records: 
                    current_batch_records.append(record)
                    current_batch_bp += record_len
                    record_idx += 1
                    continue 
                is_final_created_batch_loop = (batch_num_idx == MAX_ALLOWED_BATCHES - 1)
                if not is_final_created_batch_loop: 
                    if bp_still_to_distribute - current_batch_bp <= MAX_BP_DONT_SPLIT_THRESHOLD * 1.05 and \
                       len(all_records) - record_idx <= THRESHOLD_DONT_SPLIT_MAX_SEQS * 1.5: 
                        pass
                if batch_num_idx == MAX_ALLOWED_BATCHES - 1: 
                    current_batch_records.append(record)
                    current_batch_bp += record_len
                    record_idx += 1
                    continue
                else: 
                    SMALL_BATCH_FILL_FACTOR = 0.75 
                    add_current_record = False
                    if target_bp_for_this_batch == 0: 
                        if record_len == 0: 
                            add_current_record = True
                    elif current_batch_bp + record_len <= target_bp_for_this_batch: 
                        add_current_record = True
                    elif current_batch_bp < target_bp_for_this_batch * SMALL_BATCH_FILL_FACTOR and \
                         record_len <= target_bp_for_this_batch: 
                        add_current_record = True
                    elif current_batch_bp < target_bp_for_this_batch * 0.3 and \
                         current_batch_bp + record_len <= target_bp_for_this_batch * 1.02 : 
                         add_current_record = True
                    if add_current_record:
                        current_batch_records.append(record)
                        current_batch_bp += record_len
                        record_idx += 1
                    else:
                        break 
            if current_batch_records:
                output_path = os.path.join(output_dir_parts, "query_part_{}.fasta".format(batch_num_for_file))
                SeqIO.write(current_batch_records, output_path, "fasta")
                output_files_created.append(output_path)
        if not output_files_created and total_num_sequences > 0:
             if 'log' in globals() or 'log' in locals() and hasattr(log, 'write'):
                log.write("\nNo batch files were created by splitting despite having {} sequences. Returning original file.".format(total_num_sequences)); log.flush()
             if output_dir_parts and os.path.exists(output_dir_parts): shutil.rmtree(output_dir_parts) 
             return [input_file]
        if 'log' in globals() or 'log' in locals() and hasattr(log, 'write'):
            log.write("\nSplit complete: {} batch files created in '{}' from {} total sequences ({} bp).".format(
                len(output_files_created), output_dir_parts, total_num_sequences, total_bp_all_sequences))
            log.flush()
        return output_files_created
    except IOError as e:
        if 'log' in globals() or 'log' in locals() and hasattr(log, 'write'):
            log.write("\nIOERROR during FASTA splitting for file '{}': {}".format(input_file, str(e)))
            log.write(traceback.format_exc())
            log.flush()
        raise 
    except Exception as e:
        if 'log' in globals() or 'log' in locals() and hasattr(log, 'write'):
            log.write("\nERROR: UNEXPECTED ERROR during FASTA splitting for file '{}': {}".format(input_file, str(e)))
            log.write(traceback.format_exc())
            log.flush()
        if output_dir_parts and os.path.exists(output_dir_parts):
            try:
                shutil.rmtree(output_dir_parts)
                if 'log' in globals() or 'log' in locals() and hasattr(log, 'write'):
                    log.write("\nCleaned up query_parts directory '{}' due to error.".format(output_dir_parts))
                    log.flush()
            except OSError as ose: 
                if 'log' in globals() or 'log' in locals() and hasattr(log, 'write'):
                    log.write("\nERROR: Error during cleanup of directory '{}' after splitting error: {}. Directory might persist.".format(output_dir_parts, str(ose)))
                    log.flush()
        raise
def count_unique_queries_from_blast_outputs(list_of_blast_output_files):
    processed_query_ids = set()
    if not isinstance(list_of_blast_output_files, list):
        list_of_blast_output_files = [list_of_blast_output_files]
    for blast_file_path in list_of_blast_output_files:
        if blast_file_path and os.path.exists(blast_file_path):
            try:
                with open(blast_file_path, "r") as f:
                    for line in f:
                        if line.startswith("# Query:"): 
                            try:
                                query_id = line.split(":", 1)[1].strip() 
                                if query_id:
                                    processed_query_ids.add(query_id)
                            except IndexError:
                                if 'log' in globals() or 'log' in locals() and hasattr(log, 'write'):
                                    log.write("\nWARNING: Malformed '# Query:' line in {}: {}".format(blast_file_path, line.strip()))
                                    log.flush()
            except IOError as e_io:
                if 'log' in globals() or 'log' in locals() and hasattr(log, 'write'):
                    log.write("\nERROR: Could not read BLAST output file {} for progress counting: {}".format(blast_file_path, str(e_io)))
                    log.flush()
            except Exception as e_gen: 
                if 'log' in globals() or 'log' in locals() and hasattr(log, 'write'):
                    log.write("\nERROR: Unexpected error processing file {} for progress: {}".format(blast_file_path, str(e_gen)))
                    log.write(traceback.format_exc())
                    log.flush()
    return len(processed_query_ids)
def run_blast_batch(batch_file, batch_index, param, override_output=None):
    is_original_query = os.path.abspath(batch_file) == os.path.abspath(param.get('q', ''))
    batch_suffix = "original" if is_original_query else "batch_%d" % (batch_index + 1)
    if not batch_suffix == 'original':
        if hasattr(log, 'write'): 
            log.write("\nRunning BLAST for %s (file: %s)..." % (batch_suffix, batch_file))
            log.flush()
    batch_output_path = ""
    if override_output:
        batch_output_path = override_output
    else:
        blastn_parts_output_dir = os.path.join(param["out"], "blastn_parts")
        if not os.path.exists(blastn_parts_output_dir):
            try:
                os.makedirs(blastn_parts_output_dir)
            except OSError as e_mkdir:
                if hasattr(log, 'write'):
                    log.write("\nERROR: Could not create blastn_parts directory '%s': %s. Batch output may fail." % (blastn_parts_output_dir, str(e_mkdir)))
                    log.flush()
        batch_output_path = os.path.join(blastn_parts_output_dir, "blast_batch_%d.tab" % (batch_index + 1))
    if param.get('run') == 'web':
        max_attempts_web = int(param.get('max_web_attempts', 10)) 
        initial_wait_time_web = int(param.get('web_initial_wait', 60)) 
        max_wait_time_web = int(param.get('web_max_wait', 600)) 
        for attempt in range(max_attempts_web):
            try:
                blastn_cline_obj = NcbiblastnCommandline(
                    query=batch_file,
                    db='nt',
                    outfmt="'7 qseqid sseqid qcovs qlen slen qstart qend evalue'",
                    out=batch_output_path,
                    remote=True,
                    task='megablast',
                    max_target_seqs=100 
                )
                if 'org' in param and param['org']: 
                    blastn_cline_obj.entrez_query = "'%s'" % param['org']
                cmd_str = str(blastn_cline_obj)
                if hasattr(log, 'write'):
                    log.write("\n[%s] Running Web BLAST (attempt %d/%d): %s" % (
                        batch_suffix, attempt + 1, max_attempts_web, cmd_str))
                    log.flush()
                blast_errors_log_path = os.path.join(param["out"], "blast_errors.log")
                with open(os.devnull, 'w') as dev_null_fp: 
                    with open(blast_errors_log_path, "a") as err_log_fp: 
                        err_log_fp.write("\n--- Attempt %d/%d for %s at %s ---\n" % (
                            attempt + 1, max_attempts_web, batch_suffix, time.ctime()))
                        err_log_fp.flush()
                        subprocess.check_call(cmd_str, shell=True, stdout=dev_null_fp, stderr=err_log_fp)
                if hasattr(log, 'write'):
                    log.write("\nFinished Web BLAST for %s (file: %s). Output: %s" % (
                        batch_suffix, batch_file, batch_output_path))
                    log.flush()
                return batch_output_path 
            except subprocess.CalledProcessError as e_blast:
                if hasattr(log, 'write'):
                    log.write("\nERROR: [%s] Web BLAST failed (attempt %d/%d), return code %d. See %s for details." % (
                        batch_suffix, attempt + 1, max_attempts_web, e_blast.returncode, blast_errors_log_path))
                    log.write("\n%s" % traceback.format_exc())
                    log.flush()
                if attempt < max_attempts_web - 1:
                    wait_time = min(max_wait_time_web, (initial_wait_time_web * (2 ** attempt)) + random.uniform(0, 15))
                    
                    sys.stderr.write("Retrying Web BLAST for [%s] in %.2f seconds... (Attempt %d/%d after error)\n" % (
                        batch_suffix, wait_time, attempt + 2, max_attempts_web))
                    if hasattr(log, 'write'):
                        log.write("\nRetrying [%s] in %.2f seconds..." % (batch_suffix, wait_time))
                        log.flush()
                    time.sleep(wait_time)
                else:
                    if hasattr(log, 'write'):
                        log.write("\nERROR: [%s] Web BLAST permanently failed after %d retries." % (
                            batch_suffix, max_attempts_web))
                        log.flush()
                    return None 
            except Exception as e_general:
                if hasattr(log, 'write'):
                    log.write("\nERROR: [%s] An unexpected error occurred during Web BLAST (attempt %d/%d): %s" % (
                        batch_suffix, attempt + 1, max_attempts_web, str(e_general)))
                    log.write("\n%s" % traceback.format_exc())
                    log.flush()
                if attempt < max_attempts_web - 1:
                    wait_time = min(max_wait_time_web, (initial_wait_time_web * (1.5 ** attempt)) + random.uniform(0, 10))
                    
                    sys.stderr.write("Retrying Web BLAST for [%s] due to general error in %.2f seconds... (Attempt %d/%d)\n" % (
                        batch_suffix, wait_time, attempt + 2, max_attempts_web))
                    if hasattr(log, 'write'):
                        log.write("\nRetrying [%s] due to general error in %.2f seconds..." % (batch_suffix, wait_time))
                        log.flush()
                    time.sleep(wait_time)
                else:
                    if hasattr(log, 'write'):
                        log.write("\nERROR: [%s] Web BLAST permanently failed after %d retries due to general error." % (
                            batch_suffix, max_attempts_web))
                        log.flush()
                    return None 
        if hasattr(log, 'write'): 
            log.write("\nERROR: [%s] All Web BLAST attempts ultimately failed for file %s." % (batch_suffix, batch_file))
            log.flush()
        return None
    elif param.get('run') == 'local':
        max_attempts_local = 5 
        initial_wait_time_local = 1 
        for attempt in range(max_attempts_local):
            try:
                num_threads_val = param.get('cpu')
                if num_threads_val is None:
                    try:
                        num_threads_val = multiprocessing.cpu_count()
                    except (ImportError, NotImplementedError):
                        num_threads_val = 1
                else:
                    num_threads_val = int(num_threads_val)
                blastn_cline_obj = NcbiblastnCommandline(
                    query=batch_file,
                    db=param.get('d'), 
                    outfmt="'7 qseqid sseqid qcovs qlen slen qstart qend evalue'",
                    out=batch_output_path,
                    num_threads=num_threads_val,
                    task='megablast', 
                    max_target_seqs=100
                )
                cmd_str = str(blastn_cline_obj)
                if hasattr(log, 'write'):
                    log.write("\n[%s] Running Local BLAST (attempt %d/%d): %s" % (
                        batch_suffix, attempt + 1, max_attempts_local, cmd_str))
                    log.flush()
                blast_errors_log_path = os.path.join(param["out"], "blast_errors.log")
                with open(os.devnull, 'w') as dev_null_fp, open(blast_errors_log_path, "a") as err_log_fp:
                    err_log_fp.write("\n--- Local BLAST Attempt %d/%d for %s at %s ---\n" % (
                        attempt + 1, max_attempts_local, batch_suffix, time.ctime()))
                    err_log_fp.flush()
                    subprocess.check_call(cmd_str, shell=True, stdout=dev_null_fp, stderr=err_log_fp)
                if hasattr(log, 'write'):
                    log.write("\nFinished Local BLAST for %s (file: %s). Output: %s" % (
                        batch_suffix, batch_file, batch_output_path))
                    log.flush()
                return batch_output_path 
            except subprocess.CalledProcessError as e_blast:
                if hasattr(log, 'write'):
                    log.write("\nERROR: [%s] Local BLAST failed (attempt %d/%d), return code %d. See %s." % (
                        batch_suffix, attempt + 1, max_attempts_local, e_blast.returncode, blast_errors_log_path))
                    log.write("\n%s" % traceback.format_exc())
                    log.flush()
                if attempt < max_attempts_local - 1:
                    wait_time = (initial_wait_time_local * (2 ** attempt)) + random.uniform(0, 1)
                    
                    sys.stderr.write("Retrying Local BLAST for [%s] in %.2f seconds... (Attempt %d/%d)\n" % (
                        batch_suffix, wait_time, attempt + 2, max_attempts_local))
                    if hasattr(log, 'write'):
                        log.write("\nRetrying [%s] in %.2f seconds..." % (batch_suffix, wait_time))
                        log.flush()
                    time.sleep(wait_time)
                else:
                    if hasattr(log, 'write'):
                        log.write("\nERROR: [%s] Local BLAST permanently failed after %d retries." % (
                            batch_suffix, max_attempts_local))
                        log.flush()
                    return None 
            except Exception as e_general:
                if hasattr(log, 'write'):
                    log.write("\nERROR: [%s] An unexpected error occurred during Local BLAST (attempt %d/%d): %s" % (
                        batch_suffix, attempt + 1, max_attempts_local, str(e_general)))
                    log.write("\n%s" % traceback.format_exc())
                    log.flush()
                if attempt < max_attempts_local - 1:
                    wait_time = (initial_wait_time_local * (2 ** attempt)) + random.uniform(0, 1) 
                    if hasattr(log, 'write'):
                        log.write("\nRetrying [%s] due to general error in %.2f seconds..." % (batch_suffix, wait_time))
                        log.flush()
                    time.sleep(wait_time)
                else:
                    if hasattr(log, 'write'):
                        log.write("\nERROR: [%s] Local BLAST permanently failed after %d retries due to general error." % (
                            batch_suffix, max_attempts_local))
                        log.flush()
                    return None
        if hasattr(log, 'write'): 
            log.write("\nERROR: [%s] All Local BLAST attempts failed for file %s." % (batch_suffix, batch_file))
            log.flush()
        return None
    else:
        if hasattr(log, 'write'):
            log.write("\nERROR: Invalid run type '%s' specified for batch %s." % (param.get('run'), batch_suffix))
            log.flush()
        return None
def blast(param, q_query_filepath, blast_time_start_obj=None): 
    combined_tab_filepath = os.path.join(param["out"], "blastn.tab")
    blastn_parts_dir = os.path.join(param["out"], "blastn_parts")
    if blast_time_start_obj is None:
        blast_time_start_obj = datetime.now()
    batch_files = split_fasta(q_query_filepath) 
    if not batch_files:
        if hasattr(log, 'write'):
            log.write("\nERROR: No batch files created by split_fasta. Aborting BLAST.")
            log.flush()
        return datetime.now() - blast_time_start_obj 
    total_original_sequences_in_query = 0
    try:
        all_original_records = list(SeqIO.parse(open(q_query_filepath, "r"), "fasta")) 
        total_original_sequences_in_query = len(all_original_records)
    except Exception as e_parse:
        if hasattr(log, 'write'):
            log.write("\nERROR: Could not parse original query file '%s' to count sequences for progress: %s" % (q_query_filepath, str(e_parse)))
            log.flush()
    if hasattr(log, 'write'):
        log.write("\nSplit into %d batch file(s) for BLAST. Total original sequences in query: %s." % (
            len(batch_files), str(total_original_sequences_in_query) if total_original_sequences_in_query > 0 else "N/A"))
        log.flush()
    last_console_progress_len = 0
    if len(batch_files) == 1:
        if hasattr(log, 'write'):
            log.write("\nRunning BLAST for a single batch/original file: %s" % batch_files[0])
            log.flush()
        run_blast_batch(batch_files[0], 0, param, override_output=combined_tab_filepath)
        param['tab'] = combined_tab_filepath 
        sequences_processed_so_far = count_unique_queries_from_blast_outputs([combined_tab_filepath])
        progress_message_single = ""
        if total_original_sequences_in_query > 0:
            
            display_processed = sequences_processed_so_far if sequences_processed_so_far > 0 else total_original_sequences_in_query 
            progress_message_single = "Progress: %d/%d original sequences processed. Batch 1/1 completed." % (
                display_processed, total_original_sequences_in_query)
        else:
            progress_message_single = "Batch 1/1 completed (%d queries processed)." % sequences_processed_so_far
        sys.stdout.write(progress_message_single + "\n")
        sys.stdout.flush()
        if hasattr(log, 'write'):
            log.write("\n" + progress_message_single)
            log.flush()
    else: 
        actual_max_workers = 1 
        inter_batch_submission_delay_sec = 0 
        if param.get('run') == 'web':
            max_web_workers_config = int(param.get('max_web_workers', 3))
            actual_max_workers = min(max_web_workers_config, len(batch_files))
            actual_max_workers = max(1, actual_max_workers) 

            if actual_max_workers > 1:
                inter_batch_submission_delay_sec = int(param.get('web_inter_batch_delay', 10))
            else:
                inter_batch_submission_delay_sec = 0
            
            if hasattr(log, 'write'):
                log.write("\nConfigured for Web BLAST: %d batch(es), up to %d parallel job(s), %d s delay between submissions (if parallel)." % (
                    len(batch_files), actual_max_workers, inter_batch_submission_delay_sec if actual_max_workers > 1 else 0))
        elif param.get('run') == 'local':
            try:
                cpu_count_system = multiprocessing.cpu_count()
            except (ImportError, NotImplementedError): 
                cpu_count_system = 1
            threads_per_local_job = int(param.get('cpu', 1)) 
            threads_per_local_job = max(1, threads_per_local_job) 
            if threads_per_local_job >= cpu_count_system and cpu_count_system > 0:
                 max_local_parallel_jobs = 1 
            else:
                 max_local_parallel_jobs = max(1, int(cpu_count_system / threads_per_local_job if threads_per_local_job > 0 else cpu_count_system)) 
            actual_max_workers = min(len(batch_files), int(max_local_parallel_jobs))
            actual_max_workers = max(1, actual_max_workers) 
            if hasattr(log, 'write'):
                log.write("\nConfigured for Local BLAST: %d batch(es), up to %d parallel job(s), each BLASTN using %d thread(s)." % (
                    len(batch_files), actual_max_workers, threads_per_local_job))
        else:
            if hasattr(log, 'write'):
                log.write("\nERROR: Invalid 'run' mode specified: %s" % param.get('run'))
                log.flush()
            return datetime.now() - blast_time_start_obj 
        if hasattr(log, 'write'): log.flush()
        if not os.path.exists(blastn_parts_dir):
            try:
                os.makedirs(blastn_parts_dir)
            except OSError as e_mkdir:
                if hasattr(log, 'write'):
                    log.write("\nERROR: Could not create blastn_parts directory '%s': %s." % (blastn_parts_dir, str(e_mkdir)))
                    log.flush() 
        result_files_from_futures = []
        futures_map = {} 
        executor = ProcessPoolExecutor(max_workers=actual_max_workers)
        errors_occurred_in_pool = False
        completed_batch_count_in_pool = 0
        try:
            for i, batch_f_input in enumerate(batch_files):
                if param.get('run') == 'web' and actual_max_workers > 1 and i > 0: 
                    if inter_batch_submission_delay_sec > 0:
                        if hasattr(log, 'write'):
                            log.write("\nWaiting %d s before submitting next web BLAST batch (%d/%d)..." % (
                                inter_batch_submission_delay_sec, i + 1, len(batch_files)))
                            log.flush()
                        time.sleep(inter_batch_submission_delay_sec)
                future = executor.submit(run_blast_batch, batch_f_input, i, param)
                futures_map[future] = {'file': os.path.basename(batch_f_input), 'index': i}
            for future_item in as_completed(futures_map.keys()): 
                batch_info = futures_map[future_item]
                batch_filename_for_log = batch_info['file']
                batch_idx_for_log = batch_info['index']
                completed_batch_count_in_pool += 1
                try:
                    result_blast_output_path = future_item.result() 
                    if result_blast_output_path and os.path.exists(result_blast_output_path):
                        if result_blast_output_path not in result_files_from_futures:
                             result_files_from_futures.append(result_blast_output_path)
                        if hasattr(log, 'write'):
                            log.write("\nBLAST batch %d (%s) completed successfully. Output: %s" % (
                                batch_idx_for_log + 1, batch_filename_for_log, result_blast_output_path))
                    else:
                        if hasattr(log, 'write'):
                            log.write("\nERROR: BLAST batch %d (%s) seems to have failed or produced no output file (path: '%s')." % (
                                batch_idx_for_log + 1, batch_filename_for_log, str(result_blast_output_path)))
                        errors_occurred_in_pool = True
                except Exception as e_future: 
                    if last_console_progress_len > 0: 
                        sys.stdout.write('\b' * last_console_progress_len + ' ' * last_console_progress_len + '\b' * last_console_progress_len)
                        last_console_progress_len = 0
                    sys.stdout.write("\n") 
                    sys.stdout.flush()
                    if hasattr(log, 'write'):
                        log.write("\nERROR: Processing for BLAST batch %d (%s) failed in pool: %s" % (
                            batch_idx_for_log + 1, batch_filename_for_log, str(e_future)))
                        log.write("\n%s" % traceback.format_exc()) 
                    errors_occurred_in_pool = True
                finally:
                    if hasattr(log, 'write'): log.flush()
                    current_part_files_for_progress = []
                    if os.path.exists(blastn_parts_dir): 
                        try:
                            current_part_files_for_progress = [os.path.join(blastn_parts_dir, f_name)
                                                               for f_name in os.listdir(blastn_parts_dir)
                                                               if f_name.startswith("blast_batch_") and f_name.endswith(".tab")]
                        except OSError as e_listdir:
                            if hasattr(log, 'write'):
                                log.write("\nWARNING: Could not list blastn_parts directory '%s' for progress: %s" % (blastn_parts_dir, str(e_listdir)))
                                log.flush()
                    sequences_processed_so_far_pool = count_unique_queries_from_blast_outputs(current_part_files_for_progress)
                    progress_line_pool = ""
                    if total_original_sequences_in_query > 0:
                        progress_line_pool = "Progress: %d/%d original sequences processed. Batch %d/%d jobs processed." % (
                            sequences_processed_so_far_pool, total_original_sequences_in_query,
                            completed_batch_count_in_pool, len(batch_files))
                    else:
                        progress_line_pool = "Batch %d/%d jobs processed (%d unique queries from parts)." % (
                            completed_batch_count_in_pool, len(batch_files), sequences_processed_so_far_pool)
                    if last_console_progress_len > 0: 
                        sys.stdout.write('\b' * last_console_progress_len + ' ' * last_console_progress_len + '\b' * last_console_progress_len)
                    sys.stdout.write(progress_line_pool)
                    sys.stdout.flush()
                    last_console_progress_len = len(progress_line_pool)
        finally:
            executor.shutdown(wait=True) 
            if last_console_progress_len > 0: 
                sys.stdout.write("\n")
                sys.stdout.flush()
        if errors_occurred_in_pool:
            if hasattr(log, 'write'):
                log.write("\nERROR: One or more BLAST batch jobs failed during parallel execution. Check log for details.")
                log.flush()
        if not result_files_from_futures:
            if hasattr(log, 'write'):
                log.write("\nERROR: No BLAST result files were generated from batch processing.")
                log.flush()
            param['tab'] = None 
            return datetime.now() - blast_time_start_obj
        if hasattr(log, 'write'):
            log.write("\nCombining %d BLAST result files into '%s'..." % (len(result_files_from_futures), combined_tab_filepath))
            log.flush()
        first_file_processed = False
        try:
            with open(combined_tab_filepath, "w") as outfile:
                for res_file_path in sorted(result_files_from_futures): 
                    if res_file_path and os.path.exists(res_file_path):
                        try:
                            with open(res_file_path, "r") as infile: 
                                if not first_file_processed:
                                    for line in infile:
                                        outfile.write(line)
                                    first_file_processed = True
                                else:
                                    for line in infile:
                                        if not line.startswith("#"): 
                                            outfile.write(line)
                        except IOError as e_read_part:
                            if hasattr(log, 'write'):
                                log.write("\nERROR: Could not read temporary BLAST file '%s': %s" % (res_file_path, str(e_read_part)))
                                log.flush()
                        try:
                            os.remove(res_file_path)
                        except OSError as e_os_rm:
                            if hasattr(log, 'write'):
                                log.write("\nWARNING: Could not remove temporary BLAST file '%s': %s" % (res_file_path, str(e_os_rm)))
                                log.flush()
        except IOError as e_io_combine:
            if hasattr(log, 'write'):
                log.write("\nERROR: Could not open/write to combined BLAST output file '%s': %s" % (combined_tab_filepath, str(e_io_combine)))
                log.flush()
            param['tab'] = None 
            return datetime.now() - blast_time_start_obj
        param['tab'] = combined_tab_filepath 
    for folder_name in ['query_parts', 'blastn_parts']:
        folder_path = os.path.join(param["out"], folder_name)
        if os.path.exists(folder_path):
            try:
                shutil.rmtree(folder_path) 
                if hasattr(log, 'write'):
                    log.write("\nCleaned up intermediate folder: %s" % folder_path)
                    log.flush()
            except OSError as e_os_rmtree:
                if hasattr(log, 'write'):
                    log.write("\nWARNING: Could not remove intermediate folder '%s': %s" % (folder_path, str(e_os_rmtree)))
                    log.flush()
    total_time_taken = datetime.now() - blast_time_start_obj
    if hasattr(log, 'write'):
        log.write("\nBLAST process completed in %s" % str(total_time_taken))
        log.flush()
    return total_time_taken
def read_fasta_sequences(fasta_path):
    sequences = []
    with open(fasta_path) as f:
        header = None
        seq_lines = []
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                if header:
                    sequences.append((header, "".join(seq_lines)))
                header = line[1:]
                seq_lines = []
            else:
                seq_lines.append(line)
        if header: 
            sequences.append((header, "".join(seq_lines)))
    return sequences
def write_fasta(sequences, out_path):
    with open(out_path, "w") as f:
        for header, seq in sequences:
            f.write(">{}\n".format(header))
            for i in range(0, len(seq), 60):
                f.write(seq[i:i+60] + "\n")
def merge_blast_results(input_files, output_file): 
    with open(output_file, "w") as outfile:
        header_written = False
        for file_path in input_files: 
            with open(file_path) as infile:
                for line in infile:
                    if line.startswith("#"): 
                        if not header_written: 
                            outfile.write(line)
                    else:
                        outfile.write(line)
            header_written = True 
def get_missing_queries(blast_tab_file, sequences): 
    hit_ids = set()
    with open(blast_tab_file) as f:
        for line in f:
            if line.startswith("#"): 
                continue
            parts = line.strip().split("\t")
            if parts: 
                hit_ids.add(parts[0]) 
    missing = [(h, s) for h, s in sequences if h not in hit_ids]
    if not missing: 
      return False
    else:
      return missing  
def write_missing_queries(missing_seqs, output_path):
    with open(output_path, "w") as f:
        for header, sequence in missing_seqs:
            f.write(">{}\n".format(header))
            for i in range(0, len(sequence), 60):
                f.write(sequence[i:i+60] + "\n")
def join_lists(list1,list2):
  joined=[]
  for i in range(len(list1)): 
    joined.append([list1[i],list2[i]])
  return joined
def assembly(reads):
  reads.sort() 
  m = [] 
  if not reads: 
      return m
  s = reads[0][0]
  max_end = reads[0][1]
  for i in range(1, len(reads)):
      a = reads[i]
      if a[0] > max_end + 1: 
          m.append([s,max_end])
          s = a[0]
          max_end = a[1]
      else: 
          if a[1] > max_end:
              max_end = a[1]
  m.append([s, max_end]) 
  return sorted(m) 
def one_block(contigs,qlen,enddist):
    start_dist=contigs[0][0]-1 
    end_dist=qlen-contigs[0][1] 
    
    
    if start_dist <= enddist and end_dist > enddist: 
      estart=contigs[0][1]+1
      eend=qlen
      elen=eend-estart+1 if eend >= estart else 0 
      return estart,eend,elen
    
    elif start_dist > enddist and end_dist <= enddist: 
      estart=1
      eend=contigs[0][0]-1
      elen=eend-estart+1 if eend >= estart else 0 
      return estart,eend,elen
    
    elif start_dist > enddist and end_dist > enddist:
      return 0,0,0 
    
    elif start_dist <= enddist and end_dist <= enddist:
      return -1,-1,-1 
    return 0,0,0 
def split_list(joined):
  list1=[]
  list2=[]
  for i in joined:
    list1.append(i[0])
    list2.append(i[1])
  return list1,list2
def find_element(start_coords,end_coords): 
  maxdist=0
  estart_final=0 
  eend_final=0   
  elen_final=0   

  if len(end_coords) <= 1 and len(start_coords) <=1 : 
      return 0,0,0


  for k in range(len(end_coords)-1): 
    distance=start_coords[k+1] - end_coords[k] -1 
    if distance > maxdist:
      maxdist = distance
      estart_final = end_coords[k]+1
      eend_final = start_coords[k+1]-1
  
  if estart_final <= eend_final : 
      elen_final = eend_final - estart_final + 1
  else: 
      return 0,0,0

  return estart_final,eend_final,elen_final

def extract(qseqs_content, query_id_to_find, estart, eend): 
  
  
  header_to_find = ">" + query_id_to_find
  seq_start_index = qseqs_content.find(header_to_find)
  if seq_start_index == -1:
      return "" 
  
  
  seq_actual_start = qseqs_content.find('\n', seq_start_index) + 1
  if seq_actual_start == 0: 
      return ""

  
  next_header_start = qseqs_content.find('>', seq_actual_start)
  
  current_seq_block = ""
  if next_header_start == -1: 
      current_seq_block = qseqs_content[seq_actual_start:]
  else:
      current_seq_block = qseqs_content[seq_actual_start:next_header_start]
  
  
  actual_seq = current_seq_block.replace('\n', '')
  
  
  return actual_seq[int(estart-1):int(eend)]


log = None

try:
  if not len(sys.argv)>1:
      print(help)
  elif args.help == True:
      print(help)
  elif args.version == True:
      print(version)
  else:
    
    elements=[]
    valid,param=validate_args(args)
    if valid==False:
      print(help) 
    else:
      print('Valid arguments!')
      
      try:
        current_working_dir = os.getcwd() 
        os.chdir(param["out"]) 
        log=open('file.log','w') 
        os.chdir(current_working_dir) 

        log.write('insertion_finder v{}\n'.format(version))
        log.flush()
      except Exception as error:
        traceback.print_exc(file=sys.stdout)
        print('Log file was not created: {}'.format(error))
        sys.exit(1) 
      
      log.write('(c) 2021. Giuliana Pola & Arthur Gruber\n') 
      log.flush()
      log.write('For more information access: https://github.com/GiulianaPola/insertion_finder\n') 
      log.flush()
      log.write('\nStart time: {}\n'.format(start_time.strftime("%d/%m/%Y, %H:%M:%S")))
      log.flush()
      log.write('\nWorking directory: {}\n'.format(call)) 
      log.flush()
      log.write('\nCommand line: {}\n'.format(' '.join(sys.argv)))
      log.flush()
      user=""
      try:
        user=os.getlogin()
      except Exception:
        try:
          user=os.environ.get('LOGNAME') 
        except Exception:
          try:
            user=os.environ.get('USER') 
          except Exception:
            pass
      if user: 
        log.write('\nUser: {}\n'.format(user))
        log.flush()
      log.write('\nParameters:\n')
      log.flush()
      for arg_name in vars(args): 
          value = getattr(args, arg_name)
          if value is not None and value is not False: 
              log.write("{}={}\n".format(arg_name,value))
              log.flush()
      log.flush() 
      
      
      qseqs_content = open(param['q'], "r").read() 

      if args.tab==None:
        log.write('\nStarting BLASTn search...\n') 
        log.flush()
        
        blast_time=blast(param, param['q']) 
        print('BLASTn search execution time: {}'.format(blast_time))
        log.write('\nBLASTn search execution time: {}\n'.format(blast_time)) 
        log.flush()
        
        if param.get('tab') and os.path.exists(param['tab']):
            param['qid'],param['columns'],param['hits']=blast_parse(param['tab'])
        else:
            log.write("\nERROR: BLAST output table not found or not generated: {}\n".format(param.get('tab')))
            log.flush()
            
            sys.exit("BLAST failed to produce an output table.")

      else: 
        original_tab_file=param['tab'] 
        sequences = read_fasta_sequences(param['q']) 
        
        
        missing_seqs_data = get_missing_queries(original_tab_file, sequences) 
        
        if missing_seqs_data: 
          missing_queries_fasta_path = os.path.join(param["out"],"missing_queries.fasta")
          write_missing_queries(missing_seqs_data, missing_queries_fasta_path)
          
          log.write("\nFound missing queries. Running BLAST for them...\n")
          log.flush()
          
          
          
          
          original_run_type = param.get('run')
          param['run']='web' 
          
          blast_time_missing = blast(param, missing_queries_fasta_path) 
          
          
          if original_run_type: param['run'] = original_run_type
          else: param.pop('run', None)


          log.write('\nBLASTn search (missing queries) execution time: {}\n'.format(blast_time_missing))
          log.flush()

          
          
          
          
          new_blast_results_path = os.path.join(param["out"], "blastn.tab") 

          
          final_combined_tab_path = os.path.join(param["out"], "blastn_FINAL.tab")
          with open(final_combined_tab_path, 'w') as dest_final:
              
              if os.path.exists(original_tab_file):
                  with open(original_tab_file, 'r') as src_orig:
                      first_orig = True
                      for line in src_orig:
                          if line.startswith("# Query:") and not first_orig: 
                              continue
                          dest_final.write(line)
                          if line.startswith("# Query:"): first_orig = False
              
              
              
              if os.path.exists(new_blast_results_path):
                  with open(new_blast_results_path, 'r') as src_new:
                      first_new = True
                      for line in src_new:
                          if line.startswith("#"): 
                              
                              
                              
                              if line.startswith("# Query:"):
                                  dest_final.write(line) 
                              elif first_new and not os.path.exists(original_tab_file): 
                                  dest_final.write(line)
                          else: 
                              dest_final.write(line)
                          if line.startswith("# Fields:"): first_new = False


          param['tab'] = final_combined_tab_path 
        
        
        if param.get('tab') and os.path.exists(param['tab']):
             param['qid'],param['columns'],param['hits']=blast_parse(param['tab'])
        else:
            log.write("\nERROR: BLAST table (original or combined) not found: {}\n".format(param.get('tab')))
            sys.exit("BLAST table processing error.")


      if 'hits' not in param:
        print("BLASTn table was not read or is empty!")
        if log: log.write("\nBLASTn table was not read or is empty!\n")
        if log: log.flush()
      else:
        print("Opening BLASTn table...")
        if log: log.write("\nOpening BLASTn table...\n")  
        if log: log.flush()
        cont=0 
        econt=0 
        
        
        elements_summary_path = os.path.join(param["out"], 'elements.txt')
        try:
          
          tabular=open(elements_summary_path,'w')
        except Exception:
          traceback.print_exc(file=sys.stdout)
          print("ERROR: Elements table was not created!\n")
          if log: log.write("\nERROR: Elements table was not created!\n")
          if log: log.flush()
        else: 
          tabular.write('insertion_finder v{}\n'.format(version))
          tabular.write("\nQuery file: {}".format(param["q"]))
          if log: log.write("\n\nQuery file: {}\n".format(param["q"])) 
          
          if args.tab is not None: 
            tabular.write("\nBlastn table file (user-provided): {}".format(args.tab))
            if log: log.write("\nBlastn table file (user-provided): {}\n".format(args.tab))
          
          
          tabular.write("\nBlastn table file (processed): {}".format(param["tab"]))
          if log: log.write("\nBlastn table file (processed): {}\n".format(param["tab"]))


          if 'd' in param: 
            tabular.write("\nDatabase file: {}".format(param["d"]))
            if log: log.write("\nDatabase file: {}\n".format(param["d"]))
          elif param.get('run') == 'web': 
            tabular.write("\nDatabase file: nt")
            if log: log.write("\nDatabase file: nt\n")
          
          
          if args.tab is None and 'blast_time' in locals(): 
            if log: log.write('\nBLASTn search execution time: {}\n'.format(blast_time))
            tabular.write('\nBLASTn search execution time: {}\n'.format(blast_time))
          
          if 'org' in param:
            tabular.write("\nTaxids: {}".format(args.org)) 
            if log: log.write("\nTaxids: {}\n".format(args.org))
          
          tabular.write("\nElement length: {}-{}".format(param["minlen"],param["maxlen"]))
          if log: log.write("\nElement length: {}-{}\n".format(param["minlen"],param["maxlen"]))
          
          tabular.write("\nQuery coverage: {}-{}".format(param["mincov"],param["maxcov"]))
          if log: log.write("\nQuery coverage: {}-{}\n".format(param["mincov"],param["maxcov"]))

          tabular.write("\nMaximum block distance: {}".format(param["enddist"]))
          if log: log.write("\nMaximum block distance: {}\n".format(param["enddist"]))

          if 'cpu' in param and param.get('run') == 'local': 
            if log: log.write("\nNumber of threads: {}\n".format(param["cpu"]))
            tabular.write("\nNumber of threads: {}".format(param["cpu"]))
          
          tabular.write("\nElement color: {}\n".format(param["color"])) 
          if log: log.write("\nElement color: {}\n".format(param["color"]))
          if log: log.flush() 

          tabular.write("\nQuery ID\tSubject ID\tElement identification\tElement 5' coordinate\tElement 3' coordinate\tElement length\tValid\n")
          
          
          df = [dict(zip(param['columns'], row)) for row in param['hits']]
          
          
          def to_numeric(value):
              try:
                  return int(value)
              except (ValueError, TypeError):
                  try:
                      return float(value)
                  except (ValueError, TypeError):
                      return value 
          
          for row_idx in range(len(df)):
              for key in df[row_idx].keys():
                  df[row_idx][key] = to_numeric(df[row_idx][key])
          
          
          seen_tuples = set()
          unique_rows = []
          for row_dict in df:
              
              row_tuple = tuple(sorted(row_dict.items()))
              if row_tuple not in seen_tuples:
                  seen_tuples.add(row_tuple)
                  unique_rows.append(row_dict)
          df = unique_rows 

          for qid_val in param['qid']: 
              
              df1 = [row for row in df if row.get('query id') == qid_val]
              
              if log: log.write("\n\nQUERY ID: {}".format(qid_val))
              if log: log.write("\n[{}] {} unique hits!".format(qid_val,len(df1)))
              if log: log.flush()

              if len(df1) == 0:
                  element_status = 'no' 
                  sid_status = 'no hits' 
                  tabular.write("{0}\t{1}\t{2}\t\t\t\t{3}\n".format(qid_val, sid_status, element_status, 'no'))
              else:
                  qlen = df1[0].get('query length', 0) 
                  if log: log.write("\n[{}] Query length: {}".format(qid_val,qlen))
                  if log: log.flush()
                  
                  
                  coverage_by_subject = {}
                  for row in df1:
                      sid = row.get('subject id')
                      coverage_val = row.get('% query coverage per subject') 
                      if coverage_val is None:
                          continue  
                      if sid not in coverage_by_subject:
                          coverage_by_subject[sid] = []
                      try: 
                          coverage_by_subject[sid].append(float(coverage_val))
                      except ValueError:
                          if log: log.write("\nWARNING: Non-numeric coverage for {} - {}: {}".format(qid_val, sid, coverage_val))
                          if log: log.flush()

                  avg_coverage_list = []
                  for sid, values in coverage_by_subject.items():
                      if values: 
                           avg_coverage_list.append({'subject id': sid, '% query coverage per subject': sum(values)/len(values)})
                  
                  avg_coverage_sorted = sorted(
                      avg_coverage_list,
                      key=lambda x: x['% query coverage per subject'],
                      reverse=True
                  )
                  
                  processed_for_qid = False 
                  for i in range(len(avg_coverage_sorted)):
                      if processed_for_qid: break

                      hit_details_log = [] 
                      current_sid = avg_coverage_sorted[i]['subject id']  
                      current_cov = avg_coverage_sorted[i]['% query coverage per subject']
                      
                      hit_details_log.append('\n[{}] Subject id: {}'.format(qid_val,current_sid))
                      hit_details_log.append('[{}] % query coverage per subject: {}'.format(qid_val,current_cov))
                      
                      
                      df2 = [row for row in df if row.get('subject id') == current_sid and row.get('query id') == qid_val]
                      
                      q_starts = [row['q. start'] for row in df2 if 'q. start' in row]
                      q_ends = [row['q. end'] for row in df2 if 'q. end' in row]

                      if not q_starts or not q_ends or len(q_starts) != len(q_ends):
                          if log: log.write("\n[{}] Missing q.start/q.end for subject {}. Skipping.".format(qid_val, current_sid))
                          if log: log.flush()
                          continue

                      
                      reads_for_assembly = sorted([[int(s), int(e)] for s, e in zip(q_starts, q_ends) if s is not None and e is not None])
                      
                      contigs = assembly(reads_for_assembly)
                      hit_details_log.append("[{}] Alignments: {}".format(qid_val,str(contigs)))
                      
                      estart, eend, elen = 0, 0, 0 
                      element_found_this_subject = False

                      if len(contigs)==1:
                        if current_cov >= param['mincov']:
                          estart,eend,elen=one_block(contigs,qlen,param['enddist'])
                          if estart > 0: 
                              element_found_this_subject = True
                          elif estart == 0 and i == len(avg_coverage_sorted)-1 and not processed_for_qid: 
                            if log: log.write("\n".join(str(v) for v in hit_details_log))
                            if log: log.write("\n[{}] Invalid element for single block (type 0: block too far from both ends)".format(qid_val))
                            if log: log.flush()
                          elif estart == -1 and i == len(avg_coverage_sorted)-1 and not processed_for_qid: 
                            if log: log.write("\n".join(str(v) for v in hit_details_log))
                            if log: log.write("\n[{}] Invalid element for single block (type -1: block too close to both ends)".format(qid_val))
                            if log: log.flush()

                        elif current_cov < param['mincov']: 
                           if log: log.write("\n".join(str(v) for v in hit_details_log))
                           if log: log.write('\n[{}] Invalid element, % query coverage less than mincov ({})!'.format(qid_val, param['mincov']))
                           if log: log.flush()
                           
                           if i == len(avg_coverage_sorted)-1 and not processed_for_qid:
                               tabular.write("{0}\t{1}\t{2}\t\t\t\t{3}\n".format(qid_val,current_sid,'no (low cov)','no'))
                               cont+=1 
                               processed_for_qid = True


                      elif len(contigs)>1:
                        if current_cov <= param['maxcov']: 
                          
                          starts_from_contigs = [c[0] for c in contigs]
                          ends_from_contigs = [c[1] for c in contigs]
                          estart,eend,elen = find_element(starts_from_contigs, ends_from_contigs)
                          if elen > 0 : 
                            element_found_this_subject = True
                        
                      
                      if element_found_this_subject:
                          if log: log.write('\n'.join(str(v) for v in hit_details_log))
                          if log: log.write("\n[{}] Element coordinates: {} - {}".format(qid_val, estart,eend))
                          if log: log.write("\n[{}] Element length: {}".format(qid_val,elen))
                          if log: log.flush()
                          
                          cont+=1 
                          processed_for_qid = True 

                          if elen>=param['minlen'] and elen<=param['maxlen']:
                            if log: log.write("\n[{}] Valid element!".format(qid_val))
                            if log: log.flush()
                            econt+=1
                            
                            query_output_dir = os.path.join(param["out"], str(qid_val))
                            if not os.path.exists(query_output_dir): os.makedirs(query_output_dir)
                            
                            tabular.write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\n".format(qid_val,current_sid,'yes',estart,eend,elen,'yes'))
                            
                            
                            ft_path = os.path.join(query_output_dir, str(qid_val+'_element.gb'))
                            try:
                              with open(ft_path,'w') as ft:
                                ft.write("     misc_feature    {0}..{1}\n".format(estart,eend))
                                ft.write("                     /label=element\n")
                                ft.write("                     /color={} {} {}\n".format(param['color'][0],param['color'][1],param['color'][2]))
                                ft.write("ORIGIN\n")
                                
                                
                                
                                
                                target_sequence_string = ""
                                header_to_find_in_qseqs = ">" + str(qid_val) 
                                seq_start_idx = qseqs_content.find(header_to_find_in_qseqs)
                                if seq_start_idx != -1:
                                    seq_actual_start_idx = qseqs_content.find('\n', seq_start_idx) + 1
                                    if seq_actual_start_idx != 0:
                                        next_header_idx = qseqs_content.find('>', seq_actual_start_idx)
                                        seq_block_for_gb = qseqs_content[seq_actual_start_idx : next_header_idx if next_header_idx!=-1 else len(qseqs_content)]
                                        target_sequence_string = seq_block_for_gb.replace('\n','')
                                
                                if not target_sequence_string:
                                    if log: log.write("\n[{}] ERROR: Could not extract sequence for {} to write GB file.\n".format(qid_val, qid_val))
                                else:
                                    
                                    seq_len_for_gb = len(target_sequence_string)
                                    for n_gb in range(0, seq_len_for_gb, 60):
                                        ft.write(str(n_gb+1).rjust(9, ' '))
                                        chunk_gb = target_sequence_string[n_gb : n_gb+60]
                                        for j_gb in range(0, len(chunk_gb), 10):
                                            ft.write(" {}".format(chunk_gb[j_gb:j_gb+10]))
                                        ft.write("\n")
                                ft.write("//\n") 
                              if log: log.write("\n[{}] Writing element's feature table...".format(qid_val))
                              if log: log.flush()
                            except Exception as e_ft:
                              traceback.print_exc(file=sys.stdout)
                              if log: log.write("\n[{}] ERROR writing element's feature table for {}: {}\n".format(qid_val, str(e_ft)))
                              if log: log.flush()

                            
                            fasta_path = os.path.join(query_output_dir, str(qid_val+'_element.fasta'))
                            try:
                              element_seq_data = extract(qseqs_content, str(qid_val), estart, eend)
                              with open(fasta_path,'w') as fasta_out:
                                fasta_out.write(">{0} - element - {1}-{2}\n".format(qid_val,estart,eend))
                                
                                for k_fasta in range(0, len(element_seq_data), 60):
                                    fasta_out.write(element_seq_data[k_fasta:k_fasta+60] + "\n")
                              if log: log.write("\n[{}] Writing element's fasta...".format(qid_val))
                              if log: log.flush()
                            except Exception as e_fasta:
                              traceback.print_exc(file=sys.stdout)
                              if log: log.write("\n[{}] ERROR writing element's fasta for {}: {}\n".format(qid_val, str(e_fasta)))
                              if log: log.flush()
                          else: 
                            if elen < param['minlen']:
                              if log: log.write("\n[{}] Invalid element, smaller than minlen ({} < {}).".format(qid_val, elen, param['minlen']))
                            if elen > param['maxlen']:
                              if log: log.write("\n[{}] Invalid element, larger than maxlen ({} > {}).".format(qid_val, elen, param['maxlen']))
                            if log: log.flush()
                            tabular.write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\n".format(qid_val,current_sid,'yes (invalid size)',estart,eend,elen,'no'))
                      
                      
                      if i == len(avg_coverage_sorted)-1 and not processed_for_qid:
                          if log: log.write("\n".join(str(v) for v in hit_details_log)) 
                          if log: log.write("\n[{}] No valid element found after checking all subjects.".format(qid_val))
                          if log: log.flush()
                          tabular.write("{0}\t{1}\t{2}\t\t\t\t{3}\n".format(qid_val, current_sid if current_sid else 'no valid subjects', 'no', 'no'))
                          cont += 1 
                          processed_for_qid = True
                  
                  if not processed_for_qid and not avg_coverage_sorted : 
                      if log: log.write("\n[{}] No subjects with sufficient data or coverage found.".format(qid_val))
                      if log: log.flush()
                      tabular.write("{0}\t{1}\t{2}\t\t\t\t{3}\n".format(qid_val, 'no valid subjects', 'no', 'no'))
                      cont += 1 
                      processed_for_qid = True


            
          tabular.write('\n\nProcessed {} queries with {} valid elements\n'.format(len(param['qid']),econt)) 
          if log: log.write('\n\nProcessed {} queries with {} valid elements\n'.format(len(param['qid']),econt))
            
          final_time = datetime.now() - start_time
          if log: log.write('\nProgram execution time: {}\n'.format(final_time))
          tabular.write('\nProgram execution time: {}\n'.format(final_time))
          tabular.close()
          if log: log.flush()
          print('Processed {} queries with {} valid elements'.format(len(param['qid']),econt))

  print('Program execution time: {}'.format(datetime.now() - start_time))
  print('\a') 
except Exception:
    traceback.print_exc(file=sys.stdout)
    if log and hasattr(log, 'write') and not log.closed: 
        log.write("\nCRITICAL ERROR: Unhandled exception occurred.\n")
        log.write(traceback.format_exc())
        log.flush()
finally:
    if log and hasattr(log, 'write') and not log.closed:
        log.write("\nEnd of script execution.\n")
        log.flush()
        log.close()
    print("End")