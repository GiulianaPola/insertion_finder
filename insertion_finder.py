#!/usr/bin/env python
# -*- coding: utf-8 -*-
import time
import traceback
import warnings
import sys

def warn_with_traceback(message, category, filename, lineno, file=None, line=None):
    log_output = file if hasattr(file,'write') else sys.stderr
    traceback.print_stack(file=log_output)
    log_output.write(warnings.formatwarning(message, category, filename, lineno, line))
    log_output.flush()

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
version="2.6.1"
help_message = 'insertion_finder v{} - element insertion finder in a genome through a BLAST search\n'.format(version)
help_message = help_message + '(c) 2021. Arthur Gruber & Giuliana Pola\n'
help_message = help_message + 'For the latest version acess: https://github.com/GiulianaPola/insertion_finder\n'
help_message = help_message + 'Usage: insertion_finder.py -q <query file> -d <database file> -run local\n'
help_message = help_message + '\tinsertion_finder.py -q <query file> -d <database file> -run local -tab <BLASTn table file>\n'
help_message = help_message + '\tinsertion_finder.py -q <query file> -run web -tab <BLASTn table file>\n'
help_message = help_message + '\tinsertion_finder.py -q <query file> -run web\n'
help_message = help_message + '\nMandatory parameters:\n'
help_message = help_message + '-q <fasta or multifasta file>\tSequence to search with\n'
help_message = help_message + '-d <multifasta file>\tDatabase to BLAST against (for -run local)\n'
help_message = help_message + '-run <local|web>\tchoice of running local or web BLAST search\n'
help_message = help_message + '\nOptional parameters:\n'
help_message = help_message + '-tab <table file>\tBLASTn search result table (fields: qseqid,sseqid,qcovs,qlen,slen,qstart,qend)\n'
help_message = help_message + '-org <int>\tTaxid(s) to restrict the database of the BLASTn search (for -run web)\n'
help_message = help_message + "-out <path>\tOutput directory (default: output_dir1)\n"
help_message = help_message + "-enddist <int>\tMaximum distance between block tip and query tip in base pairs(bp) (default: 50)\n"
help_message = help_message + "-minlen <int>\tMinimum element's length in base pairs(bp) (default: 4000)\n"
help_message = help_message + "-maxlen <int>\tMaximum element's length in base pairs(bp) (default: 50000)\n"
help_message = help_message + "-mincov <int>\tMinimum % query coverage per subject (default: 30)\n"
help_message = help_message + "-maxcov <int>\tMaximum % query coverage per subject (default: 90)\n"
help_message = help_message + "-cpu <int>\tNumber of threads to execute the local blastn search (default: 18)\n"
help_message = help_message + "-max_web_workers <int>\tMax parallel web BLAST jobs (default: 1. NCBI recommends 1-3 with significant delays. Use >1 with extreme caution.)\n"
help_message = help_message + "-web_inter_batch_delay <int>\tDelay in seconds between submitting batches for web BLAST (default: 15. Increase if penalized or using >1 worker).\n"
help_message = help_message + "-color <int>\tThe RGB color of the element that is shown by the feature table, three integers between 0 and 255 separated by commas (default: 255,0,0)"

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
parser.add_argument('-max_web_workers')
parser.add_argument('-web_inter_batch_delay')
parser.add_argument('-version', action='store_true')
parser.add_argument('-h', '--help', action='store_true')
args = parser.parse_args()

def is_fasta(filename):
    try:
        with open(filename, "r") as handle:
            fasta = SeqIO.parse(handle, "fasta")
            return any(fasta)
    except Exception:
        return False

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

def blast_parse(tab_filepath):
    hits = []
    qids_in_order_of_appearance = []
    qid_set_for_order = set()
    columns_str = ''

    if not os.path.exists(tab_filepath):
        if log and hasattr(log, 'write'):
            log.write("\nERROR: BLAST results file not found for parsing: {}\n".format(tab_filepath))
            log.flush()
        return [], [], []

    try:
        with open(tab_filepath, "r") as file_handle:
            for line in file_handle:
                line = line.replace(u"\xa0 ", " ").strip()
                if not line: continue

                if line.startswith('#'):
                    if line.startswith('# Fields:'):
                        columns_str = line.split(':', 1)[1].strip()
                    elif line.startswith('# Query:'):
                        try:
                            current_qid = line.split(':', 1)[1].strip()
                            if current_qid and current_qid not in qid_set_for_order:
                                qids_in_order_of_appearance.append(str(current_qid))
                                qid_set_for_order.add(str(current_qid))
                        except IndexError:
                            if log and hasattr(log, 'write'):
                                log.write("\nWARNING: Malformed '# Query:' line in {}: {}\n".format(tab_filepath, line))
                                log.flush()
                else:
                    hits.append(line.split('\t'))
                    if not qids_in_order_of_appearance and hits[-1] and hits[-1][0]:
                        first_col_qid = str(hits[-1][0])
                        if first_col_qid and first_col_qid not in qid_set_for_order:
                            qids_in_order_of_appearance.append(first_col_qid)
                            qid_set_for_order.add(first_col_qid)

        parsed_columns = [str(col.strip()) for col in columns_str.split(', ')] if columns_str else []


        if not qids_in_order_of_appearance and hits:
            if log and hasattr(log, 'write'):
                log.write("\nNo '# Query:' lines in {}. Inferring unique query IDs from first column of hits.\n".format(tab_filepath))
                log.flush()
            for hit_row in hits:
                if hit_row and hit_row[0] and str(hit_row[0]) not in qid_set_for_order:
                    qids_in_order_of_appearance.append(str(hit_row[0]))
                    qid_set_for_order.add(str(hit_row[0]))
        
        return qids_in_order_of_appearance, parsed_columns, hits
    except Exception as e_parse:
        if log and hasattr(log, 'write'):
            log.write("\nERROR: Failed to parse BLAST file {}: {}\n{}\n".format(tab_filepath, e_parse, traceback.format_exc()))
            log.flush()
        return [], [], []

def validate_args(args):
    valid=True
    if args.q is None:
        print("Missing query file (-q)!")
        valid=False
    elif not os.path.isfile(args.q):
        print("Query file (-q) '{}' not exist!".format(args.q))
        valid=False
    elif not is_fasta(args.q):
        print("Invalid query file (-q) '{}'! Not a valid FASTA format.".format(args.q))
        valid=False
    else:
        param['q'] = os.path.realpath(args.q)
    if args.run is None and args.tab is None:
        print("Missing the choice (-run) of BLAST search ('local' or 'web'), or a precomputed BLAST table (-tab)!")
        valid=False
    elif args.run is not None:
        run_choice = str(args.run).lower()
        if run_choice == 'web':
            param['run'] = 'web'
        elif run_choice == 'local':
            param['run'] = 'local'
        else:
            print("Invalid BLASTn run choice (-run), must be 'local' or 'web'!")
            valid=False
    if param.get('run') == 'local':
        if args.d is None:
            if args.tab is None:
                print("Missing database file (-d) for local run!")
                valid=False
        elif not os.path.isfile(args.d):
            print("Database file (-d) '{}' not exist!".format(args.d))
            valid=False
        elif not is_fasta(args.d):
            print("Invalid database file (-d) '{}', invalid FASTA formatting!".format(args.d))
            valid=False
        else:
            param['d'] = os.path.realpath(args.d)
    if param.get('run') == 'local':
        if args.cpu is None:
            param['cpu'] = 18
        else:
            try:
                cpu_val = int(args.cpu)
                if cpu_val >= 1:
                    param['cpu'] = cpu_val
                else:
                    print("Number of threads (-cpu) must be greater than or equal to 1!")
                    valid=False
            except ValueError:
                print("ERROR: Number of threads (-cpu) is not integer!\n")
                valid=False
    if param.get('run') == 'web':
        if args.max_web_workers is None:
            param['max_web_workers'] = 1
        else:
            try:
                mww = int(args.max_web_workers)
                if 1 <= mww <= 3:
                    param['max_web_workers'] = mww
                    if mww > 1:
                        print("WARNING: -max_web_workers set to {}. NCBI recommends limiting concurrent requests (e.g., 1-3 with delays). Use with caution and ensure adequate -web_inter_batch_delay.".format(mww))
                else:
                    print("Max web workers (-max_web_workers) should be between 1 and 3 for NCBI web BLAST. Defaulting to 1.")
                    param['max_web_workers'] = 1
            except ValueError:
                print("ERROR: Max web workers (-max_web_workers) must be an integer. Defaulting to 1.")
                param['max_web_workers'] = 1
        if args.web_inter_batch_delay is None:
            param['web_inter_batch_delay'] = 15
        else:
            try:
                wibd = int(args.web_inter_batch_delay)
                if wibd >= 0:
                    param['web_inter_batch_delay'] = wibd
                    if wibd < 10 and param.get('max_web_workers', 1) > 1:
                        print("WARNING: -web_inter_batch_delay is {}s. For >1 worker, NCBI recommends significant delays (e.g., 10-30s or more if penalized).".format(wibd))
                    elif wibd < 5:
                        print("WARNING: -web_inter_batch_delay is {}s. Consider increasing if NCBI warnings persist.".format(wibd))
                else:
                    print("Web inter-batch delay (-web_inter_batch_delay) cannot be negative. Defaulting to 15.")
                    param['web_inter_batch_delay'] = 15
            except ValueError:
                print("ERROR: Web inter-batch delay (-web_inter_batch_delay) must be an integer. Defaulting to 15.")
                param['web_inter_batch_delay'] = 15
    if args.tab is not None:
        if not os.path.isfile(args.tab):
            print("BLASTn table file (-tab) '{}' not exist!".format(args.tab))
            valid=False
        else:
            _, cols, _ = blast_parse(args.tab)
            if not cols:
                print("Invalid or empty BLASTn table (-tab) '{}' or required '# Fields:' line missing!".format(args.tab))
                valid=False
            else:
                param['tab'] = os.path.realpath(args.tab)
    if args.org is not None:
        try:
            taxids=[]
            for id_val in args.org.split(','):
                int(id_val)
                taxids.append("txid{}[ORGN]".format(id_val))
            param['org'] = ' OR '.join(taxids)
        except ValueError:
            print("ERROR: Taxid(s) (-org) must be integers separated by commas!\n")
            valid=False
    if args.enddist is None:
        param['enddist'] = 50
    else:
        try:
            param['enddist'] = int(args.enddist)
            if param['enddist'] < 0:
                print("ERROR: Maximum distance (-enddist) cannot be negative.\n")
                valid = False
        except ValueError:
            print("ERROR: The maximum distance (-enddist) is not an integer!\n")
            valid=False
    mincov_val, maxcov_val = args.mincov, args.maxcov
    if mincov_val is None and maxcov_val is None:
        param['mincov'], param['maxcov'] = 30, 90
    else:
        try:
            p_mincov = int(mincov_val) if mincov_val is not None else 0
            p_maxcov = int(maxcov_val) if maxcov_val is not None else 100
            if not (0 <= p_mincov <= 100):
                print("ERROR: Minimum query coverage (-mincov) must be between 0 and 100!")
                valid = False
            if not (0 <= p_maxcov <= 100):
                print("ERROR: Maximum query coverage (-maxcov) must be between 0 and 100!")
                valid = False
            if valid and p_mincov >= p_maxcov:
                print("ERROR: Minimum query coverage (-mincov) must be strictly smaller than maximum query coverage (-maxcov)!")
                valid = False
            if valid:
                param['mincov'], param['maxcov'] = p_mincov, p_maxcov
        except ValueError:
            print("ERROR: Query coverage (-mincov/-maxcov) values must be integers.")
            valid = False
    minlen_val, maxlen_val = args.minlen, args.maxlen
    if minlen_val is None and maxlen_val is None:
        param['minlen'], param['maxlen'] = 4000, 50000
    else:
        try:
            p_minlen = int(minlen_val) if minlen_val is not None else 0
            p_maxlen = int(maxlen_val) if maxlen_val is not None else float('inf')
            if p_minlen < 0:
                print("ERROR: Minimum element length (-minlen) must be greater than or equal to 0!")
                valid = False
            if p_maxlen < 0:
                print("ERROR: Maximum element length (-maxlen) must be greater than or equal to 0!")
                valid = False
            if valid and p_minlen >= p_maxlen:
                print("ERROR: Minimum element length (-minlen) must be strictly smaller than maximum element length (-maxlen)!")
                valid = False
            if valid:
                param['minlen'], param['maxlen'] = p_minlen, p_maxlen
        except ValueError:
            print("ERROR: Element length (-minlen/-maxlen) values must be integers.")
            valid = False
    if args.color is None:
        param['color'] = [255,0,0]
    else:
        try:
            colors = args.color.split(',')
            if len(colors) == 3:
                param['color'] = []
                for num_str in colors:
                    n = int(num_str)
                    if 0 <= n <= 255:
                        param['color'].append(n)
                    else:
                        print("ERROR: The RGB color of the element (-color) components must be integers between 0 and 255!")
                        valid = False
                        break
                if len(param['color']) != 3 and valid:
                    print("ERROR: The RGB color of the element (-color) must be three integers!")
                    valid = False
            else:
                print("ERROR: The RGB color of the element (-color) must be three integers separated by commas!")
                valid = False
        except ValueError:
            print("ERROR: The RGB color of the element (-color) components must be integers!\n")
            valid=False
    if valid:
        try:
            os.chdir(call)
            if args.out is None:
                param["out"] = os.path.join(call, 'output_dir')
            else:
                param["out"] = os.path.join(call, args.out) if os.path.basename(args.out) == args.out else os.path.abspath(args.out)
            original_path = param["out"]
            counter = 2
            temp_out_path = param["out"]
            if args.out is None:
                if os.path.exists(temp_out_path):
                    temp_out_path = original_path + "_" + str(counter)
                    while os.path.exists(temp_out_path) :
                        counter +=1
                        temp_out_path = original_path + "_" + str(counter)
                    param["out"] = temp_out_path
            else:
                if os.path.exists(temp_out_path):
                    base_path_for_suffix = temp_out_path
                    counter = 2
                    temp_out_path = base_path_for_suffix + "_" + str(counter)
                    while os.path.exists(temp_out_path):
                        counter +=1
                        temp_out_path = base_path_for_suffix + "_" + str(counter)
                    param["out"] = temp_out_path
            os.makedirs(param["out"])
            print("Creating output directory (-out) {}...".format(param['out']))
        except Exception as e_outdir:
            traceback.print_exc(file=sys.stdout)
            print("ERROR: Output directory (-out) not valid or could not be created: {}\n".format(e_outdir))
            valid=False
    return valid,param

def split_fasta(input_file):
    MAX_BP_PER_BATCH = 200000 
    #MAX_SEQS_PER_BATCH = 50
    output_dir_parts = None
    try:
        if log and hasattr(log, 'write'):
            log.write("\nInput file was split according to the limit parameters of {} bp.\n".format(MAX_BP_PER_BATCH)) # or {} sequences  #,MAX_SEQS_PER_BATCH
            log.flush()
        if not os.path.isfile(input_file):
            if log and hasattr(log, 'write'):
                log.write("\nERROR: Input file '{}' not found for splitting.".format(input_file)); log.flush()
            raise IOError("Input file '{}' not found for splitting.".format(input_file))
        all_records = list(SeqIO.parse(open(input_file, "r"), "fasta"))
        total_num_sequences = len(all_records)
        if total_num_sequences == 0:
            if log and hasattr(log, 'write'):
                log.write("\nInput file '{}' contains no sequences. No batches created.".format(input_file)); log.flush()
            return []
        if not param or not param.get("out"):
            err_msg = "\nERROR: param['out'] is not defined. Cannot create query_parts directory."
            if log and hasattr(log, 'write'):
                log.write(err_msg); log.flush()
            raise ValueError(err_msg.strip())
        output_dir_parts = os.path.join(param["out"], "query_parts")
        if not os.path.exists(output_dir_parts):
            try:
                os.makedirs(output_dir_parts)
            except OSError as e_mkdir:
                if log and hasattr(log, 'write'):
                    log.write("\nERROR: Could not create query_parts dir '{}': {}. Splitting aborted.".format(output_dir_parts, str(e_mkdir))); log.flush()
                raise
        output_files_created = []
        record_idx = 0
        batch_num_for_file = 0
        while record_idx < total_num_sequences:
            batch_num_for_file += 1
            current_batch_records = []
            current_batch_bp = 0
            while record_idx < total_num_sequences:
                record = all_records[record_idx]
                record_len = len(record.seq)
                if record_len > MAX_BP_PER_BATCH:
                    if log and hasattr(log, 'write'):
                        log.write("\nWARNING: Sequence {} ({} bp) exceeds MAX_BP_PER_BATCH ({} bp) and will be in its own batch.".format(record.id, record_len, MAX_BP_PER_BATCH))
                        log.flush()
                    if not current_batch_records:
                        current_batch_records.append(record)
                        record_idx += 1
                    break
                if (current_batch_bp + record_len > MAX_BP_PER_BATCH and current_batch_records): #or len(current_batch_records) >= MAX_SEQS_PER_BATCH) 
                    break
                
                current_batch_records.append(record)
                current_batch_bp += record_len
                record_idx += 1
            if current_batch_records:
                output_path = os.path.join(output_dir_parts, "query_part_{}.fasta".format(batch_num_for_file))
                log.write("\nquery_part_{} contains {} bp and {} sequences.\n".format(batch_num_for_file,current_batch_bp,len(current_batch_records)))
                SeqIO.write(current_batch_records, output_path, "fasta")
                output_files_created.append(output_path)
            else:
                if log and hasattr(log, 'write') and record_idx < total_num_sequences:
                    log.write("\nWARNING: Empty batch created during splitting. record_idx: {}, total_num_sequences: {}".format(record_idx, total_num_sequences)); log.flush()
                break
        if log and hasattr(log, 'write'):
            total_bp_all_sequences = sum(len(r.seq) for r in all_records)
            log.write("\nSplit complete: {} batch files created in '{}' from {} total sequences ({} bp).".format(
                len(output_files_created), output_dir_parts, total_num_sequences, total_bp_all_sequences))
            log.flush()
        return output_files_created
    except IOError as e_io_split:
        if log and hasattr(log, 'write'):
            log.write("\nIOERROR during FASTA splitting for file '{}': {}".format(input_file, str(e_io_split)))
            log.write(traceback.format_exc())
            log.flush()
        raise
    except Exception as e_split_generic:
        if log and hasattr(log, 'write'):
            log.write("\nERROR: UNEXPECTED ERROR during FASTA splitting for file '{}': {}".format(input_file, str(e_split_generic)))
            log.write(traceback.format_exc())
            log.flush()
        if output_dir_parts and os.path.exists(output_dir_parts):
            try:
                shutil.rmtree(output_dir_parts)
                if log and hasattr(log, 'write'):
                    log.write("\nCleaned up query_parts directory '{}' due to splitting error.".format(output_dir_parts))
                    log.flush()
            except OSError as ose_cleanup:
                if log and hasattr(log, 'write'):
                    log.write("\nERROR: Error during cleanup of directory '{}' after splitting error: {}. Directory might persist.".format(output_dir_parts, str(ose_cleanup)))
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
                                if log and hasattr(log, 'write'):
                                    log.write("\nWARNING: Malformed '# Query:' line in {}: {}".format(blast_file_path, line.strip()))
                                    log.flush()
            except IOError as e_io_count:
                if log and hasattr(log, 'write'):
                    log.write("\nERROR: Could not read BLAST output file {} for progress counting: {}".format(blast_file_path, str(e_io_count)))
                    log.flush()
            except Exception as e_gen_count:
                if log and hasattr(log, 'write'):
                    log.write("\nERROR: Unexpected error processing file {} for progress: {}".format(blast_file_path, str(e_gen_count)))
                    log.write(traceback.format_exc())
                    log.flush()
    return len(processed_query_ids)

def run_blast_batch(batch_file, batch_index, current_blast_param, override_output_path=None):
    global blastn_log_file_handle
    is_original_full_query = os.path.abspath(batch_file) == os.path.abspath(current_blast_param.get('q_original_for_this_run', ''))
    batch_log_suffix = "original_query" if is_original_full_query else "batch_{}".format(batch_index + 1)
    if log and hasattr(log, 'write') and (not batch_log_suffix == 'original_query' or (hasattr(split_fasta, 'last_split_files_global_ref') and split_fasta.last_split_files_global_ref and len(split_fasta.last_split_files_global_ref) > 1)):
        log.write("\nPreparing BLAST for {} (file: {})...".format(batch_log_suffix, batch_file))
        log.flush()
    actual_batch_output_path = ""
    if override_output_path:
        actual_batch_output_path = override_output_path
    else:
        blastn_parts_dir_for_this_run = current_blast_param.get('current_blastn_parts_dir', os.path.join(current_blast_param["out"], "blastn_parts"))
        if not os.path.exists(blastn_parts_dir_for_this_run):
            try:
                os.makedirs(blastn_parts_dir_for_this_run)
            except OSError as e_mkdir_bparts_run:
                err_msg = "\nERROR: Could not create/access blastn_parts dir '{}': {}.\n".format(blastn_parts_dir_for_this_run, e_mkdir_bparts_run)
                if log and hasattr(log, 'write'): log.write(err_msg); log.flush()
                if blastn_log_file_handle: blastn_log_file_handle.write(err_msg); blastn_log_file_handle.flush()
                return 'FAILURE', None, False
        actual_batch_output_path = os.path.join(blastn_parts_dir_for_this_run, "blast_batch_{}.tab".format(batch_index + 1))
    blastn_common_args = {
        'query': batch_file,
        "outfmt": "'7 qseqid sseqid qcovs qlen slen qstart qend evalue'",
        "out": actual_batch_output_path,
        "task": "megablast",
        "max_target_seqs": 100
    }
    if current_blast_param.get('run') == 'web':
        blastn_cline_obj = NcbiblastnCommandline(db='nt', remote=True, **blastn_common_args)
        if 'org' in current_blast_param and current_blast_param['org']:
            blastn_cline_obj.entrez_query = "'{}'".format(current_blast_param['org'])
    elif current_blast_param.get('run') == 'local':
        num_threads = max(1, int(current_blast_param.get('cpu', 1)))
        blastn_cline_obj = NcbiblastnCommandline(db=current_blast_param.get('d'), num_threads=num_threads, **blastn_common_args)
    else:
        err_msg = "\nERROR: Invalid run type '{}' for batch {}.".format(current_blast_param.get('run'), batch_log_suffix)
        if log: log.write(err_msg); log.flush()
        if blastn_log_file_handle: blastn_log_file_handle.write(err_msg + "\n"); blastn_log_file_handle.flush()
        return 'FAILURE', None, False
    max_attempts = 3
    base_backoff_delay = 15
    attempt = 0
    cpu_warning_detected = False
    while attempt < max_attempts:
        cmd_str = str(blastn_cline_obj)
        attempt_timestamp = time.ctime()
        if blastn_log_file_handle:
            blastn_log_file_handle.write("\n--- {} BLAST Attempt {}/{} for {} (file: {}) at {} ---\nCommand: {}\n".format(
                current_blast_param.get('run').upper(), attempt + 1, max_attempts,
                batch_log_suffix, batch_file, attempt_timestamp, cmd_str))
            blastn_log_file_handle.flush()
        try:
            process = subprocess.Popen(cmd_str, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            stdout_output_bytes, stderr_output_bytes = process.communicate()
            stderr_output = stderr_output_bytes.decode('utf-8', 'replace') if stderr_output_bytes is not None else ""
            return_code = process.returncode
            if stderr_output and blastn_log_file_handle:
                blastn_log_file_handle.write("STDERR Output:\n{}\n".format(stderr_output))
                blastn_log_file_handle.flush()
            if "Searches from this IP address have consumed a large amount of server CPU time" in stderr_output:
                cpu_warning_detected = True
                if blastn_log_file_handle:
                    blastn_log_file_handle.write("NCBI CPU Usage Warning Triggered.\n")
                    blastn_log_file_handle.flush()
            if return_code == 0:
                if blastn_log_file_handle:
                    blastn_log_file_handle.write("SUCCESS: {} BLAST for {} completed.\n".format(current_blast_param.get('run').upper(), batch_log_suffix))
                    blastn_log_file_handle.flush()
                return 'SUCCESS', actual_batch_output_path, cpu_warning_detected
            else:
                if "BLAST query/options error" in stderr_output or "Argument" in stderr_output:
                    if log: log.write("\nERROR: [{}] Permanent BLAST error. Not retrying. Check blastn.log.\n".format(batch_log_suffix))
                    return 'FAILURE', None, cpu_warning_detected
                raise subprocess.CalledProcessError(return_code, cmd_str, output=stderr_output)
        except (subprocess.CalledProcessError, Exception) as e:
            error_detail = e.output if hasattr(e, 'output') and e.output else str(e)
            if blastn_log_file_handle:
                blastn_log_file_handle.write("ERROR: {} BLAST for {} failed (attempt {}/{}), return code {}.\nError Detail: {}\n".format(
                    current_blast_param.get('run').upper(), batch_log_suffix, attempt + 1,
                    max_attempts, getattr(e, 'returncode', 'N/A'), error_detail.strip()))
                blastn_log_file_handle.flush()
            attempt += 1
            if attempt < max_attempts:
                wait_time = (base_backoff_delay * (2 ** (attempt - 1))) + random.uniform(0, 5)
                retry_msg = "Retrying WEB BLAST for [{}] in {:.2f} seconds... (Attempt {}/{})\n".format(
                    batch_log_suffix, wait_time, attempt + 1, max_attempts)
                sys.stderr.write(retry_msg)
                if log: log.write("\n" + retry_msg.strip())
                if blastn_log_file_handle: blastn_log_file_handle.write(retry_msg)
                log.flush()
                blastn_log_file_handle.flush()
                time.sleep(wait_time)
            else:
                final_fail_msg = "\nERROR: [{}] WEB BLAST permanently failed after {} retries.".format(batch_log_suffix, max_attempts)
                if log: log.write(final_fail_msg)
                if blastn_log_file_handle: blastn_log_file_handle.write(final_fail_msg + "\n")
                return 'FAILURE', None, cpu_warning_detected
    return 'FAILURE', None, cpu_warning_detected

split_fasta.last_split_files_global_ref = []

def blast(current_param_for_blast, query_filepath_for_this_blast_run,
          desired_final_output_filename="blastn.tab",
          blast_time_start_obj=None):
    final_combined_tab_for_this_run = os.path.join(current_param_for_blast["out"], desired_final_output_filename)
    blastn_parts_intermediate_dir = os.path.join(current_param_for_blast["out"], "blastn_parts")
    current_param_for_blast['current_blastn_parts_dir'] = blastn_parts_intermediate_dir
    if blast_time_start_obj is None:
        blast_time_start_obj = datetime.now()
    current_param_for_blast['q_original_for_this_run'] = query_filepath_for_this_blast_run
    batch_files_for_this_run = split_fasta(query_filepath_for_this_blast_run)
    split_fasta.last_split_files_global_ref = batch_files_for_this_run
    if not batch_files_for_this_run:
        if log:
            log.write("\nERROR: No batch files created for {}. Aborting BLAST.".format(query_filepath_for_this_blast_run))
            log.flush()
        return datetime.now() - blast_time_start_obj, None
    total_sequences_in_current_blast_input = 0
    try:
        total_sequences_in_current_blast_input = len(list(SeqIO.parse(open(query_filepath_for_this_blast_run, "r"), "fasta")))
    except Exception: pass
    if log:
        log.write("\nBLASTing file: {} (Split into {} compliant batch file(s)).".format(
            query_filepath_for_this_blast_run, len(batch_files_for_this_run)))
        log.flush()
    if len(batch_files_for_this_run) == 1:
        status, single_batch_output_path, _ = run_blast_batch(batch_files_for_this_run[0], 0, current_param_for_blast, override_output_path=final_combined_tab_for_this_run)
        if status != 'SUCCESS':
            if log: log.write("\nERROR: Single batch BLAST failed. See blastn.log.\n"); log.flush()
            return datetime.now() - blast_time_start_obj, None
        final_output_path = single_batch_output_path
    else:
        run_mode = current_param_for_blast.get('run', 'local')
        if log:
            log.write("\nConfigured for {} BLAST: {} batch(es), running sequentially.\n".format(
                run_mode.upper(), len(batch_files_for_this_run)))
            if run_mode == 'web': log.write("Initial inter-batch delay: {}s.".format(current_param_for_blast.get('web_inter_batch_delay', 30)))
            log.flush()
        if not os.path.exists(blastn_parts_intermediate_dir):
            os.makedirs(blastn_parts_intermediate_dir)
        pool_generated_result_files = []
        failed_batches = []
        for i, batch_file_path in enumerate(batch_files_for_this_run):
            if i > 0 and current_param_for_blast.get('run') == 'web':
                current_delay = current_param_for_blast.get('web_inter_batch_delay', 30)
                if log: log.write("\nWaiting {}s before submitting next batch (job {}/{}).".format(current_delay, i + 1, len(batch_files_for_this_run))); log.flush()
                time.sleep(current_delay)
            
            status, result_path, cpu_warning = run_blast_batch(batch_file_path, i, current_param_for_blast)

            if status == 'SUCCESS' and result_path:
                pool_generated_result_files.append(result_path)
                if log: log.write("\nBLAST batch for {} completed.".format(os.path.basename(batch_file_path))); log.flush()
                if cpu_warning:
                    current_param_for_blast['web_inter_batch_delay'] = min(300, current_param_for_blast.get('web_inter_batch_delay', 30) * 2 + 10)
                    if log: log.write("\nWARNING: CPU warning detected. Inter-batch delay increased to {}s for subsequent jobs.".format(current_param_for_blast['web_inter_batch_delay'])); log.flush()
            else:
                failed_batches.append(os.path.basename(batch_file_path))
                if log: log.write("\nERROR: BLAST batch for {} failed. See blastn.log for details. Continuing with next batch.".format(os.path.basename(batch_file_path))); log.flush()
        
        if failed_batches:
             if log:
                log.write("\nWARNING: {} of {} BLAST batches failed. Proceeding with results from successful batches.".format(len(failed_batches), len(batch_files_for_this_run)))
                log.flush()

        if not pool_generated_result_files:
            if log: log.write("\nERROR: No BLAST result files were generated as all batches failed.\n"); log.flush()
            return datetime.now() - blast_time_start_obj, None
        
        with open(final_combined_tab_for_this_run, 'wb') as outfile:
            for i, part_file in enumerate(sorted(pool_generated_result_files)):
                with open(part_file, 'rb') as infile:
                    if i == 0:
                        outfile.writelines(infile)
                    else:
                        for line in infile:
                            if not line.startswith(b'#'):
                                outfile.write(line)
        if os.path.exists(blastn_parts_intermediate_dir): shutil.rmtree(blastn_parts_intermediate_dir)
        query_parts_dir = os.path.join(current_param_for_blast["out"], "query_parts")
        if os.path.exists(query_parts_dir): shutil.rmtree(query_parts_dir)
        final_output_path = final_combined_tab_for_this_run
    time_blast_run = datetime.now() - blast_time_start_obj
    final_output_exists = os.path.exists(final_output_path)
    if log: log.write("\nBLAST process completed in {}. Final output: {}\n".format(
        str(time_blast_run),
        final_output_path if final_output_exists else "Not Generated"))
    return time_blast_run, (final_output_path if final_output_exists else None)

def get_all_query_ids_from_fasta(fasta_filepath):
    all_ids = set()
    try:
        for record in SeqIO.parse(fasta_filepath, "fasta"):
            if record.id:
                all_ids.add(str(record.id)) 
    except IOError:
        if log and hasattr(log, 'write'):
            log.write("\nERROR: Query FASTA file {} not found when trying to get all IDs.\n".format(fasta_filepath))
            log.flush()
    except Exception as e:
        if log and hasattr(log, 'write'):
            log.write("\nERROR: Could not parse query FASTA file {} to get all IDs: {}\n".format(fasta_filepath, e))
            log.flush()
    return all_ids

def get_missing_queries(blast_tab_file, all_expected_query_ids_set):
    present_in_blast_ids = set()
    if not blast_tab_file or not os.path.exists(blast_tab_file):
        if log and hasattr(log, 'write'):
            log.write("\nWARNING: BLAST tab file '{}' not found for missing query check. Assuming all {} expected queries are missing.\n".format(blast_tab_file, len(all_expected_query_ids_set)))
            log.flush()
        return list(all_expected_query_ids_set)

    try:
        with open(blast_tab_file, "r") as f:
            for line in f:
                stripped_line = line.strip()
                if stripped_line.startswith("# Query:"):
                    try:
                        query_id_from_blast = stripped_line.split(":", 1)[1].strip()
                        if query_id_from_blast:
                            present_in_blast_ids.add(str(query_id_from_blast)) 
                    except IndexError:
                        if log and hasattr(log, 'write'): log.write("\nWARNING: Malformed '# Query:' line in {}: {}\n".format(blast_tab_file, stripped_line)); log.flush()
                elif not stripped_line.startswith("#") and "\t" in stripped_line:
                    parts = stripped_line.split("\t")
                    if parts and len(parts) > 0:
                        present_in_blast_ids.add(str(parts[0])) 
    except Exception as e:
        if log and hasattr(log, 'write'):
            log.write("\nERROR reading BLAST table {} for missing query check: {}\n".format(blast_tab_file, e))
            log.flush()
        return list(all_expected_query_ids_set)

    missing_ids = list(all_expected_query_ids_set - present_in_blast_ids)

    return missing_ids if missing_ids else False


def write_missing_queries(missing_query_ids_list, original_fasta_path, output_path_for_missing_fasta):
    records_to_write = []
    if not missing_query_ids_list:
        return False

    missing_query_ids_set = set(str(qid) for qid in missing_query_ids_list) 

    try:
        for record in SeqIO.parse(original_fasta_path, "fasta"):
            if str(record.id) in missing_query_ids_set: 
                records_to_write.append(record)

        if records_to_write:
            SeqIO.write(records_to_write, output_path_for_missing_fasta, "fasta")
            if log and hasattr(log, 'write'): log.write("\nSuccessfully wrote {} missing query sequences to {}\n".format(len(records_to_write), output_path_for_missing_fasta)); log.flush()
            return True
        else:
            if log and hasattr(log, 'write'): log.write("\nWARNING: No sequences found in {} for the {} missing IDs provided. Missing queries FASTA not written.\n".format(original_fasta_path, len(missing_query_ids_list))); log.flush()
            return False
    except IOError:
        if log and hasattr(log, 'write'): log.write("\nERROR: Original query FASTA file {} not found while trying to write missing queries.\n".format(original_fasta_path)); log.flush()
        return False
    except Exception as e:
        if log and hasattr(log, 'write'): log.write("\nERROR writing missing queries FASTA from {} to {}: {}\n".format(original_fasta_path, output_path_for_missing_fasta, e)); log.flush()
        return False

def join_lists(list1,list2):
    joined=[]
    for i in range(len(list1)):
        joined.append([list1[i],list2[i]])
    return joined

def assembly(reads):
    if not reads:
        return []
    processed_reads = []
    for r_idx, r_val in enumerate(reads):
        try:
            if len(r_val) == 2:
                processed_reads.append([int(r_val[0]), int(r_val[1])])
            else:
                if log and hasattr(log, 'write'): log.write("\nWARNING: Invalid read format (not a pair) in assembly data at index {}: {}. Skipping.\n".format(r_idx, r_val)); log.flush()
        except (ValueError, TypeError):
            if log and hasattr(log, 'write'): log.write("\nWARNING: Invalid coordinate (non-integer) in assembly data at index {}: {}. Skipping.\n".format(r_idx, r_val)); log.flush()
            continue

    if not processed_reads: return []

    processed_reads.sort()

    m = []
    if not processed_reads: return m

    s = processed_reads[0][0]
    max_end = processed_reads[0][1]

    for i in range(1, len(processed_reads)):
        a = processed_reads[i]
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
    if not contigs or not contigs[0]: return 0,0,0
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

def find_element(start_coords, end_coords):
    maxdist=0
    estart_final=0
    eend_final=0
    elen_final=0
    if len(end_coords) <= 1 or len(start_coords) <=1 or len(end_coords) != len(start_coords) :
        return 0,0,0
    for k in range(len(end_coords)-1):
        try:
            current_end = int(end_coords[k])
            next_start = int(start_coords[k+1])
            distance = next_start - current_end - 1
            if distance > maxdist:
                maxdist = distance
                estart_final = current_end + 1
                eend_final = next_start - 1
        except (ValueError, TypeError):
            if log and hasattr(log, 'write'):
                log.write("\nWARNING: Non-integer coordinate in find_element: {} or {}\n".format(end_coords[k], start_coords[k+1])); log.flush()
            continue
    if estart_final <= eend_final and estart_final > 0 :
        elen_final = eend_final - estart_final + 1
    else:
        return 0,0,0
    return estart_final,eend_final,elen_final

def extract(qseqs_content_str, query_id_to_find, estart, eend):
    header_to_find = ">" + str(query_id_to_find)
    seq_start_index = qseqs_content_str.find(header_to_find)
    if seq_start_index == -1:
        if log and hasattr(log, 'write'):
            log.write("\nWARNING: Header '{}' not found in qseqs_content for extraction.\n".format(header_to_find)); log.flush()
        return ""
    seq_actual_start = qseqs_content_str.find('\n', seq_start_index) + 1
    if seq_actual_start == 0:
        if log and hasattr(log, 'write'):
            log.write("\nWARNING: Could not find sequence start for '{}' in qseqs_content.\n".format(header_to_find)); log.flush()
        return ""
    next_header_start = qseqs_content_str.find('>', seq_actual_start)
    current_seq_block = ""
    if next_header_start == -1:
        current_seq_block = qseqs_content_str[seq_actual_start:]
    else:
        current_seq_block = qseqs_content_str[seq_actual_start:next_header_start]
    actual_seq = current_seq_block.replace('\n', '').replace('\r','')
    try:
        s = int(estart) - 1
        e = int(eend)
        if s < 0 : s = 0
        if e > len(actual_seq) : e = len(actual_seq)
        if s >= e : return ""
        return actual_seq[s:e]
    except ValueError:
        if log and hasattr(log, 'write'):
            log.write("\nERROR: Invalid start/end coordinates for extraction: {}, {}\n".format(estart, eend)); log.flush()
        return ""

log = None
blastn_log_file_handle = None 

try:
    if not len(sys.argv)>1:
        print(help_message)
        sys.exit(0)
    elif args.help == True:
        print(help_message)
        sys.exit(0)
    elif args.version == True:
        print(version)
        sys.exit(0)
    else:
        elements=[]
        valid,param=validate_args(args)
        if valid==False:
            print("\nInvalid arguments provided. See help below (-h for more details):")
            print(help_message)
            sys.exit(1)
        else:
            print('Valid arguments!')

            try:
                log_file_path = os.path.join(param["out"], 'file.log')
                log=open(log_file_path,'w')
                
                blastn_log_name = os.path.join(param["out"], "blastn.log")
                blastn_log_file_handle = open(blastn_log_name, "a")
                blastn_log_file_handle.write("\n\n===== New Execution Run: {} =====\n".format(datetime.now().strftime("%Y-%m-%d %H:%M:%S"))); blastn_log_file_handle.flush()


                log.write('insertion_finder v{}\n'.format(version)); log.flush()
                log.write('(c) 2021. Giuliana Pola & Arthur Gruber\n'); log.flush()
                log.write('For more information access: https://github.com/GiulianaPola/insertion_finder\n'); log.flush()
                log.write('\nStart time: {}\n'.format(start_time.strftime("%d/%m/%Y, %H:%M:%S"))); log.flush()
                log.write('\nWorking directory (where script was called): {}\n'.format(call)); log.flush()
                log.write('\nOutput directory: {}\n'.format(param["out"])); log.flush()
                log.write('\nCommand line: {}\n'.format(' '.join(sys.argv))); log.flush()
                user=""
                try: user=os.getlogin()
                except Exception:
                    try: user=os.environ.get('LOGNAME')
                    except Exception:
                        try: user=os.environ.get('USER')
                        except Exception: pass
                if user: log.write('\nUser: {}\n'.format(user)); log.flush()

                log.write('\nEffective Parameters Used:\n')
                log.flush()
                for p_name, p_value in param.items():
                    if p_name not in ['hits','columns','qid_original_from_fasta', 'current_blastn_parts_dir']:
                        log.write("{}: {}\n".format(p_name, p_value))
                log.flush()
            except Exception as error_log_setup:
                traceback.print_exc(file=sys.stdout)
                print('Log file could not be created in {}: {}'.format(param.get("out", "unknown_output_dir"), error_log_setup))
                if log and hasattr(log, 'write') and not log.closed: log.close() 
                if blastn_log_file_handle and not blastn_log_file_handle.closed: blastn_log_file_handle.close() 
                sys.exit(1)

            qseqs_content = ""
            try:
                with open(param['q'], "r") as q_file_handle:
                    qseqs_content = q_file_handle.read()
            except Exception as e_read_q:
                err_msg = "\nCRITICAL ERROR: Could not read query FASTA content from {}: {}\n".format(param['q'], e_read_q)
                if log and hasattr(log, 'write'): log.write(err_msg); log.flush()
                sys.stderr.write(err_msg)
                traceback.print_exc(file=sys.stderr)
                sys.exit("Failed to read query file: {}".format(param['q']))

            all_query_ids_from_input_fasta = get_all_query_ids_from_fasta(param['q'])
            if not all_query_ids_from_input_fasta:
                err_msg = "\nCRITICAL ERROR: No query IDs could be read from the input FASTA file: {}. Exiting.\n".format(param['q'])
                if log and hasattr(log, 'write'): log.write(err_msg); log.flush()
                sys.stderr.write(err_msg)
                sys.exit("Could not read query IDs from {}".format(param['q']))
            log.write("\nFound {} unique query IDs in input file {}.\n".format(len(all_query_ids_from_input_fasta), param['q'])); log.flush()
            param['qid_original_from_fasta'] = all_query_ids_from_input_fasta

            initial_blast_tab_path = None
            blast_run_for_initial_table = False

            if args.tab is None:
                log.write('\nStarting initial BLASTn search for all queries...\n'); log.flush()
                initial_blast_duration, generated_initial_tab_path = blast(param, param['q'], desired_final_output_filename="blastn_0.tab")

                print('Initial BLASTn search execution time: {}'.format(initial_blast_duration))
                log.write('\nInitial BLASTn search execution time: {}\n'.format(initial_blast_duration)); log.flush()

                if not generated_initial_tab_path or not os.path.exists(generated_initial_tab_path):
                    err_msg = "\nCRITICAL ERROR: Initial BLAST did not produce an output table. Expected at {}. Exiting.\n".format(generated_initial_tab_path)
                    if log and hasattr(log, 'write'): log.write(err_msg); log.flush()
                    sys.stderr.write(err_msg)
                    sys.exit("Initial BLAST failed to produce an output table.")
                initial_blast_tab_path = generated_initial_tab_path
                blast_run_for_initial_table = True
            else:
                initial_blast_tab_path = param['tab']
                log.write("\nUsing user-provided BLASTn table as initial input: {}\n".format(initial_blast_tab_path)); log.flush()

            log.write("\nChecking for missing queries. Comparing {} input FASTA IDs against BLAST results in: {}\n".format(len(all_query_ids_from_input_fasta), os.path.basename(initial_blast_tab_path))); log.flush()
            missing_query_ids = get_missing_queries(initial_blast_tab_path, all_query_ids_from_input_fasta)

            final_blast_tab_to_analyze = initial_blast_tab_path

            if missing_query_ids:
                log.write("\nFound {} missing queries: {}{}\n".format(len(missing_query_ids), ', '.join(missing_query_ids[:5]), '...' if len(missing_query_ids) > 5 else '')); log.flush()

                missing_queries_temp_fasta = os.path.join(param["out"], "missing_queries_supplemental.fasta")
                if not write_missing_queries(missing_query_ids, param['q'], missing_queries_temp_fasta):
                    log.write("\nERROR: Failed to write missing queries to {}. Analysis will proceed with current BLAST table: {}.\n".format(missing_queries_temp_fasta, initial_blast_tab_path)); log.flush()
                else:
                    log.write("Running supplementary BLAST for {} missing queries (from {})...\n".format(len(missing_query_ids), os.path.basename(missing_queries_temp_fasta))); log.flush()

                    sup_blast_param = param.copy()

                    original_run_mode = param.get('run')
                    if original_run_mode == 'local' and 'd' in param and os.path.exists(param['d']):
                        sup_blast_param['run'] = 'local'
                        sup_blast_param['d'] = param['d']
                        if 'cpu' in param: sup_blast_param['cpu'] = param['cpu']
                    else:
                        sup_blast_param['run'] = 'web'
                        sup_blast_param['max_web_workers'] = param.get('max_web_workers', 1)
                        sup_blast_param['web_inter_batch_delay'] = param.get('web_inter_batch_delay', 15) 


                    log.write("Supplementary BLAST will run in '{}' mode.\n".format(sup_blast_param['run'])); log.flush()

                    sup_blast_duration, supplementary_results_tab_path = blast(sup_blast_param, missing_queries_temp_fasta, desired_final_output_filename="blastn_extra.tab")

                    log.write('\nSupplementary BLASTn search execution time: {}\n'.format(sup_blast_duration)); log.flush()

                    if supplementary_results_tab_path and os.path.exists(supplementary_results_tab_path):
                        final_merged_tab_path_for_analysis = os.path.join(param["out"], "blastn.tab") 
                        log.write("\nMerging initial BLAST results ({}) with supplementary results ({}) into {}\n".format(os.path.basename(initial_blast_tab_path), os.path.basename(supplementary_results_tab_path), os.path.basename(final_merged_tab_path_for_analysis))); log.flush()
                        
                        with open(final_merged_tab_path_for_analysis, 'wb+') as outfile_merged: 
                            header_lines_written = False
                            if os.path.exists(initial_blast_tab_path):
                                with open(initial_blast_tab_path, 'rb') as infile_initial: 
                                    for line_init_bytes in infile_initial:
                                        outfile_merged.write(line_init_bytes)
                                        if line_init_bytes.startswith(b"# Fields:"): header_lines_written = True
                            
                            if os.path.exists(supplementary_results_tab_path): 
                                outfile_merged.seek(0, 2) 
                                if outfile_merged.tell() > 0: 
                                    outfile_merged.seek(-1, 2) 
                                    if outfile_merged.read(1) != b'\n':
                                        outfile_merged.write(b"\n") 
                                    outfile_merged.seek(0, 2)
                                    
                                with open(supplementary_results_tab_path, 'rb') as infile_sup: 
                                    wrote_sup_specific_header = False
                                    for line_sup_bytes in infile_sup:
                                        if line_sup_bytes.startswith(b"#"):
                                            if not header_lines_written and not wrote_sup_specific_header : 
                                                outfile_merged.write(line_sup_bytes)
                                                if line_sup_bytes.startswith(b"# Fields:"): wrote_sup_specific_header = True
                                        else:
                                            outfile_merged.write(line_sup_bytes)

                        final_blast_tab_to_analyze = final_merged_tab_path_for_analysis
                        log.write("Successfully merged BLAST results. Final table for analysis: {}\n".format(final_blast_tab_to_analyze)); log.flush()

                        if os.path.exists(missing_queries_temp_fasta):
                            try: os.remove(missing_queries_temp_fasta)
                            except OSError: log.write("\nWARNING: Could not remove {}\n".format(missing_queries_temp_fasta));log.flush()

                        if supplementary_results_tab_path != final_blast_tab_to_analyze and os.path.exists(supplementary_results_tab_path):
                            try: os.remove(supplementary_results_tab_path)
                            except OSError: log.write("\nWARNING: Could not remove {}\n".format(supplementary_results_tab_path));log.flush()

                        if blast_run_for_initial_table and initial_blast_tab_path != final_blast_tab_to_analyze and \
                           initial_blast_tab_path != args.tab and os.path.exists(initial_blast_tab_path) :
                            try: os.remove(initial_blast_tab_path)
                            except OSError: log.write("\nWARNING: Could not remove {}\n".format(initial_blast_tab_path));log.flush()

                    else:
                        log.write("\nERROR: Supplementary BLAST for missing queries did not produce an output table. Analysis will use results from: {}\n".format(initial_blast_tab_path)); log.flush()
            else:
                log.write("\nNo missing queries found in the initial BLAST results. Using: {}\n".format(os.path.basename(initial_blast_tab_path))); log.flush()
                if not blast_run_for_initial_table and initial_blast_tab_path != os.path.join(param["out"], "blastn.tab"):
                    try:
                        shutil.copy(initial_blast_tab_path, os.path.join(param["out"], "blastn.tab"))
                        final_blast_tab_to_analyze = os.path.join(param["out"], "blastn.tab")
                        log.write("\nCopied user-provided table to {} for consistent analysis path.\n".format(final_blast_tab_to_analyze)); log.flush()
                    except Exception as e_copy:
                        log.write("\nWARNING: Could not copy user-provided table {} to standard name: {}. Using original path.\n".format(initial_blast_tab_path, e_copy)); log.flush()
                elif blast_run_for_initial_table and initial_blast_tab_path != os.path.join(param["out"], "blastn.tab"): 
                    try:
                        os.rename(initial_blast_tab_path, os.path.join(param["out"], "blastn.tab"))
                        final_blast_tab_to_analyze = os.path.join(param["out"], "blastn.tab")
                        log.write("\nRenamed generated BLAST table to {} for consistent analysis path.\n".format(final_blast_tab_to_analyze)); log.flush()
                    except Exception as e_rename:
                         log.write("\nWARNING: Could not rename generated BLAST table {} to standard name: {}. Using original path.\n".format(initial_blast_tab_path, e_rename)); log.flush()


            if not final_blast_tab_to_analyze or not os.path.exists(final_blast_tab_to_analyze):
                err_msg = "\nCRITICAL ERROR: Final BLAST table for parsing ('{}') not found or not generated. Exiting.\n".format(final_blast_tab_to_analyze)
                if log and hasattr(log, 'write'): log.write(err_msg); log.flush()
                sys.stderr.write(err_msg)
                sys.exit("Final BLAST table '{}' not found.".format(final_blast_tab_to_analyze))

            param['tab'] = final_blast_tab_to_analyze

            parsed_qids_from_final_blast, parsed_columns_from_final_blast, parsed_hits_from_final_blast = blast_parse(param['tab'])

            param['qid'] = parsed_qids_from_final_blast
            param['columns'] = parsed_columns_from_final_blast
            param['hits'] = parsed_hits_from_final_blast

            if not param['columns']:
                log.write("\nWARNING: No '# Fields:' line found or parsed in the final BLAST table: {}. Column names will be missing, analysis might fail.\n".format(param['tab']));log.flush()

            if not param['qid'] and param['hits']:
                log.write("\nWARNING: Could not extract any query IDs from BLAST table {} despite it having hits. Please check table format.\n".format(param['tab'])); log.flush()
            elif not param['hits']:
                log.write("\nFinal BLAST table {} contains no hits.\n".format(param['tab'])); log.flush()

            log.write("\nOpening and processing final BLASTn table: {}...\n".format(os.path.basename(param['tab']))); log.flush()
            print("Opening and processing final BLASTn table: {}...".format(os.path.basename(param['tab'])))

            cont = 0
            econt = 0
            
            target_keys = {
                'query_id': 'query id',
                'subject_id': 'subject id',
                'qcovs': '% query coverage per subject',
                'qlen': 'query length',
                'slen': 'subject length',
                'qstart': 'q. start',
                'qend': 'q. end',
                'evalue': 'evalue'
            }
            expected_numeric_cols_in_dict = [target_keys['qcovs'], target_keys['qlen'], target_keys['slen'], target_keys['qstart'], target_keys['qend'], target_keys['evalue']]


            elements_summary_path = os.path.join(param["out"], 'elements.txt')
            try:
                with open(elements_summary_path,'w') as tabular:
                    tabular.write('insertion_finder v{}\n'.format(version))
                    tabular.write("\nQuery file: {}".format(param["q"]))
                    if args.tab is not None and not blast_run_for_initial_table:
                        tabular.write("\nOriginal user-provided Blastn table file: {}".format(args.tab))
                    tabular.write("\nFinal Blastn table file processed: {}".format(param["tab"]))

                    if 'd' in param and param.get('run') == 'local':
                        tabular.write("\nDatabase file (local): {}".format(param["d"]))
                    elif param.get('run') == 'web':
                        tabular.write("\nDatabase file (web): nt")

                    if 'org' in param: tabular.write("\nTaxids: {}".format(args.org if args.org else "N/A"))
                    tabular.write("\nElement length (min-max): {}-{}".format(param["minlen"],param["maxlen"]))
                    tabular.write("\nQuery coverage (% min-max): {}-{}".format(param["mincov"],param["maxcov"]))
                    tabular.write("\nMaximum block distance (enddist): {}".format(param["enddist"]))
                    if 'cpu' in param and param.get('run') == 'local':
                        tabular.write("\nNumber of threads (local BLAST): {}".format(param["cpu"]))
                    if param.get('run') == 'web':
                        tabular.write("\nMax Web Workers: {}".format(param.get('max_web_workers', 'N/A')))
                        tabular.write("\nWeb Inter-Batch Delay (s): {}".format(param.get('web_inter_batch_delay', 'N/A')))

                    tabular.write("\nElement color (RGB): {}\n".format(param["color"]))

                    tabular.write("\nQuery ID\tElement Subject ID\tSubject Avg. % Query Cov\tNo. of Contigs in Subject\tElement identification\tElement 5' coordinate\tElement 3' coordinate\tElement length\tValid\n")

                    df_for_analysis = []
                    
                    if param['columns'] and param['hits']:
                        try:
                            df_for_analysis = []
                            for row in param['hits']:
                                row_dict = {}
                                for i, col_name_from_blast in enumerate(param['columns']):
                                    row_dict[str(col_name_from_blast)] = str(row[i]) if i < len(row) else ""
                                df_for_analysis.append(row_dict)

                            for row_dict_idx in range(len(df_for_analysis)):
                                for key_col in df_for_analysis[row_dict_idx].keys():
                                    if key_col in expected_numeric_cols_in_dict:
                                        val_to_convert = df_for_analysis[row_dict_idx][key_col]
                                        try:
                                            df_for_analysis[row_dict_idx][key_col] = float(val_to_convert)
                                        except (ValueError, TypeError):
                                            try:
                                                df_for_analysis[row_dict_idx][key_col] = int(val_to_convert)
                                            except (ValueError, TypeError):
                                                if log and hasattr(log, 'write'): log.write("\nWARNING: Could not convert field '{}' to numeric for value '{}' in qid '{}'\n".format(key_col, val_to_convert, df_for_analysis[row_dict_idx].get(target_keys['query_id'],'N/A'))); log.flush()
                                                pass
                        except Exception as e_df_create:
                            if log and hasattr(log, 'write'): log.write("\nERROR creating dictionary list from BLAST hits: {}\n".format(e_df_create)); log.flush()
                    elif not param['columns'] and param['hits']:
                        if log and hasattr(log, 'write'): log.write("\nWARNING: No column names found from BLAST table. Cannot process hits as dictionaries. Analysis will be limited.\n"); log.flush()
                    
                    
                    for qid_from_fasta in all_query_ids_from_input_fasta:
                        cont += 1
                        
                        df1 = [row_dict for row_dict in df_for_analysis if str(row_dict.get(target_keys['query_id'])) == str(qid_from_fasta)]

                        if log: log.write("\n\nQUERY ID: {}".format(qid_from_fasta)); log.flush()
                        if log: log.write("\n[{}] Found {} hits for this query in the final BLAST table.".format(qid_from_fasta, len(df1))); log.flush()

                        if not df1:
                            tabular.write("{}\t{}\t{}\t{}\t{}\t\t\t\t{}\n".format(qid_from_fasta, 'no hits in final table','','','no', 'no'))
                            if log: log.write("\n[{}] No BLAST hits found for this query in the final processed table.".format(qid_from_fasta)); log.flush()
                            continue

                        qlen_current_query = 0
                        if df1[0].get(target_keys['qlen']) is not None:
                            try:
                                qlen_current_query = int(float(df1[0][target_keys['qlen']]))
                            except (ValueError, TypeError):
                                if log: log.write("\nWARNING: Could not convert query length '{}' to int for {}\n".format(df1[0][target_keys['qlen']], qid_from_fasta)); log.flush()
                                for hit_rec_for_qlen in df1:
                                    try:
                                        qlen_current_query = int(float(hit_rec_for_qlen.get(target_keys['qlen'],0)))
                                        if qlen_current_query > 0: break
                                    except: pass

                        if qlen_current_query == 0 and log: log.write("\nWARNING: Query length for {} is 0 or could not be determined from BLAST hits.\n".format(qid_from_fasta)); log.flush()
                        if log: log.write("\n[{}] Query length (from BLAST hit '{}'): {}".format(qid_from_fasta, target_keys['qlen'], qlen_current_query)); log.flush()

                        coverage_by_subject_for_qid = {}
                        for hit_row_dict in df1:
                            sid_hit = hit_row_dict.get(target_keys['subject_id'])
                            cov_val_hit = hit_row_dict.get(target_keys['qcovs'])
                            if sid_hit is None or cov_val_hit is None: continue

                            if sid_hit not in coverage_by_subject_for_qid: coverage_by_subject_for_qid[sid_hit] = []
                            try: coverage_by_subject_for_qid[sid_hit].append(float(cov_val_hit))
                            except ValueError:
                                if log: log.write("\nWARNING: Non-numeric '{}' for {} - {}: {}\n".format(target_keys['qcovs'], qid_from_fasta, sid_hit, cov_val_hit)); log.flush()

                        avg_coverage_list_for_qid = []
                        for sid_key, values_cov in coverage_by_subject_for_qid.items():
                            if values_cov: avg_coverage_list_for_qid.append({target_keys['subject_id']: sid_key, target_keys['qcovs']: sum(values_cov)/len(values_cov)})

                        avg_coverage_sorted_for_qid = sorted(avg_coverage_list_for_qid, key=lambda x: x[target_keys['qcovs']], reverse=True)

                        processed_current_qid_for_element = False
                        best_subject_log_details_for_valid_element = [] 
                        subject_analysis_summary = {} 
                        first_subject_processed_details_log = []


                        for subject_details_idx, subject_details in enumerate(avg_coverage_sorted_for_qid):
                            if processed_current_qid_for_element and subject_analysis_summary.get(current_sid_for_analysis, "").startswith("Valid element found"):
                                break 

                            current_sid_for_analysis = subject_details[target_keys['subject_id']]
                            current_avg_cov_for_sid = subject_details[target_keys['qcovs']]
                            reason_for_decision = "N/A - Criteria not met or not processed for element"
                            
                            current_subject_temp_log = [
                                "\n[{}] Subject id: {}".format(qid_from_fasta, current_sid_for_analysis),
                                "[{}] Avg. {} for this subject: {:.2f}".format(qid_from_fasta, target_keys['qcovs'], current_avg_cov_for_sid)
                            ]

                            df2_hits_for_subject = [row for row in df1 if row.get(target_keys['subject_id']) == current_sid_for_analysis]

                            q_starts_subj = [r.get(target_keys['qstart']) for r in df2_hits_for_subject if r.get(target_keys['qstart']) is not None]
                            q_ends_subj = [r.get(target_keys['qend']) for r in df2_hits_for_subject if r.get(target_keys['qend']) is not None]


                            if not q_starts_subj or not q_ends_subj or len(q_starts_subj) != len(q_ends_subj):
                                reason_for_decision = "Missing or mismatched q.start/q.end"
                                subject_analysis_summary[current_sid_for_analysis] = reason_for_decision
                                if subject_details_idx == 0: first_subject_processed_details_log = current_subject_temp_log + ["[{}] Reason: {}".format(qid_from_fasta, reason_for_decision)]
                                continue

                            reads_for_asm_subj = []
                            assembly_coord_error = False
                            for s_asm, e_asm in zip(q_starts_subj, q_ends_subj):
                                try:
                                    s = int(float(s_asm))
                                    e = int(float(e_asm))
                                    reads_for_asm_subj.append([min(s,e), max(s,e)])
                                except (ValueError, TypeError):
                                    reason_for_decision = "Invalid coordinate for assembly"
                                    assembly_coord_error = True
                                    break
                            if assembly_coord_error:
                                subject_analysis_summary[current_sid_for_analysis] = reason_for_decision
                                if subject_details_idx == 0: first_subject_processed_details_log = current_subject_temp_log + ["[{}] Reason: {}".format(qid_from_fasta, reason_for_rejection)]
                                continue
                            
                            contigs_asm = assembly(sorted(reads_for_asm_subj))
                            current_subject_temp_log.append("[{}] Alignments (contigs) for subject {}: {}".format(qid_from_fasta, current_sid_for_analysis, str(contigs_asm)))
                            
                            estart_elem, eend_elem, elen_elem = 0, 0, 0
                            element_identified_in_subject = False

                            if len(contigs_asm) == 1:
                                if param['mincov'] <= current_avg_cov_for_sid :
                                    estart_elem, eend_elem, elen_elem = one_block(contigs_asm, qlen_current_query, param['enddist'])
                                    if estart_elem > 0:
                                        element_identified_in_subject = True
                                        reason_for_decision = "Element identified via one_block: {}-{} (len {})".format(estart_elem, eend_elem, elen_elem)
                                    elif estart_elem == 0 :
                                        reason_for_decision = "Single block did not meet 'one_block' criteria (too far from both ends)."
                                    elif estart_elem == -1 :
                                        reason_for_decision = "Single block did not meet 'one_block' criteria (too close to both ends / full coverage)."
                                else:
                                    reason_for_decision = "Single block, but coverage {:.2f}% < mincov {}%.".format(current_avg_cov_for_sid, param['mincov'])

                            elif len(contigs_asm) > 1:
                                if param['mincov'] <= current_avg_cov_for_sid <= param['maxcov']:
                                    starts_contig_list = [c[0] for c in contigs_asm]
                                    ends_contig_list = [c[1] for c in contigs_asm]
                                    estart_elem, eend_elem, elen_elem = find_element(starts_contig_list, ends_contig_list)
                                    if elen_elem > 0:
                                        element_identified_in_subject = True
                                        reason_for_decision = "Element identified via find_element: {}-{} (len {})".format(estart_elem, eend_elem, elen_elem)
                                    else:
                                        reason_for_decision = "Multiple blocks, but no element found between them."
                                else:
                                    reason_for_decision = "Multiple blocks, but avg. qcovs {:.2f}% not in range ({}-{}).".format(current_avg_cov_for_sid, param['mincov'], param['maxcov'])
                            else: 
                                reason_for_rejection = "No valid contigs formed from alignments."
                                if subject_details_idx == 0 : # Update first subject log if this is the reason
                                     first_subject_processed_details_log = current_subject_temp_log + ["[{}] Reason: {}".format(qid_from_fasta, reason_for_rejection)]


                            subject_analysis_summary[current_sid_for_analysis] = reason_for_decision
                            
                            if element_identified_in_subject:
                                if param['minlen'] <= elen_elem <= param['maxlen']:
                                    if not processed_current_qid_for_element: # First valid element found
                                        best_subject_log_details_for_valid_element = current_subject_temp_log + ["[{}] Analysis Result: {}".format(qid_from_fasta, reason_for_decision)]
                                        best_subject_log_details_for_valid_element.append("[{}] Valid element size found!".format(qid_from_fasta))
                                    
                                    processed_current_qid_for_element = True 
                                    econt += 1
                                    tabular.write("{}\t{}\t{:.2f}\t{}\tyes\t{}\t{}\t{}\tyes\n".format(
                                        qid_from_fasta,current_sid_for_analysis, current_avg_cov_for_sid, len(contigs_asm),
                                        estart_elem,eend_elem,elen_elem))


                                    query_specific_output_dir = os.path.join(param["out"], str(qid_from_fasta))
                                    if not os.path.exists(query_specific_output_dir): os.makedirs(query_specific_output_dir)

                                    ft_filepath = os.path.join(query_specific_output_dir, "{}_element.gb".format(str(qid_from_fasta)))
                                    try:
                                        with open(ft_filepath,'w') as ft_handle:
                                            ft_handle.write("     misc_feature     {}..{}\n".format(estart_elem,eend_elem))
                                            ft_handle.write("                     /label=element\n")
                                            ft_handle.write("                     /color={} {} {}\n".format(param['color'][0],param['color'][1],param['color'][2]))
                                            ft_handle.write("ORIGIN\n")
                                            original_query_seq_for_gb = ""
                                            header_to_find_gb = ">" + str(qid_from_fasta)
                                            seq_start_idx_gb = qseqs_content.find(header_to_find_gb)
                                            if seq_start_idx_gb != -1:
                                                seq_actual_start_gb = qseqs_content.find('\n', seq_start_idx_gb) + 1
                                                if seq_actual_start_gb != 0:
                                                    next_header_gb = qseqs_content.find('>', seq_actual_start_gb)
                                                    seq_block_gb = qseqs_content[seq_actual_start_gb : next_header_gb if next_header_gb!=-1 else len(qseqs_content)]
                                                    original_query_seq_for_gb = seq_block_gb.replace('\n','').replace('\r','')

                                            if original_query_seq_for_gb:
                                                for n_gb_idx in range(0, len(original_query_seq_for_gb), 60):
                                                    ft_handle.write(str(n_gb_idx+1).rjust(9, ' '))
                                                    chunk_val_gb = original_query_seq_for_gb[n_gb_idx : n_gb_idx+60]
                                                    for j_gb_idx in range(0, len(chunk_val_gb), 10):
                                                        ft_handle.write(" {}".format(chunk_val_gb[j_gb_idx:j_gb_idx+10]))
                                                    ft_handle.write("\n")
                                                ft_handle.write("//\n")
                                        if log: best_subject_log_details_for_valid_element.append("[{}] Writing element's feature table: {}".format(qid_from_fasta, ft_filepath));
                                    except Exception as e_ft_write:
                                        if log: best_subject_log_details_for_valid_element.append("[{}] ERROR writing feature table for {}: {}".format(qid_from_fasta, current_sid_for_analysis, e_ft_write));

                                    fasta_element_path = os.path.join(query_specific_output_dir, "{}_element.fasta".format(str(qid_from_fasta)))
                                    try:
                                        element_sequence_data = extract(qseqs_content, str(qid_from_fasta), estart_elem, eend_elem)
                                        if element_sequence_data:
                                            with open(fasta_element_path,'w') as fasta_h:
                                                fasta_h.write(">{0} - element - {1}-{2}\n".format(qid_from_fasta,estart_elem,eend_elem))
                                                for k_seq_idx in range(0, len(element_sequence_data), 60):
                                                    fasta_h.write(element_sequence_data[k_seq_idx:k_seq_idx+60] + "\n")
                                            if log: best_subject_log_details_for_valid_element.append("[{}] Writing element's FASTA: {}".format(qid_from_fasta, fasta_element_path));
                                        else:
                                            if log: best_subject_log_details_for_valid_element.append("[{}] WARNING: Extracted element sequence was empty for {}. FASTA not written.".format(qid_from_fasta, current_sid_for_analysis));
                                    except Exception as e_fasta_write:
                                        if log: best_subject_log_details_for_valid_element.append("[{}] ERROR writing element FASTA for {}: {}".format(qid_from_fasta, current_sid_for_analysis, e_fasta_write));
                                    
                                    if log and best_subject_log_details_for_valid_element: log.write("\n".join(best_subject_log_details_for_valid_element)); log.flush() # Log details of the valid element
                                    break # Found a valid element, move to next query_id

                                else: # Element identified, but size invalid
                                    subject_analysis_summary[current_sid_for_analysis] = "Element found but size invalid: Length {} (not within {}-{})".format(elen_elem, param['minlen'], param['maxlen'])
                                    if not processed_current_qid_for_element: # If this is the first element-like thing we see, even if size is bad
                                        best_subject_log_details = current_subject_temp_log + ["[{}] Element found but size invalid: Length {} (not within {}-{})".format(qid_from_fasta, elen_elem, param['minlen'], param['maxlen'])]
                                    
                                    tabular.write("{}\t{}\t{:.2f}\t{}\tyes (invalid size)\t{}\t{}\t{}\tno\n".format(
                                        qid_from_fasta,current_sid_for_analysis, current_avg_cov_for_sid, len(contigs_asm),
                                        estart_elem,eend_elem,elen_elem))
                                    # Do not set processed_current_qid_for_element to True here if we want to keep searching for a perfectly valid one.
                                    # However, the current logic will break if ANY element is found. This needs to be consistent with desired behavior.
                                    # For now, assume any element identification stops further search for this query.
                                    processed_current_qid_for_element = True # Treat as processed for logging summary.
                                    if log and best_subject_log_details: log.write("\n".join(best_subject_log_details)); log.flush()
                                    break 
                            elif subject_details_idx == 0: # First subject, no element identified
                                reason_for_rejection="No element identified"
                                best_subject_log_details = current_subject_temp_log + ["[{}] Analysis Result: {}".format(qid_from_fasta, reason_for_rejection)]


                        if log and hasattr(log, 'write'):
                            log.write("\n[{}] Subject analysis summary:".format(qid_from_fasta))
                            for subj_id_summary, reason_summary in subject_analysis_summary.items():
                                log.write("\n[{}]  - {}: {}".format(qid_from_fasta, subj_id_summary, reason_summary))
                            log.flush()

                        if not processed_current_qid_for_element: # No valid element found for this query_id after checking all subjects
                            if log and best_subject_log_details: # Log details of the best candidate (last one processed or first one if it failed early)
                                log.write("\n".join(best_subject_log_details)); log.flush()
                            
                            if avg_coverage_sorted_for_qid :
                                first_subject_info = avg_coverage_sorted_for_qid[0]
                                first_subject_id = first_subject_info[target_keys['subject_id']]
                                first_subject_cov = first_subject_info[target_keys['qcovs']]
                                temp_df1_for_contig_count = [r_dict for r_dict in df1 if str(r_dict.get(target_keys['subject_id'])) == str(first_subject_id)] # Filter correctly
                                first_subject_contigs_list = []
                                if temp_df1_for_contig_count:
                                    first_subject_contigs_list = assembly([ [r.get(target_keys['qstart']),r.get(target_keys['qend'])] for r in temp_df1_for_contig_count if r.get(target_keys['qstart']) is not None and r.get(target_keys['qend']) is not None])
                                first_subject_num_contigs = len(first_subject_contigs_list)

                                if log: log.write("\n[{}] No valid element identified. Best candidate subject: {}. Reason: {}".format(qid_from_fasta, first_subject_id, subject_analysis_summary.get(first_subject_id, "Criteria not met"))); log.flush()
                                tabular.write("{}\t{}\t{:.2f}\t{}\tno\t\t\t\tno\n".format(
                                    qid_from_fasta, first_subject_id, first_subject_cov, first_subject_num_contigs
                                    ))
                            else:
                                if log: log.write("\n[{}] No subjects with sufficient hits/coverage to analyze for an element.\n".format(qid_from_fasta)); log.flush()
                                tabular.write("{}\tno subjects met criteria\t\t\tno\t\t\t\tno\n".format(qid_from_fasta))


                    tabular.write('\n\nProcessed {} queries from input FASTA with {} valid elements found and written.\n'.format(cont, econt))
                    tabular.write('Total script execution time: {}\n'.format(datetime.now() - start_time))
                    if log: log.write('\n\nProcessed {} queries from input FASTA with {} valid elements found and written.\n'.format(cont, econt)); log.flush()
                    print('Processed {} queries from input FASTA with {} valid elements found and written.'.format(cont, econt))

            except Exception as e_tabular_write:
                err_msg = "\nCRITICAL ERROR writing to elements.txt or during analysis: {}\n".format(e_tabular_write)
                sys.stderr.write(err_msg)
                traceback.print_exc(file=sys.stderr)
                if log and hasattr(log, 'write'): 
                    log.write(err_msg)
                    log.write(traceback.format_exc())
                    log.flush()
            finally:
                if 'tabular' in locals() and tabular and not tabular.closed:
                    tabular.close()


    final_run_duration = datetime.now() - start_time
    if log and hasattr(log, 'write') and not log.closed:
        log.write('\nTotal script execution time: {}\n'.format(final_run_duration))
        log.flush()
    print('Total script execution time: {}'.format(final_run_duration))

except Exception as main_exception:
    tb_str = traceback.format_exc()
    sys.stderr.write("\nCRITICAL ERROR: Unhandled exception occurred in main script execution.\n")
    sys.stderr.write(tb_str)
    if log and hasattr(log, 'write') and not log.closed:
        log.write("\nCRITICAL ERROR: Unhandled exception occurred in main script execution.\n")
        log.write(tb_str)
        log.flush()
finally:
    if log and hasattr(log, 'write') and not log.closed:
        log.write("\nEnd of script execution.\n")
        log.flush()
        log.close()
    if blastn_log_file_handle and not blastn_log_file_handle.closed:
        blastn_log_file_handle.write("\n===== End of Execution Run: {} =====\n".format(datetime.now().strftime("%Y-%m-%d %H:%M:%S")))
        blastn_log_file_handle.flush()
        blastn_log_file_handle.close()
    print("End")
