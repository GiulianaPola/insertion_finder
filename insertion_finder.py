#!/usr/bin/env python
# -*- coding: utf-8 -*-
import traceback
import warnings
import sys
import os
import re
import argparse
from datetime import datetime
from Bio.Blast.Applications import NcbiblastnCommandline
from Bio import SeqIO
import io # Added for SeqIO.parse on string


def warn_with_traceback(message, category, filename, lineno, file=None, line=None):
    log = file if hasattr(file, 'write') else sys.stderr
    traceback.print_stack(file=log)
    log.write(warnings.formatwarning(message, category, filename, lineno, line))
    log.flush()


warnings.showwarning = warn_with_traceback

start_time = datetime.now()
call = os.path.abspath(os.getcwd())


param = dict()

version = "2.5.1"

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
help = help + "-out <path>\tOutput directory (default: output_dir)\n"
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


def rename(i, name, typ):
    path = ''
    if '/' in name:
        path = os.path.split(name)[0]
        name = os.path.split(name)[1]
    newname = os.path.join(path, name)
    if typ == 'dir':
        while os.path.isdir(newname):
            i += 1
            newname = os.path.join(path, "{}_{}".format(name, i))
    elif typ == 'file':
        while os.path.isfile(newname):
            i += 1
            newname = os.path.join(path, "{}_{}".format(name, i))
    return newname


def blast_parse(tab):
    hits = []
    qid = []
    columns = ''
    file = open(tab, "r")
    for line in file.readlines():
        line = line.replace("  ", " ")
        line = line.replace("\n", "")
        if not re.search('#', line):
            hits.append(line.split('\t'))
        elif re.search('# Fields: ', line):
            columns = line.split('# Fields: ')[-1]
            columns = columns.replace("\n", "").split(', ')
        if re.search('# Query: ', line):
            qid.append(line.split('# Query: ')[-1])
    return qid, columns, hits


def validate_args(args):
    valid = True

    if args.q == None:
        print("ERROR: Missing query file (-q)!")
        valid = False
    elif not args.q == None:
        if not os.path.isfile(args.q):
            print("ERROR: Query file (-q) not exist!")
            valid = False
        else:
            if is_fasta(args.q):
                param['q'] = os.path.realpath(args.q)
            else:
                print("ERROR: Invalid query file (-q)!")
                valid = False

    if args.run == None and args.tab == None:
        print("ERROR: Missing the choice (-run) of BLAST search: 'local' or 'web', or a BLASTn table file (-tab)!")
        valid = False
    elif not args.run == None:
        try:
            run_choice = str(args.run).lower()
        except Exception:
            traceback.print_exc(file=sys.stdout)
            print("ERROR: BLASTn run choice (-run) must be string: 'local' or 'web'!\n")
            valid = False
        else:
            if run_choice == 'web':
                param['run'] = 'web'
            elif run_choice == 'local':
                param['run'] = 'local'
            else:
                print("ERROR: Invalid BLASTn run choice (-run), must be 'local' or 'web'!")
                valid = False

    if 'run' in param and param['run'] == 'local':
        if args.d == None:
            print("ERROR: Missing database file (-d) for local BLAST run!")
            valid = False
        else:
            if not os.path.isfile(args.d):
                print("ERROR: Database file (-d) not exist!")
                valid = False
            else:
                if is_fasta(args.d):
                    param['d'] = os.path.realpath(args.d)
                else:
                    print("ERROR: Invalid database file (-d), invalid formatting!")
                    valid = False

    if valid == True and 'run' in param and param['run'] == 'local':
        if args.cpu == None:
            param['cpu'] = 18
        else:
            try:
                int(args.cpu)
            except Exception:
                traceback.print_exc(file=sys.stdout)
                print("ERROR: Number of threads (-cpu) is not integer!\n")
                valid = False
            else:
                if int(args.cpu) >= 1:
                    param['cpu'] = int(args.cpu)
                else:
                    print("ERROR: Number of threads (-cpu) must be greater than or equal to 1!")
                    valid = False

    if valid == True:
        if not args.tab == None:
            if not os.path.isfile(args.tab):
                print("ERROR: BLASTn table file (-tab) not exist!")
                valid = False
            else:
                missing = []
                qid, columns, hits = blast_parse(args.tab)
                if qid == [] or columns == []:
                    print("ERROR: Invalid BLASTn table (-tab)!")
                    valid = False
                else:
                    required_columns = ['query id', 'subject id', '% query coverage per subject', 'query length', 'subject length', 'q. start', 'q. end', 'evalue']
                    for col in required_columns:
                        if col not in columns:
                            valid = False
                            missing.append(col)
                    if valid == False:
                        if len(missing) == 1:
                            print("ERROR: Column {} is not in table (-tab) {}!".format(missing[0], args.tab))
                        else:
                            print("ERROR: Columns {} are not in table (-tab) {}!".format(', '.join(missing), args.tab))
                    else:
                        param['columns'] = columns
                        param['qid'] = qid
                        param['hits'] = hits
                        param['tab'] = os.path.realpath(args.tab)

    if valid == True and not args.org == None:
        try:
            taxids = []
            for id_str in args.org.split(','):
                taxid = int(id_str)
                taxids.append("txid{}[ORGN]".format(taxid))
        except Exception:
            traceback.print_exc(file=sys.stdout)
            print("ERROR: Taxid(s) must be integers separated by commas!\n")
            valid = False
        else:
            param['org'] = ' OR '.join(taxids)

    if valid == True:
        if args.enddist == None:
            param['enddist'] = 50
        else:
            try:
                param['enddist'] = int(args.enddist)
            except Exception:
                traceback.print_exc(file=sys.stdout)
                print("ERROR: The maximum distance (-enddist) is not integer!\n")
                valid = False

    if valid == True:
        if args.mincov == None and args.maxcov == None:
            param['mincov'] = 30
            param['maxcov'] = 90
        else:
            try:
                min_cov = int(args.mincov) if args.mincov else 0
                max_cov = int(args.maxcov) if args.maxcov else 100
            except Exception:
                traceback.print_exc(file=sys.stdout)
                print("ERROR: Minimum or maximum query coverage (-mincov, -maxcov) is not an integer!\n")
                valid = False
            else:
                if not (0 <= min_cov <= 100):
                    print("ERROR: Minimum query coverage (-mincov) must be between 0 and 100!")
                    valid = False
                if not (0 <= max_cov <= 100):
                    print("ERROR: Maximum query coverage (-maxcov) must be between 0 and 100!")
                    valid = False
                if valid:
                    if min_cov > max_cov:
                        print("Warning: Minimum query coverage (-mincov) is greater than maximum query coverage (-maxcov). Swapping values.")
                        param['mincov'] = max_cov
                        param['maxcov'] = min_cov
                    else:
                        param['mincov'] = min_cov
                        param['maxcov'] = max_cov

    if valid == True:
        if args.minlen == None and args.maxlen == None:
            param['minlen'] = 4000
            param['maxlen'] = 50000
        else:
            try:
                min_len = int(args.minlen) if args.minlen else 0
                max_len = int(args.maxlen) if args.maxlen else float('inf')
            except Exception:
                traceback.print_exc(file=sys.stdout)
                print("ERROR: Minimum or maximum element length (-minlen, -maxlen) is not an integer!\n")
                valid = False
            else:
                if min_len < 0:
                    print("ERROR: Minimum element length (-minlen) must be greater than or equal to 0!")
                    valid = False
                if max_len < 0:
                    print("ERROR: Maximum element length (-maxlen) must be greater than or equal to 0!")
                    valid = False
                if valid:
                    if min_len > max_len:
                        print("Warning: Minimum element length (-minlen) is greater than maximum element length (-maxlen). Swapping values.")
                        param['minlen'] = max_len
                        param['maxlen'] = min_len
                    else:
                        param['minlen'] = min_len
                        param['maxlen'] = max_len

    if valid == True:
        if args.color == None:
            param['color'] = [255, 0, 0]
        else:
            try:
                color_parts = args.color.split(",")
                if len(color_parts) == 3:
                    colors = []
                    for num in color_parts:
                        n = int(num)
                        if 0 <= n <= 255:
                            colors.append(n)
                        else:
                            print("ERROR: The RGB color of the element (-color) must be three integers between 0 and 255!")
                            valid = False
                            break
                    if valid:
                        param['color'] = colors
                else:
                    print("ERROR: The RGB color of the element (-color) must be three integers between 0 and 255 separated by commas!")
                    valid = False
            except Exception:
                traceback.print_exc(file=sys.stdout)
                print("ERROR: The RGB color of the element (-color) must be three integers between 0 and 255!\n")
                valid = False

    if valid == True:
        try:
            os.chdir(call)
            if args.out is None:
                param["out"] = os.path.join(call, 'output_dir')
            else:
                param["out"] = os.path.join(call, args.out) if os.path.basename(args.out) == args.out else os.path.abspath(args.out)

            original_path = param["out"]
            counter = 1
            while os.path.exists(param["out"]):
                param["out"] = original_path + "_{}".format(counter)
                counter += 1

            os.makedirs(param["out"])
            os.chdir(param["out"])
        except Exception:
            traceback.print_exc(file=sys.stdout)
            print("ERROR: Output directory (-out) not valid!\n")
            valid = False
        else:
            print("Creating output directory (-out) {}...".format(param["out"]))

    return valid, param


def blast(param, q, blast_time='', seqs=[]):
    tab = os.path.join(param["out"], 'blastn2.tab')
    blast_start = datetime.now()

    blastn_cmd_args = {
        'query': q,
        'outfmt': "'7 qseqid sseqid qcovs qlen slen qstart qend evalue'",
        'out': tab
    }

    if param['run'] == 'local':
        print("Starting local BLAST search...")
        blastn_cmd_args['db'] = param['d']
        blastn_cmd_args['num_threads'] = param['cpu']
    elif param['run'] == 'web':
        print("Starting web BLAST search...")
        blastn_cmd_args['db'] = "nt"
        blastn_cmd_args['remote'] = True
        blastn_cmd_args['task'] = 'megablast'
        if 'org' in param:
            blastn_cmd_args['entrez_query'] = "'{}'".format(param['org'])

    comando_blastn = NcbiblastnCommandline(**blastn_cmd_args)
    stdout, stderr = comando_blastn()

    if blast_time == '':
        blast_time = (datetime.now() - blast_start)
    else:
        blast_time += (datetime.now() - blast_start)

    if seqs == []:
        seqs = missingquery(tab=tab, qry=q)
    else:
        seqs = missingquery(tab=tab, qry=q, seqs=seqs)

    if not seqs == False:
        blast(param, os.path.join(param["out"], 'newquery.fasta'), blast_time, seqs)
    else:
        if os.path.isfile(os.path.join(param["out"], 'newquery.fasta')):
            os.remove(os.path.join(param["out"], 'newquery.fasta'))

        if os.path.isfile(os.path.join(param["out"], 'blastn.tab')):
            os.remove(os.path.join(param["out"], 'blastn.tab'))
        os.rename(os.path.join(param["out"], 'blastn2.tab'), os.path.join(param["out"], 'blastn.tab'))
    return blast_time


def missingquery(tab, qry, seqs=[]):
    qfile_content = open(qry, "r").read()
    tabfile_content = open(tab, "r").read()

    tabqs_raw = tabfile_content.split("# BLAST")
    tabqs = [t for t in tabqs_raw if t.strip()]

    hits_found = {}
    misseqs = []
    newq_content = ''

    query_ids_in_fasta = []
    for record in SeqIO.parse(qry, "fasta"):
        query_ids_in_fasta.append(record.id)

    for query_id in query_ids_in_fasta:
        found_in_tab = False
        for block in tabqs:
            if "# Query: {}".format(query_id) in block:
                found_in_tab = True
                break

        if not found_in_tab:
            for record in SeqIO.parse(io.StringIO(qfile_content), "fasta"):
                if record.id == query_id:
                    misseqs.append(query_id)
                    newq_content += ">{}\n{}\n".format(record.id, str(record.seq))
                    break

    if newq_content == '':
        for line in tabfile_content.split('\n'):
            if not re.search('#', line) and line.strip():
                h = line.split('\t')
                if len(h) > 0:
                    query_id = h[0]
                    if query_id not in hits_found:
                        hits_found[query_id] = []
                    hits_found[query_id].append(h)

    if len(misseqs) >= 1 and len(seqs) < len(misseqs):
        seqs = misseqs
        with open(os.path.join(param["out"], 'newquery.fasta'), "w") as file:
            file.write(newq_content)
        return seqs
    elif len(misseqs) >= 1 and len(seqs) == len(misseqs):
        return False
    else:
        return False


def extract(seqs_content, qid_to_find, start, end):
    for record in SeqIO.parse(io.StringIO(seqs_content), "fasta"):
        if record.id == qid_to_find:
            return str(record.seq)[start - 1:end]
    return ""


try:
    if args.version:
        print("insertion_finder v{}".format(version))
        sys.exit(0)

    if args.help:
        print(help)
        sys.exit(0)

    valid, param = validate_args(args)

    if valid:
        if 'tab' not in param:
            blast_time = blast(param, param['q'])
            param['tab'] = os.path.join(param["out"], 'blastn.tab')
            qid, columns, hits = blast_parse(param['tab'])
            param['columns'] = columns
            param['qid'] = qid
            param['hits'] = hits
        else:
            blast_time = 'Not performed (table provided)'

        log = open(os.path.join(param["out"], 'log.txt'), "w")
        tabular = open(os.path.join(param["out"], 'elements_table.tab'), "w")

        log.write('insertion_finder v{}\n'.format(version))
        log.write('(c) 2021. Arthur Gruber & Giuliana Pola\n')
        log.write('BLASTn search performed on: {}\n'.format(datetime.now()))
        log.write('Parameters:\n')
        for key in param:
            log.write('{}: {}\n'.format(key, param[key]))
        log.write('\n')
        tabular.write('Query ID\tSubject ID\tIs element valid?\tElement start\tElement end\tElement length\tIs insertion found?\n')
        log.write("Extracting possible insertion elements from BLASTn table: {}\n".format(param['tab']))
        log.flush()

        file = open(param['tab'], "r")
        tab_content = file.read()
        file.close()

        qfile_content = open(param['q'], "r").read()

        econt = 0
        cont = 0

        blast_reports = re.split(r'# BLAST', tab_content)
        blast_reports = [report.strip() for report in blast_reports if report.strip()]

        for qid_from_param in param['qid']:
            current_query_report = None
            for report in blast_reports:
                if "# Query: {}".format(qid_from_param) in report:
                    current_query_report = report
                    break

            if not current_query_report:
                log.write("[{}] ERROR: No BLAST hits found for this query in the table.\n".format(qid_from_param))
                log.flush()
                continue

            report_lines = current_query_report.split('\n')
            report_header_lines = []
            report_hit_lines = []

            for line in report_lines:
                if line.startswith('#'):
                    report_header_lines.append(line)
                elif line.strip():
                    report_hit_lines.append(line)

            current_columns = param['columns']

            query_len_from_fasta = None
            for record in SeqIO.parse(io.StringIO(qfile_content), "fasta"):
                if record.id == qid_from_param:
                    query_len_from_fasta = len(record.seq)
                    break

            if query_len_from_fasta is None:
                log.write("[{}] ERROR: Query sequence not found in query file!\n".format(qid_from_param))
                log.flush()
                continue


            for hit_line in report_hit_lines:
                h = hit_line.split('\t')

                if h[current_columns.index('query id')] != qid_from_param:
                    continue

                try:
                    qstart = int(h[current_columns.index('q. start')])
                    qend = int(h[current_columns.index('q. end')])
                    qlen = int(h[current_columns.index('query length')])
                    qcovs = float(h[current_columns.index('% query coverage per subject')])
                    sid = h[current_columns.index('subject id')]
                    slen = int(h[current_columns.index('subject length')])
                    evalue = float(h[current_columns.index('evalue')])
                except (ValueError, IndexError) as e:
                    log.write("[{}] ERROR: Malformed BLAST hit line or missing column: {} - {}\n".format(qid_from_param, hit_line, e))
                    log.flush()
                    continue

                is_valid_element = 'no'
                is_insertion_found = 'no'
                elen = abs(qend - qstart) + 1

                if (qstart <= param['enddist'] or query_len_from_fasta - qend <= param['enddist']) and \
                   (param['mincov'] <= qcovs <= param['maxcov']) and \
                   (param['minlen'] <= elen <= param['maxlen']):

                    is_valid_element = 'yes'
                    is_insertion_found = 'yes'
                    econt += 1

                    fasta_filename = os.path.join(param["out"], 'element_{}.fasta'.format(qid_from_param))
                    file_idx = 0
                    temp_fasta_filename = fasta_filename
                    while os.path.exists(temp_fasta_filename):
                        file_idx += 1
                        temp_fasta_filename = os.path.join(param["out"], 'element_{}_{}.fasta'.format(qid_from_param, file_idx))

                    try:
                        with open(temp_fasta_filename, 'w') as fasta_file:
                            fasta_file.write(">{0} - element - {1}-{2}\n".format(qid_from_param, qstart, qend))
                            fasta_file.write(extract(qfile_content, qid_from_param, qstart, qend) + "\n")
                    except Exception as e:
                        traceback.print_exc(file=sys.stdout)
                        log.write("\n[{}] ERROR: Element's fasta was not written for {}! Error: {}\n".format(qid_from_param, temp_fasta_filename, e))
                        log.flush()
                else:
                    cont += 1
                    if not (qstart <= param['enddist'] or query_len_from_fasta - qend <= param['enddist']):
                        log.write("[{}] ERROR: Invalid element, distance to query tip ({}) is greater than allowed ({})!\n".format(qid_from_param, min(qstart, query_len_from_fasta - qend), param['enddist']))
                        log.flush()
                    if not (param['mincov'] <= qcovs <= param['maxcov']):
                        log.write("[{}] ERROR: Invalid element, query coverage ({}%) is not within {}-{}%!\n".format(qid_from_param, qcovs, param['mincov'], param['maxcov']))
                        log.flush()
                    if not (param['minlen'] <= elen <= param['maxlen']):
                        log.write("[{}] ERROR: Invalid element, length ({} bp) is not within {}-{} bp!\n".format(qid_from_param, elen, param['minlen'], param['maxlen']))
                        log.flush()

                tabular.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(qid_from_param, sid, is_valid_element, qstart, qend, elen, is_insertion_found))

        tabular.write('\n\nProcessed {} queries with {} valid elements\n'.format(len(param['qid']), econt))
        log.write('\n\nProcessed {} queries with {} valid elements\n'.format(len(param['qid']), econt))
        log.flush()
        print('Processed {} queries with {} valid elements'.format(len(param['qid']), econt))

        program_execution_time = datetime.now() - start_time
        log.write('\nProgram execution time: {}\n'.format(program_execution_time))
        log.flush()
        tabular.write('\nProgram execution time: {}\n'.format(program_execution_time))

        log.close()
        tabular.close()
        print('Program execution time: {}'.format(program_execution_time))
        print('\a')
    else:
        print("ERROR: Validation failed. Please check the error messages above.")
        sys.exit(1)

except Exception as e:
    traceback.print_exc(file=sys.stdout)
    print("ERROR: An unhandled error occurred: {}".format(e))
    sys.exit(1)