#!/usr/bin/python
import traceback
import warnings
import sys

def warn_with_traceback(message, category, filename, lineno, file=None, line=None):

    log = file if hasattr(file,'write') else sys.stderr
    traceback.print_stack(file=log)
    log.write(warnings.formatwarning(message, category, filename, lineno, line))

warnings.showwarning = warn_with_traceback

#python insertion_finder.py -query "NZ_EQ973357.fasta" -db "db/genomas.fasta"
#python insertion_finder.py -query "PVBG01000001_1.fasta" -db "patric/Rhodo_Rhizo.fasta"
#python insertion_finder.py -q "elements.fasta" -d "patric/Rhodo_Rhizo.fasta" -o "teste_Rhodobacterales"

import argparse
import pandas as pd
from Bio import SeqIO
from Bio.Blast.Applications import NcbiblastnCommandline
import sys
import os
log=[]

ajuda = 'insertion_finder - element insertion finder in a genome through a BLAST search\n'
ajuda = ajuda + '(c) 2021. Arthur Gruber & Giuliana Pola\n'
ajuda = ajuda + 'Usage: insertion_finder.py -q <query file> -d <database file>\n'
ajuda = ajuda + '\tinsertion_finder.py -q <query file> -tab <BLASTn table file>\n'
ajuda = ajuda + '\nMandatory parameters:\n'
ajuda = ajuda + '-q <fasta or multifasta file>\tSequence to search with\n'
ajuda = ajuda + '-d <multifasta file>\tDatabase to BLAST against\n'
ajuda = ajuda + '-tab <table file>\tBLASTn search result table (fields: qseqid,sseqid,qcovs,qlen,slen,qstart,qend)\n'
ajuda = ajuda + '\nOptional parameters:\n'
ajuda = ajuda + "-o <path>\tOutput directory (default: output_dir1)\n"
ajuda = ajuda + "-enddist <int>\tMaximum distance between block tip and query tip in base pairs(bp) (default: 1000)\n"
ajuda = ajuda + "-minlen <int>\tMinimum element's length in base pairs(bp) (default: 5000)\n"
ajuda = ajuda + "-maxlen <int>\tMaximum element's length in base pairs(bp) (default: 50000)\n"
ajuda = ajuda + "-mincov <int>\tMinimum % query coverage per subject (default: 30)\n"
ajuda = ajuda + "-maxcov <int>\tMaximum % query coverage per subject (default: 90)\n"
ajuda = ajuda + "-c <int>\tElement RGB color that is shown by the feature table, three integers between 0 and 255 separated by commas (default: 255,0,0)"

parser = argparse.ArgumentParser(add_help=False)
parser.add_argument('-q')
parser.add_argument('-d')
parser.add_argument('-o')
parser.add_argument('-c')
parser.add_argument('-tab')
parser.add_argument('-enddist')
parser.add_argument('-minlen')
parser.add_argument('-maxlen')
parser.add_argument('-mincov')
parser.add_argument('-maxcov')
parser.add_argument('-h', '--help', action='store_true')
args = parser.parse_args()

def isfasta(filename):
    with open(filename, "r") as handle:
        fasta = SeqIO.parse(handle, "fasta")
        return any(fasta)

def validateargs(args):
  valid=True
  param=dict()
  
  if args.q==None:
    print("Missing query file!")
    valid=False
  elif not args.q==None:
   if not os.path.isfile(args.q):
     print("Query file not exist!")
     valid=False
   else:
     if isfasta(args.q):
       param['q']=args.q
     else:
       print("Invalid query file, invalid fasta formatting!")
       valid=False
  
  if args.tab==None and args.d==None:
    print("Missing BLASTn table or database file!")
    valid=False
  elif not args.tab==None:
    if not os.path.isfile(args.tab):
      print("BLASTn table file not exist!")
      valid=False
    else:
      qid,colunas,hits,d=opentable(args.tab)
      if qid==[] or colunas==[]:
        valid=False
      else:
        for col in ['query id', 'subject id', '% query coverage per subject', 'query length', 'subject length', 'q. start', 'q. end']:
          if not col in colunas:
            valid=False
            print("Column {} is not in table!".format(col))
      if valid==True:
        param['qid']=qid
        param['colunas']=colunas
        param['hits']=hits 
        param['d']=d
  elif not args.d==None:
    if not os.path.isfile(args.d):
      print("Database file not exist!")
      valid=False
    else:
      if isfasta(args.d):
        param['d']=args.d
      else:
        print("Invalid database file, invalid fasta formatting!")
        valid=False
  
  if args.enddist==None:
    param['enddist']=1000
  else:
    try:
      param['enddist']=int(args.enddist)
    except:
      print("Maximum distance is not integer!")
      valid=False
  
  if args.minlen==None:
    param['minlen']=4000
  else:
    try:
      param['minlen']=int(args.minlen)
    except:
      print("Minimum element's length is not integer!")
      valid=False
  
  if args.mincov==None:
    param['mincov']=30
  else:
    try:
      param['mincov']=int(args.mincov)
    except:
      print("Minimum query coverage is not integer!")
      valid=False
  
  if args.maxlen==None:
    param['maxlen']=50000
  else:
    try:
      param['maxlen']=int(args.maxlen)
    except:
      print("Maximum element's length is not integer!")
      valid=False

  if args.maxcov==None:
    param['maxcov']=90
  else:
    try:
      param['maxcov']=int(args.maxcov)
    except:
      print("Maximum query coverage is not integer!")
      valid=False
  
  if args.c==None:
    param['c']=[255,0,0]
  else:
      try:
        if args.c.count(",")==2:
          param['c']=[]
          for num in args.c.split(","):
            n=int(num)
            if n>=0 and n<=255:
              param['c'].append(int(num))
            else:
              print("Invalid element RGB color!\n")
              valid=False
        else:
          print("Invalid element RGB color!\n")
          valid=False
      except:
        print("Invalid element RGB color!\n")
        valid=False
    
  if valid==True:
   if not args.o==None:
     if os.path.isdir(args.o):
       print("Output directory {} exists!".format(args.o))
       i=1
       out=args.o
       while os.path.isdir(out):
         i+=1
         out=args.o+"_"+str(i)
       os.mkdir(out)
       print("Creating output directory {}...".format(out))
       param['o']=out
     else:
       try:
         os.mkdir(args.o)
         print("Creating output directory {}...".format(args.o))
         param['o']=args.o
       except:
         print("Output directory not valid\n")
         valid=False
   else:
     i=1
     out='output_dir1'
     while os.path.isdir(out):
       i+=1
       out='output_dir'+str(i)
     os.mkdir(out)
     print("Creating output directory {}...".format(out))
     param['o']=out
    
  return valid,param

def opentable(tab):
  blast_result = open(tab,"r")
  qid=[]
  hits=[]
  colunas=[]
  d=''
  for linha in set(blast_result.readlines()):
    linha=linha.replace('\n','')
    if 'Query:' in linha:
      qid.append(linha.split(': ')[1])
    if 'Fields:' in linha:
      colunas=linha.split(': ')[1].split(', ')
    if 'Database:' in linha:
      d=linha.split(': ')[1]
    if not '#' in linha:
      hits.append(linha.split('\t'))
  return qid,colunas,hits,d

def joinlists(list1,list2):
  joined=[]
  for i in range(len(list1)):
    joined.append([list1[i],list2[i]])
  return joined

def assembly(reads):
  reads.sort()
  m = []
  s = -10000
  max = -100000
  for i in range(len(reads)):
      a = reads[i]
      if a[0] > max:
          if i != 0:
              m.append([s,max])
          max = a[1]
          s = a[0]
      else:
          if a[1] >= max:
              max = a[1]
  if max != -100000 and [s, max] not in m:
      m.append([s, max])
  return sorted(m)

def oneblock(contigs,qlen,enddist):
    start=contigs[0][0]-1
    end=qlen-contigs[0][1]
    if start<=enddist and end>=enddist:
      estart=contigs[0][1]+1
      eend=qlen
      elen=eend-estart
      return estart,eend,elen
    elif start>=enddist and end<=enddist:
      estart=1
      eend=contigs[0][0]-1
      elen=eend-estart
      return estart,eend,elen
    else:
      return 0,0,0

def splitlist(joined):
  list1=[]
  list2=[]
  for i in joined:
    list1.append(i[0])
    list2.append(i[1])
  return list1,list2

def findelement(start,end):
  maxdist=0
  for k in range(len(end)-1): 
    distance=start[k+1]-end[k]
    if distance>maxdist:
      #print("maxdist=",distance)
      maxdist=distance
      #print("estart=",end[k]+1)
      estart=end[k]+1
      #print("eend=",start[k+1]-1)
      eend=start[k+1]-1
  elen=eend-estart
  return estart,eend,elen

def extract(qseqs, qid, estart, eend):
  lin=qseqs[0:qseqs.find(qid)].count("\n")+1
  seq=qseqs.split("\n")[lin]
  return seq[int(estart-1):eend]

if not len(sys.argv)>1:
    print(ajuda)
elif args.help == True:
    print(ajuda)
elif args.help == False:
  valid,param=validateargs(args)
  #print(param)
  if valid==False:
    print(ajuda)
  else:
    qseqs=open(param['q'], "r").read()
    try:
      if args.tab ==None:
        comando_blastn = NcbiblastnCommandline( \
        query=param['q'], db=param['d'], \
        outfmt="'7 qseqid sseqid qcovs qlen slen qstart qend'", out=os.path.join(param['o'], "blastn.tab"),num_threads=10)
        stdout, stderr = comando_blastn()
        tab=os.path.join(param['o'], "blastn.tab")
        qid,colunas,hits,d=opentable(tab)
        param['qid']=qid
        param['colunas']=colunas
        param['hits']=hits 
        print("Opening BLASTn table...")   
    except:
      print("The BLASTn search was not completed successfully!")
    else:
      cont=0
      econt=0
      tabular=open(os.path.join(param['o'], 'elements.txt'),'w')
      tabular.write("Query file: {}\n".format(param['q']))
      log.append(("Query file: {}".format(param['q'])))
      tabular.write("Database file: {}\n".format(param['d']))
      log.append("Database file: {}".format(param['d']))
      tabular.write("Element valid length: {0}-{1}\n".format(param['minlen'],param['maxlen']))
      log.append("Element valid length: {0}-{1}".format(param['minlen'],param['maxlen']))
      tabular.write("Valid % query coverage per subject: {0}-{1}\n".format(param['mincov'],param['maxcov']))
      log.append("Valid % query coverage per subject: {0}-{1}".format(param['mincov'],param['maxcov']))
      tabular.write("query ID\tsubject ID\telement identification\telement 5' coordinate\telement 3' coordinate\telement length\tvalid\n")
      df=pd.DataFrame(columns=param['colunas'],data=param['hits'])
      df=df.apply(pd.to_numeric, errors='ignore').drop_duplicates()
      for qid in param['qid']:
        df1=df.loc[df['query id'] == qid]
        log.append("\nQuery id\t{}".format(qid))
        log.append("{} hits!".format(len(df1)))
        if len(df1)==0:
          element='no'
          sid='no hits'
          tabular.write("{0}\t{1}\t{2}\t\t\t\t{3}\n".format(qid,sid,element,'no'))
        if len(df1)>0:
          qlen=df1['query length'].tolist()[0]
          log.append("Query length\t{}".format(qlen))
          df1=df1.groupby('subject id').agg({'subject length':'mean','% query coverage per subject':'mean'}).reset_index()
          df1=df1.loc[(df1['% query coverage per subject'] >= param['mincov']) & (df1['% query coverage per subject'] <= param['maxcov'])]
          df1=df1.sort_values(by=['% query coverage per subject'],ascending=False)
          if len(df1)==0:
            element='no'
            sid='no hits'
            tabular.write("{0}\t{1}\t{2}\t\t\t\t{3}\n".format(qid,sid,element,'no'))
            log.append('No hits with % query coverage per subject between {} and {}!'.format(param['mincov'],param['maxcov']))
          else:
            #df1.to_csv(os.path.join(param['o'],'{}_blastn.tab'.format(qid)),sep="\t",index=False)
            sid=df1['subject id'].tolist()[0]
            log.extend(str(df1.iloc[1]).split('\n')[:-1])
            df1=df.loc[(df['subject id'] == sid)&(df['query id'] == qid)]
            reads=sorted(joinlists(df1['q. start'].tolist(),df1['q. end'].tolist()))
            contigs=assembly(reads)
            log.append("Alignments\t{}".format(str(reads)))
            if len(contigs)==1:
              log.append("One block alignment!")
              estart,eend,elen=oneblock(contigs,qlen,param['enddist'])
              #print(estart,eend,elen)
              if estart>0:
                element='yes'
                log.append("Query has element!")
                log.append("Element coodinates\t{} - {}".format(estart,eend))
                log.append("Element length\t{}".format(elen))
                if elen>=param['minlen'] and elen<=param['maxlen']:
                  log.append("Valid element!")
                  econt+=1
                  os.mkdir(os.path.join(param['o'], qid))
                  tabular.write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\n".format(qid,sid,element,estart,eend,elen,'yes'))
                  try:
                    log.append("Writing element's feature table...")
                    ft=open(os.path.join(param['o'], qid, str(qid+'_element.gb')),'w')
                    ft.write("     misc_feature    {0}..{1}\n".format(estart,eend))
                    ft.write("                     /label=element\n")
                    ft.write("                     /color={} {} {}\n".format(param['c'][0],param['c'][1],param['c'][2]))
                    ft.close()
                  except:
                    log.append("Element's feature table was't writen!")
                  try:
                    log.append("Writing element's fasta...")
                    fasta=open(os.path.join(param['o'], qid, str(qid+'_element.fasta')),'w')
                    fasta.write(">{0} - element - {1}-{2}\n".format(qid,estart,eend))
                    fasta.write(extract(qseqs, qid, estart, eend))
                    fasta.close()
                  except:
                    log.append("Element's fasta was't writen!")
                else:
                  log.append("Invalid element!")
                  tabular.write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\n".format(qid,sid,element,estart,eend,elen,'no'))
              else:
                log.append("No element!")
                element='no'
                tabular.write("{0}\t{1}\t{2}\t\t\t\t{3}\n".format(qid,sid,element,'no'))    
            elif len(contigs)>1:
              element='yes'
              start,end=splitlist(contigs)
              estart,eend,elen=findelement(start,end)
              log.append("Query has element!")
              log.append("Element coodinates\t{} - {}".format(estart,eend))
              log.append("Element length\t{}".format(elen))
              if elen>=param['minlen'] and elen<=param['maxlen']:
                log.append("Valid element!")
                econt+=1
                os.mkdir(os.path.join(param['o'], qid))
                tabular.write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\n".format(qid,sid,element,estart,eend,elen,'yes'))
                try:
                  log.append("Writing element's feature table...")
                  ft=open(os.path.join(param['o'], qid, str(qid+'_element.gb')),'w')
                  ft.write("     misc_feature    {0}..{1}\n".format(estart,eend))
                  ft.write("                     /label=element\n")
                  ft.write("                     /color={} {} {}\n".format(param['c'][0],param['c'][1],param['c'][2]))
                  ft.close()
                except:
                  log.append("Element's feature table was't writen!")
                try:
                  log.append("Writing element's fasta...")
                  fasta=open(os.path.join(param['o'], qid, str(qid+'_element.fasta')),'w')
                  fasta.write(">{0} - element - {1}-{2}\n".format(qid,estart,eend))
                  fasta.write(extract(qseqs, qid, estart, eend))
                  fasta.close()
                except:
                  log.append("Element's fasta was't writen!")
              else:
                log.append("Invalid element!")
                tabular.write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\n".format(qid,sid,element,estart,eend,elen,'no'))
      log.append('\nProcessed {} queries with {} valid elements'.format(len(param['qid']),econt))
      tabular.close()
      erros=open(os.path.join(param['o'], 'file.log'),'w')
      erros.write('\n'.join(str(v) for v in log))
