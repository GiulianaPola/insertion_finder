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
ajuda = ajuda + 'Usage\tinsertion_finder.py -q <query file> -d <database file>\n'
ajuda = ajuda + '\nMandatory parameters:\n'
ajuda = ajuda + '-q <file>\tSequence to search with\n'
ajuda = ajuda + '-d <file>\tDatabase to BLAST against\n'
ajuda = ajuda + '\nOptional parameters:\n'
ajuda = ajuda + "-o <path>\tOutput directory\n"
ajuda = ajuda + "-min <int>\tMinimum element's length in base pairs(bp) (default: >=5000)\n"
ajuda = ajuda + "-max <int>\tMaximum element's length in base pairs(bp) (default: <=150000)\n"
ajuda = ajuda + "-c <int>\tElement RGB color that is shown by the feature table (default: 255,0,0)"

parser = argparse.ArgumentParser(add_help=False)
parser.add_argument('-q')
parser.add_argument('-d')
parser.add_argument('-o')
parser.add_argument('-c')
parser.add_argument('-min')
parser.add_argument('-max')
parser.add_argument('-h', '--help', action='store_true')
args = parser.parse_args()

def isfasta(filename):
    with open(filename, "r") as handle:
        fasta = SeqIO.parse(handle, "fasta")
        return any(fasta)
  
def extract(qseqs, qid, estart, eend):
  lin=qseqs[0:qseqs.find(qid)].count("\n")+1
  seq=qseqs.split("\n")[lin]
  return seq[int(estart-1):eend]

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
  return m
  
def joinlists(list1,list2):
  joined=[]
  for i in range(len(list1)):
    joined.append([list1[i],list2[i]])
  return joined

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

def validateargs(args):
  valid=True
  param=dict()

  if not args.o==None:
    if os.path.isdir(args.o):
      param['o']=args.o
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
  
  if args.q==None:
    print("Missing query file!\n")
    valid=False
  elif not args.q==None:
   if not os.path.isfile(args.q):
     print("Query file not exist!\n")
     valid=False
   else:
     if isfasta(args.q):
       param['q']=args.q
     else:
       print("Invalid query file, invalid fasta formatting!\n")
       valid=False
  
  if args.d==None:
    print("Missing database file!\n")
    valid=False
  elif not args.d==None:
   if not os.path.isfile(args.d):
     print("Database file not exist!\n")
     valid=False
   else:
     if isfasta(args.d):
       param['d']=args.d
     else:
       print("Invalid query file, invalid fasta formatting!\n")
       valid=False
  
  if args.min==None:
    param['min']=5000
  else:
    try:
      param['min']=int(args.min)
    except:
      print("Minimum element's length is not integer\n")
      valid=False
  
  if args.max==None:
    param['max']=150000
  else:
    try:
      param['max']=int(args.max)
    except:
      print("Maximum element's length is not integer\n")
      valid=False
  
  if args.c==None:
    param['c']=[255,0,0]
  else:
      try:
        if args.c.count(",")==2:
          param['c']=[]
          for num in args.c.split(","):
            if int(num)>=0 & int(num)<=255:
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
    
  return valid,param

if not len(sys.argv)>1:
    print("Missing all arguments!\n")
    print(ajuda)
elif args.help == True:
    print(ajuda)
elif args.help == False:
  valid,param=validateargs(args)
  if valid==False:
    print(ajuda)
    
  else:
    qseqs=open(param['q'], "r").read()
    
    try:
      comando_blastn = NcbiblastnCommandline( \
      query=param['q'], db=param['d'], \
      outfmt="'7 qseqid sseqid qcovs sstrand score qlen slen qstart qend sstart send pident'", out=os.path.join(param['o'], "blastn.tab"),num_threads=10)
      # Executando
      stdout, stderr = comando_blastn()
     
      # Abrindo resultado
      blast_result = open(os.path.join(param['o'], "blastn.tab"),"r")
      linhas=[]
      colunas=[]
      cont=0
      econt=0
      qid=''
      element=''
      print("Opening BLASTn table...")   
    except:
      print("The BLASTn search was not completed successfully!")
    else:
      tabular=open(os.path.join(param['o'], 'elements.txt'),'w')
      tabular.write("Query file: {}\n".format(param['q']))
      log.append(("Query file: {}".format(param['q'])))
      tabular.write("Database file: {}\n".format(param['d']))
      log.append("Database file: {}".format(param['d']))
      tabular.write("Element valid length: {0}-{1}\n".format(param['min'],param['max']))
      log.append("Element valid length: {0}-{1}".format(param['min'],param['max']))
      tabular.write("query ID\tsubject ID\telement identification\telement 5' coordinate\telement 3' coordinate\telement length\tvalid\n")
      
      for linha in blast_result.readlines():
        linha=linha.replace('\n','')
        if 'Query:' in linha:
          qid=linha.split(': ')[1]
          log.append("\nQuery id\t{}".format(qid))
          cont+=1
        if 'Fields:' in linha:
          colunas=linha.split(': ')[1]
          colunas=colunas.split(', ')
        if 'hits' in linha:
          if int(linha.split(' ')[1])==0:
            element='no'
            sid='no hits'
            tabular.write("{0}\t{1}\t{2}\t\t\t\t{3}\n".format(qid,sid,element,'no'))
            log.append('No hits!')
          else:
            log.append('{} hits!'.format(linha.split(' ')[1]))
        if not '#' in linha:
          linhas.append(linha.split('\t'))
        if '# BLAST' in linha and not qid=='':
          estart=0
          eend=0
          elen=0
          if not linhas==[]:
            df=pd.DataFrame(columns=colunas,data=linhas)
            df=df.apply(pd.to_numeric, errors='ignore')
            df1=df[['query id','query length','subject id','subject length','% query coverage per subject','score','% identity']].drop_duplicates()
            qlen=df1['query length'].tolist()[0]
            log.append("Query length\t{}".format(qlen))
            filtro='score'
            df1=df1.groupby('subject id').agg({'subject length':'mean',filtro:'sum','% query coverage per subject':'mean','% identity':'mean'}).reset_index()
            df1=df1.sort_values(by=[filtro],ascending=False)
            sid=df1['subject id'].tolist()[0]
            log.extend(str(df1.iloc[1]).split('\n')[:-1])
            df1=df.loc[(df['subject id'] == sid)&(df['query id'] == qid)]
            if len(df1)==1:
              log.append("Query has no elements!")
              element='no'
              tabular.write("{0}\t{1}\t{2}\t\t\t\t{3}\n".format(qid,sid,element,'no'))
            else:
              reads=sorted(joinlists(df1['q. start'].tolist(),df1['q. end'].tolist()))
              contigs=assembly(reads)
              log.append("Alignments\t{}".format(str(reads)))
              if len(contigs)>1:
                element='yes'
                start,end=splitlist(contigs)
                estart,eend,elen=findelement(start,end)
                log.append("Query has element!")
                log.append("Element coodinates\t{} - {}".format(estart,eend))
                log.append("Element length\t{}".format(elen))
                if elen>=param['min'] and elen<=param['max']:
                  log.append("Valid element!")
                  os.mkdir(os.path.join(param['o'], qid))
                  tabular.write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\n".format(qid,sid,element,estart,eend,elen,'yes'))
                  try:
                    log.append("Writing element's feature table...")
                    ft=open(os.path.join(param['o'], qid, str(qid+'_element.gb')),'w')
                    ft.write("     misc_feature    {0}..{1}\n".format(eend,estart))
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
                  econt+=1
                else:
                  log.append("Invalid element!")
                  tabular.write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\n".format(qid,sid,element,estart,eend,elen,'no'))
              else:
                log.append("Query has no elements!")
                element='no'
                tabular.write("{0}\t{1}\t{2}\t\t\t\t{3}\n".format(qid,sid,element,'no'))
          linhas=[]
          colunas=[]
          qid=''
          element=''
      log.append('\nProcessed {} queries with {} valid elements'.format(cont,econt))
      tabular.close()
      erros=open(os.path.join(param['o'], 'file.log'),'w')
      erros.write('\n'.join(str(v) for v in log))
