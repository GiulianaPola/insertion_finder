#!/usr/bin/python
#import traceback
#import warnings
#import sys

#def warn_with_traceback(message, category, filename, lineno, file=None, line=None):

    #log = file if hasattr(file,'write') else sys.stderr
    #traceback.print_stack(file=log)
    #log.write(warnings.formatwarning(message, category, filename, lineno, line))

#warnings.showwarning = warn_with_traceback

#python insertion_finder.py -query "NZ_EQ973357.fasta" -db "db/genomas.fasta"
#python insertion_finder.py -query "PVBG01000001_1.fasta" -db "patric/Rhodo_Rhizo.fasta"
#python insertion_finder.py -query "elements.fasta" -db "patric/Rhodo_Rhizo.fasta"

import argparse
import pandas as pd
from Bio.Blast.Applications import NcbiblastnCommandline
import sys
import os

ajuda = 'InsertionFinder - element insertion locator in a genome through a BLAST search\n'
ajuda = ajuda + '(c) 2021. Arthur Gruber & Giuliana Pola\n'
ajuda = ajuda + 'Usage: insertion_finder.py -query <sequence to search with> -db <database to BLAST against> -out <output file> \n'
ajuda = ajuda + '\nMandatory parameters:\n'
ajuda = ajuda + '-q <file>\tSequence to search with\n'
ajuda = ajuda + '-d <file>\tDatabase to BLAST against\n'
ajuda = ajuda + '\nOptional parameters:\n'
ajuda = ajuda + "-o <path>\tOutput directory\n"
ajuda = ajuda + "-min <int>\tMinimum element's length in base pairs(bp) (default: >=5000)\n"
ajuda = ajuda + "-max <int>\tMaximum element's length in base pairs(bp) (default: <=150.000)\n"

parser = argparse.ArgumentParser(add_help=False)
parser.add_argument('-q')
parser.add_argument('-d')
parser.add_argument('-o')
parser.add_argument('-min')
parser.add_argument('-max')
parser.add_argument('-h', '--help', action='store_true')
args = parser.parse_args()

def splitfile(query):
  qfastas=dict()
  arquivo=open(query,"r")
  linhas=arquivo.readlines()
  for i in range(len(linhas)):
    if ">" in linhas[i]:
      head=linhas[i][1:]
      if " " in head:
        qid=head.split(" ")[0]
      else:
        qid=head[0:-1]
      qfastas[qid]=linhas[i+1]
  return qfastas
  
def extract(seq, estart, eend):
  return seq[int(estart-1):eend]

def assembly(reads):
  reads.sort(key = lambda x: x[0])
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
  print("Calculating element's coordinates...")
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

if not len(sys.argv)>1:
    print(ajuda)
elif args.help == True:
    print(ajuda)
elif not args.q==None and not args.d==None:
  if not os.path.isfile(args.q):
    print ("Query file not exist")
  elif not os.path.isfile(args.d):
    print ("Database file not exist")
  else:  
    print ("Database and query file exist")
    if not args.o==None:
     if not os.path.isdir(args.o):
       print ("Output directory not exist")
       print ("Creating directory...")
       os.mkdir(args.o)
     else:
       print ("Output directory exist")
  
    qfastas=splitfile(args.q)
    #print(sorted(qfastas.keys()))
    
    if not args.o==None:
      comando_blastn = NcbiblastnCommandline( \
      query=args.q, db=args.d, \
      outfmt="'7 qseqid sseqid qcovs sstrand score qlen slen qstart qend sstart send pident'", out=os.path.join(args.o, "blastn.tab"),num_threads=10)
    else:
      comando_blastn = NcbiblastnCommandline( \
      query=args.q, db=args.d, \
      outfmt="'7 qseqid sseqid qcovs sstrand score qlen slen qstart qend sstart send pident'", out="blastn.tab",num_threads=10)
   
    # Executando
    stdout, stderr = comando_blastn()
   
    # Abrindo resultado
    if not args.o==None:
      blast_result = open(os.path.join(args.o, "blastn.tab"),"r")
    else:
      blast_result = open("blastn.tab","r")
    linhas=[]
    colunas=[]
    cont=0
    econt=0
    qry=''
    element=''
    
    print("Opening BLAST table...")    
    
    if not args.o==None:
      tabular=open(os.path.join(args.o, 'elements.txt'),'w')
    else:
      tabular=open('elements.txt','w')
    tabular.write("Query file: {}\n".format(args.q))
    tabular.write("Database file: {}\n".format(args.d))
    if args.min==None:
      min=5000
    else:
      min=int(args.min)
    if args.max==None:
      max=150000
    else:
      max=int(args.max)
    tabular.write("Element length: {0}-{1}\n".format(min,max))
    tabular.write("Query ID\tSubject ID\telement identification\telement 5' coordinate\telement 3' coordinate\telement length\n")
    
    for linha in blast_result.readlines():
      linha=linha.replace('\n','')
      if 'Query:' in linha:
        qry=linha.split(': ')[1]
        cont+=1
      if 'Fields:' in linha:
        colunas=linha.split(': ')[1]
        colunas=colunas.split(', ')
      if 'hits' in linha:
        if int(linha.split(' ')[1])==0:
          element='no'
          sbj='no hits'
          tabular.write("{0}\t{1}\t{2}\t\t\t\n".format(qry,sbj,element))
          print('Query '+qry+' has no hits!')
        else:
          print('Query '+qry+' has '+linha.split(' ')[1]+' hits!')
      if not '#' in linha:
        linhas.append(linha.split('\t'))
      if '# BLAST' in linha and not qry=='':
        estart=0
        eend=0
        elen=0
        if not linhas==[]:
          df=pd.DataFrame(columns=colunas,data=linhas)
          print("Searching element for "+qry+"...")
          df1=df[['query id','query length','subject id','% query coverage per subject','score','% identity']].drop_duplicates()
          df1[['query length','% query coverage per subject','score','% identity']] = df1[['query length','% query coverage per subject','score','% identity']].apply(pd.to_numeric)
          qlen=df1['query length'].tolist()[0]
          filtro='score'
          df1=df1.groupby('subject id').agg({filtro:'sum'}).reset_index()
          df1=df1.sort_values(by=[filtro],ascending=False)
          sbj=df1['subject id'].tolist()[0]
          df1=df.loc[(df['subject id'] == sbj)&(df['query id'] == qry)]
          if len(df1)==1:
            print("Query "+qry+" has no elements!")
            element='no'
            tabular.write("{0}\t{1}\t{2}\t\t\t\n".format(qry,sbj,element))
          else:
            df1[['% query coverage per subject','score','q. start','q. end']]=df1[['% query coverage per subject','score','q. start','q. end']].apply(pd.to_numeric)
            contigs=assembly(sorted(joinlists(df1['q. start'].tolist(),df1['q. end'].tolist())))
            if len(contigs)>1:
              element='yes'
              start,end=splitlist(contigs)
              estart,eend,elen=findelement(start,end)
              if elen>=min and elen<=max:
                econt+=1
                print("Query "+qry+" has element!")
                if not args.o==None:
                  os.mkdir(os.path.join(args.o, qry))
                  ft=open(os.path.join(args.o, qry, str(qry+'_element.gb')),'w')
                  fasta=open(os.path.join(args.o, qry, str(qry+'_element.fasta')),'w')
                else:
                  os.mkdir(qry)
                  ft=open(os.path.join(qry, str(qry+'_element.gb')),'w')
                  fasta=open(os.path.join(qry, str(qry+'_element.fasta')),'w')
                fasta.write(">{0} - element - {1}-{2}\n".format(qry,estart,eend))
                fasta.write(extract(qfastas[qry],estart,eend))
                fasta.close()
                print("Writing element's feature table...")
                ft.write("     misc_feature    {0}..{1}\n".format(eend,estart))
                ft.write("                     /label=element\n")
                ft.write("                     /color=255 8 234\n")
                ft.close()
                tabular.write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\n".format(qry,sbj,element,estart,eend,elen))
              else:
                print("Query "+qry+" has no elements!")
                element='no'
                tabular.write("{0}\t{1}\t{2}\t\t\t\n".format(qry,sbj,element))
            else:
              element='no'
              print("Query "+qry+" has no elements!")
              tabular.write("{0}\t{1}\t{2}\t\t\t\n".format(qry,sbj,element))
        linhas=[]
        colunas=[]
        qry=''
        element=''
    print('Processed '+str(cont)+' queries with '+str(econt)+' elements')
    tabular.close()
elif args.help == False:
  if args.q==None:
    print('ERROR: Missing query\n')
  if args.d == None:
    print('ERROR: Missing database file\n')
  print(ajuda)
    
