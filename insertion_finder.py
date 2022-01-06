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

ajuda = 'InsertionFinder - element insertion locator in a genome through a BLAST search\n'
ajuda = ajuda + '(c) 2021. Arthur Gruber & Giuliana Pola\n'
ajuda = ajuda + 'Usage: insertion_finder.py -query <sequence to search with> -db <database to BLAST against> -out <output file> \n'
ajuda = ajuda + '\nMandatory parameters:\n'
ajuda = ajuda + '-query\tSequence to search with\n'
ajuda = ajuda + '-db\tDatabase to BLAST against\n'
ajuda = ajuda + '\nOptional parameters:\n'
ajuda = ajuda + "-out\tOutput file (default: 'BLASTn_'+query filename+'.txt', 'elements_'+query filename+'.txt', query id+'_element.gb', query id+'_element.fasta')\n"
ajuda = ajuda + "-minlen\tMinimum element's length (default: >=5000)\n"
ajuda = ajuda + "-maxlen\tMaximum element's length (default: <=150.000)\n"

parser = argparse.ArgumentParser(add_help=False)
parser.add_argument('-query')
parser.add_argument('-db')
parser.add_argument('-out')
parser.add_argument('-minlen')
parser.add_argument('-maxlen')
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

def assembly(contigs):
  contigs.sort(key = lambda x: x[0])
  m = []
  s = -10000
  max = -100000
  for i in range(len(contigs)):
      a = contigs[i]
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

def findelement(start,end,qlen):
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

if args.help == True:
    print(ajuda)
elif not args.query==None and not args.db==None:
  if args.out == None:
    if not args.query==None:
      if args.query.count("/")>0:
        nome=args.query.split("/")[-1]
        nome=nome.split(".")[0]
        saida='BLASTn_'+nome+".txt"
      else:
        nome=args.query.split(".")[0]
        saida='BLASTn_'+nome+".txt"
  else:
    if args.out.count(".")>0:
      nome=nome.split(".")[0]
    else:
      nome=args.out
    saida=args.out
    
  if args.db==None:
    db="nt"
  else:
    db=args.db

  qfastas=splitfile(args.query)
  #print(sorted(qfastas.keys()))

  comando_blastn = NcbiblastnCommandline( \
  query=args.query, db=db, \
  outfmt="'7 qseqid sseqid qcovs sstrand score qlen slen qstart qend sstart send pident'", out=str(saida),num_threads=10)
 
  # Executando
  stdout, stderr = comando_blastn()
 
  # Abrindo resultado
  blast_result = open(str(saida),"r")
  linhas=[]
  colunas=[]
  cont=0
  econt=0
  qry=''
  element=''
  
  print("Opening BLAST table...")    
  
  tabular=open(str('elements.txt'),'w')
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
        df1[['% query coverage per subject','score','q. start','q. end']]=df1[['% query coverage per subject','score','q. start','q. end']].apply(pd.to_numeric)
        df1=df1.sort_values(by=['q. start','q. end'])
        contigs=assembly(joinlists(df1['q. start'].tolist(),df1['q. end'].tolist()))
        if len(contigs)>1:
          element='yes'
          start,end=splitlist(contigs)
          estart,eend,elen=findelement(start,end,qlen)
          if args.minlen==None:
            minlen=5000
          else:
            minlen=int(args.minlen)
          if args.maxlen==None:
            maxlen=150000
          else:
            maxlen=int(args.maxlen)
          if elen>=minlen and elen<=maxlen:
            econt+=1
            print("Query "+qry+" has element!")
            fasta=open(str(qry+'_element.fasta'),'w')
            fasta.write(">{0} - element - {1}-{2}\n".format(qry,estart,eend))
            fasta.write(extract(qfastas[qry],estart,eend))
            fasta.close()
            ft=open(str(qry+'_element.gb'),'w')
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
  if args.query==None:
    print('ERROR: Missing query or BLASTn table file\n')
  if args.db == None:
    print('ERROR: Missing database file\n')
  print(ajuda)
    
