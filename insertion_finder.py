#!/usr/bin/python
import traceback
import warnings
import sys

from datetime import datetime
start_time = datetime.now()

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
elements=[]

ajuda = 'insertion_finder - element insertion finder in a genome through a BLAST search\n'
ajuda = ajuda + '(c) 2021. Arthur Gruber & Giuliana Pola\n'
ajuda = ajuda + 'Usage: insertion_finder.py -q <query file> -d <database file>\n'
ajuda = ajuda + '\tinsertion_finder.py -q <query file> -tab <BLASTn table file>\n'
ajuda = ajuda + '\nMandatory parameters:\n'
ajuda = ajuda + '-q <fasta or multifasta file>\tSequence to search with\n'
ajuda = ajuda + '-d <multifasta file>\tDatabase to BLAST against\n'
ajuda = ajuda + '-tab <table file>\tBLASTn search result table (fields: qseqid,sseqid,qcovs,qlen,slen,qstart,qend)\n'
ajuda = ajuda + '-run <local|web>\tchoice of running local or web BLAST search\n'
ajuda = ajuda + '\nOptional parameters:\n'
ajuda = ajuda + "-o <path>\tOutput directory (default: output_dir1)\n"
ajuda = ajuda + "-enddist <int>\tMaximum distance between block tip and query tip in base pairs(bp) (default: 50)\n"
ajuda = ajuda + "-minlen <int>\tMinimum element's length in base pairs(bp) (default: 5000)\n"
ajuda = ajuda + "-maxlen <int>\tMaximum element's length in base pairs(bp) (default: 50000)\n"
ajuda = ajuda + "-mincov <int>\tMinimum % query coverage per subject (default: 30)\n"
ajuda = ajuda + "-maxcov <int>\tMaximum % query coverage per subject (default: 90)\n"
ajuda = ajuda + "-cpu <int>\tNumber of threads to execute the blastn search (default: 10)\n"
ajuda = ajuda + "-color <int>\tThe RGB color of the element that is shown by the feature table, three integers between 0 and 255 separated by commas (default: 255,0,0)"

parser = argparse.ArgumentParser(add_help=False)
parser.add_argument('-q')
parser.add_argument('-d')
parser.add_argument('-o')
parser.add_argument('-color')
parser.add_argument('-cpu')
parser.add_argument('-tab')
parser.add_argument('-enddist')
parser.add_argument('-minlen')
parser.add_argument('-maxlen')
parser.add_argument('-mincov')
parser.add_argument('-maxcov')
parser.add_argument('-run')
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
    print("Missing query file (-q)!")
    valid=False
  elif not args.q==None:
   if not os.path.isfile(args.q):
     print("Query file (-q) not exist!")
     valid=False
   else:
     if isfasta(args.q):
       param['q']=args.q
     else:
       print("Invalid query file (-q)!")
       valid=False
  
  if valid==True:
    if args.run==None:
      if args.tab==None:
        print("Missing BLASTn table file (-tab)!")
        valid=False
      else:
        if not os.path.isfile(args.tab):
          print("BLASTn table file (-tab) not exist!")
          valid=False
        else:
          qid,colunas,hits,d=opentable(args.tab)
          if qid==[] or colunas==[]:
            print("Invalid BLASTn table (-tab)!")
            valid=False
          else:
            for col in ['query id', 'subject id', '% query coverage per subject', 'query length','q. start', 'q. end']:
              if not col in colunas:
                valid=False
                print("Column {} is not in table (-tab) {}!".format(col,args.tab))
          if valid==True:
            param['qid']=qid
            param['colunas']=colunas
            param['hits']=hits 
            param['d']=d
            param['tab']=args.tab
    elif args.run.lower()=='local':
      if args.d==None:
        print("Missing database file (-d)!")
        valid=False
      else:
        if not os.path.isfile(args.d):
          print("Database file (-d) not exist!")
          valid=False
        else:
          if isfasta(args.d):
            param['d']=args.d
          else:
            print("Invalid database file (-d), invalid formatting!")
            valid=False
  
  if valid==True:
    if args.tab==None:
      if args.run==None:
        print("Missing the choice (-run) of BLAST search: 'local' or 'web'!")
        valid=False
      else:
        try:
          args.run.lower()
        except:
          print("BLASTn run choice (-run) must be string: 'local' or 'web'!")
          valid=False
        else:
          if args.run.lower()=='web':
            param['run']='web'
          elif args.run.lower()=='local':
            param['run']='local'
          else:
            print("Invalid BLASTn run choice (-run), must be 'local' or 'web'!")
            valid=False
  


  if valid==True:
    if args.enddist==None:
      param['enddist']=50
    else:
      try:
        int(args.enddist)
      except:
        print("The maximum distance (-enddist) is not integer!")
        valid=False
      else:
        param['enddist']=int(args.enddist)
  
  if valid==True:
    if args.cpu==None:
      param['cpu']=10
    else:
      try:
        int(args.cpu)
      except:
        print("Number of threads (-cpu) is not integer!")
        valid=False
      else:
        if int(args.cpu)>=1:
          param['cpu']=int(args.cpu)
        else:
          print("Number of threads (-cpu) must be greater than or equal to 1!")
          valid=False
          
  
  if valid==True:
    if args.mincov==None and args.maxcov==None:
      param['mincov']=30
      param['maxcov']=90
    elif not args.mincov==None and not args.maxcov==None:
      try:
        int(args.mincov)
      except:
        print("Minimum query coverage (-mincov) is not integer!")
        valid=False
      try:
        int(args.maxcov)
      except:
        print("Maximum query coverage (-maxcov) is not integer!")
        valid=False
      if valid==True:
        if int(args.mincov)>100 or int(args.mincov)<0:
          print("Minimum query coverage (-mincov) must be between 0 and 100!")
          valid=False
        if int(args.maxcov)>100 or int(args.maxcov)<0:
          print("Maximum query coverage (-maxcov) must be between 0 and 100!")
          valid=False
        if valid==True and int(args.mincov)<int(args.maxcov):
          param['mincov']=args.mincov
          param['maxcov']=args.maxcov
        elif int(args.mincov)>int(args.maxcov):
          param['maxcov']=args.mincov
          param['mincov']=args.maxcov
        elif int(args.mincov)==int(args.maxcov):
          print("Minimum query coverage (-mincov) parameter must be smaller than the maximum query coverage (-maxcov)!")
          valid=False
    elif args.mincov==None and not args.maxcov==None:
      try:
        int(args.maxcov)
      except:
        print("Maximum query coverage (-maxcov) is not integer!")
        valid=False
      else:
        if int(args.maxcov)<=100 and int(args.maxcov)>=0:
          param['maxcov']=args.maxcov
        else:
          print("Maximum query coverage (-maxcov) must be between 0 and 100!")
          valid=False
    elif args.maxcov==None and not args.mincov==None:
      try:
        int(args.mincov)
      except:
        print("Minimum query coverage (-mincov) is not integer!")
        valid=False
      else:
        if int(args.mincov)<=100 and int(args.mincov)>=0:
          param['mincov']=args.mincov
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
      except:
        print("Minimum element length (-minlen) is not integer!")
        valid=False
      try:
        int(args.maxlen)
      except:
        print("Maximum element length (-maxlen) is not integer!")
        valid=False
      if valid==True:
        if int(args.minlen)<0:
          print("Minimum element length (-minlen) must be greater than or equal to 0!")
          valid=False
        if int(args.maxlen)<0:
          print("Maximum element length (-maxlen) must be greater than or equal to 0!")
          valid=False
        if valid==True and int(args.minlen)<int(args.maxlen):
          param['minlen']=args.minlen
          param['maxlen']=args.maxlen
        elif int(args.minlen)>int(args.maxlen):
          param['maxlen']=args.minlen
          param['minlen']=args.maxlen
        elif int(args.minlen)==int(args.maxlen):
          print("Minimum element length (-minlen) must be smaller than the maximum element length (-maxlen)!")
          valid=False
    elif args.minlen==None and not args.maxlen==None:
      try:
        int(args.maxlen)
      except:
        print("Maximum element length (-maxlen) is not integer!")
        valid=False
      else:
        if int(args.maxlen)>=0:
          param['maxlen']=args.maxlen
        else:
          print("Maximum element length (-maxlen) must be greater than or equal to 0!")
          valid=False
    elif args.maxlen==None and not args.minlen==None:
      try:
        int(args.minlen)
      except:
        print("Minimum element length (-minlen) is not integer!")
        valid=False
      else:
        if int(args.minlen)>=0:
          param['minlen']=args.minlen
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
        except:
          print("The RGB color of the element (-color) must be three integers between 0 and 255!")
          valid=False
    
  if valid==True:
   if not args.o==None:
     if os.path.isdir(args.o):
       print("Output directory (-o) {} exists!".format(args.o))
       i=1
       out=args.o
       while os.path.isdir(out):
         i+=1
         out=args.o+"_"+str(i)
       os.mkdir(out)
       print("Creating output directory (-o) {}...".format(out))
       param['o']=out
     else:
       try:
         os.mkdir(args.o)
         print("Creating output directory (-o) {}...".format(args.o))
         param['o']=args.o
       except:
         print("Output directory (-o) not valid\n")
         valid=False
   else:
     i=1
     out='output_dir'
     while os.path.isdir(out):
       i+=1
       out='output_dir'+str(i)
     os.mkdir(out)
     print("Creating output directory (-o) {}...".format(out))
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
    elif start>=enddist and end>=enddist:
      return 0,0,0
    elif start<=enddist and end<=enddist:
      return -1,-1,-1

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
        if param['run']=='local':
          comando_blastn = NcbiblastnCommandline(query=param['q'], db=param['d'],outfmt="'7 qseqid sseqid qcovs qlen slen qstart qend'", out=os.path.join(param['o'], "blastn.tab"),num_threads=param['cpu'])
        elif param['run']=='web':
          comando_blastn = NcbiblastnCommandline(query=param['q'], db="nt", outfmt="'7 qseqid sseqid qcovs qlen slen qstart qend'", out=os.path.join(param['o'], "blastn.tab"),remote=True)
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
      log=open(os.path.join(param['o'], 'file.log'),'w')
      hits=param.pop('hits')
      colunas=param.pop('colunas')
      qids=param.pop('qid')
      tabular.write("Query file: {}".format(param["q"]))
      log.write("Query file: {}".format(param["q"]))
      if not args.tab==None:
        tabular.write("\nBlastn table file: {}".format(args.tab))
        log.write("\nBlastn table file: {}".format(args.tab))
      if 'd' in param:
        tabular.write("\nDatabase file: {}".format(param["d"]))
        log.write("\nDatabase file: {}".format(param["d"]))
      tabular.write("\nElement length: {}-{}".format(param["minlen"],param["maxlen"]))
      log.write("\nElement length: {}-{}".format(param["minlen"],param["maxlen"]))
      tabular.write("\nQuery coverage: {}-{}".format(param["mincov"],param["maxcov"]))
      log.write("\nQuery coverage: {}-{}".format(param["mincov"],param["maxcov"]))
      tabular.write("\nMaximum block distance: {}".format(param["enddist"]))
      log.write("\nMaximum block distance: {}".format(param["enddist"]))
      log.write("\nNumber of threads: {}".format(param["cpu"]))
      tabular.write("\nNumber of threads: {}".format(param["cpu"]))
      tabular.write("\nElement color: {}\n".format(param["color"]))
      log.write("\nElement color: {}".format(param["color"]))
      tabular.write("\nQuery ID\tSubject ID\tElement identification\tElement 5' coordinate\tElement 3' coordinate\tElement length\tValid")
      df=pd.DataFrame(columns=colunas,data=hits)
      df=df.apply(pd.to_numeric, errors='ignore').drop_duplicates()
      for qid in qids:
        df1=df.loc[df['query id'] == qid]
        log.write("\n\nQuery id: {}".format(qid))
        log.write("\n{} hits!".format(len(df1)))
        qlen=df1['query length'].tolist()[0]
        log.write("\nQuery length: {}\n".format(qlen))
        df1=df1.groupby('subject id').agg({'% query coverage per subject':'mean'}).reset_index()
        df1=df1.sort_values(by=['% query coverage per subject'],ascending=False)
        if len(df1)==0:
          element='no'
          sid='no hits'
          tabular.write("\n{0}\t{1}\t{2}\t\t\t\t{3}".format(qid,sid,element,'no'))
        elif len(df1)>0:
          i=0
          nsid=len(df1['subject id'].tolist())
          while not i==-1:
            hit=[]
            sid=df1['subject id'].tolist()[i]
            cov=df1['% query coverage per subject'].tolist()[i]
            hit.append('Subject id: {}'.format(sid))
            hit.append('% query coverage per subject: {}'.format(cov))
            df2=df.loc[(df['subject id'] == sid)&(df['query id'] == qid)]
            reads=sorted(joinlists(df2['q. start'].tolist(),df2['q. end'].tolist()))
            contigs=assembly(reads)
            hit.append("Alignments: {}".format(str(reads)))
            if len(contigs)==1:
              if cov>=param['mincov'] and cov<=param['maxcov']:
                estart,eend,elen=oneblock(contigs,qlen,param['enddist'])
                #print(estart,eend,elen)
                if estart>0:
                  log.write('\n'.join(str(v) for v in hit))
                  i=-1                  
                  element='yes'
                  log.write("\nElement coodinates: {} - {}".format(estart,eend))
                  log.write("\nElement length: {}".format(elen))
                  if elen>=param['minlen'] and elen<=param['maxlen']:
                    log.write("\nValid element!")
                    econt+=1
                    os.mkdir(os.path.join(param['o'], qid))
                    tabular.write("\n{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}".format(qid,sid,element,estart,eend,elen,'yes'))
                    try:
                      ft=open(os.path.join(param['o'], qid, str(qid+'_element.gb')),'w')
                      ft.write("     misc_feature    {0}..{1}\n".format(estart,eend))
                      ft.write("                     /label=element\n")
                      ft.write("                     /color={} {} {}\n".format(param['color'][0],param['color'][1],param['color'][2]))
                      ft.close()
                    except:
                      log.write("\nElement's feature table was't writen!")
                    else:
                      log.write("\nWriting element's feature table...")
                    try:
                      fasta=open(os.path.join(param['o'], qid, str(qid+'_element.fasta')),'w')
                      fasta.write(">{0} - element - {1}-{2}\n".format(qid,estart,eend))
                      fasta.write(extract(qseqs, qid, estart, eend))
                      fasta.close()
                    except:
                      log.write("\nElement's fasta was't writen!")
                    else:
                      log.write("\nWriting element's fasta...")
                elif estart==0:
                  log.write('\n'.join(str(v) for v in hit))
                  log.write("\nInvalid element!")
                  tabular.write("\n{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}".format(qid,sid,element,estart,eend,elen,'no'))
                  i=-1
              elif cov<=param['mincov']:
                log.write('\n'.join(str(v) for v in hit))
                element='no'
                tabular.write("\n{0}\t{1}\t{2}\t\t\t\t{3}".format(qid,sid,element,'no'))
                i=-1
            elif len(contigs)>1:
              if cov>=param['mincov'] and cov<=param['maxcov']:
                element='yes'
                start,end=splitlist(contigs)
                estart,eend,elen=findelement(start,end)
                log.write('\n'.join(str(v) for v in hit))
                log.write("\nElement coodinates: {} - {}".format(estart,eend))
                log.write("\nElement length: {}".format(elen))
                i=-1
                if elen>=param['minlen'] and elen<=param['maxlen']:
                  log.write("\nValid element!")
                  econt+=1
                  os.mkdir(os.path.join(param['o'], qid))
                  tabular.write("\n{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}".format(qid,sid,element,estart,eend,elen,'yes'))
                  try:
                    log.write("\nWriting element's feature table...")
                    ft=open(os.path.join(param['o'], qid, str(qid+'_element.gb')),'w')
                    ft.write("     misc_feature    {0}..{1}\n".format(estart,eend))
                    ft.write("                     /label=element\n")
                    ft.write("                     /color={} {} {}\n".format(param['color'][0],param['color'][1],param['color'][2]))
                    ft.close()
                  except:
                    log.write("\nElement's feature table was't writen!")
                  try:
                    log.write("\nWriting element's fasta...")
                    fasta=open(os.path.join(param['o'], qid, str(qid+'_element.fasta')),'w')
                    fasta.write(">{0} - element - {1}-{2}\n".format(qid,estart,eend))
                    fasta.write(extract(qseqs, qid, estart, eend))
                    fasta.close()
                  except:
                    log.write("\nElement's fasta was't writen!")
                else:
                  log.write("\nInvalid element!")
                  tabular.write("\n{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}".format(qid,sid,element,estart,eend,elen,'no'))
            if i==nsid and not i==-1:
              i=-1
              element='no'
              tabular.write("\n{0}\t{1}\t{2}\t\t\t\t{3}".format(qid,sid,element,'no'))
            elif not i==-1:
              i+=1
      tabular.write('\n\nProcessed {} queries with {} valid elements'.format(len(qids),econt))
      log.write('\n\nProcessed {} queries with {} valid elements'.format(len(qids),econt))
      print('Processed {} queries with {} valid elements'.format(len(qids),econt))
      end_time = datetime.now()
      print('Duration: {}'.format(end_time - start_time))
      log.write('\nDuration: {}'.format(end_time - start_time))
      tabular.write('\nDuration: {}'.format(end_time - start_time))
      log.close()
      tabular.close()