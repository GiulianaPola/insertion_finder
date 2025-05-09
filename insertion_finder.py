#!/usr/bin/env python3
import traceback
import warnings
def warn_with_traceback(message, category, filename, lineno, file=None, line=None):

    log = file if hasattr(file,'write') else sys.stderr
    traceback.print_stack(file=log)
    log.write(warnings.formatwarning(message, category, filename, lineno, line))

warnings.showwarning = warn_with_traceback

import argparse
import sys
from datetime import datetime
start_time = datetime.now()
import os
import re
param=dict()

start_time = datetime.now()
call=os.path.abspath(os.getcwd())

version="2.3.1"

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

def isfasta(filename):
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
      newname=os.path.join(path, str(name+str(i)))
  elif typ=='file':
    while os.path.isfile(newname):
      i+=1
      newname=os.path.join(path, str(name+str(i)))
  return newname

def blastparse(tab):
  hits=[]
  qid=[]
  columns=''
  file = open(tab, "r")
  for line in file.readlines():
    line=line.replace("  "," ")
    line=line.replace("\n","")
    if not re.search('#', line):
      hits.append(line.split('\t'))
    elif re.search('# Fields: ', line):
      columns=line.split('# Fields: ')[-1]
      columns=columns.replace("\n","").split(', ')
    if re.search('# Query: ', line):
      qid.append(line.split('# Query: ')[-1])
  return qid,columns,hits

def validateargs(args):
  valid=True
  
  if args.q==None:
    print("Missing query file (-q)!")
    valid=False
  elif not args.q==None:
   if not os.path.isfile(args.q):
     print("Query file (-q) not exist!")
     valid=False
   else:
     if isfasta(args.q):
       param['q']=os.path.realpath(args.q)
     else:
       print("Invalid query file (-q)!")
       valid=False

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

  if args.run.lower()=='local':
    if args.d==None:
      print("Missing database file (-d)!")
      valid=False
    else:
      if not os.path.isfile(args.d):
        print("Database file (-d) not exist!")
        valid=False
      else:
        if isfasta(args.d):
          param['d']=os.path.realpath(args.d)
        else:
          print("Invalid database file (-d), invalid formatting!")
          valid=False

  if valid==True and not args.run==None:
    if args.run.lower()=='local':
      if args.cpu==None:
        param['cpu']=18
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
    if not args.tab==None:
      if not os.path.isfile(args.tab):
        print("BLASTn table file (-tab) not exist!")
        valid=False
      else:
        missing=[]
        qid,columns,hits=blastparse(args.tab)
        if qid==[] or columns==[]:
          print("Invalid BLASTn table (-tab)!")
          valid=False
        else:
          for col in ['query id', 'subject id', '% query coverage per subject', 'query length','q. start', 'q. end']:
            if not col in columns:
              valid=False
              missing.append(col)                
        if valid==False:
          if len(missing)==1:
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
      for id in args.org.split(','):
        int(id)
        taxids.append("txid{}[ORGN]".format(id))
    except:
      print("Taxid(s) must be integers separated by commas!")
      valid=False
    else:
      param['org']=' OR '.join(taxids)

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
          param['mincov']=int(args.mincov)
          param['maxcov']=int(args.maxcov)
        elif int(args.mincov)>int(args.maxcov):
          param['maxcov']=int(args.mincov)
          param['mincov']=int(args.maxcov)
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
          param['maxcov']=int(args.maxcov)
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
          param['minlen']=int(args.minlen)
          param['maxlen']=int(args.maxlen)
        elif int(args.minlen)>int(args.maxlen):
          param['maxlen']=int(args.minlen)
          param['minlen']=int(args.maxlen)
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
          param['maxlen']=int(args.maxlen)
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
        except:
          print("The RGB color of the element (-color) must be three integers between 0 and 255!")
          valid=False
    
  if valid==True:
   if not args.out==None:
     if os.path.isdir(args.out):
       print("Output directory (-out) {} exists!".format(args.out))
       out=rename(1,args.out,'dir')
       os.mkdir(out)
       print("Creating output directory (-out) {}...".format(out))
       param["out"]=out
     else:
       try:
         os.mkdir(args.out)
         print("Creating output directory (-out) {}...".format(args.out))
         param["out"]=args.out
       except:
         print("Output directory (-out) not valid\n")
         valid=False
   else:
     out=rename(1,'output_dir','dir')
     os.mkdir(out)
     print("Creating output directory (-out) {}...".format(out))
     param["out"]=out
    
  return valid,param

def blast(param, q,blast_time='',seqs=[]):
  tab=os.path.join(param["out"], 'blastn2.tab')
  blast_start = datetime.now()
  from Bio.Blast.Applications import NcbiblastnCommandline
  if param['run']=='local':
    print("Starting BLAST search...")
    comando_blastn = NcbiblastnCommandline(query=q, db=param['d'],outfmt="7 qseqid sseqid qcovs qlen slen qstart qend", out=tab,num_threads=param['cpu'])
    stdout, stderr = comando_blastn()
  elif param['run']=='web':
    if not 'org' in param:
      print("Starting BLAST search...")
      comando_blastn = NcbiblastnCommandline(query=q, db="nt", outfmt="7 qseqid sseqid qcovs qlen slen qstart qend", out=tab,remote=True,task='megablast')
      stdout, stderr = comando_blastn()
    else:
      print("Starting BLAST search...")
      comando_blastn = NcbiblastnCommandline(query=q, db="nt", outfmt="7 qseqid sseqid qcovs qlen slen qstart qend", out=tab,remote=True,entrez_query="'{}'".format(param['org']),task='megablast')
      stdout, stderr = comando_blastn()
  if blast_time=='':
    blast_time=(datetime.now() - blast_start)
  else:
    blast_time+=(datetime.now() - blast_start)
  if seqs==[]:
    seqs=missingquery(tab=tab,qry=q)
  else:
    seqs=missingquery(tab=tab,qry=q,seqs=seqs)
  if not seqs==False:
    blast(param, os.path.join(param["out"], 'newquery.fasta'),blast_time,seqs)
  else:
   if os.path.isfile(os.path.join(param["out"], 'newquery.fasta')):
     os.remove(os.path.join(param["out"], 'newquery.fasta'))
   if os.path.isfile(os.path.join(param["out"], 'blastn.tab')):
     os.remove(os.path.join(param["out"], 'blastn2.tab'))
   else:
     os.rename(os.path.join(param["out"], 'blastn2.tab'), os.path.join(param["out"], 'blastn.tab'))
   return blast_time

def missingquery(tab,qry,seqs=[]):
  qfile=open(qry,"r").read()
  tabfile=open(tab,"r").read()
  tabqs=tabfile.split("# BLAST")
  hits=dict()
  misseqs=[]
  newq=''
  table=''
  tabqs=tabqs[:-1]
  tabfile="# BLAST".join(tabqs)
  for q in qfile.split(">"):
    if not q=='':
     qid=q.split("\n")[0]
     qseqs='\n'.join(q.split("\n")[1:])
     if qid not in tabfile:
       misseqs.append(qid)
       #print(qid,qseqs)
       newq+=">{}\n{}".format(qid,qseqs)
     else:
       seqs.append(qid)
  if misseqs==[]:
    if os.path.isfile(os.path.join(param["out"], 'blastn.tab')):
      with open(os.path.join(param["out"], 'blastn.tab'),'a+') as f:
        if "# BLAST processed " not in f.read():
          print("Complete blast table!")
          f.write(tabfile)
          f.write('\n# BLAST processed {} queries'.format(len(seqs)))
    return False
  else:
    if len(misseqs)==1:
      print("BLAST search missing query: {}!".format(misseqs[0]))
    else:
      print("BLAST search missing queries: {}!".format(", ".join(misseqs)))
    newqfile=(os.path.join(param["out"], 'newquery.fasta'),'w')
    newqfile.write(newq)
    newqfile.close()
    if os.path.isfile(os.path.join(param["out"], 'blastn.tab')):
      newtabfile=open(os.path.join(param["out"], 'blastn.tab'),'a')
    else:
      newtabfile=open(os.path.join(param["out"], 'blastn.tab'),'w')
    newtabfile.write(tabfile)
    return seqs    

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
    print(help)
elif args.help == True:
    print(help)
elif args.version == True:
    print(version)
else:
  import pandas as pd
  from Bio import SeqIO
  import os
  log=[]
  elements=[]
  valid,param=validateargs(args)
  #print(param)
  if valid==False:
    print(help)
  else:
    print('Valid arguments!')
    try:
      log=open(os.path.join(param["out"], 'file.log'),'w')
      log.write('insertion_finder v{}'.format(version))
    except:
      print('Log file was not created: {}'.format(error))
    else:  
      log.write('(c) 2021. Giuliana Pola & Arthur Gruber\n')
      log.write('For more information access: https://github.com/GiulianaPola/insertion_finder')
      log.write('\nStart time: {}\n'.format(start_time.strftime("%d/%m/%Y, %H:%M:%S")))
      log.write('\nWorking directory: {}\n'.format(call))
      log.write('\nCommand line: {}\n'.format(' '.join(sys.argv)))
      user=""
      try:
        user=os.getlogin()
      except Exception as e:
        try:
          user=os.environ['LOGNAME']
        except Exception as e:
          try:
            user=os.environ['USER']
          except Exception as e:
            pass
          else:
            pass
        else:
          pass
      else:
        pass
      if not user=="":
        log.write('\nUser: {}\n'.format(user))
      log.write('\nParameters:\n')
      for arg in vars(args):
          value = getattr(args, arg)
          if value is not None and value is not False:
              log.write("{}={}\n".format(arg,value))
      log.flush() 
      qseqs=open(param['q'], "r").read()
      if args.tab==None:
        log.write('\nStarting BLASTn search...')
        blast_time=blast(param, param['q'])
        print('BLASTn search execution time: {}'.format(blast_time))
        log.write('\nBLASTn search execution time: {}'.format(blast_time))
        param['tab']=os.path.join(param["out"], "blastn.tab")
        param['qid'],param['columns'],param['hits']=blastparse(param['tab'])
      else:
        seqs=missingquery(tab=param['tab'],qry=param['q'])
        if not seqs==False:
          blast_time=blast(param, os.path.join(param["out"], 'newquery.fasta'),seqs=seqs)
          param['tab']=os.path.join(param["out"], "blastn.tab")
          param['qid'],param['columns'],param['hits']=blastparse(param['tab'])
      if 'hits' not in param:
        print("Table BLASTn was not read!")
        log.write("\nTable BLASTn was not read!")
      else:
        print("Opening BLASTn table...")
        log.write("\nOpening BLASTn table...")  
        cont=0
        econt=0
        try:
          tabular=open(os.path.join(param["out"], 'elements.txt'),'w')
        except:
          print('Elements table was not created!')
          log.write('\nElements table was not created!')
        else:
          tabular.write('insertion_finder v{}'.format(version))
          tabular.write("\nQuery file: {}".format(param["q"]))
          log.write("\n\nQuery file: {}".format(param["q"]))
          if not args.tab==None:
            tabular.write("\nBlastn table file: {}".format(args.tab))
            log.write("\nBlastn table file: {}".format(args.tab))
          if 'd' in param:
            tabular.write("\nDatabase file: {}".format(param["d"]))
            log.write("\nDatabase file: {}".format(param["d"]))
          else:
            tabular.write("\nDatabase file: nt")
            log.write("\nDatabase file: nt")
          if args.tab==None and not args.run==None:
            log.write('\nBLASTn search execution time: {}'.format(blast_time))
            tabular.write('\nBLASTn search execution time: {}'.format(blast_time))
          if 'org' in param:
            tabular.write("\nTaxids: {}".format(args.org))
            log.write("\nTaxids: {}".format(args.org))
          tabular.write("\nElement length: {}-{}".format(param["minlen"],param["maxlen"]))
          log.write("\nElement length: {}-{}".format(param["minlen"],param["maxlen"]))
          tabular.write("\nQuery coverage: {}-{}".format(param["mincov"],param["maxcov"]))
          log.write("\nQuery coverage: {}-{}".format(param["mincov"],param["maxcov"]))
          tabular.write("\nMaximum block distance: {}".format(param["enddist"]))
          log.write("\nMaximum block distance: {}".format(param["enddist"]))
          if 'cpu' in param:
            log.write("\nNumber of threads: {}".format(param["cpu"]))
            tabular.write("\nNumber of threads: {}".format(param["cpu"]))
          tabular.write("\nElement color: {}\n".format(param["color"]))
          log.write("\nElement color: {}".format(param["color"]))
          tabular.write("\nQuery ID\tSubject ID\tElement identification\tElement 5' coordinate\tElement 3' coordinate\tElement length\tValid")
          #print("Columns {},{}".format(len(param['columns']),param['columns'][0]))
          #print("Hits {},{}".format(len(param['hits']),param['hits'][0]))
          df=pd.DataFrame(columns=param['columns'],data=param['hits'])
          df = df.apply(pd.to_numeric, errors='coerce').drop_duplicates()
          for qid in param['qid']:
            df1=df.loc[df['query id'] == qid]
            log.write("\n\nQuery id: {}".format(qid))
            log.write("\n{} hits!".format(len(df1)))
            if len(df1)==0:
              element='no'
              sid='no hits'
              tabular.write("\n{0}\t{1}\t{2}\t\t\t\t{3}".format(qid,sid,element,'no'))
            elif len(df1)>0:
              qlen=df1['query length'].tolist()[0]
              log.write("\nQuery length: {}".format(qlen))
              df1=df1.groupby('subject id').agg({'% query coverage per subject':'mean'}).reset_index()
              df1=df1.sort_values(by=['% query coverage per subject'],ascending=False)
              i=0
              nsid=len(df1['subject id'].tolist())
              while i>-1 and i<nsid:
                hit=[]
                sid=df1['subject id'].tolist()[i]
                cov=df1['% query coverage per subject'].tolist()[i]
                hit.append('\nSubject id: {}'.format(sid))
                hit.append('% query coverage per subject: {}'.format(cov))
                df2=df.loc[(df['subject id'] == sid)&(df['query id'] == qid)]
                reads=sorted(joinlists(df2['q. start'].tolist(),df2['q. end'].tolist()))
                contigs=assembly(reads)
                hit.append("Alignments: {}".format(str(contigs)))
                if len(contigs)==1:
                  if cov>=param['mincov']:
                    estart,eend,elen=oneblock(contigs,qlen,param['enddist'])
                    #print(estart,eend,elen)
                    if estart==0 and i==nsid-1:
                      log.write("\nInvalid element, block distance greater than the maximum distance on both sides!")
                      log.write("\nNo valid hits!")
                      cont+=1
                      tabular.write("\n{0}\t{1}\t{2}\t\t\t\t{3}".format(qid,'no valid hits','no','no'))
                    if estart==-1 and i==nsid-1:
                      log.write("\nInvalid element, block distance smaller than the maximum distance on both sides!")
                      log.write("\nNo valid hits!")
                      cont+=1
                      tabular.write("\n{0}\t{1}\t{2}\t\t\t\t{3}".format(qid,'no valid hits','no','no'))
                    if estart>0:
                      log.write('\n'.join(str(v) for v in hit))
                      i=-1      
                      cont+=1            
                      log.write("\nElement coodinates: {} - {}".format(estart,eend))
                      log.write("\nElement length: {}".format(elen))
                      if elen>=param['minlen'] and elen<=param['maxlen']:
                        log.write("\nValid element!")
                        econt+=1
                        os.mkdir(os.path.join(param["out"], qid))
                        tabular.write("\n{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}".format(qid,sid,'yes',estart,eend,elen,'yes'))
                        try:
                          ft=open(os.path.join(param["out"], qid, str(qid+'_element.gb')),'w')
                          ft.write("     misc_feature    {0}..{1}\n".format(estart,eend))
                          ft.write("                     /label=element\n")
                          ft.write("                     /color={} {} {}\n".format(param['color'][0],param['color'][1],param['color'][2]))
                          ft.close()
                        except:
                          log.write("\nElement's feature table was't writen!")
                        else:
                          log.write("\nWriting element's feature table...")
                        try:
                          fasta=open(os.path.join(param["out"], qid, str(qid+'_element.fasta')),'w')
                          fasta.write(">{0} - element - {1}-{2}\n".format(qid,estart,eend))
                          fasta.write(extract(qseqs, qid, estart, eend))
                          fasta.close()
                        except:
                          log.write("\nElement's fasta was't writen!")
                        else:
                          log.write("\nWriting element's fasta...")
                      else:
                        if elen<=param['minlen']:
                          log.write("\nInvalid element, smaller than valid size!")
                        if elen>=param['maxlen']:
                          log.write("\nInvalid element, larger than valid size!")
                        tabular.write("\n{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}".format(qid,sid,'yes',estart,eend,elen,'no'))
                  elif cov<=param['mincov']:
                    log.write('\n'.join(str(v) for v in hit))
                    log.write('\nInvalid element, % query coverage less than valid coverage!')
                    tabular.write("\n{0}\t{1}\t{2}\t\t\t\t{3}".format(qid,sid,'no','no'))
                    i=-1
                    cont+=1
                elif len(contigs)>1:
                  if cov<=param['maxcov']:
                    start,end=splitlist(contigs)
                    estart,eend,elen=findelement(start,end)
                    log.write('\n'.join(str(v) for v in hit))
                    log.write("\nElement coodinates: {} - {}".format(estart,eend))
                    log.write("\nElement length: {}".format(elen))
                    i=-1
                    cont+=1
                    if elen>=param['minlen'] and elen<=param['maxlen']:
                      log.write("\nValid element!")
                      econt+=1
                      os.mkdir(os.path.join(param["out"], qid))
                      tabular.write("\n{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}".format(qid,sid,'yes',estart,eend,elen,'yes'))
                      try:
                        log.write("\nWriting element's feature table...")
                        ft=open(os.path.join(param["out"], qid, str(qid+'_element.gb')),'w')
                        ft.write("     misc_feature    {0}..{1}\n".format(estart,eend))
                        ft.write("                     /label=element\n")
                        ft.write("                     /color={} {} {}\n".format(param['color'][0],param['color'][1],param['color'][2]))
                        ft.close()
                      except:
                        log.write("\nElement's feature table was't writen!")
                      try:
                        log.write("\nWriting element's fasta...")
                        fasta=open(os.path.join(param["out"], qid, str(qid+'_element.fasta')),'w')
                        fasta.write(">{0} - element - {1}-{2}\n".format(qid,estart,eend))
                        fasta.write(extract(qseqs, qid, estart, eend))
                        fasta.close()
                      except:
                        log.write("\nElement's fasta was't writen!")
                    else:
                      cont+=1
                      i=-1
                      if elen<=param['minlen']:
                        log.write("\nInvalid element, smaller than valid size!")
                      if elen>=param['maxlen']:
                        log.write("\nInvalid element, larger than valid size!")
                      tabular.write("\n{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}".format(qid,sid,'yes',estart,eend,elen,'no'))
                if not i==-1:
                  i+=1
          tabular.write('\n\nProcessed {} queries with {} valid elements'.format(len(param['qid']),econt))
          log.write('\n\nProcessed {} queries with {} valid elements'.format(len(param['qid']),econt))
          print('Processed {} queries with {} valid elements'.format(len(param['qid']),econt))
          log.write('\nProgram execution time: {}'.format(datetime.now() - start_time))
          tabular.write('\nProgram execution time: {}'.format(datetime.now() - start_time))
          log.close()
          tabular.close()
print('Program execution time: {}'.format(datetime.now() - start_time))
print('\a')
