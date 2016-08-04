#!/usr/bin/env python
""" MAGeCK count module
Copyright (c) 2014 Wei Li, Han Xu, Xiaole Liu lab 
This code is free software; you can redistribute it and/or modify it
under the terms of the BSD License (see the file COPYING included with
the distribution).
@status:  experimental
@version: $Revision$
@author:  Wei Li 
@contact: li.david.wei AT gmail.com
"""
from __future__ import print_function

import sys
import argparse
import math
import logging
import string
from mageck.testVisualCount import *

def mageckcount_parseargs():
  """
  Parse arguments. Only used when mageckCount.py is executed directly.
  """
  parser=argparse.ArgumentParser(description='Collecting read counts for multiple samples.')
  
  parser.add_argument('-l','--list-seq',required=True,help='A file containing the list of sgRNA names, their sequences and associated genes. Support file format: csv and txt.')
  parser.add_argument('--sample-label',default='',help='Sample labels, separated by comma (,). Must be equal to the number of samples provided. Default "sample1,sample2,...".')
  parser.add_argument('-n','--output-prefix',default='sample1',help='The prefix of the output file(s). Default sample1.')
  parser.add_argument('--trim-5',type=int,default=0,help='Length of trimming the 5\' of the reads. Default 0')
  parser.add_argument('--sgrna-len',type=int,default=20,help='Length of the sgRNA. Default 20')
  parser.add_argument('--count-n',action='store_true',help='Count sgRNAs with Ns. By default, sgRNAs containing N will be discarded.')
  parser.add_argument('--fastq',nargs='+',help='Sample fastq files, separated by space; use comma (,) to indicate technical replicates of the same sample. For example, "--fastq sample1_replicate1.fastq,sample1_replicate2.fastq sample2_replicate1.fastq,sample2_replicate2.fastq" indicates two samples with 2 technical replicates for each sample.')

  
  args=parser.parse_args()
  
  
  return args

def mageckcount_checkargs(args):
  """
  Check args
  """
  if args.sample_label!='':
    nlabel=args.sample_label.split(',')
    #nfq=args.fastq.split(',')
    nfq=(args.fastq)
    if len(nlabel)!=len(nfq):
      logging.error('The number of labels ('+str(nlabel)+') must be equal to the number of fastq files provided.')
      sys.exit(-1)
  return 0

def mageckcount_gettotalnormfactor(ctable):
  """
  Get the factor by total normalization
  """
  n=len(ctable[list(ctable.keys())[0]]) # samples
  m=len(ctable) # sgRNAs
  # calculate the sum
  sumsample=[0]*n
  for (k,v) in ctable.items():
    sumsample=[sumsample[i]+v[i] for i in range(n)]
  # normalizing factor
  avgsample=sum(sumsample)/float(n)
  samplefactor=[avgsample/k for k in sumsample]
  return samplefactor

def mageckcount_getmediannormfactor(ctable):
  """
  Get the factor by median normalization
  """
  n=len(ctable[list(ctable.keys())[0]]) # samples
  m=len(ctable) # sgRNAs
  meanval={k:math.exp( (sum( [ math.log(v2+1.0) for v2 in v])*1.0/n) ) for (k,v) in ctable.items() if sum(v)>0}  # geometric mean
  meanval={k:(lambda x: x if x>0 else 1)(v) for (k,v) in meanval.items()} # delete those with all 0 read counts
  #samplefactor=[0]*n
  medianfactor=[0.0]*n
  for ni in range(n):
    meanfactor=[ v[ni]/meanval[k] for (k,v) in ctable.items() if k in meanval]
    #print(str(sorted(meanfactor)))
    xfactor=sorted(meanfactor)[len(meanfactor)//2] # corrected 
    if xfactor>0.0:
      medianfactor[ni]=1.0/xfactor
      #logging.debug('xfactor:'+str(xfactor))
  return medianfactor



def normalizeCounts(ctable,method='median',returnfactor=False,reversefactor=False,controlsgfile=None):
  """
  Central function for normalizing read counts
  Parameters:
  --------------
  ctable
    A dictionary of read counts: {sgrnaID:[count0,count1,...]}
  method
    Normalization methods: none,total,median,control
  returnfactor
    Whether to normalize read counts, or just return factors
  reversefactor
    Whether the factor should be reversed (1/factor)
  controlsgfile
    A file name containing control sgRNAs
 
  Return value: 
  --------------
  {sgRNA:[read counts]} if returnfactor == False, or [size_factor] if returnfactor == True
  
  By default, for higher read depths, the size factor is <1. If reversefactor is set, the factor is >1 (or 1/factor) 
  """
  # sums
  if len(ctable)==0:
    return ctable.copy()
  n=len(ctable[list(ctable.keys())[0]]) # samples
  m=len(ctable) # sgRNAs
  # calculate the total factor
  if method=='control':
    # use control sgRNAs 
    if controlsgfile == None:
      logging.error('Error: a list of control sgRNAs should be specified when using control normalization.')
      sys.exit(-1)
    controlsglist=[line.strip() for line in open(controlsgfile)]
    logging.info('Loaded '+str(len(controlsglist))+' control sgRNAs from '+controlsgfile)
    # get the tables for control sgRNAs
    ctable_nm={k:v for (k,v) in ctable.items() if k in controlsglist}
    logging.info('control sgRNAs for normalization:'+str(len(ctable_nm)))
    if len(ctable_nm) ==0:
      logging.error('Error: cannot find records of control sgRNAs in the read count table.')
      sys.exit(-1) 
    method='median' # force to use median of controls
  else:
    ctable_nm=ctable
  samplefactor=mageckcount_gettotalnormfactor(ctable_nm)
  logging.debug('Initial (total) size factor: '+' '.join([str(x) for x in samplefactor]))
  if method=='median':
    # calculate the medianfactor
    medianfactor=mageckcount_getmediannormfactor(ctable_nm)
    usetotalnorm=False
    for ni in range(n):
      if medianfactor[ni]==0.0:
        logging.warning('Sample '+str(ni)+' has zero median count, so median normalization is not possible. Switch to total read count normalization.')
        usetotalnorm=True
    if usetotalnorm == False:
      samplefactor=medianfactor
      logging.debug('Median factor: '+' '.join([str(x) for x in samplefactor]))
    # end median normalization
  elif method=='none': 
    # no normalization
    samplefactor=[1.0]*n
  logging.info('Final size factor: '+' '.join([str(x) for x in samplefactor]))
  
  #
  if returnfactor:
    if reversefactor:
      return [1.0/x for x in samplefactor]
    else:
      return samplefactor
  # normalize the table
  ntable={ k: [ samplefactor[i]*v[i] for i in range(n)] for (k,v) in ctable.items()}
  return ntable
  

def mageckcount_revcomp(x):
  '''
  Reverse complement
  '''
  return x.translate(string.maketrans("ACGT","TGCA"))[::-1]

def mageckcount_gini(x):
  '''
  Return the Gini index of an array
  Calculation is based on http://en.wikipedia.org/wiki/Gini_coefficient
  '''
  xs=sorted(x)
  n=len(xs)
  gssum=sum([ (i+1.0)*xs[i] for i in range(n)])
  ysum=sum(xs)
  if ysum==0.0:
    ysum=1.0
  gs=1.0-2.0*(n-gssum/ysum)/(n-1)
  return gs

def mageckcount_processonefile(filename,args,ctab,genedict,datastat):
  '''
  Go through one fastq file
  Parameters
  ----------
  filename
    Fastq filename to be sequence
  args
    Arguments
  ctab
    A dictionary of sgRNA sequence and count
  genedict
    {sequence:(sgRNA_id,gene_id)} dictionary
  datastat
    Statistics of datasets ({key:value})

  Return value
  -----------
  datastat
    a dictionary structure of statistics
  '''
  # ctab={}
  nline=0
  logging.info('Parsing file '+filename+'...')
  nreadcount=0
  # checking possible sgRNA length
  lengthpool={}
  for k in genedict.keys():
    if len(k) not in lengthpool:
      lengthpool[len(k)]=0
    lengthpool[len(k)]+=1
  lengthpoolkeys=sorted(lengthpool.keys())
  logging.info('Possible gRNA lengths:'+','.join([str(t) for t in lengthpoolkeys]))
  if filename.upper().endswith('.GZ'):
    import gzip
    openobj=gzip.open(filename,'rt')
  else:
    openobj=open(filename)
  for line in openobj:
    # line=line.encode('utf-8')
    nline=nline+1
    if nline%1000000==1:
      logging.info('Processing '+str(round(nline/1000000))+ 'M lines..')
    if nline%4 == 2:
      nreadcount+=1
      fseq=line.strip()
      if args.trim_5 >0:
        fseq=fseq[args.trim_5:]
      # check length
      # for l in lengthpool.keys():
      if len(genedict)==0:
        if len(fseq)<args.sgrna_len:
          continue
        fseq=fseq[:args.sgrna_len]
        if fseq.count('N')>0 and args.count_n==False:
          continue
        if fseq not in ctab:
          ctab[fseq]=0
        ctab[fseq]=ctab[fseq]+1
      else:
        findrecord=False
        for l in lengthpoolkeys: # iterate all possible lengths
          testl=l
          if len(fseq)<testl:
            continue
          fseqc=fseq[:testl]
          if fseqc.count('N')>0 and args.count_n==False:
            continue
          if fseqc not in genedict:
            continue
          else:
            if fseqc not in ctab:
              ctab[fseqc]=0
            ctab[fseqc]=ctab[fseqc]+1
            findrecord=True
            break
        # save unmapped file
        if args.unmapped_to_file and findrecord==False: 
          if len(fseq)<args.sgrna_len:
            continue
          fseqc=fseq[:args.sgrna_len]
          if fseqc.count('N')>0 and args.count_n==False:
            continue
          if fseqc not in ctab:
            ctab[fseqc]=0
          ctab[fseqc]=ctab[fseqc]+1
  # 
  openobj.close()
  # calculate statistics
  datastat['reads']=nreadcount
  # check if a library is provided
  if len(genedict)==0:
    datastat['mappedreads']=0
    datastat['totalsgrnas']=0
    datastat['zerosgrnas']=0
    datastat['giniindex']=1
  else:
    nmapped=0
    nrdcnt=[]
    for (k,v) in ctab.items():
      if k in genedict:
        nmapped+=v
        nrdcnt+=[math.log(v+1.0)]
    nzerosg=0
    for (k,v) in genedict.items():
      if k not in ctab:
        nzerosg+=1
        nrdcnt+=[math.log(0.0+1.0)]
    logging.info('mapped:'+str(nmapped))
    datastat['mappedreads']=nmapped
    datastat['totalsgrnas']=len(genedict);
    datastat['zerosgrnas']=nzerosg
    datastat['giniindex']=mageckcount_gini(nrdcnt)
  #return ctab
  return 0

def mageckcount_mergedict(dict0,dict1):
  '''
  Merge all items in dict1 to dict0.
  '''
  nsample=0
  if len(dict0)>0:
    nsample=len(dict0[dict0.keys()[0]])
  for (k,v) in dict0.items():
    if k in dict1:
      v+=[dict1[k]]
    else:
      v+=[0]
  for (k,v) in dict1.items():
    if k not in dict0:
      if nsample>0:
        dict0[k]=[0]*nsample
      else:
        dict0[k]=[]
      dict0[k]+=[v]
  # return dict0

def mageckcount_printdict(dict0,args,ofile,ounmappedfile,sgdict,datastat,sep='\t'):
  '''
  Write the table count to file
  '''
  allfastq=args.fastq
  nsample=len(allfastq)
  slabel=[datastat[f.split(',')[0]]['label'] for f in allfastq]
  # print header
  print('sgRNA'+sep+'Gene'+sep+sep.join(slabel),file=ofile)
  # print items
  if len(sgdict)==0:
    for (k,v) in dict0.items():
      print(k+sep+'None'+sep+sep.join([str(x) for x in v]),file=ofile)
  else:
    for (k,v) in dict0.items():
      if k not in sgdict: # only print those in the genedict
        if ounmappedfile != None:
          print(sep.join([k,k])+sep+sep.join([str(x) for x in v]),file=ounmappedfile)
        continue
      sx=sgdict[k]
      print(sep.join([sx[0],sx[1]])+sep+sep.join([str(x) for x in v]),file=ofile)
    # print the remaining counts, fill with 0
    for (k,v) in sgdict.items():
      if k not in dict0:
        print(sep.join([v[0],v[1]])+sep+sep.join(["0"]*nsample),file=ofile)

def mageck_printdict(dict0,args,sgdict,sampledict,sampleids):
  """Write the normalized read counts to file
  
  Parameters
  ----------
  dict0 : dict
    a {sgRNA: [read counts]} structure
  args : class
    a argparse class
  sgdict: dict
    a {sgrna:gene} dictionary
  sampledict: dict
    a {sample name: index} dict
  sampleids: list
    a list of sample index. Should include control+treatment
  
  """
  # print header
  # print items
  dfmt="{:.5g}"
  ofile=open(args.output_prefix+'.normalized.txt','w')
  # headers
  mapres_list=['']*len(sampledict)
  for (k,v) in sampledict.items():
    mapres_list[v]=k
  if len(sampledict)>0:
    cntheader=[mapres_list[x] for x in sampleids]
  else:
    cntheader=None
  logging.info('Writing normalized read counts to '+args.output_prefix+'.normalized.txt')
  if cntheader !=None:
    print('sgRNA\tGene\t'+'\t'.join(cntheader),file=ofile)
  if len(sgdict)==0:
    for (k,v) in dict0.items():
      print(k+'\t'+'None'+'\t'+'\t'.join([str(x) for x in v]),file=ofile)
  else:
    for (k,v) in dict0.items():
      if k not in sgdict: # only print those in the genedict
        logging.warning(k+' not in the sgRNA list')
        continue
      print('\t'.join([k,sgdict[k]])+'\t'+'\t'.join([str(x) for x in v]),file=ofile)
    # print the remaining counts, fill with 0
  ofile.close()




def mageckcount_checklists(args):
  """
  Read sgRNA library file
  Including sgRNAs and associated sequences and lists in csv or txt file
  format: sgRNAid  seq  geneid
  """
  genedict={}
  hascsv=False
  if args.list_seq.upper().endswith('CSV'):
    hascsv=True
  n=0
  seqdict={}
  for line in open(args.list_seq):
    if hascsv:
      field=line.strip().split(',')
    else:
      field=line.strip().split()
    n+=1
    if field[0] in genedict:
      logging.warning('Duplicated sgRNA label '+field[0]+' in line '+str(n)+'. Skip this record.')
      continue
    if len(field)<3:
      logging.warning('Not enough field in line '+str(n)+'. Skip this record.')
      continue
    sgrnaseq=field[1].upper()
    if n==1:
      import re
      if re.search('[^ATCG]',sgrnaseq) is not None:
        logging.info('Header line of the library file detected; skip the first line ...')
        continue
    if hasattr(args,'reverse_complement') and args.reverse_complement:
      sgrnaseq=mageckcount_revcomp(sgrnaseq)
    if sgrnaseq in seqdict:
      logging.warning('Duplicated sgRNA sequence '+field[1]+' in line '+str(n)+'. Skip this record.')
      continue
    genedict[field[0]]=(sgrnaseq,field[2])
  logging.info('Loading '+str(len(genedict))+' predefined sgRNAs.')
  return genedict

def mageckcount_printstat(args,datastat):
  '''
  Write data statistics to PDF file 
  '''
  for (k,v) in datastat.items():
    logging.info('Summary of file '+k+':')
    for (v1,v2) in v.items():
      logging.info(str(v1)+'\t'+str(v2))
  # write to table
  crv=VisualRCount()
  crv.setPrefix(args.output_prefix)
  for (fq, fqstat) in datastat.items():
    crv.fastqfile+=[fq]
    if 'label' in fqstat:
      crv.fastqlabels+=[fqstat['label']]
    else:
      crv.fastqlabels+=['NA']
    if 'reads' in fqstat:
      crv.reads+=[fqstat['reads']]
    else:
      crv.reads+=[0]
    if 'mappedreads' in fqstat:
      crv.mappedreads+=[fqstat['mappedreads']]
    else:
      crv.mappedreads+=[0]
    if 'totalsgrnas' in fqstat:
      crv.totalsgrnas+=[fqstat['totalsgrnas']]
    else:
      crv.totalsgrnas+=[0]
    if 'zerosgrnas' in fqstat:
      crv.zerocounts+=[fqstat['zerosgrnas']]
    else:
      crv.zerocounts+=[0]
    if 'giniindex' in fqstat:
      crv.gini+=[fqstat['giniindex']]
    else:
      crv.gini+=[0.0]
  #
  crv.startRTemplate()
  crv.writeCountSummary()
  # write to TXT file
  crv.writeCountSummaryToTxt(args.output_prefix+'.countsummary.txt')
  # write to PDF file
  outcsvfile=args.output_prefix+'.count_normalized.txt'
  crv.insertReadCountBoxPlot(os.path.basename(outcsvfile))
  if len(args.fastq) > 1:
    crv.insertPCAPlot(os.path.basename(outcsvfile))
    crv.insertClusteringPlot(os.path.basename(outcsvfile))
  crv.closeRTemplate()
  if hasattr(args,"pdf_report") and args.pdf_report:
    if hasattr(args,"keep_tmp") :
      crv.generatePDF(keeptmp=args.keep_tmp)
    else:
      crv.generatePDF()

def mageckcount_main(args):
  """
  Main entry for mageck count module
  """
  # check arguments
  mageckcount_checkargs(args)
  # check the listed files
  # listfq=args.fastq.split(',')
  listfq=[[z for z in x.split(',')] for x in args.fastq]
  nsample=len(listfq)
  datastat={}
  # check labels
  alllabel=args.sample_label
  if alllabel=='':
    slabel=['sample'+str(x) for x in range(1,nsample+1)]
  else:
    slabel=alllabel.split(',')
  for i in range(nsample):
    for fi in listfq[i]:
      datastat[fi]={}
      datastat[fi]['label']=slabel[i]
  # process gene dicts
  genedict={}
  if args.list_seq is not None:
    genedict=mageckcount_checklists(args)
  # save sgRNA ID and gene name
  sgdict={} #
  for (k,v) in genedict.items():
    sgdict[v[0]]=(k,v[1])
  alldict={}
  # go through the fastq files
  for filenamelist in listfq:
    dict0={}
    for filename in filenamelist: # technical replicates; should be merged together
      dict00={}
      mageckcount_processonefile(filename,args,dict00,sgdict,datastat[filename])
      for (k,v) in dict00.items():
        if k not in dict0:
          dict0[k]=0
        dict0[k]+=v
    mageckcount_mergedict(alldict,dict0)
  # write to file
  ofilel=open(args.output_prefix+'.count.txt','w')
  if hasattr(args,'unmapped_to_file') and args.unmapped_to_file:
    ounmappedfilel=open(args.output_prefix+'.unmapped.txt','w')
  else:
    ounmappedfilel=None
  mageckcount_printdict(alldict,args,ofilel,ounmappedfilel,sgdict,datastat)
  ofilel.close()
  if hasattr(args,'unmapped_to_file') and args.unmapped_to_file:
    ounmappedfilel.close()
  # write the median normalized read counts to csv file
  ofilel=open(args.output_prefix+'.count_normalized.txt','w')
  if len(sgdict)>0:
    allmappeddict={k:v for (k,v) in alldict.items() if k in sgdict} # only keep those with known sgRNAs
  else:
    allmappeddict=alldict
  if hasattr(args,"norm_method"):
    normmethod=args.norm_method
  else:
    normmethod="median"
  if hasattr(args,"control_sgrna"):
    ctrlsg=args.control_sgrna
  else:
    ctrlsg=None
  medalldict=normalizeCounts(allmappeddict,method=normmethod,controlsgfile=ctrlsg)
  mageckcount_printdict(medalldict,args,ofilel,None,sgdict,datastat,sep='\t')
  ofilel.close()
  # print statistics
  mageckcount_printstat(args,datastat)
  return 0


def getcounttablefromfile(filename):
  """
  read count table from file
  Returns:
  ---------------
  x: dict
    {sgrna:[read counts]} 
  y: dict
    {sgrna:gene}
  z: dict
    z={sample_id:index}
  """
  gtab={}
  mapptab={}
  sampleids={}
  nline=0
  nfield=-1
  # if it is CSV file
  hascsv=False
  if filename.upper().endswith('.CSV'):
    hascsv=True
  logging.info('Loading count table from '+filename+' ')
  for line in open(filename):
    nline+=1
    if nline % 100000 == 1:
      logging.info('Processing '+str(nline)+' lines..')
    try:
      if hascsv==False:
        field=line.strip().split()
      else:
        field=line.strip().split(',')
      if len(field)<3:
        logging.warning('Line '+str(nline)+' of the read count table has fewer than 3 columns. Skip this line ...')
      sgid=field[0]
      geneid=field[1]
      # check if duplicate sgRNA IDs are detected
      if sgid in gtab:
        logging.warning('Duplicated sgRNA IDs: '+sgid+' in line '+str(nline)+'. Skip this record.')
        continue
      mapptab[sgid]=geneid
      sgrecs=[float(x) for x in field[2:]]
      # check the number of fields
      if nfield!=-1 and len(sgrecs)!=nfield:
        logging.error('Error: incorrect number of dimensions in line '+str(nline)+'. Please double-check your read count table file.')
        sys.exit(-1)
      nfield=len(sgrecs)
      gtab[sgid]=sgrecs
    except ValueError:
      if nline!=1:
        logging.warning('Parsing error in line '+str(nline)+'. Skip this line.')
      else:
        logging.debug('Parsing error in line '+str(nline)+' (usually the header line). Skip this line.')
        ids=field[2:]
        for i in range(len(ids)):
          sampleids[ids[i]]=i
      continue
  logging.info('Loaded '+str(len(gtab))+' records.')
  return (gtab,mapptab,sampleids)



if __name__ == '__main__':
  try:
    args=mageckcount_parseargs()
    mageckcount_main(args)
  except KeyboardInterrupt:
    sys.stderr.write("Interrupted.\n")
    sys.exit(0)



