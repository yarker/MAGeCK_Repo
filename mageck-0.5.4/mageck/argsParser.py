"""MAGeCK argument parser
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
import logging

from mageck.mlemageck import *

def crisprseq_parseargs():
  """
  Parsing mageck arguments.
  """
  parser=argparse.ArgumentParser(description='mageck: performs sgRNA, gene and pathway analysis on CRISPR-Cas9 screening data.')
  # definition of sub commands
  subparser=parser.add_subparsers(help='commands to run mageck',dest='subcmd')
  
  parser.add_argument('-v', '--version',action='version',version='%(prog)s 0.5.4')
  
  # run: execute the whole program
  # disabled since 0.5.4
  # subp_run=subparser.add_parser('run',help='Run the whole MAGeCK program from fastq files, with limited options.')
  # subp_run.add_argument('-l','--list-seq',help='A file containing the list of sgRNA names, their sequences and associated genes. Support file format: csv and txt.')
  # subp_run.add_argument('-n','--output-prefix',default='sample1',help='The prefix of the output file(s). Default sample1.')
  # subp_run.add_argument('--sample-label',default='',help='Sample labels, separated by comma. Must be equal to the number of samples provided (in --fastq option). Default "sample1,sample2,...".')
  # subp_run.add_argument('--trim-5',type=int,default=0,help='Length of trimming the 5\' of the reads. Default 0')
  # subp_run.add_argument('--sgrna-len',type=int,default=20,help='Length of the sgRNA. Default 20  (after v0.5.3, the program will automatically determine the sgRNA length from library file)')
  # subp_run.add_argument('--count-n',action='store_true',help='Count sgRNAs with Ns. By default, sgRNAs containing N will be discarded.')
  # #
  # subp_run.add_argument('--gene-test-fdr-threshold',type=float,default=0.25,help='FDR threshold for gene test, default 0.25.')
  # subp_run.add_argument('-t','--treatment-id',required=True,action='append',help='Sample label or sample index (integer, 0 as the first sample)  as treatment experiments, separated by comma (,). If sample label is provided, the labels must match the labels in the first line of the count table; for example, "HL60.final,KBM7.final". For sample index, "0,2" means the 1st and 3rd samples are treatment experiments.')
  # subp_run.add_argument('-c','--control-id',action='append',help='Sample label or sample index in the count table as control experiments, separated by comma (,). Default is all the samples not specified in treatment experiments.')
  # subp_run.add_argument('--adjust-method',choices=['fdr','holm','pounds'],default='fdr',help='Method for sgrna-level p-value adjustment, including false discovery rate (fdr), holm\'s method (holm), or pounds\'s method (pounds).')
  # subp_run.add_argument('--variance-from-all-samples',action='store_true',help='Estimate the variance from all samples, instead of from only control samples. Use this option only if you believe there are relatively few essential sgRNAs or genes between control and treatment samples.')
  # subp_run.add_argument('--keep-tmp',action='store_true',help='Keep intermediate files.')
  # subp_run.add_argument('--fastq',nargs='+',help='Sample fastq files, separated by space; use comma (,) to indicate technical replicates of the same sample. For example, "--fastq sample1_replicate1.fastq,sample1_replicate2.fastq sample2_replicate1.fastq,sample2_replicate2.fastq" indicates two samples with 2 technical replicates for each sample.')
  # subp_run.add_argument('--pdf-report',action='store_true',help='Generate pdf report of the analysis.')
  
  # countonly: only collect counts
  subp_count=subparser.add_parser('count',help='Collecting read counts from fastq files.')
  subp_count.add_argument('-l','--list-seq',required=True,help='A file containing the list of sgRNA names, their sequences and associated genes. Support file format: csv and txt. Provide an empty file for collecting all possible sgRNA counts.')
  subp_count.add_argument('--sample-label',default='',help='Sample labels, separated by comma (,). Must be equal to the number of samples provided (in --fastq option). Default "sample1,sample2,...".')
  subp_count.add_argument('-n','--output-prefix',default='sample1',help='The prefix of the output file(s). Default sample1.')
  subp_count.add_argument('--trim-5',type=int,default=0,help='Length of trimming the 5\' of the reads. Default 0')
  subp_count.add_argument('--sgrna-len',type=int,default=20,help='Length of the sgRNA. Default 20. ATTENTION: after v0.5.3, the program will automatically determine the sgRNA length from library file; so only use this if you turn on the --unmapped-to-file option.')
  subp_count.add_argument('--count-n',action='store_true',help='Count sgRNAs with Ns. By default, sgRNAs containing N will be discarded.')
  subp_count.add_argument('--unmapped-to-file',action='store_true',help='Save unmapped reads to file, with sgRNA lengths specified by --sgrna-len option.')
  subp_count.add_argument('--reverse-complement',action='store_true',help='Reverse complement the sequences in library for read mapping.')
  subp_count.add_argument('--fastq',nargs='+',help='Sample fastq files (or fastq.gz files), separated by space; use comma (,) to indicate technical replicates of the same sample. For example, "--fastq sample1_replicate1.fastq,sample1_replicate2.fastq sample2_replicate1.fastq,sample2_replicate2.fastq" indicates two samples with 2 technical replicates for each sample.')
  subp_count.add_argument('--pdf-report',action='store_true',help='Generate pdf report of the fastq files.')
  subp_count.add_argument('--keep-tmp',action='store_true',help='Keep intermediate files.')
  subp_count.add_argument('--norm-method',choices=['none','median','total','control'],default='median',help='Method for normalization, including "none" (no normalization), "median" (median normalization, default), "total" (normalization by total read counts), "control" (normalization by control sgRNAs specified by the --control-sgrna option).')
  subp_count.add_argument('--control-sgrna',help='A list of control sgRNAs for normalization and for generating the null distribution of RRA.')
  
  # stat test: only do the statistical test
  subn_stattest=subparser.add_parser('test',help='Perform statistical test from a given count table (generated by count command).')
  # required parameters
  reqgroup=subn_stattest.add_argument_group(title='Required arguments',description='')
  reqgroup.add_argument('-k','--count-table',required=True,help='Provide a tab-separated count table instead of sam files. Each line in the table should include sgRNA name (1st column), gene name (2nd column) and read counts in each sample.')
  # this parameter is depreciated
  # subn_stattest.add_argument('--gene-test',help='Perform rank association analysis. A tab-separated, sgRNA to gene mapping file is required. Each line in the file should include two columns, the sgRNA name and the gene name.')
  reqgroup.add_argument('-t','--treatment-id',required=True,action='append',help='Sample label or sample index (0 as the first sample) in the count table as treatment experiments, separated by comma (,). If sample label is provided, the labels must match the labels in the first line of the count table; for example, "HL60.final,KBM7.final". For sample index, "0,2" means the 1st and 3rd samples are treatment experiments.')
  ## Optional general arguments
  gengroup=subn_stattest.add_argument_group(title='Optional general arguments',description='')
  gengroup.add_argument('-c','--control-id',action='append',help='Sample label or sample index in the count table as control experiments, separated by comma (,). Default is all the samples not specified in treatment experiments.')
  gengroup.add_argument('--norm-method',choices=['none','median','total','control'],default='median',help='Method for normalization, including "none" (no normalization), "median" (median normalization, default), "total" (normalization by total read counts), "control" (normalization by control sgRNAs specified by the --control-sgrna option).')
  gengroup.add_argument('--gene-test-fdr-threshold',type=float,default=0.25,help='FDR threshold for gene test, default 0.25.')
  gengroup.add_argument('--adjust-method',choices=['fdr','holm','pounds'],default='fdr',help='Method for sgrna-level p-value adjustment, including false discovery rate (fdr), holm\'s method (holm), or pounds\'s method (pounds).')
  gengroup.add_argument('--variance-from-all-samples',action='store_true',help='Estimate the variance from all samples, instead of from only control samples. Use this option only if you believe there are relatively few essential sgRNAs or genes between control and treatment samples.')
  gengroup.add_argument('--sort-criteria',choices=['neg','pos'],default='neg',help='Sorting criteria, either by negative selection (neg) or positive selection (pos). Default negative selection.')
  gengroup.add_argument('--remove-zero',choices=['none','control','treatment','both'],default='none',help='Whether to remove zero-count sgRNAs in control and/or treatment experiments. Default: none (do not remove those zero-count sgRNAs).')
  gengroup.add_argument('--pdf-report',action='store_true',help='Generate pdf report of the analysis.')
  ## Optional IO arguments
  iogroup=subn_stattest.add_argument_group(title='Optional arguments for input and output',description='')
  iogroup.add_argument('-n','--output-prefix',default='sample1',help='The prefix of the output file(s). Default sample1.')
  #parser.add_argument('--count-control-index',help='If -k/--counts option is given, this parameter defines the control experiments in the table.'); 
  iogroup.add_argument('--control-sgrna',help='A list of control sgRNAs for normalization and for generating the null distribution of RRA.')
  iogroup.add_argument('--normcounts-to-file',action='store_true',help='Write normalized read counts to file ([output-prefix].normalized.txt).')
  iogroup.add_argument('--keep-tmp',action='store_true',help='Keep intermediate files.')
  
  # pathway test
  subw_pathway=subparser.add_parser('pathway',help='Perform significant pathway analysis from gene rankings generated by the test command.')
  subw_pathway.add_argument('--gene-ranking',required=True,help='The gene summary file (containing both positive and negative selection tests) generated by the gene test step. Pathway enrichment will be performed in both directions.')
  subw_pathway.add_argument('--single-ranking',action='store_true',help='The provided file is a (single) gene ranking file, either positive or negative selection. Only one enrichment comparison will be performed.')
  # subw_pathway.add_argument('--gene-ranking-2',help='An optional gene ranking file of opposite direction, generated by the gene test step.')
  subw_pathway.add_argument('--gmt-file',required=True,help='The pathway file in GMT format.')
  subw_pathway.add_argument('-n','--output-prefix',default='sample1',help='The prefix of the output file(s). Default sample1.')
  subw_pathway.add_argument('--sort-criteria',choices=['neg','pos'],default='neg',help='Sorting criteria, either by negative selection (neg) or positive selection (pos). Default negative selection.')
  subw_pathway.add_argument('--keep-tmp',action='store_true',help='Keep intermediate files.')
  subw_pathway.add_argument('--ranking-column',default='2',help='Column number or label in gene summary file for gene ranking; can be either an integer of column number, or a string of column label. Default "2" (the 3rd column).')
  subw_pathway.add_argument('--ranking-column-2',default='7',help='Column number or label in gene summary file for gene ranking; can be either an integer of column number, or a string of column label. This option is used to determine the column for positive selections and is disabled if --single-ranking is specified. Default "5" (the 6th column).')
  
  # plot
  subp_plot=subparser.add_parser('plot',help='Generating graphics for selected genes.')
  subp_plot.add_argument('-k','--count-table',required=True,help='Provide a tab-separated count table instead of sam files. Each line in the table should include sgRNA name (1st column), gene name (2nd column) and read counts in each sample.')
  subp_plot.add_argument('-g','--gene-summary',required=True,help='The gene summary file generated by the test command.')
  subp_plot.add_argument('--genes',help='A list of genes to be plotted, separated by comma. Default: none.')
  subp_plot.add_argument('-s','--samples',help='A list of samples to be plotted, separated by comma. Default: using all samples in the count table.')
  subp_plot.add_argument('-n','--output-prefix',default='sample1',help='The prefix of the output file(s). Default sample1.')
  subp_plot.add_argument('--norm-method',choices=['none','median','total','control'],default='median',help='Method for normalization, including "none" (no normalization), "median" (median normalization, default), "total" (normalization by total read counts), "control" (normalization by control sgRNAs specified by the --control-sgrna option).')
  subp_plot.add_argument('--control-sgrna',help='A list of control sgRNAs.')
  subp_plot.add_argument('--keep-tmp',action='store_true',help='Keep intermediate files.')
  
  subm_mle=subparser.add_parser('mle',help='Perform MLE estimation of gene essentiality.')
  # required parameters
  reqgroup=subm_mle.add_argument_group(title='Required arguments',description='')
  reqgroup.add_argument('-k','--count-table',required=True,help='Provide a tab-separated count table. Each line in the table should include sgRNA name (1st column), target gene (2nd column) and read counts in each sample.')
  reqgroup.add_argument('-d','--design-matrix',required=True,help='Provide a design matrix, either a file name or a quoted string of the design matrix. For example, "1,1;1,0". The row of the design matrix must match the order of the samples in the count table (if --include-samples is not specified), or the order of the samples by the --include-samples option.')
  # optional arguments
  #opgroup=subm_mle.add_argument_group(title='Optional arguments',description='Optional arguments')
  ## input and output
  iogroup=subm_mle.add_argument_group(title='Optional arguments for input and output',description='')
  iogroup.add_argument('-n','--output-prefix',default='sample1',help='The prefix of the output file(s). Default sample1.')
  iogroup.add_argument('-i', '--include-samples', help='Specify the sample labels if the design matrix is not given by file in the --design-matrix option. Sample labels are separated by ",", and must match the labels in the count table.')
  iogroup.add_argument('-b', '--beta-labels', help='Specify the labels of the variables (i.e., beta), if the design matrix is not given by file in the --design-matrix option. Should be separated by ",", and the number of labels must equal to (# columns of design matrix), including baseline labels. Default value: "bata_0,beta_1,beta_2,...".')
  iogroup.add_argument('--control-sgrna',help='A list of control sgRNAs.')
  ## General options
  mlegroup=subm_mle.add_argument_group(title='Optional arguments for MLE module',description='')
  mlegroup.add_argument('--debug',action='store_true',help='Debug mode to output detailed information of the running.')
  mlegroup.add_argument('--debug-gene',help='Debug mode to only run one gene with specified ID.')
  mlegroup.add_argument('--norm-method',choices=['none','median','total','control'],default='median',help='Method for normalization, including "none" (no normalization), "median" (median normalization, default), "total" (normalization by total read counts), "control" (normalization by control sgRNAs specified by the --control-sgrna option).')
  mlegroup.add_argument('--genes-varmodeling',type=int,default=0,help='The number of genes for mean-variance modeling. Default 0.')
  mlegroup.add_argument('--permutation-round',type=int,default=10,help='The rounds for permutation (interger). The permutation time is (# genes)*x for x rounds of permutation. Suggested value: 100 (may take longer time). Default 10.')
  mlegroup.add_argument('--remove-outliers', action='store_true', help='Try to remove outliers. Turning this option on will slow the algorithm.')
  mlegroup.add_argument('--threads', type=int,default=1, help='Using multiple threads to run the algorithm. Default using only 1 thread.')
  mlegroup.add_argument('--adjust-method',choices=['fdr','holm','pounds'],default='fdr',help='Method for sgrna-level p-value adjustment, including false discovery rate (fdr), holm\'s method (holm), or pounds\'s method (pounds).')
  ## EM options
  emgroup=subm_mle.add_argument_group(title='Optional arguments for the EM iteration',description='')
  emgroup.add_argument('--sgrna-efficiency',help='An optional file of sgRNA efficiency prediction. The efficiency prediction will be used as an initial guess of the probability an sgRNA is efficient. Must contain at least two columns, one containing sgRNA ID, the other containing sgRNA efficiency prediction.')
  emgroup.add_argument('--sgrna-eff-name-column',type=int,default=0,help='The sgRNA ID column in sgRNA efficiency prediction file (specified by the --sgrna-efficiency option). Default is 0 (the first column).')
  emgroup.add_argument('--sgrna-eff-score-column',type=int,default=1,help='The sgRNA efficiency prediction column in sgRNA efficiency prediction file (specified by the --sgrna-efficiency option). Default is 1 (the second column).')
  emgroup.add_argument('--update-efficiency',action='store_true',help='Iteratively update sgRNA efficiency during EM iteration.')
  ## Bayes option
  bayesgroup=subm_mle.add_argument_group(title='Optional arguments for the Bayes estimation of gene essentiality (experimental)',description='')
  bayesgroup.add_argument("--bayes",action='store_true',help="Use the experimental Bayes module to estimate gene essentiality")
  bayesgroup.add_argument("-p","--PPI-prior",action='store_true',help="Specify whether you want to incorporate PPI as prior")
  bayesgroup.add_argument("-w","--PPI-weighting",type=float,help="The weighting used to calculate PPI prior. If not provided, iterations will be used.",default=None)
  bayesgroup.add_argument("-e","--negative-control",help="The gene name of negative controls. The corresponding sgRNA will be viewed independently.",default=None)

  #
  # subm_mle.add_argument('args', nargs=argparse.REMAINDER)
  
  args=parser.parse_args()
  
  if args.subcmd == None:
    parser.print_help()
    sys.exit(0)
  
  if args.subcmd == 'mle':
    mageckmle_main(parsedargs=args); # ignoring the script path, and the sub command
    sys.exit(0)
  
  # logging status
  if args.subcmd=='pathway':
    logmode="a"
  else:
    logmode="w"
  
  logging.basicConfig(level=10,
    format='%(levelname)-5s @ %(asctime)s: %(message)s ',
    datefmt='%a, %d %b %Y %H:%M:%S',
    # stream=sys.stderr,
    filename=args.output_prefix+'.log',
    filemode=logmode
  )
  console = logging.StreamHandler()
  console.setLevel(logging.INFO)
  # set a format which is simpler for console use
  formatter = logging.Formatter('%(levelname)-5s @ %(asctime)s: %(message)s ','%a, %d %b %Y %H:%M:%S')
  #formatter.formatTime('%a, %d %b %Y %H:%M:%S')
  # tell the handler to use this format
  console.setFormatter(formatter)
  # add the handler to the root logger
  logging.getLogger('').addHandler(console)
  
  # add paramters
  logging.info('Parameters: '+' '.join(sys.argv))
  
  
  return args


