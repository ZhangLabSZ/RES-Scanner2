#!/usr/bin/env python3               
# -*- coding: utf-8 -*-              
#                                    
# Author: Ji Li, liji1@genomics.cn   

import os
import sys
import gzip
import pysam
import getopt
import logging
import subprocess

import shell
from scanner import candidateRES
from replacegeo import replace_genome

# 
SCRIPT = os.path.basename(__file__)
PROGRAM = os.path.splitext(SCRIPT)[0]
PATH = os.path.abspath(__file__)
DIR = os.path.dirname(PATH)

# software path
HISAT , HISAT_BUILD , BWA = '' , '' , ''
SAM2BASE = os.path.join(DIR, 'sam2base.py')
if os.path.isfile(SAM2BASE) == False:
    shell.eprint('['+PROGRAM+'] Error: there is no sam2base.py in the directory '+DIR)
    sys.exit(1)
SCANNER = os.path.join(DIR, 'scanner.py')
if os.path.isfile(SCANNER) == False:
    shell.eprint('['+PROGRAM+'] Error: there is no scanner.py in the directory '+DIR)
    sys.exit(1)

# output directory
OUTDIR = os.getcwd()
GEO = ''

# global variable
BASE = ('A', 'C', 'G', 'T')
BT = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 'N': 'N'}
LINE = 1000000  # RNA reads per iteration

# hyper editing type
TYPE_ALL = ('AC', 'AG', 'AT', 'CA', 'CG', 'CT', 'GA', 'GC', 'GT', 'TA', 'TC', 'TG')
TYPE_SUB = ('AC', 'AG', 'TC', 'TG', 'AT', 'CG')
TYPE_GEO = ('AC', 'AG', 'AT', 'CG')
try:
    TT = {tp: BT[tp[0]]+BT[tp[1]] for tp in TYPE_ALL}
except:
    TT = {'AC': 'TG', 'AG': 'TC', 'AT': 'TA', 'CA': 'GT', 'CG': 'GC', 'CT': 'GA', 'GA': 'CT', 'GC': 'CG', 'GT': 'CA', 'TA': 'AT', 'TC': 'AG', 'TG': 'AC'}
'''
            TYPE      geo(+)  geo(-)  fq1   fq2
alignment:  AC(CA)    A->C    T->G    T->G  A->C
            AG(GA)    A->G    T->C    T->C  A->G
            TC(CT)    T->C    A->G    A->G  T->C
            TG(GT)    T->G    A->C    A->C  T->G
            AT(TA)    A->T    T->A    T->A  A->T 
            CG(GC)    C->G    G->C    G->C  C->G
'''

def extract_unmapped_reads(bamfile, fq1, fq2, mapQ=20, suffix='.fq.gz', filter=True):
    '''
        get unmapped reads name from bamfile, and extract it from fastq file.
    '''
    def _filter(filename):
        def __IsGoodReads(st):
            '''
                Good: N <= 0.1, 0.1 <= [ACGT] <= 0.6
            '''
            n = len(st)
            if st.count('N') / n > 0.1:
                return True
            for b in BASE:
                if (0.1 <= st.count(b) / n <= 0.6) == False:
                    return True
            return False
        s1 = set()
        f = shell.IsGzipFile(filename)
        for idx,line in enumerate(f, 1):
            if idx % 4 == 1:
                ## deal with seq name endswiths /1 or /2
                key = line.strip().split()[0][1:-2] if line.strip().split()[0][-2:] in ('/1', '/2') else line.strip().split()[0][1:]
            elif idx % 4 == 2:
                if __IsGoodReads(line.strip().upper()):
                    s1.add(key)
        f.close()
        return s1
    def _write(infile, outfile):
        f = shell.IsGzipFile(infile)
        for idx,line in enumerate(f, 1):
            if idx % 4 == 1:
                key = line.strip().split()[0][1:-2] if line.strip().split()[0][-2:] in ('/1', '/2') else line.strip().split()[0][1:]
                flag = 1 if key in s else 0
            if flag == 1:
                outfile.write(line)
        f.close()
        return 0
    # check and open bamfile.
    if bamfile.endswith('.bam'):
        opentp = 'rb'
    elif bamfile.endswith('.sam'):
        opentp = 'r'
    else:
        shell.eprint('['+PROGRAM+'] Error: result of alignment should be saved on file which ends with .bam or .sam')
        sys.exit(1)
    s = set() # Global variable only in this function.
    with pysam.AlignmentFile(bamfile, opentp) as f:
        for r in f.fetch(until_eof=True):
            if r.flag & 2 and r.mapping_quality > mapQ:
                continue
            if r.query_name[-2:] in ('/1', '/2'):
                s.add(r.query_name[:-2])
            else:
                s.add(r.query_name)
    # filter low quality reads name.
    if filter:
        for filename in fq1.split(',') + fq2.split(','):
            s = s - _filter(filename)
    # get filtered fastq file.
    with gzip.open('Nmap_1'+suffix, 'wt') as fw1:
        for filename in fq1.split(','):
            _write(filename, fw1)
    with gzip.open('Nmap_2'+suffix, 'wt') as fw2:
        for filename in fq2.split(','):
            _write(filename, fw2)

    return 0

def replace_fq(fq1, fq2, tp):
    for idx,filename in enumerate([fq1, fq2], 1):
        out = gzip.open('fq/'+tp+'_'+str(idx)+'.fq.gz', 'wt')
        with gzip.open(filename, 'rt') as f:
            for n,line in enumerate(f, 1):
                if n % 4 == 2:
                    line = line.replace(BT[tp[0]], BT[tp[1]]) if idx == 1 else line.replace(tp[0], tp[1])
                out.write(line)
        out.close()

    return 0

def alignment(tp, geotp, strand, output='bam/'):
    phred = shell.checkFqQuality('fq/'+tp+'_1.fq.gz')
    if phred == 33:
        p1 = subprocess.Popen([BWA, 'aln', 'geo/'+geotp+'.fa', 'fq/'+tp+'_1.fq.gz', '-n', '0.02', '-o', '0', '-t', '2', '-f', output+tp+'.'+strand+'_1.sai'], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        p2 = subprocess.Popen([BWA, 'aln', 'geo/'+geotp+'.fa', 'fq/'+tp+'_2.fq.gz', '-n', '0.02', '-o', '0', '-t', '2', '-f', output+tp+'.'+strand+'_2.sai'], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    elif phred == 64:
        p1 = subprocess.Popen([BWA, 'aln', 'geo/'+geotp+'.fa', 'fq/'+tp+'_1.fq.gz', '-I', '-n', '0.02', '-o', '0', '-t', '2', '-f', output+tp+'.'+strand+'_1.sai'], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        p2 = subprocess.Popen([BWA, 'aln', 'geo/'+geotp+'.fa', 'fq/'+tp+'_2.fq.gz', '-I', '-n', '0.02', '-o', '0', '-t', '2', '-f', output+tp+'.'+strand+'_2.sai'], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    p1.wait()
    p2.wait()
    if p1.returncode != 0 or p2.returncode != 0:
        shell.eprint('['+PROGRAM+'] Error: Bwa aln wrong')
        sys.exit(1)
    p3 = subprocess.Popen([BWA, 'sampe', '-a', '10000', 'geo/'+geotp+'.fa', output+tp+'.'+strand+'_1.sai', output+tp+'.'+strand+'_2.sai', 'fq/'+tp+'_1.fq.gz', 'fq/'+tp+'_2.fq.gz'], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    infile = pysam.AlignmentFile(p3.stdout, 'r')
    oufile = pysam.AlignmentFile(output+tp+'.'+strand+'.bam', 'wb', template=infile)
    for r in infile:
        oufile.write(r)
    p3.wait()
    for fd in [oufile, infile, p3.stdout]:
        fd.close()
    if p3.returncode != 0:
        shell.eprint('['+PROGRAM+'] Error: Bwa sampe wrong')
        sys.exit(1)

    return 0

def _check(r, d):
    if r.mapping_quality >= 20 and r.cigarstring != '*' and r.flag < 256:
        qual = r.qual
        if r.flag & 128:
            try:
                r.query_sequence = ''.join([BT[x] for x in d[r.query_name][1][::-1]]) if r.flag & 16 else d[r.query_name][1]
                r.qual = qual
            except KeyError:
                shell.eprint('['+PROGRAM+'] Error: Check bam file\n'+r.tostring())
                sys.exit(1)
        else:
            try:
                r.query_sequence = ''.join([BT[x] for x in d[r.query_name][0][::-1]]) if r.flag & 16 else d[r.query_name][0]
                r.qual = qual
            except KeyError:
                shell.eprint('['+PROGRAM+'] Error: Check bam file.\n'+r.tostring())
                sys.exit(1)
        return r
    else:
        return 0
            
def rebuild(f1,f2, p, n, pw, nw):
    d , flag = {} , 1
    for i in range(1, LINE*4+1):
        if i % 4 == 1:
            try:
                line1 , line2 = f1.readline().strip().split()[0][1:] , f2.readline().strip().split()[0][1:]
                key1 = line1[:-2] if line1[-2:] in ('/1', '/2') else line1
                key2 = line2[:-2] if line2[-2:] in ('/1', '/2') else line2
            except IndexError:
                flag = 0
                break
        elif i % 4 == 2:
            if key1 == key2:
                d[key1] = [f1.readline().strip(), f2.readline().strip()]
            else:
                shell.eprint('['+PROGRAM+'] Error: Check fastq file')
                sys.exit(1)
        else:
            f1.readline()
            f2.readline()
    for i in range(1, LINE*2+1):
        try:
            r1 , r2 = _check(next(p), d) , _check(next(n), d)
            if r1:
                pw.write(r1)
            if r2:
                nw.write(r2)
        except StopIteration:
            break
    if flag == 1:
        return rebuild(f1, f2, p, n, pw, nw)
    else:
        return 0

def get_genome_sequence(genomefile):
    dg , f = {} , shell.IsGzipFile(genomefile)
    for line in f:
        if line[0] == '>':
            try:
                dg[key] = st
            except UnboundLocalError:
                pass
            key = line.strip().split()[0][1:]
            st = ''
        else:
            st += line.strip().upper()
    dg[key] = st
    f.close()
    return dg

def recovery(dg, pos, outfile1, neg, outfile2, tp):
    '''
        for fr-firststrand data.
    '''
    def _diff(refer, query, length, trim):
        d , n = {} , 0
        for i in range(trim-1, length-trim):
            if refer[i] != query[i]:
                n += 1
                try:
                    d[(refer[i], query[i])] += 1
                except KeyError:
                    d[(refer[i], query[i])] = 1
        return d,n
    def _quality(seq):
        return sum([ord(i) for i in seq]) / len(seq)
    def _judge(r, tp, dg):
        if len(r.cigar) != 1 or r.cigar[0][0] != 0:
            return -1
        if _quality(r.qual) < 36:
            return -1
        length = len(r.query_sequence)
        t = dg[r.reference_name][r.reference_start: r.reference_start+length]
        if t == r.query_sequence or len(t) != length:
            return -1
        trim = 1 if length < 10 else int(length * 0.1)
        dif , al = _diff(t, r.query_sequence, length, trim)
        if (tp[0], tp[1]) in dif and dif[(tp[0], tp[1])] / al > 0.6 and dif[(tp[0], tp[1])] / (length-2*trim) > 0.05:
            return 1
        elif (tp[1], tp[0]) in dif and dif[(tp[1], tp[0])] / al > 0.6 and dif[(tp[1], tp[0])] / (length-2*trim) > 0.05:
            return 1
        else:
            return 0
    infile = pysam.AlignmentFile(pos, 'rb')
    oufile = pysam.AlignmentFile(outfile1, 'wb', template=infile)
    for r in infile.fetch(until_eof=True):
        # 82 = 2 + 16 + 64; 162 = 2 + 32 + 128
        if (r.flag & 82 == 82 and (r.flag & 32) == 0) or (r.flag & 162 == 162 and (r.flag & 16) == 0):
            if _judge(r, tp, dg) == 1:
                oufile.write(r)
    infile.close()
    oufile.close()
    infile = pysam.AlignmentFile(neg, 'rb')
    oufile = pysam.AlignmentFile(outfile2, 'wb', template=infile)
    for r in infile.fetch(until_eof=True):
        # 98 = 2 + 32 + 64; 146 = 2 + 16 + 128
        if ((r.flag & 98 == 98) and (r.flag & 16) == 0) or ((r.flag & 146 == 146) and (r.flag & 32) == 0):
            if _judge(r, TT[tp], dg) == 1:
                oufile.write(r)
    infile.close()
    oufile.close()

    return 0

def filtersite(fn, outname=''):
    fw = gzip.open('.'.join(fn.split('.')[:-1]) + '.filter.gz', 'wt')
    with gzip.open(fn, 'rt') as f:
        fw.write(f.readline())
        for line in f:
            fd = line.strip().split()
            etp = fd[10].split('->')
            geo , rna = {BASE[i]: int(v) for i,v in enumerate(fd[5].split(','))} , {BASE[i]: int(v) for i,v in enumerate(fd[9].split(','))}
            if fd[2] == '+':
                rna[etp[0]] += 1
                rna[etp[1]] += 1
            elif fd[2] == '-':
                rna[BT[etp[0]]] += 1
                rna[BT[etp[1]]] += 1
            j = 0
            for k in BASE:
                if geo[k] > 0 and rna[k] > 0:
                    j += 1
            if j <= 1:
                fw.write(line)

    return 0

def print_help():
    print('''
Usage:   ires hyper [options] <genomefile> <RNAbamfile> <DNAsbfile> <fq1> <fq2>
         
         <genomefile>  the geonome file, should have the hisat index file in the same directory.
         <RNAbamfile>  RNA bamfile, should be the strand-specific at now.
         <DNAsbfile>   DNA single base file.
         <fq1>         fastq file, read1 in paired-end. 
         <fq2>         fastq file, read2 in paired-end.
         For more than one fastq file, use comma as separator, and gzip file is allowed.

Options: --outdir,-o    DIR    the output directory, '<pwd>/hyper/' is recommended. [./]
         --mapQ         INT    the threshold of mapping quality of bam. [20]    
         --hisat        PATH   the absolute path of hisat2 program or system env. (force)
         --hisat-build  PATH   the absolute path of hisat2-build program or system env. (force)
         --bwa          PATH   the absolute path of bwa alignment program or system env. (force)

         --geo          DIR    the directory of replaced and indexed genome.
         --coverage     INT    the minium depth in sb file. [2]
         --Method       STR    method for detecting SNPs. (Bayesian | Binomial | Frequency). [Bayesian]
         --rate         INT    rate of transitions over transversions of the genome. (force --Method Bayesian) [4]
         --HomoPrior    FLOAT  prior probability of homozygous genomic positions. (force --Method Bayesian) [0.99]
         --Ploid        INT    ploidy level, 1 for monoploid , 2 for diploid, 3 for triploid and so on. [2]
         --Qual-cutoff  INT    quality cutoff for BWA alignment. [30]
         --Phred-DNA    INT    phred base quality for query quality of DNA. (33 | 64). [33]
         --Phred-RNA    INT    phred base quality for query quality of RNA. (33 | 64). [33]

         --help,-h             print help information
''')

def get_config():
    try:
        optlist , args = getopt.getopt(sys.argv[1:], 'ho:m:', ['help', 'outdir=', 'mapQ=', 'hisat=', 'hisat-build=', 'bwa=', 'coverage=', 'Method=', 'rate=', 'HomoPrior=', 'Ploid=', 'Qual-cutoff=', 'Phred-DNA=', 'Phred-RNA=', 'geo='])
    except getopt.GetoptError as e:
        shell.eprint('['+PROGRAM+'] Error: '+str(e))
        sys.exit(2)
    
    if optlist == [] and args == []:
        print_help()
        sys.exit(0)
    
    tobool , config = {'T': True, 'F': False} , {'extract_unmapped_reads': {}, 'candidateRES': {}}
    global OUTDIR , GEO , HISAT , HISAT_BUILD , BWA
    for opt , value in optlist:
        if opt in ('--help', '-h'):
            print_help()
            sys.exit(0)
        elif opt in ('--outdir', '-o'):
            OUTDIR = os.path.abspath(value) + '/'
        elif opt == '--geo':
            if os.path.exists(value) == True:
                GEO = os.path.abspath(value) + '/'
        elif opt in ('--mapQ', 'm'):
            try:
                config['extract_unmapped_reads']['mapQ'] = int(value)
            except ValueError as e:
                shell.eprint('['+PROGRAM+'] Error: --mapQ should be non-negative integer')
                sys.exit(1)
        elif opt == '--hisat':
            if os.path.exists(value):
                HISAT = os.path.abspath(value)
            else:
                shell.eprint('['+PROGRAM+'] Error: hisat2 path does not exist')
                sys.exit(1)
        elif opt == '--hisat-build':
            if os.path.exists(value):
                HISAT_BUILD = os.path.abspath(value)
            else:
                shell.eprint('['+PROGRAM+'] Error: hisat2-build path does not exist')
                sys.exit(1)
        elif opt == '--bwa':
            if os.path.exists(value):
                BWA = os.path.abspath(value)
            else:
                shell.eprint('['+PROGRAM+'] Error: bwa path does not exist')
                sys.exit(1)
        elif opt == '--coverage':
            try:
                config['candidateRES']['Coverage'] = int(value)
            except ValueError as e:
                shell.eprint('['+PROGRAM+'] Error: --coverage should be non-negative integer')
                sys.exit(1)
        elif opt == '--Method':
            config['candidateRES']['Method'] = value
        elif opt == '--rate':
            try:
                config['candidateRES']['Ratio'] = int(value)
            except ValueError as e:
                shell.eprint('['+PROGRAM+'] Error: --rate should be non-negative integer')
                sys.exit(1)
        elif opt == '--HomoPrior':
            try:
                config['candidateRES']['HomoPrior'] = float(value)
            except ValueError as e:
                shell.eprint('['+PROGRAM+'] Error: --HomoPrior should be float')
                sys.exit(1)
        elif opt == '--Ploid':
            try:
                config['candidateRES']['Ploid'] = int(value)
            except ValueError as e:
                shell.eprint('['+PROGRAM+'] Error: --Ploid should be non-negative integer')
                sys.exit(1)
        elif opt == '--Qual-cutoff':
            try:
                config['candidateRES']['Qual_cutoff'] = int(value)
            except ValueError as e:
                shell.eprint('['+PROGRAM+'] Error: --Qual-cutoff should be non-negative integer')
                sys.exit(1)
        elif opt == '--Phred-DNA':
            try:
                config['candidateRES']['Phred_DNA'] = int(value)
            except ValueError as e:
                shell.eprint('['+PROGRAM+'] Error: --Phred-DNA should be 33 or 64')
                sys.exit(1)
            if int(value) != 33 or int(value) != 64:
                shell.eprint('['+PROGRAM+'] Error: --Phred-DNA should be 33 or 64')
                sys.exit(1)
        elif opt == '--Phred-RNA':
            try:
                config['candidateRES']['Phred_RNA'] = int(value)
            except ValueError as e:
                shell.eprint('['+PROGRAM+'] Error: --Phred-RNA should be 33 or 64')
                sys.exit(1)
            if int(value) != 33 or int(value) != 64:
                shell.eprint('['+PROGRAM+'] Error: --Phred-RNA should be 33 or 64')
                sys.exit(1)
        else:
            assert False, 'unhandled option'

    try:
        genomefile , bamfile , DNAsbfile , fq1 , fq2 = args
    except ValueError as e:
        shell.eprint('['+PROGRAM+'] Error: '+str(e))
        sys.exit(1)
    # check input file
    genomefile = os.path.abspath(genomefile)
    bamfile = os.path.abspath(bamfile)
    DNAsbfile = os.path.abspath(DNAsbfile)
    fq1 = ';'.join([os.path.abspath(fn) for fn in fq1.split(';')])
    fq2 = ';'.join([os.path.abspath(fn) for fn in fq2.split(';')])
    for fn in [genomefile, bamfile, DNAsbfile] + fq1.split(';') + fq2.split(';'):
        if os.path.isfile(fn) == False:
            shell.eprint('['+PROGRAM+'] Error: please check the input file')
            sys.exit(1)

    if HISAT != '':
        temp = os.path.join(os.path.dirname(HISAT), 'hisat2-build')
        HISAT_BUILD = temp if os.path.isfile(temp) else ''
    if HISAT == '':
        r = subprocess.getstatusoutput('which hisat2')
        if r[0] == 1:
            shell.eprint('['+PROGRAM+'] Warning: lack hisat2 program')
            sys.exit(1)
        else:
            HISAT = r[1]
    if HISAT_BUILD == '':
        r = subprocess.getstatusoutput('which hisat2-build')
        if r[0] == 1:
            shell.eprint('['+PROGRAM+'] Warning: lack hisat2-build program')
            sys.exit(1)
        else:
            HISAT_BUILD = r[1]
    if BWA == '':
        r = subprocess.getstatusoutput('which bwa')
        if r[0] == 1:
            shell.eprint('['+PROGRAM+'] Warning: lack bwa program')
            sys.exit(1)
        else:
            BWA = r[1]

    if os.path.exists(OUTDIR) == False:
        try:
            os.makedirs(OUTDIR)
        except:
            shell.eprint('['+PROGRAM+'] Error: outdir could not be created, please check')
            sys.exit(1)

    os.chdir(OUTDIR)

    logging.basicConfig(level=logging.INFO,
                        filename=OUTDIR+PROGRAM+'.log',
                        filemode='w',
                        format='%(asctime)s : %(message)s',
                        datefmt='%m/%d/%Y %I:%M:%S %p')

    return genomefile , bamfile , DNAsbfile , fq1 , fq2 , config

def main():
    genomefile , bamfile , DNAsbfile , fq1 , fq2 , config = get_config()
    # check file
    FLAG = 0    # check hisat index
    logging.info('Program Start')
    logging.info('make work diretory done')
    # mkdir in the output directory.
    dirs = ['geo', 'fq', 'bam', 'singlebase', 'table']
    for i in range(len(dirs)):
        if os.path.exists(dirs[i]) == False:
            try:
                os.makedirs(dirs[i])
            except FileExistsError:
                pass
            except OSError as e:
                shell.eprint('['+PROGRAM+'] Error: '+str(e))
                for j in range(i):
                    os.system('rm -rf '+dirs[j])
                sys.exit(1)
    if GEO == '':
        logging.info('do genome base change and bwa index, this step may use large time if your genome is huge.')
        for tp in TYPE_GEO:
            replace_genome(genomefile, tp, 'geo/'+tp+'.fa', 'geo/'+TT[tp]+'.fa', BWA)
        logging.info('done')
    else:
        l = shell.isbwaindex(GEO, TYPE_GEO, TT)
        if l != []:
            try:
                for i in range(len(l)):
                    os.symlink(GEO+l[i], 'geo/'+l[i])
            except FileExistsError:
                pass
            except:
                shell.eprint('['+PROGRAM+'] Error: can not make soft link in geo/')
                for j in range(i):
                    os.system('rm -f '+'geo/'+l[j])
                sys.exit(1)
        else:
            shell.eprint('['+PROGRAM+'] Error: --geo directory do not have all replaced genome index')
            os.system('rm -rf geo/')
            sys.exit(1)
    logging.info('make sub dirs done.')
    # check hisat index
    logging.info('check hisat index file')
    if shell.ishisatindex(os.path.dirname(genomefile)) == 0:
        logging.info('do hisat-build')
        FLAG = 1
        try:
            os.symlink(genomefile, 'geo/'+os.path.basename(genomefile))
        except:
            shell.eprint('['+PROGRAM+'] Error: can not make soft link of genome file')
            os.system('rm -rf geo/')
            sys.exit(1)
        p = subprocess.Popen([HISAT_BUILD, 'geo/'+os.path.basename(genomefile), 'geo/'+os.path.basename(genomefile)], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        p.wait()
        logging.info('done')
        if p.returncode != 0:
            shell.eprint('['+PROGRAM+'] Error: hisat index failed')
            sys.exit(1)
    # extract unmapped reads.
    logging.info('extract unmapped reads from RNA alignment result')
    extract_unmapped_reads(bamfile, fq1, fq2, **config['extract_unmapped_reads'])
    logging.info('done')
    # Hisat realign.
    logging.info('do hisat alignment')
    if FLAG == 0:
        if 'Phred_RNA' in config['candidateRES'].keys() and config['candidateRES']['Phred_RNA'] == 64:
            p = subprocess.Popen([HISAT, '--phred64', '--rna-strandness', 'RF', '-x', genomefile, '-1', 'Nmap_1.fq.gz', '-2', 'Nmap_2.fq.gz'], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        else:
            p = subprocess.Popen([HISAT, '--rna-strandness', 'RF', '-x', genomefile, '-1', 'Nmap_1.fq.gz', '-2', 'Nmap_2.fq.gz'], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    elif FLAG == 1:
        if 'Phred_RNA' in config['candidateRES'].keys() and config['candidateRES']['Phred_RNA'] == 64:
            p = subprocess.Popen([HISAT, '--phred64', '--rna-strandness', 'RF', '-x', 'geo/'+os.path.basename(genomefile), '-1', 'Nmap_1.fq.gz', '-2', 'Nmap_2.fq.gz'], stdout=subprocess.PIPE, stderr=subprocess.PIPE)            
        else:
            p = subprocess.Popen([HISAT, '--rna-strandness', 'RF', '-x', 'geo/'+os.path.basename(genomefile), '-1', 'Nmap_1.fq.gz', '-2', 'Nmap_2.fq.gz'], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    infile = pysam.AlignmentFile(p.stdout, 'r')
    oufile = pysam.AlignmentFile('hisat.bam', 'wb', template=infile)
    for r in infile:
        oufile.write(r)
    infile.close()
    oufile.close()
    p.wait()
    if p.returncode != 0:
        shell.eprint(p.stderr.read().decode())
    logging.info('done')
    # extract unmapped reads from hisat align.
    logging.info('extract unmapped reads from hisat alignment result')
    extract_unmapped_reads('hisat.bam', 'Nmap_1.fq.gz', 'Nmap_2.fq.gz', mapQ=-1, suffix='.clean.fq.gz', filter=False)
    logging.info('done')
    # base masking.
    logging.info('do base change for unmapped reads')
    for tp in TYPE_SUB:
        replace_fq('Nmap_1.clean.fq.gz', 'Nmap_2.clean.fq.gz', tp)
    logging.info('done')
    # then bwa align.
    logging.info('bwa alignment for base-changed unmapped reads and rebuild alignment result')
    for tp in TYPE_SUB:
        alignment(tp, tp, 'pos')
        alignment(tp, TT[tp], 'neg')
        fq1 , fq2 , p , n = gzip.open('Nmap_1.clean.fq.gz', 'rt') , gzip.open('Nmap_2.clean.fq.gz', 'rt') , pysam.AlignmentFile('bam/'+tp+'.pos.bam', 'rb') , pysam.AlignmentFile('bam/'+tp+'.neg.bam', 'rb')
        pw , nw = pysam.AlignmentFile('bam/'+tp+'.pos.rebuild.bam', 'wb', template=p) , pysam.AlignmentFile('bam/'+tp+'.neg.rebuild.bam', 'wb', template=n)
        rebuild(fq1, fq2, p, n, pw, nw)
        for fd in [nw, pw, n, p, fq2, fq1]:
            fd.close()
    logging.info('done')
    # filter low mutation.
    logging.info('do some filter')
    dg = get_genome_sequence(genomefile)
    for tp in TYPE_SUB:
        recovery(dg, 'bam/'+tp+'.pos.rebuild.bam', 'bam/'+tp+'.RNA.pos.bam', 'bam/'+tp+'.neg.rebuild.bam', 'bam/'+tp+'.RNA.neg.bam', tp)
        for strand in ['pos', 'neg']:
            pysam.sort('bam/'+tp+'.RNA.'+strand+'.bam', '-o', 'bam/'+tp+'.RNA.'+strand+'.sort.bam')
            pysam.rmdup('bam/'+tp+'.RNA.'+strand+'.sort.bam', 'bam/'+tp+'.RNA.'+strand+'.sort.rmdup.bam')
            pysam.index('bam/'+tp+'.RNA.'+strand+'.sort.rmdup.bam')
            os.system('rm -f bam/'+tp+'.RNA.'+strand+'.sort.bam')
            os.system('rm -f bam/'+tp+'.RNA.'+strand+'.bam')
            os.system('rm -f bam/'+tp+'.'+strand+'.rebuild.bam')
            os.system('rm -f bam/'+tp+'.'+strand+'_1.sai')
            os.system('rm -f bam/'+tp+'.'+strand+'_2.sai')
            os.system('rm -f bam/'+tp+'.'+strand+'.bam')
    logging.info('done')
    # form single base file.
    logging.info('form single base file')
    for tp in TYPE_SUB:
        p1 = subprocess.Popen(['python', SAM2BASE, '--trim', '10,10', '--prefix', 'singlebase/'+tp+'.pos', genomefile, 'bam/'+tp+'.RNA.pos.sort.rmdup.bam'])
        p2 = subprocess.Popen(['python', SAM2BASE, '--trim', '10,10', '--prefix', 'singlebase/'+tp+'.neg', genomefile, 'bam/'+tp+'.RNA.neg.sort.rmdup.bam'])
        p1.wait()
        p2.wait()
    logging.info('done')
    # do scanner.
    logging.info('do scanner')
    for tp in TYPE_SUB:
        candidateRES('singlebase/'+tp+'.pos.sb.gz', DNAsbfile, genomefile, '+', Output=OUTDIR+'table/', **config['candidateRES'])
        candidateRES('singlebase/'+tp+'.neg.sb.gz', DNAsbfile, genomefile, '-', Output=OUTDIR+'table/', **config['candidateRES'])
    logging.info('done')
    # filter sites.
    logging.info('do editing site filter')
    for tp in TYPE_SUB:
        filtersite('table/'+tp+'.pos.sb.gz.homo.gz')
        filtersite('table/'+tp+'.neg.sb.gz.homo.gz')
    logging.info('done')
    logging.info('All things had been done. Have a good day!')

    return 0




if __name__ == '__main__':
    main()
