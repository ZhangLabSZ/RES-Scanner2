#!/usr/bin/python3
# -*- coding: utf-8 -*-
#
# Author: Ji Li, liji1@genomics.cn

import os
import sys
import gzip
import getopt
import logging
import subprocess

from time import time,asctime,localtime

D_BASE_COMPLIMENT = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 'N': 'N'}
T_BASE_ORDER = ('A', 'C', 'G', 'T')

def log(fn):
    def decorator(func):
        logging.basicConfig(level=logging.INFO,
                            filename=fn,
                            filemode='w',
                            format='%(asctime)s : %(message)s',
                            datefmt='%m/%d/%Y %I:%M:%S %p')
        def wrapper(*args, **kw):
            return func(*args, **kw)
        return wrapper
    return decorator

def log2(st):
    return asctime(localtime(time())) + ' : ' + st + '\n'

def reads_split(fq1, fq2=None, outdir='./', size='3G'):
    '''
        Split file into sized compressed files.
    '''
    def _filetype(filename):
        '''
            Check the filetype of Input and return the handle of file.
        '''
        file_suffix_1 , file_suffix_2 = filename.split('.')[-2:]
        if file_suffix_2 == 'gz':
            if file_suffix_1 == 'fq' or file_suffix_1 == 'fastq':
                f = gzip.open(filename, 'rb')
            else:
                eprint("Please check the filename, which is supposed to end with '.fq.gz' "
                      "or '.fastq.gz' if it's compressed or with '.fq' or '.fastq'.")
                sys.exit(1)
        elif file_suffix_2 == 'fq' or file_suffix_2 == 'fastq':
            f = open(filename, 'rb')
        else:
            eprint("Please check the filename, which is supposed to end with '.fq.gz' "
                  "or '.fastq.gz' if it's compressed or with '.fq' or '.fastq'.")
            sys.exit(1)
        return f
    #Size of every file.
    size_suffix = {'G': 3, 'M': 2, 'K': 1}          ##the value of power
    if size[-1] not in size_suffix.keys():
        eprint("The format of size for splitted files are supposed to be like 3G, 6M or 9K. "
              "Default is 3G.")
        sys.exit(1)
    else:
        try:
            size_value = 1024 ** size_suffix[size[-1]] * int(size[:-1]) ##the maxium size of uncompressed file.
        except ValueError:
            eprint("The number before 'G\M\K' should be integer.")
            sys.exit(1)
    #Split
    f1 = _filetype(fq1)
    if fq2:
        f2 = _filetype(fq2)
    dir_num , size_add = 1 , 0
    try:
        os.system('mkdir ' + outdir + '/fq' + str(dir_num))
    except:
        pass
    fw1 = open(outdir + '/fq' + str(dir_num) + '/1.fq', 'wb')
    if fq2:
        fw2 = open(outdir + '/fq' + str(dir_num) + '/2.fq', 'wb')
    line = f1.readline()
    while line:
        t = line
        for i in range(3):
            t += f1.readline()
        size_add += len(t)
        fw1.write(t)
        if fq2:
            t = ''
            for i in range(4):
                t += f2.readline()
            size_add += len(t)
            fw2.write(t)
        if size_add >= size_value:
            fw1.close()
            if fq2:
                fw2.close()
            subprocess.Popen(['gzip', outdir + '/fq' + str(dir_num) + '/1.fq'])
            subprocess.Popen(['gzip', outdir + '/fq' + str(dir_num) + '/2.fq'])
            dir_num += 1
            size_add = 0
            os.system('mkdir ' + outdir + '/fq' + str(dir_num))
            fw1 = open(outdir + '/fq' + str(dir_num) + '/1.fq', 'wb')
            if fq2:
                fw2 = open(outdir + '/fq' + str(dir_num) + '/2.fq', 'wb')
        line = f1.readline()
    fw1.close()
    if fq2:
        fw2.close()
    subprocess.Popen(['gzip', outdir + '/fq' + str(dir_num) + '/1.fq'])
    subprocess.Popen(['gzip', outdir + '/fq' + str(dir_num) + '/2.fq'])

    return 0

def IsGzipFile(filename):
    '''
        If it's a GIP file then return the file handle.
    '''
    with open(filename, 'rb') as f:
        if f.read(2) == b'\037\213':
            fg = gzip.open(filename, 'rt')
        else:
            fg = open(filename, 'r')
    return fg

def Fa2Geno(filename, chromosome=None):
    '''
        From fasta format file handle to Geno = {}.
    '''
    def _iterator(f):
        '''
            Generate the iterator of fastq file.
            '''
        seq = ''
        for line in f:
            fd = line.strip().split()
            if fd[0][0] == '>':
                if seq:
                    yield (key, seq, len(seq))
                key = fd[0][1:]
                seq = ''
            else:
                seq += fd[0]
        yield (key, seq, len(seq))

    f = IsGzipFile(filename)
    if chromosome is None:
        return _iterator(f)
    else:
        for line in _iterator(f):
            if line[0] == chromosome:
                return line
        return ()

def splitFa(fn, outdir, bases=50000000):
    def _wr(dt, order):
        fw1 , fw2 = open(outdir+os.path.basename(fn)+'.'+str(order)+'.fa', 'w'), open(outdir+os.path.basename(fn)+'.'+str(order)+'.bed', 'w')
        for k,v in dt.items():
            fw1.write('>'+k+'\n')
            fw1.write(v[1])
            fw2.write(k.split()[0]+'\t0\t'+str(v[0])+'\n')
        fw1.close()
        fw2.close()
        return {} ,  order+1
    f = IsGzipFile(fn)
    seq , dt = '' , {}
    for line in f:
        if line[0] == '>':
            if seq:
                dt[key] = [leng, seq]
            key = line[1:-1]
            seq , leng = '' , 0
        else:
            seq += line
            leng += len(line) - 1
    if seq:
        dt[key] = [leng, seq]
    f.close()
    base , order , dtt = 0 , 1 , {}
    for k,v in sorted(dt.items(), key=lambda x: x[1][0], reverse=True):
        base += v[0]
        dtt.update({k: v})
        if base >= bases:
            dtt , order = _wr(dtt, order)
            base = 0
    if dtt:
        dtt,order = _wr(dtt, order)
    
    return 0

def checkFqQuality(fn):
    f = IsGzipFile(fn)
    minx = 128
    for idx,line in enumerate(f, 1):
        if idx >= 40000:
            break
        if idx % 4 == 0:
            minx = min([x for x in map(ord, list(line.strip()))]+[minx])
    f.close()
    if minx < 64:
        return 33
    else:
        return 64

def QualityProbability(Qual, Phred=33, Type="Illumina"):
    '''
        
    '''
    if Type.upper() == 'ILLUMINA':
        return 10 ** (-(ord(Qual) - Phred) / 10)
    elif Type.upper() == 'SOLEXA':
        p = 10 ** (-(ord(Qual) - Phred) / 10)
        return p / (1 + p)
    else:
        eprint('not supported')
        sys.exit(1)

def isbwaindex(dirname, TYPE_GEO, TT):
    suffix , l = ['.amb', '.ann', '.bwt', '.pac', '.sa'] , []
    t = os.listdir(dirname)
    for i in list(TYPE_GEO)+[TT[i] for i in TYPE_GEO]:
        for j in suffix:
            if i+'.fa'+j in t:
                l.append(i+'.fa'+j)
            else:
                logging.info('do not have'+i+'.fa'+j+'file in --geo')
                return []
    logging.info('--geo directory have all necessary index file')
    return l

def ishisatindex(dirs):
    if dirs == '':
        dirs = './'
    for fn in os.listdir(dirs):
        if os.path.splitext(fn)[1] == '.ht2':
            logging.info('hisat index file exist')
            return 1
    logging.info('hisat index file not exist')
    return 0

def issamtoolsindex(tp):
    if tp == 'geo':
        ['.fai']
    elif tp == 'bam':
        ['.bai']
    else:
        eprint('')
        sys.exit(1)

def eprint(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)

def _BasePercent(filename):
    BaseContent = {'A': 0, 'C': 0, 'G': 0, 'T': 0}
    a = ('A', 'C', 'G', 'T')
    for line in IsGzipFile(filename):
        line = line.strip().upper()
        if line[0] != '>':
            for b in a:
                BaseContent[b] += line.count(b)
    num = 0
    for key in BaseContent:
        num += BaseContent[key]
    for key in BaseContent:
        BaseContent[key] = '%.4f'%(BaseContent[key]/num)
    return BaseContent

def _DealBasequality(Bases, Quals, Phred):
    bases_new , quals_new = [] , []
    for i in range(len(Bases)):
        score = ord(Quals[i]) - Phred
        if score >= Qual_cutoff and Bases[i] != 'N':
            bases_new.append(Bases[i])
            quals_new.append(Quals[i])
    return ''.join(bases_new) , ''.join(quals_new)

def _FormatP(value):
    if value > 0.0001:
        return '%.8f'%value
    else:
        return '%.8e'%value

def _Fdr(ar):
    '''
        FDR examination
    '''
    fdr_result , n = [] , len(ar)
    for i, pvalue in enumerate(ar):
        if i == 0:
            temp = n
        elif pvalue < pre:
            temp = n - i
        pre = pvalue
        fdr_result.append(pvalue * n / temp)
    return fdr_result

def _SNPvalue_Bayesian(Bases, Quals, Ploid, FixedKey, FixedError, BaseContent):
    '''
        Bayesian, p(G|D) = p(D|G)*p(G)/p(D)
    '''
    ## likeilyhood
    p_allele = {'A': [], 'C': [], 'G': [], 'T': []}
    for idx , base in enumerate(Bases):
        if base in T_BASE_ORDER:
            score = QualityProbability(Quals[idx])
            for b in T_BASE_ORDER:
                if b == base:
                    p_allele[b].append(1 - score)
                else:
                    p_allele[b].append(score * FixedError[b][base])
        else:
            continue
    ## prior homo
    prior , typeprior = {} , {}
    for base in T_BASE_ORDER:
        genotype = base * Ploid
        prior[genotype] = 1
        for score in p_allele[base]:
            prior[genotype] += log(score)
        typeprior[genotype] = HomoPrior * float(BaseContent[base])
    heteNum = (Ploid - 1) * 4 + (Ploid - 1) * 2 * Ratio
    ## prior hete
    ### Warning: There needs some check!!!
    if Ploid > 1:
        for key , ke in FixedKey:
            for i in range(1, Ploid):
                value = 1
                genotype = key * i + ke * (Ploid - i)
                for idx , j in enumerate(p_allele[key]):
                    value += log((i*j + (Ploid-i)*p_allele[ke][idx]) / Ploid)
                prior[genotype] = value
                if (key == 'G' and ke == 'A') \
                    or (key == 'A' and ke == 'G') \
                    or (key == 'T' and ke == 'C') \
                    or (key == 'C' and ke == 'T'):
                    typeprior[genotype] = HetePrior/heteNum*Ratio
                else:
                    typeprior[genotype] = HetePrior/heteNum
    ## caculate posterior
    ### p(D)
    p_pileup = 0
    greatestmax = max([value for key , value in prior.items()])
    for gt in prior.keys():
        prior[gt] -= greatestmax    # ?
        try:
            p_pileup += typeprior[gt] * exp(prior[gt])
        except KeyError:
            eprint('['+PROGRAM+'] Error: there is not the prior probability for the '+genotype+'\n')
            sys.exit(1)
    ### p(G|D)
    posterior = {}
    for gt , value in prior.items():
        try:
            posterior[gt] = typeprior[gt] * exp(value) / p_pileup
        except KeyError:
            eprint('['+PROGRAM+'] Error: there is not the prior probability for the '+gt+'\n')
            sys.exit(1)
    return sorted(posterior.items(), key=lambda x:x[1], reverse=True)[0]

def _SNVvalue_Binomial(Info, Refbase, RNA=False):
    '''
        Binomial method.
    '''
    def __dbinom(x, n, pie):
        '''
            C(m,n)*p(m)q(n-m)
        '''
        if x > n/2:
            a , b = x , n - x
        else:
            a , b = n - x , x
        log_fac1 = 0
        i = n
        while i > a:
            log_fac1 += log(i)
            i -= 1
        log_fac2 = 0
        i = b
        while i >= 1:
            log_fac2 += log(i)
            i -= 1
        return exp((log_fac1-log_fac2) + x*log(pie) + (n-x)*log(1-pie))
    def __pbinom(x, n, pie, tail):
        '''
            Sum
        '''
        tail_p = 0
        if tail == 'lower':
            if x <= n/2:
                for i in range(x+1):
                    tail_p += __dbinom(i, n, pie)
            else:
                for i in range(x+1, n+1):
                    tail_p += __dbinom(i, n, pie)
                tail_p = 1 - tail_p
        elif tail == 'upper':
            if x <= n/2:
                for i in range(x):
                    tail_p += __dbinom(i, n ,pie)
                tail_p = 1 - tail_p
            else:
                for i in range(x, n+1):
                    tail_p += __dbinom(i, n, pie)
        else:
            eprint('['+PROGRAM+'] Error: tail is not set correctly, set it to lower or upper\n')
            sys.exit(1)
        if tail_p <= 0:
            tail_p = __dbinom(x, n, pie)
        return tail_p
    ## core
    if not RNA:
        diffdep , totaldep = 0 , 0
        dt = {0: 'A', 1: 'C', 2: 'G', 3: 'T'}
        for i in range(4):
            totaldep += Info[i]
            diffdep += Info[i] if Refbase != dt[i] else 0
        return __pbinom(diffdep, totaldep, 1/Ploid, "lower")
    else:
        diffdep , totaldep = Info
        pie = Refbase
        return __pbinom(diffdep, totaldep, pie, "upper")

def _SNPvalue_Frequency(Bases, Refbase):
    x , y = len(Bases) , 0
    for b in Bases:
        y += 1 if b != Refbase else 0
    return y ,'%.4f'%(y/x)
