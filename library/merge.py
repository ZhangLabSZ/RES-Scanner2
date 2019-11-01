#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
# Author: guoqunfei, guoqunfei@genomics.cn

import os
import re
import sys
import gzip
import subprocess
import argparse

import shell


def merge_RES(Config, Genomefile, Phred_DNA=33, Phred_RNA=33, 
    Qual_cutoff=30, HomoPrior=0.99, Rate=2, Method='Bayesian', 
    Ploidy=2, Intron=None, DNAdepth=10, RNAdepth=3, Bayesian_Posterior_Probability=0.95, 
    FDR_DNA_Heterozygosis=0.05, Non_Ref_BaseCount=0, Paralogous_D=1, Intronic=6, 
    Homopolymer=1,out_path='./'):
    '''
        merge the result of scanner and check
    '''

    BT = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 'N': 'N'}
    BASE = ('A', 'C', 'G', 'T')

    def read(f, s, table, several_tp):
        for line in f:
            if line[0] == '#':
                continue
            fd = line.strip().split()
            key1,key2,key3,key4,key5 = (fd[0], fd[1], fd[2]),(s, 'DNA'),(s,'RNA'),(s,'SNPvalue'),(s,'Editvalue')
            try:
                table[key1].update({key2: fd[5]})
            except KeyError:
                table[key1] = {key2: fd[5]}
            table[key1][key3] = fd[9]
            if Method == 'Bayesian':
                table[key1][key4] = (fd[6], fd[7])
            else:
                table[key1][key4] = fd[6]
            table[key1][key5] = fd[12]
            table[key1]['gbase'] = fd[3].upper()
            try:
                if table[key1]['type'] != fd[10]:
                    several_tp.add(key1)
            except KeyError:
                table[key1]['type'] = fd[10]
        return 0


    def deleteintron(sca, l, table):
        for i in l:
            try:
                del table[(sca, str(i), '+')]
            except KeyError:
                pass
            try:
                del table[(sca, str(i), '-')]
            except KeyError:
                pass
            try:
                del table[(sca, str(i), '.')]
            except KeyError:
                pass
        return 0


    def caculus(f, table, errorRate, st):
        for line in f:
            fd = line.strip().split()
            ref = fd[2]
            for strand in st:
                key1 = (fd[0], fd[1], strand)
                if key1 in table:
                    key2 = (s, 'RNA')
                    if key2 not in table[key1]:
                        editV = (s, 'Editvalue')
                        if fd[5] == 0:
                            table[key1][key2] = '0,0,0,0'
                            table[key1][editV] = '1'
                        else:
                            seq , qual = shell._DealBasequality(fd[3][4:4+int(fd[-1])].upper(), fd[4][4:4+int(fd[-1])], Phred_RNA)
                            l = [seq.count(i) for b in BASE]
                            table[key1][key2] = '{},{},{},{}'.format(*l)
                            tol , editdep = 0 , 0
                            for i in range(4):
                                tol += l[i]
                                if BASE[i] != ref and l[i] > editdep:
                                    editdep = l[i]
                            table[key1][editV] = shell._SNVvalue_Binomial((editdep, tol), errorRate, RNA=True)
        return 0


    if Method not in ('Bayesian', 'Binomial', 'Frequency'):
        print_help()
        sys.exit(1)
    elif Method == 'Bayesian' and Genomefile == None:
        print_help()
        sys.exit(1)
    HetePrior = 1 - HomoPrior
    ## read config file
    sample , order = {} , []
    with open(Config, 'r') as f:
        for line in f:
            fd = line.strip().split()
            order.append(fd[0])
            if os.path.isfile(fd[1]+'.stat') == True:
                sample[fd[0]] = {3: fd[1]+'.stat'}
            else:
                shell.eprint('[bigtable] Error: DNA single base file do not have .stat file in the same directory')
                sys.exit(1)
            sample[fd[0]][0] = fd[1]
            if len(fd) < 4:
                shell.eprint('[bigtable] Error: config file format error')
                sys.exit(1)
            sample[fd[0]][1.1] = fd[2]
            sample[fd[0]][1.2] = fd[3]
            if len(fd) == 6:
                sample[fd[0]][2.1] = fd[4]
                sample[fd[0]][2.2] = fd[5]
    ## read table
    several_tp , table , peak_dep = set() , {} , {}
    for s,v in sample.items():
        f = shell.IsGzipFile(v[1.2])
        f.readline()
        read(f, s, table, several_tp)
        f.close()
        with open(v[3], 'r') as f:
            fd = f.readlines()[1].strip().split()
            peak_dep[s] = float(fd[4]) if float(fd[4]) > float(fd[3]) else float(fd[3])
        try:
            f = shell.IsGzipFile(v[2.2])
            f.readline()
            read(f, s, table, several_tp)
            f.close()
        except:
            pass
    ## delete sites have different editing type in different sample
    for key in several_tp:
        del table[key]
    ## filter sites locating near junctions
    if Intron != None:
        f = shell.IsGzipFile(Intron)
        for line in f:
            fd = line.strip().split()
            if int(fd[3]) < int(fd[4]):
                beg , end = int(fd[3]) , int(fd[4])
            else:
                beg , end = int(fd[4]) , int(fd[3])
            deleteintron(fd[1], range(beg, beg+Intronic), table)
            deleteintron(fd[1], range(end-Intronic+1, end+1), table)
        f.close()
    ## filter homopolymer
    if Homopolymer:
        seq , leng , st = {} , {} , ''
        m = [re.compile(i*5) for i in BASE]
        f = shell.IsGzipFile(Genomefile)
        for line in f:
            if line[0] == '>':
                if st:
                    seq[key] = st
                    leng[key] = len(st)
                st = ''
                key = line.strip().split()[0][1:]
            else:
                st += line.strip().upper()
        if st:
            seq[key] = st
            leng[key] = len(st)
        delkey = set()
        for key in table:
            sca , pos , strand = key
            beg = int(pos)-4 if int(pos)-4 > 1 else 1
            end = int(pos)+4 if int(pos)+4 < leng[sca] else leng[sca]
            nt = seq[sca][beg-1: end]
            for m1 in m:
                if m1.match(nt):
                    delkey.add(key)
                    break
        for key in delkey:
            del table[key]
        del seq
        del leng
        del delkey

    if Method == 'Bayesian':
        basecontent = shell._BasePercent(Genomefile)
        weight = 0.5
        weight_other = (1 - weight) / 2
        ## adjust substitution rate for illumina
        FixedError = {'A': {'C': weight, 'T': weight_other, 'G': weight_other},
                    'C': {'A': weight, 'T': weight_other, 'G': weight_other},
                    'G': {'T': weight, 'A': weight_other, 'C': weight_other},
                    'T': {'G': weight, 'A': weight_other, 'C': weight_other}}
        FixedKey = []
        for key in FixedError.keys():
            for ke in FixedError[key].keys():
                if (key, ke) not in FixedKey and (ke, key) not in FixedKey:
                    FixedKey.append((key, ke))
    errorRate = 10**(-1*Qual_cutoff/10)
    for s,v in sample.items():
        with gzip.open(v[0], 'rt') as f:
            for line in f:
                fd = line.strip().split()
                ref = fd[2].upper()
                for strand in ('+', '-', '.'):
                    key1 = (fd[0], fd[1], strand)
                    if key1 in table:
                        key2 = (s, 'DNA')
                        if key2 not in table[key1]:
                            snpV = (s, 'SNPvalue')
                            if fd[5] == 0:
                                table[key1][snpV] = 'NA'
                                table[key1][key2] = '0,0,0,0'
                            else:
                                seq , qual = shell._DealBasequality(fd[3][4:4+int(fd[-1])].upper(), fd[4][4:4+int(fd[-1])], Phred_DNA)
                                table[key1][key2] = ','.join([str(seq.count(i)) for b in BASE])
                                if seq == '':
                                    table[key1][snpV] = 'NA'
                                    continue
                                if Method == 'Bayesian':
                                    result = shell._SNPvalue_Bayesian(seq, qual, Ploidy, FixedKey, FixedError, basecontent)
                                    table[key1][snpV] = (result[0], result[1])
                                elif Method == 'Binomial':
                                    table[key1][snpV] = shell._SNVvalue_Binomial(','.join([str(seq.count(i)) for b in BASE]), ref)
                                elif Method == 'Frequency':
                                    table[key1][snpV] = shell._SNPvalue_Frequency(seq, ref)[1]
                                else:
                                    shell.eprint('[bigtable] Error: --method not recognize')
                                    sys.exit(1)
        with gzip.open(v[1.1], 'rt') as f:
            caculus(f, table, errorRate, ['+', '.'])
        if 2.1 not in v.keys():
            continue
        with gzip.open(v[2.1], 'rt') as f:
            caculus(f, table, errorRate, ['-'])
    ## Complement DNA and RNA information for all sites in the table
    for k,v in table.items():
        for spl in order:
            keyDNA,keysnpV,keyRNA,keyeditV = (spl,'DNA'),(spl,'SNPvalue'),(spl,'RNA'),(spl, 'Editvalue')
            if keyDNA not in v:
                table[k][keyDNA] = '0,0,0,0'
                table[k][keysnpV] = 'NA'
            if keyRNA not in v:
                table[k][keyRNA] = '0,0,0,0'
                table[k][keyeditV] = '1'
    ## Remove sites with high DNA depth and multiple RNA editing types
    delkey = set()
    for key1 in table:
        for s in order:
            fd = [int(i) for i in table[key1][(s, 'DNA')].split(',')]
            if Paralogous_D:
                dep = sum(fd)
                if dep > 2*peak_dep[s]:
                    delkey.add(key1)
                    break
            fd = [int(i) for i in table[key1][(s, 'RNA')].split(',')]
            rna_count = {}
            for i in range(4):
                if BASE[i] != table[key1]['gbase']:
                    rna_count[BASE[i]] = fd[i]
            key3 = sorted(rna_count.keys(), key=lambda x: rna_count[x], reverse=True)
            if len(key3) != 3:
                shell.eprint('[bigtable] Error: rna depth error')
                sys.exit(1)
            if rna_count[key3[0]] > 0 and (rna_count[key3[1]] / rna_count[key3[0]]) > 0.01:
                delkey.add(key1)
                break
    for key in delkey:
        del table[key]

    delkey = set()
    for key1 in table:
        kplus , kminus = (k[0], k[1], '+') , (k[0], k[1], '-')
        if kplus in table and kminus in table:
            plusdep , miusdep = 0 , 0
            for s in order:
                plusdep += sum([int(i) for i in table[kplus][(s, 'RNA')].split(',')])
                miusdep += sum([int(i) for i in table[kminus][(s, 'RNA')].split(',')])
            if plusdep > miusdep:
                delkey.add(kminus)
            elif plusdep < miusdep:
                delkey.add(kplus)
    for key in delkey:
        del table[key]
    del delkey

    fdr , binomial = {} , {}
    for key1 in table:
        for s in order:
            try:
                fdr[s].append([key1, table[key1][(s, 'Editvalue')]])
            except KeyError:
                fdr[s] = [[key1, table[key1][(s, 'Editvalue')]]]
            if Method == 'Binomial':
                try:
                    binomial[s].append([key1, table[key1][(s, 'SNPvalue')]])
                except KeyError:
                    binomial[s] = [[key1, table[key1][(s, 'SNPvalue')]]]
    for s,v in fdr.items():
        fd = sorted(v, key=lambda x: float(x[1]), reverse=True)
        fd_fdr = shell._Fdr([float(i[1]) for i in fd])
        for i in range(len(fd_fdr)):
            table[fd[i][0]][(s, 'Editvalue')] = fd_fdr[i]
    least_dep = 1
    if Method == 'Binomial':
        for s,v in binomial.items():
            fd = sorted(v, key=lambda x: float(x[1]), reverse=True)
            fd_fdr = shell._Fdr([float(i[1]) for i in fd])
            for i in range(len(fd_fdr)):
                table[fd[i][0]][(s, 'SNPvalue')] = shell._FormatP(fd_fdr[i])
        p = sorted(peak_dep.keys(), key=lambda x: peak_dep[x])
        ratio = 1 / Ploidy
        while least_dep < p[0]:
            if shell._SNVvalue_Binomial((0, least_dep), ratio, RNA=True) < 0.05:
                break
            least_dep += 1

    fw = open(out_path + '/RES.txt', 'w')
    title = '#1.Chromosome\t2.Coordinate\t3.Strand\t4.Gbase\t5.EditType'
    for idx,s in enumerate(order):
        title += '\t'+str(6+idx*2)+'.'+s+'.DNA_baseCount[A,C,G,T]\t'+str(7+idx*2)+'.'+s+'.RNA_basecount[A,C,G,T];P_value'
    title += '\n'
    fw.write(title)
    for key1 in sorted(table.keys(), key=lambda x: (x[0], int(x[1]))):
        info = []
        flag = 1
        for s in order:
            rnainfo = table[key1][(s, 'RNA')].split(',')
            rna_dep = sum([int(i) for i in rnainfo])
            dna_dep = sum([int(i) for i in table[key1][(s, 'DNA')].split(',')])
            if table[key1][(s, 'SNPvalue')] != 'NA':
                if Method == 'Bayesian':
                    if len(set(table[key1][(s, 'SNPvalue')][0])-set(table[key1]['gbase'])) != 0 or float(table[key1][(s, 'SNPvalue')][1]) < Bayesian_Posterior_Probability: flag = 0
                elif Method == 'Binomial':
                    if dna_dep < least_dep:
                        n = 0
                        for i in table[key1][(s, 'RNA')].split(','):
                            if int(i) == 0: n += 1
                        if n < 3: flag = 0
                    else:
                        if table[key1][(s, 'SNPvalue')] < FDR_DNA_Heterozygosis: flag = 0
                elif Method == 'Frequency':
                    if table[key1][(s, 'SNPvalue')] < Non_Ref_BaseCount: flag = 0
            if dna_dep >= DNAdepth and rna_dep >= RNAdepth and float(table[key1][(s, 'Editvalue')]) < 0.05:
                alt_base = table[key1]['type'].split('->')[1]
                if key1[2] == '-':
                    alt_base = BT[alt_base]
                base_dep = {}
                for i in range(4):
                    base_dep[BASE[i]] = int(rnainfo[i])
                if base_dep[alt_base] > 0:
                    info.append(table[key1][(s, 'DNA')]+'\t'+table[key1][(s, 'RNA')]+';'+str(table[key1][(s, 'Editvalue')])+'*')
                else:
                    info.append(table[key1][(s, 'DNA')]+'\t'+table[key1][(s, 'RNA')]+';'+str(table[key1][(s, 'Editvalue')]))
            else:
                info.append(table[key1][(s, 'DNA')]+'\t'+table[key1][(s, 'RNA')]+';'+str(table[key1][(s, 'Editvalue')]))
        if flag == 1:
            fw.write(('{}\t'*5).format(*(list(key1)+[table[key1]['gbase'], table[key1]['type']]))+'\t'.join(info)+'\n')

    return 0


def merge_RHES(table_path,RNA_pos_sbd,RNA_neg_sbd,out_path='./'):

    def getting_RHES(RHES_files_lst):

        RHES_dt,title = [{},[]]

        for RHES_file in RHES_files_lst:
            with gzip.open(RHES_file,'rb') as f_in:
                for line in f_in:
                    l = line.decode().split()
                    if l[0] == '1.Chromosome_ID':
                        title = line.decode().rstrip('\n')
                    else:
                        RHES = '_'.join(l[0:3])
                        try:
                            pop_obj = RHES_dt.pop(RHES)
                        except KeyError:
                            RHES_dt[RHES] = l[3:]
        return RHES_dt,title


    def getting_handle(sbd_lst):

        A,C,G,T = [0,0,0,0]
        base_str = sbd_lst[0].split(':')[1]
        if base_str[-1] == ']':
            base_str = base_str[0:-4]
        base_lst = list(base_str)
        qual_lst = list(map(ord,list(sbd_lst[1][4:4+len(base_lst)])))
        for i in range(len(qual_lst)):
            if qual_lst[i] >= 63:
                if base_lst[i] == 'A':
                    A += 1
                elif base_lst[i] == 'C':
                    C += 1
                elif base_lst[i] == 'G':
                    G += 1
                elif base_lst[i] == 'T':
                    T += 1
        return [A,C,G,T]


    def getting_RHES_sbd_dt(RNA_sbd,RHES_sbd_dt,strand,RHES_dt,complement_dt):

        with gzip.open(RNA_sbd,'rb') as f_in:
            for line in f_in:
                l = line.decode().split()
                key = '_'.join([l[0],l[1],strand])
                try:
                    pop_obj = RHES_sbd_dt.pop(key)
                    RHES_sbd_dt[key] = getting_handle(l[-3:])
                except KeyError:
                    pass

        for k,v in RHES_sbd_dt.items():
            try:
                lst = RHES_dt[k]
            except KeyError:
                continue
            else:
                TYPE,A_,C_,G_,T_ = [lst[7]] + lst[6].split(',')
                ref_base = complement_dt[TYPE.split('>')[-1]]
                A,C,G,T = [v[0]+int(A_),v[1]+int(C_),v[2]+int(G_),v[3]+int(T_)]
                base_dt = {'A':A,'C':C,'G':G,'T':T}
                num = str(A+C+G+T)
                RATE = '%.4f'%(float(base_dt[ref_base])/float(num))
                RHES_sbd_dt[k] = lst[0:5] + [num,','.join([str(A),str(C),str(G),str(T)]),TYPE,RATE] + lst[9:]
        return RHES_sbd_dt


    def getting_sbd(RHES_dt,RNA_pos_sbd,RNA_neg_sbd):

        RHES_pos_dt,RHES_neg_dt = [{},{}]

        for RHES in RHES_dt.keys():
            strand = RHES.split('_')[-1]
            if strand == '+':
                RHES_pos_dt[RHES] = [0,0,0,0]
            elif strand == '-':
                RHES_neg_dt[RHES] = [0,0,0,0]

        pos_complement_dt = {'A':'A','C':'C','G':'G','T':'T'}
        neg_complement_dt = {'A':'T','C':'G','G':'C','T':'A'}
        RHES_pos_dt = getting_RHES_sbd_dt(RNA_pos_sbd,RHES_pos_dt,'+',RHES_dt,pos_complement_dt)
        RHES_neg_dt = getting_RHES_sbd_dt(RNA_neg_sbd,RHES_neg_dt,'-',RHES_dt,neg_complement_dt)

        return RHES_pos_dt,RHES_neg_dt


    def getting_RHES_merge_dt(RHES_pos_dt,RHES_neg_dt):

        RHES_merge_dt = {}

        for k,v in RHES_pos_dt.items():
            key = '_'.join(k.split('_')[:-1])
            RHES_merge_dt[key] = ['+'] + v

        for k,v in RHES_neg_dt.items():
            key = '_'.join(k.split('_')[:-1])
            try:
                lst = RHES_merge_dt[key]
            except KeyError:
                RHES_merge_dt[key] = ['-'] + v
            else:
                if int(v[5]) > int(lst[6]):
                    RHES_merge_dt[key] = ['-'] + v
                elif int(v[5]) == int(lst[6]):
                    pop_obj = RHES_merge_dt.pop(key)
                else:
                    continue
        return RHES_merge_dt


    RHES_files_lst = []
    for i in os.listdir(table_path):
        if 'homo.filter.gz' in i:
            RHES_files_lst.append(table_path+'/'+i)

    RHES_dt,title = getting_RHES(RHES_files_lst)
    RHES_pos_dt,RHES_neg_dt = getting_sbd(RHES_dt,RNA_pos_sbd,RNA_neg_sbd)
    RHES_merge_dt = getting_RHES_merge_dt(RHES_pos_dt,RHES_neg_dt)

    out = open(os.getcwd() + '/RHES.txt','w')
    print(title,file=out)
    for k in sorted(RHES_merge_dt.keys()):
        chro,pos = ['_'.join(k.split('_')[:-1]),k.split('_')[-1]]
        print('\t'.join([chro,pos] + RHES_merge_dt[k]),file=out)

    return 0


def help():

    def res(typ):

        add_dt = {}
        if typ == 'res':
            add_dt[typ] = res_parser.add_argument
        else:
            add_dt[typ] = aes_parser.add_argument
        add_dt[typ]('--config',metavar='FILE',type=str,required=True,help="Each line of the tab-delimited configuration table represents a sample and contains six columns:\nColumn 1: Sample name that used for RES-Scanner2 pipeline.\nColumn 2: The file path of the DNA single base file generated by preprocess pipeline for the corresponding sample.\nColumn 3: The file path of the RNA single base depth file on the positive strand generated by preprocess pipeline for the corresponding sample.\nColumn 4: The file path of RNA editing sites file on the positive strand generated by scanner pipeline for the corresponding sample.\nColumn 5: The file path of RNA single base depth file on the negative strand generated by preprocess pipeline for the corresponding sample.\nColumn 6: The file path of RNA editing sites file on the negative strand generated by scanner pipeline for the corresponding sample.")
        add_dt[typ]('--phred',metavar='',help="The Phred base quality for query QUALity of DNA.bam and RNA.bam files respectively, default DNA_ASCII-33,RNA_ASCII-33. [33,33]")
        add_dt[typ]('--qual_cutoff',metavar='',type=int,default=30,help="Quality cutoff for BWA alignment. [30]")
        add_dt[typ]('--method',metavar='',type=str,default='Bayesian',help="Method for detecting SNPs.'Bayesian' or 'Binomial' or 'Frequency'. [Bayesian]")
        add_dt[typ]('--HomoPrior',metavar='',type=float,default=0.99,help="The prior probability of homozygous genomic positions. (force --method Bayesian) [0.99]")
        add_dt[typ]('--rate',metavar='',type=int,default=2,help="The rate of transitions over transversions of the genome. (force --method Bayesian) [2]")
        add_dt[typ]('--Bayesian_P',metavar='',type=float,default=0.95,help="The minimun Bayesian Posterior Probability cut off for corresponding genotype, range from 0 to 1. (force --method Bayesian) [0.95]")
        add_dt[typ]('--Binomial_FDR',metavar='',type=float,default=0.05,help="The maximun FDR cutoff of Binomial test for DNA Heterozygosis, range from 0 to 1. (force --method Binomial) [0.05]")
        add_dt[typ]('--Frequency_N',metavar='',type=int,default=0,help="The maximun non-reference base count cutoff. (force --method Frequency) [0]")
        add_dt[typ]('--genome',metavar='FILE',type=str,required=True,help="The refference genome.")
        add_dt[typ]('--ploidy',metavar='',type=int,default=2,help="The ploidy level, 1 for monoploid , 2 for diploid, 3 for triploid, 4 for tetraploid, and so on. [2].")
        add_dt[typ]('--intron',metavar='FILE',type=str,default='null',help="The intron region file with POS format. This option is called for filter editing sites locating near junctions. [null]")
        add_dt[typ]('--intronic',metavar='NUM',type=int,default=6,help="Remove intronic candidate editing sites occurring within (the number of) bases of a splice site. [6]")
        add_dt[typ]('--DNAdepth',metavar='',type=int,default=10,help="Coverage of DNA base, the site which is less than the coverage will not be marked as edit site in corresponding sample.[default 10]")
        add_dt[typ]('--RNAdepth',metavar='',type=int,default=3,help="Coverage of RNA base, the site which is less than the coverage will not be marked as edit site in corresponding sample.[default 3]")
        add_dt[typ]('--homopolymer',metavar='NUM',type=int,default=1,help="Whether remove candidate editing sites in homopolymer runs of >= 5 base pairs. 1 for yes, 0 for not. [default 1]")
        add_dt[typ]('--paralogous_D',metavar='NUM',type=int,default=1,help="Whether discard candidate editing sites with DNA reads depth of more than twice the genome-wide peak or mean depth.\t1 for yes, 0 for not. [default 1]")
        return 0


    def rhes(typ):

        add_dt = {}
        if typ == 'rhes':
            add_dt[typ] = rhes_parser.add_argument
        else:
            add_dt[typ] = aes_parser.add_argument
        add_dt[typ]('--table',dest='t',metavar='table_path',type=str,required=True,help="The path to the files that generate the hyper editing site.")
        add_dt[typ]('--pos',dest='p',metavar='RNA_pos_sbd',type=str,required=True,help="The RNA positive single base depth file.")
        add_dt[typ]('--neg',dest='n',metavar='RNA_neg_sbd',type=str,required=True,help="The RNA negative single base depth file.")
        return 0


    parser = argparse.ArgumentParser(prog='ires merge',usage='%(prog)s <cmd> [option] arg...',formatter_class=argparse.RawDescriptionHelpFormatter,description='Where <cmd> is one of : res ,rhes , aes .')
    parser.add_argument('-v','--version',action='version',version='RES-Scanner2 0.1.0')
    commands = parser.add_subparsers(help='sub-command help')
    res_parser = commands.add_parser('res',prog='ires merge res',usage='%(prog)s [option] arg...')
    rhes_parser = commands.add_parser('rhes',prog='ires merge rhes',usage ='%(prog)s [option] arg...')
    aes_parser = commands.add_parser('aes',prog='ires merge aes',usage='%(prog)s [option] arg...')

    if len(sys.argv) < 2:
        print("Please enter the correct command or -h/--help for help")
        sys.exit(1)
    elif sys.argv[1] == "res":
        res('res')
        if len(sys.argv) < 3:
            print("Please enter the correct parameters or -h/--help for help")
            sys.exit(1)
    elif sys.argv[1] == "rhes":
        rhes('rhes')
        if len(sys.argv) < 3:
            print("Please enter the correct parameters or -h/--help for help")
            sys.exit(1)
    elif sys.argv[1] == "aes":
        rhes('aes')
        res('aes')
        if len(sys.argv) < 3:
            print("Please enter the correct parameters or -h/--help for help")
            sys.exit(1)
    elif sys.argv[1] == "-h" or sys.argv[1] == "--help":
        pass
    else:
        print("Please enter the correct command or -h/--help for help")
    args = parser.parse_args()

    return args

def main():

    args = help()
    if sys.argv[1] == "res" or sys.argv[1] == "aes":
        Config,Genome = [args.config,args.genome]
        merge_RES(Config,Genome)
    elif sys.argv[1] == "rhes" or sys.argv[1] == "aes":
        table_path,RNA_pos_sbd,RNA_neg_sbd = [args.t,args.p,args.n]
        merge_RHES(table_path,RNA_pos_sbd,RNA_neg_sbd)


if __name__ == "__main__":
    main()
