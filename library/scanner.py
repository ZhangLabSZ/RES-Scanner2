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

from copy import deepcopy
from math import log, exp
from itertools import combinations

import shell

# basic env
SCRIPT = os.path.basename(__file__)
PROGRAM = os.path.splitext(SCRIPT)[0]
PATH = os.path.abspath(__file__)
DIR = os.path.dirname(PATH)

# output directory
OUTDIR = os.getcwd()

# software path
BLAT = ''

D_BASE_COMPLIMENT = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 'N': 'N'}
T_BASE_ORDER = ('A', 'C', 'G', 'T')

def candidateRES(RNA_singleBase, DNA_singleBase, Genomefile, Strand='Unknown', Coverage=2, Method='Bayesian', Ratio=4, Ploid=2, Phred_RNA=33, Phred_DNA=33, Qual_cutoff=30, HomoPrior=0.99, Output='./', Suffix='.homo.gz'):
    '''

    '''
    def _DealBasequality(Bases, Quals, Phred):
        bases_new , quals_new = [] , []
        for i in range(len(Bases)):
            score = ord(Quals[i]) - Phred
            if score >= Qual_cutoff and Bases[i] != 'N':
                bases_new.append(Bases[i])
                quals_new.append(Quals[i])
        return ''.join(bases_new) , ''.join(quals_new)
    def _BasePercent(filename):
        BaseContent = {'A': 0, 'C': 0, 'G': 0, 'T': 0}
        a = ('A', 'C', 'G', 'T')
        for line in shell.IsGzipFile(filename):
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
    def _SNPvalue_Bayesian(Bases, Quals):
        '''
            Bayesian, p(G|D) = p(D|G)*p(G)/p(D)
        '''
        ## likeilyhood
        p_allele = {'A': [], 'C': [], 'G': [], 'T': []}
        for idx , base in enumerate(Bases):
            if base in T_BASE_ORDER:
                score = shell.QualityProbability(Quals[idx])
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
                shell.eprint('['+PROGRAM+'] Error: there is not the prior probability for the '+genotype+'\n', file=sys.stderr)
                sys.exit(1)
        ### p(G|D)
        posterior = {}
        for gt , value in prior.items():
            try:
                posterior[gt] = typeprior[gt] * exp(value) / p_pileup
            except KeyError:
                shell.eprint('['+PROGRAM+'] Error: there is not the prior probability for the '+gt+'\n', file=sys.stderr)
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
                shell.eprint('['+PROGRAM+'] Error: tail is not set correctly, set it to lower or upper\n')
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
    # main
    ## Paramete tackle
    if Method == 'Bayesian':
        title1 , title2 = 'Bayesian_Genotype(DNA)' , 'Bayesian_Posterior_Probability(DNA)'
    elif Method == 'Binomial':
        title1 , title2 = 'P_value(DNA_Heterozygosis)' , 'FDR(DNA_Heterozygosis)'
    elif Method == 'Frequency':
        title1 , title2 = 'Non_Ref_BaseCount(DNA)' , 'Non_Ref_BaseRatio(DNA)'
    else:
        shell.eprint('['+PROGRAM+'] Error: Unknown --Method')
        sys.exit(1)
    head = "1.Chromosome_ID\t2.Coordinate\t3.strand\t4.Ref_base\t"\
        "5.coverage(DNA)\t6.DNA_BaseCount[A,C,G,T]\t7."+title1+"\t"\
        "8."+title2+"\t9.coverage(RNA)\t10.RNA_BaseCount[A,C,G,T]\t"\
        "11.RNA_Editing_Type\t12.RNA_Editing_Level\t13.P_value(RNA_editing)\t"\
        "14.FDR(RNA_editing)\n"
    HetePrior = 1 - HomoPrior
    if Method == 'Bayesian':
        if not os.access(Genomefile, os.R_OK):
            shell.eprint('['+PROGRAM+'] Error: genomefile is not available')
            sys.exit(1)
    ## read RNA sam2base file
    rna_editing = {}
    with gzip.open(RNA_singleBase, 'rt') as f:
        for line in f:
            fd = line.strip().split()
            if int(fd[-1]) >= Coverage and fd[3][:3] == '[M]' and fd[4][:3] == '[M]':
                bases = fd[3][4:4+int(fd[-1])].upper()
                quals = fd[4][4:4+int(fd[-1])]
                try:
                    bases_new , quals_new = _DealBasequality(bases, quals, Phred_RNA)
                except IndexError:
                    shell.eprint('['+PROGRAM+'] Error: please check the RNA single base file')
                    sys.exit(1)
                cov_new = len(bases_new)
                if cov_new < Coverage:
                    continue
                refbase = fd[2].upper()
                for b in T_BASE_ORDER:
                    num = bases_new.count(b) if b != refbase else 0
                    if num > 0:
                        try:
                            rna_editing[fd[0]].update({fd[1]: [refbase, bases_new, quals_new, cov_new]})
                        except KeyError:
                            rna_editing[fd[0]] = {fd[1]: [refbase, bases_new, quals_new, cov_new]}
                        break
    logging.info('rna single base done')
    ## BaseContent
    if Method == 'Bayesian':
        BaseContent = _BasePercent(Genomefile)
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
    ## read DNA sam2base file
    homo = {}
    with gzip.open(DNA_singleBase, 'rt') as f:
        for line in f:
            fd = line.strip().split()
            try:
                ref_base = rna_editing[fd[0]][fd[1]][0]
            except KeyError:
                continue
            if ref_base != fd[2].upper():
                shell.eprint('['+PROGRAM+'] Error: RNA ref_base does not equal to DNA ref_base!')
                sys.exit(1)
            if fd[3][:3] != '[M]' or fd[4][:3] != '[M]':
                continue
            bases = fd[3][4:4+int(fd[-1])].upper()
            quals = fd[4][4:4+int(fd[-1])]
            try:
                bases_new , quals_new = _DealBasequality(bases, quals, Phred_DNA)
            except IndexError:
                shell.eprint('['+PROGRAM+'] Error: please check the DNA single base file')
                sys.exit(1)
            cov_new = len(bases_new)
            if cov_new == 0:
                continue
            dna_info = [bases_new.count('A'), bases_new.count('C'), bases_new.count('G'), bases_new.count('T')]
            if Method == 'Bayesian':
                genotype , probability = _SNPvalue_Bayesian(bases_new, quals_new)
                try:
                    homo[fd[0]].update({fd[1]: [cov_new, ref_base, dna_info, genotype, '%.4f'%probability]})
                except KeyError:
                    homo[fd[0]] = {fd[1]: [cov_new, ref_base, dna_info, genotype, '%.4f'%probability]}
            elif Method == 'Binomial':
                pvalue = _SNVvalue_Binomial(dna_info, refbase)
                try:
                    homo[fd[0]].update({fd[1]: [cov_new, refbase, dna_info, pvalue]})
                except KeyError:
                    homo[fd[0]] = {fd[1]: [cov_new, refbase, dna_info, pvalue]}
            elif Method == 'Frequency':
                diff , probability = _SNPvalue_Frequency(bases_new, refbase)
                try:
                    homo[fd[0]].update({fd[1]: [cov_new, refbase, dna_info, diff, probability]})
                except KeyError:
                    homo[fd[0]] = {fd[1]: [cov_new, refbase, dna_info, diff, probability]}
            else:
                shell.eprint('['+PROGRAM+'] Error: unknown --method')
    logging.info('dna single base done')
    ## FDR for Binomial
    if Method == 'Binomial':
        barray , parray = [] , []
        for sca in homo.keys():
            for site in homo[sca].keys():
                barray.append(homo[sca][site]+[sca, site])
                parray.append(homo[sca][site][3])
        barray = sorted(barray, key=lambda x: x[0][3], reverse=True)
        parray = sorted(parray, reverse=True)
        parray_cor = _Fdr(parray)
        for i , line in enumerate(barray):
            barray[i][3] = _FormatP(line[3])
            parray_cor[i] = _FormatP(1) if parray_cor[i] > 1 else _FormatP(parray_cor[i])
            homo[line[4]][line[5]] = line[:4] + [parray_cor[i]]
    ## 
    readInfo_result = []
    for sca in sorted(rna_editing.keys()):
        for si in sorted(rna_editing[sca].keys(), key=lambda x:int(x)):
            try:
                cov, homobase, dna_info, homo_tag, posterior = homo[sca][si]
            except KeyError:
                continue
            if homobase == 'N':
                continue
            bases , quals = rna_editing[sca][si][1:3]
            editHash = {}
            for b in bases:
                if b != homobase:
                    try:
                        editHash[b] += 1
                    except KeyError:
                        editHash[b] = 1
            rna_info = []
            for i in range(len(T_BASE_ORDER)):
                rna_info.append(str(bases.count(T_BASE_ORDER[i])))
                dna_info[i] = str(dna_info[i])
            tag = []
            for idx , k in enumerate(sorted(editHash.keys(), key=lambda x:editHash[x], reverse=True)):
                tag.append([k, editHash[k]])
                if idx == 0:
                    best = [k, editHash[k]]
                    RNA_degree = '%.4f'%(editHash[k]/rna_editing[sca][si][-1])
            ###
            if best[0] == 'N' or best[1] < 2:
                continue
            if len(tag) > 2:
                num_array = []
                for tag_ele in tag:
                    num_array.append(tag_ele[1])
                num_array = sorted(num_array, reverse=True)
                if num_array[1]/num_array[0] > 0.01:
                    continue
            ###
            if Strand == '+':
                editingtype = homobase + '->' + best[0]
            elif Strand == '-':
                editingtype = D_BASE_COMPLIMENT[homobase] + '->' + D_BASE_COMPLIMENT[best[0]]
            elif Strand == 'Unknown':
                editingtype = homobase + '->' + best[0] + '|' + D_BASE_COMPLIMENT[homobase] + '->' + D_BASE_COMPLIMENT[best[0]]
            ###
            if Strand == '+' or Strand == '-':
                readInfo_result.append([sca, str(si), Strand, homobase, str(cov), ','.join(dna_info), homo_tag, str(posterior), 
                                    str(rna_editing[sca][si][-1]), ','.join(rna_info), editingtype, best[0], str(best[1]), str(RNA_degree)])
            else:
                readInfo_result.append([sca, str(si), '.', homobase, str(cov), ','.join(dna_info), homo_tag, str(posterior), 
                                    str(rna_editing[sca][si][-1]), ','.join(rna_info), editingtype, best[0], str(best[1]), str(RNA_degree)])
    del rna_editing
    logging.info('editing sites done')
    logging.info('add fdr')
    ##
    errorRate = 10**(-1*Qual_cutoff/10)     # Q=-10log(10)P
    for i in range(len(readInfo_result)):
        readInfo_result[i].append(_SNVvalue_Binomial((int(readInfo_result[i][12]),
                                                    int(readInfo_result[i][8])), 
                                                    errorRate, 
                                                    RNA=True))
    ##
    readInfo_result = sorted(readInfo_result, key=lambda x:x[-1], reverse=True)
    ##
    p = []
    for line in readInfo_result:
        p.append(line[-1])
    p_cor = _Fdr(p)
    ##
    for i in range(len(readInfo_result)):
        readInfo_result[i][-1] = _FormatP(readInfo_result[i][-1])
        readInfo_result[i].append(_FormatP(p_cor[i]))
    ##
    readInfo_result = sorted(readInfo_result, key=lambda x:(x[0],int(x[1])))
    for i in range(len(readInfo_result)):
        readInfo_result[i] = '\t'.join(readInfo_result[i][:11])+'\t'+'\t'.join(readInfo_result[i][-3:])+'\n'
    # write to file
    logging.info('write file')
    with gzip.open(Output+RNA_singleBase.split('/')[-1]+Suffix, 'wt') as f:
        f.write(head)
        f.writelines(readInfo_result)

    return 0

def filterRES(CRes, TDNAdepth=7, TRNAdepth=3, Edl=0.05, ExtremeLevel=0, Method='Bayesian', Bposterior=0.95, BPhete=0.05, Bfdr=0.05, FreNonRefCount=0, FreNonRefRatio=0, EdlP=0.05, Edep=3):
    '''

    '''
    # 
    f = shell.IsGzipFile(CRes)
    wr = [f.readline()]   # pass the headline
    for line in f:
        (sca, site, strand, gbase, covDNA, dna_BaseCount, lab1, lab2,
        covRNA, rna_BaseCount, editType, editDegree, pvalue, fdr) = line.strip().split()
        if int(covDNA)<TDNAdepth or int(covRNA)<TRNAdepth or float(editDegree)<Edl or float(fdr)>=EdlP:
            continue
        eH = {}
        eH['A'], eH['C'], eH['G'], eH['T'] = rna_BaseCount.split(',')
        if strand == '-':
            eB = D_BASE_COMPLIMENT[editType.split('->')[-1]]
        else:
            eB = editType.split('->')[-1]
        ##
        if int(eH[eB]) < Edep:
            continue
        ##
        if ExtremeLevel:
            if int(eH[gbase]) == 0 or editDegree == 1:
                continue
        ##
        if Method == 'Bayesian':
            str1 = set(list(lab1)) - set(gbase)
            if str1:
                continue
            if float(lab2) < Bposterior:
                continue
        elif Method == 'Binomial':
            if float(lab1) >= BPhete:
                continue
            if float(lba2) >= Bfdr:
                continue
        elif Method == 'Frequency':
            if int(lab1) > FreNonRefCount:
                continue
            if float(lab2) > FreNonRefRatio:
                continue
        else:
            shell.eprint('['+PROGRAM+'] Error: --method is not supported')
            sys.exit(1)
        wr.append(line)
    with gzip.open('.'.join(CRes.split('.')[:-1])+'.filter.gz', 'wt') as f:
        f.writelines(wr)

    return 0

def classifyRNA(BigTable, RNAbam, Readtype=3, Refined=1, RefinedDepth=1, Phred_DNA=33, Phred_RNA=33, Qual_cutoff=30, Trim=(6,6), ParalogousR=1):
    '''
        根据过滤后的RES位点从bam中得到符合条件的read并统计
    '''
    def _trim(cigartuples, t, start, seq, qual):
        '''
            Trim the reads' 5' and 3' end.
        '''
        newcigarlist , flag , tr = [] , 1 , 0   ## tr
        for operation , num in cigartuples:
            if flag:
                if operation in (0, 1, 4):  ## 0 for M, 1 for I, 4 for S
                    if num + tr > t[0]:
                        #shell.eprint(num, tr, t[0])
                        newcigarlist.append((operation, num + tr - t[0]))
                        if operation == 0 and len(t) == 2:
                            start += t[0] - tr
                        flag = 0
                    else:
                        tr += num
                        if operation == 0 and len(t) == 2:
                            start += num
                        #shell.eprint(newcigarlist, 'ok')
                elif operation in (2, 3):    ## 2 for D, 3 for N
                    if len(t) == 2:
                        start += num
                else:
                    shell.eprint("The cigar string may have unsupported char. Please check.")
            else:
                newcigarlist.append((operation, num))
            #shell.eprint(newcigarlist, t)
        if len(t) == 2:   
            return _trim(newcigarlist[::-1], [t[1]], start, seq[t[0]:-t[1]], qual[t[0]:-t[1]])
        elif len(t) == 1:
            ##
            if newcigarlist[0][0] == 1:
                newcigarlist[0] = (4, newcigarlist[0][1])
            elif newcigarlist[0][0] in (2, 3):
                newcigarlist = newcigarlist[1:]
            ##
            if newcigarlist[-1][0] == 1:
                newcigarlist[-1] = (4, newcigarlist[-1][1])
            elif newcigarlist[-1][0] in (2, 3):
                newcigarlist = newcigarlist[:-1]
            return newcigarlist[::-1], start , seq, qual
    def _offset(cigartuples, start):
        '''
            得到read在genome上的终止位置
        '''
        for operation , num in cigartuples:
            if operation == 0:
                start += num - 1
            elif operation in (2, 3):
                start += num + 1
            elif operation == 1:
                start += 1
            elif operation == 4:
                pass
            else:
                shell.eprint('The cigar string should be [MISDN]. Please check.')
        return start
    def _locateOnread(span, cigartuples):
        '''
            得到编辑位点具体在某条read上的位置
        '''
        offset , count = 0 , 0
        while span > 0:
            info = cigartuples.pop()
            count += 1
            if info[0] == 0:    # M
                if info[1] >= span:
                    offset += span
                    span = 0
                else:
                    offset += info[1]
                    span -= info[1]
            elif info[0] == 1:  # I
                offset += info[1]
            elif info[0] == 2:  # D
                if info[1] >= span:
                    return 0
                else:
                    span -= info[1]
            elif info[0] == 3:  # N
                if info >= span:
                    return 0
                else:
                    span -= info[1]
            elif info[0] == 4:  # S
                if count == 1:
                    offset += info[1]
                else:
                    return 0
            else:
                shell.eprint("The cigar string may have unsupported char. Please check.")                
        return offset
    ## Parameter check
    if RNAbam is None:
        shell.eprint('['+PROGRAM+'] Error: RNA bam file does not provided!')
        sys.exit(1)
    elif RNAbam.split('.')[-1] != 'bam':
        shell.eprint('['+PROGRAM+'] Error: do not support this format of alignment file! It needed to be compressed as bam if the offers is sam.')
        sys.exit(1)
    if BigTable is None:
        shell.eprint('['+PROGRAM+'] Error: the table of RNA candidate editing sites does not provided!')
        sys.exit(1)
    if Refined == 0 and ParalogousR == 0:
        shell.eprint('['+PROGRAM+'] Warning: do nothing for --ParalogousR 0 && --Refined 0')
        return 1
    ### Refined and ParalogousR
    if Refined:
        fw1 = gzip.open(BigTable+'.temp.gz', 'wt')
    if ParalogousR:
        fw2 = open(RNAbam.split('/')[-1]+'.editRead.fa', 'wt')
    ## Core
    ### Read Bigtable
    site , site_sort = {} , {}
    with gzip.open(BigTable, 'rt') as f:
        for line in f:
            line = line.strip()
            if line[:len('1.Chromosome_ID')] == '1.Chromosome_ID' and Refined:
                fw1.write('\t'.join([line, '15.Quadrant1/4_Read_Num(RNA)', '16.Quadrant2/3_Read_Num(RNA)', '17.Read_Type_Num(RNA)']) + '\n')
            else:
                fd = line.split()
                ## Ignore the strand ?
                try:
                    site[fd[0]].update({int(fd[1]): {'info': line, 'end': 0, 'mid': 0, 'readnum': 0, 'lastbeg': 0, 'lastend': 0}})
                except KeyError:
                    site[fd[0]] = {int(fd[1]): {'info': line, 'end': 0, 'mid': 0, 'readnum': 0, 'lastbeg': 0, 'lastend': 0}}
    ### sort the site
    for sca in site.keys():
        site_sort[sca] = sorted(site[sca].keys())
    ### Read RNAbam
    with pysam.AlignmentFile(RNAbam, 'rb') as f:
        for sca in sorted(site.keys()):
            for r in f.fetch(contig=sca):
                if len(site_sort[sca]) <= 0:
                    break
                if len(r.query_sequence) != len(r.qual):
                    continue
                if r.flag & 16:
                    Trim = Trim[::-1]
                newcigarlist , newstart , qread , qqual = _trim(list(r.cigartuples), Trim, r.reference_start+1, r.query_sequence, r.qual)
                newend = _offset(newcigarlist, newstart)
                l , cutmax , tailEndlen = deepcopy(site_sort[sca]) , max(Trim) , len(qread) // 4
                for pos in l:
                    if pos + cutmax + 1 < newstart and pos + cutmax < r.reference_start:
                        del site_sort[sca][0]
                    if pos > newend or l[-1] < newstart:
                        break
                    if newstart <= pos <= newend:
                        info = site[sca][pos]['info'].split()
                        editType = info[10].split('->')[1].upper() if info[2] == '+' else D_BASE_COMPLIMENT[info[10].split('->')[1].upper()]
                        posOnread = _locateOnread(pos-newstart+1, deepcopy(newcigarlist)[::-1])
                        if posOnread == 0 or qread[posOnread-1].upper() != editType or ord(qqual[posOnread-1]) - Phred_RNA < Qual_cutoff:
                            continue
                        if ParalogousR:
                            fw2.write('>{}\n{}\n'.format(r.query_name+'/'+('2' if r.flag & 128 else '1')+';'+r.reference_name+':'+str(pos)+':'+str(r.reference_start+1),
                                                        ''.join([D_BASE_COMPLIMENT[i] for i in r.query_sequence.upper()[::-1]]) if r.flag & 16 else r.query_sequence.upper()))
                        if Refined:
                            if tailEndlen <= posOnread <= tailEndlen * 3:
                                site[sca][pos]['mid'] += 1
                            else:
                                site[sca][pos]['end'] += 1
                            if site[sca][pos]['lastbeg'] != newstart or site[sca][pos]['lastend'] != newend:
                                site[sca][pos]['readnum'] += 1
                            site[sca][pos]['lastbeg'] = newstart
                            site[sca][pos]['lastend'] = newend
    if Refined:
        for sca in sorted(site.keys()):
            wr = []
            for pos in sorted(site[sca].keys()):
                if site[sca][pos]['mid'] < RefinedDepth or site[sca][pos]['readnum'] < Readtype:
                    continue
                wr.append('{}\t{}\t{}\t{}\n'.format(site[sca][pos]['info'], site[sca][pos]['end'], site[sca][pos]['mid'], site[sca][pos]['readnum']))
            fw1.writelines(wr)
        fw1.close()
        os.system('mv -f '+BigTable+'.temp.gz '+BigTable)

    if ParalogousR:
        fw2.close()

    return 0

def src_utils_pslScore_pslScore(*args):
    '''
        pslScore.pl - a reproduction of the C library calculations from src/lib/psl.c
        Source code of this script is from http://genome-source.cse.ucsc.edu/gitweb/?p=kent.git;a=blob_plain;f=src/utils/pslScore/pslScore.pl
        This is a python version.
    '''
    def _pslIsProtein(*subargs):
        (blockCount, strand, tStart, tEnd, tSize, tStarts, blockSizes) = subargs
        starts, sizes, lastBlock, answer = tStarts.split(','), blockSizes.split(','), blockCount - 1, 1
        if strand[0] == '+':
            answer = 3 if (tEnd == (int(starts[lastBlock]) + (3 * int(sizes[lastBlock])))) else 1
        elif strand[0] == '-':
            answer = 3 if (tStart == (tSize - (int(starts[lastBlock]) + (3 * int(sizes[lastBlock]))))) else 1
        return answer
    def _pslCalcMilliBad(*subargs):
        (sizeMul, qEnd, qStart, tEnd, tStart, qNumInsert, 
        tNumInsert, matches, repMatches, misMatches, isMrna) = subargs
        milliBad, qAliSize, tAliSize = 0, sizeMul*(qEnd-qStart), tEnd-tStart
        aliSize = tAliSize if tAliSize < qAliSize else qAliSize
        if aliSize <= 0:
            return milliBad
        sizeDif = qAliSize - tAliSize
        if sizeDif < 0:
            sizeDif = 0 if isMrna else (-sizeDif)
        insertFactor = (qNumInsert + tNumInsert) if isMrna == 0 else qNumInsert
        total = sizeMul * (matches + repMatches + misMatches)
        if total != 0:
            roundAwayFromZero = 3 * log(1+sizeDif)
            if roundAwayFromZero < 0:
                roundAwayFromZero = int(roundAwayFromZero - 0.5)
            else:
                roundAwayFromZero = int(roundAwayFromZero + 0.5)
            milliBad = (1000 * (misMatches * sizeMul + insertFactor + roundAwayFromZero)) / total
        return milliBad
    wr = []
    for fn in args:
        for line in shell.IsGzipFile(fn):
            (matches, misMatches, repMatches, nCount, qNumInsert, qBaseInsert, 
            tNumInsert, tBaseInsert, strand, qName, qSize, qStart, qEnd,
            tName, tSize, tStart, tEnd, blockCount, blockSizes, qStarts, 
            tStarts) = line.strip().split('\t')
            sizeMul = _pslIsProtein(int(blockCount), strand, int(tStart), int(tEnd), int(tSize), tStarts, blockSizes)
            pslScore = sizeMul * (int(matches) + (int(repMatches) >> 1)) - sizeMul * int(misMatches) - int(qNumInsert) - int(tNumInsert)
            milliBad = int(_pslCalcMilliBad(sizeMul, int(qEnd), int(qStart), int(tEnd), int(tStart), int(qNumInsert), int(tNumInsert), int(matches), int(repMatches), int(misMatches), 1))
            percentIdentity = 100.0 - milliBad * 0.1
            wr.append('{}\t{}\t{}\t{}:{}-{}\t{}\t{:.2f}\n'.format(tName, tStart, tEnd, qName, qStart, qEnd, pslScore, percentIdentity))
    with open(args[0]+'.score', 'w') as f:
        f.writelines(wr)

    return 0

def pslScore2editSite(PslScore, EditTable, BestHitRatio=0.5, JunctionCoordinate=None):
    '''
        Read .psl to form editSite
    '''
    # Parameter tickle
    def _assess_pslSocre(l):
        bestScore = int(sorted(l, key=lambda x: int(x[4]), reverse=True)[0][4])
        if bestScore == 0: 
            return 0
        for line in l:
            tID, tbeg, tend, qID, score, editChr, editPos, bamBeg = line
            if (tID == editChr and int(editPos) >= int(tbeg) and int(editPos) <= int(tend)) == False and int(score) / bestScore >= 0.95:
                return 0
        return 1
    if PslScore is None or EditTable is None:
        shell.eprint()
        sys.exit(1)
    junc = {}
    if JunctionCoordinate:
        for line in IsGzipFile(JunctionCoordinate):
            fd = line.strip().split()
            junc[fd[0]] = {'chr': fd[1], 'junction': fd[2]}
    psl , array = {} , []
    with open(PslScore, 'r') as f:
        for line in f:
            fd = line.strip().split()
            ## handle junc
            if fd[0] in junc:
                chrID = junc[fd[0]]['chr']
                l = []
                for b in junc[fd[0]]['junction'].split(','):
                    l += list(range(int(b.split('-')[0]), int(b.split('-')[1])+1))
                fd[0] , fd[1] , fd[2] = junc[fd[0]]['chr'] , l[fd[1]] , l[fd[2]-1]

            queryID = ':'.join(fd[3].split(':')[:-1])
            editChr , editPos , bamBeg = queryID.split(';')[1].split(':')
            try:
                if lastID != queryID or lastKey != (editChr, editPos):
                    try:
                        psl[lastKey] += _assess_pslSocre(array)
                    except KeyError:
                        psl[lastKey] = _assess_pslSocre(array)
                    lastID, lastKey = queryID, (editChr, editPos)
                    array = []
            except NameError:
                lastID, lastKey = queryID, (editChr, editPos)
            array.append([fd[0], fd[1], fd[2], queryID, fd[4], editChr, editPos, bamBeg])
    if array:
        try:
            psl[lastKey] += _assess_pslSocre(array)
        except KeyError:
            psl[lastKey] = _assess_pslSocre(array)
    wr = []
    with gzip.open(EditTable, 'rt') as f:
        for line in f:
            if line[:len('1.Chromosome_ID')] == '1.Chromosome_ID':
                wr.append(line.strip()+'\t18.Extracted_Read_Num(Blat)\t19.BestHit_Read_Num(Blat)\t20.BestHit_Read_Ratio(Blat)\n')
                continue
            fd = line.strip().split()
            base = {}
            for idx, num in enumerate(fd[9].split(',')):
                base[T_BASE_ORDER[idx]] = int(num)
            editBase = fd[10].split('->')[1] if fd[2] == '+' else D_BASE_COMPLIMENT[fd[10].split('->')[1]]
            totalEdit = base[editBase]
            try:
                rate = psl[(fd[0], fd[1])] / base[editBase]
            except KeyError:
                rate = 0 / base[editBase]
            if rate >= BestHitRatio:
                wr.append('{}\t{}\t{}\t{:.4f}\n'.format(line.strip(), base[editBase], psl[(fd[0], fd[1])], rate))
    with gzip.open(EditTable+'.temp', 'wt') as f:
        f.writelines(wr)
    os.system('mv -f '+EditTable+'.temp '+ EditTable)

    return 0

def print_help():
    print('''
Usage:   ires scanner [options] <genomefile> <DNAsbfile> <RNAsbfile> <RNAbamfile>

         <genomefile>     the genomefile file.
         <DNAsbfile>      the DNA single base file.
         <RNAsbfile>      the RNA single base file.
         <RNAbamfile>     the bam file to form RNAsbfile.
         Notes: --blat is not necessary for our program if you set --ParalogousR F.
                --strand should be provided correctly if you use strand-specific data.

Options: --outdir,-o      DIR    the output directory, '<pwd>/scanner/' is recommended. [./]
         --blat           PATH   the installed blat software. (force --ParalogousR T)
         --strand         STR    the strand infomation of the input RNAsbfile. (+ | - | unknown). [unknown]
         --coverage       INT    the minium depth in sb file. [2]
         --Method         STR    method for detecting SNPs. (Bayesian | Binomial | Frequency). [Bayesian]
         --rate           INT    rate of transitions over transversions of the genome. (force --Method Bayesian) [4]
         --HomoPrior      FLOAT  prior probability of homozygous genomic positions. (force --Method Bayesian) [0.99]
         --Ploid          INT    ploidy level, 1 for monoploid , 2 for diploid, 3 for triploid and so on. [2]
         --Qual-cutoff    INT    quality cutoff for BWA alignment result. [30]
         --Phred-DNA      INT    phred base quality for query quality of DNA. (33 | 64). [33]
         --Phred-RNA      INT    phred base quality for query quality of RNA. (33 | 64). [33]

         --DNA-depth      INT    homozygous DNA depth. [7]
         --RNA-depth      INT    minimum cutoff for RNA reads coverage. [3]
         --edl            FLOAT  minimum cutoff for editing level. [0.05]
         --extremelev     BOOL   exclude polymorphic sites with extreme degree of variation (100%) or not. [F]
         --Bposterior     FLOAT  minimun Bayesian Posterior Probability cutoff for corresponding genotype. (force --method Bayesian) [0.95]
         --Binophete      FLOAT  maximun P value cutoff of Binomial test for DNA Heterozygosis. (force --method Binomial) [0.05]
         --Binofdr        FLOAT  maximun FDR cutoff of Binomial test for DNA Heterozygosis. (force --method Binomial) [0.05]
         --Frequency-N    INT    maximun non-refference base count cutoff. (force --method Frequency) [0]
         --Frequency-R    FLOAT  maximun non-refference base ratio cutoff. (force --method Frequency) [0]
         --edit-Pvalue    FLOAT  cutoff of FDR RNA editing. [0.05]
         --edit-Depth     INT    minimum number of RNA reads supporting editing for a candidate editing site. [3]

         --trim           INT,INT  eg: 6,3 means trim 6 base from 5'end and 3 from 3'end in reads. [6,6]
         --goodread       INT      minimum number of unique RNA reads supporting editing for a candidate editing site. [3]
         --refined        BOOL     whether refined the number of RNA reads supporting candidate editing sites. [T]
         --refined-depth  INT      minimum number of RNA reads in the middle of its length supporting editing for a candidate editing site. [1]
                                   eg: from positions 23~68 of a 90-bp read
         --ParalogousR    BOOL     Remove candidate editing sites from those regions that are similar to other parts of the genome by BLAT alignment. [T]

         --help,-h               print help information
''')

def get_config():
    tobool , toint = {'T': True, 'F': False} , {'T': 1, 'F': 0}

    shortopts = 'ho:'
    longopts = ['help', 'outdir=', 'blat=', 'strand=', 'coverage=', 'Method=', 'rate=', 'HomoPrior=', 'Ploid=', 'Qual-cutoff=', 'Phred-DNA=', 'Phred-RNA=',
            'DNA-depth=', 'RNA-depth=', 'edl=', 'extremelev=', 'Bposterior=', 'Binophete=', 'Binofdr=', 'Frequency-N=', 'Frequency-R=', 'edit-Pvalue=', 'edit-Depth=',
            'trim=', 'goodread=', 'refined=', 'refined-depth=', 'ParalogousR=']

    try:
        optlist , args = getopt.getopt(sys.argv[1:], shortopts, longopts)
    except getopt.GetoptError as e:
        shell.eprint('['+PROGRAM+'] Error: '+str(e))
        sys.exit(2)
    
    if optlist == [] and args == []:
        print_help()
        sys.exit(0)
    
    try:
        genomefile , DNAsbfile , RNAsbfile , RNAbamfile = args
        genomefile = os.path.abspath(genomefile)
        DNAsbfile = os.path.abspath(DNAsbfile)
        RNAsbfile = os.path.abspath(RNAsbfile)
        RNAbamfile = os.path.abspath(RNAbamfile)
    except ValueError as e:
        print_help()
        sys.exit(1)

    config = {'candidateRES': {}, 'filterRES': {}, 'classifyRNA':{}, 'pslScore2editSite': {}}
    global OUTDIR , BLAT
    for opt , value in optlist:
        if opt in ('-h', '--help'):
            print_help()
            sys.exit(0)
        elif opt in ('-o', '--outdir'):
            OUTDIR = os.path.abspath(value) + '/'
        elif opt == '--blat':
            if os.path.isfile(value) == False:
                shell.eprint('['+PROGRAM+'] Error: --blat could be run')
                sys.exit(1)
            else:
                BLAT = os.path.abspath(value)
        elif opt == '--strand':
            if value in ('+', '-'):
                config['candidateRES']['Strand'] = value
            elif value.lower() == 'unknown':
                config['candidateRES']['Strand'] = 'Unknown'
            else:
                shell.eprint('['+PROGRAM+'] Error: --strand should be +, - or unknown')
                sys.exit(1)
        elif opt == '--coverage':
            try:
                config['candidateRES']['Coverage'] = int(value)
            except ValueError:
                shell.eprint('['+PROGRAM+'] Error: --coverage should be integer')
                sys.exit(1)
        elif opt == '--Method':
            if value in ('Bayesian', 'Binomial', 'Frequency'):
                config['candidateRES']['Method'] = value
                config['filterRES']['Method'] = value
            else:
                shell.eprint('['+PROGRAM+'] Error: --Method should be integer')
                sys.exit(1)
        elif opt == '--rate':
            try:
                config['candidateRES']['Ratio'] = int(value)
            except ValueError:
                shell.eprint('['+PROGRAM+'] Error: --Ratio should be integer')
                sys.exit(1)
        elif opt == '--HomoPrior':
            try:
                config['candidateRES']['HomoPrior'] = float(value)
            except ValueError:
                shell.eprint('['+PROGRAM+'] Error: --HomoPrior should be float')
                sys.exit(1)
        elif opt == '--Ploid':
            try:
                config['candidateRES']['Ploid'] = int(value)
            except ValueError:
                shell.eprint('['+PROGRAM+'] Error: --Ploid should be 1,2,3 and so on')
                sys.exit(1)
        elif opt == '--Qual-cutoff':
            try:
                config['candidateRES']['Qual_cutoff'] = int(value)
                config['classifyRNA']['Qual_cutoff'] = int(value)
            except ValueError:
                shell.eprint('['+PROGRAM+'] Error: --Qual-cutoff should be integer')
                sys.exit(1)
        elif opt == '--Phred-DNA':
            try:
                config['candidateRES']['Phred_DNA'] = int(value)
                config['classifyRNA']['Phred_DNA'] = int(value)
            except ValueError:
                shell.eprint('['+PROGRAM+'] Error: --Phred_DNA should be 33 or 64')
                sys.exit(1)
        elif opt == '--Phred-RNA':
            try:
                config['candidateRES']['Phred_RNA'] = int(value)
                config['classifyRNA']['Phred_RNA'] = int(value)
            except ValueError:
                shell.eprint('['+PROGRAM+'] Error: --Phred_RNA should be integer')
                sys.exit(1)
        ## fileterRES
        elif opt == '--DNA-depth':
            try:
                config['filterRES']['TDNAdepth'] = int(value)
            except ValueError:
                shell.eprint('['+PROGRAM+'] Error: --DNA-depth should be integer')
                sys.exit(1)
        elif opt == '--RNA-depth':
            try:
                config['filterRES']['TRNAdepth'] = int(value)
            except ValueError:
                shell.eprint('['+PROGRAM+'] Error: --RNA-depth should be integer')
                sys.exit(1)
        elif opt == '--edl':
            try:
                config['filterRES']['Edl'] = float(value)
            except ValueError:
                shell.eprint('['+PROGRAM+'] Error: --edl should be float')
                sys.exit(1)
        elif opt == '--extremelev':
            try:
                config['filterRES']['ExtremeLevel'] = toint[value]
            except KeyError:
                shell.eprint('['+PROGRAM+'] Error: --extremelev should be T or F')
                sys.exit(1)
        elif opt == '--Bposterior':
            try:
                config['filterRES']['Bposterior'] = float(value)
            except ValueError:
                shell.eprint('['+PROGRAM+'] Error: --Bposterior should be float')
                sys.exit(1)
        elif opt == '--Binophete':
            try:
                config['filterRES']['BPhete'] = float(value)
            except ValueError:
                shell.eprint('['+PROGRAM+'] Error: --Binophete should be float')
                sys.exit(1)
        elif opt == '--Binofdr':
            try:
                config['filterRES']['Bfdr'] = float(value)
            except ValueError:
                shell.eprint('['+PROGRAM+'] Error: --Binofdr should be float')
                sys.exit(1)
        elif opt == '--Frequency-N':
            try:
                config['filterRES']['FreNonRefCount'] = int(value)
            except ValueError:
                shell.eprint('['+PROGRAM+'] Error: --Frequency-N should be integer')
                sys.exit(1)
        elif opt == '--Frequency-R':
            try:
                config['filterRES']['FreNonRefRatio'] = float(value)
            except ValueError:
                shell.eprint('['+PROGRAM+'] Error: --Frequency-R should be float')
                sys.exit(1)
        elif opt == '--edit-Pvalue':
            try:
                config['filterRES']['EdlP'] = float(value)
            except ValueError:
                shell.eprint('['+PROGRAM+'] Error: --edit-Pvalue should be float')
                sys.exit(1)
        elif opt == '--edit-Depth':
            try:
                config['filterRES']['Edep'] = int(value)
            except ValueError:
                shell.eprint('['+PROGRAM+'] Error: --edit-Depth should be integer')
                sys.exit(1)
        elif opt == '--trim':
            try:
                t = [int(i) for i in value.split(',')]
                config['classifyRNA']['trim'] = t
            except ValueError:
                shell.eprint('['+PROGRAM+'] Error: --trim should be like int,int')
                sys.exit(1)
        elif opt == '--goodread':
            try:
                config['classifyRNA']['Readtype'] = int(value)
            except ValueError:
                shell.eprint('['+PROGRAM+'] Error: --goodread should be integer')
                sys.exit(1)
        elif opt == '--refined':
            try:
                config['classifyRNA']['Refined'] = toint[value]
            except KeyError:
                shell.eprint('['+PROGRAM+'] Error: --refined should be T or F')
        elif opt == '--refined-depth':
            try:
                config['classifyRNA']['RefinedDepth'] = int(value)
            except ValueError:
                shell.eprint('['+PROGRAM+'] Error: --refined-depth should be integer')
                sys.exit(1)
        elif opt == '--ParalogousR':
            try:
                config['classifyRNA']['ParalogousR'] = toint[value]
            except KeyError:
                shell.eprint('['+PROGRAM+'] Error: --ParalogousR should be T or F')
                sys.exit(1)
    try:
        flag = 1 if config['classifyRNA']['ParalogousR'] == 1 else 0
    except KeyError:
        flag = 1
    if flag == 1 and BLAT == '':
        r = subprocess.getstatusoutput('which blat')
        if r[0] == 1:
            shell.eprint('['+PROGRAM+'] Error: lack blat program')
            sys.exit(1)
        else:
            BLAT = r[1]

    if os.path.exists(OUTDIR) == False:
        try:
            os.makedirs(OUTDIR)
        except FileExistsError:
            shell.eprint('['+PROGRAM+'] Warning: outdir have existed, may have some conflict')
        except:
            shell.eprint('['+PROGRAM+'] Error: outdir could not be created, please check')
            sys.exit(1)

    os.chdir(OUTDIR)

    logging.basicConfig(level=logging.INFO,
                        filename=OUTDIR+PROGRAM+'.log',
                        filemode='w',
                        format='%(asctime)s : %(message)s',
                        datefmt='%m/%d/%Y %I:%M:%S %p')

    return config , genomefile , DNAsbfile , RNAsbfile , RNAbamfile

def main():
    config , genomefile , DNAsbfile , RNAsbfile , RNAbamfile = get_config()
    logging.info('Program start')
    logging.info('check input')
    # check input genomefile, DNAbamfile, RNAbamfile
    if not os.path.isfile(genomefile) or not os.path.isfile(DNAsbfile) or not os.path.isfile(RNAsbfile):
        shell.eprint('['+PROGRAM+'] Error: Please check the input file')
        sys.exit(1)
    logging.info('do identification of RNA editing sites')
    candidateRES(RNAsbfile, DNAsbfile, genomefile, **config['candidateRES'])
    logging.info('done')
    logging.info('do filter')
    filterRES(RNAsbfile.split('/')[-1]+'.homo.gz', **config['filterRES'])
    logging.info('done')
    logging.info('do extract reads')
    classifyRNA(RNAsbfile.split('/')[-1]+'.homo.filter.gz', RNAbamfile, **config['classifyRNA'])
    logging.info('done')
    logging.info('do blat')
    p = subprocess.Popen([BLAT, genomefile, RNAbamfile.split('/')[-1]+'.editRead.fa', RNAbamfile.split('/')[-1]+'.editRead.fa.psl', '-t=dna', '-q=dna', '-minIdentity=95', '-noHead'], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    p.wait()
    if p.returncode != 0:
        shell.eprint('['+PROGRAM+'] Error: blat does not complete, need check')
        shell.eprint(p.stderr.read().decode())
        sys.exit(1)
    logging.info('done')
    logging.info('do add blat result to table')
    src_utils_pslScore_pslScore(RNAbamfile.split('/')[-1]+'.editRead.fa.psl')
    pslScore2editSite(RNAbamfile.split('/')[-1]+'.editRead.fa.psl.score', RNAsbfile.split('/')[-1]+'.homo.filter.gz')
    logging.info('All things have been done. Have a good day!')

    return 0




if __name__ == '__main__':
    main()
