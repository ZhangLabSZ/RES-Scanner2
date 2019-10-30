#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
# Author: Ji Li, liji1@genomics.cn

import os
import sys
import gzip
import pysam
import getopt

import shell

# program path
SCRIPT = os.path.basename(__file__)
PROGRAM = os.path.splitext(SCRIPT)[0]
#PATH = os.path.abspath(__file__)
#DIR = os.path.dirname(PATH)

D_BASE_COMPLIMENT = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 'N': 'N'}
T_BASE_ORDER = ('A', 'C', 'G', 'T')

def sam2base(genomefile, sortbamfile, trim=(0,0), output='', suffix='.sb.gz'):
    '''
        Bam to sb file. sb -> singlebase, for the future, sbz for the sb specific compressed file,
        sbi for the sb index file.
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
                    shell.eprint('['+PROGRAM+'] Error: the cigar string may have unsupported char, should only have [MISDN]')
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
    def _transfer_cigar(cigartuples, start, query, qual):
        '''
           Judge the extend end of scaffold by read if out, and return  
        '''
        #shell.eprint(cigartuples)
        d = {}  # record handle of each cigar
        offset = 0  # the offset on read
        i = start  # i for Delete position in genome
        #shell.eprint(query)
        for operation, num in cigartuples:
            #shell.eprint(offset)
            #shell.eprint(query[offset: offset+num])
            if operation == 0:  ## M
                try:
                    d['M'].append((start, start+num, query[offset: offset+num], qual[offset: offset+num]))
                except KeyError:
                    d['M'] = [(start, start+num, query[offset: offset+num], qual[offset: offset+num])]
                #shell.eprint(d)
                start += num - 1
                i += num
                offset += num
            elif operation == 1:    ## I
                try:
                    d['I'].append((start, query[offset: offset+num], qual[offset: offset+num]))
                except KeyError:
                    d['I'] = [(start, query[offset: offset+num], qual[offset: offset+num])]
                start += 1
                offset += num
            elif operation == 2:    ## D
                try:
                    d['D'].append((i, num))
                except KeyError:
                    d['D'] = [(i, num)]
                start += num + 1
                i += num
            elif operation == 3:    ## N
                start += num + 1
                i += num
                offset += num
            elif operation == 4:    ## S
                offset += num
            else:
                shell.eprint('['+PROGRAM+'] Error: the cigar string may have unsupported char, should only have [MISDN]')
        return start, d
    def _update(seq, dw, **kw):
        if 'M' in kw:
            for line in kw['M']:
                #shell.eprint(line)
                for j , i in enumerate(range(line[0], line[1])):
                    if i not in dw:
                        dw[i] = {'M': [line[2][j], line[3][j]]}
                    else:
                        try:
                            dw[i]['M'][0] += line[2][j]
                            dw[i]['M'][1] += line[3][j]
                        except KeyError:
                            #shell.eprint(j)
                            dw[i].update({'M': [line[2][j], line[3][j]]})
        if 'I' in kw:
            for line in kw['I']:
                if line[0] not in dw:
                    dw[line[0]] = {'I': [line[1], line[2]]}
                else:
                    try:
                        dw[line[0]]['I'][0] += ',' + line[1]
                        dw[line[0]]['I'][1] += ',' + line[2]
                    except KeyError:
                        dw[line[0]].update({'I': [line[1], line[2]]})
        if 'D' in kw:
            for line in kw['D']:
                st = ''
                for i in range(line[1]):
                    st += seq[line[0]+i-1]
                if line[0] not in dw:
                    dw[line[0]] = {'D': [st, '0'*line[1]]}
                else:
                    try:
                        dw[line[0]]['D'][0] += ',' + st
                        dw[line[0]]['D'][1] += ',' + '0'*line[1]
                    except KeyError:
                        dw[line[0]].update({'D': [st, '0'*line[1]]})
        return dw
    def _write2gzipfile(dw, fw):
        wr = []
        for key in sorted(dw.keys()):
            if 'M' in dw[key]:
                if 'I' in dw[key]:
                    if 'D' in dw[key]:
                        wr.append(sca+'\t'+str(key)+'\t'+seq[key-1]+'\t'
                                +'[M]:'+dw[key]['M'][0]
                                +';[I]:'+dw[key]['I'][0]
                                +';[D]:'+dw[key]['D'][0]+'\t'
                                +'[M]:'+dw[key]['M'][1]
                                +';[I]:'+dw[key]['I'][1]
                                +';[D]:'+dw[key]['D'][1]+'\t'
                                +str(len(dw[key]['M'][0]))+'\n')
                    else:
                        wr.append(sca+'\t'+str(key)+'\t'+seq[key-1]+'\t'
                                +'[M]:'+dw[key]['M'][0]
                                +';[I]:'+dw[key]['I'][0]+'\t'
                                +'[M]:'+dw[key]['M'][1]
                                +';[I]:'+dw[key]['I'][1]+'\t'
                                +str(len(dw[key]['M'][0]))+'\n')
                elif 'D' in dw[key]:
                    wr.append(sca+'\t'+str(key)+'\t'+seq[key-1]+'\t'
                            +'[M]:'+dw[key]['M'][0]
                            +';[D]:'+dw[key]['D'][0]+'\t'
                            +'[M]:'+dw[key]['M'][1]
                            +';[D]:'+dw[key]['D'][1]+'\t'
                            +str(len(dw[key]['M'][0]))+'\n')
                else:
                    wr.append(sca+'\t'+str(key)+'\t'+seq[key-1]+'\t'
                            +'[M]:'+dw[key]['M'][0]+'\t'
                            +'[M]:'+dw[key]['M'][1]+'\t'
                            +str(len(dw[key]['M'][0]))+'\n')
        fw.write(''.join(wr))
        return {}
    if not sortbamfile.endswith('.bam'):
        shell.eprint('['+PROGRAM+'] Error: input should be sorted and indexed bam')
        sys.exit(1)
    if trim[0] < 0 or trim[1] < 0:
        shell.eprint('Are you serious ? --trim should not be negative integer.')
        sys.exit(255)
    if output == '':
        output = sortbamfile
    # Get a chromosome or scaffold from genome then do with bam
    ## If scaffold has no mapping reads, could it raise error?
    fw = gzip.open(output + suffix, 'wt')
    with pysam.AlignmentFile(sortbamfile, 'rb') as f:
        for sca, seq, seqlen in shell.Fa2Geno(genomefile):
            dw , newend = {} , 1
            try:
                for read in f.fetch(contig=sca):
                    if 1 < newend < read.reference_start + 1:
                        dw = _write2gzipfile(dw, fw)
                    if len(read.query_sequence) != len(read.qual):
                        continue
                    #shell.eprint(read.query_sequence)
                    if trim == (0,0):
                        newstart , d = _transfer_cigar(read.cigartuples, read.reference_start + 1, read.query_sequence, read.qual)
                    else:
                        if read.flag & 16:
                            trim = trim[::-1]
                        newstart , d =  _transfer_cigar(*_trim(list(read.cigartuples), trim, read.reference_start + 1, read.query_sequence, read.qual))
                    #
                    if seqlen < newstart:
                        continue
                    if newstart > newend:
                        newend = newstart
                    dw = _update(seq, dw, **d)
                _write2gzipfile(dw, fw)
            except ValueError as e:
                shell.eprint('['+PROGRAM+'] Error: input should be sorted and indexed bam')
                os.system('rm -f '+output+suffix)
                sys.exit(1)
    fw.close()

    return output + suffix

def stat(genomefile, sbf):
    f = shell.IsGzipFile(sbf)
    bases , cov , d , sca = 0 , 0 , {} , set()
    for line in f:
        fd = line.strip().split()
        sca.add(fd[0])
        if int(fd[-1]) > 0:
            cov += 1
        bases += int(fd[-1])
        try:
            d[int(fd[-1])] += 1
        except KeyError:
            d[int(fd[-1])] = 1
    f.close()
    f = shell.IsGzipFile(genomefile)
    st , lenth = '' , 0
    for line in f:
        if line[0] == '>':
            if st and key in sca:
                lenth += len(st)
            key = line.strip().split()[0][1:]
            st = ''
        else:
            st += line.strip()
    if st and key in sca:
        lenth += len(st)
    try:
        peak_depth = sorted(d.items(), key=lambda x: x[1], reverse=True)[0][0]
    except IndexError:
        peak_depth = 0
    ave_depth = 0 if cov == 0 else '%.2f'%(bases / cov)
    cov_rate = '%.2f'%(cov / lenth * 100)
    with open(sbf+'.stat', 'w') as f:
        f.write('#total_bases\t#covered_bases\tcoveage(%)\tcoverage(X)\tpeak_coverage(X)\n')
        f.write('{}\t{}\t{}\t{}\t{}\n'.format(lenth, cov, cov_rate, ave_depth, peak_depth))

    return 0

def print_help():
    print('''
Usage:   ires sam2base [options] <genomefile> <bamfile-with-index>

Options: --trim     INT,INT  eg: 6,3 means trim 6 base from 5'end and 3 from 3'end in reads. [0,0]
                                 0,0 is recommended for DNA; 6,6 is recommended for RNA
         --prefix   STR      the prefix of output file. [same as bam name]
         --stat     BOOL     stat single base file. [T]
         --help,-h           print help information.
''')

def get_config():
    try:
        optlist , args = getopt.getopt(sys.argv[1:], 'hp:t:', ['help', 'trim=', 'prefix='])
    except getopt.GetoptError as e:
        shell.eprint('['+PROGRAM+'] Error: '+str(e))
        sys.exit(2)

    if optlist == [] and args == []:
        print_help()
        sys.exit(0)

    config = {'sam2base': {}}
    for opt , value in optlist:
        if opt in ('-h', '--help'):
            print_help()
            sys.exit(0)
        elif opt in ('-t', '--trim'):
            try:
                config['sam2base']['trim'] = (int(value.split(',')[0]), int(value.split(',')[1]))
            except ValueError:
                shell.eprint('['+PROGRAM+'] Error: --trim parameter should be integer')
                sys.exit(1)
        elif opt in ('p', '--prefix'):
            if value.endswith('.'):
                shell.eprint('['+PROGRAM+'] Error: --prefix parameter should not end with \'.\'')
                sys.exit(1)
            else:
                config['sam2base']['output'] = value
        else:
            assert False, 'unhandled option'
    try:
        genomefile , bamfile = args
    except ValueError as e:
        shell.eprint('['+PROGRAM+'] Error: only two input files should be provided')
        sys.exit(1)

    return config , genomefile , bamfile

def main():
    config , genomefile , bamfile = get_config()
    sbfile = sam2base(genomefile, bamfile, **config['sam2base'])
    try:
        stat(genomefile, sbfile)
    except ZeroDivisionError:
        pass

    return 0




if __name__ == '__main__':
    main()
