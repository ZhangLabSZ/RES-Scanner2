#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
# Author: Ji Li, liji1@genomics.cn

import os
import sys
import getopt
import subprocess

import shell

#
SCRIPT = os.path.basename(__file__)
PROGRAM = os.path.splitext(SCRIPT)[0]
BWA = ''
OUTDIR = ''
BT = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 'N': 'N'}
TYPE_ALL = ('AC', 'AG', 'AT', 'CA', 'CG', 'CT', 'GA', 'GC', 'GT', 'TA', 'TC', 'TG')
TYPE_GEO = ('AC', 'AG', 'AT', 'CG')
try:
    TT = {tp: BT[tp[0]]+BT[tp[1]] for tp in TYPE_ALL}
except:
    TT = {'AC': 'TG', 'AG': 'TC', 'AT': 'TA', 'CA': 'GT', 'CG': 'GC', 'CT': 'GA', 'GA': 'CT', 'GC': 'CG', 'GT': 'CA', 'TA': 'AT', 'TC': 'AG', 'TG': 'AC'}

def replace_genome(filename, tp, outfile1, outfile2, BWA):
    f , out1 , out2 = shell.IsGzipFile(filename) , open(outfile1, 'wt') , open(outfile2, 'wt')
    for line in f:
        if line[0] == '>':
            out1.write(line)
            out2.write(line)
        else:
            out1.write(line.upper().replace(tp[0], tp[1]))
            out2.write(line.upper().replace(BT[tp[0]], BT[tp[1]]))
    for fd in [f, out1, out2]:
        fd.close()
    # index and rm
    p1 = subprocess.Popen([BWA, 'index', outfile1], stderr=subprocess.PIPE, stdout=subprocess.PIPE)
    p2 = subprocess.Popen([BWA, 'index', outfile2], stderr=subprocess.PIPE, stdout=subprocess.PIPE)
    p1.wait()
    p2.wait()
    if p1.returncode != 0 or p2.returncode != 0:
        shell.eprint('['+PROGRAM+'] Error: bwa index failed')
        sys.exit(1)
    os.system('rm -f '+outfile1)
    os.system('rm -f '+outfile2)

    return 0

def print_help():
    print('''
Usage:   ires replacegeo [options] <genomefile>

Options: --outdir,-o   DIR   default is '<pwd>/regeo/'
         --bwa         PATH  the absolute path of hisat2 program or system env. (force)
''')

def get_config():
    try:
        optlist , args = getopt.getopt(sys.argv[1:], 'ho:', ['help', 'outdir=', 'bwa='])
    except getopt.GetoptError as e:
        shell.eprint('['+PROGRAM+'] Error: '+str(e))
        sys.exit(2)
    
    if optlist == [] and args == []:
        print_help()
        sys.exit(0)

    global BWA, OUTDIR
    for opt , value in optlist:
        if opt in ('-h', '--help'):
            print_help()
            sys.exit(1)
        elif opt in ('-o', '--outdir'):
            OUTDIR = os.path.abspath(value) + '/'
        elif opt == '--bwa':
            if os.path.exists(value):
                BWA = os.path.abspath(value)
            else:
                shell.eprint('['+PROGRAM+'] Error: bwa path does not exist')
                sys.exit(1)
    try:
        genomefile = args[0]
        if os.path.isfile(genomefile):
            genomefile = os.path.abspath(genomefile)
        else:
            shell.eprint('['+PROGRAM+'] Error: genome file does not exist')
            sys.exit(1)
    except ValueError as e:
        shell.eprint('['+PROGRAM+'] Error: '+str(e))
        sys.exit(1)
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
    if OUTDIR == '':
        OUTDIR = os.getcwd()+'/regeo'
        try:
            os.makedirs(OUTDIR)
        except FileExistsError:
            pass
        except:
            shell.eprint('['+PROGRAM+'] Error: outdir could not be created, please check')
            sys.exit(1)
    os.chdir(OUTDIR)

    return genomefile

def main():
    genomefile = get_config()
    for tp in TYPE_GEO:
        replace_genome(genomefile, tp, tp+'.fa', TT[tp]+'.fa', BWA)

    return 0




if __name__ == '__main__':
    main()
