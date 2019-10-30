#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
# Author: Ji Li, liji1@genomics.cn

import os
import sys
import gzip
import json
import pysam
import getopt
import logging
import subprocess

import shell

SCRIPT = os.path.basename(__file__)
PROGRAM = os.path.splitext(SCRIPT)[0]

SOAPNUKE , JAVA , PILON , BWA = '' , '' , '' , ''

#
OUTDIR = os.getcwd() + '/'
# 0 means neither 2 means soapnuke yes, 4 means pilon yes.
FLAG = 0
# to solve different number of input files.
DBAM , RBAM , DFQ1 , DFQ2 , RFQ1 , RFQ2 = '' , '' , '' , '' , '' , ''
PDFQ1 , PDFQ2 , PRFQ1 , PRFQ2 = None , None , None , None

D_BASE_COMPLIMENT = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 'N': 'N'}
T_BASE_ORDER = ('A', 'C', 'G', 'T')

def soapnuke(tp, fq1, fq2, phred, config):
    if phred == 33:
        args = ['-Q', '2']
    elif phred == 64:
        args = ['-Q', '1']
    else:
        shell.eprint('['+PROGRAM+'] Error: phred error')
        sys.exit(1)
    if fq2 == '':
        args += ['-1', fq1, '-o', 'soapnuke/'+tp, '-C', os.path.basename(fq1)+'.clean.fq.gz']
        fq1 = 'soapnuke/'+tp+'/'+os.path.basename(fq1)+'.clean.fq.gz'
    else:
        args += ['-1', fq1, '-2', fq2, '-o', 'soapnuke/'+tp, '-C', os.path.basename(fq1)+'_1.clean.fq.gz', '-D', os.path.basename(fq2)+'_2.clean.fq.gz']
        fq1 , fq2 = 'soapnuke/'+tp+'/'+os.path.basename(fq1)+'_1.clean.fq.gz' , 'soapnuke/'+tp+'/'+os.path.basename(fq2)+'_2.clean.fq.gz'
    for k,v in config['soapnuke']['filter'][tp].items():
        if v != '':
            args += [k, v]
    p = subprocess.Popen([SOAPNUKE, 'filter']+args+['-G', '-5', '1'], stdout=open('log', 'w'), stderr=subprocess.STDOUT)
    p.wait()
    if p.returncode != 0:
        shell.eprint('['+PROGRAM+'] Error: soapnuke run error')
        shell.eprint(''.join(open('log', 'r').read()))
        sys.exit(1)

    return fq1, fq2 , 33

def bwaaln(tp, genomefile, fq1, fq2, phred, outdir, config):
    if phred == 33:
        args = []
    elif phred == 64:
        args = ['-I']
    else:
        shell.eprint('['+PROGRAM+'] Error: phred error')
        sys.exit(1)
    for k,v in config['bwa']['aln'][tp].items():
        if v != '':
            args += [k, v]
    p1 = subprocess.Popen([BWA, 'aln']+args+[genomefile, fq1, '-f', outdir+os.path.basename(fq1)+'_1.sai'], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    if fq2 != '':
        p2 = subprocess.Popen([BWA, 'aln']+args+[genomefile, fq2, '-f', outdir+os.path.basename(fq2)+'_2.sai'], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        p2.wait()
    p1.wait()
    if p1.returncode != 0 and p2.returncode != 0:
        shell.eprint('['+PROGRAM+'] Error: bwa aln error')
        sys.exit(1)
    args = []
    for k,v in config['bwa']['sampe'][tp].items():
        args += [k, v]
    if fq2 != '':
        p = subprocess.Popen([BWA, 'sampe']+args+[genomefile, outdir+os.path.basename(fq1)+'_1.sai', outdir+os.path.basename(fq2)+'_2.sai', fq1, fq2], stdout=subprocess.PIPE, stderr=subprocess.DEVNULL)
    else:
        p = subprocess.Popen([BWA, 'samse', genomefile, outdir+os.path.basename(fq1)+'_1.sai', fq1], stdout=subprocess.PIPE, stderr=subprocess.DEVNULL)
    infile = pysam.AlignmentFile(p.stdout, 'r')
    oufile = pysam.AlignmentFile(outdir+tp+'_aln.bam', 'wb', template=infile)
    for r in infile:
        oufile.write(r)
    p.wait()
    oufile.close()
    infile.close()
    os.system('rm -f '+outdir+os.path.basename(fq1)+'_1.sai')
    if fq2:
        os.system('rm -f '+outdir+os.path.basename(fq2)+'_2.sai')

    return outdir+tp+'_aln.bam'

def bestuniqbam(bamfile, mapQ=20, DNA=None, RNA=None, SS=False, RmDup=True, Uniq=True, outdir='./'):
    '''
        Filter bam to best uniq bam.
        ##eg: bestuniqbam(bamfile, RNA-mapQ=20, RNA=True, SS=True, RmDup=True, Uniq=True, outdir='./scanner')
    '''
    def _flag_split(g):
        '''
            According to sam format, 115 = 64 + 32 + 16 + 2 + 1.
        '''
        return set([ g & 2**i for i in range(0,12)])
    # check bamfile
    if bamfile.split('.')[-1] == 'bam':
        file_type_tag = 'rb'
    elif bamfile.split('.')[-1] == 'sam':
        file_type_tag = 'r'
    else:
        shell.eprint('['+PROGRAM+'] Error: the input bamfile/samfile should be *.bam or *.sam.')
        sys.exit(1)
    # parameter tickle
    if (RNA is None and DNA is None) or (RNA != None and DNA != None):
        shell.eprint('['+PROGRAM+'] Error: --DNA or --RNA is needed!')
        sys.exit(1)
    if DNA == True and SS == True:
        shell.eprint('['+PROGRAM+'] Warning: DNA do not have --ss, but it\'s ok to run this program.')
    if os.path.exists(outdir) == False:
        try:
            os.makedirs(outdir)
        except:
            sys.exit(1)
    # core
    try:
        f = pysam.AlignmentFile(bamfile, file_type_tag)
    except:
        shell.eprint('['+PROGRAM+'] Error: check the bam file '+os.path.basename(bamfile)+' please!')
        sys.exit(1)
    ## DNA and RNA's bamfile name should be different by yourself.
    if SS:
        fn1 = outdir+'.'.join(os.path.basename(bamfile).split('.')[:-1])+'.negative.bam'
        fn2 = outdir+'.'.join(os.path.basename(bamfile).split('.')[:-1])+'.positive.bam'
        fw1 = pysam.AlignmentFile(fn1, 'wb', template=f)
        fw2 = pysam.AlignmentFile(fn2, 'wb', template=f)
    else:
        fn3 = outdir+'.'.join(os.path.basename(bamfile).split('.')[:-1])+'.best.bam'
        fw3 = pysam.AlignmentFile(fn3, 'wb', template=f)
    for line in f.fetch(until_eof=True):
        if int(line.mapping_quality) < mapQ:
            continue
        flag = int(line.flag)
        tags = line.get_tags()
        if 1024 & flag and RmDup:
            continue
        if RNA:
            if Uniq and (('XT', 'U') not in tags or ('X0', 1) not in tags or ('X1', 0) not in tags):
                continue
            if 256 & flag:
                continue
            if flag in (67, 131, 115, 179):
                continue
            if 32 & flag and 16 & flag:
                continue
            if SS:
                if 64 & flag and 32 & flag and line.template_length >= 0:
                    fw1.write(line)
                elif 128 & flag and 16 & flag and line.template_length <= 0:
                    fw1.write(line)
                elif 64 & flag and 16 & flag and line.template_length <= 0:
                    fw2.write(line)
                elif 128 & flag and 32 & flag and line.template_length >= 0:
                    fw2.write(line)
                elif 64 & flag and 16 & flag != 16 and 32 & flag != 32:
                    fw1.write(line)
                elif 128 & flag and 16 & flag != 16 and 32 & flag != 32:
                    fw2.write(line)
                elif flag == 16:
                    fw2.write(line)
                elif flag == 0:
                    fw1.write(line)
            else:
                fw3.write(line)
        elif DNA:
            fw3.write(line)
        else:
            pass
    if SS:
        fw1.close()
        fw2.close()
        return fn1 , fn2
    else:
        fw3.close()
        return fn3

def bamsortindex(bamfile):
    prefix = '.'.join(bamfile.split('.')[:-1])
    try:
        pysam.sort('-o', prefix+'.sorted.bam', bamfile)
    except:
        shell.eprint('['+PROGRAM+'] Error: '+bamfile+' sort error')
        sys.exit(1)
    try:
        pysam.index(prefix+'.sorted.bam')
    except:
        shell.eprint('['+PROGRAM+'] Error: '+bamfile+' index error')
    os.system('rm -f '+bamfile)
    return prefix+'.sorted.bam'

def print_help():
    print('''
Usage:   ires preprocess [options] <genomefile>

         <genomefile>   the genome file   
         Notes: fastq file or bam should be provided only one.
                Our program will not check the parameters in the config file, so make sure that 
                your parameters will be processed correctly when using these sofeware alone.

Options: --outdir,-o  DIR      the output directory, '<pwd>/preprocess/' is recommended. [./]
         --DNA-fq1    FILE     DNA fastq file, read1 in paired-end or single-end. (force)
         --DNA-fq2    FILE     DNA fastq file, read2 in paired-end.
         --RNA-fq1    FILE     RNA fastq file, read1 in paired-end or single-end. (force)
         --RNA-fq2    FILE     RNA fastq file, read2 in paired-end.
         --DNA-bam    FILE     bam file, formed by bwa aln (force)
         --RNA-bam    FILE     bam file, formed by bwa aln (force)
         --bwa        PATH     absolute path of bwa. (force)
         --soapnuke   PATH     absolute path of soapnuke, will be added -G by force.
         --java       PATH     absolute path of java
         --pilon      PATH     absolute path of pilon. (force --java)
         --config     FILE     config of three software above, default file is in testdata/
         --DNA-mapQ   INT      the minium mapping quality of best uniq bam. [20]
         --RNA-mapQ   INT      the minium mapping quality of best uniq bam. [20]
         --rmdup      BOOL     remove PCR duplication for RNA. [T]
         --uniq       BOOL     get the uniq mapping for RNA. [T]
         --ss         BOOL     if RNA library strategy is strand-specific. [F]

         --help,-h            print help information
''')

def get_config():
    shortopts = 'ho:'
    longopts = ['help', 'outdir=', 'DNA-fq1=', 'DNA-fq2=', 'RNA-fq1=', 'RNA-fq2=', 'DNA-bam=', 'RNA-bam=', 
                'bwa=', 'soapnuke=', 'java=', 'pilon=', 'config=', 'DNA-mapQ=', 'RNA-mapQ=', 'rmdup=', 'uniq=', 
                'ss=', 'DNA-I=', 'RNA-I=']
    try:
        optlist , args = getopt.getopt(sys.argv[1:], shortopts, longopts)
    except getopt.GetoptError as e:
        shell.eprint('['+PROGRAM+'] Error: '+str(e))
        sys.exit(2)
    
    if optlist == [] and args == []:
        print_help()
        sys.exit(0)
 
    config = {'bestuniqbam': {'DNA': {}, 'RNA': {}}}
    tobool = {'T': True, 'F': False}
    global DBAM , DFQ1 , DFQ2 , RBAM , RFQ1 , RFQ2
    global PDFQ1 , PDFQ2 , PRFQ1 , PRFQ2
    global BWA , SOAPNUKE , JAVA , PILON
    global FLAG , OUTDIR
    for opt , value in optlist:
        if opt in ('-h', '--help'):
            print_help()
            sys.exit(0)
        elif opt in ('-o', '--outdir'):
            OUTDIR = os.path.abspath(value) + '/'
        elif opt == '--DNA-fq1':
            if os.path.exists(value):
                DFQ1 = os.path.abspath(value)
            else:
                shell.eprint('['+PROGRAM+'] Error: dna fq1 file does not exist')
                sys.exit(1)
        elif opt == '--DNA-fq2':
            if os.path.exists(value):
                DFQ2 = os.path.abspath(value)
            else:
                shell.eprint('['+PROGRAM+'] Error: dna fq2 file does not exist')
                sys.exit(1)
        elif opt == '--RNA-fq1':
            if os.path.exists(value):
                RFQ1 = os.path.abspath(value)
            else:
                shell.eprint('['+PROGRAM+'] Error: rna fq1 file does not exist')
                sys.exit(1)
        elif opt == '--RNA-fq2':
            if os.path.exists(value):
                RFQ2 = os.path.abspath(value)
            else:
                shell.eprint('['+PROGRAM+'] Error: rna fq2 file does not exist')
                sys.exit(1)
        elif opt == '--DNA-bam':
            if os.path.exists(value):
                DBAM = os.path.abspath(value)
            else:
                shell.eprint('['+PROGRAM+'] Error: dna bam file does not exist')
                sys.exit(1)
        elif opt == '--RNA-bam':
            if os.path.exists(value):
                RBAM = os.path.abspath(value)
            else:
                shell.eprint('['+PROGRAM+'] Error: rna bam file does not exist')
                sys.exit(1)
        elif opt == '--bwa':
            if os.path.exists(value):
                BWA = os.path.abspath(value)
            else:
                shell.eprint('['+PROGRAM+'] Error: bwa path does not exist')
                sys.exit(1)
        elif opt == '--soapnuke':
            if os.path.exists(value):
                SOAPNUKE = os.path.abspath(value)
            else:
                shell.eprint('['+PROGRAM+'] Error: soapnuke path does not exist')
                sys.exit(1)
        elif opt == '--java':
            if os.path.exists(value):
                JAVA = os.path.abspath(value)
            else:
                shell.eprint('['+PROGRAM+'] Error: java path does not exist')
                sys.exit(1)
        elif opt == '--pilon':
            if os.path.exists(value):
                PILON = os.path.abspath(value)
            else:
                shell.eprint('['+PROGRAM+'] Error: pilon path does not exist')
                sys.exit(1)
        elif opt == '--config':
            if os.path.exists(value):
                with open(value, 'r') as f:
                    config.update(json.load(f))
            else:
                shell.eprint('['+PROGRAM+'] Error: config file does not exist')
                sys.exit(1)
            if len(set(['-C', '-D', '-o', '-1', '-2']) & config['soapnuke']['filter'].keys()) > 0:
                shell.eprint('['+PROGRAM+'] Error: soapnuke parameters -C -D -o -1 -2 do not needed to provided')
                sys.exit(1)
        elif opt == '--DNA-mapQ':
            try:
                config['bestuniqbam']['DNA']['mapQ'] = int(value)
            except ValueError:
                shell.eprint('['+PROGRAM+'] Error: --DNA-mapQ should be integer')
                sys.exit(1)
        elif opt == '--RNA-mapQ':
            try:
                config['bestuniqbam']['RNA']['mapQ'] = int(value)
            except ValueError:
                shell.eprint('['+PROGRAM+'] Error: --RNA-mapQ should be integer')
                sys.exit(1)
        elif opt == '--rmdup':
            try:
                config['bestuniqbam']['RNA']['RmDup'] = tobool[value]
            except KeyError:
                shell.eprint('['+PROGRAM+'] Error: --rmdup should be T or F')
                sys.exit(1)
        elif opt == '--uniq':
            try:
                config['bestuniqbam']['RNA']['Uniq'] = tobool[value]
            except KeyError:
                shell.eprint('['+PROGRAM+'] Error: --uniq should be T or F')
                sys.exit(1)
        elif opt == '--ss':
            try:
                config['bestuniqbam']['RNA']['SS'] = tobool[value]
            except KeyError:
                shell.eprint('['+PROGRAM+'] Error: --ss should be T or F')
                sys.exit(1)
        else:
            assert False, 'unhandled option'

    if BWA == '':
        shell.eprint('['+PROGRAM+'] Error: bwa is necessary')
        sys.exit(1)

    if PILON != '' and JAVA == '':
        shell.eprint('['+PROGRAM+'] Error: pilon process needs java')
        sys.exit(1)
    
    if os.path.exists(OUTDIR) == False:
        try:
            os.makedirs(OUTDIR)
        except FileExistsError:
            shell.eprint('['+PROGRAM+'] Warning: outdir have existed, may have some conflict')
        except:
            shell.eprint('['+PROGRAM+'] Error: outdir could not be created, please check')
            sys.exit(1)

    try:
        genomefile = args[0]
        genomefile = os.path.abspath(genomefile)
        try:
            os.symlink(genomefile, OUTDIR+os.path.basename(genomefile))
        except FileExistsError:
            pass
        except:
            shell.eprint('['+PROGRAM+'] Error: could not create soft link of genomefile')
            sys.exit(1)
        genomefile = os.path.basename(genomefile)
    except ValueError:
        shell.eprint('['+PROGRAM+'] Error: '+str(e))
        sys.exit(1)

    if (DFQ1 != '' and DBAM != '') or (DFQ1 == '' and DBAM == ''):
        shell.eprint('['+PROGRAM+'] Error: DNA fastq file or bamfile should be provided only one')
        sys.exit(1)

    if (RFQ1 != '' and RBAM != '') or (RFQ1 == '' and RBAM == ''):
        shell.eprint('['+PROGRAM+'] Error: RNA fastq file or bamfile should be provided only one')
        sys.exit(1)

    if (DBAM != '' and RBAM != '') and (SOAPNUKE != '' or PILON != ''):
        shell.eprint('['+PROGRAM+'] Error: there are some logical error. Offering bam file means there is no need to do soapnuke or pilon')
        sys.exit(1)

    if DFQ1 != '':
        PDFQ1 = shell.checkFqQuality(DFQ1)
    if RFQ1 != '':
        PRFQ1 = shell.checkFqQuality(RFQ1)
    
    if SOAPNUKE != '':
        FLAG += 2
    
    if PILON != '':
        FLAG += 4

    os.chdir(OUTDIR)

    logging.basicConfig(level=logging.INFO,
                        filename=OUTDIR+PROGRAM+'.log',
                        filemode='w',
                        format='%(asctime)s : %(message)s',
                        datefmt='%m/%d/%Y %I:%M:%S %p')

    return config , genomefile


def main():
    global DBAM , RBAM , DFQ1 , DFQ2 , RFQ1 , RFQ2 , PDFQ1 , PRFQ1
    config , genomefile = get_config()
    logging.info('Program Start')

    if FLAG & 2:
        logging.info('do soapnuke')
        if DFQ1 != '':
            DFQ1 , DFQ2 , PDFQ1 = soapnuke('DNA', DFQ1, DFQ2, PDFQ1, config)
        if RFQ1 != '':
            RFQ1 , RFQ2 , PRFQ1 = soapnuke('RNA', RFQ1, RFQ2, PRFQ1, config)
        logging.info('soapnuke has done')

    if FLAG & 4:
        logging.info('do pilon')
        logging.info('check genome index')
        for suffix in ['.amb', '.ann', '.bwt', '.pac', '.sa']:
            if os.path.isfile(genomefile+suffix) == False:
                logging.info('genome index not found')
                p = subprocess.Popen([BWA, 'index', genomefile], stdout=open('log', 'w'), stderr=subprocess.STDOUT)
                p.wait()
                if p.returncode != 0:
                    shell.eprint('['+PROGRAM+'] Error: bwa index run error')
                    shell.eprint(''.join(open('log', 'r').read()))
                    sys.exit(1)
                logging.info('genome index done')
                break
        try:
            os.makedirs('pilon/')
        except FileExistsError:
            pass
        except:
            shell.eprint('['+PROGRAM+'] Error: mkdir pilon/ error')
            sys.exit(1)
        args = []
        for k,v in config['bwa']['mem'].items():
            args += (k, v)
        args += [genomefile, DFQ1] if DFQ2 == '' else [genomefile, DFQ1, DFQ2]
        p = subprocess.Popen([BWA, 'mem']+args, stdout=subprocess.PIPE, stderr=open('log', 'w'))
        infile = pysam.AlignmentFile(p.stdout, 'r')
        oufile = pysam.AlignmentFile('pilon/mem_just.bam', 'wb', template=infile)
        for r in infile:
            oufile.write(r)
        infile.close()
        oufile.close()
        p.wait()
        if p.returncode != 0:
            shell.eprint('['+PROGRAM+'] Error: bwa mem run error')
            shell.eprint(''.join(open('log', 'r').read()))
            sys.exit(1)
        DBAM_mem = bamsortindex('pilon/mem_just.bam')
        try:
            os.makedirs('pilon/sub/')
        except FileExistsError:
            pass
        except:
            shell.eprint('['+PROGRAM+'] Error: mkdir pilon/sub/ error')
            sys.exit(1)
        shell.splitFa(genomefile, 'pilon/sub/')
        f = pysam.AlignmentFile(DBAM_mem, 'rb')
        for fn in ['pilon/sub/'+bed for bed in os.listdir('pilon/sub/') if bed.endswith('.bed')]:
            fw = pysam.AlignmentFile('.'.join(fn.split('.')[:-1])+'.bam', 'wb', template=f)
            with open(fn, 'r') as fg:
                for line in fg:
                    for r in f.fetch(contig=line.strip().split()[0]):
                        if r.flag & 256 or r.flag & 2048:
                            continue
                        fw.write(r)
            fw.close()
            try:
                pysam.index('.'.join(fn.split('.')[:-1])+'.bam')
            except:
                shell.eprint('['+PROGRAM+'] Error: '+'.'.join(fn.split('.')[:-1])+'.bam'+' index error')
        f.close()

        for fn in ['pilon/sub/'+fa for fa in os.listdir('pilon/sub/') if fa.endswith('.fa')]:
            maxn = int(subprocess.getstatusoutput('awk \'BEGIN{x};{x+=$3};END{print x}\' ' + '.'.join(fn.split('.')[:-1])+'.bed')[1])
            maxmem = '-Xmx8g' if maxn <= 100000000 else '-Xmx16g'
            p =  subprocess.Popen([JAVA, maxmem, '-jar', PILON, '--fix', 'snps', '--genome', fn, '--bam', '.'.join(fn.split('.')[:-1])+'.bam', '--output', fn+'_fix_snps'], stdout=open('log', 'w'), stderr=subprocess.STDOUT)
            p.wait()
            if p.returncode != 0:
                shell.eprint('['+PROGRAM+'] Error: pilon run error')
                shell.eprint(''.join(open('log', 'r').read()))
                sys.exit(1)
        os.system('cat pilon/sub/*_fix_snps.fasta > '+genomefile+'.fix_snps.fa')
        os.system('sed -i s/_pilon// '+genomefile+'.fix_snps.fa')
        os.system('rm -rf pilon/sub/*')
        genomefile = genomefile + '.fix_snps.fa'
        logging.info('pilon has done')
    # check and do bwa aln
    logging.info('do bwa aln')
    logging.info('check fixed genome index')
    for suffix in ['.amb', '.ann', '.bwt', '.pac', '.sa']:
        if os.path.isfile(genomefile+suffix) == False:
            logging.info('fixed genome not found')
            p = subprocess.Popen([BWA, 'index', genomefile], stdout=open('log', 'w'), stderr=subprocess.STDOUT)
            p.wait()
            if p.returncode != 0:
                shell.eprint('['+PROGRAM+'] Error: bwa index for fixed snps genome run error')
                shell.eprint(''.join(open('log', 'r').read()))
                sys.exit(1)
            logging.info('fixed genome index done')
            break
    try:
        os.makedirs('aln/')
    except FileExistsError:
        pass
    except:
        shell.eprint('['+PROGRAM+'] Error: make aln directory error')
        sys.exit(1)
    if DFQ1 != '' and DBAM == '':
        DBAM = bwaaln('DNA', genomefile, DFQ1, DFQ2, PDFQ1, 'aln/', config)
    if RFQ1 != '' and RBAM == '':
        RBAM = bwaaln('RNA', genomefile, RFQ1, RFQ2, PRFQ1, 'aln/', config)
    logging.info('bwa aln has done')
    logging.info('do get best bam')
    if DBAM:
        bestuniqbam(DBAM, DNA=True, **config['bestuniqbam']['DNA'])
        DBAM = bamsortindex('.'.join(os.path.basename(DBAM).split('.')[:-1])+'.best.bam')
        logging.info('The output file of DNA alignment is '+DBAM)
    if RBAM:
        bestuniqbam(RBAM, RNA=True, **config['bestuniqbam']['RNA'])
        if 'SS' in config['bestuniqbam']['RNA']:
            bamsortindex('.'.join(os.path.basename(RBAM).split('.')[:-1])+'.negative.bam')
            bamsortindex('.'.join(os.path.basename(RBAM).split('.')[:-1])+'.positive.bam')
            logging.info('The output files of strand-specific RNA alignment are '+'.'.join(os.path.basename(RBAM).split('.')[:-1])+'.negative.bam and '+'.'.join(os.path.basename(RBAM).split('.')[:-1])+'.positive.bam')
        else:
            bamsortindex('.'.join(os.path.basename(RBAM).split('.')[:-1])+'.best.bam')
            logging.info('The output file of RNA alignment is '+'.'.join(os.path.basename(RBAM).split('.')[:-1])+'.best.bam')
    os.system('rm -f log')
    logging.info('All things have been done! Have a good day!')

    return 0




if __name__ == '__main__':
    main()
