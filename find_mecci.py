#!/usr/bin/env python

# This program aims to identify mecciRNAs from RNA-seq data
# I started to write this program on 2020/11/17, and I got all things done on 2020/12/15
# School of Life Sciences, University of Science & Technology of China, Yang Yan
# wxlyy@mail.ustc.edu.cn
# 2020/11/18

# This program is accomplished at a huge cost.

import os
import re
import gzip
import argparse
import copy
from collections import defaultdict
import sys
import logging
import subprocess
import math
import signal
import time
BIG = 100000

logging.basicConfig(
    format='[%(asctime)s %(levelname)s] %(message)s',
    stream=sys.stdout
)
log = logging.getLogger(__name__)


# This function is from ENCODE project to call Linux in Python
def run_shell_cmd(cmd):
    p = subprocess.Popen(
        ['/bin/bash', '-o', 'pipefail'],
        stdin=subprocess.PIPE,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        universal_newlines=True,
        preexec_fn=os.setsid
    )  # to make a new process with a new PGID
    pid = p.pid
    pgid = os.getpgid(pid)
    log.info('run_shell_cmd: PID={},PGID={},CMD={}'.format(pid, pgid, cmd))
    stdout, stderr = p.communicate(cmd)
    rc = p.returncode
    err_str = 'PID={}, PGID={}, RC={}\nSTDERR={}\nSTDOUT={}'.format(
        pid, pgid, rc, stderr.strip(), stdout.strip()
    )
    if rc:
        # kill all child processes
        try:
            os.killpg(pgid, signal.SIGKILL)
        except:
            pass
        finally:
            raise Exception(err_str)
    else:
        log.info(err_str)
    return stdout.strip('\n')


# Function stripping the extension of a fasta file
def strip_ext_fasta(fasta):
    return re.sub(r'\.(fasta|fa|Fasta|Fa)\$', '',
                  str(fasta))


# Function making a directory
def mkdir_p(dirname):
    if dirname and not os.path.exists(dirname):
        os.makedirs(dirname)


# Function removing a file forcely
def rm_f(files):
    if files:
        if type(files) == list:
            run_shell_cmd('rm -f {}'.format(' '.join(files)))
        else:
            run_shell_cmd('rm -f {}'.format(files))


# list all the files in a directory
def ls_l(d):
    cmd = 'ls -l {}'.format(d)
    run_shell_cmd(cmd)


# parsing arguments
def args_parser():
    parser = argparse.ArgumentParser(
        description='Identifying mecciRNAs from directional paired-end RNA-seq data')
    parser.add_argument('--sample', type=str, required=True,
                        help='The sample name of the RNA-seq dataset')
    parser.add_argument('--fastq1', type=str, required=True,
                        help='Files with #1 mates, paired with files in <m2>.')
    parser.add_argument('--fastq2', type=str, required=True,
                        help='Files with #2 mates, paired with files in <m1>.')
    parser.add_argument('--platform', type=str, choices=[
                        'Novo', 'BGI'], default='Novo', help='The platform the RNA-seq dataset was generated')
    parser.add_argument('--mito', type=str, default='chrM',
                        help='The name of the mitochondrial chromosome')
    parser.add_argument('--genome', type=str, required=True,
                        help='The reference fasta file')
    parser.add_argument('--bwa', type=str, default='bwa',
                        help='The path to BWA')
    parser.add_argument('--bowtie', type=str,
                        default='bowtie', help='The path to bowtie')
    parser.add_argument('--bowtie_build', type=str,
                        default='bowtie-build', help='The path to bowtie-build')
    parser.add_argument('--bwaindex', type=str, required=True,
                        help='The path to the bwa index of reference genome')
    parser.add_argument('--mitoidx_dir', type=str, required=True,
                        help='The directory used to store the index of mitochondrial genome during the pipeline')
    parser.add_argument('--seqtk', type=str, default='seqtk',
                        help='The path to seqtk (https://github.com/lh3/seqtk)')
    parser.add_argument('--mismatch_bwa', type=int, default=2,
                        help='The number of maximum mismatches allowed during mapping reads against reference genome')
    parser.add_argument('--mismatch_anchor', type=int, default=2,
                        help='The number of maximum mismatches allowed during mapping anchors against reference genome')
    parser.add_argument('--mismatch_internal', type=int, default=2,
                        help='The number of maximum mismatches allowed during mapping internal sequences against reference sequences')
    parser.add_argument('--threads', type=int, default=8,
                        help='The number of threads to map reads using BWA')
    parser.add_argument('--anchor_size', type=int, default=15,
                        help='The length of anchors extracted from both ends of unmapped reads to find the breakpoint of the sequence')
    parser.add_argument('--tmpdir', type=str, default='tmp',
                        help='The name of the temporary directory')
    parser.add_argument('--result', type=str, default='result',
                        help='The name of the directory storing the final results')
    parser.add_argument('--log-level', default='INFO',
                        choices=['NOTSET', 'DEBUG', 'INFO',
                                 'WARNING', 'CRITICAL', 'ERROR',
                                 'CRITICAL'],
                        help='Log level')
    args = parser.parse_args()
    log.setLevel(args.log_level)
    log.info(sys.argv)
    return args


# define a class representing a read in the fastq file
class Fastq(object):
    def __init__(self, name='', seq='', mark='+', quality=''):
        self.name = name.split(' ')[0].split('/')[0][1:]
        self.seq = seq
        self.mark = mark
        self.qual = quality
        self.rc_seq = rev_comp(self.seq)
        self.rc_qual = self.qual[::-1]

    def write_to_file(self, fh):
        fh.write("@" + self.name + "\n" + self.seq +
                 "\n" + self.mark + "\n" + self.qual + "\n")

    def write_rc_to_file(self, fh):
        fh.write("@" + self.name + "\n" + self.rc_seq +
                 "\n" + self.mark + "\n" + self.rc_qual + "\n")


# The dict indicating complement
COMPLEMENT = {
    'a': 't',
    't': 'a',
    'c': 'g',
    'g': 'c',
    'k': 'm',
    'm': 'k',
    'r': 'y',
    'y': 'r',
    's': 's',
    'w': 'w',
    'b': 'v',
    'v': 'b',
    'h': 'd',
    'd': 'h',
    'n': 'n',
    'A': 'T',
    'T': 'A',
    'C': 'G',
    'G': 'C',
    'K': 'M',
    'M': 'K',
    'R': 'Y',
    'Y': 'R',
    'S': 'S',
    'W': 'W',
    'B': 'V',
    'V': 'B',
    'H': 'D',
    'D': 'H',
    'N': 'N',
}


def complement(s):
    return "".join([COMPLEMENT[x] for x in s])


# return the reversed complement sequences
def rev_comp(seq):
    return complement(seq)[::-1]


# Function judging whether interval [start1, end1] is the subset of interval [start2, end2]
def Is_subset(start1, end1, start2, end2):
    if start1 >= start2 and end1 <= end2:
        return True
    return False


# KMP class
class KMP:
    def partial(self, pattern):
        """ Calculate partial match table: String -> [Int]"""
        ret = [0]
        for i in range(1, len(pattern)):
            j = ret[i - 1]
            while j > 0 and pattern[j] != pattern[i]:
                j = ret[j - 1]
            ret.append(j + 1 if pattern[j] == pattern[i] else j)
        return ret

    def search(self, T, P):
        """
        KMP search main algorithm: String -> String -> [Int]
        Return all the matching position of pattern string P in T
        """
        partial, ret, j = self.partial(P), [], 0
        for i in range(len(T)):
            while j > 0 and T[i] != P[j]:
                j = partial[j - 1]
            if T[i] == P[j]:
                j += 1
            if j == len(P):
                ret.append(i - (j - 1))
                j = partial[j - 1]
        return ret


# KMP function to exact matching
def match_str_kmp(T, P):
    kmp = KMP()
    return kmp.search(T, P)


# Function calculating the edit distance between two strings
def minDistance(str1, str2):
    len1 = len(str1)
    len2 = len(str2)
    dp = []
    for row in range(len1 + 1):
        tmp = []
        for col in range(len2 + 1):
            if row == 0:
                tmp.append(col)
            elif col == 0:
                tmp.append(row)
            else:
                tmp.append(False)
        dp.append(tmp)
    for row in range(1, len1 + 1):
        for col in range(1, len2 + 1):
            if str1[row - 1] == str2[col - 1]:
                dp[row][col] = dp[row - 1][col - 1]
            else:
                dp[row][col] = min(dp[row - 1][col], dp[row - 1]
                                   [col - 1], dp[row][col - 1]) + 1
    return dp[len1][len2]


# define a class to represent a read in the fasta file
class Fasta(object):
    def __init__(self, name='', seq=''):
        self.name = name.split(' ')[0].lstrip('>')
        self.seq = seq

    def get_len(self):
        return len(self.seq)

    def write_to_file(self, fh):
        fh.write(">" + self.name + "\n" + self.seq + "\n")


class Circ(object):
    def __init__(self, chr='chrM', start=0, end=0, name='', reads_num=0, strand='+', reads='', visited=False):
        self.chr = chr
        self.start = int(start)
        self.end = int(end)
        self.name = name
        self.reads_num = int(reads_num)
        self.strand = strand
        self.reads = reads
        self.visited = visited


class Circ_final(object):
    def __init__(self, chr='', start=0, end=0, name='', reads_num=0, strand='+', reads_list=[]):
        self.chr = chr
        self.start = int(start)
        self.end = int(end)
        self.name = name
        self.reads_num = int(reads_num)
        if strand == 1:
            self.strand = '+'
        elif strand == 2:
            self.strand = '-'
        else:
            self.strand = '*'
        self.reads_list = reads_list
        self.visited = False

    def length(self):
        return self.end - self.start

    def write_to_file(self, fh):
        fh.write('\t'.join([self.chr, str(self.start), str(self.end), self.name, str(
            self.reads_num), self.strand, ','.join(self.reads_list)]) + "\n")


# Class representing a mecciRNA to deal with the repeat within mecciRNA
class Mecci(object):
    def __init__(self, chr, start, end, name, read_num=0, strand='+', reads_list={}, group=0, visited=False):
        self.chr = chr
        self.start = int(start)
        self.end = int(end)
        self.name = name
        self.read_num = int(read_num)
        self.strand = strand
        self.reads_list = reads_list
        self.group = int(group)
        self.visited = visited


# Function parsing the fastq file
def parse_fastq(fqfile):
    name2fq = {}  # This is bad because the use of dictionary will occupy a lot of memory ,and a better method may be the use of generator
    with open_file(fqfile) as fq:
        line_num = 0
        name = ''
        seq = ''
        mark = ''
        quality = ''
        for line in fq:
            if line_num % 4 == 0:
                name = line.rstrip('\n').split(" ")[0].split('/')[0][1:]
            elif line_num % 4 == 1:
                seq = line.rstrip('\n')
            elif line_num % 4 == 2:
                mark = line.rstrip('\n')
            else:
                quality = line.strip()
            name2fq[name] = Fastq(name, seq, mark, quality)
            line_num += 1
    return name2fq


# Function parsing a fasta file and returning a dict
def parse_fasta(fasta):
    # I use the dict to map seqname to sequence,
    # so all the sequences will be read into the memory. It's unwise.
    # Maybe there are better methods, such as to use a generator.
    name2fa = dict()
    seq_num = 0
    name_list = []
    seq = ''
    with open_file(fasta) as fa:
        for line in fa:
            line = line.rstrip('\n')
            if re.match('^>', line):
                seqname = re.sub('^>', '', line).split(
                    ' ')[0]  # seqname = line.lstrip('>')
                print("Reading the sequence of {} ...\n".format(seqname))
                name_list.append(seqname)
                if seq_num > 0:
                    name2fa[name_list[seq_num - 1]
                            ] = Fasta(name_list[seq_num - 1], seq)
                seq = ''  # clear seq to be ready to store the next sequence
                seq_num += 1
            else:
                seq += line
        # the last sequence in the fasta file
        name2fa[name_list[seq_num - 1]] = Fasta(name_list[seq_num - 1], seq)
        return name2fa


# class representing a read of a fasta file
class Read(object):
    def __init__(self, name='', seq=''):
        self.name = name
        self.seq = seq

    def length(self):
        return len(self.seq)

    def rc(self):
        return rev_comp(self.seq)

    def get_shortname(self):
        return ''.join(self.name.split('/')[0].split(' ')[0])


# class representing a paired-end, including its name, sequence, quality and details about mapping
class PE(object):
    def __init__(self, name='', seq1='', seq2='', mark='+', qual1='', qual2='', anchor_size=15):
        self.name = name
        self.mark = mark
        self.fwseq = seq1
        self.fwseq_rc = rev_comp(seq1)
        self.fwseq_left = seq1[:anchor_size]
        self.fwseq_right = seq1[-anchor_size:]
        self.fwseq_rc_left = rev_comp(seq1)[:anchor_size]
        self.fwseq_rc_right = rev_comp(seq1)[-anchor_size:]
        self.rvseq = seq2
        self.rvseq_rc = rev_comp(seq2)
        self.rvseq_left = seq2[:anchor_size]
        self.rvseq_right = seq2[-anchor_size:]
        self.rvseq_rc_left = rev_comp(seq2)[:anchor_size]
        self.rvseq_rc_right = rev_comp(seq2)[-anchor_size:]
        self.fwqual = qual1
        self.fwqual_rc = qual1[::-1]
        self.fwqual_left = qual1[:anchor_size]
        self.fwqual_right = qual1[-anchor_size:]
        self.fwqual_rc_left = qual1[-anchor_size:]
        self.fwqual_rc_right = qual1[:anchor_size]
        self.rvqual = qual2
        self.rvqual_rc = qual2[::-1]
        self.rvqual_left = qual2[:anchor_size]
        self.rvqual_right = qual2[-anchor_size:]
        self.rvqual_rc_left = qual2[-anchor_size:]
        self.rvqual_rc_right = qual2[:anchor_size]

    def write_R1_rc_to_file(self, fh):
        fh.write(
            '\n'.join(["@"+self.name, self.fwseq_rc, self.mark, self.fwqual_rc]) + '\n')

    def write_R1_left_to_file(self, fh):
        fh.write('\n'.join(["@"+self.name, self.fwseq_left,
                            self.mark, self.fwqual_left]) + '\n')

    def write_R1_right_to_file(self, fh):
        fh.write('\n'.join(["@"+self.name, self.fwseq_right,
                            self.mark, self.fwqual_right]) + '\n')

    def write_R1_RC_left_to_file(self, fh):
        fh.write('\n'.join(["@"+self.name, self.fwseq_rc_left,
                            self.mark, self.fwqual_rc_left]) + '\n')

    def write_R1_RC_right_to_file(self, fh):
        fh.write('\n'.join(["@"+self.name, self.fwseq_rc_right,
                            self.mark, self.fwqual_rc_right]) + '\n')

    def write_R2_rc_to_file(self, fh):
        fh.write(
            '\n'.join(["@"+self.name, self.rvseq_rc, self.mark, self.rvqual_rc]) + '\n')

    def write_R2_left_to_file(self, fh):
        fh.write('\n'.join(["@"+self.name, self.rvseq_left,
                            self.mark, self.rvqual_left]) + '\n')

    def write_R2_right_to_file(self, fh):
        fh.write('\n'.join(["@"+self.name, self.rvseq_right,
                            self.mark, self.rvqual_right]) + '\n')

    def write_R2_RC_left_to_file(self, fh):
        fh.write('\n'.join(["@"+self.name, self.rvseq_rc_left,
                            self.mark, self.rvqual_rc_left]) + '\n')

    def write_R2_RC_right_to_file(self, fh):
        fh.write('\n'.join(["@"+self.name, self.rvseq_rc_right,
                            self.mark, self.rvqual_rc_right]) + '\n')

    def get_R1_map_pos(self, map_pos_dict=defaultdict(dict)):
        return map_pos_dict[self.name]['R1']

    def get_R1_rc_map_pos(self, map_pos_dict=defaultdict(dict)):
        return map_pos_dict[self.name]['R1_RC']

    def get_R1_left_map_pos(self, map_pos_dict=defaultdict(dict)):
        return map_pos_dict[self.name]['R1_left']

    def get_R1_right_map_pos(self, map_pos_dict=defaultdict(dict)):
        return map_pos_dict[self.name]['R1_right']

    def get_R1_rc_left_map_pos(self, map_pos_dict=defaultdict(dict)):
        return map_pos_dict[self.name]['R1_RC_left']

    def get_R1_rc_right_map_pos(self, map_pos_dict=defaultdict(dict)):
        return map_pos_dict[self.name]['R1_RC_right']

    def get_R2_map_pos(self, map_pos_dict=defaultdict(dict)):
        return map_pos_dict[self.name]['R2']

    def get_R2_rc_map_pos(self, map_pos_dict=defaultdict(dict)):
        return map_pos_dict[self.name]['R2_RC']

    def get_R2_left_map_pos(self, map_pos_dict=defaultdict(dict)):
        return map_pos_dict[self.name]['R2_left']

    def get_R2_right_map_pos(self, map_pos_dict=defaultdict(dict)):
        return map_pos_dict[self.name]['R2_right']

    def get_R2_rc_left_map_pos(self, map_pos_dict=defaultdict(dict)):
        return map_pos_dict[self.name]['R2_RC_left']

    def get_R2_rc_right_map_pos(self, map_pos_dict=defaultdict(dict)):
        return map_pos_dict[self.name]['R2_RC_right']


# Get all the position of associated sequences of a pair of reads
def get_PE_MAP(mismatch_anchor=2, fastq1='', fastq2='', fastq1_rc='', fastq2_rc='', fastq1_left='', fastq1_right='', fastq2_left='',
               fastq2_right='', fastq1_rc_left='', fastq1_rc_right='', fastq2_rc_left='', fastq2_rc_right='',
               anchor_size=15, ref_fasta='', nthread=1, bowtie='bowtie', bowtie_build='bowtie-build', idxdir='', outdir=''):
    name2fq1 = parse_fastq(fastq1)
    name2fq2 = parse_fastq(fastq2)
    fq1_rc_fh = open(fastq1_rc, 'w')
    fastq1_left_fh = open(fastq1_left, 'w')
    fastq1_right_fh = open(fastq1_right, 'w')
    fastq1_rc_left_fh = open(fastq1_rc_left, 'w')
    fastq1_rc_right_fh = open(fastq1_rc_right, 'w')
    fq2_rc_fh = open(fastq2_rc, 'w')
    fastq2_left_fh = open(fastq2_left, 'w')
    fastq2_right_fh = open(fastq2_right, 'w')
    fastq2_rc_left_fh = open(fastq2_rc_left, 'w')
    fastq2_rc_right_fh = open(fastq2_rc_right, 'w')
    pe_list = []

    # write all the associated sequences of Read1 and Read2 to new file
    for name in name2fq1:
        pe_read = PE(name, name2fq1[name].seq, name2fq2[name].seq, name2fq1[name].mark,
                     name2fq1[name].qual, name2fq2[name].qual, anchor_size)
        pe_list.append(pe_read)
        pe_read.write_R1_rc_to_file(fq1_rc_fh)
        pe_read.write_R1_left_to_file(fastq1_left_fh)
        pe_read.write_R1_right_to_file(fastq1_right_fh)
        pe_read.write_R1_RC_left_to_file(fastq1_rc_left_fh)
        pe_read.write_R1_RC_right_to_file(fastq1_rc_right_fh)
        pe_read.write_R2_rc_to_file(fq2_rc_fh)
        pe_read.write_R2_left_to_file(fastq2_left_fh)
        pe_read.write_R2_right_to_file(fastq2_right_fh)
        pe_read.write_R2_RC_left_to_file(fastq2_rc_left_fh)
        pe_read.write_R2_RC_right_to_file(fastq2_rc_right_fh)
    fq1_rc_fh.close()
    fastq1_left_fh.close()
    fastq1_right_fh.close()
    fastq1_rc_left_fh.close()
    fastq1_rc_right_fh.close()
    fq2_rc_fh.close()
    fastq2_left_fh.close()
    fastq2_right_fh.close()
    fastq2_rc_left_fh.close()
    fastq2_rc_right_fh.close()

    # Mapping R1, R2, R1_RC, R2_RC to the reference genome
    # Building the bowtie index of reference fasta file
    btindex = bowtie_idx(bowtie_build, ref_fasta, nthread, idxdir)
    # Mapping
    R1_pos = get_map_pos(bowtie, fastq1, mismatch_anchor,
                         btindex, nthread, outdir)
    R1_rc_pos = get_map_pos(
        bowtie, fastq1_rc, mismatch_anchor, btindex, nthread, outdir)
    R1_left_pos = get_map_pos(
        bowtie, fastq1_left, mismatch_anchor, btindex, nthread, outdir)
    R1_right_pos = get_map_pos(
        bowtie, fastq1_right, mismatch_anchor, btindex, nthread, outdir)
    R1_rc_left_pos = get_map_pos(
        bowtie, fastq1_rc_left, mismatch_anchor, btindex, nthread, outdir)
    R1_rc_right_pos = get_map_pos(
        bowtie, fastq1_rc_right, mismatch_anchor, btindex, nthread, outdir)
    R2_pos = get_map_pos(bowtie, fastq2, mismatch_anchor,
                         btindex, nthread, outdir)
    R2_rc_pos = get_map_pos(
        bowtie, fastq2_rc, mismatch_anchor, btindex, nthread, outdir)
    R2_left_pos = get_map_pos(
        bowtie, fastq2_left, mismatch_anchor, btindex, nthread, outdir)
    R2_right_pos = get_map_pos(
        bowtie, fastq2_right, mismatch_anchor, btindex, nthread, outdir)
    R2_rc_left_pos = get_map_pos(
        bowtie, fastq2_rc_left, mismatch_anchor, btindex, nthread, outdir)
    R2_rc_right_pos = get_map_pos(
        bowtie, fastq2_rc_right, mismatch_anchor, btindex, nthread, outdir)

    # Get the map position
    name2pos = defaultdict(dict)
    for name in name2fq1:
        name2pos[name]['R1'] = R1_pos[name]
        name2pos[name]['R1_RC'] = R1_rc_pos[name]
        name2pos[name]['R1_left'] = R1_left_pos[name]
        name2pos[name]['R1_right'] = R1_right_pos[name]
        name2pos[name]['R1_RC_left'] = R1_rc_left_pos[name]
        name2pos[name]['R1_RC_right'] = R1_rc_right_pos[name]
        name2pos[name]['R2'] = R2_pos[name]
        name2pos[name]['R2_RC'] = R2_rc_pos[name]
        name2pos[name]['R2_left'] = R2_left_pos[name]
        name2pos[name]['R2_right'] = R2_right_pos[name]
        name2pos[name]['R2_RC_left'] = R2_rc_left_pos[name]
        name2pos[name]['R2_RC_right'] = R2_rc_right_pos[name]

    return pe_list, name2pos


# Getting the matched number of bases indicated by a CIGAR in SAM-formatted file
def get_match(string):
    match_num = 0
    pattern_char = re.compile('M|I|D|N|S|H|P|=|X')
    pattern_num = re.compile(r'\d')
    list_by_char = [i for i in pattern_char.split(string) if i]
    list_by_num = [i for i in pattern_num.split(string) if i]
    for i in range(0, len(list_by_num)):
        if list_by_num[i] == 'M':
            match_num += int(list_by_char[i])
    return match_num


# Function opening gzipped file or not
def open_file(infile, mode='r'):
    if re.match('.gz$', infile) or re.match('.gzip$', infile):
        if mode == 'r':
            mode = 'rb'
        return gzip.open(infile, mode=mode)
    else:
        return open(infile, mode=mode)


# Function building bwa index
def build_bwaindex(ref_fasta, idx_dir, bwa):
    basename = os.path.basename(strip_ext_fasta(ref_fasta))
    prefix = os.path.join(idx_dir, basename)
    cmd = '{} -p {} {}'.format(bwa, prefix, ref_fasta)
    run_shell_cmd(cmd)
    return prefix


# Function mapping reads against reference sequence using bwa mem with parameters using in CIRI2 pipeline
def bwamap(sample, bwa, fastq1, fastq2, bwaindex, nthread, outdir):
    prefix = os.path.join(outdir, sample)
    sam_file = '{}.bwa_mem.sam'.format(prefix)
    log_file = '{}.bwa_mem.log'.format(prefix)
    cmd = '{} mem -t {} -T 19 {} {} {} 1>{} 2>{}'
    cmd = cmd.format(bwa, nthread, bwaindex, fastq1,
                     fastq2, sam_file, log_file)
    run_shell_cmd(cmd)
    return sam_file


# NOTE: I use seqtk to get sequences form original fastq files very quickly
# BUT! Sometimes there may be some problems because the format of sequence name of different
# Illuimina platform may be different so the name outputed by this function may not be concordant
# with the raw fastq files and seqtk will get nothing so I write two functions
def unmappedreads2fastq_novo(sample, seqtk, fastq1, fastq2, reads_list, outdir):
    prefix = os.path.join(outdir, sample)
    reads_out = '{}.bwa_unmap.txt'.format(prefix)
    fastq1_out = '{}.bwa_unmap_1.fastq'.format(prefix)
    fastq2_out = '{}.bwa_unmap_2.fastq'.format(prefix)
    reads_out_fh = open(reads_out, 'w')
    reads_out_fh.write('\n'.join(reads_list))
    reads_out_fh.close()

    cmd1 = '{} subseq {} {} >{}'.format(seqtk, fastq1, reads_out, fastq1_out)
    cmd2 = '{} subseq {} {} >{}'.format(seqtk, fastq2, reads_out, fastq2_out)
    run_shell_cmd(cmd1)
    run_shell_cmd(cmd2)

    return fastq1_out, fastq2_out


# NOTE: I use seqtk to get sequences form original fastq files very quickly
# BUT! Sometimes there may be some problems because the format of sequence name of different
# Illuimina platform may be different so the name outputed by this function may not be concordant
# with the raw fastq files and seqtk will get nothing so I write two functions
def unmappedreads2fastq_BGI(sample, seqtk, fastq1, fastq2, reads_list, outdir):
    prefix = os.path.join(outdir, sample)
    reads_out = '{}.bwa_unmap.txt'.format(prefix)
    fastq1_out = '{}.bwa_unmap_1.fastq'.format(prefix)
    fastq2_out = '{}.bwa_unmap_2.fastq'.format(prefix)
    reads1_out_fh = open(fastq1_out, 'w')
    reads2_out_fh = open(fastq2_out, 'w')
    reads1_out_fh.write('\n'.join([i + '/1' for i in reads_list]))
    reads2_out_fh.write('\n'.join([i + '/2' for i in reads_list]))
    reads1_out_fh.close()
    reads2_out_fh.close()

    cmd1 = '{} subseq {} {} >{}'.format(seqtk, fastq1, reads_out, fastq1_out)
    cmd2 = '{} subseq {} {} >{}'.format(seqtk, fastq2, reads_out, fastq2_out)
    run_shell_cmd(cmd1)
    run_shell_cmd(cmd2)

    return fastq1_out, fastq2_out


# Estimating the length of reads in a fastq file
def get_read_length(fastq):
    line_num = 0
    sample_reads_num = 100000
    analyzed_reads_num = 0
    readLen2count = defaultdict(int)
    with open_file(fastq) as fh:
        for line in fh:
            if line_num % 4 == 1:
                readLen2count[len(line.strip())] += 1
                analyzed_reads_num += 1
            if analyzed_reads_num >= sample_reads_num:
                break
            line_num += 1
    max_occurrence_num = 0
    read_len = 0
    for length in readLen2count:
        if readLen2count[length] > max_occurrence_num:
            max_occurrence_num = readLen2count[length]
            read_len = length
    return read_len


# Function judging whether a CIGAR in a samfile outputed by bwa mem indicating a read is matched or not
def Is_cigar_match(CIGAR, read_len, mismatch=2):
    match_num = get_match(CIGAR)
    if match_num >= read_len - mismatch:
        return True
    return False


# remove redundant elements in a list
def remove_repeat(list_a):
    list_b = [i for i in set(list_a)]
    return list_b


# Function estimating whether two lists have any common element
def Is_overlap(list_a, list_b):
    for a in list_a:
        for b in list_b:
            if a == b:
                return True
    return False


# Function getting the common elements in two lists
def common_elements(list_a, list_b):
    common_eles = []
    list_a = remove_repeat(list_a)
    list_b = remove_repeat(list_b)
    for i in list_a:
        for j in list_b:
            if i == j:
                common_eles.append(i)
    return common_eles


# Getting unmapped reads list from samfile outputed by the function bowtie_se
def get_unmapped_reads(mito_name='chrM', sam_file='', read_len=0, mismatch=2):
    unmapped_reads_list = []
    with open_file(sam_file) as fh:
        for line in fh:
            line = line.strip()
            if re.match('^@', line):
                continue
            lineL = line.split('\t')
            if lineL[2] == mito_name and lineL[6] == "=" and \
                    not Is_cigar_match(lineL[5], read_len, mismatch):
                unmapped_reads_list.append(lineL[0])
    unmapped_reads_list = remove_repeat(unmapped_reads_list)
    return unmapped_reads_list


# split a string trying to find a proper breakpoint to split the string into two parts
# to map the refseq according to the swapped orientation
def split_seq(string='', name='', refseq='', r_pos=[], l_pos=[], anchor_size=15, mismatch_split=2):
    str_len = len(string)
    # left_anchor = string[:anchor_size]  # r_pos: list storing the mapping position of left anchor
    # right_anchor = string[-anchor_size:]  # l_pos: list storing the mapping position of right anchor
    count = 0
    """For a string having more than one method to split the string,
    I will keep the splitting method having the least number of mismatch
    """
    split2mismatch = {}  # {NO of splitting way : number of mismatch}
    split2start = {}
    split2end = {}
    if l_pos and r_pos:
        for i in l_pos:
            for j in r_pos:
                if i < j:
                    for k in range(anchor_size, str_len - anchor_size):
                        count = k - anchor_size
                        ref_right = refseq[j:j+k]
                        ref_left = refseq[i+k+anchor_size -
                                          str_len:i+anchor_size]
                        ref_shuffle = ref_right + ref_left
                        m = minDistance(ref_shuffle, string)
                        if m < mismatch_split:
                            split2mismatch[count] = m
                            split2start[count] = i + k + anchor_size - str_len
                            split2end[count] = j + k
    min_mismatch = BIG
    output2interval = defaultdict(dict)
    for i in split2mismatch:
        if split2mismatch[i] < min_mismatch:
            min_mismatch = split2mismatch[i]
    out_num = 0
    for i in split2mismatch:
        if split2mismatch[i] == min_mismatch:
            output2interval[out_num] = {
                'start': split2start[i], 'end': split2end[i]}
            out_num += 1
    return output2interval


"""
Function mapping reads against reference using bowtie
to find all the possible alignment position allowing mismatch
"""


def bowtie_idx(bowtie_build, ref_fasta, nthread, outdir):
    basename = os.path.basename(strip_ext_fasta(ref_fasta))
    prefix = os.path.join(outdir, basename)
    cmd = '{} -f --threads {} {} {}'
    cmd = cmd.format(bowtie_build, nthread, ref_fasta, prefix)
    run_shell_cmd(cmd)
    return prefix


# Map fastq file against reference using bowtie
def bowtie_se(bowtie='', fastq='', btindex='', nthread=1, mismatch=2, outdir=''):
    prefix = os.path.basename(fastq)
    #prefix = os.path.join(outdir, prefix)
    sam_file = os.path.join(outdir, '{}.bowtie.sam'.format(prefix))
    cmd = '{} -q -a --best --strata --norc --no-unal --threads {}  -v {} {} {} -S {}'
    cmd = cmd.format(bowtie, nthread, mismatch, btindex, fastq, sam_file)
    run_shell_cmd(cmd)
    return sam_file


# Parse samfile to get a dictionary :{readName : mapping position in the mitochondrial genome}
def parse_samfile(samfile):
    name2pos = defaultdict(list)
    with open_file(samfile) as fh:
        for line in fh:
            if re.match('^@', line):
                continue
            line = line.strip()
            lineL = line.split('\t')
            # because the mapping position in sam format file is 1-based
            name, pos = str(lineL[0]), int(lineL[3]) - 1
            name2pos[name].append(pos)
    return name2pos


# parse bed file which contains circular RNAs
def parse_bed(circbed):
    circ_lsit = []
    with open_file(circbed) as fh_circ:
        for line in fh_circ:
            line = line.strip()
            lineL = line.split('\t')
            chr, start, end, name, reads_num, strand, reads_list = lineL
            reads_list = reads_list.split(',')
            circ = Mecci(chr, start, end, name, reads_num, strand, reads_list)
            circ_lsit.append(circ)
    return circ_lsit


# Mapping reads against reference genome using bowtie and get the position of each read
def get_map_pos(bowtie, fastq, mismatch, bowtie_index, nthread, outdir):
    samfile = bowtie_se(bowtie, fastq, bowtie_index, nthread, mismatch, outdir)
    name2pos = parse_samfile(samfile)
    return name2pos


# This function accepts the bwa-unmapped fastq1 and fastq2 file
# This function is the master of the program
def find_mecci(mito='', sample='', bowtie='', bowtie_build='', fastq1='', fastq2='', mito_fa='', idxdir='', outdir='', nthread=1, anchor_size=15, mismatch_split=2, mismatch_consecutive=2):
    """
    :param mito: the name of mitochondrial chromosome
    :param sample: the name of sample
    :param bowtie: the path of bowtie
    :param fastq1: the path of fastq1
    :param fastq2: the path of fastq2
    :param mito_fa: the fasta file storing the sequences of mitochondrial genome
    :param idxdir: the directory which is used to store the bowtie index
    :param outdir: the directory which is used to store the results
    :param anchor_size: the length of anchor(default: 15nt)
    :param mismatch_split: the number of mismatch which is allowed in the alignment of anchors
    :param mismatch_consecutive: the number of mismatch which is allowed in the alignment of internal sequences
    :return:mecciRNAs
    """
    ref_fa = parse_fasta(mito_fa)[
        mito]  # The fasta object of mitochondrial genome
    ref_name = ref_fa.name
    ref_seq = ref_fa.seq
    prefix = os.path.join(outdir, sample)
    fastq1_rc = '{}_1_RC.fastq'.format(prefix)
    fastq1_left = '{}_1_left.fastq'.format(prefix)
    fastq1_right = '{}_1_right.fastq'.format(prefix)
    fastq1_rc_left = '{}_1_RC_left.fastq'.format(prefix)
    fastq1_rc_right = '{}_1_RC_right.fastq'.format(prefix)
    fastq2_rc = '{}_2_RC.fastq'.format(prefix)
    fastq2_left = '{}_2_left.fastq'.format(prefix)
    fastq2_right = '{}_2_right.fastq'.format(prefix)
    fastq2_rc_left = '{}_2_RC_left.fastq'.format(prefix)
    fastq2_rc_right = '{}_2_RC_right.fastq'.format(prefix)

    circ_list = []  # The list saving circular RNAs
    pe_list, name2pos = get_PE_MAP(mismatch_split, fastq1, fastq2, fastq1_rc, fastq2_rc, fastq1_left, fastq1_right, fastq2_left, fastq2_right,
                                   fastq1_rc_left, fastq1_rc_right, fastq2_rc_left, fastq2_rc_right, anchor_size, mito_fa, nthread, bowtie, bowtie_build, idxdir, outdir)
    count = 0
    for read in pe_list:
        count += 1
        sys.stderr.write('parsing {}th read {}\n'.format(count, read.name))
        read1_len = len(read.fwseq)
        read2_len = len(read.rvseq)
        read1_map_start = name2pos[read.name]['R1']
        read1_rc_map_start = name2pos[read.name]['R1_RC']
        read2_map_start = name2pos[read.name]['R2']
        read2_rc_map_satrt = name2pos[read.name]['R2_RC']
        read1_interval = split_seq(read.fwseq, read.name, ref_seq, read.get_R1_left_map_pos(
            name2pos), read.get_R1_right_map_pos(name2pos), anchor_size, mismatch_split)
        read1_RC_interval = split_seq(read.fwseq_rc, read.name, ref_seq, read.get_R1_rc_left_map_pos(
            name2pos), read.get_R1_rc_right_map_pos(name2pos), anchor_size, mismatch_split)
        read2_interval = split_seq(read.rvseq, read.name, ref_seq, read.get_R2_left_map_pos(
            name2pos), read.get_R2_right_map_pos(name2pos), anchor_size, mismatch_split)
        read2_RC_interval = split_seq(read.rvseq_rc, read.name, ref_seq, read.get_R2_rc_left_map_pos(
            name2pos), read.get_R2_rc_right_map_pos(name2pos), anchor_size, mismatch_split)
        strand = 0  # The variable indicating the strand mecciRNAs originated from; 1:plus 2:minus
        """
        Case1: read1 and read2 are both junction-spanning and therefore have the same split point
        """
        for i in read1_interval:
            start1 = read1_interval[i]['start']
            end1 = read1_interval[i]['end']
            for j in read2_RC_interval:
                start2 = read2_RC_interval[j]['start']
                end2 = read2_RC_interval[j]['end']
                if start1 == start2 and end1 == end2:
                    strand = 2
                    print("mecci found!")
                    circ_list.append(
                        Circ(ref_name, start1, end1, 'circ', 1, strand, read.name))
        for i in read2_interval:
            start2 = read2_interval[i]['start']
            end2 = read2_interval[i]['end']
            for j in read1_RC_interval:
                start1 = read1_RC_interval[j]['start']
                end1 = read1_RC_interval[j]['end']
                if start1 == start2 and end1 == end2:
                    strand = 1
                    print("mecci found!")
                    circ_list.append(
                        Circ(ref_name, start1, end1, 'circ', 1, strand, read.name))
        """
        Case2: one of the read is junction-spanning and the other is within 
        the interval which is defined by the split point of the previous read
        """
        for i in read1_interval:
            for start in read2_rc_map_satrt:
                end = start + read2_len
                if Is_subset(start, end, read1_interval[i]['start'], read1_interval[i]['end']):
                    strand = 2
                    print("mecci found!")
                    circ_list.append(Circ(
                        ref_name, read1_interval[i]['start'], read1_interval[i]['end'], 'circ', 1, strand, read.name))
        for i in read1_RC_interval:
            for start in read2_map_start:
                end = start + read2_len
                if Is_subset(start, end, read1_RC_interval[i]['start'], read1_RC_interval[i]['end']):
                    strand = 1
                    print("mecci found!")
                    circ_list.append(Circ(
                        ref_name, read1_RC_interval[i]['start'], read1_RC_interval[i]['end'], 'circ', 1, strand,
                        read.name))
        for i in read2_interval:
            for start in read1_rc_map_start:
                end = start + read1_len
                if Is_subset(start, end, read2_interval[i]['start'], read2_interval[i]['end']):
                    strand = 1
                    print("mecci found!")
                    circ_list.append(Circ(
                        ref_name, read2_interval[i]['start'], read2_interval[i]['end'], 'circ', 1, strand, read.name))
        for i in read2_RC_interval:
            for start in read1_map_start:
                end = start + read1_len
                if Is_subset(start, end, read2_RC_interval[i]['start'], read2_RC_interval[i]['end']):
                    strand = 2
                    print("mecci found!")
                    circ_list.append(Circ(
                        ref_name, read2_RC_interval[i]['start'], read2_RC_interval[i]['end'], 'circ', 1, strand,
                        read.name))
    return circ_list


# Function judging whether two circRNAs are the same
def circ_equal(circ1=Circ(), circ2=Circ()):
    if circ1.chr == circ2.chr and circ1.start == circ2.start and circ1.end == circ2.end and circ1.strand == circ2.strand:
        return True
    return False


# Function traversing a list containing all circRNAs and merge the same circRNA into a single one
def summary_circ(circ_list=[]):
    circ_num = len(circ_list)
    circ_list_summary = []
    for i in range(0, circ_num):
        if not circ_list[i].visited:
            name = circ_list[i].chr + ':' + \
                str(circ_list[i].start) + '|' + str(circ_list[i].end)
            reads_list = []
            reads_list.append(circ_list[i].reads)
            circrna = Circ_final(circ_list[i].chr, circ_list[i].start,
                                 circ_list[i].end, name, 1, circ_list[i].strand, reads_list)
            for j in range(i+1, circ_num):
                if not circ_list[j].visited and circ_equal(circ_list[i], circ_list[j]):
                    circrna.reads_num += 1
                    circrna.reads_list.append(circ_list[j].reads)
                    circrna.reads_list = remove_repeat(circrna.reads_list)
                    circ_list[j].visited = True
            circ_list_summary.append(circrna)
    return circ_list_summary


# The function is also the heart of the program
# I thought out the idea that used three "pointers"
# to indicate three RNAs before dawn on 2020.1.17
# which is the first day of the winter holiday in 2020
# when I was in warm Lab 702 and there's nobody there.
def group_mecci(meccirnas):
    # Accepts a list composed of mecciRNAs which is represented by class Mecci
    i = 0
    j = 0
    k = 0
    group2mecci = {}
    meccinum = len(meccirnas)
    for i in range(0, meccinum):
        # I learnt the technology of setting visit flag from Data Structure class taught by Prof.MingZhu when I was an undergraduate.
        if not meccirnas[i].visited:
            meccirnas[i].group = i
            group2mecci[i] = []
            group2mecci[i].append(meccirnas[i])
        else:
            continue

        # the core of this function
        k = i
        j = i + 1
        while j < meccinum:
            start_dis = meccirnas[j].start - meccirnas[k].start
            end_dis = meccirnas[j].end - meccirnas[k].end
            # The heart
            if not meccirnas[j].visited and start_dis == 1 and end_dis == 1 and meccirnas[j].strand == meccirnas[
                    k].strand:
                meccirnas[j].group = i
                group2mecci[i].append(meccirnas[j])
                meccirnas[j].visited = True
                k = j
            j += 1

    mecci_group_list = []
    for groupnum in group2mecci:
        # Initializing of the variables
        chr = ''
        strand = '+'
        read_num = 0
        mecciID = "mecci" + str(groupnum)
        start_list = []
        end_list = []
        reads_list = []
        for mecci in group2mecci[groupnum]:
            chr, strand = mecci.chr, mecci.strand
            start_list.append(mecci.start)
            end_list.append(mecci.end)
            reads_list += mecci.reads_list
        # because of the bug of CIRI2_256 or something about algorothm I cannot understand
        start_list = remove_repeat(start_list)
        end_list = remove_repeat(end_list)
        reads_list = remove_repeat(reads_list)
        read_num = len(reads_list)
        mecci_group_list.append('\t'.join([chr, ','.join([str(i) for i in start_list]), ','.join(
            [str(j) for j in end_list]), mecciID, str(read_num), str(strand), ','.join(reads_list)]))
    return mecci_group_list


# The main function
def main():

    # parse the arguments from command line
    args = args_parser()
    sample = args.sample
    fastq1 = args.fastq1
    fastq2 = args.fastq2
    platform = args.platform
    genome = args.genome
    mito_name = args.mito
    threads = args.threads
    bwa = args.bwa
    bowtie = args.bowtie
    bowtie_build = args.bowtie_build
    bwaindex = args.bwaindex
    seqtk = args.seqtk
    bwa_mismatch = args.mismatch_bwa
    anchor_mismatch = args.mismatch_anchor
    internal_mismatch = args.mismatch_internal
    mitoidx_dir = args.mitoidx_dir
    anchor_size = args.anchor_size
    tmpdir = args.tmpdir
    result_dir = args.result

    # Initialization
    log.info('Initializing and making output directory...')
    mkdir_p(tmpdir)
    mkdir_p(result_dir)
    mkdir_p(mitoidx_dir)

    # Determine the reads length of the RNA-seq dataset
    log.info('Getting the length of reads from fastq file {}'.format(fastq1))
    read_len = get_read_length(fastq1)
    print('The length of the RNA-seq reads is ' + str(read_len))

    # mapping reads using BWA MEM
    log.info('Mapping reads {}, {} against reference genome {} using {}'.format(
        fastq1, fastq2, genome, bwa))
    bwa_samfile = bwamap(sample, bwa, fastq1, fastq2,
                         bwaindex, threads, tmpdir)

    # Getting reads name unmapped from mitochondria
    log.info('Getting unmapped reads name from mitochondria')
    mito_unmapped_reads_list = get_unmapped_reads(
        mito_name, bwa_samfile, read_len, bwa_mismatch)

    # Getting reads from fastq file and then write them into the intermediate files
    log.info('Getting unmapped reads from original fastq file and write them into the temporary fastq files')
    if platform == 'Novo':
        candidate_circ_fastq1, candidate_circ_fastq2 = unmappedreads2fastq_novo(
            sample, seqtk, fastq1, fastq2, mito_unmapped_reads_list, tmpdir)
    
    elif platform == 'BGI':
        candidate_circ_fastq1, candidate_circ_fastq2 = unmappedreads2fastq_BGI(
            sample, seqtk, fastq1, fastq2, mito_unmapped_reads_list, tmpdir)

    else:
        sys.stderr.write('The data generated in this platform is not supported for the time being and I will quit.')
        sys.exit(0)

    # Parsing the reference genome
    log.info('Parsing the reference genome file {}'.format(genome))
    chr2fa = parse_fasta(genome)

    # Prepare the mitochondrial genome
    mito_fasta_file = '{}.fasta'.format(os.path.join(result_dir, mito_name))
    mito_fasta_fh = open(mito_fasta_file, 'w')
    log.info('Getting the sequence of mitochondrial genome and write it to new file {}'.format(
        mito_fasta_file))
    chr2fa[mito_name].write_to_file(mito_fasta_fh)
    mito_fasta_fh.close()

    # Now, let's find mitochondrial-encoded circular RNAs
    log.info('Mapping reads against the mitochondrial genome')
    tmpbed = os.path.join(tmpdir, '{}.tmp.bed'.format(sample))
    tmpbed_sort = os.path.join(tmpdir, '{}.tmp.sort.bed'.format(sample))
    fout = open(tmpbed, 'w')
    mecci_list = find_mecci(mito_name, sample, bowtie, bowtie_build, candidate_circ_fastq1,
                            candidate_circ_fastq2, mito_fasta_file, mitoidx_dir, result_dir, threads, anchor_size, anchor_mismatch, internal_mismatch)

    mecci_list_summary = summary_circ(mecci_list)
    for circ in mecci_list_summary:
        circ.write_to_file(fout)
    fout.close()

    cmd = 'sort -k1,1 -k2,2n {} >{}'.format(tmpbed, tmpbed_sort)
    run_shell_cmd(cmd)

    meccirnas_tmp_list = parse_bed(tmpbed_sort)
    out_prefix = os.path.join(result_dir, sample)
    mecci_output_file = '{}.meccirna.txt'.format(out_prefix)
    mecci_output_fh = open(mecci_output_file, 'w')
    meccirnas = meccirnas_tmp_list
    meccirnas_group = group_mecci(meccirnas)
    for mecci in meccirnas_group:
        mecci_output_fh.write(mecci + '\n')
    mecci_output_fh.close()

    log.info('Congratulations! Pipeline done.')


if __name__ == "__main__":
    main()
