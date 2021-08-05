#!/usr/bin/env python3

import os
import sys
from attrdict import AttrDict
import re
from configparser import ConfigParser
from natsort import natsorted
import subprocess
import shlex

import logging

#============================= Function =================================
##logging info
DEBUG="" #change it when debugging
logFormat = "%(asctime)s [%(levelname)s] %(message)s"
level = "DEBUG" if DEBUG != "" else "INFO"
logging.basicConfig( stream=sys.stderr, level=level, format=logFormat )
#========================================================================

class VCF2rGFAHelper:
    def __init__(self, fasta, vcf, outDir, forceUseAlt=False):
        try:
            import pysam

            if forceUseAlt: raise ModuleNotFoundError()
            self.usePysam = True
        except ModuleNotFoundError:
            from pyfaidx import Fasta

            self.usePysam = False

        self.ENV = {}

        if not self.usePysam:
            if sys.platform == 'win32':
                toolsPath = os.path.join(os.path.dirname(os.path.realpath(__file__)), '..', '..', 'panGraphViewerApp', 'thirdPartyTools', 'SAMtools')
                self.ENV['SAMTOOLS'] = os.path.join(toolsPath, 'samtools.exe')
                self.ENV['BCFTOOLS'] = os.path.join(toolsPath, 'bcftools.exe')
                self.ENV['BGZIP'] = os.path.join(toolsPath, 'bgzip.exe')

            for name in self.ENV:
                if not os.path.isfile(self.ENV[name]):
                    sys.exit(f'ERROR: {self.ENV[name]} not found')

        self.func = self.pysamFunc if self.usePysam else self.altFunc

        self.outDir = outDir
        os.makedirs(self.outDir, mode=0o777, exist_ok=True)

        self.fasta = fasta
        self.fasta_fai = f'{self.fasta}.fai'
        self.vcf = vcf
        self.vcf_sorted = os.path.join(self.outDir, f'{os.path.basename(self.vcf)}_sorted')
        self.vcf_gz = os.path.join(self.outDir, f'{os.path.basename(self.vcf)}.gz')
        self.vcf_gz_tbi = f'{self.vcf_gz}.tbi'

        self.tempFiles = [self.fasta_fai, self.vcf_sorted, self.vcf_gz, self.vcf_gz_tbi]

    def preprocess(self, force=True):
        if force:
            self.removeTempFiles()

        if self.fasta:
            if not os.path.isfile(self.fasta_fai) or os.path.getmtime(self.fasta)>os.path.getmtime(self.fasta_fai):
                self.indexFasta()
        if not os.path.isfile(self.vcf_gz) or os.path.getmtime(self.vcf)>os.path.getmtime(self.vcf_gz):
            self.sortAndCompressVCF()
        if not os.path.isfile(self.vcf_gz_tbi) or os.path.getmtime(self.vcf_gz)>os.path.getmtime(self.vcf_gz_tbi):
            self.indexVCF()

    def getVcfSamples(self):
        vcfSamples = []

        with open(self.vcf) as f:
            for line in f:
                if line[0] != '#': break

                if line.startswith('#CHROM'):
                    fields = line.strip().split('\t')
                    vcfSamples = fields[8:]
                    break

        return vcfSamples

    def runCmdLine(self, cmd, shell=False, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE):
        cmd_list = [str.strip('"').strip("'") for str in shlex.split(cmd, posix=False)]

        ret = subprocess.run(cmd_list, shell=False, stdout=stdout, stderr=stderr, stdin=stdin, text=True)

        if ret.returncode:
            logging.error(f'cmd: {cmd}')
            logging.error(f'stderr: {ret.stderr.strip()}')

        return ret.stdout.strip().split('\n') if ret.stdout else ret.stdout

    def pysamFunc(self, param):
        import pysam

        action = param['action']
        output = ''

        if action == 'sortAndCompressVCF':
            self.sortVcf(self.vcf, self.vcf_sorted)
            pysam.tabix_compress(self.vcf_sorted, self.vcf_gz, force=True)
        elif action == 'indexFasta':
            pysam.faidx(self.fasta, force=True)
        elif action == 'indexVCF':
            pysam.tabix_index(self.vcf_gz, preset='vcf', force=True)
        elif action == 'fetchFasta':
            chr, start, end, env = param['chr'], param['start'], param['end'], param['env']

            if 'pysam' not in env: env['pysam'] = pysam.FastaFile(self.fasta)
            output = env['pysam'].fetch(chr, start-1, end-1)
        elif action == 'fetchVCF':
            chr, start, end = param['chr'], param['start'], param['end']

            vcf = pysam.VariantFile(self.vcf_gz)
            output = vcf.fetch(chr, start, end)
        else:
            print(f'ERROR: invalid action: {action}')
            return None

        return output

    def altFunc(self, param):
        from pyfaidx import Fasta

        action = param['action']
        bgzip, samtools, bcftools = self.ENV['BGZIP'], self.ENV['SAMTOOLS'], self.ENV['BCFTOOLS']
        cmd, output = '', ''

        if action == 'sortAndCompressVCF':
            self.sortVcf(self.vcf, self.vcf_sorted)
            with open(self.vcf_gz, 'wb') as f_vcf_gz:
                cmd = f'{bgzip} -c "{self.vcf_sorted}"'
                self.runCmdLine(cmd, stdout=f_vcf_gz)
        elif action == 'indexFasta':
            cmd = f'{samtools} faidx "{self.fasta}"'
            output = self.runCmdLine(cmd)
        elif action == 'indexVCF':
            cmd = f'{bcftools} index -t "{self.vcf_gz}"'
            output = self.runCmdLine(cmd)
        elif action == 'fetchFasta':
            chr, start, end, env = param['chr'], param['start'], param['end'], param['env']

            if not 'Fasta' in param: param['Fasta'] = Fasta(self.fasta)
            output = param['Fasta'][chr][start-1:end-1]
            output = str(output)
        elif action == 'fetchVCF':
            chr, start, end = param['chr'], param['start'], param['end']

            cmd = f'{bcftools} view "{self.vcf_gz}" {chr}:{start}-{end} -U -H'
            output = self.runCmdLine(cmd)
        else:
            print(f'ERROR: invalid action: {action}')
            return None

        if action == 'fetchVCF':
            vcfColNames = self.getVcfColNames()

            records = []
            for vcfLine in output:
                record = self.formatVcfLine(vcfLine.split('\t'), vcfColNames)
                records.append(record)

            output = records

        return output

    def formatVcfLine(self, rec_raw, colNames):
        row  = [field if i<7 else dict(x.split("=") if '=' in x else [x,1] for x in field.split(";")) if i==7
                     else field.split(':') if i==8 else AttrDict(zip(rec_raw[8].split(':'), field.split(':')))
                     for i,field in enumerate(rec_raw)]

        record = AttrDict()
        for i in range(9):
            record[colNames[i].lower()] = row[i]

        record['pos'] = int(record['pos'])
        record['alts'] = record['alt'].split(',')
        if 'END' in record['info']: record['stop'] = int(record['info']['END'])
        if 'STRLEN' in record['info']: record['info']['STRLEN'] = int(record['info']['STRLEN'])

        record.samples = {colNames[9+i]:field for i,field in enumerate(row[9:])}
        for sampleName in record.samples:
            sample = record.samples[sampleName]
            gt = re.split('/|\|', sample['GT'])
            sample['allele_indices'] = (None,None) if gt == ['.','.'] else (int(min(gt)), int(max(gt)))

        return record

    def sortAndCompressVCF(self):
        param = {'action':'sortAndCompressVCF'}

        return self.func(param)

    def indexFasta(self):
        param = {'action':'indexFasta'}

        return self.func(param)

    def indexVCF(self):
        param = {'action':'indexVCF'}

        return self.func(param)

    # 1-based, exclusive
    def fetchFasta(self, chr, start, end, env):
        param = {'action':'fetchFasta','chr':chr,'start':start,'end':end,'env':env}

        return self.func(param)

    # 1-based, exclusive
    def fetchVCF(self, chr, start, end, n_row = None):
        param = {'action':'fetchVCF','chr':chr,'start':start,'end':end,'n_row':n_row}

        return self.func(param)

    # 1-based, exclusive
    def getRef(self, chr, start_orig, end, env):
        if end < 1 or end < start_orig:
            return ''

        padding = 'N' * (-start_orig)
        start = 1 if start_orig < 1 else start_orig

        if not self.fasta:
            return f'{padding}{"N"*(end-start)}'
        else:
            return f'{padding}{self.fetchFasta(chr, start, end, env).upper()}'

    def getVcfChroms(self):
        vcfChroms = {}
        with open(self.vcf) as f:
            for line in f:
                if line[0] == '#':
                    continue

                fields = line.split('\t')
                vcfChroms[fields[0]] = 1

        return sorted(list(vcfChroms.keys()))

    def getVcfColNames(self):
        found = None
        with open(self.vcf) as f:
            for line in f:
               if line[0] != '#' :
                   break

               if line.startswith('#CHROM'):
                   found = line.strip()
                   break

        return found[1:].split('\t') if found else None

    def sortVcf(self, vcf, vcf_sorted):
        vcfLines = {}
        with open(vcf) as f_in, open(vcf_sorted, 'w') as f_out:
            for line in f_in:
                if line[0] == '#':
                    f_out.write(line)
                else:
                    fields = line.strip().split('\t')
                    vcfLines[f'{fields[0]}_{fields[1]}'] = line
            for entry in natsorted(vcfLines):
                 f_out.write(vcfLines[entry])

    def getRefChromsFromVcf(self):
        refChroms = {}
        with open(self.vcf) as f:
            for line in f:
               if line[0] != '#':
                   break
               #m = re.search('^##contig=<(.*)>', line.strip())
               m = re.search('^##contig=<ID=(.*),length=(.*)>', line.strip())
               if m and m.group(1) and m.group(2):
                   refChroms[m.group(1)] = {'chrom':m.group(1), 'length':int(m.group(2))}

        return refChroms

    def getRefChroms(self):
        fasta_idx = f"{self.fasta}.fai"

        refChroms = {}
        with open(fasta_idx) as f:
            for line in f:
                 row = line.strip().split('\t')
                 refChroms[row[0]] = {'chrom':row[0], 'length':int(row[1])}

        return refChroms

    def removeTempFiles(self):
        for file in self.tempFiles:
            try:
                print(f'removing {file}')
                os.remove(file)
            except:
                pass

    def test(self):
        self.removeTempFiles()
        self.preprocess(force=True)
        for file in self.tempFiles:
            print(file)
            filesize = os.path.getsize(file)
            print(filesize)

if __name__=="__main__":
    from argparse import ArgumentParser

    parser = ArgumentParser(description='helper functions for vcf2rGFA')
    parser.add_argument('-v', dest='vcf', help='the vcf file', type=str)
    parser.add_argument('-o', dest='outDir', help='the output directory', type=str)
    # optional
    parser.add_argument('-f', dest='fasta', help='the fasta file', type=str)

    args = parser.parse_args()

    if None not in [args.vcf, args.outDir]:
        VCF2rGFAHelper(args.fasta, args.vcf, args.outDir, forceUseAlt=True).test()
    else:
        print('\n%s\n' % parser.print_help())
