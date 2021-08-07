#!/usr/bin/env python3

import os
import sys
import shlex
import subprocess 
from attrdict import AttrDict
import re
from configparser import ConfigParser

class VCF2rGFAHelper:
    def __init__(self, fasta, vcf, outDir):
        try:
            import pysam

            self.usePysam = True
        except ModuleNotFoundError:
            from pyfaidx import Fasta

            self.usePysam = False

        self.ENV = {}

        if not self.usePysam:
            if sys.platform == 'win32':
                self.ENV['SAMTOOLS'] = os.path.join(os.path.dirname(os.path.realpath(__file__)), '..', 'thirdPartyTools', 'SAMtools', 'samtools.exe')
                self.ENV['BCFTOOLS'] = os.path.join(os.path.dirname(os.path.realpath(__file__)),'..', 'thirdPartyTools', 'SAMtools', 'bcftools.exe')
                self.ENV['BGZIP'] = os.path.join(os.path.dirname(os.path.realpath(__file__)), '..', 'thirdPartyTools', 'SAMtools','bgzip.exe')
                self.samtools = self.ENV['SAMTOOLS']
                self.bcftools = self.ENV['BCFTOOLS']
                self.bgzip = self.ENV['BGZIP']
            for name in self.ENV:
                if not os.path.isfile(self.ENV[name]):
                    sys.exit(f'ERROR: {self.ENV[name]} not found')

        self.func = self.pysamFunc if self.usePysam else self.altFunc

        self.outDir = outDir
        try:
            os.makedirs(self.outDir, mode=0o777, exist_ok=True)
        except:
            pass

        self.fasta = fasta
        #self.fasta_fai = f'{self.outDir}/{os.path.basename(self.fasta)}.fai'
        self.fasta_fai = f'{self.fasta}.fai'
        self.vcf = vcf
        self.vcf_sorted = os.path.join(self.outDir, f'{os.path.basename(self.vcf)}_sorted')
        self.vcf_gz = os.path.join(self.outDir, f'{os.path.basename(self.vcf)}.gz')
        self.vcf_gz_tbi = os.path.join(self.outDir, f'{self.vcf_gz}.tbi')
        self.tempFiles = []

        if not self.usePysam:
            cmd = f'"{self.bcftools}" view "{self.vcf}" -h | grep ^#CHROM'
            self.vcfHeader = ['CHROM']+self.runCmdLine(shlex.split(cmd))[-1].split('\t')[1:]

    def preprocess(self, force=False):
        if self.fasta:
            if not os.path.isfile(self.fasta_fai) or os.path.getmtime(self.fasta)>os.path.getmtime(self.fasta_fai):
                self.tempFiles.append(self.fasta_fai)
                self.indexFasta()
        if not os.path.isfile(self.vcf_gz) or os.path.getmtime(self.vcf)>os.path.getmtime(self.vcf_gz):
            self.tempFiles.append(self.vcf_sorted)
            self.tempFiles.append(self.vcf_gz)
            self.sortAndCompressVCF()
        if not os.path.isfile(self.vcf_gz_tbi) or os.path.getmtime(self.vcf_gz)>os.path.getmtime(self.vcf_gz_tbi):
            self.tempFiles.append(self.vcf_gz_tbi)
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

    def runCmdLine(self, cmd, shell=True, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE):
        ret = subprocess.run(cmd, shell=True, stdout=stdout, stderr=stderr, stdin=stdin, text=True)
        return ret.stdout.strip().split('\n') if ret.stdout else ret.stdout

    def pysamFunc(self, param):
        import pysam

        action = param['action']
        output = ''

        if action == 'sortAndCompressVCF':
            # should run both sort and compress together. to-be-fixed

            cmd = f'grep "#" "{self.vcf}" > "{self.vcf_sorted}"'
            os.system(cmd)
            cmd = f'grep -v "#" "{self.vcf}" | sort -k1V -k2 >> "{self.vcf_sorted}"'
            os.system(cmd)
           
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
            cmd = f'(grep # "{self.vcf}" && grep -v # "{self.vcf}" | sort -k1V -k2) | "{bgzip}" > "{self.vcf_gz}"'
            output = subprocess.run(shlex.split(cmd), shell=True)

        elif action == 'indexFasta':
            cmd = f'"{samtools}" faidx "{self.fasta}"'
            output = subprocess.run(shlex.split(cmd), shell=True)
        elif action == 'indexVCF':
            cmd = f'"{bcftools}" index -t "{self.vcf_gz}"'
            output = subprocess.run(shlex.split(cmd), shell=True)
        elif action == 'fetchFasta':
            chr, start, end, env = param['chr'], param['start'], param['end'], param['env']

            if not 'Fasta' in param: param['Fasta'] = Fasta(self.fasta)
            output = param['Fasta'][chr][start-1:end-1]
            output = str(output)
        elif action == 'fetchVCF':
            chr, start, end = param['chr'], param['start'], param['end']

            cmd = f'"{bcftools}" view "{self.vcf_gz}" "{chr}:{start}-{end}" -U -H'
            output = self.runCmdLine(shlex.split(cmd))
        else:
            print(f'ERROR: invalid action: {action}')
            return None

        if action == 'fetchVCF':
            records = []
            for vcfLine in output:
                record = self.formatVcfLine(vcfLine.split('\t'))
                records.append(record)

            output = records

        return output

    def formatVcfLine(self, rec_raw):
        row  = [field if i<7 else dict(x.split("=") if '=' in x else [x,1] for x in field.split(";")) if i==7
                     else field.split(':') if i==8 else AttrDict(zip(rec_raw[8].split(':'), field.split(':')))
                     for i,field in enumerate(rec_raw)]

        record = AttrDict()
        for i in range(9):
            record[self.vcfHeader[i].lower()] = row[i]

        record['pos'] = int(record['pos'])
        record['alts'] = record['alt'].split(',')
        if 'END' in record['info']: record['stop'] = int(record['info']['END'])
        if 'STRLEN' in record['info']: record['info']['STRLEN'] = int(record['info']['STRLEN'])

        record.samples = {self.vcfHeader[9+i]:field for i,field in enumerate(row[9:])}
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
        if sys.platform == "win32":
            cmd = f'cat "{self.vcf}" | grep -v # | awk "{{print $1}}" | sort | uniq'
            vcfChroms = self.runCmdLine(shlex.split(cmd))
        else:
            cmd = f'cat "{self.vcf}" | grep -v "#" | awk \'{{print $1}}\' | sort | uniq > "{self.vcf}.chrs"'
            os.system(cmd)
            vcfChroms = []
            with open(f"{self.vcf}.chrs", 'r') as chrs:
                for line in chrs:
                    line = line.strip()
                    vcfChroms.append(line)
            os.remove(f"{self.vcf}.chrs")
        return vcfChroms

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
                os.remove(file)
            except:
                pass
