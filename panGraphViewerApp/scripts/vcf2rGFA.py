#!/usr/bin/env python3

import os
import sys
import logging
import re

from argparse import ArgumentParser
from multiprocessing import Process, Pool, Manager
from natsort import natsorted, ns

try:
    from scripts.vcf2rGFAHelper import *
except ModuleNotFoundError:
    from vcf2rGFAHelper import *

import json

#============================= Function =================================
##logging info
DEBUG="" #change it when debugging
logFormat = "%(asctime)s [%(levelname)s] %(message)s"
level = "DEBUG" if DEBUG != "" else "INFO"
logging.basicConfig( stream=sys.stderr, level=level, format=logFormat )
#========================================================================

config = None
script_directory = os.path.dirname(os.path.realpath(__file__))
copied = os.path.join(script_directory, '..',  "config.ini")

def get_var(cfg, group, var, must_have=False):
    global config

    if not config:
        config = ConfigParser()
        config.read(cfg)

    try:
        return config.get(group, var)
    except:
        if must_have:
            raise
        else:
            logging.info(f'Config value [{group}][{var}] not found')
        return None

class VCF2rGFA:
    def __init__(self, options):
        self.vcf = options.vcf
        self.fasta = options.fasta
        self.backbone = options.backbone
        self.outDir = options.outdir
        self.nthread = options.nthread
        self.convertAllChr = 0 if options.chr else 1
        self.prefix = options.prefix if options.prefix else f'{options.backbone}_vcf2rGFA'

        try:
            os.makedirs(self.outDir, mode=0o777, exist_ok=True)
        except:
            pass

        self.nodeID = {}

        self.vcf2rGFAHelper = VCF2rGFAHelper(self.fasta, self.vcf, self.outDir)
        self.vcfSamples = self.vcf2rGFAHelper.getVcfSamples()
        if not self.vcfSamples:
            logging.error('No sample columns found in VCF. Program aborted')
            sys.exit(1)

        if not self.fasta:
            logging.info(f"As reference fasta is not provided, sequences will not be available in the output rGFA")

        logging.info(f'Preprocessing ...')
        self.vcf2rGFAHelper.preprocess(force=True)
        logging.info(f'Preprocessing done')

        self.vcfChroms = self.vcf2rGFAHelper.getVcfChroms()
        self.refChroms = self.vcf2rGFAHelper.getRefChroms() if self.fasta else self.vcf2rGFAHelper.getRefChromsFromVcf()
        missChroms = sorted(list(set(self.vcfChroms) - set(self.refChroms)))
        if missChroms:
            if self.fasta:
                logging.error(f"Some chromosomes in VCF are missing in REF: {', '.join(missChroms)}")
                sys.exit(2)
            else:
                logging.error(f"Please provide fasta reference as some contig info are missing in VCF headers: {', '.join(missChroms)}")
                sys.exit(3)

        commChroms = sorted(list(set(self.refChroms).intersection(self.vcfChroms)))
        logging.info(f"Chromosomes in VCF to be converted: {', '.join(commChroms)}")

        self.chroms = options.chr if options.chr else self.vcfChroms
        self.SN_delim = get_var(copied, 'nodes','SN_delim', must_have=True)

    def getNodeID(self, chr, nodeType):
        try:
            nodeID = self.nodeID[chr][nodeType]
        except:
            if chr not in self.nodeID:
                self.nodeID[chr] = {}
            if nodeType not in self.nodeID[chr]:
                self.nodeID[chr][nodeType] = 0

            nodeID = self.nodeID[chr][nodeType]

        self.nodeID[chr][nodeType] += 1

        return f'::{chr}_{nodeType}_{str(nodeID)}:::'

    def getNodeIDMap(self, workerResults):
        results2 = {}
        sampleNameDict = {}
        nodeIDMap = {}

        currNodeID = 1
        for result in workerResults:
            chr = result['chr']
            results2[chr] = result

            for sampleName in result['nodeIDInfo']:
                sampleNameDict[sampleName] = 1

        sampleNameList = sorted(sampleNameDict.keys())

        for sampleName in sampleNameList:
            for chr in self.chroms:
                if sampleName in results2[chr]['nodeIDInfo']:
                    for i in range(results2[chr]['nodeIDInfo'][sampleName]):
                       nodeIDMap[f'::{chr}_{sampleName}_{i}:::'] = f's{currNodeID}'
                       currNodeID += 1

        return nodeIDMap

    def save2GFA(self, workerResults):
        logging.info(f"Generating rGFA ...")

        gfa = os.path.join(self.outDir, f'{self.prefix}.gfa')
        if os.path.isfile(gfa):
            os.remove(gfa)

        nodeIDMap = self.getNodeIDMap(workerResults)

        with open(gfa, 'w') as fd:
            for result in workerResults:
                chr = result['chr']
                nodes = result['nodes']
                edges = result['edges']
                nodeIDInfo = result['nodeIDInfo']

                for nodeID in nodes:
                    # get updated nodeID
                    nodeID2 = nodeIDMap[nodeID]
                    node = nodes[nodeID]
                    node['INF'] = self.updateNodeID(nodeIDMap, node['INF'])
                    if 'DEL_SAMPLE' in node and node['DEL_SAMPLE']:
                       node['INF'] = f"INF:Z:SV||DEL||{'||'.join(sorted(node['DEL_SAMPLE']))}"
                       print('check 3', node['INF'])
                    output = ['S',nodeID2,node['seq'],node['LN'],node['SN'],node['SO'],node['SR'],node['INF']]
                    fd.write('%s\n' % "\t".join(output))

                for i in edges:
                    fromNodeID = nodeIDMap[i[0]]
                    toNodeID = nodeIDMap[i[1]]
                    fd.write(f'L\t{fromNodeID}\t{edges[i][0]}\t{toNodeID}\t{edges[i][1]}\t0M\n')

        os.chmod(gfa, 0o666)

        logging.info(f"Complete generating rGFA {gfa}")

    def worker(self, dummy):
        results = []
        # per-worker env (if any)
        env = {}

        try:
           while True:
              chr = self.d.pop(0)
              results.append(self.convertVCF2rGFA(chr, env=env))
        except IndexError as indexError:
           pass

        return results

    def rev_comp(self, seq):
        trans = str.maketrans('ACGT', 'TGCA')
        return seq.translate(trans)[::-1]

    def updateNodeID(self, nodeIDMap, input):
        output = input
        p = re.compile('::(.*?):::')
        found = p.findall(input)
        for s in found:
           updatedNodeID = nodeIDMap[f'::{s}:::']
           output = output.replace(f'::{s}:::', updatedNodeID)

        return output

    def convertVCF2rGFA(self, chr, env):
        logging.info(f"Start to convert the VCF file to rGFA file at chr: {chr} ...")

        try:
            chrLen = self.refChroms[chr]['length']
        except KeyError:
            logging.error(f"Illegal chromosome name '{chr}'. Abort!")
            return -3

        nodes, edges, backbone = {}, {}, self.backbone

        # svList[sampleName][varStart] = variants
        svList = {}

        # anchors in backbone for sample nodes
        # (e.g. posDict[8]['end'] mean anchor ends at 8, posDict[9]['start'] means anchor atarts at 9)
        posDict = {}

        for rec in self.vcf2rGFAHelper.fetchVCF(chr, 1, chrLen + 1):
            varType, varStart, varEnd, varSeq, refSeq, varLen, posShift1, isReverse = None, None, None, None, None, None, False, False

            # varStart, varEnd: exclusive
            # e.g. if pos = 8, varStart = 7 (before the start of var)
            varId = rec.id
            varStart = rec.pos - 1
            if 'SVTYPE' in rec.info:
                varType = rec.info['SVTYPE']

                # don't know why it sometimes returns tuple. to-be-fixed
                #tmp = rec.info['SVLEN'][0] if type(rec.info['SVLEN']) is tuple else rec.info['SVLEN']
                #varLen = int(abs(tmp))

                # cannot use info['END'] in pysam as 'END' is reserved attribute
                if varType != 'BND':
                    varEnd = rec.stop
                    varLen = varEnd - varStart - 1

                if varType in ['DUP']:
                    refSeq = ''
                    varSeq = self.vcf2rGFAHelper.getRef(chr, varStart+1, varEnd, env)

                    # assume the SV is inserted at the end of the duplicated region
                    #varStart = varStart + varLen
                elif varType in ['DEL']:
                    refSeq = self.vcf2rGFAHelper.getRef(chr, varStart+1, varEnd, env)
                    varSeq = refSeq
                elif varType in ['INV']:
                    refSeq = self.vcf2rGFAHelper.getRef(chr, varStart+1, varEnd, env)
                    varSeq = self.rev_comp(refSeq)
                    isReverse = True
                elif varType in ['INS']:
                    refSeq = ''
                    varSeq = rec.alts[0]
                    varLen = len(varSeq)
                elif varType in ['BND']:
                    refSeq = ''

                    alt = rec.alts[0]
                    p = re.compile('](.*?)]|\[(.*?)\[')
                    found = p.findall(alt)
                    if found[0][0] or found[0][1]:
                        bpChr, bpPos = found[0][0].split(':') if found[0][0] else found[0][1].split(':')
                        bpPos = int(bpPos)
                        varLen = -1

                        seqLen = 100
                        if found[0][0]:
                            if alt.startswith(']'): # ]p]t: piece extending to the left of p is joined before t
                                varStart = None
                                varEnd = bpPos
                                varSeq = self.vcf2rGFAHelper.getRef(bpChr, varEnd - seqLen, varEnd, env)
                            else:                   # t]p]: reverse comp piece extending left of p is joined after t
                                varStart = bpPos
                                varEnd = None
                                varSeq = self.rev_comp(self.vcf2rGFAHelper.getRef(bpChr, varStart - seqLen, varStart, env))
                                isReverse = True
                        else:
                            if alt.startswith('['): # [p[t: reverse comp piece extending right of p is joined before t
                                varStart = None
                                varEnd = bpPos
                                varSeq = self.rev_comp(self.vcf2rGFAHelper.getRef(bpChr, varEnd, varEnd + seqLen, env))
                                isReverse = True
                            else:                   # t[p[: piece extending to the right of p is joined after t
                                varStart = bpPos
                                varEnd = None
                                varSeq = self.vcf2rGFAHelper.getRef(bpChr, varStart, varStart + seqLen, env)
                    else:
                        logging.warning(f"Invalid format for BND: {rec.alts[0]}. Variant ignored")
                        continue
                else:
                    logging.warning(f"SVTYPE ignored: {varType}")
                    continue
            else:
                ref = rec.ref
                alt = rec.alts[0]
                varType = 'SNP' if len(ref)==len(alt) else 'INS' if len(ref)<len(alt) else 'DEL'
                varSeq = None

                if rec.ref == 'N':
                    refSeq = ''
                elif len(ref)>=1 and len(alt)>=1 and ref[0] == alt[0][0]:
                    posShift1 = True
                    varStart += 1
                    refSeq = rec.ref[1:]
                else:
                    refSeq = rec.ref

            """
            for sampleName in rec.samples:
                sample = rec.samples[sampleName]
                #if sample.allele_indices != (0,0) and sample.allele_indices != (None, None):
                if any(sample.allele_indices):
                    if varSeq is None:
                        varSeq = rec.alts[max(sample.allele_indices)-1]
            """

            sampleList = []
            for sampleName in rec.samples:
                if any(rec.samples[sampleName].allele_indices):
                    sampleList.append(sampleName)

            if True:
                if sampleList:
                    #sampleName = sampleList[0]
                    sampleName = '~'.join(sampleList)

                    if varSeq is None:
                        varSeq = rec.alts[0]

                    if posShift1:
                        varSeq = varSeq[1:]

                    if varEnd is None and varType != 'BND':
                        varEnd = varStart + len(refSeq) + 1

                    if varType != 'BND' and ((varStart and varStart+1 > chrLen) or (varEnd and varEnd-1 > chrLen)):
                        tempStart  = varStart + 1 - varLen if varType == 'DUP_old' else varStart + 1
                        tempEnd = varEnd - 1
                        logging.warning(f"{varType} ({tempStart}-{tempEnd}) for {sampleName} is outside the range of the chr. Variant is ignored")
                        continue

                    if varType == 'DUP_old':
                        varDesc = f'SV||{varType}||{varStart+1-len(varSeq)}||{varStart}'

                        # split duplicated backbone node
                        if varStart-len(varSeq) not in posDict: posDict[varStart-len(varSeq)] = {}
                        posDict[varStart-len(varSeq)]['end'] = 1
                    elif varType == 'INS':
                        #varDesc = f'SV||{varType}||{varStart+1}'
                        varDesc = f'SV||{varType}||{chr}||{varStart+1}||{chr}||{varEnd}||{varId}||{sampleName}'
                    elif varType == 'BND':
                        varDesc = f"SV||{varType}||{rec.alts[0].replace(':','||')}"
                    else:
                        #varDesc = f'SV||{varType}||{varStart+1}||{varEnd-1}'
                        varDesc = f'SV||{varType}||{chr}||{varStart+1}||{chr}||{varEnd}||{varId}||{sampleName}'

                    if sampleName not in svList: svList[sampleName] = {}
                    if varStart not in svList[sampleName]: svList[sampleName][varStart] = []
                    svList[sampleName][varStart].append({'sampleName':sampleName,'varDesc':varDesc,'varType':varType,'varStart':varStart,'varEnd':varEnd,'varLen':varLen,'varSeq':varSeq,'refSeq':refSeq,'isReverse':isReverse})

                    if varStart:
                        if varStart not in posDict: posDict[varStart] = {}
                        posDict[varStart]['end'] = 1
                    if varEnd:
                        if varEnd not in posDict: posDict[varEnd] = {}
                        posDict[varEnd]['start'] = 1

        # init backbone node
        nodeID = self.getNodeID(chr=chr, nodeType='0backbone')
        start = 1

        # create backbone nodes and add to global list
        pos2NodeID = {}
        prevNodeID = None
        for pos in sorted(posDict.keys()):
            if 'end' in posDict[pos] and not 'start' in posDict[pos]:
                end = pos
                seq = self.vcf2rGFAHelper.getRef(chr, start, end+1, env)
                nodes[nodeID] = {'seq':seq,'LN':f'LN:i:{len(seq)}','SN':f'SN:Z:{backbone}{self.SN_delim}{chr}','SO':f'SO:i:{start-1}','SR':'SR:i:0','INF':''}

                pos2NodeID[start] = nodeID
                pos2NodeID[end] = nodeID

                posDict[pos]['end'] = nodeID
                prevNodeID = nodeID

                nodeID = self.getNodeID(chr=chr, nodeType='0backbone')
                start = pos + 1

                edges[(prevNodeID, nodeID)] = ('+', '+')
            elif not 'end' in posDict[pos] and 'start' in posDict[pos]:
                if start != pos:
                    end = pos -1
                    seq = self.vcf2rGFAHelper.getRef(chr, start, end+1, env)
                    nodes[nodeID] = {'seq':seq,'LN':f'LN:i:{len(seq)}','SN':f'SN:Z:{backbone}{self.SN_delim}{chr}','SO':f'SO:i:{start-1}','SR':'SR:i:0','INF':''}

                    pos2NodeID[start] = nodeID
                    pos2NodeID[end] = nodeID

                    prevNodeID = nodeID

                    nodeID = self.getNodeID(chr=chr, nodeType='0backbone')
                    start = pos
                    posDict[pos]['start'] = nodeID

                    edges[(prevNodeID, nodeID)] = ('+', '+')
                else:
                    posDict[pos]['start'] = nodeID
            elif 'end' in posDict[pos] and 'start' in posDict[pos]:
                if start == pos:
                    posDict[pos]['start'] = nodeID
                    posDict[pos]['end'] = nodeID
                    nodes[nodeID] = {'seq':seq,'LN':f'LN:i:{len(seq)}','SN':f'SN:Z:{backbone}{self.SN_delim}{chr}','SO':f'SO:i:{start-1}','SR':'SR:i:0','INF':''}

                    pos2NodeID[start] = nodeID
                    pos2NodeID[end] = nodeID

                    prevNodeID = nodeID

                    nodeID = self.getNodeID(chr=chr, nodeType='0backbone')
                    start = pos + 1

                    edges[(prevNodeID, nodeID)] = ('+', '+')
                else:
                    end = pos - 1
                    seq = self.vcf2rGFAHelper.getRef(chr, start, end+1, env)
                    nodes[nodeID] = {'seq':seq,'LN':f'LN:i:{len(seq)}','SN':f'SN:Z:{backbone}{self.SN_delim}{chr}','SO':f'SO:i:{start-1}','SR':'SR:i:0','INF':''}

                    pos2NodeID[start] = nodeID
                    pos2NodeID[end] = nodeID

                    prevNodeID = nodeID

                    nodeID = self.getNodeID(chr=chr, nodeType='0backbone')
                    start = pos
                    end = pos
                    seq = self.vcf2rGFAHelper.getRef(chr, start, end+1, env)
                    nodes[nodeID] = {'seq':seq,'LN':f'LN:i:{len(seq)}','SN':f'SN:Z:{backbone}{self.SN_delim}{chr}','SO':f'SO:i:{start-1}','SR':'SR:i:0','INF':''}

                    pos2NodeID[start] = nodeID
                    pos2NodeID[end] = nodeID

                    posDict[pos]['end'] = nodeID
                    posDict[pos]['start'] = nodeID

                    edges[(prevNodeID, nodeID)] = ('+', '+')

                    prevNodeID = nodeID

                    nodeID = self.getNodeID(chr=chr, nodeType='0backbone')

                    start = pos + 1

                    edges[(prevNodeID, nodeID)] = ('+', '+')

        end = self.refChroms[chr]['length']+1
        seq = self.vcf2rGFAHelper.getRef(chr, start, end, env)
        nodes[nodeID] = {'seq':seq,'LN':f'LN:i:{len(seq)}','SN':f'SN:Z:{backbone}{self.SN_delim}{chr}','SO':f'SO:i:{start-1}','SR':'SR:i:0','INF':''}

        pos2NodeID[start] = nodeID
        pos2NodeID[end] = nodeID

        # create sample nodes and add to global list (without edges)
        sampleNames = sorted(list(svList.keys()))
        for sampleName in svList:
            for pos in svList[sampleName]:
                for sv in svList[sampleName][pos]:
                    svType = sv['varType']
                    pos1 = sv['varStart']
                    pos2 = sv['varEnd']
                    seq = sv['varSeq']
                    sampleName = sv['sampleName']
                    varDesc = sv['varDesc']
                    isReverse = sv['isReverse']

                    sv['edges'] = []
                    if svType == 'DEL_old':
                        sv['edges'].append((posDict[pos1]['end'], posDict[pos2]['start'], '+', '+'))
                        #nodes[pos2NodeID[pos1+1]][5].append(sampleName)
                        if 'DEL_SAMPLE' not in nodes[pos2NodeID[pos1+1]]:
                            nodes[pos2NodeID[pos1+1]]['DEL_SAMPLE'] = []
                        nodes[pos2NodeID[pos1+1]]['DEL_SAMPLE'].append(sampleName)
                    else:
                        nodeID = self.getNodeID(chr=chr, nodeType=sv['sampleName'])
                        # sample node

                        #SR_val = sampleNames.index(sampleName) + 1
                        SR_val = 1

                        size = len(seq)
                        #size = 123

                        #nodes[nodeID] = {'seq':seq,'LN':f'LN:i:{len(seq)}','SN':f'SN:Z:{sampleName}{self.SN_delim}{chr}','SO':f'SO:i:{len(seq)}','SR':f'SR:i:{SR_val}','INF':''}
                        nodes[nodeID] = {'seq':seq,'LN':f'LN:i:{size}','SN':f'SN:Z:{backbone}{self.SN_delim}{chr}','SO':f'SO:i:{pos1}','SR':f'SR:i:{SR_val}','INF':''}

                        if svType != 'INS':
                            # add ref nodeID info
                            if svType == 'DUP_old':
                                startRefNodeID = pos2NodeID[pos1+1-len(seq)]
                                endRefNodeID = pos2NodeID[pos1]
                                nodes[nodeID]['INF'] = f'INF:Z:{varDesc}_REF_{startRefNodeID}_{endRefNodeID}'
                            elif svType == 'BND':
                                   nodes[nodeID]['INF'] = f'INF:Z:{varDesc}'
                            else:
                                startRefNodeID = pos2NodeID[pos1+1]
                                endRefNodeID = pos2NodeID[pos2-1]
                                #nodes[nodeID]['INF'] = f'INF:Z:{varDesc}_REF_{startRefNodeID}_{endRefNodeID}'
                                nodes[nodeID]['INF'] = f'INF:Z:{varDesc}'
                        else:
                            nodes[nodeID]['INF'] = f'INF:Z:{varDesc}'

                        sv['nodeID'] = nodeID
                        if pos1: sv['edges'].append((posDict[pos1]['end'], nodeID, '+', '+' if not isReverse else '-'))
                        if pos2: sv['edges'].append((nodeID, posDict[pos2]['start'], '+' if not isReverse else '-', '+'))

        # update edges
        for sampleName in svList:
            for pos in svList[sampleName]:
                for sv in svList[sampleName][pos]:
                    if sv['varType'] == 'SNP' and len(sv['varSeq']) == 1 and pos-1 in svList[sampleName]:
                        for sv2 in svList[sampleName][pos-1]:
                            if sv2['varType'] == 'SNP' and len(sv2['varSeq']) == 1:
                                sv2['edges'][1] = (sv2['nodeID'], sv['nodeID'], '+', '+')
                                sv['edges'][0] = None

        # add edges to global list
        for sampleName in svList:
            for pos in svList[sampleName]:
                svs = svList[sampleName][pos]
                for sv in svs:
                    for edge in sv['edges']:
                        if edge:
                            edges[edge[0:2]] = edge[2:4]

        logging.info(f"Finish to convert the VCF file to rGFA file at chr: {chr} ...")

        output = {'chr': chr, 'nodes': nodes, 'edges': edges, 'nodeIDInfo': self.nodeID[chr]}
        return output

    def run(self):
        manager = Manager()
        self.d = manager.list()
        for chrom in self.chroms:
           self.d.append(chrom)

        pool = Pool(self.nthread)
        worker_output = pool.map(self.worker, [''] * self.nthread)
        pool.close()
        pool.join()

        results = []
        for output in worker_output:
           results.extend(output)
        try:
            results = natsorted(results, key=lambda x:x['chr'])
        except TypeError:
            sys.exit(4)

        # hardcode for testing
        #results = [self.convertVCF2rGFA(self.chroms[0], env={})]

        self.save2GFA(results)
        self.vcf2rGFAHelper.removeTempFiles()

if __name__=="__main__":
    parser = ArgumentParser(description='Convert a vcf file to an rGFA file')
    parser.add_argument('--version', action='version', version='1.0')
    parser.add_argument('-b', dest='backbone', help='the name of the backbone sample', type = str)
    parser.add_argument('-v', dest='vcf', help='the vcf file', type=str)
    parser.add_argument('-o', dest='outdir', help='the output directory', type=str)
    # optional
    parser.add_argument('-p', dest='prefix', help='output filename prefix', type = str)
    parser.add_argument('-f', dest='fasta', help='a fasta format file that from the backbone sample', type = str)
    parser.add_argument('-c', dest='chr', nargs='*', help='the name of the chromosome(s) [default: all chroms]', type=str)
    parser.add_argument('-n', dest='nthread', help='number of threads [default: 4]', type=int, default = 4)

    args = parser.parse_args()

    if None not in [args.backbone, args.vcf, args.outdir]:
        VCF2rGFA(args).run()
    else:
        print('\n%s\n' % parser.print_help())
