from argparse import ArgumentParser
import networkx as nx
import sys
import math
import logging
from subprocess import Popen, PIPE

from os import path
import os
from pathlib import Path
from configparser import ConfigParser

from random import choice
import json
from natsort import natsorted, ns

from bokeh.embed import file_html
from bokeh.resources import CDN
from bokeh.plotting._tools import process_tools_arg
from dna_features_viewer import GraphicFeature, GraphicRecord

from shutil import copyfile

import datetime
from bisect import bisect_left
import shlex

import subprocess
import re

import zipfile

try:
    from scripts.gfa2rGFA import *
    from scripts.utilities import *
except ModuleNotFoundError:
    try:
        from gfa2rGFA import *
        from utilities import *
    except ModuleNotFoundError:
        from pangraphviewer.utilities import *

from enum import IntEnum
class NODE(IntEnum):
    seqDesc = 0
    seqLastDesc = 1
    len = 2
    seq = 3
    sample = 4
    chr = 5
    lenBefore = 6
    rank = 7
    inf = 8

    posStart = 9
    posEnd = 10

class EDGE(IntEnum):
    #type = 0
    fromNodeId = 0
    fromStrand = 1
    toNodeId = 2
    toStrand = 3

class RAWNODEDATA(IntEnum):
    fileOffset = 0
    lineLen = 1
    #sample = 2
    rank = 2
    contig = 3
    lenBefore = 4
    seqLen = 5

class DrawGraphResult:
    @staticmethod
    def toJson(panGraph, drawGraphResult):
        output = {}
        output['plotErrorSignal'] = panGraph.plotErrorSignal
        output['emptyGraphSignal'] = panGraph.emptyGraphSignal
        output['overNodeLimitSignal'] = panGraph.overNodeLimitSignal
        output['subNodesCount'] = len(panGraph.subNodes)
        output['outHtml'] = drawGraphResult['outHtml']

        return json.dumps(output)

    def fromJson(panGraph, resultJson):
        obj = json.loads(resultJson)
        panGraph.plotErrorSignal = obj['plotErrorSignal']
        panGraph.emptyGraphSignal = obj['emptyGraphSignal']
        panGraph.overNodeLimitSignal = obj['overNodeLimitSignal']
        panGraph.subNodesCount = obj['subNodesCount']
        panGraph.drawGraphResult = {'outHtml':obj['outHtml']}

class FormatException(Exception):
    def __init__(self, message="Format error occurred"):
        super().__init__(message)
        self.custom_data = "Optional additional data"

class PanGraph:
    seqDescLen = 10
    SN_delim = getVar(copied, 'nodes', 'SN_delim', mustHave=True)
    maxNodesLimit = int(getVar(copied, 'nodes', 'maxNodesLimit', mustHave=True))
    maxNodesDisplay = int(getVar(copied, "nodes", "maxNodesDisplay"))

    def __init__(self, gfa, outdir, parseRGFA=True, nodeIdDict=None):
        self.gfa = gfa
        self.outdir = outdir

        try:
            os.makedirs(self.outdir, mode=0o777, exist_ok=True)
        except:
            pass

        self.nodes = {}
        self.edges = []
        self.inf = None
        self.nameCols = {}
        self.config = None
        self.backbone = None
        self.nodesInfo = None
        self.G = None

        self.illegalrGFA = 0
        self.plotErrorSignal = 0
        self.overNodeLimitSignal = 0
        self.emptyGraphSignal = 0
        self.noOverlap = False

        self.error_unknown = 0
        self.drawgraph_by_nodeids = False

        self.subnodes = -1
        self.exceedNodeIdCount = False
        self.maxNodeCount = 20000

        # init
        if parseRGFA:
            logging.info("Parsing rGFA")
            self.inf = self.parseRGFA(nodeIdDict)
            self.neededGFA = self.inf['neededGFA']
            self.backbone_info = self.inf['backbone']
            self.colorPalettes()

        self.bigNodeSeqLen = int(getVar(copied, 'nodes', 'maxNodeLenDisplay'))
        self.bigNodeList = []

    def runCmdline(self, cmd, shell=True, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE):
        ret = subprocess.run(cmd, shell=True, stdout=stdout, stderr=stderr, stdin=stdin, text=True)

        return ret.stdout.strip().split('\n') if ret.stdout else ret.stdout

    # assign color
    def colorPalettes(self):
        self.nameCols = {s:None for s in self.inf['samples']}
        self.num = len(self.nameCols)
        cols = [f'#{s}' for s in getVar(copied, 'nodes', 'colors').replace('\n','').split(',')]
        #seed(1)
        sequence = [i for i in range(len(cols))]
        self.cols = list()
        for _ in range(self.num):
            select = choice(sequence)
            self.cols.append(cols[select])

        i = 0
        for key in self.nameCols.keys():
            self.nameCols[key] = self.cols[i]
            i += 1

    # get backbone info and sample names
    def parseRGFA(self, nodeIdDict):
        samples = {}
        neededGFA = False
        #backbone = {'name':None, 'contigs':{}}
        backbone = {}
        nodeIdDDlist = []

        rawNodeData = {}
        rawEdgeData = []
        fileOffset = 0

        with open(self.gfa) as f:
            try:
                for lineNum, line in enumerate(f):
                    if line[0] == 'S':
                        row = line.strip().split('\t')

                        nodeId = row[1]

                        if nodeIdDict and nodeId not in nodeIdDict:
                            fileOffset += len(line)
                            continue

                        if lineNum <= self.maxNodeCount:
                            nodeIdDDlist.append(nodeId)
                        else:
                            self.exceedNodeIdCount = True

                        tags = {}
                        for val in row[3:]:
                            lst = val.split(':')
                            tags[lst[0]] = lst[2]

                        seqLen = int(tags['LN']) if 'LN' in tags else len(row[2])
                        rank = tags['SR'] if 'SR' in tags else ''
                        res = tags['SN'] if 'SN' in tags else ''
                        if len(res.split(self.SN_delim)) >= 2:
                            sample = res.split(self.SN_delim)[0]
                            contig = self.SN_delim.join(res.split(self.SN_delim)[1:])
                        else:
                            #neededGFA = True
                            #sample = rank
                            #contig = res
                            raise FormatException()
                        lenBefore = int(tags['SO']) if 'SO' in tags else 0

                        samples[sample] = 1
                        if rank == '0':
                            if sample not in backbone:
                                backbone[sample] = {'name':'', 'contigs':{}}
                            backbone[sample]['contigs'][contig] = 1
                            backbone[sample]['name'] = sample

                        # not needed (TBC)
                        #rawNodeData[nodeId] = [fileOffset, len(line), rank, contig, lenBefore, seqLen]
                    elif line[0] == 'L':
                        row = line.split()
                        fromNodeId, fromNodeStrand, toNodeId, toNodeStrand = row[1:5]

                        if nodeIdDict and (fromNodeId not in nodeIdDict or toNodeId not in nodeIdDict):
                            continue

                        # not needed (TBC)
                        #str = f'{fromNodeId}|{fromNodeStrand}|{toNodeId}|{toNodeStrand}'
                        #rawEdgeData.append(str)

                    fileOffset += len(line)
            except FormatException as e:
                neededGFA = 1
            except Exception as e:
                logging.error(f'!!!!! Parsing aborted: invalid format at line: {lineNum}')
                neededGFA = -1
                print(e)

        if neededGFA == 1:
            if self.illegalrGFA == 0:
                self.illegalrGFA += 1
                if self.illegalrGFA < 2:
                    logging.warning("!!!!! Not the wanted rGFA format. Please refer to 'Manual' to generate a proper one")
                else:
                    self.illegalrGFA == 0

        #backbone['contigs'] = list(backbone['contigs'].keys())
        for bb in backbone:
             backbone[bb]['contigs'] = list(backbone[bb]['contigs'].keys())
        #nodeIdDDlist = natsorted(nodeIdDDlist)
        if len(rawNodeData) > self.maxNodeCount:
            self.exceedNodeIdCount = True
        nodeIdDDlist = natsorted(list(rawNodeData.keys())[:self.maxNodeCount])

        self.rawNodeData = rawNodeData
        self.rawEdgeData = rawEdgeData

        return {'NodeID':nodeIdDDlist,'samples':samples,'neededGFA':neededGFA,'backbone':backbone}

    def buildGraph(self, targetChr, targetStart=None, targetEnd=None, sampleList=None, nodeIdDict=None):
        G = nx.DiGraph()
        count = 0
        rawNodeData = self.rawNodeData
        rawEdgaData = self.rawEdgeData

        #backbone = {'name':None, 'contigs':{}, 'nodes':{}}
        backbone = {}

        for nodeId in rawNodeData:
            if nodeIdDict and nodeId not in nodeIdDict:
                continue

            field = rawNodeData[nodeId]
            rank, contig, lenBefore, seqLen = field[RAWNODEDATA.rank], field[RAWNODEDATA.contig], field[RAWNODEDATA.lenBefore], field[RAWNODEDATA.seqLen]

            if rank == '0':
                count += 1
                #if contig != targetChr or \
                #   (targetStart and targetStart > lenBefore + seqLen) or \
                #   (targetEnd and targetEnd < lenBefore+1): continue

            G.add_node(nodeId)

        backbone['nodes'] = [nodeId for nodeId in rawNodeData if rawNodeData[nodeId][RAWNODEDATA.rank] == '0']

        for edge in self.rawEdgeData:
            fromNodeId, fromStrand, toNodeId, toStrand = edge.split('|')

            if fromNodeId not in G.nodes or toNodeId not in G.nodes:
                 continue

            if (fromStrand == '-' and toStrand == '+' and fromNodeId in backbone['nodes']) or \
               (toStrand == '-' and fromStrand == '+' and toNodeId in backbone['nodes']):
                fromNodeId, toNodeId = toNodeId, fromNodeId
            elif fromStrand == '-' and toStrand == '-':
                fromNodeId, toNodeId = toNodeId, fromNodeId
                fromStrand, toStrand = '+', '+'

            """
            if not (fromStrand == '-' and toStrand == '-'):
                if fromStrand == '-':
                    fromNodeId = f'{fromNodeId}*'
                if fromNodeId not in G.nodes:
                    G.add_node(fromNodeId)

                if toStrand == '-':
                    toNodeId = f'{toNodeId}*'
                if toNodeId not in G.nodes:
                    G.add_node(toNodeId)
            """

            G.add_edge(fromNodeId, toNodeId, strands=f'{fromStrand}{toStrand}')

        #self.firstNodeId = {}
        #self.firstNodeId[targetChr] = min([id for id in G.nodes if '*' not in id and rawNodeData[id][RAWNODEDATA.rank] == '0'],
        #                                   key=lambda id:rawNodeData[id][RAWNODEDATA.lenBefore])

        #G = G.subgraph(list(nx.node_connected_component(G.to_undirected(), self.firstNodeId[targetChr])))

        self.nodes = {}
        with open(self.gfa, 'rb') as f:
            for nodeId in G.nodes:

                rawNodeId = nodeId if nodeId[-1] != '*' else nodeId[:-1]

                node = rawNodeData[rawNodeId]
                fileOffset, lineLen = node[0], node[1]
                f.seek(node[RAWNODEDATA.fileOffset])
                line = f.read(node[RAWNODEDATA.lineLen])
                row = line.decode('ascii').split()
                #node[RAWNODEDATA.seq] = field[2]

                # to-be-refined
                seq = row[2]
                tags = {}
                for val in row[3:]:
                    lst = val.split(':')
                    tags[lst[0]] = lst[2]

                seqLen = int(tags['LN']) if 'LN' in tags else len(seq)
                seqDesc = seq[0:self.seqDescLen]
                seqLastDesc = seq[-self.seqDescLen:]
                rank = tags['SR'] if 'SR' in tags else ''
                res = tags['SN'] if 'SN' in tags else ''
                if len(res.split(self.SN_delim)) == 2:
                    sample = res.split(self.SN_delim)[0]
                    contig = res.split(self.SN_delim)[1]
                else:
                    sample = tags['SR'] if 'SR' in tags else ''
                    contig = res
                lenBefore = int(tags['SO']) if 'SO' in tags else 0
                inf = {'sv_type':tags['INF'].split(self.SN_delim)[0],'raw':tags['INF']} if 'INF' in tags else {}
                if inf and inf['sv_type'] == 'SV':
                    inf['sv_type'] = tags['INF'].split(self.SN_delim)[1]

                seq = None
                self.nodes[nodeId] = [seqDesc, seqLastDesc, seqLen, seq, sample, contig, lenBefore, rank, inf, None, None]

        self.G = G

    def parseVCF(self, vcf, backbone):
        samples = {}
        neededGFA = False
        backbone = {'name':backbone, 'contigs':{}}
        nodeID = []
        warning = ''
        error = ''

        refChroms = {}

        with open(vcf) as f:
            try:
                for lineNum, line in enumerate(f):
                    line = line.strip()
                    if line[0] == '#':
                        if line.startswith('#CHROM\t'):
                            fields = line.split('\t')
                            samples = fields[9:]
                        else:
                            m = re.search('^##contig=<ID=(.*),length=(.*)>', line.strip())
                            if m and m.group(1) and m.group(2):
                                refChroms[m.group(1)] = {'chrom':m.group(1), 'length':int(m.group(2))}

                        continue

                    fields = line.split('\t')
                    backbone['contigs'][fields[0]] = 1
            except Exception as e:
                logging.error(f'!!!!! Parsing aborted: invalid format at line: {lineNum+1}')
                neededGFA = -1
                #raise

        backbone['contigs'] = list(backbone['contigs'].keys())

        results = {'NodeID':nodeID,'samples':samples,'neededGFA':neededGFA,'backbone':backbone,'warning':warning,'error':error,'refChroms':refChroms}

        return results

    # load nodes and edges
    def loadRGFA(self, targetBb, targetChr=None, targetStart=None, targetEnd=None, sampleList=None, nodeIdList=None):
        G = nx.DiGraph()
        self.nodes = {}
        edges = []
        firstNodeId, firstLenBefore = {}, {}

        # data same as parseRGFA
        samples = {}
        neededGFA = False
        backbone = {'name':None, 'contigs':{}, 'nodes':{}}
        #nodeID = []

        nodeFileInfo = {}

        # find target node
        with open(self.gfa) as f:
            try:
                for lineNum, line in enumerate(f):
                    if line[0] == 'S':
                        row = line.strip().split('\t')

                        nodeId = row[1]

                        if nodeIdList and nodeId not in nodeIdList:
                            continue

                        seq = row[2]

                        # get tag values
                        tags = {}
                        for val in row[3:]:
                            lst = val.split(':')
                            tags[lst[0]] = lst[2]

                        seqLen = int(tags['LN']) if 'LN' in tags else len(seq)
                        #seqDesc = seq[0:self.seqDescLen]
                        #seqLastDesc = seq[-self.seqDescLen:]
                        rank = tags['SR'] if 'SR' in tags else ''
                        res = tags['SN'] if 'SN' in tags else ''

                        """
                        if len(res.split(self.SN_delim)) == 2:
                            sample = res.split(self.SN_delim)[0]
                            contig = res.split(self.SN_delim)[1]
                        else:
                            sample = tags['SR'] if 'SR' in tags else ''
                            contig = res
                        """
                        if len(res.split(self.SN_delim)) >= 2:
                            if 'SR' in tags:
                                sample = res.split(self.SN_delim)[0] if tags['SR'] == '0' else 'Samples'
                            else:
                                sample = res.split(self.SN_delim)[0]
                            contig = res.split(self.SN_delim)[1]
                        else:
                            sample = res.split(self.SN_delim)[0]
                            contig = res

                        lenBefore = int(tags['SO']) if 'SO' in tags else 0


                        # check backbone and sample name
                        if (rank == '0' and sample != targetBb) or \
                           (rank != '0' and sampleList and sample not in sampleList):
                            continue

                        # check chr and pos
                        if (targetChr and contig != targetChr) or \
                           (targetStart and targetStart > lenBefore+seqLen) or \
                           (targetEnd and targetEnd < lenBefore+1):
                            continue

                        # update backbone info
                        if rank == '0':
                            backbone['contigs'][contig] = 1
                            backbone['name'] = sample

                            if contig not in firstNodeId or lenBefore < firstLenBefore[contig]:
                                firstNodeId[contig] = nodeId
                                firstLenBefore[contig] = lenBefore

                            # record backbone node for use in reverse strand
                            backbone['nodes'][nodeId] = 1

                        samples[sample] = 1
                        G.add_node(nodeId)
                    elif line[0] == 'L':
                        row = line.strip().split('\t')
                        fromNodeId, fromStrand, toNodeId, toStrand = row[1:5]
                        tags = {val.split(':')[0]:val.split(':')[2] for val in row[6:]}

                        if fromNodeId not in G.nodes or toNodeId not in G.nodes:
                            continue

                        if (fromStrand == '-' and toStrand == '+' and fromNodeId in backbone['nodes']) or \
                           (toStrand == '-' and fromStrand == '+' and toNodeId in backbone['nodes']):
                            fromNodeId, toNodeId = toNodeId, fromNodeId
                        elif fromStrand == '-' and toStrand == '-':
                            fromNodeId, toNodeId = toNodeId, fromNodeId
                            fromStrand, toStrand = '+', '+'

                        """
                        if not (fromStrand == '-' and toStrand == '-'):
                            if fromStrand == '-':
                                fromNodeId = f'{fromNodeId}*'
                                if fromNodeId not in G.nodes:
                                    G.add_node(fromNodeId)

                            if toStrand == '-':
                                toNodeId = f'{toNodeId}*'
                                if toNodeId not in G.nodes:
                                    G.add_node(toNodeId)
                        """

                        G.add_edge(fromNodeId, toNodeId, strands=f'{fromStrand}{toStrand}', tags=tags)

            except Exception as e:
                logging.error(f'!!!!! Loading aborted: invalid format at line: {lineNum}: {line}')
                neededGFA = -1
                raise
                #print(e)

        H = G.to_undirected()
        subNodes = {}
        for nodeId in nx.node_connected_component(H, firstNodeId[targetChr]):
            subNodes[nodeId] = 1

        self.G = G.subgraph(list(subNodes.keys()))
        self.firstNodeId = firstNodeId

        backbone['contigs'] = list(backbone['contigs'].keys())

        self.inf = {'samples':samples,'neededGFA':neededGFA,'backbone':backbone}
        self.neededGFA = neededGFA
        self.backbone = backbone

        self.colorPalettes()

        if neededGFA == 1:
            if self.illegalrGFA == 0:
                self.illegalrGFA += 1
                if self.illegalrGFA < 2:
                    logging.warning("!!!!! Not the wanted rGFA format. Please refer to 'Manual' to generate a proper one")
                else:
                    self.illegalrGFA == 0

        # load node detail ONLY
        with open(self.gfa) as f:
            try:
                for lineNum, line in enumerate(f):
                    if line[0] == 'S':
                        row = line.strip().split('\t')

                        nodeId = row[1]

                        if nodeId not in subNodes and f'{nodeId}*' not in subNodes: continue

                        seq = row[2]

                        # get tag values
                        tags = {}
                        for val in row[3:]:
                            lst = val.split(':')
                            tags[lst[0]] = lst[2]

                        seqLen = int(tags['LN']) if 'LN' in tags else len(seq)
                        seqDesc = seq[0:self.seqDescLen]
                        seqLastDesc = seq[-self.seqDescLen:]
                        rank = tags['SR'] if 'SR' in tags else ''
                        res = tags['SN'] if 'SN' in tags else ''

                        """
                        if len(res.split(self.SN_delim)) == 2:
                            sample = res.split(self.SN_delim)[0]
                            contig = res.split(self.SN_delim)[1]
                        else:
                            sample = tags['SR'] if 'SR' in tags else ''
                            contig = res
                        """
                        if len(res.split(self.SN_delim)) >= 2:
                            if 'SR' in tags:
                                sample = res.split(self.SN_delim)[0] if tags['SR'] == '0' else 'Samples'
                            else:
                                sample = res.split(self.SN_delim)[0]
                            contig = res.split(self.SN_delim)[1]
                        else:
                            sample = res.split(self.SN_delim)[0]
                            contig = res

                        lenBefore = int(tags['SO']) if 'SO' in tags else 0
                        inf = {'sv_type':tags['INF'].split(self.SN_delim)[0],'raw':tags['INF']} if 'INF' in tags else {}
                        if inf and inf['sv_type'] == 'SV':
                            inf['sv_type'] = tags['INF'].split(self.SN_delim)[1]

                        seq = None
                        if nodeId in subNodes:
                            self.nodes[nodeId] = [seqDesc,seqLastDesc,seqLen,seq,sample,contig,lenBefore,rank,inf,None,None]
                        if f'{nodeId}*' in subNodes:
                            self.nodes[f'{nodeId}*'] = [self.revComp(seqLastDesc),self.revComp(seqDesc),seqLen,seq,sample,contig,lenBefore,rank,inf,None,None]
            except Exception as e:
                logging.error(f'!!!!! Loading aborted: invalid format at line: {lineNum}: {line}')
                neededGFA = -1
                raise
                #print(e)

    def loadLenBeforeDict(self, chromList):
        lenBeforeDict = {}
        with open(self.gfa) as f:
            try:
                for lineNum, line in enumerate(f):
                    if line[0] == 'S':
                        row = line.strip().split('\t')

                        nodeId = row[1]
                        tags = {}
                        for val in row[3:]:
                            lst = val.split(':')
                            tags[lst[0]] = lst[2]

                        rank = tags['SR'] if 'SR' in tags else ''
                        res = tags['SN'] if 'SN' in tags else ''
                        if len(res.split(self.SN_delim)) == 2:
                            sample = res.split(self.SN_delim)[0]
                            contig = res.split(self.SN_delim)[1]
                        else:
                            sample = tags['SR'] if 'SR' in tags else ''
                            contig = res
                        lenBefore = int(tags['SO']) if 'SO' in tags else 0

                        if rank != '0' or contig not in chromList:  continue
                        if contig not in lenBeforeDict:
                            lenBeforeDict[contig] = {}

                        lenBeforeDict[contig][nodeId] = lenBefore

            except Exception as e:
                raise
                #print(e)

        self.lenBeforeDict = lenBeforeDict

    def drawGraph(self, backbone, sampleList, targetChr, targetStart, targetEnd, isGenHtml=True, nodeIdDict=None):
        self.error_unknown = 0
        self.drawgraph_by_nodeids = False

        try:
            logging.info("Loading rGFA ...")
            self.loadRGFA(targetBb=backbone, targetChr=targetChr, targetStart=targetStart, targetEnd=targetEnd, sampleList=sampleList, nodeIdList=nodeIdDict)

            logging.info("Updating nodes ...")
            self.updateNodes()

            logging.info("Generating subgraph ...")
            sampleList = list(self.inf['samples'].keys())
            posDict = {targetChr:{'posFrom':targetStart,'posTo':targetEnd}}
            self.genSubGraph(sampleList, posDict)

            logging.info(f"No. of nodes: {len(self.subGraph.nodes)}, no. of edges: {len(self.subGraph.edges)}")

            if self.emptyGraphSignal == 1:
                return -1

            logging.info("Generating graph results")
            self.genDrawGraphResult(self.subGraph, posDict)
            if isGenHtml:
                outHtml = self.genHtml(self.subGraph, posDict)
                self.drawGraphResult['outHtml'] = outHtml

            logging.info("The graph is generated")

            return self.drawGraphResult
        except Exception as e:
            self.error_unknown = 1

            raise e
            return None

    def checkNodeIds(self, nodeIdDict, targetChr):
        count = {'total':0, 'targetChr':0}

        rawNodeData = self.rawNodeData
        for nodeId in nodeIdDict:
            field = rawNodeData[nodeId]

            contig = field[RAWNODEDATA.contig]
            if contig == targetChr:
                count['targetChr'] += 1

            count['total'] += 1

        return count

    def drawGraphByNodeId(self, sampleList, targetChr, nodeIdDict, isGenHtml=True):
        self.error_unknown = 0
        self.drawgraph_by_nodeids = True

        try:
            logging.info("Building graph ...")
            self.buildGraph(targetChr, nodeIdDict=nodeIdDict)

            #logging.info("updating nodes ...")
            #self.updateNodes()

            """
            logging.info("Generating subgraph ...")
            sampleList = list(self.inf['samples'].keys())
            posDict = {targetChr:{'posFrom':None,'posTo':None}}
            #self.genSubGraph(sampleList, posDict)
            """

            posDict = {}
            self.subGraph = self.G

            logging.info(f"No. of nodes: {len(self.subGraph.nodes)}, no. of edges: {len(self.subGraph.edges)}")

            if self.emptyGraphSignal == 1:
                return -1

            logging.info("Generating graph results")
            self.genDrawGraphResult(self.subGraph, posDict)
            if isGenHtml:
                outHtml = self.genHtml(self.subGraph, posDict)
                self.drawGraphResult['outHtml'] = outHtml

            logging.info("The graph is generated")

            self.subNodesCount = len(self.subGraph.nodes)

            return self.drawGraphResult
        except Exception as e:
            self.error_unknown = 1

            #raise e
            return None

    def drawGraphCmdline(self, backbone, sampleList, targetChr, targetStart, targetEnd, isGenHtml=True, nodeIdFile=None):
        self.error_unknown = 0
        self.drawgraph_by_nodeids = False

        try:
            script = os.path.realpath(__file__)

            cmd = f'"{sys.executable}" "{script}" -b "{backbone}" -g "{self.gfa}" -o "{self.outdir}" -c "{targetChr}" -a drawGraph'
            if targetStart:
                cmd = f'{cmd} -s "{targetStart}"'
            if targetEnd:
                cmd = f'{cmd} -e "{targetEnd}"'
            if sampleList:
                str = ' '.join(['"'+s+'"' for s in sampleList])
                cmd = f'{cmd} -l {str}'
            #print('cmd', cmd)

            stdout = self.runCmdline(cmd)
            DrawGraphResult.fromJson(self, stdout[-1])
        except:
            self.error_unknown = 1
            return None

    def saveNodesInfo(self, nodeIDList, prefix):
        outfile = os.path.join(self.outdir, f'{prefix}.extractedSeqInfo.fasta')
        with open (outfile, 'w') as extract:
            if self.nodesInfo == None:
                self.searchNodes(nodeIDList)
            extract.write('%s\n' % self.nodesInfo)

        logging.info(f"The information of selected node(s) has been saved to {outfile}")

    def getBedGffFormat(self, file):
        filename, fileext = os.path.splitext(file)
        fileext = fileext.upper()
        if not fileext.startswith('.GFF') and not fileext.startswith('.BED') and not fileext.startswith('.GTF'):
            return None

        fields = None
        with open(file) as f:
            for line in f:
                if line[0] == '#': continue
                fields = line.strip().split('\t')
                break

        if (fileext.startswith('.GFF') or fileext.startswith('.GTF')) \
            and fields[3].isdigit() and fields[4].isdigit():
            # check for compulsory field in attr
            # ...

            if fileext.startswith('.GFF'): return 'GFF'
            if fileext.startswith('.GTF'): return 'GTF'
        elif fileext.startswith('.BED') and fields[1].isdigit() and fields[2].isdigit():
            # check for compulsory field in attr
            # ...
            return 'BED'
        else:
            return None

    def loadBedGff(self, file):
        self.bed = {}
        self.warning = ''
        self.error = ''

        format = self.getBedGffFormat(file)
        if not format:
            self.error = f'!!!!! loading Bed/Gff aborted: unknown format'
            logging.error(self.error)
            return

        self.warning = f"The file is in {format} format"
        logging.info(f"Start to load {format.upper()} file: '{file}'")
        with open(file, 'r') as f:
            try:
                for line in f:
                    if line[0] == '#': continue
                    fields = line.strip().split('\t')
                    if format == 'BED':
                        geneId = fields[3]
                        self.bed[geneId] = {'Chr':fields[0],'Start':int(fields[1]),'End':int(fields[2]),
                                            'GeneID':fields[3], 'Orientation':fields[5]}
                    elif format == 'GFF' or format == 'GTF':
                        if fields[2] != "gene": continue

                        if format == 'GFF':
                            attr = dict(x.split('=') for x in fields[8].split(';') if x.count('=') == 1)
                        else:
                            attr = dict(shlex.split(x.strip()) for x in fields[8].split(';') if x.strip().count(' ') == 1)

                        geneId = None
                        for entry in ['ID', 'Name', 'gene_id', 'gene_name']:
                            if entry in attr:
                                geneId = attr[entry]
                                break
                        if not geneId:
                            self.error = f'!!!!! Missing ID/Name/gene_id/gene_name in {format}'
                            logging.error(self.error)
                            break

                        self.bed[geneId] = {'Chr':fields[0],'Start':int(fields[3]),'End':int(fields[4]),
                                            'GeneID':geneId, 'Orientation':fields[6]}
            except:
                self.error = 'Error in parsing the file'
                logging.error(self.error)
                return

        logging.info(f"Complete {format.upper()} file parsing ...")

    def overlapGenes(self, targetBb, geneId):
        if geneId not in self.bed:
            logging.info(f'GeneID {geneId} not found. Action abort')
            return

        geneInfo = self.bed[geneId]

        geneChr, geneStart, geneEnd, geneOri = geneInfo['Chr'], geneInfo['Start'], geneInfo['End'], geneInfo['Orientation']
        self.loadRGFA(targetBb, geneChr, geneStart, geneEnd)
        self.updateNodes()

        posDict = {}
        posDict[geneChr] = {'posFrom':geneStart,'posTo':geneEnd}
        self.genSubGraph(list(self.inf['samples'].keys()), posDict)

        if not self.subGraph.nodes:
            logging.info(f'No overlapped nodes with Gene {geneId}')
            self.noOverlap = True
            return

        newList = natsorted(self.subGraph.nodes.keys())
        startNode = newList[0]
        endNode = newList[-1]

        logging.info(f'Nodes between {startNode} --> {endNode} overlap with Gene: {geneId}')

        self.features=[]
        plotNodes=[]
        if startNode == endNode:
            node = self.subGraph.nodes[startNode]
            self.features.append(GraphicFeature(start=node['posStart'], end=node['posEnd'], strand=0,
                                                color=self.nameCols[node['sample']], label=startNode))

        else:
            connectedNodes = {}
            for nodeId in self.subGraph.nodes:
                node = self.nodes[nodeId]
                if nodeId not in plotNodes:
                    plotNodes.append(nodeId)
                    if len(node[NODE.inf]) != 0 and node[NODE.inf]['sv_type'] == 'DUP':
                        sPos = int(node[NODE.inf]['raw'].split(self.SN_delim)[1])
                        ePos = int(node[NODE.inf]['raw'].split(self.SN_delim)[2])
                        self.features.append(GraphicFeature(start=sPos, end=ePos, strand=0,
                                                color=self.nameCols[node[NODE.sample]], label="This is node '%s'" % nodeId))
                    else:
                        self.features.append(GraphicFeature(start=node[NODE.posStart], end=node[NODE.posEnd], strand=0,
                                                color=self.nameCols[node[NODE.sample]], label="This is node '%s'" % nodeId))

        if geneOri == "+":
            self.features.append(GraphicFeature(start=geneStart, end=geneEnd, strand=+1, color="#cffccc",
                                                label="This is gene '%s @%s: %d-%d'" % (geneId, geneChr, geneStart, geneEnd)))
        if geneOri == "-":
            self.features.append(GraphicFeature(start=geneStart, end=geneEnd, strand=-1, color="#cffccc",
                                                label="This is gene '%s @%s: %d-%d'" % (geneId, geneChr, geneStart, geneEnd)))
        self.record = GraphicRecord(sequence_length=geneEnd*2, features=self.features)

    def drawOverlapGenes(self, geneId):
        #self.record.plot(figure_width=20)
        bokeh_figure = self.record.plot_with_bokeh(figure_width=13, figure_height=8)
        #bokeh.plotting.show(bokeh_figure)

        tool_objs, tool_map = process_tools_arg(bokeh_figure, "xpan,xwheel_zoom,reset,tap")
        bokeh_figure.tools = tool_objs

        drawOverlap = os.path.join(self.outdir, 'drawOverlap_with_Gene-%s.html' % geneId)

        with open(drawOverlap, 'w') as f:
            f.write(file_html(bokeh_figure, CDN, "Nodes overlap with Gene: %s" % geneId))

    def getNodeFromRGFA(self, nodeIdList=None, maxNodeCount=None, getSeq=True):
        nodeIdDict = {nodeId.strip():1 for nodeId in nodeIdList}
        nodes = {}

        with open(self.gfa) as f:
            try:
                for lineNum, line in enumerate(f):
                    row = line.strip().split('\t')
                    if row[0] != 'S': continue

                    nodeId = row[1]
                    if nodeId not in nodeIdDict:
                        continue

                    seq = row[2]
                    tags = {}
                    for val in row[3:]:
                        lst = val.split(':')
                        tags[lst[0]] = lst[2]

                    seqLen = int(tags['LN']) if 'LN' in tags else len(seq)
                    seqDesc = seq[0:self.seqDescLen]
                    seqLastDesc = seq[-self.seqDescLen:]
                    rank = tags['SR'] if 'SR' in tags else ''
                    res = tags['SN'] if 'SN' in tags else ''
                    if len(res.split(self.SN_delim)) >= 2:
                        if 'SR' in tags:
                            sample = res.split(self.SN_delim)[0] if tags['SR'] == '0' else 'Samples'
                        else:
                            sample = res.split(self.SN_delim)[0]
                        contig = res.split(self.SN_delim)[1]
                    else:
                        sample = res.split(self.SN_delim)[0]
                        contig = res

                    lenBefore = int(tags['SO']) if 'SO' in tags else 0
                    inf = {'sv_type':tags['INF'].split(self.SN_delim)[0],'raw':tags['INF']} if 'INF' in tags else {}
                    if inf and inf['sv_type'] == 'SV':
                        inf['sv_type'] = tags['INF'].split(self.SN_delim)[1]

                    if not getSeq:
                        seq = None

                    nodes[nodeId] = [seqDesc,seqLastDesc,seqLen,seq,sample,contig,lenBefore,rank,inf,None,None]

                    del nodeIdDict[nodeId]
                    if not nodeIdDict:
                        break
            except Exception as e:
                raise

                logging.error(f'!!!!! Parsing aborted: invalid format at line: {lineNum}')
                #print(e)

        return nodes

    def searchNodes(self, nodeIdList):
        nodesInfo = []

        nodes = self.getNodeFromRGFA(nodeIdList)

        for nodeId in nodeIdList:
            if nodeId in nodes:
                node = nodes[nodeId]
                info = self.formatNodeOutput(nodeId, node, showSeqDesc=False)
                nodesInfo.append(f">{nodeId}\t{info['title']}\n{node[NODE.seq]}")

                if len(node[NODE.seq]) > self.bigNodeSeqLen:
                    self.bigNodeList.append(nodeId)

        self.nodesInfo = '\n'.join(nodesInfo)

    def revComp(self, seq):
        trans = str.maketrans('ACGTN*', 'TGCAN*')
        return seq.translate(trans)[::-1]

    def updateNodes(self):
        G = self.G

        for contig in self.firstNodeId:
            nodeId = self.firstNodeId[contig]

            # init firstNodeId
            node = self.nodes[nodeId]
            if node[NODE.sample] == self.backbone['name'] and node[NODE.chr] == contig:
                node[NODE.posStart] = node[NODE.lenBefore] + 1
                node[NODE.posEnd] = node[NODE.posStart] + node[NODE.len] - 1

            nextPos = 0
            lastNode = None

            for tuple in list(nx.bfs_edges(G, source=nodeId)):
                fromNode = self.nodes[tuple[0]]
                toNode = self.nodes[tuple[1]]

                if fromNode[NODE.sample] == self.backbone['name'] and fromNode[NODE.chr] == contig:
                    fromNode[NODE.posStart] = fromNode[NODE.lenBefore] + 1
                    fromNode[NODE.posEnd] = fromNode[NODE.posStart] + fromNode[NODE.len] - 1
                    nextPos = fromNode[NODE.posEnd] + 1

                if toNode[NODE.sample] == self.backbone['name'] and toNode[NODE.chr] == contig:
                    toNode[NODE.posStart] = toNode[NODE.lenBefore] + 1
                    toNode[NODE.posEnd] = toNode[NODE.posStart] + toNode[NODE.len] - 1
                    nextPos = toNode[NODE.posEnd] + 1
                else:
                    toNode[NODE.posStart] = nextPos
                    toNode[NODE.posEnd] = toNode[NODE.posStart]

                lastNode = toNode

            if lastNode and lastNode[NODE.sample] == self.backbone['name']:
                lastNode[NODE.posStart] = lastNode[NODE.lenBefore] + 1

        self.G = G

    # e.g. posDict = {'Chr01':{'posFrom':1,'posTo':2}}
    def genSubGraph(self, sampleList, posDict):
        if not self.G:
            logging.error(f'genSubGraph(): need to run genGraph() first')
            return nx.DiGraph()

        G = self.G
        H = self.G.to_undirected()

        subNodes = []
        notConnectCount = 0
        connectCount = 0

        nodeCount = 0
        for contig in posDict:
            if contig != 'all':
                posFrom = posDict[contig]['posFrom']
                posTo = posDict[contig]['posTo']
                anyNodeId = self.firstNodeId[contig]

                nodeIdDict = {}
                for nodeId in nx.node_connected_component(H, anyNodeId):
                    nodeIdDict[nodeId] = self.nodes[nodeId][NODE.lenBefore]
                sortedDict = {k: v for k, v in sorted(nodeIdDict.items(), key=lambda item: item[1])}

                for nodeId in sortedDict:
                    node = self.nodes[nodeId]

                    # note: links from out-of-region nodes
                    if node[NODE.posStart] is None:
                        node[NODE.posStart] = 0
                        node[NODE.posEnd] = 0

                    if sampleList and (node[NODE.rank] != '0' and node[NODE.sample] != self.backbone and node[NODE.sample] not in sampleList):
                        continue

                    if (not posTo or node[NODE.posStart] <= posTo) and (not posFrom or node[NODE.posEnd] >= posFrom):
                        subNodes.append(nodeId)
            else:
                for contig in self.firstNodeId:
                    anyNodeId = self.firstNodeId[contig]
                    for nodeId in nx.node_connected_component(H, anyNodeId):
                        subNodes.append(nodeId)

        subGraph = G.subgraph(subNodes)
        self.subNodes = subGraph.nodes
        if len(subGraph.nodes) > self.maxNodesDisplay:
            self.overNodeLimitSignal = 1
        if len(subGraph.nodes) == 0:
            self.emptyGraphSignal = 1
        logging.info(f'subGraph: number of nodes: {len(subGraph.nodes)}, number of edges: {len(subGraph.edges)}')

        self.nodes = {nodeId:self.nodes[nodeId] for nodeId in subGraph.nodes}
        self.subGraph = subGraph
        self.G = None

    def formatNodeOutput(self, nodeId, node, showSeqDesc=True):
        shape = getVar(copied, 'nodes',f"{node[NODE.inf]['sv_type']}_shape") if 'sv_type' in node[NODE.inf] else getVar(copied, 'nodes', 'BB_shape')
        shape_cy = getVar(copied, 'cytoscape',f"{node[NODE.inf]['sv_type']}_shape") if 'sv_type' in node[NODE.inf] else getVar(copied, 'cytoscape', 'BB_shape')
        sv_type = node[NODE.inf]['sv_type'] if 'sv_type' in node[NODE.inf] else ''

        size = 1 if not node[NODE.len] else float(f"{math.log(abs(node[NODE.len]), 10)*8 + 1:.1f}")
        color = self.nameCols[node[NODE.sample]] if node[NODE.sample] in self.nameCols else '#A2A2A2'
        sample = node[NODE.sample]
        pos = node[NODE.lenBefore] if node[NODE.lenBefore] else 0
        inf = node[NODE.inf]['raw'] if 'raw' in node[NODE.inf] else ''
        seqDesc = f"{node[NODE.seqDesc]}..." if node[NODE.len] > len(node[NODE.seqDesc]) and node[NODE.seqDesc] != '*' else node[NODE.seqDesc]

        title = f"NodeId: {nodeId}; Resource: {node[NODE.sample]}_{node[NODE.chr]}; Len: {node[NODE.len]}"
        if sample == self.backbone['name']: title += f"; Pos: {pos} - {pos + node[NODE.len] - 1}"
        if inf: title += f"; Info: {inf}"
        if showSeqDesc: title += f"; Seq: {seqDesc}"

        return {'color':color,'id':nodeId,'label':nodeId,'shape':shape,'size':size,'title':title,'shape_cy':shape_cy,'sample':sample,'sv_type':sv_type}

    def formatEdgeOutput(self, edge):
        sourceLabel = '*' if 'strands' in edge[2] and edge[2]['strands'][0] == '-' else ''
        targetLabel = '*' if 'strands' in edge[2] and edge[2]['strands'][1] == '-' else ''

        label = f'({sourceLabel},{targetLabel})' if sourceLabel or targetLabel else ''

        return {'from':edge[0],'to':edge[1],'arrows':'to','sourceLabel':sourceLabel,'targetLabel':targetLabel,'label':label,'color':'red','labelHighlightBold':'true'}

    def genDrawGraphResult(self, graph, posDict):
        self.colorPalettes()

        nodes = [self.formatNodeOutput(nodeId, self.nodes[nodeId]) for nodeId in graph.nodes]
        edges = [self.formatEdgeOutput(edge) for edge in graph.edges(data=True)]

        inNodeIdList = [n for n,d in graph.in_degree() if d == 0]
        outNodeIdList = [n for n,d in graph.out_degree() if d == 0]

        startNodeId, endNodeId = None, None
        try:
            startNodeId = min([id for id in graph.nodes if self.nodes[id][NODE.sample] == self.backbone['name'] and self.nodes[id][NODE.posStart]], key=lambda id:self.nodes[id][NODE.posStart])
            endNodeId = max([id for id in graph.nodes if self.nodes[id][NODE.sample] == self.backbone['name'] and self.nodes[id][NODE.posEnd]], key=lambda id:self.nodes[id][NODE.posEnd])
        except:
            pass

        if startNodeId:
            nodes.append({'color':'green','id':f'start','label':'start','shape':'star','size':20,'title':'start','shape_cy':'star'})
            edges.append({'from':f'start','to':startNodeId,'arrows':'to'})

        if endNodeId:
            nodes.append({'color':'red','id':f'end','label':'end','shape':'star','size':20,'title':'end','shape_cy':'star'})
            edges.append({'from':endNodeId,'to':f'end','arrows':'to'})

        for idx, nodeId in enumerate(inNodeIdList):
            if nodeId == startNodeId: continue
            labelNodeId = f'other_start_{idx+1}'
            nodes.append({'color':'gray','id':labelNodeId,'label':labelNodeId,'shape':'star','size':20,'title':labelNodeId,'shape_cy':'star'})
            edges.append({'from':labelNodeId,'to':nodeId,'arrows':'to'})

        for idx, nodeId in enumerate(outNodeIdList):
            if nodeId == endNodeId: continue
            labelNodeId = f'other_end_{idx+1}'
            nodes.append({'color':'gray','id':labelNodeId,'label':labelNodeId,'shape':'star','size':20,'title':labelNodeId,'shape_cy':'star'})
            edges.append({'from':nodeId,'to':labelNodeId,'arrows':'to'})

        self.drawGraphResult = {'error':False, 'nodes_data':nodes,'edges_data':edges}

    def genCyDataFromDrawGraphResult(self, drawGraphResult):
        cyData = []
        outputNodeInfo = {}

        nodes = drawGraphResult['nodes_data']
        edges = self.drawGraphResult['edges_data']

        for node in nodes:
            size = node['size']
            ele = {'data':{'id':node['id'],'name':node['id'],'weight':1,'size':size,'color':node['color'],'shape':node['shape_cy'],'title':node['title']}}
            cyData.append(ele)

            outputNodeInfo[node['id']] = ele

        for edge in edges:
            sourceLabel = edge['sourceLabel'] if 'sourceLabel' in edge else ''
            targetLabel = edge['targetLabel'] if 'targetLabel' in edge else ''
            color = outputNodeInfo[edge['from']]['data']['color']
            # avoid too light color for edges
            color = f'#{int(color[1:3],16)//2:02x}{int(color[3:5],16)//2:02x}{int(color[5:],16)//2:02x}' if color[0] == '#' else color

            if sourceLabel or targetLabel:
                color = 'red'

            arrow = 'triangle-backcurve'
            ele = {'data':{'source':edge['from'],'target':edge['to'],'weight':1,'color':color,'arrow':arrow,'sourceLabel':sourceLabel,'targetLabel':targetLabel}}
            cyData.append(ele)

        return cyData

    def genHtml(self, graph, posDict={}, outHtmlPrefix=None):
        genHtmlInfo = {'vis':{'template':os.path.join(os.path.dirname(os.path.realpath(__file__)),'template','htmlTemplate.html')}, \
                       'cytoscape':{'template':os.path.join(os.path.dirname(os.path.realpath(__file__)),'template','htmlTemplate_cytoscape.html')}}

        outFile = {}
        interaction = getVar(copied, 'enable', 'interaction')
        graphLayoutModification = getVar(copied, 'enable', 'graphLayoutModification')
        addRemoveNodes = getVar(copied, 'enable', 'addRemoveNodes')
        canvas_height = getVar(copied, 'canvas', 'height')
        canvas_width = getVar(copied, 'canvas', 'width')
        #maxNodesDisplay = int(getVar(copied, "nodes", "maxNodesDisplay"))

        nodes = self.drawGraphResult['nodes_data']
        edges = self.drawGraphResult['edges_data']

        colors = {node['sample']:node['color'] for node in nodes if 'sample' in node}
        shapes = {node['sv_type']:node['shape'] for node in nodes if 'sv_type' in node and node['sv_type']}
        shapes_cy = {node['sv_type']:node['shape_cy'] for node in nodes if 'sv_type' in node and node['sv_type']}
        hasReversed = True if [node for node in nodes if node['id'][-1] == '*'] else False

        if self.drawgraph_by_nodeids:
            title = f"PanGraph -- by node IDs"
        else:
            tmp = []
            for contig in posDict:
                if posDict[contig]['posFrom'] == None and posDict[contig]['posTo'] == None:
                    tmp = [f"{contig}"]
                if posDict[contig]['posFrom'] != None and posDict[contig]['posTo'] != None:
                    tmp = [f"{contig}: {posDict[contig]['posFrom']} - {posDict[contig]['posTo']}"]
                if posDict[contig]['posFrom'] == None and posDict[contig]['posTo'] != None:
                    tmp = [f"{contig}: 1 - {posDict[contig]['posTo']}"]
                if posDict[contig]['posFrom'] != None and posDict[contig]['posTo'] == None:
                    tmp = [f"{contig}: {posDict[contig]['posFrom']} -"]
            title = f"PanGraph -- {', '.join(tmp)}"


        # vis.js specific
        filterInfo = 'false'
        if 'Yes' in [interaction, graphLayoutModification, addRemoveNodes]:
            modification = []
            if interaction == 'Yes':
                modification.append('interaction')
            if graphLayoutModification == 'Yes':
                modification.append('physics')
            if addRemoveNodes == 'Yes':
                #modification.append('nodes')
                #modification.append('edges')
                modification.append('manipulation')
            filterInfo = f'true, "filter": {modification}'

        for libraryName in genHtmlInfo:
            htmlInfo = genHtmlInfo[libraryName]

            templateHtml = os.path.join(os.path.dirname(os.path.realpath(__file__)), htmlInfo['template'])
            if not outHtmlPrefix:
                now = datetime.datetime.now()
                self.name = now.strftime("%Y-%m-%d-%H_%M_%S")
                outHtml = os.path.join(self.outdir, f'{self.name}_{libraryName}.html')
            else:
                outHtml = f"{outHtmlPrefix}_{libraryName}.html"

            outFile[libraryName] = outHtml

            #logging.info(f'The template HTML file is: {templateHtml}')

            zipFile = 'static.zip'
            zipFile = os.path.join(os.path.dirname(os.path.realpath(__file__)), 'template', zipFile)
            with zipfile.ZipFile(zipFile, 'r') as zip:
                zip.extractall(self.outdir)

            if libraryName == 'vis' and len(graph.nodes) <= self.maxNodesDisplay:
                with open(templateHtml) as f_in, open(outHtml,'w') as f_out:
                    data = f_in.read()
                    data = data.replace('{%title%}', title)
                    data = data.replace('{%INFO%}', filterInfo)
                    data = data.replace('{%canvas_height%}', canvas_height)
                    data = data.replace('{%canvas_width%}', canvas_width)
                    data = data.replace('{%nodes%}', json.dumps(nodes))
                    data = data.replace('{%edges%}', json.dumps(edges))
                    data = data.replace('{%colors%}', json.dumps(colors))
                    data = data.replace('{%shapes%}', json.dumps(shapes))
                    data = data.replace('{%has_reversed%}', json.dumps(hasReversed))
                    f_out.write(data)
                    logging.info(f'The output HTML file is: {outHtml}')

            if libraryName == 'cytoscape' and len(graph.nodes) > self.maxNodesDisplay:
                with open(templateHtml) as f_in, open(outHtml,'w') as f_out:
                    data = f_in.read()
                    cyData = self.genCyDataFromDrawGraphResult(self.drawGraphResult)

                    data = data.replace('{%title%}', title)
                    data = data.replace('{%data%}', json.dumps(cyData, indent=4, sort_keys=True))
                    data = data.replace('{%colors%}', json.dumps(colors))
                    data = data.replace('{%shapes%}', json.dumps(shapes_cy))
                    data = data.replace('{%has_reversed%}', json.dumps(hasReversed))
                    f_out.write(data)
                    logging.info(f'The output HTML file is: {outHtml}')

        #return outFile
        if len(graph.nodes) <= self.maxNodesDisplay:
            return outFile['vis']
        else:
            return outFile['cytoscape']

    def test(self):
        pass

    def binSearch(self, a, x):
        i = bisect_left(a, x)
        if i :
            return i-1
        else:
            return -1

    def getGeneNodeOverlapThreadhold(self):
        num = int(getVar(copied, "nodes", "geneNodeOverlapCntThreshold"))
        return num

    # note: bed info have been loaded into self.bed
    def nodeGeneOverlap(self, overlapNodeCountThreshold):
        logging.info(f'finding genes which overlap more than {overlapNodeCountThreshold} nodes ...')

        results = {}

        bedInfo = self.bed

        # get sorted list of backbone node position
        posStartDict = {}
        for contig in self.lenBeforeDict:
            posStartDict[contig] = sorted(list(self.lenBeforeDict[contig].values()))

        for geneId in bedInfo:
            gene = bedInfo[geneId]
            contig = gene['Chr']

            if contig not in posStartDict:
                continue

            startIdx = self.binSearch(posStartDict[contig], gene['Start'])
            endIdx = self.binSearch(posStartDict[contig], gene['End'])

            if endIdx - startIdx + 1 > overlapNodeCountThreshold:
                results[geneId] = {'gene_chr':contig,'gene_start':posStartDict[contig][startIdx],'gene_end':posStartDict[contig][endIdx]}

        logging.info(f'finding genes done')

        return results

if __name__=="__main__":
    parser = ArgumentParser(description='Generate graph html')
    parser.add_argument('-g', dest='gfa', help='the gfa file', type = str)
    parser.add_argument('-o', dest='outdir', help='the output directory', type=str)
    parser.add_argument('-c', dest='chr', help='the name of the chromosome(s)', type=str)
    parser.add_argument('-s', dest='start', help='start pos', type=int)
    parser.add_argument('-e', dest='end', help='end pos', type=int)
    parser.add_argument('-l', dest='samplelist', nargs='*', help='sample list', type=str)
    parser.add_argument('-n', dest='nodeidlist', nargs='*', help='nodeID list', type=str)
    parser.add_argument('-b', dest='backbone', help='backbone', type=str)

    parser.add_argument('-a', dest='action', help='action [parseRGFA, drawGraph]', type = str)

    args = parser.parse_args()

    if None not in [args.gfa, args.outdir]:
        if args.action == 'drawGraph':
            panGraph = PanGraph(args.gfa, args.outdir, parseRGFA=False)
            #panGraph = PanGraph(args.gfa, args.outdir)
            drawGraphResult = panGraph.drawGraph(args.backbone, args.samplelist, args.chr, args.start, args.end)
            print(DrawGraphResult.toJson(panGraph, drawGraphResult))
        elif args.action == 'drawGraphCmdline': # for debug only
            panGraph = PanGraph(args.gfa, args.outdir, parseRGFA=False)
            panGraph.drawGraphCmdline(args.backbone, args.samplelist, args.chr, args.start, args.end)
            print(panGraph.drawGraphResult['outHtml'])
            print(panGraph.subNodesCount)
        elif args.action == 'test':
            nodeIdDict = {nodeId:1 for nodeId in args.nodeidlist}
            panGraph = PanGraph(args.gfa, args.outdir, nodeIdDict=nodeIdDict)
            panGraph.drawGraphByNodeId(sampleList=args.samplelist, targetChr=args.chr, nodeIdDict=nodeIdDict)
        elif args.action == 'parseRGFA':
            panGraph = PanGraph(args.gfa, args.outdir, parseRGFA=True)
        else:
            panGraph = PanGraph(args.gfa, args.outdir, parseRGFA=False)
            panGraph.drawGraph(args.backbone, args.samplelist, args.chr, args.start, args.end)
    else:
        print('\n%s\n' % parser.print_help())
