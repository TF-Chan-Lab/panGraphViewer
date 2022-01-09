from argparse import ArgumentParser
import networkx as nx
import sys
import math
import logging
from subprocess import Popen, PIPE
import shlex

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
import re

from pangraphviewer.utilities import *

class PanGraph:
    seqDescLen = 10
    SN_delim = getVar(copied, 'nodes', 'SN_delim', mustHave=True)
    maxNodesLimit = int(getVar(copied, 'nodes', 'maxNodesLimit', mustHave=True))

    def __init__(self, gfa, outdir, parseRGFA=True):
        self.gfa = gfa
        self.outdir = outdir

        prev = os.umask(0)
        os.makedirs(self.outdir, mode=0o777, exist_ok=True)
        os.umask(prev)

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

        # init
        if parseRGFA:
            logging.info("parsing rGFA ...")
            self.inf = self.parseRGFA()
            self.neededGFA = self.inf['neededGFA']
            self.backbone = self.inf['backbone']
            self.colorPalettes()

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
    def parseRGFA(self):
        samples = {}
        neededGFA = False
        backbone = {'name':None, 'contigs':{}}
        nodeID = []
        warning = ''
        error = ''

        with open(self.gfa) as f:
            try:
                for lineNum, line in enumerate(f):
                    if line[0] != 'S': continue

                    row = line.strip().split('\t')

                    nodeId = row[1]
                    nodeID.append(nodeId)

                    tags = {}
                    for val in row[3:]:
                        lst = val.split(':')
                        tags[lst[0]] = lst[2]

                    rank = tags['SR'] if 'SR' in tags else ''
                    res = tags['SN'] if 'SN' in tags else ''
                    if len(res.split(self.SN_delim)) >= 2:
                        sample = res.split(self.SN_delim)[0]
                        contig = self.SN_delim.join(res.split(self.SN_delim)[1:])
                    else:
                        neededGFA = True
                        sample = rank
                        contig = res

                    samples[sample] = 1
                    if rank == '0':
                        backbone['contigs'][contig] = 1
                        backbone['name'] = sample
            except Exception as e:
                logging.error(f'!!!!! Parsing aborted: invalid format at line: {lineNum+1}')
                error = f'!!!!! Parsing aborted: invalid format at line: {lineNum+1}'
                neededGFA = -1
                #print(e)

        if neededGFA == 1:
            if self.illegalrGFA == 0:
                self.illegalrGFA += 1
                if self.illegalrGFA < 2:
                    logging.warning("!!!!! Not the wanted rGFA format. Please refer to 'Manual' to generate a proper one")
                else:
                    self.illegalrGFA == 0

        backbone['contigs'] = list(backbone['contigs'].keys())

        return {'NodeID':nodeID,'samples':samples,'neededGFA':neededGFA,'backbone':backbone, 'warning':warning, 'error':error}

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
                    else:
                        fields = line.split('\t')
                        if len(fields) < 8:
                            raise ValueError('Invalid format')
                        backbone['contigs'][fields[0]] = 1
            except Exception as e:
                message = f'Parsing aborted: invalid format at line: {lineNum+1}'
                logging.error(message)
                error = message

        backbone['contigs'] = list(backbone['contigs'].keys())

        results = {'NodeID':nodeID,'samples':samples,'neededGFA':neededGFA,'backbone':backbone,'warning':warning,'error':error,'refChroms':refChroms}

        return results

    def parseBedGff(self, file):
        output = {'gene':{},'warning':'','error':'','bed_contigs':{}}

        logging.info("Start to parse bed/gff/gtf file: '%s'" % file)

        format = self.getBedGffFormat(file)
        if not format:
            output['error'] = 'Unknown format'
            return output

        warning = f'The file is in {format} format'
        gene, bedContigs, warning, error = {},{},'',''
        with open(file, 'r') as f:
            try:
                for line in f:
                    if line[0] == '#': continue
                    fields = line.strip().split('\t')
                    bedContigs[fields[0]] = 1

                    if format == 'BED':
                        geneId = fields[3]
                        gene[geneId] = 1
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
                            error = f'Missing ID/Name/gene_id/gene_name in {format}'
                            logging.error(error)
                            break

                        gene[geneId] = 1
            except:
                #raise
                output['error'] = 'Error in parsing the file'
                return output

        return {'gene':gene,'warning':warning,'error':error,'bed_contigs':list(bedContigs.keys())}

    # load nodes and edges
    def loadRGFA(self, targetChr = None):
        nodes = {}
        edges = []
        firstNodeId, firstLenBefore = {}, {}

        # data same as parseRGFA
        samples = {}
        neededGFA = False
        backbone = {'name':None, 'contigs':{}}
        nodeID = []

        with open(self.gfa) as f:
            try:
                for lineNum, line in enumerate(f):
                    if line[0] == 'S':
                        #if targetChr and targetChr not in line: continue

                        row = line.strip().split('\t')

                        nodeId = row[1]
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
                        if len(res.split(self.SN_delim)) == 2:
                            sample = res.split(self.SN_delim)[0]
                            contig = res.split(self.SN_delim)[1]
                        else:
                            sample = tags['SR'] if 'SR' in tags else ''
                            contig = res
                        lenBefore = int(tags['SO']) if 'SO' in tags else 0

                        # update backbone info
                        if rank == '0':
                            if targetChr and contig != targetChr: continue

                            backbone['contigs'][contig] = 1
                            backbone['name'] = sample

                            if contig not in firstNodeId or lenBefore < firstLenBefore[contig]:
                                firstNodeId[contig] = nodeId
                                firstLenBefore[contig] = lenBefore

                        # update nodes
                        samples[sample] = 1
                        inf = {'sv_type':tags['INF'].split('_')[0],'raw':tags['INF']} if 'INF' in tags else {}
                        nodes[nodeId] = {'nodeId':nodeId,'seqDesc':seqDesc,'seqLastDesc':seqLastDesc,'len':seqLen,
                                         'sample':sample,'chr':contig,'lenBefore':lenBefore,'rank':rank,'inf':inf}
                    elif line[0] == 'L':
                        row = line.strip().split('\t')
                        type, fromNodeId, fromStrand, toNodeId, toStrand = row[:5]

                        if fromNodeId not in nodes or toNodeId not in nodes: continue

                        if (nodes[fromNodeId]['rank'] == '0' and fromStrand == '-' and nodes[toNodeId]['rank'] != '0' and toStrand == '+' ) or \
                           (nodes[fromNodeId]['rank'] != '0' and fromStrand == '+' and nodes[toNodeId]['rank'] == '0' and toStrand == '-' ):
                           fromNodeId, toNodeId = toNodeId, fromNodeId

                        edges.append({'type':type,'fromNodeId':fromNodeId,'fromStrand':fromStrand,
                                           'toNodeId':toNodeId,'toStrand':toStrand})
            except Exception as e:
                logging.error(f'!!!!! Loading aborted: invalid format at line: {lineNum+1}')
                neededGFA = -1
                #print(e)

        if neededGFA == 1:
            if self.illegalrGFA == 0:
                self.illegalrGFA += 1
                if self.illegalrGFA < 2:
                    logging.warning("!!!!! Not the wanted rGFA format. Please refer to 'Manual' to generate a proper one")
                else:
                    self.illegalrGFA == 0

        self.nodes = nodes
        self.edges = edges
        self.firstNodeId = firstNodeId

        backbone['contigs'] = list(backbone['contigs'].keys())

        self.inf = {'samples':samples,'neededGFA':neededGFA,'backbone':backbone}
        self.neededGFA = neededGFA
        self.backbone = backbone

        self.colorPalettes()

    def drawGraph(self, sampleList, targetChr, targetStart, targetEnd, isGenHtml=True):
        logging.info("loading rGFA ...")
        self.loadRGFA(targetChr=targetChr)

        logging.info("generating graph ...")
        self.genGraph()

        logging.info("updating nodes ...")
        self.updateNodes()

        logging.info("generating subgraph ...")
        posDict = {targetChr:{'posFrom':targetStart,'posTo':targetEnd}}
        subGraph = self.genSubGraph(sampleList, posDict)

        if self.emptyGraphSignal == 1:
            return -1

        logging.info("generating graph results")
        self.genDrawGraphResult(subGraph, posDict)
        if isGenHtml:
            outHtml = self.genHtml(subGraph, posDict)
            self.drawGraphResult['outHtml'] = outHtml

        logging.info("The graph is generated")

        return self.drawGraphResult

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

    def overlapGenes(self, geneId):
        if geneId not in self.bed:
            logging.error(f'Gene {geneId} not found in bed. Gene ignored')
            return

        geneInfo = self.bed[geneId]

        self.loadRGFA()
        self.genGraph()
        self.updateNodes()

        if geneId not in self.bed:
            logging.info(f'GeneID {geneId} not found. Action abort')
            return

        geneChr, geneStart, geneEnd, geneOri = geneInfo['Chr'], geneInfo['Start'], geneInfo['End'], geneInfo['Orientation']

        posDict = {}
        posDict[geneChr] = {'posFrom':geneStart,'posTo':geneEnd}
        subGraph = self.genSubGraph(list(self.inf['samples'].keys()), posDict)

        if not subGraph.nodes:
            logging.info(f'No overlapped nodes with Gene {geneId}')
            self.noOverlap = True
            return

        newList = natsorted(subGraph.nodes.keys())
        startNode = newList[0]
        endNode = newList[-1]

        logging.info(f'Nodes between {startNode} --> {endNode} overlap with Gene: {geneId}')

        self.features=[]
        plotNodes=[]
        if startNode == endNode:
            node = subGraph.nodes[startNode]
            self.features.append(GraphicFeature(start=node['posStart'], end=node['posEnd'], strand=0,
                                                color=self.nameCols[node['sample']], label=startNode))

        else:
            connectedNodes = {}
            for nodeId in subGraph.nodes:
                node = subGraph.nodes[nodeId]
                if nodeId not in plotNodes:
                    plotNodes.append(nodeId)
                    if len(node['inf']) != 0 and node['inf']['sv_type'] == 'DUP':
                        sPos = int(node['inf']['raw'].split('_')[1])
                        ePos = int(node['inf']['raw'].split('_')[2])
                        self.features.append(GraphicFeature(start=sPos, end=ePos, strand=0, 
                                                color=self.nameCols[node['sample']], label="This is node '%s'" % nodeId))
                    else:
                        self.features.append(GraphicFeature(start=node['posStart'], end=node['posEnd'], strand=0,
                                                color=self.nameCols[node['sample']], label="This is node '%s'" % nodeId))

        if geneOri == "+":
            self.features.append(GraphicFeature(start=geneStart, end=geneEnd, strand=+1, color="#cffccc",
                                                label="This is gene '%s @%s: %d-%d'" % (geneId, geneChr, geneStart, geneEnd)))
        if geneOri == "-":
            self.features.append(GraphicFeature(start=geneStart, end=geneEnd, strand=-1, color="#cffccc",
                                                label="This is gene '%s @%s: %d-%d'" % (geneId, geneChr, geneStart, geneEnd)))
        self.record = GraphicRecord(sequence_length=geneEnd*2, features=self.features)

    def drawOverlapGenes(self, geneId):
        self.record.plot(figure_width=20)
        bokeh_figure = self.record.plot_with_bokeh(figure_width=13, figure_height=8)
        #bokeh.plotting.show(bokeh_figure)

        tool_objs, tool_map = process_tools_arg(bokeh_figure, "xpan,xwheel_zoom,reset,tap")
        bokeh_figure.tools = tool_objs

        drawOverlap = os.path.join(self.outdir, 'drawOverlap_with_Gene-%s.html' % geneId)

        with open(drawOverlap, 'w') as f:
            f.write(file_html(bokeh_figure, CDN, "Nodes overlap with Gene: %s" % geneId))

    def getNodeFromRGFA(self, nodeIdList):
        nodeIdDict = {nodeId:1 for nodeId in nodeIdList}
        nodes = {}

        with open(self.gfa) as f:
            try:
                for lineNum, line in enumerate(f):
                    row = line.strip().split('\t')
                    if row[0] != 'S': continue

                    nodeId = row[1]
                    if nodeId not in nodeIdDict and f'{nodeId}*' not in nodeIdDict:
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
                        sample = res.split(self.SN_delim)[0]
                        contig = self.SN_delim.join(res.split(self.SN_delim)[1:])
                    else:
                        sample = rank
                        contig = res
                    lenBefore = int(tags['SO']) if 'SO' in tags else 0
                    inf = {'sv_type':tags['INF'].split('_')[0],'raw':tags['INF']} if 'INF' in tags else {}

                    nodes[nodeId] = {'nodeId':nodeId,'seqDesc':seqDesc,'seqLastDesc':seqLastDesc,'len':seqLen,'seq':seq,
                                     'sample':sample,'chr':contig,'lenBefore':lenBefore,'rank':rank,'inf':inf}
            except Exception as e:
                raise

                logging.error(f'!!!!! Parsing aborted: invalid format at line: {lineNum+1}')
                #print(e)

        return nodes

    def searchNodes(self, nodeIdList):
        self.bigNodeSeqLen = int(getVar(copied, 'nodes', 'maxNodeLenDisplay'))
        self.bigNodeList = []
        nodesInfo = []

        nodes = self.getNodeFromRGFA(nodeIdList)

        for nodeId in nodeIdList:
            if nodeId in nodes:
                node = nodes[nodeId]
                info = self.formatNodeOutput(nodeId, node, showSeq=False)
                nodesInfo.append(f">{nodeId}\t{info['title'].rsplit(';',1)[0]}\n{node['seq']}")

                if len(node['seq']) > self.bigNodeSeqLen:
                    self.bigNodeList.append(nodeId)

        self.nodesInfo = '\n'.join(nodesInfo)

    def revComp(self, seq):
        trans = str.maketrans('ACGTN*', 'TGCAN*')
        return seq.translate(trans)[::-1]

    def genGraph(self):
        G = nx.DiGraph()
        for nodeId in self.nodes:
            G.add_node(nodeId, **self.nodes[nodeId])

        for edge in self.edges:
            fromNodeId = edge['fromNodeId']
            fromStrand = edge['fromStrand']
            toNodeId = edge['toNodeId']
            toStrand = edge['toStrand']

            if fromNodeId not in G.nodes or toNodeId not in G.nodes:
                continue

            if fromStrand == '-':
                newNodeId = f'{fromNodeId}*'
                if newNodeId not in G.nodes:
                    G.add_node(newNodeId, **self.nodes[fromNodeId])
                    G.nodes[newNodeId]['nodeId'] = newNodeId
                    G.nodes[newNodeId]['seqDesc'] = self.revComp(self.nodes[fromNodeId]['seqLastDesc'])
                    G.nodes[newNodeId]['seqLastDesc'] = self.revComp(self.nodes[fromNodeId]['seqDesc'])
                fromNodeId = newNodeId

            if toStrand == '-':
                newNodeId = f'{toNodeId}*'
                if newNodeId not in G.nodes:
                    G.add_node(newNodeId, **self.nodes[toNodeId])
                    G.nodes[newNodeId]['nodeId'] = newNodeId
                    G.nodes[newNodeId]['seqDesc'] = self.revComp(self.nodes[toNodeId]['seqLastDesc'])
                    G.nodes[newNodeId]['seqLastDesc'] = self.revComp(self.nodes[toNodeId]['seqDesc'])
                toNodeId = newNodeId

            G.add_edge(fromNodeId, toNodeId)

        # fixing for added reverse
        for edge in self.edges:
            fromNodeId, fromStrand = edge['fromNodeId'], edge['fromStrand']
            toNodeId, toStrand = edge['toNodeId'], edge['toStrand']

            if fromStrand == '+' and toStrand == '+' and \
               f'{fromNodeId}*' in G.nodes and f'{toNodeId}*' in G.nodes:
               G.add_edge(f'{toNodeId}*', f'{fromNodeId}*')
            elif fromStrand == '+' and toStrand == '-' and \
               f'{fromNodeId}*' in G.nodes and f'{toNodeId}' in G.nodes:
               G.add_edge(f'{toNodeId}', f'{fromNodeId}*')
            elif fromStrand == '-' and toStrand == '+' and \
               f'{toNodeId}*' in G.nodes and f'{fromNodeId}' in G.nodes:
               G.add_edge(f'{toNodeId}*', f'{fromNodeId}')

        self.G = G

    def updateNodes(self):
        G = self.G

        for contig in self.firstNodeId:
            nodeId = self.firstNodeId[contig]

            nextPos = 0
            lastNode = None

            for tuple in list(nx.bfs_edges(G, source=nodeId)):
                fromNode = G.nodes[tuple[0]]
                toNode = G.nodes[tuple[1]]

                if fromNode['sample'] == self.backbone['name'] and fromNode['chr'] == contig:
                    fromNode['posStart'] = fromNode['lenBefore'] + 1
                    fromNode['posEnd'] = fromNode['posStart'] + fromNode['len'] - 1
                    nextPos = fromNode['posEnd'] + 1

                if toNode['sample'] == self.backbone['name'] and toNode['chr'] == contig:
                    toNode['posStart'] = toNode['lenBefore'] + 1
                    toNode['posEnd'] = toNode['posStart'] + toNode['len'] - 1
                    nextPos = toNode['posEnd'] + 1
                else:
                    toNode['posStart'] = nextPos
                    toNode['posEnd'] = toNode['posStart']

                lastNode = toNode

            if lastNode and lastNode['sample'] == self.backbone['name']:
                lastNode['posStart'] = lastNode['lenBefore'] + 1

        self.G = G
        self.H = self.G.to_undirected()

        #logging.info(f'whole graph: number of nodes: {len(G.nodes)}, number of edges: {len(G.edges)}')

    # e.g. posDict = {'Chr01':{'posFrom':1,'posTo':2}}
    def genSubGraph(self, sampleList, posDict):
        if not self.G:
            logging.error(f'genSubGraph(): need to run genGraph() first')
            return nx.DiGraph()

        G = self.G

        subNodes = []
        notConnectCount = 0
        connectCount = 0

        for contig in posDict:
            if contig != 'all':
                if contig not in self.firstNodeId:
                   logging.error(f'Contig {contig} not found in GFA. Contig ignored')
                   continue

                posFrom = posDict[contig]['posFrom']
                posTo = posDict[contig]['posTo']
                anyNodeId = self.firstNodeId[contig]

                for nodeId in nx.node_connected_component(self.H, anyNodeId):
                    node = G.nodes[nodeId]

                    # why needed? to-be-fix
                    if 'posStart' not in node:
                        continue

                    if sampleList and (node['rank'] != '0' and node['sample'] != self.backbone and node['sample'] not in sampleList):
                        continue

                    if (not posTo or node['posStart'] <= posTo) and (not posFrom or node['posEnd'] >= posFrom):
                        subNodes.append(nodeId)
            else:
                for contig in self.firstNodeId:
                    anyNodeId = self.firstNodeId[contig]
                    for nodeId in nx.node_connected_component(self.H, anyNodeId):
                        subNodes.append(nodeId)

        subGraph = G.subgraph(subNodes)
        self.maxNodesDisplay = int(getVar(copied, "nodes", "maxNodesDisplay"))
        if len(subGraph.nodes) > self.maxNodesDisplay:
            self.overNodeLimitSignal = 1
        if len(subGraph.nodes) == 0:
            self.emptyGraphSignal = 1
        logging.info(f'subGraph: number of nodes: {len(subGraph.nodes)}, number of edges: {len(subGraph.edges)}')

        return subGraph

    def formatNodeOutput(self, nodeId, node, showSeq=True):
        shape = getVar(copied, 'nodes',f"{node['inf']['sv_type']}_shape") if 'sv_type' in node['inf'] else getVar(copied, 'nodes', 'BB_shape')
        shape_cy = getVar(copied, 'cytoscape',f"{node['inf']['sv_type']}_shape") if 'sv_type' in node['inf'] else getVar(copied, 'cytoscape', 'BB_shape')
        sv_type = node['inf']['sv_type'] if 'sv_type' in node['inf'] else ''

        size = 1 if not node['len'] else float(f"{math.log(abs(node['len']), 10)*8 + 1:.1f}") 
        color = self.nameCols[node['sample']] if node['sample'] in self.nameCols else '#A2A2A2'
        sample = node['sample']
        pos = node['lenBefore'] if 'lenBefore' in node else 0
        inf = node['inf']['raw'] if 'raw' in node['inf'] else ''
        seqDesc = f"{node['seqDesc']}..." if node['len'] > len(node['seqDesc']) and node['seqDesc'] != '*' else node['seqDesc']

        #title = f"Resource: {node['sample']}_{node['chr']}; len: {node['len']}"
        title = f"NodeId: {nodeId}; Resource: {node['sample']}_{node['chr']}; Len: {node['len']}"
        if sample == self.backbone['name']: title += f"; Pos: {pos} - {pos + node['len'] - 1}"
        if inf: title += f"; Info: {inf}"
        if showSeq: title += f"; Seq: {seqDesc}"

        return {'color':color,'id':nodeId,'label':nodeId,'shape':shape,'size':size,'title':title,'shape_cy':shape_cy,'sv_type':sv_type,'sample':sample}

    def formatEdgeOutput(self, edge):
        return {'from':edge[0],'to':edge[1],'arrows':'to'}

    def genDrawGraphResult(self, graph, posDict):
        self.colorPalettes()

        nodes = [self.formatNodeOutput(nodeId, graph.nodes[nodeId]) for nodeId in graph.nodes]
        edges = [self.formatEdgeOutput(edge) for edge in graph.edges]

        startNodeIdList = [n for n,d in graph.in_degree() if d == 0]
        endNodeIdList = [n for n,d in graph.out_degree() if d == 0]
        for idx, nodeId in enumerate(startNodeIdList):
            if graph.nodes[nodeId]['sample'] != self.backbone['name']: continue
            if 'all' not in posDict and graph.nodes[nodeId]['chr'] not in posDict: continue
            nodes.append({'color':'green','id':f'start_{idx+1}','label':'start','shape':'star','size':20,'title':'start','shape_cy':'star'})
            edges.append({'from':f'start_{idx+1}','to':nodeId,'arrows':'to'})

        for idx, nodeId in enumerate(endNodeIdList):
            if graph.nodes[nodeId]['sample'] != self.backbone['name']: continue
            if 'all' not in posDict and graph.nodes[nodeId]['chr'] not in posDict: continue
            nodes.append({'color':'red','id':f'end_{idx+1}','label':'end','shape':'star','size':20,'title':'end','shape_cy':'star'})
            edges.append({'from':nodeId,'to':f'end_{idx+1}','arrows':'to'})

        self.drawGraphResult = {'error':False, 'nodes_data':nodes,'edges_data':edges}

    def genCyDataFromDrawGraphResult(self, drawGraphResult):
        cyData = []
        outputNodeInfo = {}

        nodes = drawGraphResult['nodes_data']
        edges = self.drawGraphResult['edges_data']

        for node in nodes:
            #size = 25 if len(node['title'].split(";")) <3 else node['title'].split(";")[2].split(" ")[2]
            #print(f"node is: {node}")
            #size = float(f"{math.log(abs(int(size)), 10)*8 + 1:.1f}") if int(size) != 1 else 1
            size = node['size']
            ele = {'data':{'type':'','id':node['id'],'name':node['id'],'weight':1,'size':size,'color':node['color'],'shape':node['shape_cy'],'title':node['title']}}
            cyData.append(ele)

            outputNodeInfo[node['id']] = ele

        for edge in edges:
            color = outputNodeInfo[edge['from']]['data']['color']
            # avoid too light color for edges
            color = f'#{int(color[1:3],16)//2:02x}{int(color[3:5],16)//2:02x}{int(color[5:],16)//2:02x}' if color[0] == '#' else color
            ele = {'data':{'type':'','source':edge['from'],'target':edge['to'],'weight':1,'color':color}}
            cyData.append(ele)

        return cyData

    def genHtml(self, graph, posDict, outHtmlPrefix=None):
        genHtmlInfo = {'vis':{'template':os.path.join(os.path.dirname(os.path.realpath(__file__)),'template','htmlTemplate.html')}, \
                       'cytoscape':{'template':os.path.join(os.path.dirname(os.path.realpath(__file__)),'template','htmlTemplate_cytoscape.html')}}

        outFile = {}
        interaction = getVar(copied, 'enable', 'interaction')
        graphLayoutModification = getVar(copied, 'enable', 'graphLayoutModification')
        addRemoveNodes = getVar(copied, 'enable', 'addRemoveNodes')
        canvas_height = getVar(copied, 'canvas', 'height')
        canvas_width = getVar(copied, 'canvas', 'width')
        maxNodesDisplay = int(getVar(copied, "nodes", "maxNodesDisplay"))

        nodes = self.drawGraphResult['nodes_data']
        edges = self.drawGraphResult['edges_data']

        #if not posDict or 'all' in posDict:
        #    title = 'All graphs'
        #else:
        #    tmp = [f"{contig}: {posDict[contig]['posFrom']} - {posDict[contig]['posTo']}" for contig in posDict]
        #    title = f"PanGraph -- {', '.join(tmp)}"

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


        # vis.js specific (?)
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

            if libraryName == 'vis' and len(graph.nodes) <= maxNodesDisplay:
                with open(templateHtml) as f_in, open(outHtml,'w') as f_out:
                    data = f_in.read()
                    data = data.replace('{%title%}', title)
                    data = data.replace('{%INFO%}', filterInfo)
                    data = data.replace('{%canvas_height%}', canvas_height)
                    data = data.replace('{%canvas_width%}', canvas_width)
                    data = data.replace('{%nodes%}', json.dumps(nodes))
                    data = data.replace('{%edges%}', json.dumps(edges))
                    f_out.write(data)
                    logging.info(f'The output HTML file is: {outHtml}')

            if libraryName == 'cytoscape' and len(graph.nodes) > maxNodesDisplay:
                with open(templateHtml) as f_in, open(outHtml,'w') as f_out:
                    data = f_in.read()
                    cyData = self.genCyDataFromDrawGraphResult(self.drawGraphResult)

                    suppFiles = ['loader.gif','images.js']
                    for suppFile in suppFiles:
                        src = os.path.join(os.path.dirname(os.path.realpath(__file__)), 'template', suppFile)
                        dest = os.path.join(self.outdir, suppFile)
                        copyfile(src, dest)

                    data = data.replace('{%title%}', title)
                    data = data.replace('{%data%}', json.dumps(cyData, indent=4, sort_keys=True))
                    f_out.write(data)
                    logging.info(f'The output HTML file is: {outHtml}')

        #return outFile
        if len(graph.nodes) <= maxNodesDisplay:
            return outFile['vis']
        else:
            return outFile['cytoscape']

    # for testing
    def run(self, targetChr, targetStart, targetEnd):
        logging.info('loading GFA ...')

        """
        self.inf = self.loadRGFA(targetChr=targetChr)
        self.neededGFA = self.inf['neededGFA']
        self.backbone = self.inf['backbone']
        self.colorPalettes()
        """

        self.loadRGFA(targetChr=targetChr)

        logging.info(f'generating graph ...')
        self.genGraph()
        #firstNodeId = self.firstNodeId[list(self.firstNodeId.keys())[0]]
        self.updateNodes()

        logging.info('generating subgraph ...')
        sampleList = list(self.inf['samples'].keys())
        posDict = {targetChr:{'posFrom':targetStart,'posTo':targetEnd}}
        subGraph = self.genSubGraph(sampleList, posDict)
        logging.info(f"number of nodes in subgraph: {len(subGraph.nodes)}")
        logging.info(f"number of edges in subgraph: {len(subGraph.edges)}")

        logging.info('generating html ...')
        self.genDrawGraphResult(subGraph, posDict)
        self.genHtml(subGraph, posDict)

        #self.nodeGeneOverlap('../../test/1-GFA/demo.bed', overlapNodeCountThreshold=1)

    def test(self):
        vcf = '../../test/2-VCF/test.vcf'
        backbone = 'demo'

        logging.info(f'Parsing vcf ...')
        self.parseVCF(vcf, backbone)
        logging.info(f'Done')

    def binSearch(self, a, x):
        i = bisect_left(a, x)
        if i :
            return i-1
        else:
            return -1

    def getGeneNodeOverlapThreadhold(self):
        num = int(getVar(copied, "nodes", "geneNodeOverlapCntThreshold"))
        return num

    def loadBedGff(self, file):
        self.bed = {}

        format = self.getBedGffFormat(file)
        logging.info(f"File format: '{format}'")
        if not format:
            logging.error(f'!!!!! loading Bed/Gff aborted: unknown format')
            return

        logging.info(f"Start to load {format.upper()} file: '{file}'")
        with open(file, 'r') as f:
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
                        logging.error(f'!!!!! Missing ID/Name/gene_id/gene_name in {format}')
                        break

                    self.bed[geneId] = {'Chr':fields[0],'Start':int(fields[3]),'End':int(fields[4]),
                                       'GeneID':geneId, 'Orientation':fields[6]}

        logging.info(f"Complete {format.upper()} file parsing ...")

    def nodeGeneOverlap(self, bed, overlapNodeCountThreshold):
        logging.info(f'finding genes which overlap more than {overlapNodeCountThreshold} nodes ...')

        results = {}

        bedInfo = {}
        # load bed
        self.loadBedGff(bed)
        bedInfo = self.bed

        # get sorted list of backbone node position
        posStartList = sorted([self.G.nodes[nodeId]['lenBefore'] for nodeId in self.G.nodes
                       if self.G.nodes[nodeId]['sample'] == self.backbone['name']])

        for geneId in bedInfo:
            gene = bedInfo[geneId]
            startIdx = self.binSearch(posStartList, gene['Start'])
            endIdx = self.binSearch(posStartList, gene['End'])

            if endIdx - startIdx + 1 > overlapNodeCountThreshold:
                 gene_chr = self.bed[geneId]['Chr']
                 gene_start = self.bed[geneId]['Start']
                 gene_end = self.bed[geneId]['End']
                 results[geneId] = {'gene_chr':gene_chr,'gene_start':gene_start,'gene_end':gene_end,
                                    'startPos':posStartList[startIdx],'endPos':posStartList[endIdx]}

        #print(results)

        logging.info(f'finding genes done')

        return results

if __name__=="__main__":
    parser = ArgumentParser(description='Generate graph html')
    parser.add_argument('-g', dest='gfa', help='the gfa file', type = str)
    parser.add_argument('-o', dest='outdir', help='the output directory', type=str)
    #parser.add_argument('-c', dest='chrlist', nargs='*', help='the name of the chromosome(s)', type=str)
    parser.add_argument('-c', dest='chr', help='the name of the chromosome(s)', type=str)
    parser.add_argument('-s', dest='start', help='start pos', type=int)
    parser.add_argument('-e', dest='end', help='end pos', type=int)
    parser.add_argument('-l', dest='samplelist', nargs='*', help='sample list', type=str)

    args = parser.parse_args()

    if None not in [args.gfa, args.outdir, args.chr]:
        #PanGraph(args.gfa, args.outdir).drawGraph(args.samplelist, args.chr, args.start, args.end)
        PanGraph(args.gfa, args.outdir, parseRGFA=False).run(args.chr, args.start, args.end)
        #PanGraph(args.gfa, args.outdir, parseRGFA=False).test()
    else:
        print('\n%s\n' % parser.print_help())
