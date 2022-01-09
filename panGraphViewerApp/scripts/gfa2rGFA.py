#!/usr/bin/env python3

import os
import sys
import logging
import re

import networkx as nx

from argparse import ArgumentParser
import json

# to_be_fix
try:
    from scripts.utilities import *
except ModuleNotFoundError:
    try:
        from utilities import *
    except ModuleNotFoundError:
        from pangraphviewer.utilities import *

#============================= Function =================================
##logging info
DEBUG="" #change it when debugging
logFormat = "%(asctime)s [%(levelname)s] %(message)s"
level = "DEBUG" if DEBUG != "" else "INFO"
logging.basicConfig( stream=sys.stderr, level=level, format=logFormat )
#========================================================================

class GFA2rGFA:
    def __init__(self, in_gfa, out_rgfa):
        self.in_gfa = in_gfa
        self.out_rgfa = out_rgfa
        self.SN_delim = getVar(copied, 'nodes', 'SN_delim', mustHave=True)
        self.checknLines = getVar(copied, 'nodes', 'checknLines', mustHave=True)

    #@staticmethod
    def checkGfaFormat(self):
        formats = {'RGFA':{'headers':'SLW','tags':['LN','SN','SO','SR','INF'],'must_headers':'S','must_tags':['SN','SO','SR']},
                   'GFA1':{'headers':'#HSLCPW','tags':['LN','RC','FC','KC','SH','UR'],'must_headers':'S','must_tags':[]}}
        headers, tags = {}, {}

        try:
            nLines = max(int(self.checknLines), 10)
        except:
            nLines = 10

        with open(self.in_gfa) as f:
            n = 0
            for line in f:
                if not line.strip() or line[0] == '#':
                    continue

                fields = line.strip().split()

                headers[fields[0]] = 1
                if line[0] == 'S':
                    tags.update({val.split(':')[0]:1 for val in fields[3:]})

                n += 1
                if n >= nLines:
                    break

            for name in formats:
                format = formats[name]
                if set(format['must_headers']).issubset(set(headers.keys())) and \
                   set(format['must_tags']).issubset(tags.keys()) and \
                   set(headers.keys()).issubset(set(format['headers'])) and \
                   set(tags.keys()).issubset(set(format['tags'])):
                    if name == 'RGFA':
                        return 1
                    if name == 'GFA1':
                        return 2
            return 3

    #@staticmethod
    def convert(self):
        return_code = {'error':''}

        defaultName = {'backbone':'backbone','non-backbone':'non-backbone','contig':'contig'}
        S, L, P  = {}, [], {}

        rect = 0
        with open(self.in_gfa) as f:
            lineNum = 0
            try:
                for line in f:
                    lineNum += 1

                    if not line.strip() or line[0] == '#':
                        continue
                    row = line.strip().split()
                    if not row or row[0] not in ['L','S','P']:
                        continue

                    if row[0] == 'S':
                        if len(row) < 3:
                            logging.error(f'missing field in line {lineNum}. Abort!')
                            rect = 4

                        nodeId, seq = row[1], row[2]

                        S[nodeId] = {'seq':seq}
                        tag = {}
                        for col in row[3:]:
                            tag[col.split(':')[0]] = col.split(':')[2]

                        if seq == '*' and 'LN' not in tag:
                            logging.error(f'both seq and LN tag are missing in line {lineNum}. Abort!')
                            rect = 5

                        S[nodeId]['len'] = int(tag['LN']) if 'LN' in tag else len(S[nodeId]['seq'])
                    elif row[0] == 'L':
                        try:
                            cigar = row[5]
                        except:
                            cigar = ''
                        L.append({'fromNodeId':f'{row[1]}{row[2]}','toNodeId':f'{row[3]}{row[4]}','cigar':cigar})
                    elif row[0] == 'P':
                        P[row[1]] = {'path':row[2].split(','),'cigar':row[3].split(',')}
            except:
                logging.error(f'The file is in GFA v1, but error occurs during conversion at line {lineNum}. Abort!')
                rect = 6 

        # Path is the must line
        if not P:
            logging.error('Path info is missing from the input GFA1 file, which is necessary for conversion. Abort!')
            rect = 7

        if rect in [4, 5, 6, 7]:
            sys.exit(rect) 

        # update S
        for pathName in P:
            for idx, nodeId2 in enumerate(P[pathName]['path']):
                nodeId, strand = nodeId2[:-1], nodeId2[-1]
                if strand == '-':
                    newId = f'{nodeId}*'
                    S[newId] = {'seq':rev_comp(S[nodeId]['seq']),'len':S[nodeId]['len']}
                    del S[nodeId]

                    P[pathName]['path'][idx] = f'{newId}+'

        # update L
        for idx, edge in enumerate(L):
            fromNodeId, toNodeId = edge['fromNodeId'], edge['toNodeId']
            # eg. 12- to 12*+, 12+ to 12*-

            if fromNodeId not in S:
                if fromNodeId[-1:] == '-' and f'{fromNodeId[:-1]}*' in S:
                    fromNodeId = f'{fromNodeId[:-1]}*+'
                elif fromNodeId[-1:] == '+' and f'{fromNodeId[:-1]}*' in S:
                    fromNodeId = f'{fromNodeId[:-1]}*-'
                edge['fromNodeId'] = fromNodeId
            if toNodeId not in S:
                if toNodeId[-1:] == '-' and f'{toNodeId[:-1]}*' in S:
                    toNodeId = f'{toNodeId[:-1]}*+'
                elif toNodeId[-1:] == '+' and f'{toNodeId[:-1]}*' in S:
                    toNodeId = f'{toNodeId[:-1]}*-'
                edge['toNodeId'] = toNodeId

        G = nx.DiGraph()
        for nodeId in S:
            G.add_node(nodeId, seq=S[nodeId])

        for edge in L:
            fromNodeId, toNodeId = edge['fromNodeId'][:-1], edge['toNodeId'][:-1]
            G.add_edge(fromNodeId, toNodeId)

        # update nodes
        H = G.to_undirected()
        for pathName in P:
            ele = P[pathName]
            lenBefore = 0
            firstNodeId = None
            for pathNodeId in ele['path']:
                node = S[pathNodeId[:-1]]
                if not firstNodeId:
                    firstNodeId = nodeId
                node['sample'] = defaultName['backbone']
                node['contig'] = pathName
                node['lenBefore'] = lenBefore
                node['rank'] = 0
                lenBefore += node['len']

            for pathNodeId in nx.node_connected_component(H, firstNodeId):
                node = S[pathNodeId]
                if 'sample' in node: continue

                node['sample'] = defaultName['non-backbone']
                node['contig'] = pathName
                node['lenBefore'] = 0
                node['rank'] = 1

        # handle unconnected nodes
        for nodeId in S:
            node = S[nodeId]
            if 'sample' in node: continue

            node['sample'] = defaultName['non-backbone']
            node['contig'] = defaultName['contig']
            node['lenBefore'] = 0
            node['rank'] = 1

        with open(self.out_rgfa, 'w') as f:
            for nodeId in S:
                node = S[nodeId]

                field = []
                field.append('S')
                field.append(nodeId)
                field.append(node['seq'])
                field.append(f"LN:i:{node['len']}")
                field.append(f"SN:Z:{node['sample']}{self.SN_delim}{node['contig']}")
                field.append(f"SO:i:{node['lenBefore']}")
                field.append(f"SR:i:{node['rank']}")
                print('\t'.join(field),file=f)

            for edge in L:
                fromNodeId, fromStrand, toNodeId, toStrand, cigar = \
                    edge['fromNodeId'][:-1], edge['fromNodeId'][-1], edge['toNodeId'][:-1], edge['toNodeId'][-1], edge['cigar']

                field = ['L', fromNodeId, fromStrand, toNodeId, toStrand, cigar]
                print('\t'.join(field), file=f)

        return {}

if __name__=="__main__":
    parser = ArgumentParser(description='Convert a GFA(v1) file to an rGFA file')
    parser.add_argument('--version', action='version', version='1.0')
    parser.add_argument('-in_gfa', dest='in_gfa', help='input gfa', type = str)
    parser.add_argument('-out_rgfa', dest='out_rgfa', help='output rgfa', type=str)

    args = parser.parse_args()
    if None not in [args.in_gfa, args.out_rgfa]:
        run = GFA2rGFA(args.in_gfa, args.out_rgfa)
        if run.checkGfaFormat() == 1:
            logging.info("The input file is in rGFA format.")
        if run.checkGfaFormat() == 2:
            logging.info("The input file is in GFA1 format and we are converting it to rGFA. Please wait!")
            run.convert()
        if run.checkGfaFormat() == 3:
            logging.error("The input file is in an unknown format. Please check!")
    else:
         print('\n%s\n' % parser.print_help())
