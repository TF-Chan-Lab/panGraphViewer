#!/usr/bin/env python3

import os
import sys
import logging
import re

import csv
import networkx as nx

from argparse import ArgumentParser
from natsort import natsorted, ns

from configparser import ConfigParser

import json
import traceback

#sys.tracebacklimit = 0

DEBUG=False

from bisect import bisect_left

import threading

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

class ManageRGFA:
    backbone = 'bk1'
    sample_delim = '~'
    sample_delim2 = ':'

    def __init__(self, options):
        self.options = options

        self.nodes = {}
        self.edges = {}
        self.bb_chr = {}
        self.node_id_max = {'node':0,'gene':0}
        self.G = {}
        self.outDir = options.outdir
        self.field_delim = get_var(copied, 'nodes','SN_delim', must_have=True)

        try:
            os.makedirs(self.outDir, mode=0o777, exist_ok=True)
        except:
            pass

        self.refChroms = {}

        self.sample_list = []

    def _check_gfa(self, gfa_file):
        logging.info(f'Checking rGFA {gfa_file}')

        bb_list = {}
        gtf_list = {}
        vcf_list = {}

        with open(gfa_file) as f:
            for line in f:
                if line[0] != 'S':
                    continue

                row = line.strip().split('\t')
                nodeId = row[1]
                seq = row[2]
                tags = {val.split(':')[0]:val.split(':')[2] for val in row[3:]}
                sample, chr = tags['SN'].split(self.field_delim)[0:2]
                inf = tags['INF'] if 'INF' in tags else ''
                node_type = inf.split(self.field_delim)[0]  # GENE / SV / REF

        return {'ref_list':[], 'vcf_list':[], 'gene_list':[]}

    def _load_gfa(self, gfa_file, id_offset=0):
        logging.info(f'Loading rGFA file {gfa_file}')

        if self.node_id_max['node']:
            id_offset = self.node_id_max['node']

        # *** need to add checking

        nodes = {}
        edges = {}
        nodeIdToChr = {}
        nodeIdMap = {}
        sample_dict = {}

        with open(gfa_file) as f:
            for line in f:
                if line[0] == 'S':
                    row = line.strip().split('\t')

                    nodeId = row[1]
                    if id_offset:
                        match = re.match(r'(\D+)(\d+)', nodeId)
                        if match:
                            nodeIdMap[nodeId] = f'{match.group(1)}{int(match.group(2))+id_offset}'
                            nodeId = nodeIdMap[nodeId]

                    seq = row[2]
                    tags = {val.split(':')[0]:val.split(':')[2] for val in row[3:]}
                    # set tags default value if not found
                    if 'LN' not in tags:
                        tags['LN'] = len(seq)
                    if 'SN' not in tags:
                        tags['SN'] = ''
                    if 'SO' not in tags:
                        tags['SO'] = 0
                    #if 'SR' not in tags:
                    #    tags['SR'] = ''

                    #rank = int(tags['SR'])
                    rank = int(tags['SR']) if 'SR' in tags else 0
                    tags['SR'] = rank
                    if self.field_delim in tags['SN']:
                        bb, chr = tags['SN'].split(self.field_delim)[0:2]
                    else:
                        bb, chr = 'backbone', tags['SN']

                    sample = ''
                    sv_id = ''
                    if rank != '0' and 'INF' in tags and len(tags['INF'].split('||')) >= 7:
                        #sample = tags['INF'].split('||')[6].split('~')
                        sv_id = tags['INF'].split('||')[6]
                        sample_dict.update({s:1 for s in tags['INF'].split('||')[7].split('~')})

                    start, end = int(tags['SO']), int(tags['SO'])+int(tags['LN'])-1

                    if 'INF' in tags:
                        inf = tags['INF']
                        if tags['INF'].split(self.field_delim)[0] == 'ANNOT':
                           node_type = 'GENE'
                        else:
                           node_type = 'SV'
                    else:
                        node_type = 'REF'
                        inf = ''

                    if 'LN' in tags and int(tags['LN']) != len(seq) and node_type not in ['GENE']:
                        #print(f"ERROR: LN != len(seq) ({tags['LN']}, {len(seq)}), {tags}")
                        logging.error(f"LN != len(seq) ({tags['LN']}, {len(seq)}), {tags}")

                    try:
                       size = int(tags['LN']) if node_type != 'GENE' else int(tags['INF'].split(self.field_delim)[2])
                    except:
                       size = 0

                    node_info = {'nodeType':node_type,'nodeId':nodeId,'nodeId2':nodeId,'sv_id':sv_id,'start':start,'end':end,'seq':seq,'LN':int(tags['LN']),'SN':tags['SN'],'SO':int(tags['SO']),'SR':int(tags['SR']), 'INF':inf, 'size':size}
                    bb_chr = f'{bb}{self.field_delim}{chr}'

                    try:
                       nodes[bb_chr]['all'][nodeId] = node_info
                    except:
                       nodes[bb_chr] = {'all':{},'REF':{},'SV':{},'GENE':{}}
                       nodes[bb_chr]['all'][nodeId] = node_info

                    #nodes[bb_chr][node_type] = nodeId

                    # update self.node_id_max
                    nodeId_int = int(re.findall(r'\d+', nodeId)[0])
                    if nodeId_int > self.node_id_max['node']:
                        self.node_id_max['node'] = nodeId_int

                    nodeIdToChr[nodeId] = bb_chr

                    # added tentatively
                    if bb_chr not in edges:
                        edges[bb_chr] = {}

                    #if node_type == 'SV' and sample not in sample_dict:
                    #    sample_dict[sample] = 1

                elif line[0] == 'L':
                    row = line.strip().split('\t')
                    fromNodeId, fromStrand, toNodeId, toStrand, cigar = row[1:6]

                    if id_offset:
                      fromNodeId = nodeIdMap[fromNodeId]
                      toNodeId = nodeIdMap[toNodeId]

                    tags = {'from_strand':fromStrand, 'to_strand':toStrand, 'cigar':cigar}
                    tags.update({val.split(':')[0]:int(val.split(':')[2]) for val in row[6:]})

                    #edges.append({'from':fromNodeId,'fromStrand':'+','to':toNodeId,'toStrand':'+','tag':tags})
                    edge_info = {'from':fromNodeId,'fromStrand':'+','to':toNodeId,'toStrand':'+','tag':tags}
                    fromNodeChr = nodeIdToChr[fromNodeId]
                    try:
                        edges[fromNodeChr][(fromNodeId, toNodeId)] = edge_info
                    except:
                        edges[fromNodeChr] = {}
                        edges[fromNodeChr][(fromNodeId, toNodeId)] = edge_info

        # backbone_node_list (any need?)
        #self.backbone_nodes = {nodes[id]['SO']:nodes[id]['nodeId'] for id in nodes if nodes[id]['SR'] == self.backbone_SR}

        #self.results = {}
        #for bb_chr in nodes.keys():
        #    self.results[bb_chr] = {'bb_chr':bb_chr,'nodes':nodes[bb_chr],'edges':edges[bb_chr]}

        for bb_chr in nodes.keys():
            for node_id in nodes[bb_chr]['all']:
                node_type = nodes[bb_chr]['all'][node_id]['nodeType']
                nodes[bb_chr][node_type][node_id] = 1

            self.nodes[bb_chr] = nodes[bb_chr]
            self.edges[bb_chr] = edges[bb_chr]

            self.G[bb_chr] = nx.MultiDiGraph()
            for node_id in nodes[bb_chr]['all']:
                self.G[bb_chr].add_node(node_id, **nodes[bb_chr]['all'][node_id])
            for edge in edges[bb_chr]:
                self.G[bb_chr].add_edge(edge[0], edge[1])

        """
        # create network
        self.G[bb_chr] = nx.MultiDiGraph()
        for node_id in nodes[bb_chr]['all']:
            self.G[bb_chr].add_node(node_id, **nodes[bb_chr]['all'][node_id])
        for edge in edges[bb_chr]:
            self.G[bb_chr].add_edge(edge[0], edge[1])
        """

        self.sample_list = natsorted(sample_dict.keys())

    def _save_gfa(self, gfa=None):
        gfa = os.path.join(self.outDir, gfa if gfa else f'{self.prefix}.gfa')

        logging.info(f'Generating rGFA file {gfa}')

        if os.path.isfile(gfa):
            os.remove(gfa)

        with open(gfa, 'w') as fd:
            for bb_chr in self.nodes.keys():
                nodes = self.nodes[bb_chr]
                if not nodes['all']:
                    continue

                if bb_chr not in self.G:
                    continue

                edges = self.G[bb_chr].edges
                bb, chr = bb_chr.split(self.field_delim)

                for nodeId in nodes['all']:
                    node = nodes['all'][nodeId]

                    if node['nodeType'] == 'REF':
                        #ref = self.refChroms[bb]['fasta']
                        #if not node['seq']:
                        #    node['seq'] = self.get_seq(node['nodeType'],chr,node['start'],node['end'],env={},nodeId=nodeId,ref=ref)
                        output = ['S',node['nodeId2'],node['seq'],f"LN:i:{node['LN']}",f"SN:Z:{node['SN']}",f"SO:i:{node['SO']}",f"SR:i:{node['SR']}"]
                    elif node['nodeType'] == 'SV':
                        #ref = self.refChroms[bb]['fasta']
                        #if not node['seq']:
                        #    node['seq'] = self.get_seq(node['nodeType'],chr,node['start'],node['end'],env={},nodeId=nodeId,ref=ref) if not node['seq'] else node['seq']
                        output = ['S',node['nodeId2'],node['seq'],f"LN:i:{node['LN']}",f"SN:Z:{node['SN']}",f"SO:i:{node['SO']}",f"SR:i:{node['SR']}",f"INF:Z:{node['INF']}"]
                    elif node['nodeType'] == 'GENE':
                        #output = ['S', node['nodeId2'],'',f"LN:i:{node['LN']}",f"SN:Z:{node['SN']}",f"SO:i:{node['SO']}",f"SR:i:{node['SR']}",f"INF:Z:ANNOT{self.field_delim}GENE{self.field_delim}{node['size']}"]
                        output = ['S', node['nodeId2'],'',f"LN:i:{node['LN']}",f"SN:Z:{node['SN']}",f"SO:i:{node['SO']}",f"SR:i:{node['SR']}",f"INF:Z:ANNOT{self.field_delim}GENE{self.field_delim}{node['end']}"]

                    fd.write('%s\n' % "\t".join([str(x) for x in output]))

            for bb_chr in self.nodes.keys():
                if bb_chr not in self.G:
                    continue

                edges = self.G[bb_chr].edges
                for edge in edges:
                    if edge[0] not in self.nodes[bb_chr]['all'] or edge[1] not in self.nodes[bb_chr]['all']:
                        continue

                    node_from = self.nodes[bb_chr]['all'][edge[0]]
                    node_to = self.nodes[bb_chr]['all'][edge[1]]
                    row = ['L',node_from['nodeId2'],'+',node_to['nodeId2'],'+','0M']
                    fd.write('\t'.join(row)+'\n')

        os.chmod(gfa, 0o666)

    def _filter_sv(self, bb_chr):
        # need to check if self.nodes[bb_chr]['GENE'] is empty

        sv_filtered = {bb_chr:[]}
        node_all = self.nodes[bb_chr]['all']
        node_gene = self.nodes[bb_chr]['GENE']
        sv_info = self.nodes[bb_chr]['SV']
        bed_info = {node_all[node_id]['start']:node_all[node_id]['end'] for node_id in node_gene}
        bed_start_sorted = sorted(bed_info.keys())
        for node_id in sv_info:
            sv_start = node_all[node_id]['start']
            sv_end = node_all[node_id]['end']
            sv_start_idx = bisect_left(bed_start_sorted, sv_start)
            sv_end_idx = bisect_left(bed_start_sorted, sv_end)

            if sv_start_idx != sv_end_idx:
                sv_filtered[bb_chr].append(node_id)
            else:
                if sv_start_idx >= len(bed_start_sorted):
                    sv_start_idx -= 1
                bed_start = bed_start_sorted[sv_start_idx]
                bed_end = bed_info[bed_start]
                if sv_start <= bed_end and sv_end >= bed_start:
                    sv_filtered[bb_chr].append(node_id)

        return sv_filtered

    def _save_vcf(self, bb, nodes=None, vcf=None, overlap_gene_only=False):
        vcf = os.path.join(self.outDir, vcf if vcf else f'{self.prefix}.vcf')

        logging.info(f"Generating vcf file {vcf}")
        if overlap_gene_only:
            if not [bb_chr for bb_chr in self.nodes if self.nodes[bb_chr]['GENE']]:
                logging.warning(f'Annotation is not found. SV are not filtered by Genes')
            else:
                logging.info(f'SV are filtered by genes [-overlap_gene_only 1]')
        else:
            logging.info(f'SV are NOT filtered by genes')

        with open(vcf, 'w') as out_csv:
            out_fieldnames = ['#CHROM','POS','ID','REF','ALT','QUAL','FILTER','INFO','FORMAT'] + self.sample_list

            writer = csv.DictWriter(out_csv, fieldnames=out_fieldnames, delimiter='\t', extrasaction='ignore', dialect='unix', quoting=csv.QUOTE_NONE)
            writer.writeheader()

            if not nodes:
                nodes = {}
                for bb_chr in self.nodes:
                    if bb_chr.startswith(bb):
                        if overlap_gene_only and self.nodes[bb_chr]['GENE']:
                            nodes[bb_chr] = self._filter_sv(bb_chr)[bb_chr]
                        else:
                            nodes[bb_chr] = self.nodes[bb_chr]['SV']

            count = 0
            for bb_chr in nodes:
                for node_id in sorted(nodes[bb_chr], key=lambda x: self.nodes[bb_chr]['all'][x]['start']):
                    sv_node = self.nodes[bb_chr]['all'][node_id]

                    if (sv_node['SN'].split(self.field_delim)[0] != bb):
                        continue

                    nodeId = sv_node['nodeId']
                    nodeId2 = sv_node['nodeId2']
                    varSeq = sv_node['seq']
                    fields = sv_node['INF'].split(self.field_delim)

                    sv_type = fields[1]
                    sv_chr = fields[2]
                    sv_pos = int(fields[3])
                    sv_end = int(fields[5])
                    sv_len = len(sv_node['seq']) if sv_type == 'INS' else int(sv_end)-int(sv_pos)
                    sv_id = fields[6]
                    sv_sample_dict = {entry:'0/1' for entry in fields[7].split(self.sample_delim) if entry}
                    # if haplotype is stored
                    #sv_sample_dict = {entry.split(self.sample_delim2)[0]:entry.split(self.sample_delim2)[1] for entry in fields[7].split(self.sample_delim) if entry}

                    results = {}
                    results['#CHROM'] = sv_chr
                    results['POS'] = sv_pos
                    results['ID'] = sv_id
                    results['REF'] = 'N'
                    results['ALT'] = sv_node['seq'] if sv_type == 'INS' else f'<{sv_type}>'
                    results['QUAL'] = '.'
                    results['FILTER'] = 'PASS'
                    results['INFO'] = f'SVTYPE={sv_type};SVLEN={sv_len};END={sv_end}'
                    results['FORMAT'] = 'GT'
                    results.update({sample:sv_sample_dict[sample] if sample in sv_sample_dict.keys() else '0/0' for sample in self.sample_list})

                    writer.writerow(results)
                    count += 1

        logging.info(f"{count} SV are extracted from {bb}")
        os.chmod(vcf, 0o666)

    def _save_bed(self, bb=None, bed=None):
        bed = os.path.join(self.outDir, bed if bed else f'{self.prefix}.bed')

        logging.info(f"Generating bed file {bed} from backbone {bb}")

        with open(bed, 'w') as out_bed:
            out_fieldnames = ['chr', 'start', 'end', 'node_id', 'ref']
            writer = csv.DictWriter(out_bed, fieldnames=out_fieldnames, delimiter='\t', extrasaction='ignore', dialect='unix', quoting=csv.QUOTE_NONE)
            writer.writeheader()

            count = 0
            for bb_chr in self.nodes:
                if bb and not bb_chr.startswith(bb):
                    continue

                backbone, chr = bb_chr.split(self.field_delim)
                for node_id in self.nodes[bb_chr]['GENE']:
                    bed_node = self.nodes[bb_chr]['all'][node_id]

                    nodeId = bed_node['nodeId']

                    results = {}
                    results['chr'] = chr
                    results['start'] = bed_node['start']
                    results['end'] = int(results['start']) + bed_node['size']
                    results['node_id'] = nodeId
                    results['ref'] = backbone

                    writer.writerow(results)
                    count += 1

        logging.info(f"{count} genes are extracted from backbone {bb}")

        os.chmod(bed, 0o666)

    def get_nodes(self, bb_chr=None):
        return self.nodes

    def get_edges(self, bb_chr=None):
        return self.edges

    def _get_node_id(self, node_type='node'):
        assert node_type in self.node_id_max

        prefix = 'g' if node_type == 'gene' else 's'

        self.node_id_max[node_type] += 1

        return f'{prefix}{self.node_id_max[node_type]}'

    def _add_annot(self, annot_file, bb):
        logging.info(f'Adding gff file {annot_file} to backbone {bb}')

        not_found_bb_chr = {}
        with open(annot_file) as f:
            for line in f:
                if not line.strip() or line[0] == '#':
                    continue
                fields = line.strip().split('\t')
                if len(fields) < 3 or fields[2] != 'gene':
                    continue
                chr,name,start,end,info = fields[0],fields[1],int(fields[3]),int(fields[4]),fields[8]
                new_node_id = self._get_node_id('gene')
                bb_chr = f'{bb}{self.field_delim}{chr}'

                if bb_chr not in self.nodes:
                    not_found_bb_chr[bb_chr] = 1
                    continue
                    #self.nodes[bb_chr] = {'all':{},'REF':{},'SV':{},'GENE':{}}

                #node = {'nodeType':'GENE','end':end,'size':end-start+1,'nodeId':new_node_id,'nodeId2':new_node_id,'LN':end-start+1,'SN':f"{bb}{self.field_delim}{chr}",'SO':start,'SR':0,'INF':''}
                node = {'nodeType':'GENE','start':start,'end':end,'size':end-start+1,'nodeId':new_node_id,'nodeId2':new_node_id,'LN':end-start+1,'SN':f"{bb}{self.field_delim}{chr}",'SO':start,'SR':0,'INF':''}
                self.nodes[bb_chr]['all'][new_node_id] = node
                self.nodes[bb_chr]['GENE'][new_node_id] = 1

        entries_added = 0
        for bb_chr in self.nodes:
            if not bb_chr.startswith(bb):
                continue
            if self.nodes[bb_chr]['GENE']:
                entries_added += len(self.nodes[bb_chr]['GENE'])

        if not_found_bb_chr:
            logging.warning(f'Annotation of {len(not_found_bb_chr)} chroms/contigs not used in GFA have been ignored')

        logging.info(f'{entries_added} gene entries added')

    def _remove_annot(self, backbone):
        logging.info(f"Removing annotation from backbone {backbone}")

        for bb_chr in self.nodes:
            if not bb_chr.startswith(backbone):
                continue

            for node_id in self.nodes[bb_chr]['GENE']:
                del self.nodes[bb_chr]['all'][node_id]
            del self.nodes[bb_chr]['GENE']

    def _remove_backbone(self, backbone):
        logging.info(f"Removing backbone {backbone}")

        temp = []
        for bb_chr in self.nodes:
            if bb_chr.startswith(backbone):
                temp.append(bb_chr)

        for bb_chr in temp:
            del self.G[bb_chr]
            del self.nodes[bb_chr]

    def test(self):
        self._load_gfa('by_release_fix_test/W05_vcf2rGFA.gfa')
        self._load_gfa('by_release_fix_test/ZH13_vcf2rGFA.gfa')
        self._load_gfa('by_release_fix_test/Wm82_vcf2rGFA.gfa')

        self._add_annot('gff/Gsoja.W05.gene.gff','W05')
        self._add_annot('gff/Gmax_ZH13v2.ChrRenamed_modified.gff','ZH13')
        self._add_annot('gff/Gmax_508_Wm82.a4.v1.gene_exons.gff3','Wm82')

        self._save_gfa('merged_W05_ZH13_Wm82.gfa')

        self._save_bed(bb='W05', bed='W05_working_save_annot.bed')
        self._save_bed(bb='ZH13', bed='ZH13_working_save_annot.bed')
        self._save_bed(bb='Wm82', bed='Wm82_working_save_annot.bed')

        self._save_vcf(bb='W05', vcf='W05_working_save_vcf.vcf')
        self._save_vcf(bb='ZH13', vcf='ZH13_working_save_vcf.vcf')
        self._save_vcf(bb='Wm82', vcf='Wm82_working_save_vcf.vcf')

        self._save_vcf(bb='W05', vcf='W05_working_save_vcf_filtered.vcf', overlap_gene_only=True)
        self._save_vcf(bb='ZH13', vcf='ZH13_working_save_vcf_filtered.vcf', overlap_gene_only=True)
        self._save_vcf(bb='Wm82', vcf='Wm82_working_save_vcf_filtered.vcf', overlap_gene_only=True)

    # need to pass through check_param() first
    def run(self):
        logging.info(f'Action: {self.options.action}')
        if self.options.action == 'merge_gfa':
            for gfa in self.options.in_gfa:
                self._load_gfa(gfa)
            self._save_gfa(self.options.out_gfa)
        elif self.options.action == 'add_annot':
            for gfa in self.options.in_gfa:
                self._load_gfa(gfa)
            for entry in self.options.in_gff:
                temp = entry.split(',')
                gff, backbone = temp[0], temp[1]
                self._add_annot(gff, backbone)
            self._save_gfa(self.options.out_gfa)
        elif self.options.action == 'get_sv':
            for gfa in self.options.in_gfa:
                self._load_gfa(gfa)
            for entry in self.options.out_vcf:
                temp = entry.split(',')
                vcf, backbone = temp[0], temp[1]
                self._save_vcf(bb=backbone, vcf=vcf, overlap_gene_only=self.options.overlap_gene_only)
        elif self.options.action == 'get_annot':
            for gfa in self.options.in_gfa:
                self._load_gfa(gfa)
            for entry in self.options.out_bed:
                temp = entry.split(',')
                bed, backbone = temp[0], temp[1]
                self._save_bed(bb=backbone, bed=bed)
        elif self.options.action == 'remove_annot':
            for gfa in self.options.in_gfa:
                self._load_gfa(gfa)
            for backbone in self.options.backbone:
                self._remove_annot(backbone=backbone)
            self._save_gfa(self.options.out_gfa)
        elif self.options.action == 'remove_backbone':
            for gfa in self.options.in_gfa:
                self._load_gfa(gfa)
            for backbone in self.options.backbone:
                self._remove_backbone(backbone=backbone)
            self._save_gfa(self.options.out_gfa)

def check_param(args):
    error = None
    """
    if self.options.action == 'merge_gfa':
        if not self.options.in_gfa:
            return False
    elif self.options.action == 'add_annot':
        return False
    elif self.options.action == 'get_sv':
        return False
    elif self.options.action == 'get_annot':
        return False
    elif self.options.action == 'remove_annot':
        return False
    elif self.options.action == 'remove_backbone':
        return False
    """

    return error

if __name__=="__main__":
    parser = ArgumentParser(description='Convert a vcf file to an rGFA file')
    parser.add_argument('--version', action='version', version='1.0')
    parser.add_argument('-outdir', dest='outdir', help='the output directory', type=str)

    action_choices=['merge_gfa','add_annot','get_sv','get_annot','remove_sv','remove_annot','remove_backbone']
    parser.add_argument('-action', dest='action', help='Action to be done on GFA. Correponding options needed', choices=action_choices, type=str, default='add_vcf')
    parser.add_argument('-in_gfa', dest='in_gfa', nargs='+', help='input GFA', type=str)
    parser.add_argument('-in_gff', dest='in_gff', nargs='*', help='add GFF of specified backbone [-in_gff gff1,backbone1 gff2,backbone2 ...]', type=str)
    parser.add_argument('-out_gfa', dest='out_gfa', help='output GFA', type=str)
    parser.add_argument('-out_vcf', dest='out_vcf', nargs='*', help='output VCF of specified backbone [-out_vcf vcf1,backbone1 vcf2,backbone2]', type=str)
    parser.add_argument('-out_bed', dest='out_bed', nargs='*', help='output BED of specified backbone [-out_bed bed1,backbone1 bed2,backbone2]', type=str)

    parser.add_argument('-backbone', dest='backbone', nargs='+', help='the name of the backbone sample', type=str)
    parser.add_argument('-overlap_gene_only', dest='overlap_gene_only', help='output SV that overlap with genes only [1]', type=int, default=1)

    args = parser.parse_args()
    error = check_param(args)
    if None not in [args.outdir] and not error:
        if args.action == 'test':
            ManageRGFA(args).test()
        else:
            ManageRGFA(args).run()
    else:
        print('\n%s\n' % parser.print_help())
        if error:
            print(error)
