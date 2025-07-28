from django.shortcuts import render

from django.http import JsonResponse, HttpResponse
from rest_framework.decorators import api_view
import json

from .forms import UploadForm

import time

from django.conf import settings

from pangraphviewer.utilities import *
#from pangraphviewer.panGraph import *
from pangraphviewer.panGraphWeb import *
from pangraphviewer.gfa2rGFA import *

from .models import *

from django.contrib.auth.decorators import login_required
import datetime


from os import listdir
from os.path import isfile, isdir, join
from django.db import IntegrityError, transaction

work_base_dir = None

upload_file_ext = {'gfa':['.gfa','.rgfa'],
                   'vcf':['.vcf'],
                   'fasta':['.fa','.fasta','.fna'],
                   'bed':['.bed','.gtf','.gff','.gff3']}

@login_required
def get_work_dir(request, folder_list=[], file_name=''):
    global work_base_dir

    if not work_base_dir:
        profile = Profile.objects.get_or_create(username=request.user.username)
        work_base_dir = profile[0].work_base_dir if profile[0].work_base_dir else getVar(copied, 'web', 'work_dir', mustHave=True)

    if sys.platform == 'win32' and work_base_dir[0] == '/':
        drive = os.path.dirname(os.path.realpath(__file__)).split(':')[0]
        work_base_dir = f'{drive}:{work_base_dir}'.replace('/','\\')
        profile = Profile.objects.get_or_create(username=request.user.username)
        profile[0].work_base_dir = work_base_dir
        profile[0].save()

    work_dir = os.path.join(work_base_dir, request.user.username, *folder_list)
    if os.path.isdir(work_dir):
        if not os.access(work_dir, os.W_OK):
            raise Exception(f'Working directory {work_dir} exists but is inaccessible')
    else:
        makedirs(work_dir)

    return os.path.join(work_dir, file_name)

@login_required
def get_uploaded_list(request):
    file_type = request.POST.get('file_type','')

    work_dir = get_work_dir(request, [file_type])
    files = sorted([name for name in listdir(work_dir) if isfile(join(work_dir, name)) and f".{name.split('.')[-1].lower()}" in upload_file_ext[file_type]])
    status = 200

    return JsonResponse({'files':files,'work_dir':work_dir}, safe=False, status=status)

@login_required
def getdata(request):
    input_type = request.POST.get('input_type','')
    gfa = request.POST.get('gfa','')
    vcf = request.POST.get('vcf','')
    fasta = request.POST.get('fasta','')
    backbone = request.POST.get('backbone','')
    chr = request.POST.get('chr','')
    start = int(request.POST.get('start',0)) if request.POST.get('start') else None
    end = int(request.POST.get('end',0)) if request.POST.get('end') else None
    node_ids = request.POST.get('node_ids','')
    #sampleList = request.POST.getlist('sample_list[]','')
    sampleList = None

    data = None
    error = ''
    warning = ''

    # check for input
    if input_type == 'gfa':
        if not gfa or not chr:
            data = {'error':True, 'msg': 'missing value', 'nodes_data':[], 'edges_data':[]}
    elif input_type == 'vcf':
        if not vcf or not backbone or not chr:
            data = {'error':True, 'msg': 'missing value', 'nodes_data':[], 'edges_data':[]}
    if data: return JsonResponse(data, safe=False, status=400 if data['error'] else 200)

    if input_type == 'vcf':
        if fasta: fasta = get_work_dir(request, ['fasta'], fasta)
        out_dir = get_work_dir(request, ['gfa'])
        prefix = f'{backbone}_{os.path.splitext(os.path.basename(vcf))[0]}_{chr}'
        vcf = get_work_dir(request, ['vcf'], vcf)
        results = convert_vcf_to_gfa(vcf, chr, backbone, out_dir, prefix, fasta)
        gfa = get_work_dir(request, ['gfa'], f'{prefix}.gfa')
    else:
        gfa = get_work_dir(request, ['gfa'], gfa)

    if node_ids:
        nodeIdDict = {id.strip():1 for id in node_ids.split()}
        graph = PanGraph(gfa, outdir=get_work_dir(request), nodeIdDict=nodeIdDict)
        drawGraphResult = graph.drawGraphByNodeId(sampleList, targetChr=chr, nodeIdDict=nodeIdDict, isGenHtml=False)
    else:
        graph = PanGraph(gfa, outdir=get_work_dir(request), parseRGFA=False)
        drawGraphResult = graph.drawGraph(sampleList, chr, start, end, isGenHtml=False)

    nodeCount = len(drawGraphResult['nodes_data'])
    if nodeCount > graph.maxNodesLimit:
        error = f"The num. of nodes in the selected region is '{nodeCount}', which is bigger than the limit '{graph.maxNodesLimit}' in settings. Please modify Start/End Position to limit the number of nodes."

    if input_type == 'vcf':
        warning = f"The vcf file has been converted to rGFA format as {prefix}.gfa, and can be found in 'uploaded (r)GFA file' list"

    status = 400 if drawGraphResult['error'] else 200

    cyData = graph.genCyDataFromDrawGraphResult(drawGraphResult)
    colors = {node['sample']:node['color'] for node in drawGraphResult['nodes_data'] if 'sample' in node and ',' not in node['sample']}
    shapes = {node['sv_type']:node['shape_cy'] for node in drawGraphResult['nodes_data'] if 'sv_type' in node and node['sv_type']}
    hasReversed = True if [node for node in drawGraphResult['nodes_data'] if node['id'][-1] == '*'] else False
    legend = {'colors':colors, 'shapes':shapes, 'has_reversed': hasReversed}

    return JsonResponse({'error':error,'warning':warning,'cyData':cyData,'gfa':os.path.basename(gfa),'legend':legend}, safe=False, status=status)

@login_required
def parse_gfa(request):
    start_time = time.time()

    gfa = request.POST.get('gfa','')
    gfa = get_work_dir(request,['gfa'], gfa)

    now = datetime.datetime.now()
    graph = PanGraph(gfa, outdir=get_work_dir(request))

    warning, error = '', ''
    if graph.illegalrGFA:
        warning = "The format of 'SN:Z' in the input rGRF file is not the one we want. Although we can display the graph, many features would not show. Please refer to the 'Manual' in Help if you want a better experience in using this tool."
    elif graph.neededGFA:
        error = "Unknown GFA format"

    results = {}
    results['status'] = 0
    results['backbone_info'] = graph.backbone
    #results['backbone'] = list(graph.backbone.keys())
    results['chr'] =  []
    results['gfa'] = gfa
    results['time'] = time.time() - start_time
    results['warning'] = warning
    results['error'] = error

    return JsonResponse(results, status=200)

@login_required
def parse_vcf(request):
    start_time = time.time()
    vcf = request.POST.get('vcf','')
    vcf_backbone = request.POST.get('vcf_backbone','')
    vcf = get_work_dir(request, ['vcf'], vcf)

    data = PanGraph(None, outdir=get_work_dir(request), parseRGFA=False).parseVCF(vcf, vcf_backbone)
    error = data['error']
    warning = data['warning']

    needFasta = 0 if data['refChroms'] else 1
    results = {}
    results['status'] = 0
    results['backbone'] = [data['backbone']['name']]
    results['chr'] = data['backbone']['contigs']
    results['vcf'] = vcf
    results['time'] = time.time() - start_time
    results['warning'] = warning
    results['error'] = error
    results['needFasta'] = needFasta

    return JsonResponse(results, status=200)

@login_required
def parse_bed(request):
    start_time = time.time()
    gfa = request.POST.get('gfa','')
    bed = request.POST.get('bed','')
    data = {}
    warning = ''
    error = ''

    if not gfa or not bed:
        return JsonResponse({'error':True, 'msg': 'missing value'}, safe=False, status=400)

    gfa = get_work_dir(request, ['gfa'], gfa)
    bed = get_work_dir(request, ['bed'], bed)

    #graph = PanGraph(gfa, outdir=get_work_dir(request), parseRGFA=False)
    """
    graph = PanGraph(gfa, outdir=get_work_dir(request))
    results = graph.parseBedGff(bed)
    bed_contigs = results['bed_contigs']
    gfa_contigs = graph.inf['backbone']['contigs']
    common_contigs = [value for value in bed_contigs if value in gfa_contigs]

    if not common_contigs and not results['error']:
        results['error'] = 'The contig/chr names in the selected file are different from those in the rGFA file. Please check'

    if results['error']:
        return JsonResponse({'error':results['error']}, status=200)

    graph.loadRGFA()
    graph.genGraph()
    #graph.updateNodes()
    data = graph.nodeGeneOverlap(bed, 2)

    results = {}
    results['status'] = 0
    results['gene'] = sorted(list(data.keys()))
    results['gene_info'] = data
    results['time'] = time.time() - start_time
    results['warning'] = warning
    results['error'] = error
    """

    graph = PanGraph(gfa, outdir=get_work_dir(request), parseRGFA=False)
    graph.loadBedGff(bed)
    chromDict = {graph.bed[geneId]['Chr']:1  for geneId in graph.bed }
    graph.loadLenBeforeDict(chromDict)

    data = graph.nodeGeneOverlap(overlapNodeCountThreshold=2)

    results = {}
    results['status'] = 0
    results['gene'] = sorted(list(data.keys()))
    results['gene_info'] = data
    results['time'] = time.time() - start_time
    results['warning'] = warning
    results['error'] = error

    return JsonResponse(results, status=200)

@login_required
def graph(request):
    template_name = 'pangraphviewer/graph.html'

    gfa = request.GET.get('gfa','')
    backbone = request.GET.get('backbone','')
    chr = request.GET.get('chr','')
    start = request.GET.get('start','')
    end = request.GET.get('end','')
    plot = request.GET.get('plot','')

    context = {
        'graph_page': 'active',
        'gfa': gfa,
        'backbone': backbone,
        'chr': chr,
        'start': start,
        'end': end,
        'plot': plot,
        'upload_file_ext': {type:','.join(upload_file_ext[type]) for type in upload_file_ext}
    }

    return render(request, template_name, context=context)

@transaction.atomic
@login_required
def config(request):
    global work_base_dir

    template_name = 'pangraphviewer/config.html'

    if request.is_ajax():
        try:
            with transaction.atomic():
                work_base_dir = request.POST.get('work_base_dir','')
                profile = Profile.objects.get_or_create(username=request.user)
                profile[0].work_base_dir = work_base_dir
                profile[0].save()

                work_base_dir = None
                get_work_dir(request)

                status = 200
        except:
            status = 400

        return JsonResponse(None, safe=False, status=status)
    else:
        try:
            profile = Profile.objects.get_or_create(username=request.user)
            if not profile[0].work_base_dir:
                profile[0].work_base_dir = getVar(copied, 'web', 'work_dir', mustHave=True)
                profile[0].save()

            get_work_dir(request)
        except:
            profile[0].work_base_dir = getVar(copied, 'web', 'work_dir', mustHave=True)
            profile[0].save()

            work_base_dir = None
            get_work_dir(request)

        context = {
            'config_page': 'active',
            'work_base_dir': work_base_dir,
        }

        return render(request, template_name, context=context)

def getnodes(request):
    results = []

    gfa = request.POST.get('gfa','')
    seqname = request.POST.get('seqname','').replace(',',' ').split()

    gfa = get_work_dir(request, ['gfa'], gfa)
    panGraph = PanGraph(gfa, os.path.dirname(gfa), parseRGFA=False)
    nodes = panGraph.getNodeFromRGFA(seqname)

    for nodeId in nodes:
        node = nodes[nodeId]
        title = f"NodeId: {nodeId}; Resource: {node[NODE.sample]}_{node[NODE.chr]}; Len: {node[NODE.len]}"
        if node[NODE.rank] == '0':
            title += f"; Pos: {node[NODE.lenBefore]} - {node[NODE.lenBefore]+node[NODE.len]-1}"
        else:
            title += f"; SV-Info: {node[NODE.inf]['raw']}"
        results.append({'seqname':nodeId, 'title':title, 'seq':node[NODE.seq]})

    return results

@login_required
def viewseq(request):
    nodes = getnodes(request)
    text = '\n'.join([f">{node['seqname']}\t{node['title']}\n{node['seq']}" for node in nodes])

    mime_type = "text/plain"
    response = HttpResponse(content=text, content_type=mime_type)

    return response

@login_required
def downloadseq(request):
    nodes = getnodes(request)
    text = '\n'.join([f">{node['seqname']}\t{node['title']}\n{node['seq']}" for node in nodes])

    mime_type = "text/plain"
    filename = 'sequence.txt'

    response = HttpResponse(content=text, content_type=mime_type)
    response['Content-Disposition'] = "attachment; filename=%s" % filename

    return response

@login_required
def download_rgfa(request):
    rgfa = request.GET.get('rgfa','')
    rgfa = get_work_dir(request, ['gfa'], rgfa)

    status = 200 if os.path.isfile(rgfa) else 400
    if status == 200:
        mime_type = "text/plain"

        filename = rgfa
        f = open(filename, 'r')

        response = HttpResponse(content=f, content_type=mime_type)
        response['Content-Disposition'] = "attachment; filename=%s" % os.path.basename(rgfa)

        return response
    else:
        return JsonResponse('', safe=False, status=400)

@login_required
def convert_vcf(request):
    output = None

    vcf = request.POST.get('vcf','')
    fasta = request.POST.get('fasta','')
    backbone = request.POST.get('backbone','')
    chr = request.POST.get('chr','')

    if vcf: vcf = get_work_dir(request, ['vcf'], vcf)
    if fasta: fasta = get_work_dir(request, ['fasta'], fasta)

    if not vcf or not fasta or not backbone:
        output = {'error':True, 'msg': 'Missing value'}
    elif not os.path.isfile(vcf):
        output = {'error':True, 'msg': 'VCF file not exists'}
    elif not os.path.isfile(fasta):
        output = {'error':True, 'msg': 'Fasta file not exists'}

    if output:
        return JsonResponse(output, safe=False, status=400 if output['error'] else 200)

    out_dir = get_work_dir(request, ['gfa'])
    prefix = f'{backbone}_{os.path.splitext(os.path.basename(vcf))[0]}'
    if chr:
        prefix += f'_{chr}'
    convert_vcf_to_gfa(vcf, chr, backbone, out_dir, prefix, fasta)
    outfile = f'{prefix}.gfa'

    output = {'error':'','rgfa':os.path.basename(outfile)}
    status = 200 if os.path.isfile(get_work_dir(request, ['gfa'], output['rgfa'])) else 400

    return JsonResponse(output, safe=False, status=200)

@login_required
def vcf_to_gfa(request):
    template_name = 'pangraphviewer/vcf_to_gfa.html'
    context = {
        'utilities_page': 'active',
        'convert_vcf_page': 'active',
    }

    return render(request, template_name, context=context)

@login_required
def upload_file(request):
    if request.is_ajax():
        error, warning = '', ''
        work_dir = get_work_dir(request)
        upload = request.FILES['image']
        file_type = request.POST.get('file_type','')
        file_path = get_work_dir(request, [file_type], upload.name)
        file_path_orig = file_path
        with open(file_path,'wb') as f:
            for chunk in upload.chunks():
                f.write(chunk)
        os.chmod(file_path, 0o666)

        if file_type == 'gfa':
            retcode = GFA2rGFA(file_path, None).checkGfaFormat()
            if retcode == 2: # GFA1
                file = os.path.split(file_path)[1]
                converted = f'{os.path.splitext(file)[0]}_converted{os.path.splitext(file)[1]}'
                warning = f'The file is in GFA v1. Conversion has been taken, and the converted file {converted} will be used instead'
                converted = get_work_dir(request,['gfa'], converted)
                results = convert_gfa_to_rgfa(file_path, converted)
                returncode = results['returncode']
                if returncode == 4:
                    error = 'Missing GFA field in the input file. Abort!'
                elif returncode == 5:
                    error = 'Both seq and LN tag are missing in the input file. Abort!'
                elif returncode == 6:
                    error = 'The file is in GFA v1, but error occurs during conversion. Abort!'
                elif returncode == 7:
                    error = 'Path info is missing from the input GFA1 file, which is necessary for conversion. Abort!'
                else:
                    file_path = converted
            elif retcode == 3: # UNKNOWN
                error = 'The input file is in an unknown format. Abort!'

            if retcode != 1:
                try:
                    os.remove(file_path_orig)
                except OSError:
                    pass

        return JsonResponse({'filepath': file_path,'error':error,'warning':warning}, status=200)

    return JsonResponse({}, status=400)

def manage_file(request):
    action = request.POST.get('action')
    file_name = request.POST.get('file_name','')
    file_type = request.POST.get('file_type','')

    results, status = {}, 200

    if action == 'upload':
        work_dir = get_work_dir(request)
        upload = request.FILES['image']
        file_path = get_work_dir(request, [file_type], upload.name)
        with open(file_path,'wb') as f:
            for chunk in upload.chunks():
                f.write(chunk)
        os.chmod(file_path, 0o666)

        results, status = {'filepath': file_path}, 200
    elif action == 'download':
        file_path = get_work_dir(request, [file_type], file_name)
        status = 200 if os.path.isfile(file_path) else 400
        if status == 200:
            mime_type = "text/plain"
            f = open(file_path, 'r')
            response = HttpResponse(content=f, content_type=mime_type)
            response['Content-Disposition'] = "attachment; filename=%s" % os.path.basename(file_path)

            return response
        else:
            results, status = {'File cannot be downloaded'}, 400
    elif action == 'delete':
        file_path = get_work_dir(request, [file_type], file_name)
        results = {}
        try:
            os.remove(file_path)
            status = 200
        except:
            status = 400

    return JsonResponse(results, safe=False, status=status)

def draw_overlap_gene(request):
    gfa = request.GET.get('gfa','')
    bed = request.GET.get('bed','')
    gene_id = request.GET.get('gene_id','')

    if not gfa or not bed or not gene_id:
        return JsonResponse({'error':True, 'msg': 'missing value'}, safe=False, status=400)

    gfa = get_work_dir(request, ['gfa'], gfa)
    bed = get_work_dir(request, ['bed'], bed)

    graph = PanGraph(gfa, outdir=get_work_dir(request, ['gene']), parseRGFA=False)
    graph.loadBedGff(bed)
    graph.overlapGenes(gene_id)
    graph.drawOverlapGenes(gene_id)
    path = get_work_dir(request, ['gene'])
    file = os.path.join(path, f'drawOverlap_with_Gene-{gene_id}.html')
    test_file = open(file)
    response = HttpResponse(content=test_file)

    return response

def check_node_id(request):
    start_time = time.time()

    gfa = request.POST.get('gfa','')
    input_ids = request.POST.get('input_ids','')

    if not gfa or not input_ids:
        return JsonResponse({'error':True, 'msg': 'missing value'}, safe=False, status=400)

    gfa = get_work_dir(request, ['gfa'], gfa)
    pangraph = PanGraph(gfa, outdir=get_work_dir(request, ['gene']), parseRGFA=False)
    node_ids = [id.strip() for id in input_ids.split() if id.strip()]
    checkedNodeIds = pangraph.getNodeFromRGFA(node_ids)

    missing = list(set(node_ids) - set(list(checkedNodeIds.keys())))
    if missing:
        return JsonResponse({'error':True, 'msg': f'Some input node ids cannot be found: {missing[0]} ...'}, safe=False, status=400)

    # keep the order of node IDs in input_ids
    outputNodeIds = [nodeId for nodeId in node_ids if nodeId in checkedNodeIds]

    results = {}
    results['status'] = 0
    results['checked_node_ids'] = '\n'.join(outputNodeIds)
    #results['time'] = time.time() - start_time

    return JsonResponse(results, status=200)
