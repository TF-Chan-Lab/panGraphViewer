const gfa_upload = document.getElementById('gfa_upload')

const gfa = document.getElementById('gfa_path')
const alertBox = document.getElementById('alert-box')
const plotBtn = document.getElementById('plot-btn')
const csrf = document.getElementsByName('csrfmiddlewaretoken')
const upload_file_url = document.getElementById('upload_file_url').value
const parse_gfa_url = document.getElementById('parse_gfa_url').value
const parse_vcf_url = document.getElementById('parse_vcf_url').value
const parse_bed_url = document.getElementById('parse_bed_url').value
const getdata_url = document.getElementById('getdata_url').value
var sample_list = '';
var gene_info;

function gen_graph(chr, start, end, tabHeader, by_node_id=false) {
    var url = getdata_url;
    var input_type = document.getElementById('input_type').value;
    var gfa = document.getElementById('gfa_path').value;
    var vcf = document.getElementById('vcf_path').value;
    var fasta = document.getElementById('fasta_path').value;
    var backbone = document.getElementById('backbone').value;
    if (!chr) chr = document.getElementById('chr').value;
    if (!start) start = document.getElementById('start').value;
    if (!end) end = document.getElementById('end').value;

    node_ids = null;
    if (by_node_id) {
        node_ids = $('#extract_node_checked_node_id').val();
    }

    if (input_type == 'vcf') {
        str = 'While a fasta file is missing, we can still generate a graph. However, node sequence checking will not be available. Do you want to continue ?';
        if (!fasta && !confirm(str)) {
            update_alert_box('Plotting by VCF is cancelled', 'alert-info');
            return;
        }
    }

    start_str = start;
    end_str = end;
    if (!start_str && !end_str) title = `${chr}`
    else {
        if (!start_str) start_str = '1';
        title = `${chr}: ${start_str} - ${end_str}`
    }

    if (!tabHeader) tabHeader = title;

    update_alert_box('Generating data ...', 'alert-info')

    data = {'csrfmiddlewaretoken': csrf[0].value,'input_type':input_type,'gfa':gfa,'vcf':vcf,'fasta':fasta,
            'backbone':backbone,'chr':chr,'start':start,'end':end,'sample_list':sample_list,'node_ids':node_ids}

    $.ajax({
        type:'POST',
        url: url,
        data: data,
        dataType: 'json',
        success: function(result) {
            if (result.error) {
                alert(result.error);
                update_alert_box(result.error, 'alert-danger')
                return;
            }
            if (result.warning) {
                alert(result.warning);
            }
            update_alert_box('Generation done', 'alert-success')

            cyData = result.cyData;
            legend = result.legend;
            result_gfa = result.gfa
            addOutputTab(title, cyData, result_gfa, legend, tabHeader);
        },
        error: function(result) {
            obj = result.responseJSON;
            str = 'Generation failed';
            if (obj && 'msg' in obj) str += ': '+ obj.msg;
            update_alert_box(str, 'alert-danger')
        }
    });
}

$( "#input-form" ).submit(function( event ) {
    event.preventDefault();

    gen_graph();
});

$('#parse-btn').click(function() {
    event.preventDefault();

    var gfa = document.getElementById('gfa_path').value;

    if (!gfa) {
        update_alert_box('Please select uploaded (r)GFA file, or upload new local (r)GFA file', 'alert-danger')
        document.getElementById('gfa_path').focus();
        return;
    }

    var post_data = {
        'csrfmiddlewaretoken': csrf[0].value,
        'gfa': gfa
    }

    update_alert_box('Parsing gfa file ...', 'alert-info')
    $.ajax({
        type:'POST',
        url: parse_gfa_url,
        data: post_data,
        success: function(result) {
            sample_list = result.backbone;

            backbone_obj = $('#backbone').find('option:not(:first)').remove();
            chr_obj = $('#chr').find('option:not(:first)').remove();
            $('#start').val('');
            $('#end').val('');

            if (result.warning) update_alert_box(result.warning, 'alert-info')
            else if (result.error) {
                $("#tab-action a").addClass('disabled')
                $("#tab-action a").removeClass('active')
                $("#tab-content-action :input").prop('disabled', true)

                update_alert_box(result.error, 'alert-danger')
                return
            } else {
                update_alert_box('Parsing done', 'alert-success')
            }

            // update backbone dropdown
            $.each(result.backbone, function(index, value) {
                backbone_obj.end().append('<option value="'+value+'">'+value+'</option>');
            });

            // update chr dropdown
            $.each(result.chr, function(index, value) {
                chr_obj.end().append('<option value="'+value+'">'+value+'</option>');
            });

            // enable tabs
            $("#tab-action a").removeClass('disabled')
            $("#plot_tab").trigger('click');
            // enable inputs in tab-content
            $("#tab-content-action :input").prop('disabled', false)
            $("#bed_path").trigger('change');

            // disable specific controls
            btn = ['extract_node_view_btn',
                   'extract_node_plot_btn',
                   'extract_node_download_btn'];
            $.each(btn, function(index, value) {
              $(`#${value}`).prop('disabled', true);
            });
        },
        error: function(result) {
            /*
            $("#backbone").prop('disabled', true);
            $("#chr").prop('disabled', true);
            $("#start").prop('disabled', true);
            $("#end").prop('disabled', true);
            $("#plot-btn").prop('disabled', true);
            */
            control_plot_input('gfa', 'reset');

            update_alert_box('Parsing failed', 'alert-danger')
        }
    });
});

$('#plot-gene-btn').click(function() {
    event.preventDefault();

    var gfa = document.getElementById('gfa_path').value;
    var bed = document.getElementById('bed_path').value;
    var gene = document.getElementById('gene').value;

    if (!gfa) {
        update_alert_box('Please select uploaded (r)GFA file, or upload new local (r)GFA file', 'alert-danger');
        document.getElementById('gfa_path').focus();
        return;
    }

    if (!bed) {
        update_alert_box('Please select uploaded BED/GTF/GFF file, or upload new local BED/GTF/GFF file', 'alert-danger');
        document.getElementById('bed_path').focus();
        return;
    }

    if (!gene) {
        update_alert_box('Please select gene', 'alert-danger');
        document.getElementById('gene').focus();
        return;
    }

    plot_gene(gfa, bed, gene)
});

$('#parse-vcf-btn').click(function() {
    event.preventDefault();

    var vcf = document.getElementById('vcf_path').value;
    var vcf_backbone = document.getElementById('vcf_backbone').value;

    if (!vcf_backbone) {
        update_alert_box('Please input vcf backbone name', 'alert-danger')
        document.getElementById('vcf_backbone').focus();
        return;
    }

    if (!vcf) {
        update_alert_box('Please select uploaded VCF file, or upload new local VCF file', 'alert-danger');
        document.getElementById('vcf_path').focus();
        return;
    }

    var post_data = {
        'csrfmiddlewaretoken': csrf[0].value,
        'vcf': vcf,
        'vcf_backbone': vcf_backbone,
    }

    update_alert_box('Parsing vcf file ...', 'alert-info')
    $.ajax({
        type:'POST',
        url: parse_vcf_url,
        data: post_data,
        success: function(result) {
            sample_list = result.backbone;

            backbone_obj = $('#backbone').find('option:not(:first)').remove();
            chr_obj = $('#chr').find('option:not(:first)').remove();
            $('#start').val('');
            $('#end').val('');

            if (result.error) {
                alert(result.error);
                update_alert_box(result.error, 'alert-danger')
                control_plot_input('vcf', 'reset');
                return;
            } else if (result.warning) {
                alert(result.warning);
                update_alert_box(result.warning, 'alert-info')
            } else {
                update_alert_box('Parsing done', 'alert-success')
            }

            if (result.needFasta == 1){
                $('#fasta_path').attr('required', true);
                if (!document.getElementById('fasta_path').value) {
                   alert('As VCF contains no contig length, fasta file is needed to do conversion');
                }
            } else {
                $('#fasta_path').attr('required', false);
            }

            // update backbone dropdown
            $.each(result.backbone, function(index, value) {
                backbone_obj.end().append('<option value="'+value+'">'+value+'</option>');
            });

            // update chr dropdown
            $.each(result.chr, function(index, value) {
                chr_obj.end().append('<option value="'+value+'">'+value+'</option>');
            });

            // enable tabs
            $("#plot_tab").removeClass('disabled')
            $("#plot_tab").trigger('click');
            // enable inputs in tab-content
            $("#plot :input").prop('disabled', false)
        },
        error: function(result) {
            /*
            $("#backbone").prop('disabled', true);
            $("#chr").prop('disabled', true);
            $("#start").prop('disabled', true);
            $("#end").prop('disabled', true);
            $("#plot-btn").prop('disabled', true);
            */
            control_plot_input('vcf', 'reset');

            update_alert_box('Parsing failed', 'alert-danger')
        }
    });
});

$('#parse-bed-btn').click(function() {
    event.preventDefault();

    var gfa = document.getElementById('gfa_path').value;
    var bed = document.getElementById('bed_path').value;

    if (!gfa) {
        update_alert_box('Please select uploaded (r)GFA file, or upload new local (r)GFA file', 'alert-danger');
        document.getElementById('gfa_path').focus();
        return;
    }

    if (!bed) {
        update_alert_box('Please select uploaded BED file, or upload new local BED file', 'alert-danger');
        document.getElementById('bed_path').focus();
        return;
    }

    var post_data = {
        'csrfmiddlewaretoken': csrf[0].value,
        'gfa': gfa,
        'bed': bed,
    }

    update_alert_box('Parsing bed file ...', 'alert-info')
    $.ajax({
        type:'POST',
        url: parse_bed_url,
        data: post_data,
        success: function(result) {
            gene_obj = $('#gene').find('option:not(:first)').remove();
            gene_list = result.gene;
            gene_info = result.gene_info;

            if (result.error) {
                update_alert_box(result.error, 'alert-danger');
                $("#plot-gene-btn").prop('disabled', true);
                return
            } else if (result.warning) {
                update_alert_box(result.warning, 'alert-info');
            } else {
                update_alert_box('Parsing done', 'alert-success')
            }

            // update gene dropdown
            $.each(result.gene, function(index, value) {
                gene_obj.end().append('<option value="'+value+'">'+value+'</option>');
            });
            $("#gene").prop('disabled', false);
            //$("#plot-gene-btn").prop('disabled', false);
        },
        error: function(result) {
            $("#gene").prop('disabled', true);
            update_alert_box('Parsing failed', 'alert-danger')
        }
    });
});

function bed_path_onchange() {
    $('#gene').find('option:not(:first)').remove();
    $("#plot-gene-btn").prop('disabled', true);
    if ($('#bed_path').val()) {
        $("#parse-bed-btn").prop('disabled', false);
    } else {
        $("#parse-bed-btn").prop('disabled', true);
    }
}

function gene_onchange() {
    if ($('#gene').val()) {
        $("#plot-gene-btn").prop('disabled', false);
    } else {
        $("#plot-gene-btn").prop('disabled', true);
    }
}

function download_sequence(ids, gfa_path) {
    var form = $('<form></form>').attr('action', 'downloadseq').attr('method', 'post');
    form.append($("<input></input>").attr('type', 'hidden').attr('name', 'csrfmiddlewaretoken').attr('value', csrf[0].value));
    form.append($("<input></input>").attr('type', 'hidden').attr('name', 'gfa').attr('value', gfa_path));
    form.append($("<input></input>").attr('type', 'hidden').attr('name', 'seqname').attr('value', ids));
    form.appendTo('body').submit().remove();
}

function view_sequence(ids, gfa_path) {
    var form = $('<form></form>').attr('action', 'viewseq').attr('method', 'post').attr('target','_blank');
    form.append($("<input></input>").attr('type', 'hidden').attr('name', 'csrfmiddlewaretoken').attr('value', csrf[0].value));
    form.append($("<input></input>").attr('type', 'hidden').attr('name', 'gfa').attr('value', gfa_path));
    form.append($("<input></input>").attr('type', 'hidden').attr('name', 'seqname').attr('value', ids));
    form.appendTo('body').submit().remove();
}

function drawGraph2(loadStatusId, cyId, cyData) {

    function makePopper(ele) {
        let ref = ele.popperRef();

        ele.tippy = tippy(ref, {
            content: () => {
                let content = document.createElement('div');

                content.innerHTML = ele.data('title');

                return content;
            },
            trigger: 'manual'
        });
    }

    var cy = window.cy = cytoscape({
        container: document.getElementById(cyId),

        layout: {
            //name: 'euler',
            name: 'fcose',
            randomize: true,
            animate: false,

            stop: function() {
                $('#' + loadStatusId).val(1);
                $('#loading').hide();
            },
        },

        style: [
            {
                selector: 'node',
                style: {
                    'content': 'data(name)',
                    'background-color': 'data(color)',
                    'shape': 'data(shape)',
                    'width': 'data(size)',
                    'height': 'data(size)',
                    'border-width': '0.09px',
                    'border-style': 'solid',
                    'border-color': '#E5E5E5',
                    'font-size':'8px'
                }
            },

            {
                selector: 'edge',
                style: {
                    'source-label': 'data(sourceLabel)',
                    'target-label': 'data(targetLabel)',
                    'width': 'data(weight)',
                    'opacity': 0.5,
                    'line-color':'data(color)',
                    'arrow-scale': 1.1,
                    'curve-style': 'unbundled-bezier',
                    'control-point-distance': 7,
                    'control-point-weight': '0.5',
                    'target-arrow-shape':'data(arrow)',
                    'target-arrow-color':'data(color)',
                    'font-size':'20px',
                    'color':'red',
                    'source-text-offset':'7px',
                    'target-text-offset':'7px',
                }
            },

            {
                selector: 'node:selected',
                style: {
                    'border-color': 'red',
                    'border-width': '2px',
                }
            },

            {
                selector: 'edge:selected',
                style: {
                    'line-color': 'red',
                    'target-arrow-color':'red',
                    'width': '2px',
                }
            },
        ],

        elements: cyData
    });

    cy.contextMenus({
      menuItems: [
        {
          id: 'download',
          content: 'download sequence',
          tooltipText: 'download sequence',
          selector: '*',
          onClickFunction: function (event) {
              curr_id = (event.target || event.cyTarget).data('id');
              selected = cy.$(':selected').jsons();
              let ids = selected.map(a => a.data.id);
              if ($.inArray(curr_id, ids)<0) ids.push(curr_id);
              gfaPath = $(`#${cyId.replace('cy','gfa_path')}`).val();
              download_sequence(ids, gfaPath);
          },
          disabled: false
        },
        {
          id: 'view',
          content: 'view sequence',
          tooltipText: 'view sequence',
          selector: '*',
          onClickFunction: function (event) {
              curr_id = (event.target || event.cyTarget).data('id');
              selected = cy.$(':selected').jsons();
              let ids = selected.map(a => a.data.id);
              if ($.inArray(curr_id, ids)<0) ids.push(curr_id);
              gfaPath = $(`#${cyId.replace('cy','gfa_path')}`).val();
              view_sequence(ids, gfaPath);
          },
          disabled: false,
          hasTrailingDivider: true
        },
        {
          id: 'hide_selected',
          content: 'hide selected',
          tooltipText: 'hide selected',
          selector: '*',
          onClickFunction: function (event) {
              cy.$(':selected').hide();
              curr_node = (event.target || event.cyTarget);
              curr_node.hide();
          },
          disabled: false
        },
        {
          id: 'show_selected_only',
          content: 'show selected (only)',
          tooltipText: 'show unselected (only)',
          selector: '*',
          onClickFunction: function (event) {
              cy.$(':unselected').hide();
              curr_node = (event.target || event.cyTarget);
              curr_node.show();
          },
          disabled: false
        },
        {
          id: 'show_all',
          content: 'show all',
          tooltipText: 'show all',
          selector: '*',
          onClickFunction: function (event) {
              cy.$(':unselected').show();
              cy.$(':selected').show();
          },
          disabled: false
        },
      ]
    });

    cy.ready(function() {
        cy.elements().forEach(function(ele) {
            makePopper(ele);
        });
    });

    cy.nodes().unbind('mouseover');
    cy.nodes().bind('mouseover', (event) => event.target.tippy.show());

    cy.nodes().unbind('mouseout');
    cy.nodes().bind('mouseout', (event) => event.target.tippy.hide());

    var minX=null, minY=null, maxX=null, maxY=null;
    cy.nodes().forEach(function(ele) {
        var x = ele.position().x;
        var y = ele.position().y;

        if (!minX || x < minX) minX = x;
        if (!maxX || x > maxX) maxX = x;
        if (!minY || y < minY) minY = y;
        if (!maxY || y > maxY) maxY = y;
    });
    minX -= 50; minY -= 50;
    maxX += 50; maxY += 50;
    edge_weight = (maxX - minX)/1000;
    color = '#777777';
    shape = 'ellipse';
    weight = 75;
    cy.add([{group:'nodes', data:{ id:'nodeA', name:'', weight:weight, shape:shape, size:1, color:color}, position:{x:minX, y:minY}},
            {group:'nodes', data:{ id:'nodeB', name:'', weight:weight, shape:shape, size:1, color:color}, position:{x:minX, y:maxY}},
            {group:'nodes', data:{ id:'nodeC', name:'', weight:weight, shape:shape, size:1, color:color}, position:{x:maxX, y:maxY}},
            {group:'nodes', data:{ id:'nodeD', name:'', weight:weight, shape:shape, size:1, color:color}, position:{x:maxX, y:minY}},

            {group:'edges', data:{ id:'border1', weight:edge_weight, source:'nodeA', target:'nodeB', color:color, arrow:'none', sourceLabel:'', targetLabel:''}},
            {group:'edges', data:{ id:'border2', weight:edge_weight, source:'nodeB', target:'nodeC', color:color, arrow:'none', sourceLabel:'', targetLabel:''}},
            {group:'edges', data:{ id:'border3', weight:edge_weight, source:'nodeC', target:'nodeD', color:color, arrow:'none', sourceLabel:'', targetLabel:''}},
            {group:'edges', data:{ id:'border4', weight:edge_weight, source:'nodeD', target:'nodeA', color:color, arrow:'none', sourceLabel:'', targetLabel:''}},
    ]);
}

var gId=1;

function getId() {
    return gId++;
}

function addOutputTab(title, cyData, gfa, legend, tabHeader) {
    id = getId();
    tabId = 'tab' + id;
    tabContentId = 'tabContent' + id;
    closeBtnId = 'closeBtnId' + id;
    cyId = 'cy' + id;
    loadStatusId = 'loadStatus' + id;
    gfaPathId = 'gfa_path' + id;
    captureBtnId = 'captureBtnId' + id;
    reposBtnId = 'reposBtnId' + id;

    if (!tabHeader) tabHeader = title;

    closeBtnStr = `<button type="button" class="close" aria-label="Close" id="${closeBtnId}"><span aria-hidden="true">&times;</span></button>`
    nodeStr = `<li class="nav-item"><a class="nav-link" id="${tabId}" data-toggle="tab" href="#${tabContentId}" role="tab" aria-controls="home" aria-selected="true">${tabHeader}&nbsp;${closeBtnStr}</a></li>`
    $('#myTab').children(':first').after(nodeStr);
    //$('#myTab').append(nodeStr);

    $(`#${closeBtnId}`).on('click', function() {
        // remove tab
        $($(this).parents('a').attr('href')).remove();
        // remove anchor
        $(this).parents('li').remove('li');
        // show the 1st tab
        $('#myTab a:first').tab('show');
    });

    reposBtnStr = `<button type="button" class="btn btn-default repos_btn" id="${reposBtnId}" onclick="repos_cy(${id})"><i class="fa fa-sync-alt" aria-hidden="true"></i></button>`
    captureBtnStr = `<button type="button" class="btn btn-default capture_btn" id="${captureBtnId}" onclick="capture_cy(${id})"><i class="fas fa-camera" aria-hidden="true"></i></button>`
    legendStr = `<div class="myDIV">Legend (mouseover here to show)</div>`
    legendStr += `<div class="hide"><table><tr><td style="display: table-cell;vertical-align: top"><div id='colorTable${id}'></div></td><td style="display: table-cell;vertical-align: top"><div id='shapeTable${id}'></div></td></tr></table></div>`
    body = `<center><h1>${title}${reposBtnStr}${captureBtnStr}</h1>${legendStr}</center><div class="cy" id="${cyId}"></div><input type=hidden id="${loadStatusId}" value="-1">`
    body += `<input type=hidden id="${gfaPathId}" value="${gfa}">`

    tab = `<div class="tab-pane fade show" id="${tabContentId}" role="tabpanel" aria-labelledby="${tabId}">${body}</div>`;
    $('#myTabContent').append(tab);

    $(`#loading`).show();
    $(`#${tabId}`).tab('show');
    add_legend(id, legend);

    $(`#${tabId}`).on('shown.bs.tab', function (e) {
       loadStatus = $('#' + loadStatusId).val();
       if (loadStatus == -1) {
          $('#loading').show();
          $('#' + loadStatusId).val(0);
          drawGraph2(loadStatusId, cyId, cyData);
       } else if (loadStatus == 0) {
          $('#loading').show();
       } else {
          $('#loading').hide();
       }
    })
}

function upload_add_listener_change(prefix) {
    upload = document.getElementById(prefix+'_upload')
    upload.onchange = function(e) {
        var file_type = e.target.id.split('_')[0]
        var file_data = e.target.files[0]

        var prefix = file_type;
        var fd = new FormData()
        fd.append('csrfmiddlewaretoken', csrf[0].value)
        fd.append('image', file_data)
        fd.append('file_type', prefix)

        control_plot_input(file_type, 'reset');

        if (file_type == 'gfa') $('#parse-btn').prop('disabled', true);
        else if (file_type == 'vcf') $('#parse-vcf-btn').prop('disabled', true);

        $.ajax({
            type:'POST',
            url: upload_file_url,
            enctype: 'multipart/form-data',
            data: fd,
            beforeSend: function(){
                update_alert_box('')
            },
            xhr: function(){
                const xhr = new window.XMLHttpRequest();
                xhr.upload.addEventListener('progress', e=>{
                    if (e.lengthComputable) {
                        const percent = e.loaded / e.total * 100
                        update_alert_box(`Upload status: ${percent.toFixed(1)}%`, 'alert-info')
                        if (percent == 100) {
                            update_alert_box(`Upload done. Processing ...`, 'alert-info')
                        }
                    }
                })
                return xhr
            },
            success: function(response){
                if (response.error) {
                    alert(response.error);
                    update_alert_box(response.error, 'alert-danger');

                    if (file_type == 'gfa') control_plot_input('gfa', 'reset');
                    else if (file_type == 'vcf') control_plot_input('vcf', 'reset');

                    return;
                } else if (response.warning) {
                    alert(response.warning);
                    update_alert_box(response.warning, 'alert-info');
                } else if (!response.filepath) {
                    update_alert_box('Ups... something went wrong', 'alert-danger');
                    return;
                }

                update_alert_box('Successfully uploaded the file below', 'alert-success')
                var idx = response['filepath'].lastIndexOf("/") + 1;
                var filename = response['filepath'].substr(idx);
                refresh_uploaded_list(`${prefix}_path`,prefix,filename);

                if (file_type == 'gfa') $('#parse-btn').prop('disabled', false);
                else if (file_type == 'vcf') $('#parse-vcf-btn').prop('disabled', false);
            },
            error: function(error){
                //console.log('error', error)
                if (file_type == 'gfa') control_plot_input('gfa', 'reset');
                else if (file_type == 'vcf') control_plot_input('vcf', 'reset');

                update_alert_box('Ups... something went wrong', 'alert-danger')
            },
            cache: false,
            contentType: false,
            processData: false,
        });

        $(this).prop("value", "");
    }
}

function update_alert_box(msg, alert_type) {
    const alertBox = document.getElementById('alert-box')
    if (msg) {
        alertBox.innerHTML = `<div class="alert ${alert_type}" role="alert">${msg}</div>`
    } else {
        alertBox.innerHTML = '';
    }
}

function refresh_uploaded_list(dropdown_id, file_type, selected_value) {
    get_uploaded_list_url = document.getElementById('get_uploaded_list_url').value

    if (!selected_value) selected_value = $(`#${dropdown_id}`).val();

    var post_data = {
        'csrfmiddlewaretoken': csrf[0].value,
        'file_type': file_type
    }

    $.ajax({
        type:'POST',
        url: get_uploaded_list_url,
        data: post_data,
        success: function(result) {
            var msg;
            dropdown_obj = $(`#${dropdown_id}`).find('option:not(:first)').remove();
            $.each(result.files, function(index, value) {
                if (value == selected_value) selected = 'selected'; else selected = '';
                dropdown_obj.end().append(`<option value="${value}" ${selected}>${value}</option>`);
            });
            $(`#${dropdown_id}`).trigger('change');

            if (file_type == 'bed') msg = `BED/GTF/GFF uploaded file list updated`;
            else msg = `${file_type.toUpperCase()} uploaded file list updated`;

            update_alert_box(msg, 'alert-success');
        },
        error: function(result) {
            update_alert_box(`Error in getting ${file_type.toUpperCase()} uploaded file list`, 'alert-danger');
        }
    });
}

function download_uploaded(dropdown_id, file_type) {
    action = 'download';
    file_name = $(`#${dropdown_id}`).val();
    if (!file_name) return;

    manage_file(action, file_name, file_type);
}

function delete_uploaded(dropdown_id, file_type) {
    action = 'delete';
    file_name = $(`#${dropdown_id}`).val();
    if (!file_name) return;

    manage_file(action, file_name, file_type);
}

function control_plot_input(file_type, action) {
    if (action=='reset') {
        //dropdownlist_ids = ['backbone','chr']
        dropdownlist_ids = ['backbone','chr','bed_path','gene']
        $.each(dropdownlist_ids, function(index, value) {
            $(`#${value}`).find('option:not(:first)').remove();
        });

        //input_ids = ['start','end']
        input_ids = ['start','end','extract_node_node_id','extract_node_checked_node_id']
        $.each(input_ids, function(index, value) {
            $(`#${value}`).val('');
        });

        //input_ids = ['backbone','chr','start','end','plot-btn']
        input_ids = ['backbone','chr','start','end','plot-btn',
                     'extract_node_node_id','extract_node_checked_node_id','extract_node_view_btn','extract_node_plot_btn','extract_node_download_btn',
                     'bed_path','gene','parse-bed-btn','plot-gene-btn']
        $.each(input_ids, function(index, value) {
            $(`#${value}`).prop('disabled', true);
        });

        // disable tabs
        $("#tab-action a").addClass('disabled')
        $("#tab-action a").removeClass('active')
    }
    val = $(`#${file_type}_path`).val()
    update_related_control(file_type,val);
}

$('#rgfa_tab').on('shown.bs.tab', function (e) {
    $('#gfa_path').attr('required', true);
    $('#vcf_path').attr('required', false);
    $('#vcf_backbone').attr('required', false);

    $('#input_type').val('rgfa');

    // disable tabs
    $("#tab-action a").addClass('disabled')
    $("#tab-action a").removeClass('active')
    // disable inputs in tab-content
    $("#tab-content-action :input").prop('disabled', true)
});

$('#vcf_tab').on('shown.bs.tab', function (e) {
    $('#gfa_path').attr('required', false);
    $('#vcf_path').attr('required', true);
    $('#vcf_backbone').attr('required', true);

    $('#input_type').val('vcf');

    // disable tabs
    $("#tab-action a").addClass('disabled')
    $("#tab-action a").removeClass('active')
    // disable inputs in tab-content
    $("#tab-content-action :input").prop('disabled', true)
});

$.download_file = function(url, post_data){
    var form = $('<form></form>').attr('action', url).attr('method', 'post');
    $.each(post_data, function(name, value) {
        form.append($("<input></input>").attr('type', 'hidden').attr('name', name).attr('value', value));
    });
    form.appendTo('body').submit().remove();
};

function manage_file(action, file_name, file_type, select_id) {
    url = $('#manage_file_url').val()

    var post_data = {
        'csrfmiddlewaretoken': csrf[0].value,
        'action': action,
        'file_name': file_name,
        'file_type': file_type,
    }

    if (action == 'download') {
        $.download_file(url, post_data);
    } else {
        $.ajax({
            type:'POST',
            url: url,
            data: post_data,
            success: function(result) {
                //update_alert_box(`File removed`, 'alert-success');
                if (action == 'delete') alert(`${file_name} removed`);
                refresh_uploaded_list(`${file_type}_path`,file_type);
            },
            error: function(result) {
                update_alert_box(`Error on working on ${file_name}`, 'alert-danger');
            }
        });
    }
}

function plot_gene(gfa, bed, gene_id) {
    chr = gene_info[gene_id]['gene_chr'];
    start = gene_info[gene_id]['gene_start'];
    end = gene_info[gene_id]['gene_end'];
    gen_graph(chr, start, end, `Subgraph within the region of Gene: ${gene_id}`);
    url = `draw_overlap_gene?gfa=${gfa}&bed=${bed}&gene_id=${gene_id}`
    window.open(url);
}

$('#extract_node_view_btn').click(function () {
    var gfa = $('#gfa_path').val();
    var input_ids = $('#extract_node_node_id').val();
    var checked_ids = $('#extract_node_checked_node_id').val();

    if (!gfa) {
        update_alert_box('Please select uploaded (r)GFA file, or upload new local (r)GFA file', 'alert-danger');
        document.getElementById('gfa_path').focus();
        return;
    }

    if (!checked_ids) {
        if (input_ids) {
            update_alert_box('Checked node IDs empty. Node IDs not found, or missing press Check button to get checked node IDs', 'alert-danger');
            document.getElementById('extract_node_check_btn').focus();
        } else {
            update_alert_box('Please input node ids', 'alert-danger')
            document.getElementById('extract_node_node_id').focus();
        }
        return;
    }

    view_sequence(checked_ids, gfa);
});

$('#extract_node_plot_btn').click(function () {
    var gfa = $('#gfa_path').val();
    var input_ids = $('#extract_node_node_id').val();
    var checked_ids = $('#extract_node_checked_node_id').val();

    if (!gfa) {
        update_alert_box('Please select uploaded (r)GFA file, or upload new local (r)GFA file', 'alert-danger')
        document.getElementById('gfa_path').focus();
        return;
    }

    if (!checked_ids) {
        if (input_ids) {
            update_alert_box('Checked node IDs empty. Node IDs not found, or missing press Check button to get checked node IDs', 'alert-danger');
            document.getElementById('extract_node_check_btn').focus();
        } else {
            update_alert_box('Please input node ids', 'alert-danger')
            document.getElementById('extract_node_node_id').focus();
        }
        return;
    }

    gen_graph(null, null, null, 'Plot by node IDs', true);
});

$('#extract_node_download_btn').click(function () {
    var gfa = $('#gfa_path').val();
    var input_ids = $('#extract_node_node_id').val();
    var checked_ids = $('#extract_node_checked_node_id').val();

    if (!gfa) {
        update_alert_box('Please select uploaded (r)GFA file, or upload new local (r)GFA file', 'alert-danger')
        document.getElementById('gfa_path').focus();
        return;
    }

    if (!checked_ids) {
        if (input_ids) {
            update_alert_box('Checked node IDs empty. Node IDs not found, or missing press Check button to get checked node IDs', 'alert-danger');
            document.getElementById('extract_node_check_btn').focus();
        } else {
            update_alert_box('Please input node ids', 'alert-danger')
            document.getElementById('extract_node_node_id').focus();
        }
        return;
    }

    download_sequence(checked_ids, gfa);
});

function check_node_id(gfa, input_ids) {
    check_node_id_url = document.getElementById('check_node_id_url').value

    var post_data = {
        'csrfmiddlewaretoken': csrf[0].value,
        'gfa': gfa,
        'input_ids': input_ids,
    }

    $.ajax({
        type:'POST',
        url: check_node_id_url,
        data: post_data,
        success: function(result) {
            $('#extract_node_checked_node_id').val(result.checked_node_ids);

            btn = ['extract_node_view_btn',
                   'extract_node_plot_btn',
                   'extract_node_download_btn'];

            $.each(btn, function(index, value) {
              $(`#${value}`).prop('disabled', false);
            });

            update_alert_box(`Node ids checked`, 'alert-success');
        },
        error: function(result) {
            //update_alert_box(`Error in checking node IDs`, 'alert-danger');
            obj = result.responseJSON;
            str = 'Error in checking node ID';
            if (obj && 'msg' in obj) str += ': '+ obj.msg;
            update_alert_box(str, 'alert-danger');
        }
    });
}

$('#extract_node_check_btn').click(function () {
    var gfa = $('#gfa_path').val();
    var input_ids = $('#extract_node_node_id').val();

    if (!gfa) {
        update_alert_box('Please select uploaded (r)GFA file, or upload new local (r)GFA file', 'alert-danger')
        document.getElementById('gfa_path').focus();
        return;
    }

    if (!input_ids) {
        update_alert_box('Please input node ids', 'alert-danger')
        document.getElementById('extract_node_node_id').focus();
        return;
    }

    $('#extract_node_checked_node_id').val('');
    btn = ['extract_node_view_btn',
           'extract_node_plot_btn',
           'extract_node_download_btn'];
    $.each(btn, function(index, value) {
        $(`#${value}`).prop('disabled', true);
    });

    update_alert_box('Checking node IDs ...', 'alert-info')
    check_node_id(gfa, input_ids);
});

function update_related_control(update_type, update_value) {
    if (update_type == 'gfa' || update_type == 'vcf') {
        if (update_value == '') {
            // disable tabs
            $("#tab-action a").addClass('disabled')
            $("#tab-action a").removeClass('active')

            // disable inputs in tab-content
            $("#tab-content-action :input").prop('disabled', true)

            // disable buttons
            if (update_type == 'gfa') $('#parse-btn').prop('disabled', true);
            else if (update_type == 'vcf') $('#parse-vcf-btn').prop('disabled', true);
        }
        else {
            // enable buttons
            if (update_type == 'gfa') $('#parse-btn').prop('disabled', false);
            else if (update_type == 'vcf') $('#parse-vcf-btn').prop('disabled', false);
        }
    }
}

function capture_cy(id) {
    let div = document.getElementById(`cy${id}`);
    html2canvas(div).then(function (canvas) {
        var link = document.createElement('a');
        link.download = 'screencapture.png';
        link.href = canvas.toDataURL()
        link.click();
        alert('Screen is captured and downloaded');
    });
}

function repos_cy(id) {
    document.getElementById(`cy${id}`)._cyreg.cy.fit();
}

function add_legend(id, data) {
    var $table = $('<table/>');
    $.each(data['colors'], function(key, value) {
        $table.append( '<tr><td bgcolor="' + value + '">' + key + '</td></tr>' );
    });
    if (data['has_reversed']) {
        $table.append( '<tr><td style="display: table-cell">(*) nodes with seq<br>in rev-comp of orig</td></tr>' );
    }
    $('#colorTable' + id).append($table);

    var $table = $('<table/>');
    $.each(data['shapes'], function(key, value) {
        var url = '/static/images/cy/' + value + '.png';
        $table.append( '<tr><td><img src="' + url + '" width="20"></td><td>' + key + '</td></tr>' );
    });
    $('#shapeTable' + id).append($table);
}

upload_add_listener_change('gfa');
upload_add_listener_change('vcf');
upload_add_listener_change('fasta');
upload_add_listener_change('bed');
