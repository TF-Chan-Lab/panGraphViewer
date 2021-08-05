const alertBox = document.getElementById('alert-box')
const csrf = document.getElementsByName('csrfmiddlewaretoken')

var upload_file_url = document.getElementById('upload_file_url').value

upload_add_listener_change('vcf');
upload_add_listener_change('fasta');

function upload_add_listener_change(prefix) {
    upload = document.getElementById(prefix+'_upload')
    upload.onchange = function(e) {
        var prefix = e.target.id.split('_')[0]
        var file_data = e.target.files[0]

        var fd = new FormData()
        fd.append('csrfmiddlewaretoken', csrf[0].value)
        fd.append('image', file_data)
        fd.append('file_type', prefix)

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
                    }

                })
                return xhr
            },
            success: function(response){
                update_alert_box('Successfully uploaded the file below', 'alert-success')
                var idx = response['filepath'].lastIndexOf("/") + 1;
                var filename = response['filepath'].substr(idx);
                refresh_uploaded_list(`${prefix}_path`,prefix,filename);
            },
            error: function(error){
                console.log('error', error)
                update_alert_box('Ups... something went wrong', 'alert-danger')
            },
            cache: false,
            contentType: false,
            processData: false,
        });

        $(this).prop("value", "");
    }
}

$('#convert-btn').click(function() {
    event.preventDefault();

    var vcf_path = document.getElementById('vcf_path').value;
    var fasta_path = document.getElementById('fasta_path').value;
    var backbone = document.getElementById('backbone').value;
    var convert_vcf_url = document.getElementById('convert_vcf_url').value;

    var rgfa_path = document.getElementById('rgfa_path')

    if (!backbone) {
        update_alert_box('Please input backbone name', 'alert-danger')
        document.getElementById('backbone').focus();
        return;
    }

    if (!vcf_path) {
        update_alert_box('Please select uploaded VCF file, or upload new local VCF file', 'alert-danger')
        document.getElementById('vcf_path').focus();
        return;
    }

    if (!fasta_path) {
        update_alert_box('Please select uploaded fasta file, or upload new local fasta file', 'alert-danger')
        document.getElementById('fasta_path').focus();
        return;
    }

    var post_data = {
        'csrfmiddlewaretoken': csrf[0].value,
        'vcf': vcf_path,
        'fasta': fasta_path,
        'backbone': backbone,
    }

    alertBox.innerHTML = `<div class="alert alert-info" role="alert">Converting VCF ...</div>`
    $.ajax({
        type:'POST',
        url: convert_vcf_url,
        data: post_data,
        success: function(result) {
            rgfa_path.value = result['rgfa'];
            $('#download-btn').prop('disabled', false);
            alertBox.innerHTML = `<div class="alert alert-success" role="alert">Conversion done</div>`
        },
        error: function(result) {
            rgfa_path.value = '';
            alertBox.innerHTML = `<div class="alert alert-danger" role="alert">Conversion failed: ${result.responseJSON.msg}</div>`
        }
    });
});

$('#download-btn').click(function() {
    event.preventDefault();

    download_rgfa();
});

function download_rgfa() {
    var rgfa = document.getElementById('rgfa_path').value;
    var url = 'download_rgfa?rgfa=' + rgfa;
    window.open(url);
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
            dropdown_obj = $(`#${dropdown_id}`).find('option:not(:first)').remove();
            $.each(result.files, function(index, value) {
                if (value == selected_value) selected = 'selected'; else selected = '';
                dropdown_obj.end().append(`<option value="${value}" ${selected}>${value}</option>`);
            });
            update_alert_box(`${file_type.toUpperCase()} uploaded file list updated`, 'alert-success');
        },
        error: function(result) {
            update_alert_box(`Error in getting ${file_type.toUpperCase()} uploaded file list`, 'alert-danger');
        }
    });
}

function delete_uploaded(dropdown_id, file_type) {
    action = 'delete';
    file_name = $(`#${dropdown_id}`).val();
    if (!file_name) return;

    manage_file(action, file_name, file_type);
}

$.download_file = function(url, post_data){
    var form = $('<form></form>').attr('action', url).attr('method', 'post');
    $.each(post_data, function(name, value) {
        form.append($("<input></input>").attr('type', 'hidden').attr('name', name).attr('value', value));
    });
    form.appendTo('body').submit().remove();
};

function download_uploaded(dropdown_id, file_type) {
    action = 'download';
    file_name = $(`#${dropdown_id}`).val();
    if (!file_name) return;

    manage_file(action, file_name, file_type);
}

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
