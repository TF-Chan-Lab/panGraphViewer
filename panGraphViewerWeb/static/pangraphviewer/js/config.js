$( "#input-form" ).submit(function( event ) {
    event.preventDefault();

    var work_base_dir = document.getElementById('work_base_dir').value;
    var csrf = document.getElementsByName('csrfmiddlewaretoken')[0].value;

    var alertBox = document.getElementById('alert-box');
    alertBox.innerHTML = `<div class="alert alert-info" role="alert">Updating ...</div>`

    data = {'csrfmiddlewaretoken': csrf,'work_base_dir':work_base_dir}
    $.ajax({
        type:'POST',
        //url: url,
        data: data,
        dataType: 'json',
        success: function(result) {
            alertBox.innerHTML = `<div class="alert alert-success" role="alert">Update done</div>`
        },
        error: function(response) {
            obj = response.responseJSON;
            str = 'Update failed';
            if (obj && 'msg' in obj) str += ': '+ obj.msg;
            alertBox.innerHTML = `<div class="alert alert-danger" role="alert">` + str + `</div>`
        }
    });
});

