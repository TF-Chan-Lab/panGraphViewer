import base64 
import json
from shutil import copyfile

from os import linesep

shapes = {'cy':["ellipse", "triangle", "round-triangle", 
        "rectangle", "round-rectangle", "bottom-round-rectangle",
        "cut-rectangle", "barrel", "rhomboid", "diamond", 
        "round-diamond", "pentagon", "round-pentagon", "hexagon", 
        "round-hexagon", "concave-hexagon", "heptagon", "round-heptagon",
        "octagon", "round-octagon", "star", "tag", "round-tag", "vee"],\
         'vis':['dot','triangle','database','triangleDown','text','star',
        'ellipse', 'square', 'box', 'diamond']}
js = '../images.js'
js_web = '../../../../panGraphViewerWeb/static/pangraphviewer/js/images.js'

data = {}

def gen_js():
    for type in shapes:
        data[type] = {}
        for shape in shapes[type]:
            image = open(f'{type}/{shape}.png', 'rb')
            image_read = image.read()
            #image_64_encode = base64.encodestring(image_read)
            image_64_encode = base64.encodebytes(image_read)
            str = image_64_encode.decode("utf-8").replace(linesep,"")

            data[type][shape] = f'data:image/png;base64,{str}'

    with open(js, 'w') as f:
        print(f'var images = {{}};', file=f)
        for type in shapes:
             print(f"images = ", json.dumps(data), file=f)

    # copy to web
    copyfile(js, js_web)

gen_js()
