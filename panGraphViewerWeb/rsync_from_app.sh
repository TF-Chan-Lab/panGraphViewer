TEMPLATE_PATH="${APP_PATH}/scripts/template"

# css
rsync -av ${TEMPLATE_PATH}/cytoscape-context-menus.css static/pangraphviewer/css/ > /dev/null

# js
rsync -av ${TEMPLATE_PATH}/cytoscape-context-menus.js static/pangraphviewer/js/ > /dev/null
rsync -av ${TEMPLATE_PATH}/cytoscape-euler.js static/pangraphviewer/js/ > /dev/null
rsync -av ${TEMPLATE_PATH}/cytoscape-qtip.js static/pangraphviewer/js/ > /dev/null

# misc
rsync -av ${TEMPLATE_PATH}/loader.gif static/images/ > /dev/null

echo rsync done
