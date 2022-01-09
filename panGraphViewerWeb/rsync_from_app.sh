APP_PATH="../panGraphViewerApp"
TEMPLATE_PATH="${APP_PATH}/scripts/template"

# scripts
rsync -av ${APP_PATH}/scripts/gfa2rGFA.py pangraphviewer/ > /dev/null

# icons
rsync -av --delete ${TEMPLATE_PATH}/images/cy static/images/ > /dev/null

# misc
rsync -av ${TEMPLATE_PATH}/loader.gif static/images/ > /dev/null

echo rsync done
