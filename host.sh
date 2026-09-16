#!/bin/bash
# host ShinyCellPlus
source /home/ubuntu/miniforge3/etc/profile.d/mamba.sh
mamba activate epictope
Rscript -e 'shiny::runApp("./app.R", port=5041, host="0.0.0.0")'
