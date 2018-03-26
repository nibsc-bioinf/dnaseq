#!/bin/bash

project=$1
dir="/sequencing/projects/$project/"

cd $dir
echo `pwd`

bash dnaseq.bash.$project.index.sh
bash dnaseq.bash.$project.map.sh
bash dnaseq.bash.$project.report.sh
