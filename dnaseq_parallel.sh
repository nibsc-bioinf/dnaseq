#!/bin/bash

project=$1
jobs=$2
dir="/sequencing/projects/$project/"

cd $dir
echo `pwd`

parallel --jobs $jobs < dnaseq.bash.$project.index.sh
parallel --jobs $jobs < dnaseq.bash.$project.map.sh
parallel --jobs $jobs < dnaseq.bash.$project.report.sh
