#!/bin/bash
export PATH=/hwfssz1/ST_EARTH/Reference/ST_AGRIC/USER/zhongliyuan/miniconda3/bin/:$PATH

source activate st

for i in $(cat sections.txt)
do 
name1=${i##*/}
name2=${name1%%.*}
name3=${name2%%_*} 
name4=${name3%%-*}
mkdir $name2
python loom_generated.py $i $name4.gtf $name2
done
