#!/bin/bash
export PATH=/hwfssz1/ST_EARTH/Reference/ST_AGRIC/USER/zhongliyuan/miniconda3/bin/:$PATH
source activate velo

#cd /hwfssz1/ST_EARTH/Reference/ST_AGRIC/USER/zhongliyuan/miniconda3/envs/dynamo/lib/python3.8/site-packages/
#pip show dynamo-release
for i in $(cat sections.txt)
do
name1=${i##*/}
name2=${name1%%.*}
name3=${name2%%_*}
name4=${name3%%-*}
name5=${name2//_/-}
mkdir $name2
cd $name2
python \
../analysis2.py \
$name2 \
$name5 \
../../dynamo/
cd ..
done
