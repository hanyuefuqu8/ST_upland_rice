#/hwfssz1/ST_EARTH/Reference/ST_AGRIC/USER/zhongliyuan/miniconda3/envs/R4.1/bin/Rscript monocle_position.R

#sh GO.sh


for i in $(ls -d *_rice/)
do
cd $i
Rscript ../module_plot.R
Rscript ../gene_plot.R
cd ..
done
