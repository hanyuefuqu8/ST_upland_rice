for species in $(awk 'NR!=1{print $4}' ../info2|sort |uniq)
do
/hwfssz1/ST_EARTH/Reference/ST_AGRIC/USER/zhongliyuan/miniconda3/envs/R4.1/bin/Rscript \
../monocle/work1.R \
-i ../diff_gene_SCT_position_adjusting/data_pos_adj.RDS \
-f ../info2 \
-l $species

mkdir $species
cd $species
/hwfssz1/ST_EARTH/Reference/ST_AGRIC/USER/zhongliyuan/miniconda3/envs/R4.1/bin/Rscript \
../transition_gene.R \
-i ../sub.RDS \
-f ../../info2 \
-l 0.05
cd ..
rm sub.RDS
done
