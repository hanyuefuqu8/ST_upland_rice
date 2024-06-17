#input gem and interested_gene.lst and binsize
for i in $(ls *.gem)
do
name=${i%%.*}
B=40
perl TransBIN.pl $i  -size $B > TransBIN$B$name.txt
perl -ane '{print "$F[0]_$F[1]\t$F[2]\t$F[3]\n"}' TransBIN$B$name.txt > MBIN$B$name.txt
perl trans.pl MBIN$B$name.txt > BIN$B$name.csv
done

ls *.csv>file

mkdir merge
cd merge
ls ../*.csv>file
/hwfssz1/ST_EARTH/Reference/ST_AGRIC/USER/zhongliyuan/miniconda3/envs/R-2/bin/Rscript /hwfssz1/ST_EARTH/Reference/ST_AGRIC/USER/zhongliyuan/Script/Seurat_merge.R -i file -o ./ -s Gmax -t root -d 20 -r 0.6 -m 0 -k 50
cd ..

mkdir harmony
cd harmony
ls ../*.csv>file
/hwfssz1/ST_EARTH/Reference/ST_AGRIC/USER/zhongliyuan/miniconda3/envs/R-2/bin/Rscript /hwfssz1/ST_EARTH/Reference/ST_AGRIC/USER/zhongliyuan/Script/Seurat_harmony.R -i file -o ./ -s Gmax -t root -d 20 -r 0.6 -m 0 -k 50
cd ..

mkdir integrate
cd integrate
ls ../*.csv>file
/hwfssz1/ST_EARTH/Reference/ST_AGRIC/USER/zhongliyuan/miniconda3/envs/R-2/bin/Rscript /hwfssz1/ST_EARTH/Reference/ST_AGRIC/USER/zhongliyuan/Script/Seurat_SCTintegrate.R -i file -o ./ -s Gmax -t root -d 20 -r 0.7 -m 0 -k 50
cd ..

#mkdir refrence
#cd refrence
#ls ../*.csv>file
#/hwfssz1/ST_EARTH/Reference/ST_AGRIC/USER/zhongliyuan/miniconda3/envs/R-2/bin/Rscript /hwfssz1/ST_EARTH/Reference/ST_AGRIC/USER/zhongliyuan/Script/Seurat_reference.R -i file -o ./ -s Gmax -t root -d 20 -r 1 -m 200 -k 50
#cd ..

