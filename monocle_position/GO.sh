for i in $(ls -d */)
do
cd $i
module=$(awk -F',' '{print $3}' phase_module.csv|sort -n|uniq|tail -n 1)
smodule=$(awk -F',' '{print $4}' phase_module.csv|sort -n|uniq|tail -n 1)
for((j=1;j<=$module;j++))
do
awk -F',' '$3=='"$j"'{print $2}' phase_module.csv> m$j.list
sh /hwfssz1/ST_EARTH/Reference/ST_AGRIC/USER/zhongliyuan/GO/enrich_Os/OS_enrich.sh m$j.list
done

for((j=1;j<=$smodule;j++))
do
awk -F',' '$4=='"$j"'{print $2}' phase_module.csv> sm$j.list
sh /hwfssz1/ST_EARTH/Reference/ST_AGRIC/USER/zhongliyuan/GO/enrich_Os/OS_enrich.sh sm$j.list
done

ls *.GO/*.difgo3.filt > file
/hwfssz1/ST_EARTH/Reference/ST_AGRIC/USER/zhongliyuan/miniconda3/envs/R4.1/bin/Rscript /hwfssz1/ST_EARTH/Reference/ST_AGRIC/USER/zhongliyuan/GO/plot.R
cd ..
done
