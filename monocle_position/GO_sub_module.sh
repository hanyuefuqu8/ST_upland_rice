for i in $(ls -d *rice/)
do
cd $i
module=$(awk -F',' '{print $3}' sub_module7.csv|sort -n|uniq|tail -n 1)
smodule=$(awk -F',' '{print $3}' sub_module1.csv|sort -n|uniq|tail -n 1)
for((j=1;j<=$module;j++))
do
awk -F',' '$3=='"$j"'{print $2}' sub_module7.csv> m7s$j.list
sh /hwfssz1/ST_EARTH/Reference/ST_AGRIC/USER/zhongliyuan/GO/enrich_Os/OS_enrich.sh m7s$j.list
done

for((j=1;j<=$smodule;j++))
do
awk -F',' '$4=='"$j"'{print $2}' sub_module1.csv> m1s$j.list
sh /hwfssz1/ST_EARTH/Reference/ST_AGRIC/USER/zhongliyuan/GO/enrich_Os/OS_enrich.sh m1s$j.list
done

ls *.GO/*.difgo3.filt > file
/hwfssz1/ST_EARTH/Reference/ST_AGRIC/USER/zhongliyuan/miniconda3/envs/R4.1/bin/Rscript /hwfssz1/ST_EARTH/Reference/ST_AGRIC/USER/zhongliyuan/GO/plot.R
cd ..
done
