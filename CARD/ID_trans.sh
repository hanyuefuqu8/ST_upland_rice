gzip -d features.tsv.gz

awk -F'\t' 'NR==FNR{a[substr($2,1,14)]=$1}NR!=FNR{print $0"\t"a[$1]}' /hwfssz1/ST_EARTH/P20Z10200N0035/USER/zhongliyuan/zhongliyuan/reference/reference/RAP-MSU_2022-03-11.txt features.tsv > temp.tsv

mv features.tsv ori_fea.tsv

awk -F'\t' '{if($3==""||$3=="None") $3=$2; print $3"\t"$3}' temp.tsv > features.tsv
