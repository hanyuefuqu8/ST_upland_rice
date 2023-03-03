awk -F'\t' 'NR==1{print "\t"$0}$4=="\"repeatable\""{print $0}' all_pr_deg.txt |sed 's/"//g' > all_sig_pr_deg.txt
awk -F'\t' 'NR==FNR{if(NR==1){head=$0}else{a[$1]=$0}}NR!=FNR{if(FNR==1){print $0"\t"head} else if (a[$2]) {print $0"\t"a[$2]}}' total_diffgene2.info.xls all_sig_pr_deg.txt > all_sig_diff.xls
awk -F'\t' 'NR==FNR{if(NR==1){head=$0}else{a[$2]=$0}}NR!=FNR{if(FNR==1){print $0"\t"head} else if (a[$2]) {print $0"\t"a[$2]}}' /hwfssz1/ST_EARTH/P20Z10200N0035/USER/zhongliyuan/zhongliyuan/makergenes/reference/IRGSP-1.0_representative_annotation_2022-03-11.tsv all_sig_pr_deg.txt > all_sig_info.xls

