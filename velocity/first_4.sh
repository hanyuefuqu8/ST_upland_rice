#!/bin/bash
#the first parameter is annotation file,the second parameter is bam file, the third rameter is gem file, the fourth parameter is gem.gz in results
export PATH=/hwfssz1/ST_EARTH/Reference/ST_AGRIC/USER/zhongliyuan/miniconda3/bin/:$PATH

source activate velo

X=$(zcat SS200000130TR_A4.gem.gz |sed -n '5p'|awk -F'=' '{print $2}')
Y=$(zcat SS200000130TR_A4.gem.gz |sed -n '6p'|awk -F'=' '{print $2}')

samtools view -h SS200000130TR_A4_web_8.Aligned.sortedByCoord.out.merge.q10.dedup.target.bam |awk -v OFS="\t" '
BEGIN{a[0]="AA";
a[1]="AC";
a[2]="AG";
a[3]="AT";
a[4]="CA";
a[5]="CC";
a[6]="CG";
a[7]="CT";
a[8]="GA";
a[9]="GC";
a["A"]="GG";
a["B"]="GT";
a["C"]="TA";
a["D"]="TC";
a["E"]="TG";
a["F"]="TT";
}
{
if($1~/^@/)
{print $0}
else
{b=substr($18,6)
cx=substr($16,6)-"'$X'";cy=substr($17,6)-"'$Y'";d=""
for (i=1;i<=length(b);i++)
{
c[i]=a[substr(b,i,1)]
tnum="0x"substr(b,i,1)
#if (i==length(b) && strtonum(tnum)<4) c[i]=substr(c[i],2,1)
d=d c[i]
}
$18="UB:Z:"d;$16="Cx:i:"cx;$17="Cy:i:"cy
print $0
}
}' |sed 's/Cx:i:\([0-9]*\)\tCy:i:\([0-9]*\)/CB:Z:\1_\2/g'|samtools view -b -S >1.bam


awk '{print $2"_"$3}' 130TR_A4down01.gem|sort|uniq > position.txt

#samtools view -h $2|awk '{if($1~/^@/) {print $0} else {print $0"\tCB:Z:"substr($16,6)"_"substr($17,6)}}'|sed 's/UR:Z:/UB:Z:/g' |samtools view -b -S -@12 > 1.bam

samtools sort 1.bam -@ 12 -o 2.bam
#samtools view -@12 $2  |awk '{print $0"\tCB:Z:"substr($16,6)"_"substr($17,6)}'|sed 's/UR:Z:/UB:Z:/g' |samtools view -b -S -@12 > 1.bam

#awk 'NR==FNR{a[$0]=$0}NR!=FNR&&a[substr($16,6)]||/^@/{print $0}' position.txt test.sam >test2.sam

#/zfssz5/BC_PUB/Software/03.Soft_ALL/samtools-1.2/samtools view -b -S A2.sam > A2.bam

samtools sort -t CB -O BAM -@ 12 -o cellsorted_2.bam 2.bam

velocyto run -b position.txt -o ./ 2.bam CX.gff3

python /hwfssz1/ST_EARTH/Reference/ST_AGRIC/USER/zhongliyuan/Script/velo/anndata_generator.py *.loom

