#pos是紧随其后的位置发生插入缺失的坐标点
#下面输出GS106的位移坐标
#./.不做缺失值处理，当作ref
zcat all.head.vcf.gz|grep -v '#'|awk '
{split($5,alt,",");
i=substr($10,1,1);
if (i==0)
{changedl=0}
else if(i==".")
{changedl=0;
}
else{
changedl=length(alt[i])-length($4)
if(changedl > 0) pos=$2+length(alt[i])-1; else pos=$2+length($4)-1  
}
if (changedl!=0)
print $1,$2,pos,changedl
}' > GS106_position.txt
#去掉重叠区域并对位移位置进行累加
sed -n '1,/Above messages from 106/p' work.sh.e4585610 >GS106_mask.txt
awk '
NR==FNR{
split($3,m,":")
a[$3]=m[2]
}
NR!=FNR{
if ($2!=a[sprintf("%s:%s",$1,$2)]) print $0;
}' GS106_mask.txt GS106_position.txt|awk 'BEGIN{
a=0}
{
a=a+$4;
if (b!=$1) a=$4
print $0,a
b=$1
}' |sed 's/ /\t/g' > GS106_position.final.txt

#开始更改gtf坐标


awk -F'\t' 'BEGIN{
i=1
}NR==FNR{
a[i]=$1
b[i]=$3
c[i]=$4
d[i]=$5
i++
}
NR!=FNR{
i=1;
while(a[i]<$1 || (a[i]==$1 && b[i]<=$4))
{
i++
}
i--
if (a[i]==$1) 
{
if (b[i]==$4) start=$4+d[i-1]; else {
if (c[i]>0 || (b[i]-$4)<c[i]) start=$4+d[i]; else start=b[i]+d[i-1]+1 
}
} else {
start=$4
}
while(a[i]<$1 || (a[i]==$1 && b[i]<=$5))
{
i++
}
i--
if (a[i]==$1)
{
if (b[i]==$5) end=$5+d[i-1]; else{
if (c[i]>0 || (b[i]-$5)<c[i]) end=$5+d[i]; else end=b[i]+d[i-1]
}
} else {
end=$5
}
if (end-start>=0) print $1"\t"$2"\t"$3"\t"start"\t"end"\t"$6"\t"$7"\t"$8"\t"$9
}' GS106_position.final.txt Oryza_sativa.IRGSP-1.0.51.gtf |sed '/^\t/d' > GS69_test.gtf

