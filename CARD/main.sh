Rscript NIP_preprocess.R

for i in $(ls -d BIN40*/)
do
cd $i
:Rscript ../card.R
cd ..
done
