rm -f ./output/*.png
for file in ./data/*.csv; do
label=${file##*/}
Rscript visualise.r $file ${label%.csv}
done
