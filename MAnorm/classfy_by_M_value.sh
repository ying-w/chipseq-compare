

#!/bin/bash
if [ $# -ne 3 ]
then
  echo "Usage: `basename $0` fold-change_cutoff_unbiased fold-change_cutoff_biased -log10p_cutoff_biased"
  exit
fi
# sed '1d' MAnorm_result.xls | awk 'BEGIN {OFS="\t"}{print $1,$2,$3>"MAnorm_result.tmp"}'
awk -v var1=$2 -v var2=$3 'BEGIN {OFS="\t"} {if($4 == "unique_peak1" && $7 >= var1 && $9 > var2)print $1,$2,$3>"sample1_uniq_peaks.bed"}'  MAnorm_result.xls
awk -v var1=$2 -v var2=$3 'BEGIN {OFS="\t"} {if($4 == "unique_peak2" && $7 <= -var1 && $9 > var2)print $1,$2,$3>"sample2_uniq_peaks.bed"}'  MAnorm_result.xls
awk -v var=$1 'BEGIN {OFS="\t"} {if($4 ~ /common_peak*/ && $7<var && $7>-var)print $1,$2,$3>"unbiased_peaks.tmp"}'  MAnorm_result.xls
mergeBed -i unbiased_peaks.tmp > unbiased_peaks.bed
rm unbiased_peaks.tmp
