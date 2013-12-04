# Author: Ying Wu daiyingw@gmail.com
# Last major update: December 2013
# This version is heavily modified from the original version
# License for this code should be the same as the original code but license is not given for original
# Latest version: https://github.com/ying-w/chipseq-compare/tree/master/MAnorm
# Original : http://bcb.dfci.harvard.edu/~gcyuan/MAnorm/R_tutorial.html
#
# Requirements: fairly modern version of R, bedtools, and linux/unix/mac machine
#
# This script should be run with MAnorm3.R in the same directory
# StepI Reads will be filtered (remove non-chr chromosomes) and shifted
# StepII Intersecting peaks will be isolated
# StepIII Reads over peaks will be counted (both over all peaks and intersecting peaks only)
# StepIV run MAnorm3.R to use edgeR to find differentially bound peaks
# Cleanup rm tmp_*
#
# Several figures will be created to show fold change versus binding of shared 
# peaks before and after normalization and also two tab deliminated files (ending
# in .xls) will be created from the output of edgeR. This output format was made 
# to match original MAnorm output
# 
# Usage: MAnorm3.sh foldername gr_peak.bed er_peak.bed gr_rep1.bed gr_rep2.bed er_rep1.bed er_rep2.bed 120 110 95 100
# gr_peak.bed: sample 1 significant peak list
# er_peak.bed: sample 2 significant peak list
# gr_rep1.bed: sample 1 raw reads rep1
# gr_rep2.bed: sample 1 raw reads rep2
# er_rep1.bed: sample 2 raw reads rep1
# er_rep2.bed: sample 2 raw reads rep2
# shifts for:
# gr_rep1(120), gr_rep2(110), er_rep1(95), er_rep2(100)
#
# how variable names are structured 
# (terminology kept consistant with original MAnorm code):
# tmp_(merge_)(common_/unique_)peak(1/2/)_read(1a/1b/1).(bed/counts)
#
# tmp_      is used to remove all files at the end of the script
# merge_    using merged reference regions or pooled
# common_	refers to intersect of the two sets of peaks being compared (shared peaks)
# unique_	refers to setdiff of the two peaks being compared
# peak_ 	refers to current peak file of interest (1 or 2) if no number then 1+2
# read  	refers to current replicate (a/b) counted over in condition (1/2)
# bed       refers to sorted coordinates 
# counts	used for all variables that store counts (using coverageBed)
#####################################################################################################
# Todo: 
# switch for gzip
# make code work with no replicates [needs motifications to MAnorm3.R]
# resume by step# feature 
# bounds checking when shifting since can result in negative coordinates
# catch case when nothing is unique (same peaks file used)
# catch for directory already exists (overwrite)
# additional error checking (never ending)
# someday I might rewrite the whole program to run it within R without bedtools

fname=$1
if [ ! -d "$fname" ]
then #http://stackoverflow.com/questions/59838/how-to-check-if-a-directory-exists-in-a-shell-script
    mkdir $fname
elif [[ $(ls -A ./$fname/tmp* ) ]] #fix case where no tmp files found leads to error msg
then # quick error check: http://stackoverflow.com/questions/8921441/sh-test-for-existence-of-files
    echo "Please remove all tmp files before continuing"
    exit 1 
fi
echo $@ > ./$fname/${1}.log
echo `date` >> ./$fname/${1}.log  

echo "StepI: clean and sort input"
#dump and shift reads from bed (except chrM and *random*)
if [ $# -eq 7 ] # no replicates
then
    # # MAnorm3.R will fail without replicates!
    echo "This script does not currently work without replicates (will be updated in the future)"
    exit
    # echo "no replicates" >> ./$fname/${1}.log  
    # #newlines do not work in echo by default, must use echo -e which is not POSIX
    # options="$4 shift: $6\n$5 shift: $7"  
    # printf "%b\n" "$options" >> ./$fname/${1}.log # http://wiki.bash-hackers.org/commands/builtin/printf
	# #global variables must passed into awk using -v

    # #shift NOT extend
	# zcat $2 | sed 's/\s$//g' | awk -v fname="$fname" 'BEGIN {OFS="\t"}
		# {if ($1~/chr/ && $1 !="chrM" && $1 !~/random/ && $3>$2 && $2>0 && $3>0)
			# print $1,$2,$3 > "./"fname"/tmp_peak1.bed";
		# else 
			# print $0 > "./"fname"/tmp_peak1_exclude.bed"}' &
	# zcat $3 | sed 's/\s$//g' | awk -v fname="$fname" 'BEGIN {OFS="\t"}
		# {if ($$1~/chr/ && 1 !="chrM"  && $1 !~/random/ && $3>$2  && $2>0 && $3>0)
			# print $1,$2,$3 > "./"fname"/tmp_peak2.bed";
		# else 
			# print $0 > "./"fname"/tmp_peak2_exclude.bed"}' &
	# zcat $4 | sed 's/\s$//g' | awk -v fname="$fname" -v shift1=$6 'BEGIN {OFS="\t"}
		# {if ($1~/chr/ && $1 !="chrM" && $6=="+" && $1 !~/random/ && $3>$2  && $2>0 && $3>0)
			# print $1,$2+shift1,$3+shift1 > "./"fname"/tmp_read1.bed";
		# else if ($1~/chr/ && $1 !="chrM" && $6=="-" && $1 !~/random/ && $3>$2  && $2>shift1 && $3>shift1)
			# print $1,$2-shift1,$3-shift1 > "./"fname"/tmp_read1.bed";
		# else 
			# print $0 > "./"fname"/tmp_read1_exclude.bed"}' &
	# zcat $5 | sed 's/\s$//g' | awk -v fname="$fname" -v shift2=$7 'BEGIN {OFS="\t"}
		# {if ($1~/chr/ && $1 !="chrM" && $6=="+" && $1 !~/random/ && $3>$2  && $2>0 && $3>0)
			# print $1,$2+shift2,$3+shift2 > "./"fname"/tmp_read2.bed";
		# else if ($1~/chr/ && $1 !="chrM" && $6=="-" && $1 !~/random/ && $3>$2  && $2>shift2 && $3>shift2)
			# print $1,$2-shift2,$3-shift2 > "./"fname"/tmp_read2.bed";
		# else 
			# print $0 > "./"fname"/tmp_read2_exclude.bed"}' &
    # wait
    # #remember that sort is in non-ASCII order ie 1, 10, 11, 2, etc
    # sort -k1,1 -k2,3n ./$fname/tmp_peak1.bed > ./$fname/tmp_peak1s.bed &
    # sort -k1,1 -k2,3n ./$fname/tmp_peak2.bed > ./$fname/tmp_peak2s.bed &
    # sort -k1,1 -k2,3n ./$fname/tmp_read1.bed > ./$fname/tmp_read1s.bed &
    # sort -k1,1 -k2,3n ./$fname/tmp_read2.bed > ./$fname/tmp_read2s.bed &
    # wait
    # #cannot replace file that is being sorted
    # mv ./$fname/tmp_peak1s.bed ./$fname/tmp_peak1.bed
    # mv ./$fname/tmp_peak2s.bed ./$fname/tmp_peak2.bed
    # mv ./$fname/tmp_read1s.bed ./$fname/tmp_read1.bed
    # mv ./$fname/tmp_read2s.bed ./$fname/tmp_read2.bed
    # #write to log
	# echo "wc -l tmp_peak1.bed `wc -l ./$fname/tmp_peak1.bed` | cut -f 1 -d " "" >> ./$fname/${1}.log
	# echo "wc -l tmp_peak2.bed `wc -l ./$fname/tmp_peak2.bed` | cut -f 1 -d " "" >> ./$fname/${1}.log    
	# echo "wc -l tmp_read1.bed `wc -l ./$fname/tmp_read1.bed` | cut -f 1 -d " "" >> ./$fname/${1}.log
	# echo "wc -l tmp_read2.bed `wc -l ./$fname/tmp_read2.bed` | cut -f 1 -d " "" >> ./$fname/${1}.log
	# echo "wc -l tmp_read1_exclude.bed `wc -l ./$fname/tmp_read1_exclude.bed` | cut -f 1 -d " "" >> ./$fname/${1}.log
	# echo "wc -l tmp_read2_exclude.bed `wc -l ./$fname/tmp_read2_exclude.bed` | cut -f 1 -d " "" >> ./$fname/${1}.log
elif [ $# -eq 11 ] # 1 replicate
then 
	#mostly copied from above
	#add functionality to combine replicates
	#"proper" way would be to downscale replicates so that both contribute equally (not done)
    echo "duplicates" >> ./$fname/${1}.log
    #newlines do not work in echo by default, must use echo -e which is not POSIX
    #on ubuntu echo automatically adds newline
    options="$4 shift: $8 \n$5 shift: $9 \n$6 shift: ${10} \n$7 shift: ${11}"  
    printf "%b\n" "$options" >> ./$fname/${1}.log # http://wiki.bash-hackers.org/commands/builtin/printf
    
	zcat $2 | sed 's/\s$//g' | awk -v fname="$fname" 'BEGIN {OFS="\t"}
		{if ($1~/chr/ && $1 !="chrM" && $1 !~/random/ && $3>$2 && $2>0 && $3>0)
			print $1,$2,$3 > "./"fname"/tmp_peak1.bed";
		else 
			print $0 > "./"fname"/tmp_peak1_exclude.bed"}' &
	zcat $3 | sed 's/\s$//g' | awk -v fname="$fname" 'BEGIN {OFS="\t"}
		{if ($$1~/chr/ && 1 !="chrM" && $1 !~/random/ && $3>$2 && $2>0 && $3>0)
			print $1,$2,$3 > "./"fname"/tmp_peak2.bed";
		else 
			print $0 > "./"fname"/tmp_peak2_exclude.bed"}' &
	#shift NOT extend, ideally would also check for end of chr
	zcat $4 | sed 's/\s$//g' | awk -v fname="$fname" -v shift1a=$8 'BEGIN {OFS="\t"}
		{if ($1~/chr/ && $1 !="chrM" && $6=="+" && $1 !~/random/ && $3>$2 && $2>0 && $3>0)
			print $1,$2+shift1a,$3+shift1a > "./"fname"/tmp_read1a.bed";
		else if ($1~/chr/ && $1 !="chrM" && $6=="-" && $1 !~/random/ && $3>$2 && $2>shift1a && $3>shift1a)
			print $1,$2-shift1a,$3-shift1a > "./"fname"/tmp_read1a.bed";
		else 
			print $0 > "./"fname"/tmp_read1a_exclude.bed"}' &
	zcat $5 | sed 's/\s$//g' | awk -v fname="$fname" -v shift1b=$9 'BEGIN {OFS="\t"}
		{if ($1~/chr/ && $1 !="chrM" && $6=="+" && $1 !~/random/ && $3>$2 && $2>0 && $3>0)
			print $1,$2+shift1b,$3+shift1b > "./"fname"/tmp_read1b.bed";
		else if ($1~/chr/ && $1 !="chrM" && $6=="-" && $1 !~/random/ && $3>$2 && $2>shift1b && $3>shift1b)
			print $1,$2-shift1b,$3-shift1b > "./"fname"/tmp_read1b.bed";
		else 
			print $0 > "./"fname"/tmp_read1b_exclude.bed"}'	&
	zcat $6 | sed 's/\s$//g' | awk -v fname="$fname" -v shift2a=${10} 'BEGIN {OFS="\t"}
		{if ($1~/chr/ && $1 !="chrM" && $6=="+" && $1 !~/random/ && $3>$2 && $2>0 && $3>0)
			print $1,$2+shift2a,$3+shift2a > "./"fname"/tmp_read2a.bed";
		else if ($1~/chr/ && $1 !="chrM" && $6=="-" && $1 !~/random/ && $3>$2 && $2>shift2a && $3>shift2a)
			print $1,$2-shift2a,$3-shift2a > "./"fname"/tmp_read2a.bed";
		else 
			print $0 > "./"fname"/tmp_read2a_exclude.bed"}' &
	zcat $7 | sed 's/\s$//g' | awk -v fname="$fname" -v shift2b=${11} 'BEGIN {OFS="\t"}
		{if ($1~/chr/ && $1 !="chrM" && $6=="+" && $1 !~/random/ && $3>$2 && $2>0 && $3>0)
			print $1,$2+shift2b,$3+shift2b > "./"fname"/tmp_read2b.bed";
		else if ($1~/chr/ && $1 !="chrM" && $6=="-" && $1 !~/random/ && $3>$2 && $2>shift2b && $3>shift2b)
			print $1,$2-shift2b,$3-shift2b > "./"fname"/tmp_read2b.bed";
		else 
			print $0 > "./"fname"/tmp_read2b_exclude.bed"}' &
    wait
    
    # output pooled reads for MAnorm.bed step
	#cat ./$fname/tmp_read1a.bed ./$fname/tmp_read1b.bed > ./$fname/tmp_read1.bed
	#cat ./$fname/tmp_read2a.bed ./$fname/tmp_read2b.bed > ./$fname/tmp_read2.bed
    
    #sort
    sort -k1,1 -k2,3n ./$fname/tmp_peak1.bed > ./$fname/tmp_peak1s.bed &
    sort -k1,1 -k2,3n ./$fname/tmp_peak2.bed > ./$fname/tmp_peak2s.bed &
    #sort -k1,1 -k2,3n ./$fname/tmp_read1.bed > ./$fname/tmp_read1s.bed &
    #sort -k1,1 -k2,3n ./$fname/tmp_read2.bed > ./$fname/tmp_read2s.bed &
    sort -k1,1 -k2,3n ./$fname/tmp_read1a.bed > ./$fname/tmp_read1as.bed &
    sort -k1,1 -k2,3n ./$fname/tmp_read1b.bed > ./$fname/tmp_read1bs.bed &
    sort -k1,1 -k2,3n ./$fname/tmp_read2a.bed > ./$fname/tmp_read2as.bed &
    sort -k1,1 -k2,3n ./$fname/tmp_read2b.bed > ./$fname/tmp_read2bs.bed &
    wait
    
    mv ./$fname/tmp_peak1s.bed ./$fname/tmp_peak1.bed
    mv ./$fname/tmp_peak2s.bed ./$fname/tmp_peak2.bed
    mv ./$fname/tmp_read1as.bed ./$fname/tmp_read1a.bed
    mv ./$fname/tmp_read1bs.bed ./$fname/tmp_read1b.bed
    mv ./$fname/tmp_read2as.bed ./$fname/tmp_read2a.bed
    mv ./$fname/tmp_read2bs.bed ./$fname/tmp_read2b.bed

    #write to log
	echo "wc -l tmp_peak1.bed `wc -l ./$fname/tmp_peak1.bed | cut -f 1 -d " "`" >> ./$fname/${1}.log
	echo "wc -l tmp_peak2.bed `wc -l ./$fname/tmp_peak2.bed | cut -f 1 -d " "`" >> ./$fname/${1}.log    
	echo "wc -l tmp_read1a.bed `wc -l ./$fname/tmp_read1a.bed | cut -f 1 -d " "`" >> ./$fname/${1}.log
	echo "wc -l tmp_read1b.bed `wc -l ./$fname/tmp_read1b.bed | cut -f 1 -d " "`" >> ./$fname/${1}.log
	echo "wc -l tmp_read2a.bed `wc -l ./$fname/tmp_read2a.bed | cut -f 1 -d " "`" >> ./$fname/${1}.log
	echo "wc -l tmp_read2b.bed `wc -l ./$fname/tmp_read2b.bed | cut -f 1 -d " "`" >> ./$fname/${1}.log
	echo "wc -l tmp_read1a_exclude.bed `wc -l ./$fname/tmp_read1a_exclude.bed | cut -f 1 -d " "`" >> ./$fname/${1}.log
	echo "wc -l tmp_read1b_exclude.bed `wc -l ./$fname/tmp_read1b_exclude.bed | cut -f 1 -d " "`" >> ./$fname/${1}.log
	echo "wc -l tmp_read2a_exclude.bed `wc -l ./$fname/tmp_read2a_exclude.bed | cut -f 1 -d " "`" >> ./$fname/${1}.log
	echo "wc -l tmp_read2b_exclude.bed `wc -l ./$fname/tmp_read2b_exclude.bed | cut -f 1 -d " "`" >> ./$fname/${1}.log
else
    #basename command will remove the directory part of the file as well as the suffix
    echo "Error: $# arguments found"
    echo "Usage: `basename $0` foldername peak1.bed peak2.bed read1.bed read2.bed shift1 shift2"
    echo "or `basename $0` foldername peak1.bed peak2.bed read1a.bed read1b.bed read2a.bed read2b.bed shift1a shift1b shift2a shift2b"
    exit
fi

#####################################################################################################
echo "StepII: classify common or unique peaks"
#remember that b is put in memory, since input is sorted, output should also be sorted
intersectBed -a ./$fname/tmp_peak1.bed -b ./$fname/tmp_peak2.bed -u > ./$fname/tmp_common_peak1.bed &
intersectBed -a ./$fname/tmp_peak2.bed -b ./$fname/tmp_peak1.bed -u > ./$fname/tmp_common_peak2.bed &
intersectBed -a ./$fname/tmp_peak1.bed -b ./$fname/tmp_peak2.bed -v > ./$fname/tmp_unique_peak1.bed &
intersectBed -a ./$fname/tmp_peak2.bed -b ./$fname/tmp_peak1.bed -v > ./$fname/tmp_unique_peak2.bed &
wait

# WARNING!!! must sort before merging, above sort cmd is faster but this is more clear
cat ./$fname/tmp_common_peak1.bed ./$fname/tmp_common_peak2.bed | sortBed | mergeBed > ./$fname/tmp_common_peak.bed & 

#same as above, might be a faster way around this but keep running into issues
mergeBed -i ./$fname/tmp_peak1.bed > ./$fname/tmp_merge_peak1.bed & #often the same
mergeBed -i ./$fname/tmp_peak2.bed > ./$fname/tmp_merge_peak2.bed & #often the same
wait
intersectBed -a ./$fname/tmp_merge_peak1.bed -b ./$fname/tmp_merge_peak2.bed -u > ./$fname/tmp_merge_common_peak1.bed &
intersectBed -a ./$fname/tmp_merge_peak2.bed -b ./$fname/tmp_merge_peak1.bed -u > ./$fname/tmp_merge_common_peak2.bed &
intersectBed -a ./$fname/tmp_merge_peak1.bed -b ./$fname/tmp_merge_peak2.bed -v > ./$fname/tmp_merge_unique_peak1.bed &
intersectBed -a ./$fname/tmp_merge_peak2.bed -b ./$fname/tmp_merge_peak1.bed -v > ./$fname/tmp_merge_unique_peak2.bed &
wait
cat ./$fname/tmp_merge_common_peak1.bed ./$fname/tmp_merge_common_peak2.bed | sortBed | mergeBed > ./$fname/tmp_merge_common_peak.bed

echo "wc -l tmp_common_peak1.bed `wc -l ./$fname/tmp_common_peak1.bed | cut -f 1 -d " "`" >> ./$fname/${1}.log
echo "wc -l tmp_common_peak2.bed `wc -l ./$fname/tmp_common_peak2.bed | cut -f 1 -d " "`" >> ./$fname/${1}.log
echo "wc -l tmp_unique_peak1.bed `wc -l ./$fname/tmp_unique_peak1.bed | cut -f 1 -d " "`" >> ./$fname/${1}.log
echo "wc -l tmp_unique_peak2.bed `wc -l ./$fname/tmp_unique_peak2.bed | cut -f 1 -d " "`" >> ./$fname/${1}.log
echo "wc -l tmp_common_peak.bed `wc -l ./$fname/tmp_common_peak.bed | cut -f 1 -d " "`" >> ./$fname/${1}.log
echo "wc -l tmp_merge_peak1.bed `wc -l ./$fname/tmp_merge_peak1.bed | cut -f 1 -d " "`" >> ./$fname/${1}.log
echo "wc -l tmp_merge_peak2.bed `wc -l ./$fname/tmp_merge_peak2.bed | cut -f 1 -d " "`" >> ./$fname/${1}.log
echo "wc -l tmp_merge_common_peak1.bed `wc -l ./$fname/tmp_merge_common_peak1.bed | cut -f 1 -d " "`" >> ./$fname/${1}.log
echo "wc -l tmp_merge_common_peak2.bed `wc -l ./$fname/tmp_merge_common_peak2.bed | cut -f 1 -d " "`" >> ./$fname/${1}.log
echo "wc -l tmp_merge_unique_peak1.bed `wc -l ./$fname/tmp_merge_unique_peak1.bed | cut -f 1 -d " "`" >> ./$fname/${1}.log
echo "wc -l tmp_merge_unique_peak2.bed `wc -l ./$fname/tmp_merge_unique_peak2.bed | cut -f 1 -d " "`" >> ./$fname/${1}.log
echo "wc -l tmp_merge_common_peak.bed `wc -l ./$fname/tmp_merge_common_peak.bed | cut -f 1 -d " "`" >> ./$fname/${1}.log 

if [ ! -f ./$fname/tmp_common_peak1.bed ];
then
    echo "Error during StepII"
    exit 1
fi
#####################################################################################################
echo "StepIII: count peak read"
if [ -f ./$fname/tmp_MAnorm.bed ];
then
	rm ./$fname/tmp_MAnorm.bed
fi
#appends are not atomic, cannot parallel >> so split up the tasks
# coverageBed -a ./$fname/tmp_read1.bed -b ./$fname/tmp_unique_peak1.bed | awk -v fname="$fname" 'BEGIN {OFS="\t"} {print $1,$2,$3,"unique_peak1" >> "./"fname"/tmp_MAnorm.bed"; print $4 > "./"fname"/tmp_unique_peak1_count_read1"}'
# coverageBed -a ./$fname/tmp_read2.bed -b ./$fname/tmp_unique_peak1.bed | awk -v fname="$fname" '{print $4 > "./"fname"/tmp_unique_peak1_count_read2"}'
# coverageBed -a ./$fname/tmp_read1.bed -b ./$fname/tmp_common_peak1.bed | awk -v fname="$fname" 'BEGIN {OFS="\t"} {print $1,$2,$3,"common_peak1" >> "./"fname"/tmp_MAnorm.bed";print $4 > "./"fname"/tmp_common_peak1_count_read1"}'
# coverageBed -a ./$fname/tmp_read2.bed -b ./$fname/tmp_common_peak1.bed | awk -v fname="$fname" '{print $4 > "./"fname"/tmp_common_peak1_count_read2"}'
# coverageBed -a ./$fname/tmp_read1.bed -b ./$fname/tmp_common_peak2.bed | awk -v fname="$fname" 'BEGIN {OFS="\t"} {print $1,$2,$3,"common_peak2"  >> "./"fname"/tmp_MAnorm.bed";print $4 > "./"fname"/tmp_common_peak2_count_read1"}'
# coverageBed -a ./$fname/tmp_read2.bed -b ./$fname/tmp_common_peak2.bed | awk -v fname="$fname" '{print $4 > "./"fname"/tmp_common_peak2_count_read2"}'
# coverageBed -a ./$fname/tmp_read1.bed -b ./$fname/tmp_unique_peak2.bed | awk -v fname="$fname" 'BEGIN {OFS="\t"} {print $1,$2,$3,"unique_peak2">> "./"fname"/tmp_MAnorm.bed";print $4 > "./"fname"/tmp_unique_peak2_count_read1"}'
# coverageBed -a ./$fname/tmp_read2.bed -b ./$fname/tmp_unique_peak2.bed | awk -v fname="$fname" '{print $4 > "./"fname"/tmp_unique_peak2_count_read2"}'

if [ $# -eq 7 ]
then
    coverageBed -a ./$fname/tmp_read1.bed -b ./$fname/tmp_unique_peak1.bed | cut -f 1-4 > ./$fname/tmp_unique_peak1_read1.counts &
    coverageBed -a ./$fname/tmp_read2.bed -b ./$fname/tmp_unique_peak1.bed | cut -f 1-4 > ./$fname/tmp_unique_peak1_read2.counts &
    coverageBed -a ./$fname/tmp_read1.bed -b ./$fname/tmp_common_peak1.bed | cut -f 1-4 > ./$fname/tmp_common_peak1_read1.counts &
    coverageBed -a ./$fname/tmp_read2.bed -b ./$fname/tmp_common_peak1.bed | cut -f 1-4 > ./$fname/tmp_common_peak1_read2.counts &
    coverageBed -a ./$fname/tmp_read1.bed -b ./$fname/tmp_common_peak2.bed | cut -f 1-4 > ./$fname/tmp_common_peak2_read1.counts &
    coverageBed -a ./$fname/tmp_read2.bed -b ./$fname/tmp_common_peak2.bed | cut -f 1-4 > ./$fname/tmp_common_peak2_read2.counts &
    coverageBed -a ./$fname/tmp_read1.bed -b ./$fname/tmp_unique_peak2.bed | cut -f 1-4 > ./$fname/tmp_unique_peak2_read1.counts &
    coverageBed -a ./$fname/tmp_read2.bed -b ./$fname/tmp_unique_peak2.bed | cut -f 1-4 > ./$fname/tmp_unique_peak2_read2.counts &
    wait
elif [ $# -eq 11 ]
then
	coverageBed -a ./$fname/tmp_read1a.bed -b ./$fname/tmp_unique_peak1.bed | cut -f 1-4 > ./$fname/tmp_unique_peak1_read1a.counts &
	coverageBed -a ./$fname/tmp_read2a.bed -b ./$fname/tmp_unique_peak1.bed | cut -f 1-4 > ./$fname/tmp_unique_peak1_read2a.counts &
	coverageBed -a ./$fname/tmp_read1a.bed -b ./$fname/tmp_common_peak1.bed | cut -f 1-4 > ./$fname/tmp_common_peak1_read1a.counts &
	coverageBed -a ./$fname/tmp_read2a.bed -b ./$fname/tmp_common_peak1.bed | cut -f 1-4 > ./$fname/tmp_common_peak1_read2a.counts &
	coverageBed -a ./$fname/tmp_read1a.bed -b ./$fname/tmp_common_peak2.bed | cut -f 1-4 > ./$fname/tmp_common_peak2_read1a.counts &
	coverageBed -a ./$fname/tmp_read2a.bed -b ./$fname/tmp_common_peak2.bed | cut -f 1-4 > ./$fname/tmp_common_peak2_read2a.counts &
	coverageBed -a ./$fname/tmp_read1a.bed -b ./$fname/tmp_unique_peak2.bed | cut -f 1-4 > ./$fname/tmp_unique_peak2_read1a.counts &
	coverageBed -a ./$fname/tmp_read2a.bed -b ./$fname/tmp_unique_peak2.bed | cut -f 1-4 > ./$fname/tmp_unique_peak2_read2a.counts &	
    
	coverageBed -a ./$fname/tmp_read1b.bed -b ./$fname/tmp_unique_peak1.bed | cut -f 1-4 > ./$fname/tmp_unique_peak1_read1b.counts &
	coverageBed -a ./$fname/tmp_read2b.bed -b ./$fname/tmp_unique_peak1.bed | cut -f 1-4 > ./$fname/tmp_unique_peak1_read2b.counts &
	coverageBed -a ./$fname/tmp_read1b.bed -b ./$fname/tmp_common_peak1.bed | cut -f 1-4 > ./$fname/tmp_common_peak1_read1b.counts &
	coverageBed -a ./$fname/tmp_read2b.bed -b ./$fname/tmp_common_peak1.bed | cut -f 1-4 > ./$fname/tmp_common_peak1_read2b.counts &
	coverageBed -a ./$fname/tmp_read1b.bed -b ./$fname/tmp_common_peak2.bed | cut -f 1-4 > ./$fname/tmp_common_peak2_read1b.counts &
	coverageBed -a ./$fname/tmp_read2b.bed -b ./$fname/tmp_common_peak2.bed | cut -f 1-4 > ./$fname/tmp_common_peak2_read2b.counts &
	coverageBed -a ./$fname/tmp_read1b.bed -b ./$fname/tmp_unique_peak2.bed | cut -f 1-4 > ./$fname/tmp_unique_peak2_read1b.counts &
	coverageBed -a ./$fname/tmp_read2b.bed -b ./$fname/tmp_unique_peak2.bed | cut -f 1-4 > ./$fname/tmp_unique_peak2_read2b.counts &	
	wait
	# .counts have 4 col: coordinates + count
	cat ./$fname/tmp_common_peak1_read1a.counts ./$fname/tmp_common_peak2_read1a.counts > ./$fname/tmp_common_peak_read1a.counts &
	cat ./$fname/tmp_common_peak1_read2a.counts ./$fname/tmp_common_peak2_read2a.counts > ./$fname/tmp_common_peak_read2a.counts &
	cat ./$fname/tmp_common_peak1_read1b.counts ./$fname/tmp_common_peak2_read1b.counts > ./$fname/tmp_common_peak_read1b.counts &
	cat ./$fname/tmp_common_peak1_read2b.counts ./$fname/tmp_common_peak2_read2b.counts > ./$fname/tmp_common_peak_read2b.counts &
	cat ./$fname/tmp_unique_peak1_read1a.counts ./$fname/tmp_common_peak1_read1a.counts ./$fname/tmp_common_peak2_read1a.counts ./$fname/tmp_unique_peak2_read1a.counts > ./$fname/tmp_peak_read1a.counts &
	cat ./$fname/tmp_unique_peak1_read2a.counts ./$fname/tmp_common_peak1_read2a.counts ./$fname/tmp_common_peak2_read2a.counts ./$fname/tmp_unique_peak2_read2a.counts > ./$fname/tmp_peak_read2a.counts &
	cat ./$fname/tmp_unique_peak1_read1b.counts ./$fname/tmp_common_peak1_read1b.counts ./$fname/tmp_common_peak2_read1b.counts ./$fname/tmp_unique_peak2_read1b.counts > ./$fname/tmp_peak_read1b.counts &
	cat ./$fname/tmp_unique_peak1_read2b.counts ./$fname/tmp_common_peak1_read2b.counts ./$fname/tmp_common_peak2_read2b.counts ./$fname/tmp_unique_peak2_read2b.counts > ./$fname/tmp_peak_read2b.counts &	
	#http://nixtricks.wordpress.com/2009/10/24/awk-add-another-column-from-a-second-fil/
    awk -v fname="$fname" 'BEGIN {OFS="\t"} { counta=$4; getline < ("./"fname"/tmp_unique_peak1_read1b.counts"); print $1,$2,$3,$4+$counta}' ./$fname/tmp_unique_peak1_read1a.counts > ./$fname/tmp_unique_peak1_read1.counts &
    awk -v fname="$fname" 'BEGIN {OFS="\t"} { counta=$4; getline < ("./"fname"/tmp_unique_peak1_read2b.counts"); print $1,$2,$3,$4+$counta}' ./$fname/tmp_unique_peak1_read2a.counts > ./$fname/tmp_unique_peak1_read2.counts &
    awk -v fname="$fname" 'BEGIN {OFS="\t"} { counta=$4; getline < ("./"fname"/tmp_common_peak1_read1b.counts"); print $1,$2,$3,$4+$counta}' ./$fname/tmp_common_peak1_read1a.counts > ./$fname/tmp_common_peak1_read1.counts &
    awk -v fname="$fname" 'BEGIN {OFS="\t"} { counta=$4; getline < ("./"fname"/tmp_common_peak1_read2b.counts"); print $1,$2,$3,$4+$counta}' ./$fname/tmp_common_peak1_read2a.counts > ./$fname/tmp_common_peak1_read2.counts &
    awk -v fname="$fname" 'BEGIN {OFS="\t"} { counta=$4; getline < ("./"fname"/tmp_unique_peak2_read1b.counts"); print $1,$2,$3,$4+$counta}' ./$fname/tmp_unique_peak2_read1a.counts > ./$fname/tmp_common_peak2_read1.counts &
    awk -v fname="$fname" 'BEGIN {OFS="\t"} { counta=$4; getline < ("./"fname"/tmp_unique_peak2_read2b.counts"); print $1,$2,$3,$4+$counta}' ./$fname/tmp_unique_peak2_read2a.counts > ./$fname/tmp_common_peak2_read2.counts &
    awk -v fname="$fname" 'BEGIN {OFS="\t"} { counta=$4; getline < ("./"fname"/tmp_common_peak2_read1b.counts"); print $1,$2,$3,$4+$counta}' ./$fname/tmp_common_peak2_read1a.counts > ./$fname/tmp_unique_peak2_read1.counts &
    awk -v fname="$fname" 'BEGIN {OFS="\t"} { counta=$4; getline < ("./"fname"/tmp_common_peak2_read2b.counts"); print $1,$2,$3,$4+$counta}' ./$fname/tmp_common_peak2_read2a.counts > ./$fname/tmp_unique_peak2_read2.counts &
    wait
fi

#not sure how to do this w/out awk
awk 'BEGIN {OFS="\t"} {print $1,$2,$3,"unique_peak1"}' ./$fname/tmp_unique_peak1_read1.counts >>  ./$fname/tmp_MAnorm.bed
awk 'BEGIN {OFS="\t"} {print $1,$2,$3,"common_peak1"}' ./$fname/tmp_common_peak1_read1.counts >>  ./$fname/tmp_MAnorm.bed
awk 'BEGIN {OFS="\t"} {print $1,$2,$3,"common_peak2"}' ./$fname/tmp_common_peak2_read1.counts >>  ./$fname/tmp_MAnorm.bed
awk 'BEGIN {OFS="\t"} {print $1,$2,$3,"unique_peak2"}' ./$fname/tmp_unique_peak2_read1.counts >>  ./$fname/tmp_MAnorm.bed

cat ./$fname/tmp_common_peak1_read1.counts ./$fname/tmp_common_peak2_read1.counts | cut -f 4 > ./$fname/tmp_common_peak_read1.counts &
cat ./$fname/tmp_common_peak1_read2.counts ./$fname/tmp_common_peak2_read2.counts | cut -f 4 > ./$fname/tmp_common_peak_read2.counts &
cat ./$fname/tmp_unique_peak1_read1.counts ./$fname/tmp_common_peak1_read1.counts ./$fname/tmp_common_peak2_read1.counts ./$fname/tmp_unique_peak2_read1.counts | cut -f 4 > ./$fname/tmp_peak_read1.counts &
cat ./$fname/tmp_unique_peak1_read2.counts ./$fname/tmp_common_peak1_read2.counts ./$fname/tmp_common_peak2_read2.counts ./$fname/tmp_unique_peak2_read2.counts | cut -f 4 > ./$fname/tmp_peak_read2.counts &
wait

# Merge overlaps
if [ -f ./$fname/tmp_MAnorm_merge.bed ];
then
	rm ./$fname/tmp_MAnorm_merge.bed
fi

if [ $# -eq 7 ]
then
    coverageBed -a ./$fname/tmp_read1.bed -b ./$fname/tmp_merge_common_peak.bed | cut -f 1-4 > ./$fname/tmp_merge_common_read1.counts &
    coverageBed -a ./$fname/tmp_read2.bed -b ./$fname/tmp_merge_common_peak.bed | cut -f 1-4 > ./$fname/tmp_merge_common_read2.counts &
    
    coverageBed -a ./$fname/tmp_read1.bed -b ./$fname/tmp_merge_unique_peak1.bed | cut -f 1-4 > ./$fname/tmp_merge_unique_peak1_read1.counts &
    coverageBed -a ./$fname/tmp_read2.bed -b ./$fname/tmp_merge_unique_peak1.bed | cut -f 1-4 > ./$fname/tmp_merge_unique_peak1_read2.counts &
    coverageBed -a ./$fname/tmp_read1.bed -b ./$fname/tmp_merge_unique_peak2.bed | cut -f 1-4 > ./$fname/tmp_merge_unique_peak2_read1.counts &
    coverageBed -a ./$fname/tmp_read2.bed -b ./$fname/tmp_merge_unique_peak2.bed | cut -f 1-4 > ./$fname/tmp_merge_unique_peak2_read2.counts &
elif [ $# -eq 11 ]
then
	coverageBed -a ./$fname/tmp_read1a.bed -b ./$fname/tmp_merge_common_peak.bed | cut -f 1-4 > ./$fname/tmp_merge_common_read1a.counts &
	coverageBed -a ./$fname/tmp_read2a.bed -b ./$fname/tmp_merge_common_peak.bed | cut -f 1-4 > ./$fname/tmp_merge_common_read2a.counts &
	coverageBed -a ./$fname/tmp_read1b.bed -b ./$fname/tmp_merge_common_peak.bed | cut -f 1-4 > ./$fname/tmp_merge_common_read1b.counts &
	coverageBed -a ./$fname/tmp_read2b.bed -b ./$fname/tmp_merge_common_peak.bed | cut -f 1-4 > ./$fname/tmp_merge_common_read2b.counts &
    
    coverageBed -a ./$fname/tmp_read1a.bed -b ./$fname/tmp_merge_unique_peak1.bed | cut -f 1-4 > ./$fname/tmp_merge_unique_peak1_read1a.counts &
    coverageBed -a ./$fname/tmp_read1b.bed -b ./$fname/tmp_merge_unique_peak1.bed | cut -f 1-4 > ./$fname/tmp_merge_unique_peak1_read1b.counts &
    coverageBed -a ./$fname/tmp_read2a.bed -b ./$fname/tmp_merge_unique_peak1.bed | cut -f 1-4 > ./$fname/tmp_merge_unique_peak1_read2a.counts &
    coverageBed -a ./$fname/tmp_read2b.bed -b ./$fname/tmp_merge_unique_peak1.bed | cut -f 1-4 > ./$fname/tmp_merge_unique_peak1_read2b.counts &
    coverageBed -a ./$fname/tmp_read1a.bed -b ./$fname/tmp_merge_unique_peak2.bed | cut -f 1-4 > ./$fname/tmp_merge_unique_peak2_read1a.counts &
    coverageBed -a ./$fname/tmp_read1b.bed -b ./$fname/tmp_merge_unique_peak2.bed | cut -f 1-4 > ./$fname/tmp_merge_unique_peak2_read1b.counts &
    coverageBed -a ./$fname/tmp_read2a.bed -b ./$fname/tmp_merge_unique_peak2.bed | cut -f 1-4 > ./$fname/tmp_merge_unique_peak2_read2a.counts &
    coverageBed -a ./$fname/tmp_read2b.bed -b ./$fname/tmp_merge_unique_peak2.bed | cut -f 1-4 > ./$fname/tmp_merge_unique_peak2_read2b.counts &
	wait
    
	cat ./$fname/tmp_merge_unique_peak1_read1a.counts ./$fname/tmp_merge_common_read1a.counts ./$fname/tmp_merge_unique_peak2_read1a.counts > ./$fname/tmp_merge_common_peak_read1a.counts &
	cat ./$fname/tmp_merge_unique_peak1_read2a.counts ./$fname/tmp_merge_common_read2a.counts ./$fname/tmp_merge_unique_peak2_read2a.counts > ./$fname/tmp_merge_common_peak_read2a.counts &
	cat ./$fname/tmp_merge_unique_peak1_read1b.counts ./$fname/tmp_merge_common_read1b.counts ./$fname/tmp_merge_unique_peak2_read1b.counts > ./$fname/tmp_merge_common_peak_read1b.counts &
	cat ./$fname/tmp_merge_unique_peak1_read2b.counts ./$fname/tmp_merge_common_read2b.counts ./$fname/tmp_merge_unique_peak2_read2b.counts > ./$fname/tmp_merge_common_peak_read2b.counts &
    
    awk -v fname="$fname" 'BEGIN {OFS="\t"} { counta=$4; getline < ("./"fname"/tmp_merge_common_read1b.counts"); print $1,$2,$3,$4+$counta}' ./$fname/tmp_merge_common_read1a.counts > ./$fname/tmp_merge_common_read1.counts &
    awk -v fname="$fname" 'BEGIN {OFS="\t"} { counta=$4; getline < ("./"fname"/tmp_merge_common_read2b.counts"); print $1,$2,$3,$4+$counta}' ./$fname/tmp_merge_common_read2a.counts > ./$fname/tmp_merge_common_read2.counts &    
    awk -v fname="$fname" 'BEGIN {OFS="\t"} { counta=$4; getline < ("./"fname"/tmp_merge_unique_peak1_read1b.counts"); print $1,$2,$3,$4+$counta}' ./$fname/tmp_merge_unique_peak1_read1a.counts > ./$fname/tmp_merge_unique_peak1_read1.counts &
    awk -v fname="$fname" 'BEGIN {OFS="\t"} { counta=$4; getline < ("./"fname"/tmp_merge_unique_peak1_read2b.counts"); print $1,$2,$3,$4+$counta}' ./$fname/tmp_merge_unique_peak1_read2a.counts > ./$fname/tmp_merge_unique_peak1_read2.counts &    
    awk -v fname="$fname" 'BEGIN {OFS="\t"} { counta=$4; getline < ("./"fname"/tmp_merge_unique_peak2_read1b.counts"); print $1,$2,$3,$4+$counta}' ./$fname/tmp_merge_unique_peak2_read1a.counts > ./$fname/tmp_merge_unique_peak2_read1.counts &
    awk -v fname="$fname" 'BEGIN {OFS="\t"} { counta=$4; getline < ("./"fname"/tmp_merge_unique_peak2_read2b.counts"); print $1,$2,$3,$4+$counta}' ./$fname/tmp_merge_unique_peak2_read2a.counts > ./$fname/tmp_merge_unique_peak2_read2.counts &    
fi
wait

cat ./$fname/tmp_merge_unique_peak1_read1.counts | awk -v fname="$fname" 'BEGIN {OFS="\t"} {print $1,$2,$3,"unique_peak1" >> "./"fname"/tmp_MAnorm_merge.bed"}'
cat ./$fname/tmp_merge_common_read1.counts | awk 'BEGIN {OFS="\t"} {print $1,$2,$3,"merged_common_peak"}' >> ./$fname/tmp_MAnorm_merge.bed
cat ./$fname/tmp_merge_unique_peak2_read1.counts | awk -v fname="$fname" 'BEGIN {OFS="\t"} {print $1,$2,$3,"unique_peak2" >> "./"fname"/tmp_MAnorm_merge.bed"}'

cat ./$fname/tmp_merge_unique_peak1_read1.counts ./$fname/tmp_merge_common_read1.counts ./$fname/tmp_merge_unique_peak2_read1.counts | cut -f 4 > ./$fname/tmp_merge_common_peak_read1.counts &
cat ./$fname/tmp_merge_unique_peak1_read2.counts ./$fname/tmp_merge_common_read2.counts ./$fname/tmp_merge_unique_peak2_read2.counts | cut -f 4 > ./$fname/tmp_merge_common_peak_read2.counts &
wait

#####################################################################################################
echo "StepIV: normalize using common peaks"
#R --vanilla MAnorm.r >Rcommand.out 
#R CMD BATCH ../MAnorm.r Rcommand.out
cd ./$fname/
#R --no-restore --save < ../MAnorm2.r ${1}_Rcommand.out 
#doesnt load .R that might slow things down but still loads usr+global env 
#doesnt redirect nicely
R CMD BATCH ../MAnorm3.R ${1}_Rcommand.out 
#echo "R CURRENTLY DISABLED!!!"
cd ..

#This is not a proper wig, I think it is a bed (check)
#for format definition see: http://genome.ucsc.edu/goldenPath/help/wiggle.html
#This is just a bed + score
#awk -v fname="$fname" 'BEGIN{OFS="\t"}{if($4~/1/) print $1,$2,$3,$7 > "./"fname"/${1}_sample1_peaks.wig"}' MAnorm_result.xls
#awk -v fname="$fname" 'BEGIN{OFS="\t"}{if($4~/2/) print $1,$2,$3,$7 > "./"fname"/${1}_sample2_peaks.wig"}' MAnorm_result.xls

#####################################################################################################
#cleanup
#echo "CLEANUP DISABLED"
if [[ $(ls -A ./$fname/*xls ) ]] 
then
    rm ./$fname/tmp_*
else
    echo "MAnorm3.R did not create excel output file check ${fname}_Rcommand.out"
    #return 1 #use for sourcing
    exit 1
fi
