#http://bcb.dfci.harvard.edu/~gcyuan/MAnorm/R_tutorial.html
#ex. MAnorm3.sh foldername gr_peak.bed er_peak.bed gr_rep1.bed gr_rep2.bed er_rep1.bed er_rep2.bed 120 110 95 100
# peak1.bed: sample 1 peak list
# peak2.bed: sample 2 peak list
# read1.bed: sample 1 read used for mapping to peaks
# read2.bed: sample 2 read used for mapping to peaks
# peak1_dump.bed: sample 1 peaks not used
# peak2_dump.bed: sample 2 peaks not used
# read1_dump.bed: sample 1 read NOT used for mapping to peaks
# read2_dump.bed: sample 2 read NOT used for mapping to peaks
#####################################################################################################
# quick error check
# http://stackoverflow.com/questions/8921441/sh-test-for-existence-of-files
if [[ $(ls -A tmp* ) ]] 
then
    echo "Please remove all tmp files before continuing"
    return 1
fi

#dump and shift reads from bed (except chrM and *random*)
if [ $# -eq 7 ]
then
	fname=$1
	#TODO: write catch for directory creation
	mkdir $fname 
	echo $@ > ./$fname/${1}.log
    echo `date` >> ./$fname/${1}.log    
	echo "StepI: clean input"
    #newlines do not work in echo by default, must use echo -e which is not POSIX
    options="$4 shift: $6\n$5 shift: $7"  
    printf "%b\n" "$options" >> ./$fname/${1}.log # http://wiki.bash-hackers.org/commands/builtin/printf
	#global variables must passed into awk using -v
	zcat $2 | sed 's/\s$//g' | awk -v fname="$fname" 'BEGIN {OFS="\t"}
		{if ($1~/chr/ && $1 !="chrM" && $1 !~/random/ && $3>$2 && $2>0 && $3>0)
			print $1,$2,$3 > "./"fname"/tmp_peak1.bed";
		else 
			print $0 > "./"fname"/tmp_peak1_exclude.bed"}' &
	zcat $3 | sed 's/\s$//g' | awk -v fname="$fname" 'BEGIN {OFS="\t"}
		{if ($$1~/chr/ && 1 !="chrM"  && $1 !~/random/ && $3>$2  && $2>0 && $3>0)
			print $1,$2,$3 > "./"fname"/tmp_peak2.bed";
		else 
			print $0 > "./"fname"/tmp_peak2_exclude.bed"}' &
	#shift NOT extend, ideally would also bounds check since can result in negative coord
	zcat $4 | sed 's/\s$//g' | awk -v fname="$fname" -v shift1=$6 'BEGIN {OFS="\t"}
		{if ($1~/chr/ && $1 !="chrM" && $6=="+" && $1 !~/random/ && $3>$2  && $2>0 && $3>0)
			print $1,$2+shift1,$3+shift1 > "./"fname"/tmp_read1.bed";
		else if ($1~/chr/ && $1 !="chrM" && $6=="-" && $1 !~/random/ && $3>$2  && $2>shift1 && $3>shift1)
			print $1,$2-shift1,$3-shift1 > "./"fname"/tmp_read1.bed";
		else 
			print $0 > "./"fname"/tmp_read1_exclude.bed"}' &
	zcat $5 | sed 's/\s$//g' | awk -v fname="$fname" -v shift2=$7 'BEGIN {OFS="\t"}
		{if ($1~/chr/ && $1 !="chrM" && $6=="+" && $1 !~/random/ && $3>$2  && $2>0 && $3>0)
			print $1,$2+shift2,$3+shift2 > "./"fname"/tmp_read2.bed";
		else if ($1~/chr/ && $1 !="chrM" && $6=="-" && $1 !~/random/ && $3>$2  && $2>shift2 && $3>shift2)
			print $1,$2-shift2,$3-shift2 > "./"fname"/tmp_read2.bed";
		else 
			print $0 > "./"fname"/tmp_read2_exclude.bed"}' &
    wait
	echo "wc -l tmp_read1.bed `wc -l ./$fname/tmp_read1.bed` | cut -f 1 -d " "" >> ./$fname/${1}.log
	echo "wc -l tmp_read2.bed `wc -l ./$fname/tmp_read2.bed` | cut -f 1 -d " "" >> ./$fname/${1}.log
	echo "wc -l tmp_read1_exclude.bed `wc -l ./$fname/tmp_read1_exclude.bed` | cut -f 1 -d " "" >> ./$fname/${1}.log
	echo "wc -l tmp_read2_exclude.bed `wc -l ./$fname/tmp_read2_exclude.bed` | cut -f 1 -d " "" >> ./$fname/${1}.log
elif [ $# -eq 11 ]
then 
	#mostly copied from above
	#add functionality to combine replicates
	#"proper" way would be to downscale replicates so that both contribute equally (not done)
	fname=$1 
	mkdir $fname 
	echo $@ > ./$fname/${1}.log
    echo `date` >> ./$fname/${1}.log
	echo "StepI: clean input"
    #newlines do not work in echo by default, must use echo -e which is not POSIX
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
    # output combined reads to make sure further steps don't fail, a bit less efficient this way though
	cat ./$fname/tmp_read1a.bed ./$fname/tmp_read1b.bed > ./$fname/tmp_read1.bed
	cat ./$fname/tmp_read2a.bed ./$fname/tmp_read2b.bed > ./$fname/tmp_read2.bed

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
#the sort is in non-ASCII order ie 1, 10, 11, 2, etc
intersectBed -a ./$fname/tmp_peak1.bed -b ./$fname/tmp_peak2.bed -u | sort -k1,1 -k2,2n -k3,3n > ./$fname/tmp_common_peak1.bed &
intersectBed -a ./$fname/tmp_peak2.bed -b ./$fname/tmp_peak1.bed -u | sort -k1,1 -k2,2n -k3,3n > ./$fname/tmp_common_peak2.bed &
intersectBed -a ./$fname/tmp_peak1.bed -b ./$fname/tmp_peak2.bed -v | sort -k1,1 -k2,2n -k3,3n > ./$fname/tmp_unique_peak1.bed &
intersectBed -a ./$fname/tmp_peak2.bed -b ./$fname/tmp_peak1.bed -v | sort -k1,1 -k2,2n -k3,3n > ./$fname/tmp_unique_peak2.bed &
wait

# WARNING!!! must sort before merging, above sort cmd is faster but this is more clear
cat ./$fname/tmp_common_peak1.bed ./$fname/tmp_common_peak2.bed | sortBed | mergeBed > ./$fname/tmp_common_peak.bed

echo "wc -l tmp_common_peak1.bed `wc -l ./$fname/tmp_common_peak1.bed | cut -f 1 -d " "`" >> ./$fname/${1}.log
echo "wc -l tmp_common_peak2.bed `wc -l ./$fname/tmp_common_peak2.bed | cut -f 1 -d " "`" >> ./$fname/${1}.log
echo "wc -l tmp_unique_peak1.bed `wc -l ./$fname/tmp_unique_peak1.bed | cut -f 1 -d " "`" >> ./$fname/${1}.log
echo "wc -l tmp_unique_peak2.bed `wc -l ./$fname/tmp_unique_peak2.bed | cut -f 1 -d " "`" >> ./$fname/${1}.log
echo "wc -l tmp_common_peak.bed `wc -l ./$fname/tmp_common_peak.bed | cut -f 1 -d " "`" >> ./$fname/${1}.log

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

coverageBed -a ./$fname/tmp_read1.bed -b ./$fname/tmp_unique_peak1.bed | cut -f 1-4 > ./$fname/tmp_unique_peak1_count_read1 &
coverageBed -a ./$fname/tmp_read2.bed -b ./$fname/tmp_unique_peak1.bed | cut -f 1-4 > ./$fname/tmp_unique_peak1_count_read2 &
coverageBed -a ./$fname/tmp_read1.bed -b ./$fname/tmp_common_peak1.bed | cut -f 1-4 > ./$fname/tmp_common_peak1_count_read1 &
coverageBed -a ./$fname/tmp_read2.bed -b ./$fname/tmp_common_peak1.bed | cut -f 1-4 > ./$fname/tmp_common_peak1_count_read2 &
coverageBed -a ./$fname/tmp_read1.bed -b ./$fname/tmp_common_peak2.bed | cut -f 1-4 > ./$fname/tmp_common_peak2_count_read1 &
coverageBed -a ./$fname/tmp_read2.bed -b ./$fname/tmp_common_peak2.bed | cut -f 1-4 > ./$fname/tmp_common_peak2_count_read2 &
coverageBed -a ./$fname/tmp_read1.bed -b ./$fname/tmp_unique_peak2.bed | cut -f 1-4 > ./$fname/tmp_unique_peak2_count_read1 &
coverageBed -a ./$fname/tmp_read2.bed -b ./$fname/tmp_unique_peak2.bed | cut -f 1-4 > ./$fname/tmp_unique_peak2_count_read2 &
wait
if [ $# -eq 11 ]
then
	coverageBed -a ./$fname/tmp_read1a.bed -b ./$fname/tmp_unique_peak1.bed | cut -f 1-4 > ./$fname/tmp_unique_peak1_count_read1a &
	coverageBed -a ./$fname/tmp_read2a.bed -b ./$fname/tmp_unique_peak1.bed | cut -f 1-4 > ./$fname/tmp_unique_peak1_count_read2a &
	coverageBed -a ./$fname/tmp_read1a.bed -b ./$fname/tmp_common_peak1.bed | cut -f 1-4 > ./$fname/tmp_common_peak1_count_read1a &
	coverageBed -a ./$fname/tmp_read2a.bed -b ./$fname/tmp_common_peak1.bed | cut -f 1-4 > ./$fname/tmp_common_peak1_count_read2a &
	coverageBed -a ./$fname/tmp_read1a.bed -b ./$fname/tmp_common_peak2.bed | cut -f 1-4 > ./$fname/tmp_common_peak2_count_read1a &
	coverageBed -a ./$fname/tmp_read2a.bed -b ./$fname/tmp_common_peak2.bed | cut -f 1-4 > ./$fname/tmp_common_peak2_count_read2a &
	coverageBed -a ./$fname/tmp_read1a.bed -b ./$fname/tmp_unique_peak2.bed | cut -f 1-4 > ./$fname/tmp_unique_peak2_count_read1a &
	coverageBed -a ./$fname/tmp_read2a.bed -b ./$fname/tmp_unique_peak2.bed | cut -f 1-4 > ./$fname/tmp_unique_peak2_count_read2a &	
	coverageBed -a ./$fname/tmp_read1b.bed -b ./$fname/tmp_unique_peak1.bed | cut -f 1-4 > ./$fname/tmp_unique_peak1_count_read1b &
	coverageBed -a ./$fname/tmp_read2b.bed -b ./$fname/tmp_unique_peak1.bed | cut -f 1-4 > ./$fname/tmp_unique_peak1_count_read2b &
	coverageBed -a ./$fname/tmp_read1b.bed -b ./$fname/tmp_common_peak1.bed | cut -f 1-4 > ./$fname/tmp_common_peak1_count_read1b &
	coverageBed -a ./$fname/tmp_read2b.bed -b ./$fname/tmp_common_peak1.bed | cut -f 1-4 > ./$fname/tmp_common_peak1_count_read2b &
	coverageBed -a ./$fname/tmp_read1b.bed -b ./$fname/tmp_common_peak2.bed | cut -f 1-4 > ./$fname/tmp_common_peak2_count_read1b &
	coverageBed -a ./$fname/tmp_read2b.bed -b ./$fname/tmp_common_peak2.bed | cut -f 1-4 > ./$fname/tmp_common_peak2_count_read2b &
	coverageBed -a ./$fname/tmp_read1b.bed -b ./$fname/tmp_unique_peak2.bed | cut -f 1-4 > ./$fname/tmp_unique_peak2_count_read1b &
	coverageBed -a ./$fname/tmp_read2b.bed -b ./$fname/tmp_unique_peak2.bed | cut -f 1-4 > ./$fname/tmp_unique_peak2_count_read2b &	
	wait
	#all 4 col included, last col not cut
	cat ./$fname/tmp_common_peak1_count_read1a ./$fname/tmp_common_peak2_count_read1a > ./$fname/tmp_common_peak_count_read1a &
	cat ./$fname/tmp_common_peak1_count_read2a ./$fname/tmp_common_peak2_count_read2a > ./$fname/tmp_common_peak_count_read2a &
	cat ./$fname/tmp_unique_peak1_count_read1a ./$fname/tmp_common_peak1_count_read1a ./$fname/tmp_common_peak2_count_read1a ./$fname/tmp_unique_peak2_count_read1a > ./$fname/tmp_peak_count_read1a &
	cat ./$fname/tmp_unique_peak1_count_read2a ./$fname/tmp_common_peak1_count_read2a ./$fname/tmp_common_peak2_count_read2a ./$fname/tmp_unique_peak2_count_read2a > ./$fname/tmp_peak_count_read2a &
	cat ./$fname/tmp_common_peak1_count_read1b ./$fname/tmp_common_peak2_count_read1b > ./$fname/tmp_common_peak_count_read1b &
	cat ./$fname/tmp_common_peak1_count_read2b ./$fname/tmp_common_peak2_count_read2b > ./$fname/tmp_common_peak_count_read2b &
	cat ./$fname/tmp_unique_peak1_count_read1b ./$fname/tmp_common_peak1_count_read1b ./$fname/tmp_common_peak2_count_read1b ./$fname/tmp_unique_peak2_count_read1b > ./$fname/tmp_peak_count_read1b &
	cat ./$fname/tmp_unique_peak1_count_read2b ./$fname/tmp_common_peak1_count_read2b ./$fname/tmp_common_peak2_count_read2b ./$fname/tmp_unique_peak2_count_read2b > ./$fname/tmp_peak_count_read2b &	
	#wait
fi

#not sure how to do this w/out awk
awk 'BEGIN {OFS="\t"} {print $1,$2,$3,"unique_peak1"}' ./$fname/tmp_unique_peak1_count_read1 >>  ./$fname/tmp_MAnorm.bed
awk 'BEGIN {OFS="\t"} {print $1,$2,$3,"common_peak1"}' ./$fname/tmp_common_peak1_count_read1 >>  ./$fname/tmp_MAnorm.bed
awk 'BEGIN {OFS="\t"} {print $1,$2,$3,"common_peak2"}' ./$fname/tmp_common_peak2_count_read1 >>  ./$fname/tmp_MAnorm.bed
awk 'BEGIN {OFS="\t"} {print $1,$2,$3,"unique_peak2"}' ./$fname/tmp_unique_peak2_count_read1 >>  ./$fname/tmp_MAnorm.bed

cat ./$fname/tmp_common_peak1_count_read1 ./$fname/tmp_common_peak2_count_read1 | cut -f 4 > ./$fname/tmp_common_peak_count_read1 &
cat ./$fname/tmp_common_peak1_count_read2 ./$fname/tmp_common_peak2_count_read2 | cut -f 4 > ./$fname/tmp_common_peak_count_read2 &
cat ./$fname/tmp_unique_peak1_count_read1 ./$fname/tmp_common_peak1_count_read1 ./$fname/tmp_common_peak2_count_read1 ./$fname/tmp_unique_peak2_count_read1 | cut -f 4 > ./$fname/tmp_peak_count_read1 &
cat ./$fname/tmp_unique_peak1_count_read2 ./$fname/tmp_common_peak1_count_read2 ./$fname/tmp_common_peak2_count_read2 ./$fname/tmp_unique_peak2_count_read2 | cut -f 4 > ./$fname/tmp_peak_count_read2 &
wait

if [ -f ./$fname/tmp_MAnorm_merge.bed ];
then
	rm ./$fname/tmp_MAnorm_merge.bed
fi

coverageBed -a ./$fname/tmp_read1.bed -b ./$fname/tmp_common_peak.bed | cut -f 1-4 > ./$fname/tmp_merge_common_read1 &
coverageBed -a ./$fname/tmp_read2.bed -b ./$fname/tmp_common_peak.bed | cut -f 1-4 > ./$fname/tmp_merge_common_read2 &
if [ $# -eq 11 ]
then
	coverageBed -a ./$fname/tmp_read1a.bed -b ./$fname/tmp_common_peak.bed | cut -f 1-4 > ./$fname/tmp_merge_common_read1a &
	coverageBed -a ./$fname/tmp_read2a.bed -b ./$fname/tmp_common_peak.bed | cut -f 1-4 > ./$fname/tmp_merge_common_read2a &
	coverageBed -a ./$fname/tmp_read1b.bed -b ./$fname/tmp_common_peak.bed | cut -f 1-4 > ./$fname/tmp_merge_common_read1b &
	coverageBed -a ./$fname/tmp_read2b.bed -b ./$fname/tmp_common_peak.bed | cut -f 1-4 > ./$fname/tmp_merge_common_read2b &
	wait
	cat ./$fname/tmp_unique_peak1_count_read1a ./$fname/tmp_merge_common_read1a ./$fname/tmp_unique_peak2_count_read1a > ./$fname/tmp_merge_common_peak_count_read1a &
	cat ./$fname/tmp_unique_peak1_count_read2a ./$fname/tmp_merge_common_read2a ./$fname/tmp_unique_peak2_count_read2a > ./$fname/tmp_merge_common_peak_count_read2a &
	cat ./$fname/tmp_unique_peak1_count_read1b ./$fname/tmp_merge_common_read1b ./$fname/tmp_unique_peak2_count_read1b > ./$fname/tmp_merge_common_peak_count_read1b &
	cat ./$fname/tmp_unique_peak1_count_read2b ./$fname/tmp_merge_common_read2b ./$fname/tmp_unique_peak2_count_read2b > ./$fname/tmp_merge_common_peak_count_read2b &
fi
wait

cat ./$fname/tmp_unique_peak1_count_read1 | awk -v fname="$fname" 'BEGIN {OFS="\t"} {print $1,$2,$3,"tmp_unique_peak1" >> "./"fname"/tmp_MAnorm_merge.bed"}'
cat ./$fname/tmp_merge_common_read1 | awk 'BEGIN {OFS="\t"} {print $1,$2,$3,"merged_common_peak"}' >> ./$fname/tmp_MAnorm_merge.bed
cat ./$fname/tmp_unique_peak2_count_read1 | awk -v fname="$fname" 'BEGIN {OFS="\t"} {print $1,$2,$3,"tmp_unique_peak2" >> "./"fname"/tmp_MAnorm_merge.bed"}'

cat ./$fname/tmp_unique_peak1_count_read1 ./$fname/tmp_merge_common_read1 ./$fname/tmp_unique_peak2_count_read1 | cut -f 4 > ./$fname/tmp_merge_common_peak_count_read1
cat ./$fname/tmp_unique_peak1_count_read2 ./$fname/tmp_merge_common_read2 ./$fname/tmp_unique_peak2_count_read2 | cut -f 4 > ./$fname/tmp_merge_common_peak_count_read2

#####################################################################################################
echo "SetpIV: normalize using common peaks"
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
    return 1
fi
