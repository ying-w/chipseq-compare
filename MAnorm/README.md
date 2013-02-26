MAnorm {WIP}
============
My modified/improved version for testing.

Original can be found [here](http://bcb.dfci.harvard.edu/~gcyuan/MAnorm/R_tutorial.html) 
Published in [Genome Biology 2012](http://www.ncbi.nlm.nih.gov/pubmed/22424423)

Improvements
------------
* Allow for reading/combining 2 replicates (since ENCODE samples have 2 replicates)
** concatinating the read files, proper way would be to downsample so replicates have equal weight
* run bedtools run in parallel (faster if you have enough IO)
* allow to specification directory filename
* edgeR peak calling (for replicates)

Notes
-----
* There is something wrong with how the pvalues are calculated (see v2 MAnorm.r for details). 
* In v3 I have removed this feature and use edgeR instead (will not work if no replicates). 
* I run bedtools in parallel in the background for significant speedup.
* sortBed before mergeBed command has been added to v3 so the merged dataset is actually merged

Experience from running on a quad core workstation with 32gb of ram: 
* StepI is I/O bound (mostly reading from disk and writing back to disk)
* StepII & StepIII is CPU bound (lots of parallel bedtools)
* Ram usage I've seen with ENCODE datasets is 20gb+/-5gb


v2 will retain the same input/output as the original script (from supplement of paper) but will be faster
v3 will have different input/output options but the same methods will be applied
