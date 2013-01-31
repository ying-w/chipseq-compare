MAnorm {WIP}
============
My modified/improved version for testing.

Original can be found [here](http://bcb.dfci.harvard.edu/~gcyuan/MAnorm/R_tutorial.html). 
Published in [Genome Biology 2012](http://www.ncbi.nlm.nih.gov/pubmed/22424423)

Improvements
------------
* Allow for reading/combining 2 replicates (since ENCODE samples have 2 replicates)
** done simply by concatinating the files, proper way would be to downsample so they have equal weight
* option to run bedtools run in parallel (faster if you have enough IO)
* allow to specify filename for directory
