# Comparing differential transcription factor chip-seq methods

This repository holds some of my scripts and findings for comparing differential chip-seq methods published in [frontiers in genetics 2015](https://www.ncbi.nlm.nih.gov/pubmed/25972895).

* [DIME](DIME/) -- Some matlab code that I ported over to R. Method took too long to run so I gave up
* [MAnorm](MAnorm/) -- Rewrote most of this program so it runs faster, requires replicates, and uses edgeR
* [edgeR-DiffBind-DBChIP](edgeR-DiffBind-DBChIP/) -- I go over how to use these three methods and document some differences between them
