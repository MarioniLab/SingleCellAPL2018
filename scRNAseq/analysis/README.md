# scRNAseq analyses

After QC, analyses were run using the following scripts in the following order:

CD8_timecourse_analyses_20170830.Rmd
CD8_APLseparated_analyses_20170904.Rmd
CD8_APLcombined_analyses_20170904.Rmd
CD8_APLcombined_unstimulated_noise_20180226.Rmd

The diffusion pseudotime analyses in these scripts were originally performed in 2016/2017 using the dpt R package from Haghverdi et al, 2016, Nat Methods, 13:845. DPT analysis has since been incorporated into the destiny Bioconductor package, but for running the scripts above, it is essential to download and install the original dpt package.  