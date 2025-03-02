# Arnoult Lab dual CRISPR screen for DNA damage/repair pathways
This repository accumulates data processing scripts, intermediate results, and final results related to Nebenfuehr et al. (2025). Sequencing data for this experiment can be found at GEO accession number [GSE290153](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE290153).

Contents of the directory are organized as follows:
- 01_raw_data_processing contains Slurm and Python scripts used for evaluating, trimming, and mapping raw fastq files, then sorting guides and barcodes to generate files that contain read counts for each guide pair in each sample. Software versions are listed in relevant Slurm scripts, and all Python scripts used version 3.6.3. All basal screen biological samples were sequenced five times and thus have five separate associated technical replicate fastq files, while irradiated samples were only sequenced once. All fastq files were processed separately through the initial steps and initial sorting, then guide pair counts were combined within each biological replicate.

- 02_count_files contains the accumulated guide pair count files for the basal screen samples and irradiated samples (IR). For further analysis, IR reps 1 and 2 were compared to basal screen Day 0 samples reps 1 and 2, respectively.

- 03_count_processing contains custom Python scripts to process guide pair counts into phenotype value (log2(fold changes between timepoints)). Replicate 3 was dropped during this analysis due to low read depth and widely disparate count distributions compared to other samples. All scripts in this step were run using python version 3.10.12 with numpy v.1.21.2 and pandas v.1.1.5.

- 04_phenotype_results contains results output from step 4.
