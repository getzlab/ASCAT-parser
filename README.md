# ASCAT Parser

This repo contains the necessary functions to parse [ASCATNGS](https://pubmed.ncbi.nlm.nih.gov/27930809/) ([repo](https://github.com/cancerit/ascatNgs)) outputs to an allelic capseg equivalent to run ABSOLUTE.


# Requirements
- pandas
- scipy
- numpy
- gcsfs - For pandas reading from google buckets

# Test Data
- [Tumor caveman](gs://fc-0520679c-c3bc-46ab-9120-6d18617a28f3/TCGA-2G-AAEQ-01A-11D-A734-36.copynumber.caveman.csv): gs://fc-0520679c-c3bc-46ab-9120-6d18617a28f3/TCGA-2G-AAEQ-01A-11D-A734-36.copynumber.caveman.csv

- [Tumor copy number](gs://fc-0520679c-c3bc-46ab-9120-6d18617a28f3/TCGA-2G-AAEQ-01A-11D-A734-36.copynumber.txt.gz): gs://fc-0520679c-c3bc-46ab-9120-6d18617a28f3/TCGA-2G-AAEQ-01A-11D-A734-36.copynumber.txt.gz

- [Normal SNPs](gs://fc-0520679c-c3bc-46ab-9120-6d18617a28f3/TCGA-2G-AAEQ-10A-01D-A734-36.count):  gs://fc-0520679c-c3bc-46ab-9120-6d18617a28f3/TCGA-2G-AAEQ-10A-01D-A734-36.count

- [Interval list](gs://fc-0520679c-c3bc-46ab-9120-6d18617a28f3/whole_exome_illumina_coding_v1.Homo_sapiens_assembly19.baits.interval_list_reduced.csv): gs://fc-0520679c-c3bc-46ab-9120-6d18617a28f3/whole_exome_illumina_coding_v1.Homo_sapiens_assembly19.baits.interval_list_reduced.csv

## To Do

- [ ] Add requirements.txt
- [ ] Format repo.
- [ ] Better way to use test data.
