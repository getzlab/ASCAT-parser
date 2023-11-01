# ASCAT Parser

This repo contains the necessary functions to parse [ASCATNGS](https://pubmed.ncbi.nlm.nih.gov/27930809/) ([repo](https://github.com/cancerit/ascatNgs)) outputs to an allelic capseg equivalent to run ABSOLUTE.

# ASCAT Parser

Some changes were made relative to ascatNgs:
Standard ascatNGS VAF fits were suboptimal. The signature of poor het VAF fitting appears as "=" (equal signs) across ABSOLUTE allelic copy ratio profiles. Occasional "=" aCR levels can be real, but if most the genome has a small allele shift, then there is a problem. 
ASCAT-parser minimizes false shifts away from VAF=0.5 using a symmetric peaked beta distribution acround 0.5. This is done for each segment with at least one het.



# Requirements
- pandas
- scipy
- numpy
- gcsfs - For pandas reading from google buckets

# Inputs: 
see inline comments in ascatparser.py (lines 10-28) for argument information

# Output file columns (AllelicCapSeg format for input to ABSOLUTE):
- Chromosome: chromsome label for segment
- Start.bp: start coordinate for segment
- End.bp: end coordinate for segment
- n_probes: number of targeted regions used for total copy number fit
- length: genomic length of segment (End.bp-Start.bp)
- n_hets: Number of het in segment used to measure allelic imbalance. 
- f: mean het allele fraction shift (0.5 corresponds to no allelic imbalance)
- tau: total copy ratio
- sigma.tau: 1-sigma statistical error on tau
- mu.minor: allelic copy ratio (aCR) at lower level
- sigma.minor: 1-sigma statistical error on mu.minor
- mu.major: allelic copy ratio (aCR) at higher level
- sigma.major: 1-sigma statistical error on mu.major
- SegLabelCNLOH: 0 or 1 indicating false vs true copy-neutral LOH

