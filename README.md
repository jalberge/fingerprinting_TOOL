# Fingerprinting tool

This leverages GATK's `ExtractFingerprints` and `CrosscheckFingerprints` on a terra sample_set

## Current version

v0.0.1

Not working / not tested: 

- upload to bucket
- sync to terra
- sample_set instead of table of samples
- GDC/NCI bams
- LocalizeToDisk (only stream or local for now)
- scRNA (probably change haplotypeDB)
- hg19 / non human
- Heatmap and tsv export of results

I copied the code within the cloud task to allow `rm` minibam after extractfingerprints.

## Performance

Tested on:

- 150 WGS available in Google cloud without localization
- 933 SNPs from haplotype_db (default)
- *<15 min*

## Default values:

```
stream_bam_or_cram=True
haplotype_db="gs://getzlab-workflows-reference_files-oa/hg38/Homo_sapiens_assembly38.haplotype_database.txt"
ref_fa="gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.fasta"
ref_fai="gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.fasta.fai"
ref_dict="gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.dict"
```
