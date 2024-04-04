# Fingerprinting tool

wolF tasks to fingerprint many samples and compare them in a pairwise manner. This leverages GATK's `ExtractFingerprints` and `CrosscheckFingerprints` on a terra sample_set, similar to the `CrossCheckLaneFingerprints` [task](https://github.com/getzlab/picard_TOOL/blob/6ca4f895f0ddd72c9bdcaebc706b1e4d4a2d8967/wolF/tasks.py#L27-L65).

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

## Default values

```
stream_bam_or_cram=True
haplotype_db="gs://getzlab-workflows-reference_files-oa/hg38/Homo_sapiens_assembly38.haplotype_database.txt"
ref_fa="gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.fasta"
ref_fai="gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.fasta.fai"
ref_dict="gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.dict"
```

## Execution (single pass)

```
import wolf
import dalmatian
import pandas as pd
import numpy as np

def terra_na(x):
    y = ( pd.isna(x) ) | ( x=="")
    return y

fingerprinting_workflows = wolf.ImportTask("git@github.com:jalberge/fingerprinting_TOOL.git", commit="2c3d24fb", branch="master")

# extract samples from a terra workspace
WORKSPACE = "<project/workspace>"
WIC = wolf.fc.WorkspaceInputConnector(WORKSPACE)
S = WIC.samples

# subset table
S = S.loc[ ~terra_na(S['hg38_analysis_ready_bam']) ]

with wolf.Workflow(workflow=fingerprinting_workflows.fingerprint) as w:
    w.run(RUN_NAME="fingerprint_test", # fill in run name
        bam_or_cram=S["hg38_analysis_ready_bam"].tolist(),
        bai_or_crai=S["hg38_analysis_ready_bam_index"].tolist(),
        sample_id=S.index.values.tolist(),
        sample_set_id="fingerprint_test"
        )
```


