import wolf

# for cloud usage we require samtools
GATK_DOCKER = "gcr.io/broad-getzlab-workflows/gatk4_wolf:v6"


class CloudExtractFingerprints(wolf.Task):
    name = "cloud_extract_fingerprints"
    inputs = {
        "sample_id": None,
        "bam_or_cram": None,  # set to None for required arguments
        "bai_or_crai": "",
        "haplotype_db": None,
        "ref_fa": None,
        "ref_fai": None,
        "ref_dict": None,
        "gatk_path": "/gatk/gatk"
    }
    # by default, stream sequence data with gatk
    overrides = {
        "bam_or_cram": "string",
        "bai_or_crai": "string"
    }

    script = """
set -euxo pipefail

# samtools requires auth token
export GCS_OAUTH_TOKEN=$(gcloud auth application-default print-access-token)
    
# transform haplotype db to bed intervals
grep -E -v "^[@#]" "${haplotype_db}" | awk -v OFS='\t' '{print $1,$2-1,$2}' > bed

# create mini bam
samtools view -h -f 1 -F 1024 -o tmp.bam -q 30 -M -L bed -X "${bam}" "${bai}"
samtools index tmp.bam

# extract fingerprints

# we could directly use CrosscheckFingerprints, but this way we store the vcf for later use.
${gatk_path} --java-options "-Xms3G" \
  ExtractFingerprint \
  INPUT=tmp.bam \
  OUTPUT=${sample_id}.fingerprinted.vcf \
  HAPLOTYPE_MAP=${haplotype_db} \
  REFERENCE_SEQUENCE=${ref_fa}
  
rm bed
rm tmp.bam

    """
    outputs = {
        "fingerprint_vcf": "*.fingerprinted.vcf"
    }
    docker = GATK_DOCKER
    resources = {"mem": "4G"}


class ExtractFingerprints(wolf.Task):
    name = "extract_fingerprints"
    inputs = {
        "sample_id": None,
        "bam_or_cram": None,  # set to None for required arguments
        "bai_or_crai": "",
        "haplotype_db": None,
        "ref_fa": None,
        "ref_fai": None,
        "ref_dict": None,
        "gatk_path": "/gatk/gatk"
    }
    # by default, stream sequence data with gatk
    overrides = {
        "bam_or_cram": "string",
        "bai_or_crai": "string"
    }

    script = """
    export GCS_OAUTH_TOKEN=$(gcloud auth application-default print-access-token)
    
    set -euxo pipefail
    # we could directly use CrosscheckFingerprints, but this way we store the vcf for later use.
    ${gatk_path} --java-options "-Xms3G" \
      ExtractFingerprint \
      INPUT=${bam_or_cram} \
      OUTPUT=${sample_id}.fingerprinted.vcf \
      HAPLOTYPE_MAP=${haplotype_db} \
      REFERENCE_SEQUENCE=${ref_fa}
    """
    outputs = {
        "fingerprint_vcf": "*.fingerprinted.vcf"
    }
    docker = GATK_DOCKER
    resources = {"mem": "4G"}


class CrosscheckFingerprints(wolf.Task):
    name = "crosscheck_fingerprints"
    inputs = {
        "input_vcfs": None,
        "sample_set_id": None,
        "haplotype_db": None,
        "gatk_path": "/gatk/gatk"
    }
    script = """
    # does gatk accept lists of files? doc says yes https://gatk.broadinstitute.org/hc/en-us/articles/360037594711-CrosscheckFingerprints-Picard#--INPUT
    ${gatk_path} --java-options "-Xms3G" \
        CrosscheckFingerprints \
        I=${input_vcfs} \
        HAPLOTYPE_MAP=~{haplotype_db} \
        LOD_THRESHOLD=-5 \
        CROSSCHECK_MODE=CHECK_ALL_OTHERS \
        CROSSCHECK_BY=FILE \
        OUTPUT="${sample_set_id}.metrics" \
        MO="${sample_set_id}.matrix" \
        EXIT_CODE_WHEN_MISMATCH=0
    """
    outputs = {
        "crosscheck_metrics": "*.metrics",
        "crosscheck_matrix": "*.matrix"
    }
    docker = GATK_DOCKER
    resources = {"mem": "4G"}


def fingerprint(
        bam_or_cram=None,
        bai_or_crai=None,
        sample_id=None,
        sample_set_id=None,
        stream_bam_or_cram=True,
        haplotype_db="gs://getzlab-workflows-reference_files-oa/hg38/Homo_sapiens_assembly38.haplotype_database.txt",
        ref_fa="gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.fasta",
        ref_fai="gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.fasta.fai",
        ref_dict="gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.dict",
        gatk_path="gatk",
        workspace=None,
        bucket=None):
    ref_files = wolf.LocalizeToDisk(files={
        "haplotype_db": haplotype_db,
        "ref_fa": ref_fa,
        "ref_fai": ref_fai,
        "ref_dict": ref_dict
    })

    if stream_bam_or_cram:
        fingerprints = CloudExtractFingerprints(inputs={
            "bam_or_cram": bam_or_cram,
            "bai_or_crai": bai_or_crai,
            "sample_id": sample_id,
            "haplotype_db": ref_files["haplotype_db"],
            "ref_fa": ref_files["ref_fa"],
            "ref_fai": ref_files["ref_fai"],
            "ref_dict": ref_files["ref_dict"],
            "gatk_path": gatk_path
        })
    else:
        fingerprints = ExtractFingerprints(inputs={
            "bam_or_cram": bam_or_cram,
            "bai_or_crai": bai_or_crai,
            "sample_id": sample_id,
            "stream_bam_or_cram": stream_bam_or_cram,
            "haplotype_db": ref_files["haplotype_db"],
            "ref_fa": ref_files["ref_fa"],
            "ref_fai": ref_files["ref_fai"],
            "ref_dict": ref_files["ref_dict"],
            "gatk_path": gatk_path
        })

    crosschecked_fingerprints = CrosscheckFingerprints(inputs={
        "input_vcfs": [fingerprints["fingerprint_vcf"]],
        "sample_set_id": sample_set_id,
        "haplotype_db": ref_files["haplotype_db"],
        "gatk_path": gatk_path
    })

    if bucket is not None:
        upload_vcf = wolf.UploadToBucket(
            files=[
                fingerprints["fingerprint_vcf"],
                crosschecked_fingerprints["crosscheck_matrix"],
                crosschecked_fingerprints["crosscheck_metrics"]
            ],
            bucket=bucket
        )

    if workspace:
        upload_dict = {
            'crosschecked_fingerprints_matrix': crosschecked_fingerprints["crosscheck_matrix"],
            'crosschecked_fingerprints_metrics': crosschecked_fingerprints["crosscheck_metrics"]
        }
        sync_run = wolf.SyncToWorkspace(
            nameworkspace=workspace,
            entity_type="sample_set",
            entity_name=sample_set_id,
            attr_map=upload_dict
        )
