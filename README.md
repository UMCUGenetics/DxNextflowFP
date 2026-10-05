# UMCUGenetics/dxnextflowfp

[![GitHub Actions CI Status](https://github.com/UMCUGenetics/dxnextflowfp/actions/workflows/nf-test.yml/badge.svg)](https://github.com/UMCUGenetics/dxnextflowfp/actions/workflows/nf-test.yml)
[![GitHub Actions Linting Status](https://github.com/UMCUGenetics/dxnextflowfp/actions/workflows/linting.yml/badge.svg)](https://github.com/UMCUGenetics/dxnextflowfp/actions/workflows/linting.yml)
[![nf-test](https://img.shields.io/badge/unit_tests-nf--test-337ab7.svg)](https://www.nf-test.com)
[![Nextflow](https://img.shields.io/badge/version-%E2%89%A525.10.4-green?style=flat&logo=nextflow&logoColor=white&color=%230DC09D&link=https%3A%2F%2Fnextflow.io)](https://www.nextflow.io/)
[![nf-core template version](https://img.shields.io/badge/nf--core_template-4.1.0-green?style=flat&logo=nfcore&logoColor=white&color=%2324B064&link=https%3A%2F%2Fnf-co.re)](https://github.com/nf-core/tools/releases/tag/4.1.0)
[![run with docker](https://img.shields.io/badge/run%20with-docker-0db7ed?labelColor=000000&logo=docker)](https://www.docker.com/)
[![run with singularity](https://img.shields.io/badge/run%20with-singularity-1d355c.svg?labelColor=000000)](https://sylabs.io/docs/)


**UMCUGenetics/dxnextflowfp** 

Genome Diagnotics Nextflow Fingerprint. This workflow is used to perform alignment on the fingerprint data. These data are further processed by BAM_VARIANTCALLING_INTERVALS to create a vcf file for fingerprinting to detect swaps/contamination in samples.

# Usage
## HPC
On the HPC you can submit a pipeline run using the `run_dxnextflowfp.sh` script. This submits the pipeline in an sbatch and runs nextflow with a slurm profile.
```bash
./DxNextflowFP/run_dxnextflowfp.sh \
  --input fastq_dir/ \
  --outdir <path>   \
  --email <address> \
  [options]
```

## Alternatively
The pipeline can also be executed using a regular nextflow run command:
```bash
nextflow run UMCUGenetics/dxnextflowfp \
   -profile <docker/singularity/.../institute> \
   --input fastq_directory/ \
   --outdir <OUTDIR> \
   --email <address> \
   [options]
```

> [!WARNING]
> Please provide pipeline parameters via the CLI or Nextflow `-params-file` option. Custom config files including those provided by the `-c` Nextflow option can be used to provide any configuration _**except for parameters**_; see [docs](https://nf-co.re/docs/running/run-pipelines#using-parameter-files).

# Citations

An extensive list of references for the tools used by the pipeline can be found in the [`CITATIONS.md`](CITATIONS.md) file.

This pipeline uses code and infrastructure developed and maintained by the [nf-core](https://nf-co.re) community, reused here under the [MIT license](https://github.com/nf-core/tools/blob/main/LICENSE).

