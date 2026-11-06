# johannesnicolaus/nanoplasty-nf

[![GitHub Actions CI Status](https://github.com/johannesnicolaus/nanoplasty-nf/actions/workflows/nf-test.yml/badge.svg)](https://github.com/johannesnicolaus/nanoplasty-nf/actions/workflows/nf-test.yml)
[![GitHub Actions Linting Status](https://github.com/johannesnicolaus/nanoplasty-nf/actions/workflows/linting.yml/badge.svg)](https://github.com/johannesnicolaus/nanoplasty-nf/actions/workflows/linting.yml)[![Cite with Zenodo](http://img.shields.io/badge/DOI-10.5281/zenodo.XXXXXXX-1073c8?labelColor=000000)](https://doi.org/10.5281/zenodo.XXXXXXX)
[![nf-test](https://img.shields.io/badge/unit_tests-nf--test-337ab7.svg)](https://www.nf-test.com)

[![Nextflow](https://img.shields.io/badge/version-%E2%89%A525.10.4-green?style=flat&logo=nextflow&logoColor=white&color=%230DC09D&link=https%3A%2F%2Fnextflow.io)](https://www.nextflow.io/)
[![nf-core template version](https://img.shields.io/badge/nf--core_template-4.2.1-green?style=flat&logo=nfcore&logoColor=white&color=%2324B064&link=https%3A%2F%2Fnf-co.re)](https://github.com/nf-core/tools/releases/tag/4.2.1)
[![run with conda](http://img.shields.io/badge/run%20with-conda-3EB049?labelColor=000000&logo=anaconda)](https://docs.conda.io/en/latest/)
[![run with docker](https://img.shields.io/badge/run%20with-docker-0db7ed?labelColor=000000&logo=docker)](https://www.docker.com/)
[![run with singularity](https://img.shields.io/badge/run%20with-singularity-1d355c.svg?labelColor=000000)](https://sylabs.io/docs/)

## Introduction

**johannesnicolaus/nanoplasty-nf** assembles circular organelle genomes (animal mitochondria, plastids and
plant mitochondria) from Oxford Nanopore or PacBio long reads with
[NANOPlasty](https://github.com/ndierckx/NANOPlasty), then checks, annotates and draws them.

1. Read QC ([`NanoPlot`](https://github.com/wdecoster/NanoPlot), [`SeqKit`](https://bioinf.shenwei.me/seqkit/)); optional filtering ([`Chopper`](https://github.com/wdecoster/chopper))
2. Assembly ([`NANOPlasty`](https://github.com/ndierckx/NANOPlasty)), skipped for samples that already have an assembly
3. Assembly QC: reads mapped back to the circular assembly ([`minimap2`](https://github.com/lh3/minimap2)), coverage ([`mosdepth`](https://github.com/brentp/mosdepth), [`samtools stats`](https://www.htslib.org/)), reads spanning the start/end of the circle, assembly statistics ([`QUAST`](https://github.com/ablab/quast))
4. Annotation, by assembly type:
   - `mito` (animal mitochondria): [`MITOS2`](https://gitlab.com/Bernt/MITOS)
   - `chloro` (plastids): [`Chloë`](https://github.com/ian-small/Chloe.jl), or the homology annotator with `--chloro_annotator homology`
   - `plant_mito` (plant mitochondria): homology annotator — [`miniprot`](https://github.com/lh3/miniprot) against RefSeq organelle proteins, [`BLASTN`](https://blast.ncbi.nlm.nih.gov/) against RefSeq rRNAs, and [`tRNAscan-SE`](https://github.com/UCSC-LoweLab/tRNAscan-SE) in organellar mode
5. Circular map of each assembly ([`pyCirclize`](https://github.com/moshi4/pyCirclize)): genes coloured by functional class, read depth, GC content and GC skew
6. Report ([`MultiQC`](http://multiqc.info/)), including a per-sample organelle summary table

## Usage

> [!NOTE]
> If you are new to Nextflow and nf-core, please refer to [this page](https://nf-co.re/docs/get_started/environment_setup/overview) on how to set-up Nextflow. Make sure to [test your setup](https://nf-co.re/docs/get_started/run-your-first-pipeline) with `-profile test` before running the workflow on actual data.

Prepare a samplesheet with one row per sample:

`samplesheet.csv`:

```csv
sample,reads,platform,type,genome_range,seed,r10,min_read_length,chloro_reference,assembly
fly_mt,fly.fastq.gz,ont,mito,15000-22000,,true,,,
ath_cp,ath.fastq.gz,ont,chloro,140000-170000,,,,,
ath_mt,ath_hifi.fastq.gz,pacbio_hifi,plant_mito,300000-400000,,,,ath_cp.fasta,
```

Only `sample`, `reads`, `platform` and `type` are required. See the [usage documentation](docs/usage.md) for every column.

Then run:

```bash
nextflow run johannesnicolaus/nanoplasty-nf \
   -profile <docker/singularity/.../institute> \
   --input samplesheet.csv \
   --outdir <OUTDIR>
```

> [!WARNING]
> Please provide pipeline parameters via the CLI or Nextflow `-params-file` option. Custom config files including those provided by the `-c` Nextflow option can be used to provide any configuration _**except for parameters**_; see [docs](https://nf-co.re/docs/running/run-pipelines#using-parameter-files).

> [!NOTE]
> NANOPlasty runs from `ghcr.io/johannesnicolaus/nanoplasty:0.2`, built by the
> [NANOPlasty repository](https://github.com/johannesnicolaus/NANOPlasty). To run a local checkout instead, pass
> `--nanoplasty_script /path/to/NANOPlasty0.2.pl` with its dependencies (BLAST+, kalign3, Perl MCE and
> Parallel::ForkManager) from `-profile conda` or the host `PATH`. Samples with an `assembly` column skip NANOPlasty.

## Pipeline output

See the [output documentation](docs/output.md).

## Credits

johannesnicolaus/nanoplasty-nf was originally written by Johannes Nicolaus Wibisana.

We thank the following people for their extensive assistance in the development of this pipeline:

- Nicolas Dierckxsens, author of NANOPlasty

## Contributions and Support

If you would like to contribute to this pipeline, please see the [contributing guidelines](docs/CONTRIBUTING.md).

## Citations

<!-- TODO nf-core: Add citation for pipeline after first release. Uncomment lines below and update Zenodo doi and badge at the top of this file. -->
<!-- If you use johannesnicolaus/nanoplasty-nf for your analysis, please cite it using the following doi: [10.5281/zenodo.XXXXXX](https://doi.org/10.5281/zenodo.XXXXXX) -->



An extensive list of references for the tools used by the pipeline can be found in the [`CITATIONS.md`](CITATIONS.md) file.

This pipeline uses code and infrastructure developed and maintained by the [nf-core](https://nf-co.re) community, reused here under the [MIT license](https://github.com/nf-core/tools/blob/main/LICENSE).

> **The nf-core framework for community-curated bioinformatics pipelines.**
>
> Philip Ewels, Alexander Peltzer, Sven Fillinger, Harshil Patel, Johannes Alneberg, Andreas Wilm, Maxime Ulysse Garcia, Paolo Di Tommaso & Sven Nahnsen.
>
> _Nat Biotechnol._ 2020 Feb 13. doi: [10.1038/s41587-020-0439-x](https://dx.doi.org/10.1038/s41587-020-0439-x).
