# johannesnicolaus/nanoplasty-nf: Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## v0.1.0dev - [unreleased]

Initial release of johannesnicolaus/nanoplasty-nf, created with the [nf-core](https://nf-co.re/) template.

### `Added`

- NANOPlasty assembly of `mito`, `chloro` and `plant_mito` genomes from ONT or PacBio reads, or use of an existing assembly
- Read QC with NanoPlot and SeqKit, optional Chopper filtering
- Assembly QC: reads mapped back to the wrapped circular assembly (minimap2), samtools stats, mosdepth, reads spanning the start/end, QUAST
- Annotation with MITOS2 (`mito`), Chloë (`chloro`) or a miniprot + BLASTN + tRNAscan-SE homology annotator (`plant_mito`, optional for `chloro`)
- Circular maps with pyCirclize (genes by class, read depth, GC content, GC skew)
- Per-sample organelle summary table in MultiQC

### `Fixed`

### `Dependencies`

### `Deprecated`
