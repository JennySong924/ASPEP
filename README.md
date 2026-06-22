# Alternative Splicing Junctions PSI Estimation Pipeline （ASPEP）


## Overview

`ASPEP` is a R-based command line tool for calculating percent splided in(PSI) of alternative splicing junctions.

<!--![](https://github.com/JennySong924/ASPEP/blob/main/static/Pipeline-01.png)-->



## Input files

> Junction score file generaged by RegTools

The pipeline needs an input file generaged by RegTools, which includes junction coordinates, strand and scores information. Here is an example of junction score file:

|chrom|start|end|name|score|strand|splice_site|acceptors_skipped|exons_skipped|donors_skipped|anchor|known_donor|known_acceptor|known_junction|gene_names|gene_ids|transcripts|
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
|chr1|4774516|4777525|JUNC00000020|29|-|GT-AG|0|0|0|DA|1|1|1|Mrpl15|Mrpl15|NM_001177658.1,NM_025300.4,NR_033530.1|

To run the pipeline, columns of `chrom`,`start`,`end`,`strand` are required.

## Basic Usage

```bash
bash psi_calculation.sh -f [input file] -c 5 -o [output dir] -n [prefix]
```

Example:

```bash
bash psi_calculation.sh -f OHC.junc.annot.example.refseq -c 5 -o ./example -n OHC
```

### Parameters

**Required**
- `-f <input file>` File path of the junction score file
- `-c <coverage cutoff>` Junction score cutoff. Only keep junctions with score >= coverage cutoff in order to reduce noise.                 
- `-o <output dir>` Name of the output directory to create
- `-n <prefix>` Prefix of output files


## Installation & Setup

### Prerequisites
- **R**: Version 4.5.0 or higher
- bedtools: Installation struction https://bedtools.readthedocs.io/en/latest/content/installation.html

### Installation

1. **Clone the ASPEP Repository**
```bash
git clone [https://github.com/XX](https://github.com/JennySong924/ASPEP.git)
```

2. **Install R Dependencies**
```bash
Rscript r_package_preparation.r
```

## Key Output Files

- `<prefix>.psi.gz` gzipped file containing following columns:

| Column | Discription | Example |
|--------|-------------|--------|
|chr| Junction chromosome|chr10|
|start|	Junction start coordinate|60304528|
|end| Junction end coordinate|60304825|
|strand| Junction strand|-|
|gene_names| Names of gene that junction lies in|Cdh23|
|gene_ids| IDs of gene that junction lies in (usually the same as gene_names)|Cdh23|
|junction_id| ID of junction in form of chr:start-end:strand|chr10:60304528-60304825:-|	
|start_coverage	| Depth of junction start calculated based on junction file|73|
|end_coverage| Depth of junction end calculated based on junction file|73|	
|score| Junction score extracted from input file |	73|
|start_psi| PSI on junction start| 1|	
|end_psi| PSI on junction end | 1 |	
|start_psi_weight|	PSI weight on junction start calculated based on depth of start coverage and end coverage | 0.5|
|end_psi_weight	| PSI weight on junction end calculated based on depth of start coverage and end coverage | 0.5 |
|start_psi_weighted	| start_psi_weight * start_psi| 0.5|
|end_psi_weighted| end_psi_weight * end_psi | 0.5 |
|weighted_psi| Final psi of the junction : start_psi_weight + end_psi_weighted| 1 |


## Version History
* June 19, 2026: Initial release.


## Citation

If you use `ASPEP`, please cite XXX

<!-- Miao, J., Song, G., Wu, Y., Hu, J., Wu, Y., Basu, S., Andrews, J. S., Schaumberg, K., Fletcher, J. M., Schmitz, L. L., & Lu, Q. (2022). [Reimagining Gene-Environment Interaction Analysis for Human Complex Traits](https://doi.org/10.1101/2022.12.11.519973 ). bioRxiv, 2022.2012.2011.519973.  -->

## Contact

For questions and comments, please open a GitHub issue (preferred) or contact Jie Song at jsong89@sjtu.edu.cn.
