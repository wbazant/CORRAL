# Marker alignments in Nextflow

## Installation
This workflow is not containerised, but the dependencies are quite minimal:
- `bowtie2`
- [Marker alignments package](https://github.com/wbazant/marker_alignments) and its tool `marker_alignments`.

Additionally, `samtools stats` is the default and recommended for alignment stats.

By default, `bowtie2`, `marker_alignments` and `samtools` are assumed to be on `$PATH` but you can provide a path to an executable in the pipeline config.


If you want to use `--downloadMethod wget` you also need `wget`. If you want to use `--downloadMethod sra` you need the SRA EUtils, with `prefetch` and `fastq-dump` on `$PATH`.

`--unpackMethod bz2` requires `bzip2` on `$PATH`.

You also need a `bowtie2` reference database of taxonomic markers, like ChocoPhlAn or EukDetect.


## Usage

Main parameters:

| param         | value type        | description  |
| ------------- | ------------- | ------------ |
| inputPath  | path to file | TSV: sample ID,fastq URL or run ID, [second URL for paired reads] |
| downloadMethod | "wget" / "sra" | |
| libraryLayout | "single" / "paired" | |
| resultDir  | path to dir  | publish directory |
| refdb | path pattern | bowtie2 -x parameter |
| bowtie2Command | shell | Run bowtie2 |
| alignmentStatsCommand | shell | `samtools stats` or other |
| summarizeAlignmentsCommand | shell | path to `marker_alignments` optionally with filter arguments to use|

Optional parameters:

| param         | value type        | description  |
| ------------- | ------------- | ------------ |
| marker_to_taxon_path | path to file | summarize_marker_alignments --marker_to_taxon_path parameter |
| unpackMethod | "bz2" | for FTP .tar.bz2 content |

### Example 
I run it like that:

#### `run.sh`
```
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"
REF_PATH="~/eukprot"

nextflow pull wbazant/marker-alignments-nextflow -r main

nextflow run wbazant/marker-alignments-nextflow -r main \
  --inputPath $DIR/in.tsv  \
  --resultDir $DIR/results \
  --downloadMethod wget \
  --unpackMethod bz2 \
  --libraryLayout paired \
  --refdb ${REF_PATH}/ncbi_eukprot_met_arch_markers.fna \
  --marker_to_taxon_path ${REF_PATH}/busco_taxid_link.txt  \
  -c $DIR/cluster.conf \
  -with-trace -resume | tee $DIR/tee.out

```

#### `cluster.conf`

```  
process {
  executor = 'lsf'
  maxForks = 60
  
  withLabel: 'download' {
    maxForks = 5
    maxRetries = 3
  }
  withLabel: 'align' {
    errorStrategy = 'finish'
  }
}
```

#### `in.tsv`
```
SRS011061       https://downloads.hmpdacc.org/dacc/hhs/genome/microbiome/wgs/analysis/hmwgsqc/v1/SRS011061.tar.bz2
SRS011086       https://downloads.hmpdacc.org/dacc/hhs/genome/microbiome/wgs/analysis/hmwgsqc/v1/SRS011086.tar.bz2
```
