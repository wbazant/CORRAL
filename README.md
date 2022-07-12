# Marker alignments in Nextflow

### Get Started
  * Install Nextflow
    
    `curl https://get.nextflow.io | bash`

## Usage

Main parameters:

| param         | value type        | description  |
| ------------- | ------------- | ------------ |
| inputPath  | path to file | TSV: sample ID,fastq URL or run ID, [second URL for paired reads] |
| downloadMethod | "wget" / "sra" / "local" | |
| libraryLayout | "single" / "paired" | |
| resultDir  | path to dir  | publish directory |
| refdb | path pattern | bowtie2 -x parameter |
| bowtie2Command | shell | Run bowtie2 |
| alignmentStatsCommand | shell | `samtools stats` or other |
| summarizeAlignmentsCommand | shell | path to `marker_alignments` optionally with filter arguments to use|
| apiKey | string | ncbi ApiKey |

Optional parameters:

| param         | value type        | description  |
| ------------- | ------------- | ------------ |
| marker_to_taxon_path | path to file | summarize_marker_alignments --marker_to_taxon_path parameter |
| unpackMethod | "bz2" | for FTP .tar.bz2 content |

### Example 

 `nextflow run VEuPathDB/blastSimilarity -with-trace -c  <config_file> -r main`
