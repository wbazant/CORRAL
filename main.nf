nextflow.enable.dsl=2


process downloadSingleWget {
  label 'download'
  input:
  tuple val(sample), val(fastqUrl)

  output:
  file("${sample}.fastq.gz")

  script:
  """
  ${params.wgetCommand} $fastqUrl -O ${sample}.fastq.gz
  """
}

process downloadPairedWget {
  label 'download'
  input:
  tuple val(sample), val(fastqUrlR1), val(fastqUrlR2)

  output:
  tuple file("${sample}_R1.fastq.gz"), file("${sample}_R2.fastq.gz")

  script:
  """
  ${params.wgetCommand} $fastqUrlR1 -O ${sample}_R1.fastq.gz
  ${params.wgetCommand} $fastqUrlR2 -O ${sample}_R2.fastq.gz
  """
}


process downloadSingleSra {

  label 'download'

  input:
  tuple val(sample), val(runAccession)

  output:
  file("${sample}.fastq")

  script:
  """
  getFastqFromSraSingle $runAccession ${sample}.fastq
  """
}

process downloadPairedSra {

  label 'download'

  input:
  tuple val(sample), val(runAccession)

  output:
  tuple file("${sample}_R1.fastq"), file("${sample}_R2.fastq")

  script:
  """
  getFastqFromSraPaired $runAccession ${sample}_R1.fastq ${sample}_R2.fastq
  """
}
process countReads {
  label 'count'
  input:
  path(readsFastq)

  output:
  path("numReads.txt")
  

  script:
  """
  grep -c '^@' ${readsFastq} > numReads.txt
  """
}

process bowtie2Single {
  label 'align'
  input:
  path(readsFastq)

  output:
  path("alignmentsSingle.sam")

  script:
  """
  bowtie2 --omit-sec-seq --no-discordant --no-unal \
    -x ${params.refdb} \
    -U ${readsFastq} \
    -S alignmentsSingle.sam 
  """

}
process bowtie2Paired {
  label 'align'
  input:
  tuple path(readsFastqR1), path(readsFastqR2)

  output:
  path("alignmentsPaired.sam")

  script:
  """
  bowtie2 --omit-sec-seq --no-discordant --no-unal \
    -x ${params.refdb} \
    -1 ${readsFastqR1} \
    -2 ${readsFastqR2} \
    -S alignmentsPaired.sam
  """
}

process summarizeAlignmentsIntoMarkers {
  publishDir "${params.resultDir}/markers"

  label 'postAlign'

  input:
  val(sample)
  path(numReadsPath)
  path(alignmentsSam)

  output:
  path("${sample}.markers.tsv")

  script:
  """
  summarize_marker_alignments \
    --input ${alignmentsSam} \
    --refdb-marker-to-taxon-id-path ${params.marker_to_taxon_id_path} \
    --refdb-format eukprot \
    --output-type marker_all \
    --num-reads \$(cat ${numReadsPath}) \
    --output ${sample}.markers.tsv 
  """
}
process summarizeAlignmentsIntoTaxa {
  publishDir "${params.resultDir}/taxa"
  label 'postAlign'

  input:
  val(sample)
  path(numReadsPath)
  path(alignmentsSam)

  output:
  path("${sample}.taxa.tsv")

  script:
  """
  summarize_marker_alignments \
    --input ${alignmentsSam} \
    --refdb-marker-to-taxon-id-path ${params.marker_to_taxon_id_path} \
    --refdb-format eukprot \
    --output-type taxon_all \
    --num-reads \$(cat ${numReadsPath}) \
    --output ${sample}.taxa.tsv 
  """
}

process filterPostSummarize {
  publishDir "${params.resultDir}/taxon-to-cpm"
  label 'postAlign'

  input:
  val(sample)
  path(summarizedAlignmentsPath)

  output:
  path("${sample}.taxon-to-cpm.tsv")

  script:
  """
  filterAlignmentSummaries.pl ${sample} ${summarizedAlignmentsPath} \
    > ${sample}.taxon-to-cpm.tsv
  """
}

process makeTsv {
  publishDir params.resultDir, mode: 'move', overwrite: true  
  label 'postAlign'

  input:
  file("*.taxon-to-cpm.tsv")

  output:
  file("cpms.tsv")

  script:
  """
  makeTsv.pl . .taxon-to-cpm.tsv > cpms.tsv
  """
}

def postAlign(sample, numReadsPath, alignmentsSam) {
  summarizeAlignmentsIntoMarkers(sample, numReadsPath, alignmentsSam)
  sas = summarizeAlignmentsIntoTaxa(sample, numReadsPath, alignmentsSam)
  fps = filterPostSummarize(sample, sas)
  return fps
}

def singleWget(input) {
  reads = downloadSingleWget(input)
  alignments = bowtie2Single(reads)
  numReads = countReads(reads)
  return postAlign(input.map{it[0]}, numReads, alignments)
}

def pairedWget(input) {
  reads = downloadPairedWget(input)
  alignments = bowtie2Paired(reads)
  numReads = countReads(reads.map{it[0]})
  return postAlign(input.map{it[0]}, numReads, alignments)
}

def singleSra(input) {
  reads = downloadSingleSra(input)
  alignments = bowtie2Single(reads)
  numReads = countReads(reads)
  return postAlign(input.map{it[0]}, numReads, alignments)
}

def pairedSra(input) {
  reads = downloadPairedSra(input)
  alignments = bowtie2Paired(reads)
  numReads = countReads(reads.map{it[0]})
  return postAlign(input.map{it[0]}, numReads, alignments)
}

workflow {
  input = Channel.fromPath(params.inputPath).splitCsv(sep: "\t")
  if (params.downloadMethod == 'sra') {
    input = input.map{it.size() == 1 ? [it[0], it[0]] : it}
  }

  if(params.downloadMethod == 'wget' && params.libraryLayout == 'single'){
    xs = singleWget(input)
  } else if(params.downloadMethod == 'wget' && params.libraryLayout == 'paired'){
    xs = pairedWget(input)
  } else if(params.downloadMethod == 'sra' && params.libraryLayout == 'single'){
    xs = singleSra(input)
  } else if(params.downloadMethod == 'sra' && params.libraryLayout == 'paired'){
    xs = pairedSra(input)
  }

  makeTsv(xs.collect())
}
