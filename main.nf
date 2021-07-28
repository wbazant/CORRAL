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
  file("${sample}_R1.fastq.gz"), file("${sample}_R2.fastq.gz")

  script:
  """
  ${params.wgetCommand} $fastqUrlR1 -O ${sample}_R1.fastq.gz
  ${params.wgetCommand} $fastqUrlR2 -O ${sample}_R2.fastq.gz
  """
}


process downloadSingleSra {

  label 'download'

  input:
  tuple val(sample), val(stringWithRunAccession)

  output:
  file("${sample}.fastq")

  script:
  """
  getFastqFromSraSingle $stringWithRunAccession "${sample}.fastq
  """
}

process downloadPairedSra {

  label 'download'

  input:
  tuple val(sample), val(stringWithRunAccession)

  output:
  file("${sample}_R1.fastq"), file("${sample}_R2.fastq")

  script:
  """
  getFastqFromSraPaired $stringWithRunAccession ${sample}_R1.fastq ${sample}_R2.fastq
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

process summarizeAlignments {
  publishDir 'summarizedAlignments'

  label 'summarize'

  input:
  path(numReadsPath)
  path(alignmentsSam)

  output:
  path("summarizedAlignments.tsv")

  script:
  """
  summarize_marker_alignments \
    --input ${alignmentsSam} \
    --refdb-marker-to-taxon-id-path ${params.marker_to_taxon_id_path} \
    --refdb-format eukprot \
    --output-type taxon_all \
    --num-reads \$(cat ${numReadsPath}) \
    --output summarizedAlignments.tsv 
  """
}

process filterPostSummarize {
  publishDir 'filteredAlignmentSummaries'

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
  publishDir params.resultDir, mode: 'move'  

  input:
  file("*.taxon-to-cpm.tsv")

  output:
  file("cpms.tsv")

  script:
  """
  makeTsv.pl . .taxon-to-cpm.tsv > cpms.tsv
  """
}

def postAlign(sampleId, numReadsPath, alignmentsSam) {
  sas = summarizeAlignments(numReadsPath, alignmentsSam)
  fps = filterPostSummarize(sampleId, sas)
  return fps
}

def singleWget(input) {
  reads = downloadSingleWget(input)
  alignments = bowtie2Single(reads)
  numReads = countReads(reads)
  return postAlign(input.map{it[0]}, numReads, alignments)
}


workflow runSingleWget {

  main:
  input = Channel.fromPath(params.inputPath).splitCsv(sep: "\t")
  xs = singleWget(input)
  makeTsv(xs.collect())
}
