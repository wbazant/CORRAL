nextflow.enable.dsl=2


process downloadSingleWget {
  label 'download'
  input:
  tuple val(sample), val(fastqUrl)

  output:
  tuple val(sample), file("${sample}.fastq.gz")

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
  tuple val(sample), file("${sample}_R1.fastq.gz"), file("${sample}_R2.fastq.gz")

  script:
  """
  ${params.wgetCommand} $fastqUrlR1 -O ${sample}_R1.fastq.gz
  ${params.wgetCommand} $fastqUrlR2 -O ${sample}_R2.fastq.gz
  """
}

process downloadPairedWgetUnpackBz2 {
  label 'download'
  input:
  tuple val(sample), val(url)

  output:
  tuple val(sample), file("${sample}_R1.fastq"), file("${sample}_R2.fastq")

  script:
  """
  ${params.wgetCommand} $url -O ${sample}.tar.bz2
  tar -xvjf ${sample}.tar.bz2
  fastq_R1=\$(find ${sample} -name '*1.fastq')
  fastq_R2=\$(find ${sample} -name '*2.fastq')
  mv -v "\$fastq_R1" "${sample}_R1.fastq"
  mv -v "\$fastq_R2" "${sample}_R2.fastq"
  """
}


process downloadSingleSra {

  label 'download'

  input:
  tuple val(sample), val(runAccession)

  output:
  tuple val(sample), file("${sample}.fastq")

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
  tuple val(sample), file("${sample}_R1.fastq"), file("${sample}_R2.fastq")

  script:
  """
  getFastqFromSraPaired $runAccession ${sample}_R1.fastq ${sample}_R2.fastq
  """
}

process bowtie2Single {
  label 'align'
  input:
  tuple val(sample), path(readsFastq)

  output:
  tuple val(sample), path("numReads.txt"), path("alignmentsSingle.sam")

  script:
  """
  grep -c '^@' ${readsFastq} > numReads.txt

  bowtie2 --omit-sec-seq --no-discordant --no-unal \
    -x ${params.refdb} \
    -U ${readsFastq} \
    -S alignmentsSingle.sam 
  """

}
process bowtie2Paired {
  label 'align'
  input:
  tuple val(sample), path(readsFastqR1), path(readsFastqR2)

  output:
  tuple val(sample), path("numReads.txt"), path("alignmentsPaired.sam")

  script:
  """
  grep -c '^@' ${readsFastqR1} > numReads.txt

  bowtie2 --omit-sec-seq --no-discordant --no-unal \
    -x ${params.refdb} \
    -1 ${readsFastqR1} \
    -2 ${readsFastqR2} \
    -S alignmentsPaired.sam
  """
}
process alignmentStats {
  publishDir "${params.resultDir}/alignmentStats"
  label 'samtools'

  input:
  tuple val(sample), path(numReadsPath), path(alignmentsSam)

  output:
  tuple val(sample), path("${sample}.alignmentStats.txt")

  script:
  """
  samtools stats ${alignmentsSam} > ${sample}.alignmentStats.txt
  """
}

process filterAlignments {
  label 'postAlign'

  input:
  tuple val(sample), path(numReadsPath), path(alignmentsSam)

  output:
  tuple val(sample), path(numReadsPath), path("${sample}.filteredAlignments.sam")

  script:
  """
  ${params.samtoolsFilterCommand} ${alignmentsSam} -o ${sample}.filteredAlignments.sam
  """
}

process summarizeAlignmentsIntoMarkers {
  publishDir "${params.resultDir}/markers"

  label 'postAlign'

  input:
  tuple val(sample), path(numReadsPath), path(alignmentsSam)

  output:
  tuple val(sample), path("${sample}.markers.tsv")

  script:
  """
  summarize_marker_alignments \
    --input ${alignmentsSam} \
    --refdb-marker-to-taxon-path ${params.marker_to_taxon_path} \
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
  tuple val(sample), path(numReadsPath), path(alignmentsSam)

  output:
  tuple val(sample), path("${sample}.taxa.tsv")

  script:
  """
  summarize_marker_alignments \
    --input ${alignmentsSam} \
    --refdb-marker-to-taxon-path ${params.marker_to_taxon_path} \
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
  tuple val(sample), path(summarizedAlignmentsPath)

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

def postAlign(sample_numReadsPath_alignmentsSam) {
  alignmentStats(sample_numReadsPath_alignmentsSam)
  sample_numReadsPath_filteredAlignmentsSam = filterAlignments(sample_numReadsPath_alignmentsSam)
  summarizeAlignmentsIntoMarkers(sample_numReadsPath_filteredAlignmentsSam)
  sample_sas = summarizeAlignmentsIntoTaxa(sample_numReadsPath_filteredAlignmentsSam)
  fps = filterPostSummarize(sample_sas)
  return fps
}

def singleWget(input) {
  sample_reads = downloadSingleWget(input)
  sample_numReads_alignments = bowtie2Single(sample_reads)
  return postAlign(sample_numReads_alignments)
}

def pairedWget(input) {
  sample_reads = downloadPairedWget(input)
  sample_numReads_alignments = bowtie2Paired(sample_reads)
  return postAlign(sample_numReads_alignments)
}

def pairedWgetUnpackBz2(input) {
  sample_reads = downloadPairedWgetUnpackBz2(input)
  sample_numReads_alignments = bowtie2Paired(sample_reads)
  return postAlign(sample_numReads_alignments)
}

def singleSra(input) {
  sample_reads = downloadSingleSra(input)
  sample_numReads_alignments = bowtie2Single(sample_reads)
  return postAlign(sample_numReads_alignments)
}

def pairedSra(input) {
  sample_reads = downloadPairedSra(input)
  sample_numReads_alignments = bowtie2Paired(sample_reads)
  return postAlign(sample_numReads_alignments)
}

workflow {
  input = Channel.fromPath(params.inputPath).splitCsv(sep: "\t")
  if (params.downloadMethod == 'sra') {
    input = input.map{it.size() == 1 ? [it[0], it[0]] : it}
  }

  if(params.downloadMethod == 'wget' && params.libraryLayout == 'single'){
    xs = singleWget(input)
  } else if(params.downloadMethod == 'wget' && params.libraryLayout == 'paired'){
    if(params.unpackMethod == 'bz2'){
      xs = pairedWgetUnpackBz2(input)
    } else {
      xs = pairedWget(input)
    }
  } else if(params.downloadMethod == 'sra' && params.libraryLayout == 'single'){
    xs = singleSra(input)
  } else if(params.downloadMethod == 'sra' && params.libraryLayout == 'paired'){
    xs = pairedSra(input)
  }

  makeTsv(xs.collect())
}
