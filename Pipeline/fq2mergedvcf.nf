#!/usr/bin/env nextflow


// ------------------------- Reference & Software Parameters -------------------------

// Human genome fasta reference parameters
def rDir = file("/projects/pe41/Resources/resources_broad_hg38")
params.refFasta    = rDir.resolve("resources_broad_hg38_v0_Homo_sapiens_assembly38.fasta")
params.refFastaIdx = rDir.resolve("resources_broad_hg38_v0_Homo_sapiens_assembly38.fasta.fai")
params.refFastaDct = rDir.resolve("resources_broad_hg38_v0_Homo_sapiens_assembly38.dict")

// Analysis region parameters
// Pipeline designed to process data from multiple panels simultaneously,
//     maps used to dynamically dispatch panel intervals based on run/plate name
// Panel 3 CHIP-only panel directory and files
def panel03Dir = file("/projects/pe41/Resources/Archimedes/Phase_2_CHIP-only/3470281_shared")
def panel03map = [
  "panelBed"             : panel03Dir.resolve("3470281_Covered.bed"),
  "panelIntervals"       : panel03Dir.resolve("3470281_Covered.interval_list"),
  "panelBedPadded"       : panel03Dir.resolve("3470281_Covered_padded_100bp.bed"),
  "panelIntervalsPadded" : panel03Dir.resolve("3470281_Covered_padded_100bp.interval_list"),
]
// Mapping run names to associated panel files
params.runpanelmap = [
  "ARCH_run04_plate1": panel03map,
  "ARCH_run04_plate2": panel03map,
]

// Agilent AGeNT tools location
def agentDir = file("/projects/pe41/Software/AGeNT_3.0.6/agent/lib")
params.trimmerJar = agentDir.resolve("trimmer-3.0.5.jar")
params.creakJar   = agentDir.resolve("creak-1.0.5.jar")

// Conda env paths
params.condabwagatk   = "/projects/pe41/Software/miniconda3/envs/bwagatk4"    // bwa=0.7.17-r1188, gatk4=4.4.0.0-0
params.condavardict   = "/projects/pe41/Software/miniconda3/envs/vardictjava" // vardictjava=1.8.3, bcftools=1.18
params.condabcftools  = "/projects/pe41/Software/miniconda3/envs/bcftools"    // bcftools=1.18
params.condaagentjava = "/projects/pe41/Software/miniconda3/envs/openjdk-21"  // openjdk=21.0.1, AGeNT=3.0.6
params.condasamtools  = "/projects/pe41/Software/miniconda3/envs/samtools"    // samtools=1.18
params.condafastqc    = "/projects/pe41/Software/miniconda3/envs/fastqc"      // fastqc=0.12.1
params.condamultiqc   = "/projects/pe41/Software/miniconda3/envs/multiqc"     // multiqc=1.18
// Standalone binaries
params.multiqcbin     = "/projects/pe41/Software/miniconda3/envs/multiqc/bin/multiqc"

// Output directory for published workflow outputs
params.wfOutputDir = file("./pipeline_outputs")


// ------------------------- Pipeline Processes -------------------------

process trimFastqs {
  conda "${params.condaagentjava}"
  memory "12 GB"
  time "16m"
  input:
    tuple val(basename), path(fqpair), val(groupname)
  output:
    tuple val(basename), path("*_trim_R*.fastq.gz"), val(groupname)
  script:
  """
  java -Xmx${task.memory.getGiga()}g -jar ${params.trimmerJar} \
        -v2 \
        -out ./${basename}_trim \
        -fq1 ${fqpair[0]} \
        -fq2 ${fqpair[1]}
  """
}

process alignBwaSortSam {
  label "gatkshort"
  conda "${params.condabwagatk}"
  time "24m"
  input:
    tuple val(basename), path(fqpair), val(groupname)
  output:
    tuple val(basename), path("${basename}.*"), val(groupname)
  """
  set -o pipefail
  bwa mem -K 100000000 -v 3 -C -t ${task.cpus} -Y \
      -R "@RG\\tID:${basename}\\tSM:${basename}\\tLB:lib1\\tPL:ILLUMINA\\tPU:unit1" \
      ${params.refFasta} ${fqpair[0]} ${fqpair[1]} | \
  gatk SortSam \
      --INPUT /dev/stdin \
      --OUTPUT ${basename}.bam \
      --SORT_ORDER "coordinate" \
      --CREATE_INDEX true \
      --CREATE_MD5_FILE true \
      --MAX_RECORDS_IN_RAM 300000
  """
}

process dedupCreak {
  conda "${params.condaagentjava}"
  publishDir "${params.wfOutputDir}/${groupname}/bams", mode: "link"
  memory "12 GB"
  input:
    tuple val(basename), path(bamset), val(groupname)
  output:
    tuple val(basename), path("${basename}.creak.*"), val(groupname)
  script:
  def panelBed = params.runpanelmap[groupname]["panelBed"]
  """
  java -Xmx${task.memory.getGiga()}g -jar ${params.creakJar} \
      -c HYBRID -d 0 -f -mm 25 -mr 30 -F -MS 1 -MD 2 \
      -b ${panelBed} -o ${basename}.creak.bam \
      "${bamset[1]}"
  """
}

process depthSamtools {
  conda "${params.condasamtools}"
  publishDir "${params.wfOutputDir}/${groupname}/depths", mode: "link"
  cpus 1
  input:
    tuple val(basename), path(bamset), val(groupname)
  output:
    tuple val(basename), path("${basename}.depths.tsv.gz"), val(groupname)
  script:
  def panelBedPadded = params.runpanelmap[groupname]["panelBedPadded"]
  """
  samtools depth -a --threads ${task.cpus} -H \
      -b ${panelBedPadded} \
      "${bamset[1]}" | gzip > "${basename}.depths.tsv.gz"
  """
}

process collectHsMetrics {
  label "gatkshort"
  conda "${params.condabwagatk}"
  time "8m"
  input:
    tuple val(basename), path(bamset), val(groupname)
  output:
    tuple val(basename), path("${basename}.hs_metrics.txt"), val(groupname)
  script:
  def panelIntervals = params.runpanelmap[groupname]["panelIntervals"]
  """
  gatk CollectHsMetrics \
      -I "${bamset[1]}" \
      -O "${basename}.hs_metrics.txt" \
      -R "${params.refFasta}" \
      -BAIT_INTERVALS "${panelIntervals}" \
      -TARGET_INTERVALS "${panelIntervals}"
  """
}

process callVarDict {
  label "gatkshort"
  conda "${params.condavardict}"
  input:
    tuple val(basename), path(bamset), val(groupname)
  output:
    tuple val(basename), path("${basename}.vardict.vcf.gz*"), val(groupname)
  script:
  def panelBedPadded = params.runpanelmap[groupname]["panelBedPadded"]
  def af_thr = "0.005"
  """
  vardict-java \
      --nosv \
      -th 4 \
      -G "${params.refFasta}" \
      -f ${af_thr} \
      -N "${basename}" \
      -b "${bamset[1]}" \
      -c 1 -S 2 -E 3 -g 4 "${panelBedPadded}" | \
  teststrandbias.R | \
  var2vcf_valid.pl \
      -N "${basename}" \
      -E \
      -f ${af_thr} | bgzip > "${basename}.vardict.vcf.gz"
  bcftools index --threads ${task.cpus} --tbi "${basename}.vardict.vcf.gz"
  """
}

process decomposeVarDict {
  conda "${params.condabcftools}"
  publishDir "${params.wfOutputDir}/${groupname}/vcfs", mode: "link"
  time "5m"
  input:
    tuple val(basename), path(vcfset), val(groupname)
  output:
    tuple val(basename), path("${basename}.vardict.norm.vcf*"), val(groupname)
  """
  bcftools norm -a -m -any --threads ${task.cpus} \
      -f "${params.refFasta}" \
      -o "${basename}.vardict.norm.vcf.gz" \
      "${vcfset[0]}"
  bcftools index --threads ${task.cpus} --tbi "${basename}.vardict.norm.vcf.gz"
  """
}

process mergeVarDict {
  conda "${params.condabcftools}"
  publishDir "${params.wfOutputDir}/${groupname}", mode: "link"
  input:
    tuple val(basenames), val(vcfs), val(groupname)
  output:
    path "${groupname}_vardict_merged.vcf.gz*"
  script:
  vcfstr = vcfs.collect({ it[0] }).join(" ")
  """
  bcftools merge --threads ${task.cpus} -m none -o ${groupname}_vardict_merged.vcf.gz ${vcfstr}
  bcftools index --threads ${task.cpus} --tbi ${groupname}_vardict_merged.vcf.gz
  """
}

process metricsFastqc {
  conda "${params.condafastqc}"
  publishDir "${params.wfOutputDir}/${groupname}", mode: "link"
  cpus 16
  memory "16 GB"
  time "28m"
  input:
    tuple val(basenames), val(allfqpairs), val(groupname)
  output:
    tuple path("${groupname}_fqc_reports.tar.gz"), path("${groupname}_fastqc_multiqc_report*")
  script:
  """
  mkdir ${groupname}_fqc_reports
  fastqc -t ${task.cpus} -o ${groupname}_fqc_reports ${allfqpairs.flatten().join(" ")}
  tar czvf ${groupname}_fqc_reports.tar.gz ${groupname}_fqc_reports
  ${params.multiqcbin} -n ${groupname}_fastqc_multiqc_report.html -o ./ ${groupname}_fqc_reports/
  """
}

process multiqcHsmetrics {
  conda "${params.condamultiqc}"
  publishDir "${params.wfOutputDir}/${groupname}", mode: "link"
  cpus 1
  memory "4 GB"
  time "8m"
  input:
    tuple val(basenames), path(hsmetrics), val(groupname)
  output:
    path "${groupname}_hsmetrics_multiqc_report*"
  script:
  """
  multiqc -n ${groupname}_hsmetrics_multiqc_report.html -o ./ ./*.hs_metrics.txt
  """
}


// ------------------------- Main Workflow -------------------------

workflow {
  // Channel Setup
  // Get list of all fastq subdirectories
  def fqDirs = file("./fastqdirs/*", type: "dir").sort()

  // Build list of fastq channels for each directory
  // Assumes paired fastqs with pairs designated by _R1/2_ in filename
  // Results in the structure: val(basename), path(fqpair), val(groupname)
  def fqChannels = fqDirs.collect {d ->
    Channel.fromFilePairs(d.resolve("*_R{1,2}_*.fastq.gz")).merge(Channel.value(d.getName()))
  }
  // Concat individual group channels
  def allFqPairs = null
  for (fqc in fqChannels) {
    if (!allFqPairs) {
      allFqPairs = fqc
    } else {
      allFqPairs = allFqPairs.concat(fqc)
    }
  }

  // Fastqc metrics
  allFqPairs.groupTuple(by: 2) | metricsFastqc

  // Main alignment processes
  bams = allFqPairs | trimFastqs | alignBwaSortSam | dedupCreak

  // Bam metrics
  bams | depthSamtools
  hsmetrics = bams | collectHsMetrics
  hsmetrics.groupTuple(by: 2) | multiqcHsmetrics

  // VarDict variant calling + merging
  vardvcfs = bams | callVarDict | decomposeVarDict
  vardvcfs.groupTuple(by: 2) | mergeVarDict
}