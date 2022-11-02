#!/usr/bin/env nextflow

// Params from config files (system-dependent)
input=params.input
output=params.output
work=params.work
aux=params.aux

huge=params.huge
big=params.big
small=params.small

// Global Params
params.dir = null
if( !params.dir ) error "Missing dir parameter"
println "dir: $params.dir"

params.rlen = null
if( !params.rlen ) error "Missing length (average read length) parameter"
println "rlen: $params.rlen"

// flag for fastqc and multiqc (--qc)
params.qc = false


////////////////////////////////////////////////
// ** - Pull in fq files and indexed genome
////////////////////////////////////////////////

Channel.fromFilePairs(input + "/${params.dir}/*_R{1,2}.fq.gz", flat: true) //for subsampled data
//Channel.fromFilePairs(input + "/${params.dir}/*_R{1,2}_001.f[a-z]*q.gz", flat: true) //for full dataset
          .set { fqs }

bwa_indices = Channel.fromPath(input + "/Aeaegypti_ref/reference.*" )

ref_genome = file(input + "/Aeaegypti_ref/reference.fasta")

////////////////////////////////////////////////
// ** - Trim reads
////////////////////////////////////////////////

process trim_reads {

  publishDir "${output}/${params.dir}/trim_stats/", mode: 'copy', pattern: '*.{json,html}'

  cpus small
  tag { id }

  input:
    tuple val(id), file(forward), file(reverse) from fqs

  output:
    tuple id, file("${id}_R1_trim.fq.gz"), file("${id}_R2_trim.fq.gz") into trimmed_fqs
    tuple file("*.html"), file("*.json")  into trim_log

  """
    fastp -i $forward -I $reverse -w ${task.cpus} -o ${id}_R1_trim.fq.gz -O ${id}_R2_trim.fq.gz -y -l 50 -h ${id}.html -j ${id}.json
  """

}
trimmed_fqs.into { trimmed_reads_bwa; trimmed_reads_qc ; trimmed_reads_picard}



////////////////////////////////////////////////
// ** - multiQC of trimmed fqs (--qc flag)
////////////////////////////////////////////////

process fastqc {

    publishDir "${output}/${params.dir}/fastqc", mode: 'copy', pattern: '*_fastqc.{zip,html}'

    cpus small
    tag { id }

    when:
      params.qc

    input:
    tuple val(id), file(forward), file(reverse) from trimmed_reads_qc

    output:
    file "*_fastqc.{zip,html}" into fastqc_results

    script:

    """
      fastqc -q $forward $reverse -t ${task.cpus}
    """
}

process multiqc {

    publishDir "${output}/${params.dir}/fastqc", mode: 'copy', pattern: 'multiqc_report.html'

    cpus small

    when:
      params.qc

    input:
    file ('fastqc/*') from fastqc_results.collect().ifEmpty([])

    output:
    file "multiqc_report.html" into multiqc_report

    script:

    """
      multiqc .
    """
}


////////////////////////////////////////////////
// ** - bwa mapping
////////////////////////////////////////////////

process bwa_align {
    publishDir "${output}/${params.dir}/bwa_stats/", mode: 'copy', pattern: '*.flagstat.txt'
    publishDir "${output}/${params.dir}/bams", mode: 'copy', pattern: '*.bam'
    publishDir "${output}/${params.dir}/bams", mode: 'copy', pattern: '*.bam.bai'

    cpus big
    tag { id }

    input:
        tuple val(id), file(forward), file(reverse) from trimmed_reads_bwa
        file bwa_indices from bwa_indices.collect()

    output:
        file("${id}.flagstat.txt") into bwa_stats
        tuple id, file("${id}.bam") into bam_files
        file("${id}.bam.bai") into bam_indices

    script:
        index_base = bwa_indices[0].toString() - ~/.fa[.a-z]*/
        readGroup = \
          "@RG\\tID:${id}\\tLB:${id}\\tPL:illumina\\tPM:novaseq\\tSM:${id}"

        """
        bwa mem \
          -K 100000000 \
          -v 3 \
          -t ${task.cpus} \
          -M \
          -Y \
          -R \"${readGroup}\" \
          ${index_base}.fa \
          ${forward} \
          ${reverse} \
          > ${id}.sam
        samtools view -@ ${task.cpus} -bS ${id}.sam > ${id}.unsorted.bam
        rm *.sam
        samtools flagstat ${id}.unsorted.bam
        samtools sort -@ ${task.cpus} -m 16G -o ${id}.bam ${id}.unsorted.bam
        rm *.unsorted.bam
        samtools index -@ ${task.cpus} -b ${id}.bam
        samtools flagstat ${id}.bam > ${id}.flagstat.txt
        """
}

////////////////////////////////////////////////
// ** - markdups and sortsam (Picard)
////////////////////////////////////////////////

process mark_dups {
    publishDir "${output}/${params.dir}/picard_stats", mode: 'copy', pattern: '*_marked_dups_stats.txt'

    cpus big
    tag { id }

    input:
        tuple val(id), file(bam) from bam_files

    output:
        tuple id, file("${id}_marked_dups.bam") into marked_bams
        file "${id}_marked_dups_stats.txt" into picard_logs

    """
        picard -Xmx8g MarkDuplicates \
          I=${bam} \
          O=${id}_marked_dups.unsorted.bam \
          M=${id}_marked_dups_stats.txt
        
        picard -Xmx8g SortSam \
          I=${id}_marked_dups.unsorted.bam \
          O=${id}_marked_dups.bam \
          SORT_ORDER=coordinate

    """
}


////////////////////////////////////////////////
// ** - base recalibration
////////////////////////////////////////////////

process base_recalibration {
    publishDir "${output}/${params.dir}/base_recal", mode: 'copy', pattern: 'recal_data.table'

    cpus big
    tag { id }

    input:
        tuple val(id), file(bam) from marked_bams
        file("reference.fa") from ref_genome

    output:

    """

        gatk BaseRecalibrator \
          -I ${bam} \
          -R reference.fa \
          -O recal_data.table

    """
}


// ////////////////////////////////////////////////
// // ** - VARIANT CALLING PIPELINE
// ////////////////////////////////////////////////

// // GATK likes to use both aligned and unaligned reads, so we first generate a
// // uBAM from original FASTQs. NOTE: the rg variable will need to be adjusted
// // based on the @RGs from the FASTQs (the regex could be made to be more robust.)
// process picard_fastq_uBAM {

//     cpus big
//     tag { id }

//     input:
//         tuple val(id), file(forward), file(reverse) from trimmed_reads_picard

//     output:
//         tuple id, file("${id}.ubam") into sorted_ubams

//     """
//         SM=`echo ${id} | cut -c1-3 | tr -d '\n'`
//         ID=${id}
//         LB=${id}
//         gatk \
//           --java-options \
//           -Xmx6g \
//           FastqToSam \
//           -F1 ${forward} \
//           -F2 ${reverse} \
//           -O ${id}.ubam \
//           -SM "\$SM" \
//           -RG "\$ID" \
//           -LB "\$LB" \
//           -PL illumina \
//           -SO queryname \
//           -TMP_DIR ${work}
//     """
// }
//
// // uBAMs and BAMs need to be QN sorted.
// process picard_sort_bam {
//
//     cpus big
//     tag { id }
//
//     input:
//         tuple val(id), file(bam) from bam_files
//
//     output:
//         tuple val(id), file("${id}_qn_sorted.bam") into sorted_bams
//
//     """
//         gatk \
//           --java-options \
//           -Xmx6g \
//           SortSam \
//           -I ${bam} \
//           -O ${id}_qn_sorted.bam \
//           -SO queryname \
//           -TMP_DIR ${work}
//     """
// }
//
// // This joins BAM and uBAM channels by ID.
// joined_bams = sorted_bams.join(sorted_ubams)
//
// // Merge BAM and uBAM.
// process picard_merge {
//
//     cpus big
//     tag { id }
//
//     input:
//         tuple val(id), file(bam), file(ubam) from joined_bams
//
//     output:
//         tuple val(id), file("${id}_merged.bam") into merged_bams
//
//     """
//         gatk \
//           --java-options \
//           -Xmx6g \
//           MergeBamAlignment \
//           -ALIGNED ${bam} \
//           -UNMAPPED ${ubam} \
//           -O "${id}_merged.bam" \
//           -R "${genome}/reference/reference.fa" \
//           -TMP_DIR ${work}
//     """
// }
