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

ref_genome = file(input + "/Aeaegypti_ref/reference.fa")

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
trimmed_fqs.into { trimmed_reads_bwa; trimmed_reads_qc ; trimmed_reads_picard }



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
// ** - base recalibration BaseRecalibrator & ApplyBQSR  (SKIP since no known sites vcf)
////////////////////////////////////////////////

// Get a list of known Ae. aegypti LVP variants from Ensembl.
process fetch_variants {

    cpus small

    output:
        tuple file("known_variants.vcf.gz"), file("known_variants.vcf.gz.csi"), file("known_variants.vcf.gz.tbi") into known_variants

    """
        wget ftp://ftp.ensemblgenomes.org/pub/metazoa/release-55/variation/vcf/aedes_aegypti_lvpagwg/aedes_aegypti_lvpagwg.vcf.gz \
          -O known_variants.vcf.gz
        wget ftp://ftp.ensemblgenomes.org/pub/metazoa/release-55/variation/vcf/aedes_aegypti_lvpagwg/aedes_aegypti_lvpagwg.vcf.gz.csi \
          -O known_variants.vcf.gz.csi
        gatk IndexFeatureFile \
          -I known_variants.vcf.gz
    """
}

// base recalibration
process base_recalibration {

    cpus big
    tag { id }

    input:
        tuple val(id), file(bam) from marked_bams
        file("reference.fa") from ref_genome
        tuple file(vcf), file(index_csi), file(index_tbi) from known_variants

    output:
        tuple id, file("${bam}"), file("${id}_recal_data.table") into brdt

    """
        gatk CreateSequenceDictionary -R reference.fa
        samtools faidx reference.fa

        gatk BaseRecalibrator \
          -I ${bam} \
          -R reference.fa \
          --known-sites ${vcf} \
          -O ${id}_recal_data.table

    """
}

// apply recalibration
process apply_recalibration {
    //publishDir "${output}/${params.dir}/base_recal", mode: 'copy', pattern: '*.pdf'

    cpus big
    tag { id }

    input:
        tuple val(id), file(bam), file(recal_table) from brdt
        file("reference.fa") from ref_genome

    output:
        tuple id, file("${id}_recal.bam") into recal_bams
        file ("${id}_AnalyzeCovariates.pdf") into recal_logs

    """
        gatk CreateSequenceDictionary -R reference.fa
        samtools faidx reference.fa

        gatk ApplyBQSR \
          -R reference.fa \
          -I ${bam} \
          --bqsr-recal-file ${recal_table} \
          -O "${id}_recal.bam"

    """
}

 //       gatk AnalyzeCovariates \
  //        -bqsr ${recal_table} \
  //        -plots "${id}_AnalyzeCovariates.pdf"



////////////////////////////////////////////////
// ** - VARIANT CALLING PIPELINE
////////////////////////////////////////////////

// single-sample HaplotypeCaller -> GVCFs
process haplotype_caller {
   publishDir "${output}/${params.dir}/gvcf", mode: 'copy', pattern: '*.vcf.gz'

    cpus big
    tag { id }

    input:
        tuple val(id), file(bam) from recal_bams
        file ("reference.fa") from ref_genome

    output:
        tuple val(id), file("${id}.vcf.gz") into gvcfs

    """
        gatk CreateSequenceDictionary -R reference.fa
        samtools faidx reference.fa
        samtools index -@ ${task.cpus} -b ${bam}

        gatk --java-options "-Xmx4g" HaplotypeCaller  \
          -R reference.fa \
          -I ${bam} \
          -O ${id}.vcf.gz \
          -ERC GVCF

        
    """
}

gvcfs.into {gvcfs_map; gvcfs_combine}
sample_map = gvcfs_map.map { "${it[0]}\t${it[0]}vcf.gz" }.collectFile(name: "sample_map.tsv", newLine: true)


// GenomicsDBImport: import single-sample GVCFs *stuck here*
process combine_gvcfs {

    cpus big

    input:
      file (gvcfs) from gvcfs_combine.collect().ifEmpty([])
      file("sample_map.tsv") from sample_map

    output:
      file("gvcf_db.txt") into gvcf_db_success

    """
        mkdir ${work}/gvcf_db

        cat sample_map.tsv

        gatk --java-options "-Xmx4g -Xms4g" \
          GenomicsDBImport \
          --genomicsdb-workspace-path ${work}/gvcf_db \
          --sample-name-map sample_map.tsv \
          --tmp-dir . \
          --reader-threads ${task.cpus}

        echo "gvcf_db created" > gvcf_db.txt
    """
}


// Joint-Call Cohort GenotypeGVCFs
process joint_call_gvcfs {
   publishDir "${output}/${params.dir}/gvcf", mode: 'copy', pattern: 'output.vcf.gz'

    cpus big

    input:
        file ("reference.fa") from ref_genome
        file("gvcf_db.txt") from gvcf_db_success
    
    output:
        file("output.vcf.gz") into joint_vcf

    """
      gatk CreateSequenceDictionary -R reference.fa
      samtools faidx reference.fa

      gatk --java-options "-Xmx4g" GenotypeGVCFs \
        -R reference.fa \
        -V ${work}/gvcf_db \
        -O output.vcf.gz \
        --tmp-dir . \

    """
}



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
