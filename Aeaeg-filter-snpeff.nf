#!/usr/bin/env nextflow

// Params from config files (system-dependent)
input=params.input
output=params.output
work=params.work
aux=params.aux

huge=params.huge
big=params.big
med=params.med
small=params.small

// Global Params
params.dir = null
if( !params.dir ) error "Missing dir parameter"
println "dir: $params.dir"


////////////////////////////////////////////////
// ** - Pull in vcf file
////////////////////////////////////////////////
input_vcf = file(input + "/*.vcf.gz" )

//ref_genome = file(input + "/Aeaegypti_ref/reference.fa")


////////////////////////////////////////////////
// ** - hard-filter variants
////////////////////////////////////////////////

// split snps and indels
process split_snps_indels {

    cpus small
    tag { unfilt_vcf }

    input:
      file(unfilt_vcf) from input_vcf

    output:
      file "snps.vcf.gz" into snps_vcf
      file "snps.vcf.gz.tbi" into snps_vcf_index

      file "indels.vcf.gz" into indels_vcf
      file "indels.vcf.gz.tbi" into indels_vcf_index

    script:

    """
      gatk IndexFeatureFile -I ${unfilt_vcf}

      gatk SelectVariants \
        -V ${unfilt_vcf} \
        -select-type SNP \
        -O snps.vcf.gz
      gatk IndexFeatureFile -I snps.vcf.gz

      gatk SelectVariants \
        -V ${unfilt_vcf} \
        -select-type INDEL \
        -O indels.vcf.gz
      gatk IndexFeatureFile -I indels.vcf.gz
    """
}


// filter variants
process filter_snps {

    cpus small

    input:
      file(snps) from snps_vcf
      file(snps_index) from snps_vcf_index

    output:
      file "snps_filtered.vcf.gz" into snps_filtered_vcf
      file "snps_filtered.vcf.gz.tbi" into snps_filtered_vcf_index

    script:

    """
    gatk VariantFiltration \
      -V ${snps} \
      -filter "QD < 2.0" --filter-name "QD2" \
      -filter "QUAL < 30.0" --filter-name "QUAL30" \
      -filter "SOR > 3.0" --filter-name "SOR3" \
      -filter "FS > 60.0" --filter-name "FS60" \
      -filter "MQ < 40.0" --filter-name "MQ40" \
      -filter "MQRankSum < -12.5" --filter-name "MQRankSum-12.5" \
      -filter "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSum-8" \
      -O snps_filtered.vcf.gz

    gatk IndexFeatureFile -I snps_filtered.vcf.gz
    """
}

process filter_indels {

    cpus small

    input:
      file(indels) from indels_vcf
      file(indels_index) from indels_vcf_index

    output:
      file "indels_filtered.vcf.gz" into indels_filtered_vcf
      file "indels_filtered.vcf.gz.tbi" into indels_filtered_vcf_index

    script:

    """
    gatk VariantFiltration \
        -V ${indels} \
        -filter "QD < 2.0" --filter-name "QD2" \
        -filter "QUAL < 30.0" --filter-name "QUAL30" \
        -filter "FS > 200.0" --filter-name "FS200" \
        -filter "ReadPosRankSum < -20.0" --filter-name "ReadPosRankSum-20" \
        -O indels_filtered.vcf.gz

    gatk IndexFeatureFile -I indels_filtered.vcf.gz
    """
}

// merge snp and indel vcfs
process merge_vcf {

    publishDir "${output}/${params.dir}/vcf", mode: 'copy', pattern: 'filtered.vcf.gz'

    cpus small

    input:
      file (snps_filtered) from snps_filtered_vcf
      file (snps_filtered_index) from snps_filtered_vcf_index
      file (indels_filtered) from indels_filtered_vcf
      file (indels_filtered_index) from indels_filtered_vcf_index

    output:
      file "filtered.vcf.gz" into filtered_vcf

    script:

    """
      picard -Xmx8g MergeVcfs \
            I=${snps_filtered} \
            I=${indels_filtered} \
            O=filtered.vcf.gz
    """
}


////////////////////////////////////////////////
// ** - variant effect prediction (snpeff)
////////////////////////////////////////////////
