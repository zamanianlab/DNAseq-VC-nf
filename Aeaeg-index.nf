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


////////////////////////////////////////////////
// ** - Fetch genome and gene annotation files
////////////////////////////////////////////////

genome_url="https://vectorbase.org/common/downloads/Current_Release/AaegyptiLVP_AGWG/fasta/data/VectorBase-59_AaegyptiLVP_AGWG_Genome.fasta"
annot_url="https://vectorbase.org/common/downloads/Current_Release/AaegyptiLVP_AGWG/gff/data/VectorBase-59_AaegyptiLVP_AGWG.gff"

process fetch_ref {

    publishDir "${output}/${params.dir}/", mode: 'copy'

    output:
        file("reference.fa") into reference_fa

    """
        echo '${genome_url}'
        wget ${genome_url} -O reference.fa
        echo '${annot_url}'
        wget ${annot_url} -O geneset.gff
    """
}
reference_fa.into { bwa_index }


////////////////////////////////////////////////
// ** - Index Genome (bwa)
////////////////////////////////////////////////

process build_bwa_index {

    cpus huge

    publishDir "${output}/${params.dir}/", mode: 'copy'

    input:
        file("reference.fa") from bwa_index

    output:
        file "reference.*" into bwa_indices

    """
        bwa index reference.fa
    """
}