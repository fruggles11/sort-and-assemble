#!/usr/bin/env nextflow

nextflow.enable.dsl = 2



// WORKFLOW SPECIFICATION
// --------------------------------------------------------------- //
workflow {
	
	
	// input channels
    ch_fastqs = Channel
        .fromPath( "${params.fastq_dir}/*.fastq.gz" )
        .collect()

    ch_barcodes

    ch_primers
	
	
	// Workflow steps 
    MERGE_BY_BARCODE (
        ch_fastqs,
        ch_barcodes
    )

    FIND_ADAPTER_SEQS (
        MERGE_BY_BARCODE.out
    )

    SPLIT_BY_PRIMER (
        MERGE_BY_BARCODE.out,
        ch_primers
    )

    QC_TRIMMING (
        SPLIT_BY_PRIMER.out,
        FIND_ADAPTER_SEQS.out
    )

    READ_STATS (
        QC_TRIMMING.out
    )

    VISUALIZE_STATS (
        READ_STATS.out
    )

    ASSEMBLE_WITH_CANU (
        QC_TRIMMING.out
    )

    SEARCH_IGBLAST (
        ASSEMBLE_WITH_CANU.out
    )
	
	
}
// --------------------------------------------------------------- //



// DERIVATIVE PARAMETER SPECIFICATION
// --------------------------------------------------------------- //
// Additional parameters that are derived from parameters set in nextflow.config
params.merged_reads = params.results + "/1_merged_reads"
params.split_reads = params.results + "/2_split_reads"
params.read_qc = params.results + "/3_read_QC"
params.trimmed_reads = params.read_qc + "/trimmed_reads"
params.read_stats = params.read_qc + "/read_stats"
params.assembly_results  = params.results + "/4_assembly_results"
params.ig_blast = params.assembly_results + "/IgBLAST"

// --------------------------------------------------------------- //




// PROCESS SPECIFICATION 
// --------------------------------------------------------------- //

process MERGE_BY_BARCODE {
	
	/* */
	
	tag "${sample_id}"
	publishDir params.merged_reads, mode: 'copy'
	
	input:
	each path(fastqs)
    tuple val(sample_id), val(barcode)
	
	output:
	path "*.fastq.gz"
	
	script:
	"""
    seqkit grep -j ${task.cpus} -p ${barcode} -m 1 ${fastqs} -o ${sample_id}.fastq.gz
	"""

}

process FIND_ADAPTER_SEQS {
	
	/* */
	
	tag "${sample_id}"
	publishDir params.merged_reads, mode: 'copy'
	
	input:
	path demux_reads
	
	output:
	tuple path("*_adapters.fasta"), val(sample_id)

    when:
    file(demux_reads.toString()).countFastq() > 50
	
	script:
	sample_id = demux_reads.getSimpleName()
	"""
    bbmerge.sh in="${demux_reads}" outa="${sample_id}_adapters.fasta" ow reads=1m
	"""

}

process SPLIT_BY_PRIMER {
	
	/* */
	
	tag "${sample_id}"
	publishDir params.params.split_reads, mode: 'copy'
	
	input:
	each path(demux_reads)
    tuple val(primer), val(primer_id)

	output:
	tuple path("*.fastq.gz"), val(sample_id), val(primer_id)

    when:
    file(demux_reads.toString()).countFastq() > 50
	
	script:
	sample_id = demux_reads.getSimpleName()
	"""
    seqkit grep -j ${task.cpus} -p ${primer} -m 1 ${demux_reads} -o ${sample_id}_${primer_id}.fastq.gz
	"""

}

process QC_TRIMMING {
	
	/* */
	
	tag "${sample_id}"
	publishDir params.trimmed_reads, mode: 'copy'
	
	input:
	tuple path(split_reads), val(sample_id), val(primer_id)
    tuple path(adapters), val(adapter_id)
	
	output:
	tuple path("*.fastq.gz"), val(sample_id), val(primer_id)

    when:
    sample_id == adapter_id
	
	script:
	"""
	"""

}

process READ_STATS {
	
	/* */
	
	tag "${sample_id}"
	publishDir params.read_stats, mode: 'copy'
	
	input:
	tuple path(qc_reads), val(sample_id), val(primer_id)
	
	
	output:
    tuple path("*.tsv"), val(sample_id), val(primer_id)
	
	script:
	"""
	"""

}

process VISUALIZE_STATS {
	
	/* */
	
	tag "${sample_id}"
	publishDir params.read_stats, mode: 'copy'
	
	input:
    tuple path(stats_tsv), val(sample_id), val(primer_id)
	
	output:
	path "*"
	
	script:
	"""
	"""

}

process ASSEMBLE_WITH_CANU {
	
	/* */
	
	tag "${tag}"
	publishDir params.assembly_results, mode: 'copy'
	
	input:
	tuple path(qc_reads), val(sample_id), val(primer_id)
	
	output:
	tuple path("*.fasta"), val(sample_id), val(primer_id)
	
	script:
	"""
	"""

}

process SEARCH_IGBLAST {
	
	/* */
	
	tag "${sample_id}"
	publishDir params.ig_blast, mode: 'copy'
	
	input:
	tuple path(fasta), val(sample_id), val(primer_id)
	
	output:
	
	
	script:
	"""
	"""

}

// --------------------------------------------------------------- //