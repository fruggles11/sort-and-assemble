#!/usr/bin/env nextflow

nextflow.enable.dsl = 2



// WORKFLOW SPECIFICATION
// --------------------------------------------------------------- //
workflow {

	// input channels
    ch_fastq_dirs = Channel
        .fromPath( "${params.fastq_dir}/barcode*/", type: 'dir' )

    ch_primers = Channel
        .fromPath ( params.primer_table )
        .splitCsv ( header: true )
        .map { row -> tuple( row.chain, row.primer_seq, row.differentiating ) }
        .filter { it[2] == "true" }
        .groupTuple ( by: 0 )
		.map { chain, primers, diff -> tuple( chain, primers ) }
	
	// Workflow steps 
    MERGE_BY_BARCODE (
        ch_fastq_dirs
    )

    FIND_ADAPTER_SEQS (
        MERGE_BY_BARCODE.out
			.map { fastq -> tuple( file(fastq), file(fastq).countFastq() ) }
			.filter { it[1] > params.min_reads }
    )

    SPLIT_BY_PRIMER (
        MERGE_BY_BARCODE.out
			.map { fastq -> tuple( file(fastq), file(fastq).countFastq() ) }
			.filter { it[1] > params.min_reads }
			.map { fastq, count -> fastq },
        ch_primers
    )

    QC_TRIMMING (
        SPLIT_BY_PRIMER.out
			.combine( FIND_ADAPTER_SEQS.out, by: 1 )
    )

    READ_STATS (
        QC_TRIMMING.out
    )

    // VISUALIZE_STATS (
    //     READ_STATS.out
    // )

    // ASSEMBLE_WITH_CANU (
    //     QC_TRIMMING.out
    // )

    // SEARCH_IGBLAST (
    //     ASSEMBLE_WITH_CANU.out
    // )
	
	
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
	
	tag "${barcode}"
	publishDir params.merged_reads, mode: 'copy'

	cpus 4
	
	input:
	path read_dir
	
	output:
	path "*.fastq.gz"
	
	script:
	barcode = read_dir.getName()
	"""
    seqkit scat -j 4 --out-format fastq -f `realpath ${read_dir}` -o ${barcode}.fastq.gz
	"""

}

process FIND_ADAPTER_SEQS {
	
	/* */
	
	tag "${sample_id}"
	// publishDir params.merged_reads, mode: 'copy'

	errorStrategy 'ignore'
	
	input:
	tuple path(merged_reads), val(count)
	
	output:
	tuple path("*_adapters.fasta"), val(sample_id)
	
	script:
	sample_id = merged_reads.getSimpleName()
	"""
    bbmerge.sh in="${merged_reads}" outa="${sample_id}_adapters.fasta" ow qin=33
	"""

}

process SPLIT_BY_PRIMER {
	
	/* */
	
	tag "${sample_id}"
	publishDir params.split_reads, mode: 'copy'
	
	input:
	each path(merged_reads)
    tuple val(primer_id), val(primer_seqs)

	output:
	tuple path("*.fastq.gz"), val(sample_id)
	
	script:
	sample_id = merged_reads.getSimpleName()
	seq_patterns = primer_seqs
		.toString()
		.replace("[", "")
		.replace("]", "")
		.replace("'", "")
		.replace(" ", "")
		.replace(",", "\n")
	"""
	echo "${seq_patterns}" > primers.txt
    seqkit grep -j ${task.cpus} -f primers.txt -m 1 ${merged_reads} \
	-o ${sample_id}_${primer_id}.fastq.gz
	"""

}

process QC_TRIMMING {
	
	/* */
	
	tag "${sample_id}"
	publishDir params.trimmed_reads, mode: 'copy'

	cpus 4
	
	input:
	tuple val(sample_id), path(split_reads), path(adapters)
	
	output:
	tuple path("*.fastq.gz"), val(sample_id), val(primer_id)
	
	script:
	primer_id = split_reads.getSimpleName().split("_")[1]
	"""
	reformat.sh in=${split_reads} \
	out=${sample_id}_${primer_id}_filtered.fastq.gz \
	ref=${adapters} \
	forcetrimleft=30 forcetrimright2=30 \
	mincalledquality=9 qin=33 minlength=200 \
	uniquenames=t t=${task.cpus}
	"""

}

process READ_STATS {
	
	/* */
	
	tag "${sample_id}"
	publishDir params.read_stats, mode: 'copy'

	cpus 4
	
	input:
	tuple path(qc_reads), val(sample_id), val(primer_id)
	
	
	output:
    tuple path("*.tsv"), val(sample_id), val(primer_id)
	
	script:
	"""
	seqkit stats -j ${task.cpus} --tabular -a \
	${qc_reads} > ${sample_id}_${primer_id}.tsv
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
	csvtk -t plot hist -f 4 --format pdf ${stats_tsv} \
	-o ${sample_id}_${primer_id}_num_seqs.pdf
	csvtk -t plot hist -f 13 --format pdf ${stats_tsv} \
	-o ${sample_id}_${primer_id}_n50.pdf
	csvtk -t plot hist -f 15 --format pdf ${stats_tsv} \
	-o ${sample_id}_${primer_id}_q30.pdf
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