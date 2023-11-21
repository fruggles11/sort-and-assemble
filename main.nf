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
	
	ch_file_list = Channel
		.fromPath( params.url_list )

	ch_blast_files = Channel
		.fromPath ( params.mamu_database_files )
	
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

    TRIM_ADAPTERS (
        SPLIT_BY_PRIMER.out
			.combine( FIND_ADAPTER_SEQS.out, by: 1 )
			.map { id, reads, adapters -> tuple( id, file(reads), file(reads).countFastq(), file(adapters) ) }
			.filter { it[2] > params.min_reads }
			.map { id, reads, count, adapters -> tuple( id, file(reads), file(adapters) ) }
    )

	QC_VALIDATION (
		TRIM_ADAPTERS.out
	)

    READ_STATS (
        QC_VALIDATION.out
			.map { tsv, sample, primer -> tsv }
			.collect()
    )

    VISUALIZE_STATS (
        READ_STATS.out
    )

    CONVERT_TO_FASTA (
        QC_VALIDATION.out
    )

	ASSEMBLE_WITH_CANU (
		CONVERT_TO_FASTA.out
	)

	// DEDUP_CONTIGS (
	// 	ASSEMBLE_WITH_CANU.out.contigs
	// )

	// PULL_IMGT_REFS (
	// 	ch_file_list
	// )

	// PULL_MAMU_DATABASE (
	// 	ch_blast_files
	// )

	// BUILD_IGBLAST_DATABASE (
	// 	PULL_IMGT_REFS.out
	// 		.collect()
	// )

	// BUNDLE_DATABASES (
	// 	PULL_MAMU_DATABASE.out,
	// 	BUILD_IGBLAST_DATABASE.out
	// )

    // SEARCH_IGBLAST (
	// 	BUNDLE_DATABASES.out,
    //     CLUSTER_BY_IDENTITY.out
    // )
	
	
}
// --------------------------------------------------------------- //



// DERIVATIVE PARAMETER SPECIFICATION
// --------------------------------------------------------------- //
// Additional parameters that are derived from parameters set in nextflow.config
if ( params.debugmode == true ){
	errorMode = 'terminate'
} else {
	errorMode = 'ignore'
}

params.merged_reads = params.results + "/1_merged_reads"
params.split_reads = params.results + "/2_split_reads"
params.read_qc = params.results + "/3_read_QC"
params.trimmed_reads = params.read_qc + "/trimmed_reads"
params.read_stats = params.read_qc + "/read_stats"
params.assembly_results  = params.results + "/4_assembly_results"
params.clustering_results  = params.results + "/4_clustering_results"
params.ig_blast = params.assembly_results + "/IgBLAST"

// --------------------------------------------------------------- //




// PROCESS SPECIFICATION 
// --------------------------------------------------------------- //

process MERGE_BY_BARCODE {
	
	/* */
	
	tag "${barcode}"
	publishDir params.merged_reads, mode: 'copy', overwrite: true

	errorStrategy { task.attempt < 3 ? 'retry' : errorMode }
	maxRetries 2

	cpus 4
	
	input:
	path read_dir
	
	output:
	path "${barcode}.fastq.gz"
	
	script:
	barcode = read_dir.getName()
	"""
    seqkit scat -j 4 -f `realpath ${read_dir}` -o ${barcode}.fastq.gz
	"""

}

process FIND_ADAPTER_SEQS {
	
	/* */
	
	tag "${sample_id}"

	errorStrategy { task.attempt < 3 ? 'retry' : errorMode }
	maxRetries 2
	
	input:
	tuple path(merged_reads), val(count)
	
	output:
	tuple path("${sample_id}_adapters.fasta"), val(sample_id)
	
	script:
	sample_id = merged_reads.getSimpleName()
	"""
    bbmerge.sh in=`realpath ${merged_reads}` outa="${sample_id}_adapters.fasta" ow qin=33
	"""

}

process SPLIT_BY_PRIMER {
	
	/* */
	
	tag "${sample_id}"
	publishDir params.split_reads, mode: 'copy', overwrite: true

	errorStrategy { task.attempt < 3 ? 'retry' : errorMode }
	maxRetries 2
	
	input:
	each path(merged_reads)
    tuple val(primer_id), val(primer_seqs)

	output:
	tuple path("${sample_id}_${primer_id}.fastq.gz"), val(sample_id)
	
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
    seqkit grep \
	--ignore-case \
	--by-seq \
	--delete-matched \
	-f primers.txt \
	-m ${params.primer_mismatch} \
	-j ${task.cpus} \
	${merged_reads} \
	-o ${sample_id}_${primer_id}.fastq.gz
	"""

}

process TRIM_ADAPTERS {
	
	/* */
	
	tag "${sample_id}, ${primer_id}"

	errorStrategy { task.attempt < 3 ? 'retry' : errorMode }
	maxRetries 2

	cpus 4
	
	input:
	tuple val(sample_id), path(split_reads), path(adapters)
	
	output:
	tuple path("${sample_id}_${primer_id}_trimmed.fastq.gz"), val(sample_id), val(primer_id)
	
	script:
	primer_id = split_reads.getSimpleName().split("_")[1]
	"""
	reformat.sh in=`realpath ${split_reads}` \
	out=${sample_id}_${primer_id}_trimmed.fastq.gz \
	ref=`realpath ${adapters}` \
	mincalledquality=9 qin=33 \
	minlength=${params.min_len} maxlength=${params.max_len} \
	uniquenames=t overwrite=true tossbrokenreads t=${task.cpus}
	"""

}

process QC_VALIDATION {
	
	/* */
	
	tag "${sample_id}, ${primer_id}"
	publishDir params.trimmed_reads, mode: 'copy', overwrite: true

	errorStrategy 'ignore'

	cpus 4
	
	input:
	tuple path(split_reads), val(sample_id), val(primer_id)
	
	output:
	tuple path("${sample_id}_${primer_id}_filtered.fastq.gz"), val(sample_id), val(primer_id)
	
	script:
	"""
	seqkit seq \
	--min-len ${params.min_len} \
	--max-len ${params.max_len} \
	--validate-seq \
	--threads ${task.cpus} \
	--min-qual 9.0 \
	${split_reads} \
	-o ${sample_id}_${primer_id}_filtered.fastq.gz
	"""

}

process READ_STATS {
	
	/* */
	
	publishDir params.read_stats, mode: 'copy', overwrite: true

	errorStrategy { task.attempt < 3 ? 'retry' : errorMode }
	maxRetries 2

	cpus 4
	
	input:
	path fastqs 
	
	output:
    path "sequencing_run_stats.tsv"
	
	script:
	"""
	seqkit stats -j ${task.cpus} --tabular -a \
	*.fastq.* > sequencing_run_stats.tsv
	"""

}

process VISUALIZE_STATS {
	
	/* */

	publishDir params.read_stats, mode: 'copy', overwrite: true

	errorStrategy { task.attempt < 3 ? 'retry' : errorMode }
	maxRetries 2
	
	input:
    path stats_tsv
	
	output:
	path "*"
	
	script:
	"""
	csvtk -t plot hist -f 4 --format pdf ${stats_tsv} \
	-o num_seqs_histogram.pdf
	csvtk -t plot hist -f 13 --format pdf ${stats_tsv} \
	-o n50_histogram.pdf
	csvtk -t plot hist -f 15 --format pdf ${stats_tsv} \
	-o q30_histogram.pdf
	"""

}

process CONVERT_TO_FASTA {

	/* */
	
	tag "${sample_id}, ${primer_id}"

	errorStrategy { task.attempt < 3 ? 'retry' : errorMode }
	maxRetries 2

	input:
	tuple path(reads), val(sample_id), val(primer_id)

	output:
	tuple path("${sample_id}_${primer_id}.fasta"), val(sample_id), val(primer_id)

	script:
	"""
	seqkit fq2fa ${reads} -o ${sample_id}_${primer_id}.fasta
	"""

}

process CLUSTER_BY_IDENTITY {

	/* */
	
	tag "${sample_id}, ${primer_id}"
	publishDir "${params.clustering_results}/${sample_id}_${primer_id}", mode: 'copy', overwrite: true

	errorStrategy { task.attempt < 3 ? 'retry' : errorMode }
	maxRetries 2

	cpus 4
	
	input:
	tuple path(fasta), val(sample_id), val(primer_id)
	
	output:
	tuple path("${sample_id}_${primer_id}-*.fasta"), val(sample_id), val(primer_id), emit: contigs
	path "cluster_fastqs/", emit: reads

	script:
	"""
	vsearch --cluster_fast ${fasta} \
	--id ${params.id_threshold} \
	--clusters ${sample_id}_${primer_id}-cluster-seqs \
	--consout ${sample_id}_${primer_id}-consensus.fasta \
	--iddef 3 \
	--threads ${task.cpus} && \
	mkdir -p cluster_fastqs && \
	for file in ${sample_id}_${primer_id}-cluster-seqs*; do
		if [ -f "\$file" ]; then
			mv "\$file" "cluster_fastqs/\${file}.fasta"
		fi
	done
	"""

}

process ASSEMBLE_WITH_CANU {
	
	/* */
	
	tag "${sample_id}, ${primer_id}"
	publishDir "${params.assembly_results}/${sample_id}_${primer_id}", mode: 'copy'
	errorStrategy { task.attempt < 3 ? 'retry' : errorMode }
	maxRetries 2
	cpus 4
	
	input:
	tuple path(qc_reads), val(sample_id), val(primer_id)
	
	output:
	tuple path("*.fasta"), val(sample_id), val(primer_id)
	
	script:
	"""
	canu \
	-p "${sample_id}-${primer_id}" -d . \
	'genomesize=1000' \
	'minreadlength=600' \
	-trimmed \
	-nanopore `realpath ${qc_reads}`
	"""

}

process DEDUP_CONTIGS {

	/* */
	
	tag "${sample_id}, ${primer_id}"
	publishDir "${params.clustering_results}/${sample_id}_${primer_id}", mode: 'copy', overwrite: true

	errorStrategy { task.attempt < 3 ? 'retry' : errorMode }
	maxRetries 2

	input:
	tuple path(fasta), val(sample_id), val(primer_id)

	output:
	tuple path("${sample_id}_${primer_id}_deduped.fasta"), val(sample_id), val(primer_id)

	script:
	"""
	dedup_and_recal.py \
	--fasta ${fasta} \
	--output ${sample_id}_${primer_id}_deduped \
	--split_char "=" \
	--min_depth 2
	"""
}

process PULL_IMGT_REFS {

	/* */

	input:
	path url_list

	output:
	path "I*.fasta"

	script:
	"""
	goDownloadFiles \
	-http ${url_list}
	"""

}

process PULL_MAMU_DATABASE {

	/* */

	input:
	path blast_files

	output:
	path "*"

	script:
	"""
	goDownloadFiles -http ${blast_files} && \
	mkdir database && \
	tar -xvf rhesus_monkey_VJ.tar -C database && \
	mkdir -p internal_data/rhesus_monkey && \
	cp rhesus_monkey_* internal_data/rhesus_monkey/
	"""

}

process BUILD_IGBLAST_DATABASE {

	/* */

	input:
	path fastas

	output:
	path "imgt_db*"

	script:
	"""
	cat *.fasta > merged.fasta && \
	edit_imgt_file.pl merged.fasta > imgt_db && \
	makeblastdb -parse_seqids -hash_index -dbtype nucl -in imgt_db
	"""

}

process BUNDLE_DATABASES {

	/* */

	input:
	path imgt_refs
	path mamu_db

	output:
	path "databases.tar"

	script:
	"""
	tar -cvf databases.tar imgt_db* rhesus_monkey_* database/ internal_data/
	"""

}

process SEARCH_IGBLAST {
	
	/* */
	
	tag "${sample_id}, ${primer_id}"
	publishDir params.ig_blast, mode: 'copy', overwrite: true

	errorStrategy 'ignore' // { task.attempt < 3 ? 'retry' : errorMode }
	// maxRetries 2
	
	input:
	each path(databases)
	tuple path(fasta), val(sample_id), val(primer_id)
	
	output:
	path "*"
	
	script:
	"""
	tar -xvf databases.tar
	igblastn \
	-organism rhesus_monkey \
	-query ${fasta} \
	-db imgt_db
	"""

}

// --------------------------------------------------------------- //