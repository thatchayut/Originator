#!/usr/bin/env bash

main () {
	local vcf_input_file=$1
	local barcodes_file=$2
	local bam_input_file=$3
	local popscle_sif_file=$4
	local output_path=$5

	# Check inputs
	if [[ ! -f "$vcf_input_file" ]]; then
		echo "VCF input file '$vcf_input_file' does not exist."
		return 1
	fi

	if [[ ! -f "$barcodes_file" ]]; then
		echo "Barcodes tsv file '$barcodes_file' does not exist."
		return 1
	fi

	if [[ ! -f "$bam_input_file" ]]; then
		echo "BAM input file '$bam_input_file' does not exist."
		return 1
	fi

	if [[ ! -f "$popscle_sif_file" ]]; then
		echo "Popscle SIF file '$popscle_sif_file' does not exist."
		return 1
	fi

	if [[ -z "$output_path" ]]; then
		echo "Output path is not provided."
		return 1
	fi


	local filename=$(basename $bam_input_file ".bam")
	local result_path="${output_path}/${filename}"
	local sorted_vcf_file="${result_path}/sorted_as_in_bam.vcf"
	local filtered_bam_file="${result_path}/filter_bam_file_for_popscle_dsc_pileup.bam"
	local pileup_file="${result_path}/to_demultiplex.pileup"
	local freemuxlet_result="${result_path}/freemuxlet"

	mkdir -p $result_path

	echo "Working on $filename"

	# Download required scripts
	if [[ ! -d ./exec/popscle_helper_tools ]]; then
		echo "popscle_helper_tools not found, downloading"
		git clone https://github.com/aertslab/popscle_helper_tools.git ./exec/popscle_helper_tools
	fi

	# Sort vcf in bam file order
	echo "Sort vcf in bam file order"
	./exec/popscle_helper_tools/sort_vcf_same_as_bam.sh ${bam_input_file} ${vcf_input_file} v > ${sorted_vcf_file}
	if [[ ! $? -eq 0 ]]; then
		return 1
	fi

	# Filter bam for pile up
	echo "Filter bam for pile up"
	./exec/popscle_helper_tools/filter_bam_file_for_popscle_dsc_pileup.sh ${bam_input_file} ${barcodes_file} ${sorted_vcf_file} ${filtered_bam_file}
	if [[ ! $? -eq 0 ]]; then
		return 1
	fi

	# Use filtered bam file for dsc-pileup
	echo "Use filtered bam file for dsc-pileup"
	singularity run ${popscle_sif_file} "dsc-pileup \
		--sam ${filtered_bam_file} \
		--vcf ${sorted_vcf_file} \
		--group-list ${barcodes_file} \
		--out ${pileup_file}"
	if [[ ! $? -eq 0 ]]; then
		return 1
	fi

	# Freemuxlet
	echo "Freemuxlet"
	singularity run ${popscle_sif_file} "freemuxlet \
		--plp ${pileup_file} \
		--nsample 2 \
		--out ${freemuxlet_result}"
	if [[ ! $? -eq 0 ]]; then
		return 1
	fi

	# Clean up
	echo "Cleaning up"
	rm -f $sorted_vcf_file $filtered_bam_file $pileup_file

	return 0
}

main $@
exit $?
