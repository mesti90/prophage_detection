#!/usr/bin/env bash

set -euo pipefail

usage() {
	cat << EOF
Usage:
	$(basename "$0") \
		-i samples.tsv \
		-o genomes \
		-f failed_biosamples.txt \
		-c container.sif

Input TSV must contain columns:
	Strain
	Biosample

Options:
	-i  Input TSV [default: biosamples.tsv]
	-o  Output directory for genomes [default: genomes]
	-r  Output directory for reads [default: sra_reads]
	-f  Failed BioSamples file [default: failed_biosamples.txt]
	-c  Singularity container [default: /home/vasarhelyib/containers/mesti90-ncbi_edirect.24.7.20250903.sif]
	-d  Singularity container for SRA tools [default: /home/vasarhelyib/containers/ncbi-sra-tools.3.4.1.sif]
	-t  Output table with assembly paths [default: assemblies.tsv]
	-n  Number of threads for fasterq-dump [default: 8]
	-h  Show this help message
	-s  Skip assembly download phase and run only fasterq-dump on entries in failed file
	-u Skip assembly download phase and fasterq-dump, just prepare the assembly input table for the assembler with the entries in the failed file
EOF
	exit 1
}

SKIP_TO_FASTERQ=0
SKIP_TO_ASSEMBLY=0
INPUT="biosamples.tsv"
OUTDIR="genomes"
READDIR="sra_reads"
FAILED="failed_biosamples.txt"
CONTAINER="/home/vasarhelyib/containers/mesti90-ncbi_edirect.24.7.20250903.sif"
SRA_CONTAINER="/home/vasarhelyib/containers/ncbi-sra-tools.3.4.1.sif"
OUTPUT_TABLE="assemblies.tsv"
THREADS=8

while getopts ":i:o:r:f:c:d:t:n:suh" opt; do
	case "${opt}" in
		i) INPUT="${OPTARG}" ;;
		o) OUTDIR="${OPTARG}" ;;
		r) READDIR="${OPTARG}" ;;
		f) FAILED="${OPTARG}" ;;
		c) CONTAINER="${OPTARG}" ;;
		d) SRA_CONTAINER="${OPTARG}" ;;
		t) OUTPUT_TABLE="${OPTARG}" ;;
		n) THREADS="${OPTARG}" ;;
		s) SKIP_TO_FASTERQ=1 ;;
		u) SKIP_TO_ASSEMBLY=1 ;;
		h) usage ;;
		*) usage ;;
	esac
done



if [[ "$SKIP_TO_FASTERQ" -eq 0 ]]; then

	mkdir -p "${OUTDIR}" "${READDIR}"
	: > "${FAILED}"
	echo -e "Strain\tBiosample\tSpecies\tR1\tR2\tAssembly" > "${FAILED}"

	HEADER=$(head -n1 "${INPUT}")
	echo -e "${HEADER}\tAssembly" > "${OUTPUT_TABLE}"

	STRAIN_COL=$(echo "${HEADER}" | tr '\t' '\n' | nl -v1 | awk '$2=="Strain"{print $1}')
	BIOSAMPLE_COL=$(echo "${HEADER}" | tr '\t' '\n' | nl -v1 | awk '$2=="Biosample"{print $1}')
	SPECIES_COL=$(echo "${HEADER}" | tr '\t' '\n' | nl -v1 | awk '$2=="Species"{print $1}')

	[[ -z "${STRAIN_COL}" ]] && { echo "Column 'Strain' not found"; exit 1; }
	[[ -z "${BIOSAMPLE_COL}" ]] && { echo "Column 'Biosample' not found"; exit 1; }
	[[ -z "${SPECIES_COL}" ]] && { echo "Column 'Species' not found"; exit 1; }

	tail -n +2 "${INPUT}" | while IFS= read -r LINE; do

		IFS=$'\t' read -r -a FIELDS <<< "${LINE}"

		STRAIN="${FIELDS[$((STRAIN_COL - 1))]}"
		BIOSAMPLE="${FIELDS[$((BIOSAMPLE_COL - 1))]}"
		SPECIES="${FIELDS[$((SPECIES_COL - 1))]}"
		[[ -z "${STRAIN}" ]] && continue
		[[ -z "${BIOSAMPLE}" ]] && continue
		[[ -z "${SPECIES}" ]] && continue

		OUTFILE="${OUTDIR}/${STRAIN}.fna.gz"
		echo -e "${LINE}\t${OUTFILE}" >> "${OUTPUT_TABLE}"

		if [[ -s "${OUTFILE}" ]]; then
			echo "SKIP (exists): ${STRAIN}"

			continue
		fi

		echo "Processing: ${STRAIN} (${BIOSAMPLE})"

		URL=$(
			singularity exec \
				</dev/null \
				--env LC_ALL=C,LANG=C,LANGUAGE=C \
				"${CONTAINER}" \
				bash -lc "
					esearch -db biosample -query '${BIOSAMPLE}' \
					| elink -target assembly \
					| esummary \
					| xtract -pattern DocumentSummary \
						-element FtpPath_RefSeq FtpPath_GenBank \
					| head -n1
				" \
			2>/dev/null \
			| awk '
				{
					base = ($1 != "" ? $1 : $2)

					if (base != "") {
						file = gensub(/.*\//, "", "g", base)
						print base "/" file "_genomic.fna.gz"
					}
				}
			'
		)

		if [[ -z "${URL}" ]]; then
			echo "FAILED: ${STRAIN} (${BIOSAMPLE})"
			echo -e "${STRAIN}\t${BIOSAMPLE}\t${SPECIES}\t${READDIR}/${BIOSAMPLE}_1.fastq.gz\t${READDIR}/${BIOSAMPLE}_2.fastq.gz\t${OUTDIR}/${STRAIN}.fna.gz" >> "${FAILED}"
	
			continue
		fi

		if ! wget -q -O "${OUTFILE}" "${URL}"; then
			echo "FAILED DOWNLOAD: ${STRAIN} (${BIOSAMPLE})"
			echo -e "${STRAIN}\t${BIOSAMPLE}\t${SPECIES}\t${READDIR}/${BIOSAMPLE}_1.fastq.gz\t${READDIR}/${BIOSAMPLE}_2.fastq.gz\t${OUTDIR}/${STRAIN}.fna.gz" >> "${FAILED}"
			rm -f "${OUTFILE}"
			continue
		fi

		echo "OK: ${STRAIN}"

	done
fi

