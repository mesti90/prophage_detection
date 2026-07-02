#!/usr/bin/env bash

set -euo pipefail

usage() {
	cat <<'EOF'
Usage:
	./run_pipeline.sh [options]

Options:
	-i INPUT_TSV          Input TSV with Sample_ID, Biosample, and Species columns [default: biosamples.tsv]
	-o OUTDIR             Base directory for pipeline outputs [default: pipeline_outputs]
	-c CONTAINER          Singularity container for assembly lookup [default: /home/vasarhelyib/containers/mesti90-ncbi_edirect.24.7.20250903.sif]
	-d SRA_CONTAINER      Singularity container for SRA tools [default: /home/vasarhelyib/containers/ncbi-sra-tools.3.4.1.sif]
	-t THREADS            Threads for fasterq-dump [default: 8]
	--steps 1,2,3,4,5     Run only specific pipeline steps
	--dry-run             Print the commands without executing them
	-h, --help            Show this help message
EOF
}

INPUT="biosamples.tsv"
OUTDIR="pipeline_outputs"
ASSEMBLY_CONTAINER="/home/vasarhelyib/containers/mesti90-ncbi_edirect.24.7.20250903.sif"
SRA_CONTAINER="/home/vasarhelyib/containers/ncbi-sra-tools.3.4.1.sif"
THREADS=8
STEPS=(1 2 3 4 5)
DRY_RUN=0

while [[ $# -gt 0 ]]; do
	case "$1" in
		-i)
			INPUT="$2"
			shift 2
			;;
		-o)
			OUTDIR="$2"
			shift 2
			;;
		-c)
			ASSEMBLY_CONTAINER="$2"
			shift 2
			;;
		-d)
			SRA_CONTAINER="$2"
			shift 2
			;;
		-t)
			THREADS="$2"
			shift 2
			;;
		--steps)
			IFS=',' read -r -a STEPS <<< "$2"
			shift 2
			;;
		--dry-run)
			DRY_RUN=1
			shift
			;;
		-h|--help)
			usage
			exit 0
			;;
		*)
			echo "Unknown option: $1" >&2
			usage >&2
			exit 1
			;;
	esac
done

repo_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$repo_dir"

mkdir -p "$OUTDIR"
mkdir -p "$OUTDIR/genomes" "$OUTDIR/sra_reads" "$OUTDIR/genomad" "$OUTDIR/genomad_results" "$OUTDIR/clusters"

ln -sfn "$OUTDIR/genomes" genomes
ln -sfn "$OUTDIR/sra_reads" sra_reads
ln -sfn "$OUTDIR/genomad" genomad
ln -sfn "$OUTDIR/genomad_results" genomad_results
ln -sfn "$OUTDIR/clusters" clusters

if [[ "$INPUT" != "biosamples.tsv" ]]; then
	ln -sfn "$INPUT" biosamples.tsv
fi

run_step() {
	local step="$1"
	case "$step" in
		1)
			cmd=(bash ./1_download_assemblies.sh -i "$INPUT" -o "$OUTDIR/genomes" -r "$OUTDIR/sra_reads" -f "$OUTDIR/failed_biosamples.txt" -c "$ASSEMBLY_CONTAINER" -d "$SRA_CONTAINER" -t "$OUTDIR/assemblies.tsv" -n "$THREADS")
			;;
		2)
			cmd=(python3 ./2_download_sra_reads.py)
			;;
		3)
			cmd=(python3 ./3_collect_input_for_genomad.py)
			;;
		4)
			cmd=(python3 ./4_batch_genomad.py)
			;;
		5)
			cmd=(python3 ./5_genomad_postprocess.py --out_dir "$OUTDIR/phage_resistant_pipeline")
			;;
		*)
			echo "Unsupported step: $step" >&2
			exit 1
			;;
	esac

	echo "==> Step $step"
	printf '    '
	printf '%q ' "${cmd[@]}"
	printf '\n'

	if (( DRY_RUN )); then
		return
	fi

	"${cmd[@]}"
}

for step in "${STEPS[@]}"; do
	run_step "$step"
done

echo "Pipeline completed. Outputs are under $OUTDIR"
