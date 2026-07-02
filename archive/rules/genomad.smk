###############################################
# rules/genomad.smk
###############################################

rule genomad:
	"""
	Run geNomad on a genome assembly.
	"""

	input:
		assembly="results/genomes/{strain}.fna.gz"

	output:
		directory("results/genomad/{strain}")

	params:
		db=config["genomad_db"]

	threads: 16

	log:
		"logs/genomad/{strain}.log"

	container:
		config["containers"]["genomad"]

	shell:
		r"""
		set -euo pipefail

		genomad end-to-end \
			{input.assembly} \
			{output} \
			{params.db} \
			--threads {threads} \
			> {log} 2>&1
		"""
rule collect_genomad_outputs:
	"""
	Collect plasmid and provirus FASTA files into a single directory.
	"""

	input:
		expand("results/genomad/{strain}", strain=STRAINS)

	output:
		touch("results/genomad/.collection_complete")

	run:
		from pathlib import Path
		import shutil

		provirus_dir = Path("results/genomad_results/provirus_sequences/fna")
		plasmid_dir = Path("results/genomad_results/plasmid_sequences/fna")

		provirus_dir.mkdir(parents=True, exist_ok=True)
		plasmid_dir.mkdir(parents=True, exist_ok=True)

		for strain in STRAINS:

			root = Path(f"results/genomad/{strain}")

			for f in root.rglob("*provirus*.fna"):
				shutil.copy2(f, provirus_dir / f"{strain}_provirus.fna")

			for f in root.rglob("*plasmid*.fna"):
				shutil.copy2(f, plasmid_dir / f"{strain}_plasmid.fna")

		Path(output[0]).touch()


rule genomad_postprocess:
	"""
	CheckV filtering, concatenate plasmids/proviruses,
	cluster with MMseqs2, and generate presence tables.
	"""

	input:
		"results/genomad/.collection_complete"

	output:
		provirus_fasta="results/postprocess/all_prophages.fasta",
		plasmid_fasta="results/postprocess/all_plasmids.fasta",

		provirus_rep="results/postprocess/clusters/provirus.rep_sequences.fasta",
		plasmid_rep="results/postprocess/clusters/plasmid.rep_sequences.fasta",

		provirus_presence="results/postprocess/clusters/provirus.presence.tsv",
		plasmid_presence="results/postprocess/clusters/plasmid.presence.tsv"

	params:
		checkv_db=config["checkv_db"]

	threads: 20

	log:
		"logs/postprocess.log"

	container:
		config["containers"]["python"]

	shell:
		r"""
		set -euo pipefail

		python scripts/genomad_postprocess.py \
			--provirus_dir results/genomad_results/provirus_sequences/fna \
			--plasmid_dir results/genomad_results/plasmid_sequences/fna \
			--checkv_db {params.checkv_db} \
			--out_dir results/postprocess \
			--threads {threads} \
			> {log} 2>&1
		"""