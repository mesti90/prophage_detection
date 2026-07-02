###############################################
# rules/acquire_genome.smk
###############################################

import json

##########################################################
# 1. Try downloading an assembly
##########################################################

rule fetch_assembly:

	input:
		# force execution once per sample
		touch=lambda wc: config["samples"]

	output:
		assembly=temp("results/downloaded_genomes/{strain}.fna.gz"),
		status="results/download/{strain}.json"

	params:
		biosample=lambda wc: BIOSAMPLE[wc.strain],
		assembly=lambda wc: ASSEMBLY[wc.strain],
		species=lambda wc: SPECIES[wc.strain]

	log:
		"logs/download/{strain}.log"

	shell:
		r"""
		mkdir -p results/downloaded_genomes
		mkdir -p results/download

		python scripts/fetch_assembly.py \
			--biosample "{params.biosample}" \
			--assembly "{params.assembly}" \
			--species "{params.species}" \
			--output {output.assembly} \
			--status {output.status} \
			> {log} 2>&1
		"""


##########################################################
# 2. Download reads if needed
##########################################################

rule fastq_dl:

	input:
		json="results/download/{strain}.json"

	output:
		r1=temp("results/reads/{strain}_R1.fastq.gz"),
		r2=temp("results/reads/{strain}_R2.fastq.gz")

	params:
		biosample=lambda wc: BIOSAMPLE[wc.strain]

	run:

		with open(input.json) as f:
			info = json.load(f)

		if info["status"] != "missing":

			shell("""
				touch {output.r1}
				touch {output.r2}
			""")

			return

		shell(f"""
			fastq-dl \
				-a {params.biosample} \
				--group-by-sample \
				--prefix {wildcards.strain} \
				--outdir results/reads \
				--cpus {threads}
		""")


##########################################################
# 3. Assemble only missing genomes
##########################################################

rule spades:

	input:
		json="results/download/{strain}.json",
		r1="results/reads/{strain}_R1.fastq.gz",
		r2="results/reads/{strain}_R2.fastq.gz"

	output:
		assembly=temp("results/assembled_genomes/{strain}.fna.gz")

	params:
		outdir=lambda wc: f"results/spades/{wc.strain}"

	run:

		with open(input.json) as f:
			info = json.load(f)

		if info["status"] != "missing":

			shell("""
				touch {output.assembly}
			""")

			return

		shell(f"""
			spades.py \
				-1 {input.r1} \
				-2 {input.r2} \
				-o {params.outdir}

			gzip -c {params.outdir}/contigs.fasta \
				> {output.assembly}
		""")


##########################################################
# 4. Produce ONE final genome
##########################################################

rule genome:

	input:
		json="results/download/{strain}.json",
		downloaded="results/downloaded_genomes/{strain}.fna.gz",
		assembled="results/assembled_genomes/{strain}.fna.gz"

	output:
		"results/genomes/{strain}.fna.gz"

	run:

		with open(input.json) as f:
			info = json.load(f)

		if info["status"] == "downloaded":

			shell("""
				cp {input.downloaded} {output}
			""")

		elif info["status"] == "assembled":

			shell("""
				cp {input.assembled} {output}
			""")

		else:
			raise ValueError(info)