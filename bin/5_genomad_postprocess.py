#!/usr/bin/env python3

import os
import pdb
import glob
import subprocess
import pandas as pd
from pathlib import Path
from scipy.spatial.distance import pdist, squareform
import shutil
import argparse
from concurrent.futures import ThreadPoolExecutor, as_completed
from typing import NamedTuple
from dataclasses import dataclass, field
from collections import defaultdict
import csv
from Bio import SeqIO

# ========================
# --- Argument Parsing ---
# ========================

def parse_args():
	parser = argparse.ArgumentParser(description="Pipeline for plasmid/prophage clustering and Mantel analysis.")
	parser.add_argument("--table", default="assemblies.tsv", help="Input TSV file with Sample_ID and assembly_path columns")
	parser.add_argument("--genomad_path", default="genomad", help="Path to genomad results directory")
	parser.add_argument("--provirus_dir", default="genomad/provirus_sequences/fna")
	parser.add_argument("--plasmid_dir", default="genomad/plasmid_sequences/fna")
	parser.add_argument("--filtered_provirus_dir", default="genomad/filtered_provirus_sequences")
	parser.add_argument("--checkv_db", default="/node8_R10/kintses_lab/databases/checkv/checkv-db-v1.5")
	parser.add_argument("--cluster_dir", default="clusters")
	parser.add_argument("--checkv_sif", default="/home/vasarhelyib/containers/staphb-checkv.1.0.3.sif")
	parser.add_argument("--seqtk_sif", default="/home/vasarhelyib/containers/staphb-seqtk.1.4.sif")
	parser.add_argument("--mmseqs_sif", default="/home/vasarhelyib/containers/soedinglab-mmseqs2.version-13.sif")
	parser.add_argument("--r_sif", default="/home/vasarhelyib/containers/mesti90-r-vegan.1.0.sif")
	parser.add_argument("--out_prefix", default="genomad_results")
	parser.add_argument("--out_dir", default="genomad/collected_genomad_results")
	parser.add_argument("--checkv_out", default="checkv_outputs")
	parser.add_argument("--checkv_threads", type=int, default=4)
	parser.add_argument("-n", "--threads", type=int, default=20)
	return parser.parse_args()

# ========================
# --- Utility Functions --
# ========================

def concat_sequences(file_list, prefix_list, output_file):
	print(f"Concatting to {output_file}")
	with open(output_file, "w") as out:
		for fname,prefix in zip(file_list,prefix_list):
			for rcd in SeqIO.parse(fname,"fasta"):
				rcd.id = f"{prefix}|{rcd.id}"
				SeqIO.write(rcd, out, "fasta")

def build_matrix(cluster_file, cluster_type):
	if not os.path.exists(cluster_file):
		raise FileNotFoundError(f"[ERROR] Expected cluster file not found: {cluster_file}")
	df = pd.read_csv(cluster_file, sep="\t", header=None, names=["cluster", "member"])
	matrix = {}
	for _, row in df.iterrows():
		sample = row["member"].split("_")[0]
		cluster_id = row["cluster"]
		matrix.setdefault(sample, set()).add(cluster_id)

	all_clusters = sorted({c for clusters in matrix.values() for c in clusters})
	prefix = "Provirus" if cluster_type == "provirus" else "Plasmid"
	cluster_map = {old: f"{prefix}_{i+1:05d}" for i, old in enumerate(all_clusters)}
	renamed_clusters = [cluster_map[c] for c in all_clusters]

	df_binary = pd.DataFrame(0, index=matrix.keys(), columns=renamed_clusters)
	for sample, clusters in matrix.items():
		df_binary.loc[sample, [cluster_map[c] for c in clusters]] = 1
	return df_binary

def write_mantel_r_script():
	r_script = """
library(vegan)

args <- commandArgs(trailingOnly=TRUE)
dist_genetic_file <- args[1]
dist_phylo_file <- args[2]
label <- args[3]
out_file <- args[4]

dist_genetic <- as.dist(read.table(dist_genetic_file, header=TRUE, row.names=1, check.names=FALSE))
dist_phylo <- as.dist(read.table(dist_phylo_file, header=TRUE, row.names=1, check.names=FALSE))

res <- mantel.partial(dist_genetic, dist_phylo, dist_phylo, method="spearman")

sink(out_file)
cat(paste("Partial Mantel test for", label, "\n"))
cat("Rho:", res$statistic, "\n")
cat("p-value:", res$signif, "\n")
sink()
"""
	with open("run_mantel.R", "w") as f:
		f.write(r_script)

# =============================
# --- Core Pipeline Steps -----
# =============================

def run_checkv_and_filter(checkv_args, args):
	input_fna = checkv_args.input_fna
	output_dir = checkv_args.output_dir
	summary_file = checkv_args.summary_file
	sample = checkv_args.sample
	out_fna = checkv_args.out_fna
	if not os.path.exists(summary_file):
		os.makedirs(output_dir, exist_ok=True)
		cmd = [
			"singularity", "run", "-B", "/node8_R10,/node10_R10",args.checkv_sif,
			"checkv", "end_to_end", input_fna, output_dir,
			"-t", str(args.checkv_threads), "-d", args.checkv_db
		]
		subprocess.run(cmd, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
	empty = False
	if not os.path.exists(summary_file):
		empty = True
	else:
		df = pd.read_csv(summary_file, sep="\t")
		high_quality = df[df["completeness"] >= 90]["contig_id"].tolist()
		if not high_quality:
			empty = True
		else:
			tmp_ids = f"{sample}_ids.txt"
			with open(tmp_ids, "w") as f:
				f.write("\n".join(high_quality) + "\n")
			with open(out_fna, "w") as fout:
				subprocess.run([
					"singularity", "exec", args.seqtk_sif,
					"seqtk", "subseq", input_fna, tmp_ids
				], stdout=fout, stderr=subprocess.DEVNULL)
			os.remove(tmp_ids)
	
	if empty:
		Path(out_fna).touch()


@dataclass
class MmseqsArgs:
	input_fasta: str = ""
	out_prefix: str = ""
	cluster_dir: str = ""
	
	def __post_init__(self):
		self.cluster_path = os.path.join(self.cluster_dir, self.out_prefix)
		self.cluster_out = os.path.join(self.cluster_path, "clusters")
		self.tmp_dir = os.path.join(self.cluster_path, "tmp")
		self.output_file = os.path.join(self.cluster_path, "clusters_cluster.tsv")
		self.rep_file = os.path.join(self.cluster_path, "clusters_rep_seq.fasta")
		self.renamed_seq = os.path.join(self.cluster_dir, f"{self.out_prefix}.rep_sequences.fasta")
		self.presence = os.path.join(self.cluster_dir, f"{self.out_prefix}.presence.tsv")
		os.makedirs(self.cluster_path, exist_ok=True)
		os.makedirs(self.tmp_dir, exist_ok=True)
		
def run_mmseqs(mmseqs_args, args):
	if os.path.exists(mmseqs_args.output_file):
		print(f"[SKIP] MMseqs2 output already exists: {mmseqs_args.output_file}")
		return
	cmd = [
		"singularity", "exec","-B","/node8_R10,/node10_R10",
		args.mmseqs_sif,
		"mmseqs", "easy-cluster",
		mmseqs_args.input_fasta, mmseqs_args.cluster_out, mmseqs_args.tmp_dir,
		"--min-seq-id", "0.9", "-c", "0.9", "--cov-mode", "1"
	]
	print(" ".join(cmd))
	subprocess.run(cmd, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
	if not os.path.exists(mmseqs_args.output_file):
		print(f"[ERROR] MMseqs2 did not produce expected file: {mmseqs_args.output_file}")
	else:
		print(f"[OK] MMseqs2 output ready: {mmseqs_args.output_file}")

def run_r_mantel(jaccard_file, phylo_file, label, output_file, args):
	subprocess.run([
		"singularity", "exec","-B","/node8_R10,/node10_R10",
		args.r_sif,
		"Rscript", "run_mantel.R", jaccard_file, phylo_file, label, output_file,
	], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)



@dataclass
class CheckvArgs:
	input_fna: str = ""
	checkv_out: str = ""
	sample: str = field(init=False)
	output_dir: str = field(init=False)
	summary_file: str = field(init=False)
	out_fna: str = field(init=False)
	filtered_provirus_dir: str = ""
	
	def __post_init__(self):
		self.sample = Path(self.input_fna).stem
		self.output_dir = os.path.join(self.checkv_out, self.sample)
		self.summary_file = os.path.join(self.output_dir, "quality_summary.tsv")
		self.out_fna = os.path.join(self.filtered_provirus_dir, f"{self.sample}_filtered.fna")


def run_checkv(fna_files, checkv_args_list, args):
	if not fna_files:
		print("[INFO] No provirus .fna files found. Skipping CheckV and continuing with empty dataset.")
		return
	print(f"[INFO] Running CheckV on {len(fna_files)} files...",flush=False)
	with ThreadPoolExecutor(max_workers=args.threads) as executor:
		futures = {executor.submit(run_checkv_and_filter, checkv_args, args): checkv_args for checkv_args in checkv_args_list}
		total = len(futures)
		for i, future in enumerate(as_completed(futures), 1):
			result = future.result()
			print(f"\r[CheckV] Completed {i}/{total} (running: {total - i})", end="", flush=True)
	print()  # newline after status


def create_provirus_fasta(checkv_args_list, args):
	fna_list = [checkv_args.out_fna for checkv_args in checkv_args_list]
	prefix_list = [checkv_args.sample for checkv_args in checkv_args_list]
	concat_sequences(fna_list, prefix_list, args.provirus_fasta)

def create_plasmid_fasta(args):
	fna_list = glob.glob(os.path.join(args.plasmid_dir, "*.fna"))
	prefix_list = [fna.replace("_plasmid.fna","") for fna in fna_list]
	concat_sequences(fna_list, prefix_list, args.plasmid_fasta)

def build_presence_matrix(provirus_fasta: str, plasmid_fasta: str, args):
	print("\n[INFO] Clustering provirus sequences with MMseqs2...")
	provirus_clusters = run_mmseqs(provirus_fasta, "provirus", args)
	print("\n[INFO] Clustering plasmid sequences with MMseqs2...")
	plasmid_clusters = run_mmseqs(plasmid_fasta, "plasmid", args)
	#df_provirus = build_matrix(provirus_clusters, "provirus")
	#df_plasmid = build_matrix(plasmid_clusters, "plasmid")
	#df_combined = pd.concat([df_provirus, df_plasmid], axis=1).fillna(0).astype(int)
	#df_combined.to_csv(os.path.join(args.out_dir, "presence_absence_matrix.tsv"), sep="\t")
	#return df_combined

def compute_distances_and_mantel(df_combined: pd.DataFrame, args):
	jaccard = pd.DataFrame(
		squareform(pdist(df_combined, metric="jaccard")),
		index=df_combined.index,
		columns=df_combined.index
	)
	jaccard_file = os.path.join(args.out_dir, "jaccard_matrix.tsv")
	jaccard.to_csv(jaccard_file, sep="\t")
	phylo = pd.DataFrame(
		squareform(pdist(df_combined, metric="euclidean")),
		index=df_combined.index,
		columns=df_combined.index
	)
	phylo_file = os.path.join(args.out_dir, "phylogenetic_distance.tsv")
	phylo.to_csv(phylo_file, sep="\t")
	print("[INFO] Computing distance matrices and running Mantel test...")
	write_mantel_r_script()
	mantel_out = os.path.join(args.out_dir, f"mantel_result_{args.out_prefix}.txt")
	run_r_mantel(jaccard_file, phylo_file, args.out_prefix, mantel_out, args)


def mmseqs_to_prevalence(mmseqs_args):
	record_cluster = {}
	clid = 0
	clprefix = f"{mmseqs_args.out_prefix}_"
	clusters = []
	with open(mmseqs_args.renamed_seq,"w") as g:
		for rcd in SeqIO.parse(mmseqs_args.rep_file,"fasta"):
			clid += 1
			clname = f"{clprefix}{clid:05}"
			clusters.append(clname)
			record_cluster[rcd.id] = clname
			g.write(f">{clname}\n{rcd.seq}\n")
	cluster_members = defaultdict(lambda: defaultdict(lambda: 0))
	genomes = set()
	with open(mmseqs_args.output_file) as f:
		rdr = csv.reader(f, delimiter="\t")
		for row in rdr:
			cluster = record_cluster[row[0]]
			genome = row[1].split("|")[0]
			cluster_members[genome][cluster] = 1
			genomes.add(genome)
	genomes = sorted(set(genomes))
	with open(mmseqs_args.presence, "w") as g:
		wtr = csv.writer(g, delimiter="\t")
		wtr.writerow(["Genome"] + clusters + ["Total"])
		cluster_totals = [0] * len(clusters)
		for genome in genomes:
			genome_data = cluster_members[genome]
			row_vals = [genome_data[cluster] for cluster in clusters]
			row_total = sum(row_vals)
			for i, val in enumerate(row_vals):
				cluster_totals[i] += val
			wtr.writerow([genome] + row_vals + [row_total])
		# Add "Total" row
		wtr.writerow(["Total"] + cluster_totals + [sum(cluster_totals)])
	print(f"Ready : {mmseqs_args.presence}")




def run_and_process(item, mmseqs_args, args):
	print(f"\n[INFO] Clustering {item} sequences with MMseqs2...")
	run_mmseqs(mmseqs_args[item], args)
	return item

def parallel_run_mmseqs(mmseqs_args, args):
	with ThreadPoolExecutor(max_workers=args.threads) as executor:
		futures = {executor.submit(run_and_process, item, mmseqs_args, args): item for item in mmseqs_args}
		for future in as_completed(futures):
			item = futures[future]
			try:
				result = future.result()
				# Optionally: mmseqs_to_prevalence could go here if needed
				# mmseqs_to_prevalence(result)
			except Exception as e:
				print(f"[ERROR] MMseqs2 clustering failed for {item}: {e}")

def get_assemblies(args):
	if Path(args.table).exists():
		df = pd.read_csv(args.table, sep="\t")
		df['genomad_ok'] = df['Strain'].apply(lambda x: (Path(args.genomad_path) / x).exists())

		assemblies = df[df["genomad_ok"]]['Strain'].tolist()
		
		return assemblies
	return []

def collect_genomad_results(assemblies, args):
	for sample in assemblies:
		sample_dir = Path(args.genomad_path) / sample / f"{sample}_summary"
		plasmid_file = sample_dir / f"{sample}_plasmid.fna"
		virus_file = sample_dir / f"{sample}_virus.fna"
		if plasmid_file.exists():
			shutil.copy(plasmid_file, Path(args.plasmid_dir) / f"{sample}_plasmid.fna")
		if virus_file.exists():
			shutil.copy(virus_file, Path(args.provirus_dir) / f"{sample}_virus.fna")

	print(f"{len(assemblies)} Plasmid/virus sequences copied")


def main():
	args = parse_args()
	os.makedirs(args.out_dir, exist_ok=True)
	args.filtered_provirus_dir = os.path.join(args.out_dir, args.filtered_provirus_dir)
	args.checkv_out = os.path.join(args.out_dir, args.checkv_out)
	args.cluster_dir = os.path.join(args.out_dir, args.cluster_dir)
	
	args.provirus_fasta = os.path.join(args.out_dir, "all_prophages.fasta")
	args.plasmid_fasta = os.path.join(args.out_dir, "all_plasmids.fasta")
	
	for d in [args.plasmid_dir, args.provirus_dir, args.filtered_provirus_dir, args.checkv_out, args.cluster_dir]:
		Path(d).mkdir(parents=True, exist_ok=True)
	
	#0. copy virus and plasmid sequences from genomad results
	assemblies = get_assemblies(args)
	collect_genomad_results(assemblies, args)

	#1. Run checkV for provirus & filter provirus using checkV completeness & concat files
	if os.path.exists(args.provirus_fasta):
		print(f"[SKIP] Found existing {args.provirus_fasta}, skipping CheckV processing.")
	else:
		fna_files = glob.glob(os.path.join(args.provirus_dir, "*.fna"))
		checkv_args_list = [CheckvArgs(input_fna = fna, checkv_out = args.checkv_out, filtered_provirus_dir = args.filtered_provirus_dir,) for fna in fna_files]
		run_checkv(fna_files, checkv_args_list, args)
		create_provirus_fasta(checkv_args_list, args)
	#2. Concat plasmid fasta
	if os.path.exists(args.plasmid_fasta):
		print(f"[SKIP] Found existing {args.plasmid_fasta}, skipping concatenation.")
	else:
		create_plasmid_fasta(args)
	
	#3. mmseqs for provirus & plasmid
	mmseqs_args = {
		"provirus": MmseqsArgs(args.provirus_fasta, "provirus", args.cluster_dir),
		"plasmid": MmseqsArgs(args.plasmid_fasta, "plasmid", args.cluster_dir),
	}
	parallel_run_mmseqs(mmseqs_args, args)
	
	for key in mmseqs_args:
		mmseqs_to_prevalence(mmseqs_args[key])

	
	
	#if not filtered:
	#	df_matrix = pd.DataFrame()
	#else:
	#	df_matrix = build_presence_matrix(provirus_fasta, plasmid_fasta, args)
	#if not df_matrix.empty:
	#	compute_distances_and_mantel(df_matrix, args)
	#else:
	#	print("[INFO] No provirus data to process for clustering and Mantel test.")
	#print("[✅] Pipeline complete. Outputs saved.")

if __name__ == "__main__":
	main()
