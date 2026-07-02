#!/usr/bin/env python3

import pandas as pd
import subprocess
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor
import argparse

def parse_args():
	parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
	parser.add_argument("-i","--input", type=str, default="assemblies.tsv", help="Input TSV file with Sample_ID and assembly_path columns")
	parser.add_argument("--threads", type=int, default=5, help="Number of parallel geNomad runs")
	parser.add_argument("--subthreads", type=int, default=3, help="Threads per geNomad run")
	parser.add_argument("--container", type=str, default="/home/vasarhelyib/containers/antoniopcamargo-genomad.1.8.1.sif", help="Path to geNomad container")
	parser.add_argument("-a","--assemblies", type=str, default="assemblies_for_genomad.tsv", help="TSV file with existing assemblies for geNomad")
	parser.add_argument("-o","--output", type=str, default="genomad.tsv", help="Output TSV file for geNomad results")
	parser.add_argument("--db", type=str, default="/node8_R10/kintses_lab/databases/genomad/genomad_v1.7", help="Path to geNomad database")
	return parser.parse_args()


def get_existing_assemblies(table, outtable):
	if Path(table).exists():
		df = pd.read_csv(table, sep="\t")
		df['assembly_exists'] = df['Assembly'].apply(lambda x: Path(x).exists())
		existing_assemblies = df[df["assembly_exists"]]
		existing_assemblies.to_csv(outtable, sep="\t", index=False)
		print(f"Existing assemblies saved to {outtable}")
		return existing_assemblies
	return pd.DataFrame(columns=["Strain", "Assembly"])

def run_genomad(row):
	sample = row["Strain"]
	assembly_path = row["Assembly"]
	output = f"genomad/{sample}"
	if Path(output).exists():
		print(f"{sample} already processed")
		return

	if not Path(assembly_path).exists():
		print(f"{sample} assembly {assembly_path} does not exist")
		return
	cmd = f"singularity run -B /node8_R10,/node10_R10,/scratch,/home {args.container} end-to-end --conservative --splits 8 --cleanup --quiet -t {args.subthreads} {assembly_path} {output} {args.db}"

	print(f"Running {cmd}")
	subprocess.run(cmd, check=True, shell=True, text=True)
	print(f"Finished {sample}")


def main():
	global args
	args = parse_args()
	df = get_existing_assemblies(args.input, args.assemblies)
	print("Existing assemblies for geNomad:", len(df))
	rows = df.to_dict("records")
	Path("genomad").mkdir(exist_ok=True)
	with ProcessPoolExecutor(max_workers=args.threads) as executor:
		list(executor.map(run_genomad, rows))
	


if __name__ == "__main__":
	main()