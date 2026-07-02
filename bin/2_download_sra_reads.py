import shutil
import argparse
import pandas as pd
from pathlib import Path
import subprocess
from concurrent.futures import ThreadPoolExecutor, as_completed

sra_container = "/home/vasarhelyib/containers/ncbi-sra-tools.3.4.1.sif"

def parse_args():
	parser = argparse.ArgumentParser( description="Find samples without assemblies and download SRA reads.", formatter_class=argparse.ArgumentDefaultsHelpFormatter )

	parser.add_argument( "-i", "--input", default="assemblies.tsv", help="Input biosamples TSV" )

	parser.add_argument( "-m", "--missing-table", default="biosample_missing_assemblies.tsv", help="Output table of samples missing assemblies" )

	parser.add_argument( "-o", "--assembly-input", default="assembly_input_for_sra_reads.tsv", help="Output table for downstream assembly generation" )

	parser.add_argument( "--failed", default="biosamples_without_assembly_or_sra.tsv", help="Output table of samples without assemblies or downloaded reads" )
	return parser.parse_args()

def get_missing_assemblies(intable, outtable):
	df = pd.read_csv(intable, sep="\t")
	df['assembly_path'] = df['Strain'].apply(lambda x: f"genomes/{x}.fna.gz")
	df['assembly_exists'] = df['assembly_path'].apply(lambda x: Path(x).exists())
	missing_assemblies = df[~df['assembly_exists']]
	missing_assemblies.to_csv(outtable, sep="\t", index=False)
	return missing_assemblies


def fastq_dl(biosample, strain, outputdir):
	r1_path = Path(outputdir) / f"{strain}_1.fastq"
	r2_path = Path(outputdir) / f"{strain}_2.fastq"

	if r1_path.exists() and r2_path.exists() and not r1_path.stat().st_size == 0 and not r2_path.stat().st_size == 0:
		print(f"{strain} reads already exist, skipping download")
		return {
			"biosample": biosample,
			"downloaded": True,
			"r1": str(r1_path),
			"r2": str(r2_path),
			"returncode": 0,
		}
	print(f"[{strain}] downloading reads from {biosample}...")
	subprocess.run(f"singularity run -B /node8_R10,/node8_data,/node10_R10,/home,/scratch {sra_container} fasterq-dump {biosample} --outfile sra_reads/{strain} --threads 8 --temp sra_reads/{strain}_temp", shell=True, text=True,)



def create_assembly_input(intable, outtable):
	df = pd.read_csv(intable, sep="\t")
	df['assembly_exists'] = df['assembly_path'].apply(lambda x: Path(x).exists())
	missing_assemblies = df[~df['assembly_exists']].copy()
	existing_sra_reads = missing_assemblies[missing_assemblies['R1'].apply(lambda x: Path(x).exists()) & missing_assemblies['R2'].apply(lambda x: Path(x).exists())].copy()
	existing_sra_reads.drop(columns=['Assembly','assembly_exists'], inplace=True)
	existing_sra_reads.columns = ["Sample_ID", "Biosample", "Species", "Assembly_Path","R1", "R2"]
	existing_sra_reads.to_csv(outtable, sep="\t", index=False)
	
	print("Assembly input table created:", outtable)

def download_sra_reads(sample_table):
	with ThreadPoolExecutor(max_workers=4) as executor:
		futures = [executor.submit(fastq_dl, row["Biosample"], row["Strain"], "sra_reads") for _, row in sample_table.iterrows()]
		for future in as_completed(futures):
			future.result()

def check_downloaded_reads(missing, outtable, failed_outtable):	
	downloaded = []
	r1 = []
	r2 = []
	for _,row in missing.iterrows():
		sample = row['Strain']
		R1 = Path('sra_reads') / f"{sample}_1.fastq"
		R2 = Path('sra_reads') / f"{sample}_2.fastq"
		r1.append(R1)
		r2.append(R2)
		downloaded.append(R1.exists() and R2.exists())
	missing['downloaded'] = downloaded
	missing['R1'] = r1
	missing['R2'] = r2
	df = missing[missing['downloaded'] == True].copy()
	still_missing = missing[missing['downloaded'] == False].copy()
	if not df.empty:
		df.drop(columns=['downloaded', 'assembly_exists'], inplace=True)
		df.to_csv(outtable, sep="\t", index=False)
	if not still_missing.empty:
		still_missing.to_csv(failed_outtable, sep="\t", index=False)

if __name__ == "__main__":
	args = parse_args()

	missing = get_missing_assemblies( args.input, args.missing_table )
	if not missing.empty:
		download_sra_reads(missing)
		check_downloaded_reads( missing, args.assembly_input, args.failed )
		create_assembly_input(args.assembly_input, args.assembly_input)
