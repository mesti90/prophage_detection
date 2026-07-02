#!/usr/bin/python3

import pandas as pd
from pathlib import Path

df = pd.read_csv("biosamples.tsv", sep="\t")
df['assembly_path'] = df['Sample_ID'].apply(lambda x: f"genomes/{x}.fna.gz")


df2 = pd.read_csv("dorina_genomes.20260604.tsv", sep="\t")

df = pd.concat([df, df2], ignore_index=True, sort=False, axis=0, join="outer")

df['assembly_exists'] = df['assembly_path'].apply(lambda x: Path(x).exists())

df.drop(columns=['Biosample'], inplace=True)


df.to_csv("genomad_input.tsv", sep="\t", index=False)
