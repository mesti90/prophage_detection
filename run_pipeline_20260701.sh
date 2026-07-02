bash prophage_detection/bin/1_download_assemblies.sh -i samples.tsv -t assemblies.tsv -f failed_biosamples.txt
python3 prophage_detection/bin/2_download_sra_reads.py
snakemake --use-singularity -s prophage_detection/denovo_assembler/Snakefile --profile prophage_detection/denovo_assembler/profiles/server  --cores 8 --configfile prophage_detection/denovo_assembler/config.yaml --config samples=assembly_input_for_sra_reads.tsv
python3 prophage_detection/bin/4_batch_genomad.py
python3 prophage_detection/bin/5_genomad_postprocess.py