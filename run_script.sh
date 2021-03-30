#! /bin/bash

echo "Submitting snakemake jobs to SLURM cluster"
snakemake --latency-wait 120 --restart-times 1 --jobs 500 --cluster "sbatch -o {log}.slurm.out -e {log}.slurm.err -n {threads} --mem {params.memory}GB --time=01:00:00 --parsable" --cluster-status ./status.py "$@"
