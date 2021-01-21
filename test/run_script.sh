#! /bin/bash

echo "Submitting snakemake jobs to SLURM cluster"
snakemake --latency-wait 60 --restart-times 1 --jobs 10000 --cluster "sbatch -o {log}.slurm.out -e {log}.slurm.err -n {threads} --mem {params.memory}GB --time=00:01:00 --parsable" --cluster-status ./status.py "$@"
