#! /bin/bash

echo "Submitting snakemake jobs"

CLUSTER="NONE"

if [ ! -z `which sbatch` ]; then
  CLUSTER="SLURM"
fi

if [ ! -z `which bsub` ]; then
  CLUSTER="LSF"
fi

case "$CLUSTER" in
"LSF")
  echo "Submitting snakemake jobs to LSF cluster"
  snakemake --latency-wait 60 --restart-times 1 --jobs 10000 --cluster "bsub -oo {log}.bsub -n {threads} -R rusage[mem={params.memory}000] -R span[hosts=1]" "$@"
  ;;
"SLURM")
  echo "Submitting snakemake jobs to SLURM cluster"
  snakemake combine_distances --rerun-incomplete --latency-wait 60 --jobs 500 --cluster "sbatch -o {log}.slurm.out -e {log}.slurm.err -n {threads} --mem {params.memory}GB --time=04:00:00 --parsable" --cluster-status ./status.py "$@"
  ;;
*)
  snakemake --cores 8 "$@"
  ;;
esac


