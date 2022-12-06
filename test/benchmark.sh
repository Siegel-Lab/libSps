#!/bin/bash
#SBATCH -p slim16 --mem 100G -J benchmark_libSps --time=240:00:00 -o slurm_benchmark_libSps-%j.out --mail-user=markus.rainer.schmidt@gmail.com --mail-type END

#source activate libSps

python3 -u test/benchmark.py > test/benchmark_results.tsv