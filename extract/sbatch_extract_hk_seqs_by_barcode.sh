#!/bin/bash
#SBATCH --account=bgmp
#SBATCH --partition=bgmp
#SBATCH --mem=64GB
#SBATCH --cpus-per-task=40
#SBATCH --job-name=extract
#SBATCH --out=/projects/bgmp/shared/groups/2025/chimera/jonat/bgmp-2025-chimera/extract/results/slurm/slurm-%j.out
#SBATCH --err=/projects/bgmp/shared/groups/2025/chimera/jonat/bgmp-2025-chimera/extract/results/slurm/slurm-%j.err

output_dir=/projects/bgmp/shared/groups/2025/chimera/jonat/bgmp-2025-chimera/extract/results

run_iter_num=$(ls -1 -d $output_dir/*/ | sed -r 's/.+_(.+)\//\1/' | sort -n | awk '{num=$1} END{print num+1}')
results_dir="${output_dir}/results_${run_iter_num}"
mkdir $results_dir

mamba run -p "/projects/bgmp/shared/groups/2025/chimera/envs/biopython" \
    /usr/bin/time -v \
        python \
            /projects/bgmp/shared/groups/2025/chimera/jonat/bgmp-2025-chimera/extract/extract_hk_seqs_by_barcode.py \
                -i /projects/bgmp/shared/groups/2025/chimera/jonat/bgmp-2025-chimera/original_pacbio_pipeline/results/results_7/01_bedtools/segmented.fastq.gz \
                -o "${results_dir}/" \
                -t $SLURM_CPUS_PER_TASK \
                -q
