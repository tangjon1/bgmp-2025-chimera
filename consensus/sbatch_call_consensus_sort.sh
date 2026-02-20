#!/bin/bash
#SBATCH --account=bgmp
#SBATCH --partition=bgmp
#SBATCH --mem=32GB
#SBATCH --cpus-per-task=8
#SBATCH --job-name=sort
#SBATCH --out=/projects/bgmp/shared/groups/2025/chimera/jonat/bgmp-2025-chimera/consensus/sort/results/slurm/slurm-%j-%a.out
#SBATCH --err=/projects/bgmp/shared/groups/2025/chimera/jonat/bgmp-2025-chimera/consensus/sort/results/slurm/slurm-%j-%a.err
#SBATCH --array=4

input_files="/projects/bgmp/shared/groups/2025/chimera/jonat/bgmp-2025-chimera/consensus/sort/input/extract_005.fastq.gz /projects/bgmp/shared/groups/2025/chimera/jonat/bgmp-2025-chimera/consensus/sort/input/extract_010.fastq.gz /projects/bgmp/shared/groups/2025/chimera/jonat/bgmp-2025-chimera/consensus/sort/input/extract_050.fastq.gz /projects/bgmp/shared/groups/2025/chimera/shared/bgmp_pacbio/02_extract/251110_output_01/extract.fastq.gz"
input_file=$(echo $input_files | sed 's/ /\n/g' | awk "NR==$SLURM_ARRAY_TASK_ID {print \$0}")
input_file_basename=$(basename $input_file .fastq.gz)

output_dir=/projects/bgmp/shared/groups/2025/chimera/jonat/bgmp-2025-chimera/consensus/sort/results

#run_iter_num=$(ls -1 -d $output_dir/*/ | sed -r 's/.+_(.+)\//\1/' | sort -n | awk '{num=$1} END{print num+1}')
#results_dir="${output_dir}/results_${run_iter_num}"
results_dir="${output_dir}/results_1"
mkdir -p $results_dir


tmp=$(mktemp)
cut -f 1,3 /projects/bgmp/shared/groups/2025/chimera/rmurnane/pacbio_stuff/starcode_for_shared_test/03_starcode/starcode_barcodes.tsv > "$tmp"


mamba run -p "/projects/bgmp/shared/groups/2025/chimera/envs/biopython" \
    /usr/bin/time -v \
        python \
            /projects/bgmp/shared/groups/2025/chimera/jonat/bgmp-2025-chimera/consensus/call_consensus.py sort \
                --input-file $input_file \
                --output-file "${results_dir}/${input_file_basename}.sorted.fastq" \
                --cluster-file $tmp \
                -t $SLURM_CPUS_PER_TASK

/usr/bin/time -v pigz -9 -p $SLURM_CPUS_PER_TASK "${results_dir}/${input_file_basename}.sorted.fastq"


exit
