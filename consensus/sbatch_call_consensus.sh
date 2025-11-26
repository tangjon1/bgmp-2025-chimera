#!/bin/bash
#SBATCH --account=bgmp
#SBATCH --partition=bgmp
#SBATCH --mem=64GB
#SBATCH --cpus-per-task=1
#SBATCH --job-name=consen1
#SBATCH --out=/projects/bgmp/shared/groups/2025/chimera/jonat/bgmp-2025-chimera/consensus/consensus/results/slurm/slurm-%j-%a.out
#SBATCH --err=/projects/bgmp/shared/groups/2025/chimera/jonat/bgmp-2025-chimera/consensus/consensus/results/slurm/slurm-%j-%a.err
#SBATCH --array=1-4

output_dir=/projects/bgmp/shared/groups/2025/chimera/jonat/bgmp-2025-chimera/consensus/consensus/results

#run_iter_num=$(ls -1 -d $output_dir/*/ | sed -r 's/.+_(.+)\//\1/' | sort -n | awk '{num=$1} END{print num+1}')
results_dir="${output_dir}/results_2"
mkdir -p $results_dir


input_dir=/gpfs/projects/bgmp/shared/groups/2025/chimera/jonat/bgmp-2025-chimera/consensus/sort/results/results_1
input_files="${input_dir}/extract_005.sorted.fastq.gz ${input_dir}/extract_010.sorted.fastq.gz ${input_dir}/extract_050.sorted.fastq.gz ${input_dir}/extract.sorted.fastq.gz"
input_file=$(echo $input_files | sed 's/ /\n/g' | awk "NR==$SLURM_ARRAY_TASK_ID {print \$0}")
input_file_basename=$(basename $input_file .sorted.fastq.gz)


echo "${SLURM_JOB_ID}: This run only uses plurality and collision margin in determining consensus (included cluster size in tsv)" > "${results_dir}/README.md"


mamba run -p "/projects/bgmp/shared/groups/2025/chimera/envs/biopython" \
    /usr/bin/time -v \
        python \
            /projects/bgmp/shared/groups/2025/chimera/jonat/bgmp-2025-chimera/consensus/call_consensus.py consensus \
                --input-file $input_file \
                --output-file "${results_dir}/${input_file_basename}.consensus.tsv"


exit
