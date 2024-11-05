#!/bin/bash

#sbatch Lai_deduper.sh "STL96.txt" input_sam "Dedup_Outputs/filtered_sorted.sam" 
#sbatch Lai_deduper.sh "STL96.txt" "/projects/bgmp/shared/deduper/C1_SE_uniqAlign.sam" "Dedup_Outputs/filtered_sorted.sam" 

#SBATCH --job-name=Deduper_KL
#SBATCH --account=bgmp
#SBATCH --partition=bgmp
#SBATCH --cpus-per-task=8
#SBATCH --mem=100G
#SBATCH --nodes=1
#SBATCH --output=/home/kenlai/bgmp/bioinfo/Bi624/Assignments/Deduper-klai22/o_e_files/slurm-special-name-%j.out
#SBATCH --error=/home/kenlai/bgmp/bioinfo/Bi624/Assignments/Deduper-klai22/o_e_files/slurm-special-name-%j.err

#MANAGING VARIABLES / INPUTS / OUTPUTS 
known_UMIs=$1 #path to .txt of known UMIs
input_sam=$2 #path to INPUT sam. file (pre-dedup.) 
o=$3 #path to OUTPUT sorted .sam (post-dedup.)

#SORTING SAM FILE 
conda activate bgmp_samtools #activate conda envs w/ samtools 
samtools sort -o /home/kenlai/bgmp/bioinfo/Bi624/Assignments/Deduper-klai22/Dedup_Outputs/sorted_input.sam --output-fmt sam $input_sam

#NAVIGATING TO WD w/ SCRIPT(s) 
cd /home/kenlai/bgmp/bioinfo/Bi624/Assignments/Deduper-klai22

#Running Deduper Script 
/usr/bin/time -v ./Lai_deduper.py -u $known_UMIs -f "/home/kenlai/bgmp/bioinfo/Bi624/Assignments/Deduper-klai22/Dedup_Outputs/sorted_input.sam" -o $o

exit
