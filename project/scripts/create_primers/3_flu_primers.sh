#!/bin/bash
#SBATCH --job-name=create_flu_primer_files      
#SBATCH --partition=batch  
#SBATCH --ntasks=1                
#SBATCH --cpus-per-task=4      
#SBATCH --mem=40gb   
#SBATCH --time=2:00:00
#SBATCH --output=/work/gene8940/lml38336/log.%j
#SBATCH --mail-user=lml38336@uga.edu              
#SBATCH --mail-type=BEGIN,END,FAIL

#set working directory
cd /work/gene8940/lml38336/project/references/flu/primers

#1
echo ">Uni12/Inf-1" >> flu_1.fas
echo "GGGGGGAGCAAAAGCAGG" >> flu_1.fas
echo ">Uni13/Inf-1" >> flu_1.fas
echo "CGGGTTATTAGTAGAAACAAGG" >> flu_1.fas

#2
echo ">Uni12/Inf-3" >> flu_2.fas
echo "GGGGGGAGCGAAAGCAGG" >> flu_2.fas
echo ">Uni13/Inf-1" >> flu_2.fas
echo "CGGGTTATTAGTAGAAACAAGG" >> flu_2.fas