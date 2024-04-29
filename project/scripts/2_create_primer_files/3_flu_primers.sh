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
#modified original primers to remove overhang not found in reference genome 
echo ">MBTuni-12" >> flu_1.fas
echo "AGCAAAAGCAGG" >> flu_1.fas
echo ">MBTuni-13" >> flu_1.fas
echo "AGTAGAAACAAGG" >> flu_1.fas


