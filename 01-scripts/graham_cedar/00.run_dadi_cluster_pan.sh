#!/bin/bash
#SBATCH --account=your_account
#SBATCH --time=99:99:99
#SBATCH --job-name=abc
#SBATCH --output=abc-%J.out
##SBATCH --array=1-33%33
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=32
##SBATCH --gres=cpu:32
##SBATCH --mail-user=yourmail
##PBS -l nodes=1:ppn=8
##SBATCH --mail-type=EA 

#eventually load python module
#module load python/2.7.14
#activate your environment if necessary
#source /home/path_to_activate_local_env/activate
#source temp/bin/activate
# Move to directory where job was submitted
cd $SLURM_SUBMIT_DIR

#module load parallel/20170622 #useless
#List of wanted arguments
fs=$1     #dataset fs stored in 03-data
model=$2  #model names
folded=$3 #either "folded" or "unfolded"
grid_size=$4

if [[ -z "$fs" ]]
then
    echo "Error: need input file"
    exit
fi

if [[ -z "$model" ]]
then
    echo "Error: need model name ex: SC2M"
    exit
fi

if [[ -z "$folded" ]]
then
    echo "Error: need folding info"
    exit
fi

if [[ -z "$grid_size" ]]
then
    echo "Error: need gride_size for optimisation"
    exit
fi

./01-scripts/00.pan.sh "$fs" "$model" "$folded" "$grid_size"
