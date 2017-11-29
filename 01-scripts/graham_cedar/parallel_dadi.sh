#!/bin/bash
#SBATCH --account=def-audetc
#SBATCH --time=02:45:00
#SBATCH --job-name=abc
#SBATCH --output=abc-%J.out
##SBATCH --array=1-33%33
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=32
##SBATCH --gres=cpu:32
##SBATCH --mail-user=quentinrougemont@orange.fr
##PBS -l nodes=1:ppn=8
##SBATCH --mail-type=EA 

source temp/bin/activate
# Move to directory where job was submitted
cd $SLURM_SUBMIT_DIR

module load parallel/20170622

cd 01-scripts

#List of wanted arguments
fs=$1     #dataset fs stored in 03-data
model=$2  #model names
folded=$3 #either "folded" or "unfolded"

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

nrep=100
NUM_CPUS=32
#seq $bootrep |parallel -j "$NUM_CPUS" ../treemix_iterations.sh {} "$mig"
#for i in $(seq 1)
#do
seq $nrep |parallel -j "$NUM_CPUS" ./00-run_model_iteration.sh {} "$fs" "$model" "$folded"
#done
