#!/bin/bash
#
#SBATCH -p long
#SBATCH --array=0-203
#SBATCH -c 1
#SBATCH --mem=1G
#SBATCH --job-name=per_bone_jan
#SBATCH -o /mnt/shared/projects/nhm/goswamilab/BT_placental/tree1_75_80_95_all/reports/per_bone_%A_%a.out
#SBATCH -e /mnt/shared/projects/nhm/goswamilab/BT_placental/tree1_75_80_95_all/reports/per_bone_%A_%a.err


cd /mnt/shared/scratch/rfelice/apps/BayesTraitsMPI/

FILES=($(ls -1 /mnt/shared/projects/nhm/goswamilab/BT_placental/tree1_75_80_95_all/pPCscores/modules_jan2/*run.txt))
DATA_FILE=${FILES[$SLURM_ARRAY_TASK_ID]}
echo "My input file is $DATA_FILE"
CONTROLFILES=($(ls -1 /mnt/shared/projects/nhm/goswamilab/BT_placental/tree1_75_80_95_all/scripts/command_files/*.cmd))
CMD_FILE=${CONTROLFILES[$SLURM_ARRAY_TASK_ID]}
echo "My cmd file is $CMD_FILE"

TREE_FILE=/mnt/shared/projects/nhm/goswamilab/BT_placental/tree1_75_80_95_all/tree75_80_95.nex

./BayesTraitsV3 $TREE_FILE $DATA_FILE < $CMD_FILE
