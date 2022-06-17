#!/bin/bash
#
#SBATCH -p long
#SBATCH --mem=64G
#SBATCH --job-name=BMM_mods_rate_summary
#SBATCH -o /mnt/shared/projects/nhm/goswamilab/BT_placental/tree1_75_80_95_all/reports/BMM_mods_rate_summary_%A.out
#SBATCH -e /mnt/shared/projects/nhm/goswamilab/BT_placental/tree1_75_80_95_all/reports/BMM_mods_rate_summary_%A.err


cd /mnt/shared/home/anjgoswami/projects/nhm/goswamilab/BT_placental/tree1_75_80_95_all/scripts

Rscript --no-save /mnt/shared/home/agoswami/projects/nhm/goswamilab/BT_placental/tree1_75_80_95_all/scripts/bt_summarizer_mods_BMM.R
