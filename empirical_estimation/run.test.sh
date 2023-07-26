#!/bin/bash
#SBATCH -J test.sim
#SBATCH -o test.sim.o%j
#SBATCH -e test.sim.o%j
#SBATCH -p normal
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --mem=16G
#SBATCH -t 1:00:00
#SBATCH --mail-user=samuel.smith@utexas.edu
#SBATCH --mail-type=begin
#SBATCH --mail-type=end
#SBATCH --array=0-18

trait_array=(height years_of_schooling birth_weight pack_years_of_smoking hip_circumference waist_circumference household_income alcohol_intake_frequency BMI forced_vital_capacity hair_color hand_grip_strength neuroticism_score overall_health_rating pulse_rate skin_color diastolic_blood_pressure basal_metabolic_rate fluid_intelligence)
trait=${trait_array[$SLURM_ARRAY_TASK_ID]}
echo $trait
pval_array=(.00000001 .00001 .001 1)
for p_cutoff in "${pval_array[@]}"
do
python sad.py --genetic-file mh_20_stats/mh20_$trait/clumped_all_pops_1kg_phase3_autosomes_apr22_mh20_$trait"_standard.nodups.bed"  \
	--pop-gwas mh_20_stats/mh20_$trait/preproc_mh20_$trait"_standard.tsv.gz" \
	--sib-gwas mh_20_stats/mh20_$trait/preproc_mh20_$trait"_sib.tsv.gz" \
	--pvalue $p_cutoff \
	--pval-col pval \
	--ascertainment-set gwas \
	--out mh_20_stats/mh20_$trait/

# Rscript decomp_argparse.R mh_20_stats/mh20_$trait/clumped_eur_1000G_mh20_$trait"_standard.nodups" \
# 			mh_20_stats/mh20_$trait/preproc_mh20_$trait"_standard.tsv.gz" \
# 			mh_20_stats/mh20_$trait/preproc_mh20_$trait"_sib.tsv.gz" \
# 			$p_cutoff \
# 			gwas \
# 			mh_20_stats/mh20_$trait/
done



