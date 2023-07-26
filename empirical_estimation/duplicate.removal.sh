#!/bin/bash
#SBATCH -J test.sim
#SBATCH -o test.sim.o%j
#SBATCH -e test.sim.o%j
#SBATCH -p vm-small
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -t 1:00:00
#SBATCH --mail-user=samuel.smith@utexas.edu
#SBATCH --mail-type=begin
#SBATCH --mail-type=end


trait_array=(height years_of_schooling birth_weight pack_years_of_smoking hip_circumference waist_circumference household_income alcohol_intake_frequency basal_metabolic_rate BMI fluid_intelligence forced_vital_capacity hair_color hand_grip_strength neuroticism_score overall_health_rating pulse_rate skin_color diastolic_blood_pressure)

# for i in ${trait_array[@]}
# do
# 	awk '{print $2}' mh_20_stats/mh20_$i/clumped_eur_1000G_mh20_$i"_standard.bim" | uniq -c | awk '$1!=1{print $2}' > mh_20_stats/mh20_$i/duplicate.ids.txt
# 	plink --bfile mh_20_stats/mh20_$i/clumped_eur_1000G_mh20_$i"_standard" --exclude mh_20_stats/mh20_$i/duplicate.ids.txt --make-bed --out mh_20_stats/mh20_$i/clumped_eur_1000G_mh20_$i'_standard.nodups'
# done

for i in ${trait_array[@]}
do
        awk '{print $2}' mh_20_stats/mh20_$i/clumped_all_pops_1kg_phase3_autosomes_apr22_mh20_$i"_standard.bim" | uniq -c | awk '$1!=1{print $2}' > mh_20_stats/mh20_$i/duplicate.ids.allpops.txt
        plink --bfile mh_20_stats/mh20_$i/clumped_all_pops_1kg_phase3_autosomes_apr22_mh20_$i"_standard" --exclude mh_20_stats/mh20_$i/duplicate.ids.allpops.txt --make-bed --out mh_20_stats/mh20_$i/clumped_all_pops_1kg_phase3_autosomes_apr22_mh20_$i"_standard.nodups"
done