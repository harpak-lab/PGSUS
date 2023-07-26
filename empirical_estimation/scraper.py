import pandas as pd 
import numpy as np

trait_array = ['height', 'years_of_schooling','birth_weight','pack_years_of_smoking',
			'hip_circumference','waist_circumference','household_income','alcohol_intake_frequency',
			'BMI','forced_vital_capacity','hair_color',
			'hand_grip_strength','neuroticism_score','overall_health_rating','pulse_rate','skin_color',
			'diastolic_blood_pressure', 'fluid_intelligence','basal_metabolic_rate']
thresholds = ['1e-08','1e-05','0.001','1.0']

outdf = pd.DataFrame(np.zeros((len(trait_array),len(thresholds))), index = trait_array, columns = thresholds)

for trait in trait_array:
	for threshold in thresholds:
		newfile = [i.strip() for i in open('mh_20_stats/mh20_' + trait + '/pval.' + threshold + '.gamma.txt','r')]
		outdf.loc[trait,threshold] = newfile[0]

outdf.to_csv('gamma.estimate.mat.txt',sep = '\t')

outdf = pd.DataFrame(np.zeros((len(trait_array),len(thresholds))), index = trait_array, columns = thresholds)

for trait in trait_array:
	for threshold in thresholds:
		newfile = [i.strip() for i in open('mh_20_stats/mh20_' + trait + '/pval.' + threshold + '.se.txt','r')]
		outdf.loc[trait,threshold] = newfile[0]

outdf.to_csv('se.estimate.mat.txt',sep = '\t')
