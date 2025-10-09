import pandas as pd 
import numpy as np 
import sys
import matplotlib.pyplot as plt
from matplotlib import gridspec
import sys
import seaborn as sns
import os
import seaborn as sns
import statsmodels.api as sm
import warnings
import argparse
from scipy import stats
from scipy.stats import pearsonr
import matplotlib as mpl
warnings.filterwarnings('ignore')

#load classes from other scripts
from alphapillars import alphapillars
from mainfigs import main_figures
from method_scatters import method_scatters
from decomps import decomps
from decomps_nondirect import decomps_nondirect
from threshold_scatters import threshold_scatters
from stratification_comparison import stratification_comparison
from trait_component_tables import trait_components
from trait_component_tables_nondirect import trait_components_nondirect
from simulation_figures import simulation_plots
from threshold_ascertainment import thresholding_ascertainment
from external_correlates import external_correlates

parser = argparse.ArgumentParser()
parser.add_argument('-f','--figures',type=str,help='figures to make')
parser.add_argument('-a','--analyses',type=str,help='analyses to make figures')
parser.add_argument('-s','--fig-set',type=str,help='figure set to make (main or supp)',dest = 'figset')

args = parser.parse_args()
if args.figures:
	figures = [x for x in args.figures.split(",")]
if args.analyses:
	analyses = [x for x in args.analyses.split(",")]
else:
	analyses = []
figset = args.figset

label_dict = {'pack_years_of_smoking':'Pack years of smoking',\
				'pack_years_smoking':'Pack years of smoking',\
				'neuroticism_score':'Neuroticism score',\
				'alcohol_intake_frequency':'Alcohol intake frequency',\
				'alcohol_intake_freq':'Alcohol intake frequency',\
				'overall_health_rating':'Overall health rating',\
				'overall_health':'Overall health rating',\
				'years_of_schooling':'Years of schooling',\
				'years_schooling':'Years of schooling',\
				'fluid_intelligence':'Fluid intelligence',\
				'forced_vital_capacity':'Forced vital capacity',\
				'fvc':'Forced vital capacity',\
				'birth_weight':'Birth weight',\
				'diastolic_blood_pressure':'Diastolic blood pressure',\
				'dbp':'Diastolic blood pressure',\
				'bmi':'Body mass index',\
				'pulse_rate':'Pulse rate',\
				'waist_circumference':'Waist circumference',\
				'waist_circ':'Waist circumference',\
				'hand_grip_strength':'Hand grip strength',\
				'hair_color':'Hair color',\
				'hip_circumference':'Hip circumference',\
				'hip_circ':'Hip circumference',\
				'skin_color':'Skin color',\
				'household_income':'Household income',\
				'basal_metabolic_rate':'Basal metabolic rate',\
				'height':'Height',
				'giant_height':'GIANT height',
				'giant_height_rescaled':'GIANT height rescaled',
				'anthropometric':'Anthropometric',
				'behavioral':'Behavioral',
				'other':'Other',
				'all':'Median across traits',
				'null':''}

sps_traits = ['alcohol_intake_freq','birth_weight','bmi','dbp','fvc','hair_color',
		'hand_grip_strength','height','hip_circ','household_income','neuroticism_score',
		'overall_health','pack_years_smoking','pulse_rate','skin_color','waist_circ','years_schooling']

pval_array = [1e-8,1e-5,0.001,1.0]

if 'main' in figset:
	print('here')
	x = main_figures(analyses, label_dict)
	x.run()

if 'supp' in figset:
	if 'alpha_pillars' in figures:
		plots = alphapillars(analyses, label_dict)
		plots.run()

	if 'method_scatters' in figures:
		plots = method_scatters(analyses,label_dict)
		plots.run()

	if 'threshold_scatters' in figures:
		plots = threshold_scatters(analyses,label_dict)
		plots.run()

	if 'decomps' in figures:
		plots = decomps(analyses, label_dict,sps_traits)
		plots.run()

	if 'decomps_nondirect' in figures:
		plots = decomps_nondirect(analyses, label_dict,sps_traits)
		plots.run()

	if 'trait_components' in figures:
		plots = trait_components(analyses,label_dict,sps_traits)
		plots.run()

	if 'trait_components_nondirect' in figures:
		plots = trait_components_nondirect(analyses,label_dict,sps_traits)
		plots.run()

	if 'simulations' in figures:
		plots = simulation_plots(analyses,label_dict)
		plots.run()

	if 'thresholding_ascertainment' in figures:
		plots = thresholding_ascertainment(label_dict, sps_traits)
		plots.run()

	if 'external_correlates' in figures:
		plots = external_correlates(analyses, label_dict)
		plots.run()

	if 'stratification_comparison' in figures:
		plots = stratification_comparison(label_dict)
		plots.run()


