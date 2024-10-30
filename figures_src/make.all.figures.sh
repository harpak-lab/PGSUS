#!/bin/bash
#main text figure generation
python FigureDriver.py -s main -a alpha_pillars,decomps,stratification_comparison
#supplemental figures and figures for dataverse
#make the figures in the style of figure 1 for each of the analyses
python FigureDriver.py -s supp -f decomps -a wc,nopcs,bolt,giant,1kg.pcs.only,ukb.and.1kg.pcs
python FigureDriver.py -s supp -f decomps_nondirect -a wc,nopcs,bolt,giant,1kg.pcs.only,ukb.and.1kg.pcs,okbay2022,giant2022
#make the figures in the style of main text figure 2 for each analyses
python FigureDriver.py -s supp -f alpha_pillars -a wc,nopcs,bolt,1kg.pcs.only,ukb.and.1kg.pcs
#make all versions of figure 3
python FigureDriver.py -s supp -f stratification_comparison
#make thresholding ascertainment plots for each trait
python FigureDriver.py -s supp -f thresholding
#compare alpha estimates between different methods for correction
python FigureDriver.py -s supp -f method_scatters -a wc_v_nopcs,wc_v_bolt,bolt_v_bolt
#makes supplemental files 1-10 in the form of Excel sheets
python FigureDriver.py -s supp -f trait_components -a wc,nopcs,bolt,1kg.pcs.only,ukb.and.1kg.pcs
#generate other biplots of interest that are included as supplemental figures
python FigureDriver.py -s supp -f external_correlates -a pc1,pc_loadings,ses_h2,perm_v_normal_ses,component_scatters
#compare the alpha estimates in each standard GWAS approach across thresholds
python FigureDriver.py -s supp -f threshold_scatters -a wc,nopcs,bolt,1kg.pcs.only,ukb.and.1kg.pcs
#generate simulation figures
python FigureDriver.py -s supp -f simulations -a null,thresholds,all_loci,component_error