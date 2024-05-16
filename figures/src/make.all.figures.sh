#!/bin/bash
#main text figure generation
python FigureDriver.py -s main -f alpha_pillars, decomps, stratification_comparison
#supplemental figures and figures for dataverse
python FigureDriver.py -s supp -f alpha_pillars -a wc,nopcs,bolt,halves,nosingletons,shuffle
python FigureDriver.py -s supp -f method_scatters -a wc_v_nopcs,wc_v_nosingletons,wc_v_bolt,bolt_v_bolt,ascertainment_v_validation,wc_v_shuffle
python FigureDriver.py -s supp -f threshold_scatters -a wc,nopcs,bolt,halves,nosingletons,shuffle
python FigureDriver.py -s supp -f decomps -a wc,nopcs,bolt,halves,nosingletons,shuffle,giant,1kg.pcs.only,ukb.and.1kg.pcs
python FigureDriver.py -s supp -f trait_components -a wc,nopcs,bolt,halves,nosingletons,1kg.pcs.only,ukb.and.1kg.pcs
python FigureDriver.py -s supp -f stratification_comparison