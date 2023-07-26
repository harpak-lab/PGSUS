import pandas as pd 
import numpy as np 
import sys
import matplotlib.pyplot as plt
from matplotlib import gridspec
import sys
import seaborn as sns
import os


df = pd.read_csv('gamma.estimate.mat.txt',sep = '\t').set_index('Unnamed: 0')
dfplot = float(1)/df.astype(float).sort_values(by='1e-08')
dfplot['ycoordinate'] = [i for i in range(dfplot.shape[0],0,-1)]
fig, ax = plt.subplots(nrows = 1, ncols = 4, figsize = (12,4))

for i,j in enumerate(df.columns):
	
	ax[i].vlines(1, ymin = 0, ymax =20, linestyles='dashed')
	ax[i].scatter(dfplot[j],dfplot['ycoordinate'], color = 'black')
	ax[i].title.set_text('p < ' + str(j))
	ax[i].set_axisbelow(True)
	ax[i].grid(True)
	if i == 0:
		ax[i].set_yticks(dfplot['ycoordinate'])
		ax[i].set_yticklabels(dfplot.index.tolist())
	else:
		ax[i].set_yticks([x for x in range(dfplot.shape[0])])
		ax[i].set_yticklabels(['' for x in range(dfplot.shape[0])])


sns.despine()
plt.tight_layout()
plt.savefig('fig1.replicate.pdf')
