import numpy as np
import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)
warnings.simplefilter(action='ignore', category=RuntimeWarning)
import pandas as pd
import sys
import statsmodels.formula.api as smf
import statsmodels.api as sm
from permutations import block_permutation
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
import matplotlib.pyplot as plt
from scipy.odr import *

class estimate_components(object):

    def __init__(self, pc_genotypes, gwas_beta, gwas_se, sib_beta, sib_se, chr_pos, asc_p, thresh, outpath, outlabel, pos_label,
        pc_lower_bound=100, eigenvecs= None, eigenvalues = None, boot_se = 100, block_perm = True, pcs_to_test = 15, nperm = 1000):
        
        self.gwas_beta = gwas_beta.reset_index(drop = True)
        self.gwas_se = gwas_se.reset_index(drop = True)
        self.sib_beta = sib_beta.reset_index(drop = True)
        self.sib_se = sib_se.reset_index(drop = True)
        self.ascertainment_p = asc_p.reset_index(drop = True)
        self.thresh = thresh
        self.pc_lower_bound = pc_lower_bound
        self.outpath = outpath
        self.nboot = boot_se
        self.nperm = nperm
        self.pcs_to_test = pcs_to_test
        self.outlabel = outlabel
        self.plot = 0

        if eigenvecs == None:
            self.pca(pc_genotypes)
        else:
            self.eigenvalues = np.load(eigenvalues)
            self.eigenvecs = np.load(eigenvecs)

        self.alpha = 1.
        self.pc_upper_bound = self.eigenvalues.shape[0]
        self.direct_variance_component_init, self.sad_variance_component_init, self.covar_variance_component_init, self.decomp_gwas, self.decomp_sib, self.decomp_diff, self.gwas_avg_se2,self.sib_avg_se2, self.proj_gwas, self.proj_sib, self.proj_diff, self.variance_direct_vc, self.variance_sad_vc, self.variance_covar_vc, self.standard_variance_component, self.gwas_beta_threshed, self.gwas_se_threshed, self.sib_beta_threshed, self.sib_se_threshed, self.beta_sum, self.alpha, self.alpha_se, self.vartotals_init, self.tau = self.estimate_components_and_alpha(self.gwas_beta, self.gwas_se,self.sib_beta,self.sib_se,self.ascertainment_p, self.thresh,self.eigenvecs, self.eigenvalues, 1, 1, self.pc_upper_bound, pc_lower_bound = self.pc_lower_bound)
        self.alpha_report = self.alpha
        print(self.outlabel)
        print(self.thresh)
        print('first:' + str(self.alpha_report) + '(' + str(self.alpha_se) + ')')

        self.gwas_beta = gwas_beta.reset_index(drop = True)/self.alpha_report
        self.gwas_se = gwas_se.reset_index(drop = True)/self.alpha_report

        self.sib_beta = sib_beta.reset_index(drop = True)
        self.sib_se = sib_se.reset_index(drop = True)
        self.pc_upper_bound = self.eigenvalues.shape[0]
        self.ascertainment_p = asc_p.reset_index(drop = True)
        self.thresh = thresh
        self.pc_lower_bound = pc_lower_bound
        
        self.direct_variance_component, self.sad_variance_component, self.covar_variance_component, self.decomp_gwas, self.decomp_sib, self.decomp_diff, self.gwas_avg_se2,self.sib_avg_se2, self.proj_gwas, self.proj_sib, self.proj_diff, self.variance_direct_vc, self.variance_sad_vc, self.variance_covar_vc, self.standard_variance_component, self.gwas_beta_threshed, self.gwas_se_threshed, self.sib_beta_threshed, self.sib_se_threshed, self.beta_sum, self.alpha, self.alpha_se_new, self.vartotals, self.tau = self.estimate_components_and_alpha(self.gwas_beta, self.gwas_se,self.sib_beta,self.sib_se, self.ascertainment_p, self.thresh,self.eigenvecs, self.eigenvalues, 0, self.tau, self.pc_upper_bound, pc_lower_bound = self.pc_lower_bound)
        
        if block_perm:
            block_permutation(chr_pos, self.gwas_beta_threshed, self.gwas_se_threshed, self.sib_beta_threshed, self.sib_se_threshed, self.ascertainment_p, self.thresh, outlabel, self.eigenvecs, self.eigenvalues, self.direct_variance_component, self.sad_variance_component, self.covar_variance_component, self.decomp_gwas, self.decomp_sib, self.decomp_diff,self.variance_direct_vc, self.variance_sad_vc, self.variance_covar_vc, outpath, pos_label, self.pcs_to_test, nperm = nperm)
        else:
            pass

    def element_multiplier(self,x,y):
        return np.multiply(x,y)

    def prs_decomp(self,eff_sizes,eigenvalues,eigenvecs):
        eff_sizes = np.array(eff_sizes).reshape(self.eigenvecs.shape[0])
        temp = np.apply_along_axis(self.element_multiplier,0,eigenvecs,eff_sizes)
        temp = np.sum(temp, axis = 0)
        temp = np.power(temp,2)
        temp = np.multiply(eigenvalues,temp)
        return temp
        
    def error_decomp(self,ses,eigenvalues,eigenvecs):
        ses = np.power(np.array(ses).reshape(self.eigenvecs.shape[0]),2)
        temp = np.apply_along_axis(self.element_multiplier, 0, np.power(eigenvecs,2), ses)
        temp = np.sum(temp, axis = 0)
        temp = np.multiply(eigenvalues,temp)
        return temp

    def deming_error_decomp(self,errors,eigenvalues,eigenvecs):
        errors = np.array(errors).reshape(self.eigenvecs.shape[0])
        temp = np.sum(np.apply_along_axis(self.element_multiplier, 0, np.power(eigenvecs,2), errors),axis = 0)
        temp = np.power(temp,2)
        temp = 2*np.multiply(np.power(eigenvalues,2),temp)
        return temp

    def prs_decomp_unsq(self,eff_sizes,eigenvalues,eigenvecs):
        eff_sizes = np.array(eff_sizes).reshape(self.eigenvecs.shape[0])
        temp = np.apply_along_axis(self.element_multiplier, 0, eigenvecs,eff_sizes)
        temp = np.sum(temp, axis = 0)
        temp = np.multiply(eigenvalues,temp)
        return temp

    def beta_p_thresh(self,betas,p,thresh):
        filt_p = p[p >= thresh].index.tolist()
        self.ascertained_snp_index = p[p < thresh].index.tolist()
        betas.loc[betas.index.isin(filt_p)] = 0
        self.nsnp = betas.shape[0] - len(filt_p)
        return betas

    def se2_avg_p_thresh(self,ses,p,thresh):
        filt_p = p[p >= thresh].index.tolist()
        return np.mean(ses.loc[~ses.index.isin(filt_p)]**2)

    def f(self,B,x):
        return B[0]*x

    def se_bootstrapper(self,lmdf):
        linear = Model(self.f)
        estimates = []
        for strap in range(100):
            indices = np.random.choice(lmdf.index.tolist(), lmdf.shape[0])
            temp_lmdf = lmdf.loc[indices]
            mydata = RealData(x=temp_lmdf['sib_vc'],y=temp_lmdf['standard_vc'], sx = temp_lmdf['var_sib_vc'], sy = temp_lmdf['var_standard_vc'])
            myodr = ODR(mydata, linear, beta0 = [0.])
            myoutput = myodr.run()
            estimates.append(np.sqrt(np.abs(myoutput.beta[0])))
        return np.std(estimates)


    def estimate_components_and_alpha(self,gwas_beta,gwas_se,sib_beta,sib_se,ascertainment_p,thresh,eigenvecs,eigenvalues, plotter, tau,
        pc_upper_bound, pc_lower_bound = 100):
        
        #get all necessary statistics from the population gwas
        #get gwas effects that have p values less than threshold
        gwas_beta_threshed = self.beta_p_thresh(gwas_beta, ascertainment_p, thresh)
        #get gwas SEs that have p values less than threshold
        gwas_se_threshed = self.beta_p_thresh(gwas_se, ascertainment_p, thresh)
        #get gwas SEs that have p values less than threshold
        gwas_avg_se2 = self.se2_avg_p_thresh(gwas_se, ascertainment_p, thresh)
        
        #get all necessary statistics from the sib gwas
        #get sib effects that have p values less than threshold
        sib_beta_threshed = self.beta_p_thresh(sib_beta, ascertainment_p, thresh)
        #get sib effects that have p values less than threshold
        sib_se_threshed = self.beta_p_thresh(sib_se, ascertainment_p, thresh)

        #get sib SEs that have p values less than threshold
        sib_avg_se2 = self.se2_avg_p_thresh(sib_se, ascertainment_p, thresh)
        #perform the decomposition for each component of interest       
        decomp_gwas = self.prs_decomp(gwas_beta_threshed,eigenvalues,eigenvecs)
        decomp_sib = self.prs_decomp(sib_beta_threshed,eigenvalues,eigenvecs)
        decomp_diff = self.prs_decomp(gwas_beta_threshed - sib_beta_threshed, eigenvalues, eigenvecs)

        #perform the decomposition for each component of interest, unsquared
        proj_gwas = self.prs_decomp_unsq(gwas_beta_threshed,eigenvalues,eigenvecs)
        proj_sib = self.prs_decomp_unsq(sib_beta_threshed,eigenvalues,eigenvecs)
        proj_diff = self.prs_decomp_unsq(gwas_beta_threshed - sib_beta_threshed,eigenvalues,eigenvecs)
        
        #decompose the standard errors of the gwas and sibling standard errors
        decomp_gwas_se = self.error_decomp(gwas_se_threshed, eigenvalues, eigenvecs)
        decomp_sib_se = self.error_decomp(sib_se_threshed, eigenvalues, eigenvecs)
        direct_variance_component = decomp_sib - decomp_sib_se
        sad_variance_component = decomp_diff - decomp_sib_se - decomp_gwas_se
        covar_variance_component = decomp_gwas - decomp_diff - decomp_sib + 2*decomp_sib_se
        standard_variance_component = decomp_gwas - decomp_gwas_se
        
        #get variance of each variance component
        variance_direct_vc = 2*(decomp_sib_se**2)
        variance_sad_vc = 2*((decomp_sib_se+decomp_gwas_se)**2)
        variance_covar_vc = 4*(decomp_sib_se+decomp_gwas_se) + 8*(decomp_sib_se**2)
                
        variance_standard_vc = np.power(2*((decomp_gwas_se)**2),0.5)
        variance_sib_vc = np.power(2*((decomp_sib_se)**2),0.5)

        startdf = np.vstack((standard_variance_component,direct_variance_component,variance_standard_vc,variance_sib_vc))
        lmdf = pd.DataFrame(data=startdf, index = ['standard_vc','sib_vc','var_standard_vc','var_sib_vc']).T
        lmdf = lmdf.astype(float)
        lmdf = lmdf.iloc[pc_lower_bound:pc_upper_bound]
        
        linear = Model(self.f)
        mydata = RealData(x=lmdf['sib_vc'],y=lmdf['standard_vc'],sx = lmdf['var_sib_vc'], sy = lmdf['var_standard_vc'])
        myodr = ODR(mydata, linear, beta0 = [0.])
        myoutput = myodr.run()
        alpha = np.sqrt(np.abs(myoutput.beta[0]))
        alpha_se = self.se_bootstrapper(lmdf)

        if self.plot == 0:
            fig, ax = plt.subplots(1,1,figsize = (5,5))
            ax.scatter(lmdf['sib_vc'],lmdf['standard_vc'],color = 'grey',s = 3)
            ax.errorbar(lmdf['sib_vc'],lmdf['standard_vc'],xerr = lmdf['var_sib_vc'],color = 'grey',elinewidth=1)
            ax.errorbar(lmdf['sib_vc'],lmdf['standard_vc'],yerr = lmdf['var_standard_vc'],color = 'grey',elinewidth=1)
            xseq = ax.get_xlim()
            ax.plot(xseq,xseq,color = 'black',linestyle='--',alpha=0.75)
            plt.tight_layout()
            plt.savefig(self.outpath + '/' + self.outlabel +'.pval.' + str(self.thresh) + '.pdf')
            self.plot = 1
            
        total_var = np.sum([np.sum(direct_variance_component),np.sum(sad_variance_component),np.sum(covar_variance_component),np.sum(decomp_gwas_se)])
        var_props_all_pcs = np.array([np.sum(direct_variance_component),np.sum(sad_variance_component),np.sum(covar_variance_component),np.sum(decomp_gwas_se)])
        var_props_all_pcs[:3] = var_props_all_pcs[:3]/float(np.sum(var_props_all_pcs[:3]))
        var_props_bottom_pcs= [np.sum(direct_variance_component[self.pc_lower_bound:]),np.sum(sad_variance_component[self.pc_lower_bound:]),np.sum(covar_variance_component[self.pc_lower_bound:]),np.sum(decomp_gwas_se[self.pc_lower_bound:])]/np.sum(decomp_gwas[self.pc_lower_bound:])
    
        colnames = ['Direct', 'SAD', 'Direct-SAD covariance', 'Error']
        rownames = ['All_PCs', 'PC ' + str(self.pc_lower_bound) + '+']

        vartotals = pd.DataFrame(np.vstack((var_props_all_pcs,var_props_bottom_pcs)), index = rownames, columns = colnames).T

        return direct_variance_component, sad_variance_component, covar_variance_component, decomp_gwas, decomp_sib, decomp_diff, gwas_avg_se2, sib_avg_se2, proj_gwas, proj_sib, proj_diff, variance_direct_vc, variance_sad_vc, variance_covar_vc, standard_variance_component, gwas_beta_threshed, gwas_se_threshed, sib_beta_threshed, sib_se_threshed, np.sum(gwas_beta_threshed[gwas_beta_threshed != 0]), alpha, alpha_se, vartotals, tau

    def pca(self, genotype_mat):

        pca = PCA()
        pca.fit_transform(genotype_mat)
        
        self.eigenvecs = pca.components_.T
        self.eigenvalues = pca.explained_variance_.T

        if self.outlabel != None:
            np.save(self.outpath + '/' + self.outlabel + '.eigenvectors', self.eigenvecs)
            np.save(self.outpath + '/' + self.outlabel + '.eigenvalues', self.eigenvalues)
        else:
            np.save(self.outpath + '/eigenvectors', self.eigenvecs)
            np.save(self.outpath + '/eigenvalues', self.eigenvalues)

    def outputs(self):
        return {'direct_vc':self.direct_variance_component,'sad_vc':self.sad_variance_component, 'covar_vc':self.covar_variance_component, 
            'decomp_gwas':self.decomp_gwas, 'decomp_sib':self.decomp_sib, 'decomp_diff':self.decomp_diff,
            'gwas_avg_se2':self.gwas_avg_se2,'sib_avg_se2':self.sib_avg_se2,
            'proj_gwas':self.proj_gwas, 'proj_sib':self.proj_sib, 'proj_diff':self.proj_diff,
            'var_direct_vc':self.variance_direct_vc, 'var_sad_vc':self.variance_sad_vc, 'var_covar_vc':self.variance_covar_vc,
            'beta_sum':self.beta_sum, 'alpha':self.alpha_report,'alpha_se':self.alpha_se_new, 'nsnp':self.nsnp, 
            'var_totals':self.vartotals_init}


        



