import openpyxl
import pandas as pd
import brightway2 as bw
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
import re
import matplotlib.colors

def _plot_params(param_list,
                 logx=False,
                 logy=False,
                 figsize=(8,4),
                 columns=2,
                 size=1000,
                 **plot_params):
    '''
    Plot the distribution of a parameter or more
    from Thomas Gibon "10.1021:acs.est.3c03190"
    three distribution allowed: linear (uniform), triangle, lognormal (not available: beta  )
    '''
    
    if len(param_list)==1:
        fig, ax = plt.subplots(len(param_list),
                               figsize=figsize)
    else:
        fig, axes = plt.subplots(int(np.ceil(len(param_list) / columns)),
                                 columns,
                                 figsize=figsize)
    
    for i,param in enumerate(param_list):
        
        if len(param_list)>1:
            if len(param_list) == columns:
                ax=axes[i]
            else:
                ax=axes[i//columns, i%columns]
                
        if param.distrib == 'linear':
            
            sns.histplot(np.random.uniform(param.min,
                                  param.max,size=size),
                         ax=ax,kde=True,stat='density',kde_kws=dict(cut=3),
                        **plot_params)
            
        if param.distrib == 'triangle':
            
            sns.histplot(np.random.triangular(param.min,
                                              param.default,
                                              param.max,size=size),
                         ax=ax,kde=True,stat='density',kde_kws=dict(cut=3),
                        **plot_params)
            
        if param.distrib == 'lognormal':
            # Check https://en.wikipedia.org/wiki/Log-normal_distribution#Estimation_of_parameters

            sns.histplot(np.random.lognormal(param.mean,
                                             param.std,
                                             size=size),
                         ax=ax,kde=True,stat='density',kde_kws=dict(cut=3),
#                          log_scale=True,
                        **plot_params)
            
        if param.unit:
            ax.set_xlabel(f'{param.label} [{param.unit}]')
        
        ax.set_title(param.name)
    
    #for log transformation        
    if logx:  
        plt.semilogx()
    if logy:
        plt.semilogy()
    plt.tight_layout()
    plt.show()
    
