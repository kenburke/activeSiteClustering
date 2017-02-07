from .utils import Atom, Residue, ActiveSite
from .io import read_active_sites, write_clustering, write_mult_clusterings
from .cluster import individual_metrics, cluster_by_partitioning, cluster_hierarchically
import numpy as np
from itertools import combinations as combo

def runthrough():

    activeSites = read_active_sites('data')
    
    #initialize
    means = dict.fromkeys(individual_metrics(activeSites[0]))
    
    for site in activeSites:
        metrics = individual_metrics(site)
        
        for metric in means:
        
            if isinstance(site[metric],int) or isinstance(site[metric],float):
                means[metric] += site[metric]
            elif isinstance(site[metric],dict):
                
        