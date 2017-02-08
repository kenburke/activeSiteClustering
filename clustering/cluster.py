from .utils import Atom, Residue, ActiveSite
from .io import read_active_sites
import numpy as np
from random import seed, sample

def cluster_by_partitioning(active_sites,num_clusters=5, num_iters=10000):
    """
    Cluster a given set of ActiveSite instances using a partitioning method.

    Input: a list of ActiveSite instances, (OPTIONAL: number of clusters, default 5)
    Output: a clustering of ActiveSite instances
            (this is really a list of clusters, each of which is list of
            ActiveSite instances)
    """
    
    cluster_centers = k_means(active_sites,num_clusters,num_iters)
    
    return []
    
def k_means(active_sites,num_clusters,num_iters):
    """
    K-Means clustering
    
    Input: List of active site instances, number of clusters and number of iterations
    Output: List of cluster centers as numpy arrays
    """
    
    # house-keeping
    num_clusters = round(num_clusters)
    num_iters = round(num_iters)
    if num_clusters>len(active_sites) or numclusters<1:
        print("Invalid number of clusters: Default to 5")
        num_clusters = 5
    
    # initialize centroids randomly by choosing X random points
    seed()
    rand1 = round(random()*len(active_sites))
    rand2 = round(random()*len(active_sites))
    
    
    
    cluster_centers = 0
    
    return cluster_centers
    

def cluster_hierarchically(active_sites):
    """
    Cluster the given set of ActiveSite instances using a hierarchical algorithm.

    Input: a list of ActiveSite instances
    Output: a list of clusterings
            (each clustering is a list of lists of Sequence objects)
    """

    # Fill in your code here!

    return []
