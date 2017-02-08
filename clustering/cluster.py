from .utils import Atom, Residue, ActiveSite
from .io import read_active_sites
import numpy as np
from random import seed

def cluster_by_partitioning(active_sites,num_clusters=5, max_iters=10000, dist_thresh=0.01):
    """
    Cluster a given set of ActiveSite instances using a partitioning method.

    Input: a list of ActiveSite instances, (OPTIONAL: number of clusters, default 5)
    Output: a clustering of ActiveSite instances
            (this is really a list of clusters, each of which is list of
            ActiveSite instances)
    """
    
    cluster_centers = k_means(active_sites,num_clusters,max_iters)
    
    return []
    
def k_means(active_sites,num_clusters,max_iters,dist_thresh):
    """
    K-Means clustering
    
    Input: List of active site instances, number of clusters and number of iterations
    Output: List of cluster centers as numpy arrays
    """
    
    # edge cases
    num_clusters = round(num_clusters)
    max_iters = round(max_iters)
    if num_clusters>len(active_sites) or numclusters<1:
        print("Invalid number of clusters: Default to 5")
        num_clusters = 5
    
    # initialize centroids randomly by choosing X random points
    seed()
    inds = np.random.choice(len(activeSites),num_clusters,replace=False)
    cluster_centers = np.array([active_sites[i].get_norm_metrics() for i in inds])
    
    # init trackers
    iter = 0
    labels = np.zeros([len(active_sites)])    
    prev_centers = np.zeros(cluster_centers.shape)
    
    # begin algorithm
    while not iter > max_iters and not distance_cutoff(prev_centers,cluster_centers,dist_thresh):
        prev_centers = cluster_centers
        iter += 1
        
        #assign all objects to a cluster, then reevaluate cluster centers
        labels = assign_to_clusters(active_sites,cluster_centers)
        cluster_centers = evaluate_cluster_locs(active_sites, labels, num_clusters)
        
    return cluster_centers

def distance_cutoff(prev,current,thresh):
    """
    Input: numpy arrays of previous/current centroid locations, and threshold
    Output: boolean of whether centroids have moved farther than the threshold
    """
    
    
    
        
def assign_to_clusters():

def evaluate_cluster_locs():


def cluster_hierarchically(active_sites):
    """
    Cluster the given set of ActiveSite instances using a hierarchical algorithm.

    Input: a list of ActiveSite instances
    Output: a list of clusterings
            (each clustering is a list of lists of Sequence objects)
    """

    # Fill in your code here!

    return []
