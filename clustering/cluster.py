from .utils import Atom, Residue, ActiveSite
from .io import read_active_sites
import numpy as np
from random import seed

def cluster_by_partitioning(active_sites,num_clusters=3, max_iters=10000, dist_thresh=0.01):
    """
    Cluster a given set of ActiveSite instances using a partitioning method.

    Input: a list of ActiveSite instances, (OPTIONAL: number of clusters, default 5)
    Output: a clustering of ActiveSite instances
            (this is really a list of clusters, each of which is list of
            ActiveSite instances)
    """
    
    cluster_centers = k_means(active_sites,num_clusters,max_iters,dist_thresh)
    
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
    if num_clusters>len(active_sites) or num_clusters<1:
        print("Invalid number of clusters: Default to 5")
        num_clusters = 5
    
    # initialize centroids randomly by choosing X random points
    seed()
    inds = np.random.choice(len(activeSites),num_clusters,replace=False)
    cluster_centers = np.array([active_sites[i].get_norm_metrics() for i in inds])
    
    # init trackers
    iter = 0
    prev_centers = np.zeros(cluster_centers.shape)
    
    # begin algorithm
    while not iter > max_iters and not distance_cutoff(prev_centers,cluster_centers,dist_thresh):
        prev_centers = cluster_centers
        iter += 1
        if iter%1000==0:
            print("Iteration no.", iter)
        
        #assign all objects to a cluster, then reevaluate cluster centers
        labels = assign_to_clusters(active_sites,cluster_centers)
        cluster_centers = update_cluster_locs(active_sites, labels, num_clusters)
    
        
    return cluster_centers

def distance_cutoff(prev,current,thresh):
    """
    Input: numpy arrays of previous/current centroid locations, and threshold
    Output: boolean of whether centroids have moved farther than the threshold on average
    """
    
    sum = 0
    
    for i in range(prev.shape[0]):
        sum += np.linalg.norm(prev[i,:] - current[i,:])
    
    sum /= prev.shape[0]
    
    if sum < thresh:
        print("Threshold reached, mean center distance travelled is ", sum)
    
    return sum < thresh
    
        
def assign_to_clusters(active_sites,cluster_centers):
    """
    Input: List of ActiveSite instances, and 2D numpy array of cluster centers
    Output: List of labels, value is the row (cluster) that Site was assigned to.
    """
    
    labels = np.zeros([len(active_sites)])  
    
    for i in range(len(active_sites)):
        site_metrics = active_sites[i].get_norm_metrics()
        cluster_matches = np.sqrt(((cluster_centers - site_metrics)**2).sum(axis=1))
        labels[i] = np.argmax(cluster_matches)
        
    return labels

def update_cluster_locs(active_sites, labels, num_clusters):
    """
    Input: List of ActiveSite instances and labels, and the number of clusters
    Output: Numpy array of new cluster centers
    """
    
    new_cluster_centers = np.zeros(shape=(num_clusters,len(active_sites[0].get_norm_metrics())))
    
    # for each cluster, get all assigned sites and average their metrics
    for clust in range(num_clusters):
        site_inds = [ind for ind,val in enumerate(labels) if val==clust]
        all_clust_metrics = np.array([active_sites[i].get_norm_metrics() for i in site_inds])
        new_center = all_clust_metrics.mean(axis=0)
        new_cluster_centers[clust,:] = new_center

def cluster_hierarchically(active_sites):
    """
    Cluster the given set of ActiveSite instances using a hierarchical algorithm.

    Input: a list of ActiveSite instances
    Output: a list of clusterings
            (each clustering is a list of lists of Sequence objects)
    """

    # Fill in your code here!

    return []
