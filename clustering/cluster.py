from .utils import Atom, Residue, ActiveSite
from .io import read_active_sites
import numpy as np
from random import seed
from itertools import combinations as combo
from math import log

###########################
# Clustering By Partition #
###########################


def cluster_by_partitioning(active_sites,num_clusters=7, max_iters=10000, dist_thresh=0.001):
    """
    Cluster a given set of ActiveSite instances using a partitioning method.

    Input: a list of ActiveSite instances
    (OPTIONAL): number of clusters, maximum iterations, and distance threshold
    Output: a clustering of ActiveSite instances
            (this is really a list of clusters, each of which is list of
            ActiveSite instances)
    """
    
    labels = k_means(active_sites,num_clusters,max_iters,dist_thresh)
    
    clustering = []
    
    for clust in range(num_clusters):
        clustering.append([active_sites[i] for i in range(len(labels)) if labels[i]==clust])
    
    return clustering
    
def k_means(active_sites,num_clusters,max_iters,dist_thresh,printMe=False):
    """
    K-Means clustering
    
    Input: List of active site instances, number of clusters, number of iters and tresh
    Output: List of labels and cluster centers as numpy arrays
    """
    seed()

    # edge cases
    num_clusters = round(num_clusters)
    max_iters = round(max_iters)
    if printMe:
        print('------')
    if num_clusters>len(active_sites) or num_clusters<1:
        if printMe:
            print("Invalid number of clusters: Default to 7")
        num_clusters = 7
    else:
        if printMe:
            print('Number of Clusters:',num_clusters)
    
    # initialize centroids randomly by choosing X random points
    iter = 0
    inds = np.random.choice(len(active_sites),num_clusters,replace=False)
    cluster_centers = np.array([active_sites[i].get_norm_metrics() for i in inds])
    
    # init trackers (10 most recent ones)
    prev_centers = np.zeros(cluster_centers.shape)

    # begin algorithm
    if printMe:
        print('Maximum Iterations:',max_iters)
        print('Distance Threshold:',dist_thresh)
        print('------')

    while not iter > max_iters and not distance_cutoff(
        prev_centers,cluster_centers,dist_thresh,printMe
        ):
                            
        prev_centers = cluster_centers
        iter += 1
        if iter%1000==0:
            print("Iteration no.", iter)
        
        #assign all objects to a cluster, then reevaluate cluster centers
        labels = assign_to_clusters(active_sites,cluster_centers)
        cluster_centers = update_cluster_locs(active_sites, labels, num_clusters)
    
    if printMe:
        print('Total iterations:',iter,'\n')    
    
    return labels

def distance_cutoff(prev,current,thresh,printMe):
    """
    Input: numpy arrays of previous/current centroid locations, and threshold
    Output: boolean of whether centroids have moved farther than the threshold on average
    """
    
    sum = 0
    
    for i in range(prev.shape[0]):
        sum += np.linalg.norm(prev[i,:] - current[i,:])
    
    sum /= prev.shape[0]
    
    if sum < thresh and printMe:
        print("Threshold reached, mean center distance travelled this step is ", sum)
    
    return sum < thresh    
        
def assign_to_clusters(active_sites,cluster_centers):
    """
    Input: List of ActiveSite instances, and 2D numpy array of cluster centers
    Output: array of labels, value is the row (cluster) that Site was assigned to.
    """
    
    labels = np.zeros([len(active_sites)])  
    
    #for each active site
    for i in range(len(active_sites)):
        site_metrics = active_sites[i].get_norm_metrics()   #pull its metrics
        cluster_dist = np.sqrt(((cluster_centers - site_metrics)**2).sum(axis=1)) #get distances
        labels[i] = np.argmin(cluster_dist)         #label according to closest one
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
    
    return new_cluster_centers

###########################
# Hierarchical Clustering #
###########################

def cluster_hierarchically(active_sites,num_clusters=7):
    """
    Cluster the given set of ActiveSite instances using a hierarchical algorithm.

    Input: a list of ActiveSite instances
    (OPTIONAL): number of clusters (default 7)
    Output: a list of clusterings
            (each clustering is a list of lists of ActiveSite instances)
    """

    labels = centroid_linkage(active_sites,num_clusters)
    
    clustering = []
    
    for clust in np.unique(labels):
        clustering.append([active_sites[ind] for ind,val in enumerate(labels.tolist() )if val==clust])
    
    return clustering

def centroid_linkage(active_sites,min_clusters,printMe=False):
    """
    Centroid Linkage clustering
    
    This implementation builds until it hits min_clusters. Dendrograms not supported.
    
    Inputs: List of ActiveSite instances and number of clusters 
    Outputs: List of labels and cluster centers
    """
    
    seed()
        
    # edge cases
    min_clusters = round(min_clusters)
    if printMe:
        print('------')
    if min_clusters>len(active_sites) or min_clusters<1:
        if printMe:
            print("Invalid number of clusters: Default to 7")
        min_clusters = 7
    else:
        if printMe:
            print('Target Number of Clusters:',min_clusters)
         
    # initialize variables
    num_clusters = len(active_sites)
    labels = np.arange(len(active_sites))
    np.random.shuffle(labels)
       
    # begin algorithm
    while num_clusters > min_clusters:
        
        #calculate the centroid of each cluster and find the smallest distance
        clusters = shortest_centroid_dist(active_sites,labels)
        #merge the two clusters
        labels = merge(clusters,labels)
        num_clusters -= 1
        
        if num_clusters%10==0:
            print('Number of Clusters:',num_clusters)
        
    return labels

def shortest_centroid_dist(active_sites,labels):
    """
    Input: List of active_sites and labels
    Output: List of two cluster labels that are closest
    """
    
    # find unique values in labels and their indices
    
    u, u_ind = np.unique(labels,return_inverse=True)

    # go through all unique pairwise combinations of labels
    # and test the distances between centroids
    
    shortest = 1000000000
    closest_clusters = [0,0]
    
    for i,j in combo(u,2):
        # get metrics from all sites for each cluster pair
        inds_i = [ind for ind,val in enumerate(labels.tolist()) if val==i]
        metrics_i = np.array([active_sites[k].get_norm_metrics() for k in inds_i])   
         
        inds_j = [ind for ind,val in enumerate(labels.tolist()) if val==j]
        metrics_j = np.array([active_sites[m].get_norm_metrics() for m in inds_j])    
        
        # measure distance between centroids and update if necessary
        dist = np.linalg.norm(metrics_j.mean(axis=0)-metrics_i.mean(axis=0))

        if dist < shortest:
            shortest = dist
            closest_clusters = [i,j]
                        
    return closest_clusters

    
def merge(clusters, labels):
    """
    Input: List of two cluster labels and list of all labels
    Output: List of labels, with specified ones merged
    """
        
    #arbitrarily merge first INTO second (so that first label no longer exists)
    labels[labels==clusters[0]] = clusters[1]
        
    return labels
    
###########################
# Eval Clustering Quality #
###########################

def quality(clustering):
    """
    Measures quality of clustering quantitatively
    
    Input: An embedded list of clustered ActiveSite instances (generated by clustering alg)
    Output: The value Q, the average inter-object distance among clusters in the data
    """
    
    Q = 0
    
    # loop through clusters
    for n in range(len(clustering)):
        cluster_mean = 0
        num_pairs = 0
        # loop through pairs per cluster
        for i,j in combo(clustering[n],2):
            num_pairs += 1
            cluster_mean += np.linalg.norm(i.get_norm_metrics()-j.get_norm_metrics())
        # normalize to remain unbiased by cluster size
        if num_pairs==0:
            continue
        else:
            cluster_mean /= num_pairs
            Q += cluster_mean
    
    Q /= len(clustering)
    
    return Q
            
def compare(clust_a,clust_b):
    """
    Compares two culsterings. Positive = clustering B is better than A
    
    Input: Two clusterings from above clustering algorithms
    Output: log-ratio comparison of the quality of the two
    """
    
    return log(quality(clust_a)/quality(clust_b))


