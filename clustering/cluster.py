from .utils import Atom, Residue, ActiveSite
from .io import load_obj, save_obj, read_active_sites
import numpy as np
from itertools import combinations as combo
from copy import deepcopy

def compute_similarity(site_a, site_b, redo_mean_dev=False):
    """
    Compute the similarity between two given ActiveSite instances.

    Input: two ActiveSite instances
        -NOTE, optional "redo_mean_dev". Default: FALSE = normalize according to 
            the factors derived from the original dataset. Otherwise (e.g. if True), 
            normalize according to data currently in data folder
            
    Output: 
        -the similarity between them (a floating point number)
        -the two vectors used to calculate similarity (in a list)
    
    Logic behind similarity computation:
        - Average of several metrics with different biological meaning
        - Attempt to cover chemistry, morphology and homogeneity with metrics
        
    Metrics:
        Within Active Site (calculated in individual_metrics):
            - Mean and Variance of total charge/polarity of residues
            - Fraction of carbon / nitrogen / oxygen
            - Fraction of residue TYPES (by polarity/charge)     
            - Mean distance between residue centers of mass
            - Total number of atoms / residues
            - Fraction of residue TYPES (by polarity/charge)     
    
    Then:
        -Normalize the metrics vectors
        -Find the Euclidean distance between the vectors
        -Transform distance to similarity on range [0,1]
    """

    #get all of the summary metrics for each site
    site_metrics = [individual_metrics(site_a),individual_metrics(site_b)]
    
    #normalize metric vectors
    
        #if you want to recalculate mean/devs
    if redo_mean_dev:
        gen_mean_dev_normalizations()
    
    means = load_obj('means')
    devs = load_obj('devs')
    
    #normalize
    for i in range(len(site_metrics)):
        for metric, value in site_metrics[i].items():
            site_metrics[i][metric] = (value-means[metric])/float(devs[metric])
    
    #convert to arrays
    #NOTE: When converting into array, follows the metric_names list
    #so that you can consistently compare numpy arrays element-wise
    
    metric_names = ['meanChar','varChar','meanPol','varPol',
        'elemFracC','elemFracN','elemFracO','numAtoms','numRes','meanResDist','varResDist', 
        'nonpolarFrac','polarFrac','acidicFrac','basicFrac']
        
    metric_arr_a = np.ndarray(0)
    metric_arr_b = np.ndarray(0)
    
    for metric in metric_names:
        metric_arr_a = np.append(metric_arr_a,site_metrics[0][metric])
        metric_arr_b = np.append(metric_arr_b,site_metrics[1][metric])      
    
    metric_arrays = [metric_arr_a, metric_arr_b]
    
    #take the distance between them, then normalize to [0,1] range
    distance = np.linalg.norm(metric_arr_a - metric_arr_b)
    similarity = 1.0/(1.0+distance)
    
    return similarity, metric_arrays

def individual_metrics(active_site,printMe=None):
    """
    Compute metrics of chemistry, morphology and homogeneity for an ActiveSite instance.

    Input: ActiveSite instance
    Output: a dictionary of metrics
           
    Metrics:
        - Mean and Variance of total charge/polarity of residues
        - Fraction of carbon / nitrogen / oxygen
        - Total number of atoms / residues
        - Mean distance between residue centers of mass
    """

    chem_properties = {
        "GLY": "nonpolar",
        "ALA": "nonpolar",
        "VAL": "nonpolar",
        "ILE": "nonpolar",
        "LEU": "nonpolar",
        "MET": "nonpolar",
        "PHE": "nonpolar",
        "TYR": "nonpolar",
        "TRP": "nonpolar",
        "PRO": "nonpolar",
        "CYS": "nonpolar",
        "SER": "polar",
        "THR": "polar",
        "ASN": "polar",
        "GLN": "polar",
        "ARG": "basic",
        "HIS": "basic",
        "LYS": "basic",
        "ASP": "acidic",
        "GLU": "acidic",
        }
        
    charge = {"nonpolar": 0, "polar": 0, "basic": 1, "acidic": -1}
    polarity = {"nonpolar":0, "polar": 0.5, "basic": 1, "acidic": 1}
    res_type = {"nonpolar":0, "polar": 0, "basic": 0, "acidic": 0}
    
    atom_type = {"C":0,"N":0,"O":0}
    
    metric_names = ['meanChar','varChar','meanPol','varPol',
        'elemFracC','elemFracN','elemFracO','numAtoms','numRes','meanResDist','varResDist', 
        'nonpolarFrac','polarFrac','acidicFrac','basicFrac']
    
    metrics = dict.fromkeys(metric_names)
    
    #first calculate mean and variance of total charges/polarity of residues
    residue_charge = np.ndarray(0)
    residue_polarity = np.ndarray(0)
    numRes = len(active_site.residues)
    
    for residue in active_site.residues:
        res_type[chem_properties[residue.type.strip()]] += 1
        residue_charge = np.append(residue_charge, charge[chem_properties[residue.type.strip()]])
        residue_polarity = np.append(residue_polarity, polarity[chem_properties[residue.type.strip()]])
    
    for key in res_type.keys():
        res_type[key] /= numRes
        
    update_dict(metrics,
        ('meanChar','varChar','meanPol','varPol','nonpolarFrac','polarFrac','acidicFrac','basicFrac'),
        (np.mean(residue_charge),np.var(residue_charge),
            np.mean(residue_polarity),np.var(residue_polarity),
            res_type['nonpolar'],res_type['polar'],res_type['acidic'],res_type['basic'])
        )
      
    #now calculate the fraction of atoms of carbon/nitrogen/oxygen in the active zone
    numAtoms = 0
        
    for residue in active_site.residues:
        for atom in residue.atoms:
            numAtoms += 1
            type = atom.type[0]
            if type.lower()=='c':
                atom_type['C'] += 1
            elif type.lower()=='n':
                atom_type['N'] += 1
            elif type.lower()=='o':
                atom_type['O'] += 1
    
    for key in atom_type.keys():
        atom_type[key] /= numAtoms
    
    update_dict(metrics,
        ('elemFracC','elemFracN','elemFracO','numAtoms','numRes'),
        (atom_type['C'],atom_type['N'],atom_type['O'],numAtoms,numRes)
        )
    
    #now calculate the "center of mass" for each residue.
    #then calculate the mean/var of the pairwise Euclidean distance between residue COMs.  

    for residue in active_site.residues:
        numAtoms = 0
        centerOfMass = np.zeros(3)
        
        for atom in residue.atoms:
            numAtoms += 1
            centerOfMass += np.asarray(atom.coords)
        
        if not numAtoms==0:
            centerOfMass /= numAtoms
        
        residue.com = tuple(centerOfMass)
    
    residue_dist = []
    
    for pair in combo(active_site.residues,2):
        distance = np.linalg.norm(np.asarray(pair[0].com)-np.asarray(pair[1].com))
        residue_dist.append(distance)  
       
    update_dict(metrics,
        ('meanResDist','varResDist'),
        (np.mean(residue_dist),np.var(residue_dist))
        )
    
    if printMe:    
        for metric in metrics.keys():
            print('{0}\n\t\t{1}\n'.format(metric,metrics[metric]))
        
    return metrics

def gen_mean_dev_normalizations():
    """
    Used to generate Mean and Deviation dictionaries for similarity function.
        -generates values from all PDB files in data folder
        -stores them in pickled format
    cluster.compute_similarity uses 
    """
    activeSites = read_active_sites('data')
    
    #initialize
    means = deepcopy(individual_metrics(activeSites[0]))
    devs = deepcopy(individual_metrics(activeSites[0]))
    
    #sum (with scale for auto-mean)
    for i in range(len(activeSites)):
            
        site = activeSites[i]
        
        siteMetrics = individual_metrics(site)
        
        for metric, value in siteMetrics.items():
            if i==0:
                means[metric] /= float(len(activeSites))
            else:
                means[metric] += siteMetrics[metric]/float(len(activeSites))
                
    
    #now get mean absolute deviation
    for i in range(len(activeSites)):
            
        site = activeSites[i]
        
        siteMetrics = individual_metrics(site)
        
        for metric, value in siteMetrics.items():        
            if i==0:
                devs[metric] = 0
            devs[metric] += abs(siteMetrics[metric]-means[metric])/float(len(activeSites))
                
                
                    
    save_obj(means,'means')
    save_obj(devs,'devs')
        
    return means, devs

def update_dict(dictionary,keys,values):
    """
    Update dictionary given lists of keys and values
    
    Input: dictionary, keys, values
    Output: None (dictionary updated automatically)
    """

    dictionary.update(dict(zip(keys,values)))
    
def find_active_sites(activeSites, targetName):
    """
    Find an active site in a list of AS's given the name
    
    Input: list of active sites, name of target
    Output: list of indeces of matched active sites
    """
    
    return [i for i, j in enumerate(activeSites) if j.name==targetName]

def cluster_by_partitioning(active_sites,num_clusters=3):
    """
    Cluster a given set of ActiveSite instances using a partitioning method.

    Input: a list of ActiveSite instances, (OPTIONAL: number of clusters, default 3)
    Output: a clustering of ActiveSite instances
            (this is really a list of clusters, each of which is list of
            ActiveSite instances)
    """
    
    # Fill in your code here!

    return []


def cluster_hierarchically(active_sites):
    """
    Cluster the given set of ActiveSite instances using a hierarchical algorithm.

    Input: a list of ActiveSite instances
    Output: a list of clusterings
            (each clustering is a list of lists of Sequence objects)
    """

    # Fill in your code here!

    return []
