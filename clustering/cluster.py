from .utils import Atom, Residue, ActiveSite
import numpy as np
from itertools import combinations as combo

def update_dict(dictionary,keys,values):
    dictionary.update(dict(zip(keys,values)))

def compute_similarity(site_a, site_b, norm_factors):
    """
    Compute the similarity between two given ActiveSite instances.

    Input: two ActiveSite instances
        -NOTE, optional "custom_norm". Default: normalize according to 
            the factors derived from the original dataset. Otherwise
            (e.g. if True), normalize according to custom dictionary.
            (for use within large datasets only).
            
    Output: the similarity between them (a floating point number)
    
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
        -Find the cosine of the angle between the vectors
    """

    #get all of the summary metrics for each site
    metrics_site_a = individual_metrics(site_a)
    metrics_site_b = individual_metrics(site_b)
    
    #normalize metric vectors
    
    for metric in metrics_site_a.keys():
        
    
    #take the cosine of the angle between them
            
        
    
    similarity = 0.0

    return similarity

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
    
    metric_names = ['meanChar','varChar','meanPol','varPol','elemFrac',
        'numAtoms','numRes','meanResDist','varResDist', 'resTypes']
    
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
        ('meanChar','varChar','meanPol','varPol','resTypes'),
        (np.mean(residue_charge),np.var(residue_charge),
            np.mean(residue_polarity),np.var(residue_polarity),
            res_type)
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
        ('elemFrac','numAtoms','numRes'),
        (atom_type,numAtoms,numRes)
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
    
#    print(residue_dist)  
    
    update_dict(metrics,
        ('meanResDist','varResDist'),
        (np.mean(residue_dist),np.var(residue_dist))
        )
    
    if printMe:    
        for metric in metrics.keys():
            print('{0}\n\t\t{1}\n'.format(metric,metrics[metric]))
        
    return metrics


def cluster_by_partitioning(active_sites):
    """
    Cluster a given set of ActiveSite instances using a partitioning method.

    Input: a list of ActiveSite instances
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
