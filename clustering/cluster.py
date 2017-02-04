from .utils import Atom, Residue, ActiveSite
import numpy as np
from itertools import combinations as combo

def update_dict(dictionary,keys,values):
    dictionary.update(dict(zip(keys,values)))

def compute_similarity(site_a, site_b):
    """
    Compute the similarity between two given ActiveSite instances.

    Input: two ActiveSite instances
    Output: the similarity between them (a floating point number)
    
    Logic behind similarity computation:
        - Average of several metrics with different biological meaning
        - Attempt to cover chemistry, morphology and homogeneity with metrics
        
    Metrics:
        Within Active Site (calculated in individual_metrics):
            - Mean and Variance of total charge/polarity of residues
            - Fraction of carbon / nitrogen / oxygen
            - Mean distance between residue centers of mass
            - Total number of atoms / residues
        Between two sites:
            - Fraction of residues shared
            - Fraction of residue TYPES (by polarity/charge) shared            
    """


    #get all of the summary metrics for each site
    metrics_site_a = individual_metrics(site_a)
    metrics_site_b = individual_metrics(site_b)
     
    similarity = 0.0

    return similarity

def individual_metrics(active_site):
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
    
    metric_names = ['meanChar','varChar','meanPol','varPol','elemFrac',
        'numAtoms','numRes','meanResDist','varResDist']
    
    metrics = dict.fromkeys(metric_names)
    
    
    #first calculate mean and variance of total charges/polarity of residues
    residue_charge = np.ndarray(0)
    residue_polarity = np.ndarray(0)
    
    for residue in active_site.residues:
        residue_charge = np.append(residue_charge, charge[chem_properties[residue.type.strip()]])
        residue_polarity = np.append(residue_polarity, polarity[chem_properties[residue.type.strip()]])
        
    update_dict(metrics,
        ('meanChar','varChar','meanPol','varPol'),
        (np.mean(residue_charge),np.var(residue_charge),
            np.mean(residue_polarity),np.var(residue_polarity))
        )
      
    #now calculate the fraction of atoms of carbon/nitrogen/oxygen in the active zone
    fractions = np.zeros(3)
    numAtoms = 0
    numRes = len(active_site.residues)
        
    for residue in active_site.residues:
        for atom in residue.atoms:
            numAtoms += 1
            type = atom.type[0]
            if type.lower()=='c':
                fractions[0] += 1
            elif type.lower()=='n':
                fractions[1] += 1
            elif type.lower()=='o':
                fractions[2] += 1
    
    fractions /= numAtoms
    
    update_dict(metrics,
        ('elemFrac','numAtoms','numRes'),
        (fractions,numAtoms,numRes)
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
    
    print(residue_dist)  
    
    update_dict(metrics,
        ('meanResDist','varResDist'),
        (np.mean(residue_dist),np.var(residue_dist))
        )
        
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
