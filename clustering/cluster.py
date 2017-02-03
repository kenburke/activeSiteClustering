from .utils import Atom, Residue, ActiveSite
import numpy as np

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
        Within Active Site:
            - Mean and Variance of total charge/polarity of residues
            - Fraction of carbon / nitrogen / oxygen
            - Mean distance between residue centers of mass
            - Total number of atoms / residues
        Between two sites:
            - Fraction of residues shared
            - Fraction of residue TYPES (by polarity/charge) shared            
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
    
    metrics = ['meanChar','varChar','meanPol','varPol','elemFrac',
        'distCOM','numAtoms','numRes']

    metrics_site_a = dict.fromkeys(metrics)
    metrics_site_b = dict.fromkeys(metrics)
    
    #first calculate mean and variance of total charges/polarity of residues
    charges_a = np.ndarray(0)
    polarity_a = np.ndarray(0)
    charges_b = np.ndarray(0)
    polarity_b = np.ndarray(0)
    
    for residue in site_a.residues:
        charges_a = np.append(charges_a, charge[chem_properties[residue.type.strip()]])
        polarity_a = np.append(polarity_a, polarity[chem_properties[residue.type.strip()]])
        
    update_dict(metrics_site_a,
        ('meanChar','varChar','meanPol','varPol'),
        (np.mean(charges_a),np.var(charges_a),np.mean(polarity_a),np.var(polarity_a)))
    
    for residue in site_b.residues:
        charges_b = np.append(charges_b, charge[chem_properties[residue.type.strip()]])
        polarity_b = np.append(polarity_b, polarity[chem_properties[residue.type.strip()]])

    update_dict(metrics_site_b,
        ('meanChar','varChar','meanPol','varPol'),
        (np.mean(charges_b),np.var(charges_b),np.mean(polarity_b),np.var(polarity_b)))

    
    #now calculate the fraction of atoms of carbon/nitrogen/oxygen in the active zone
    
    
    similarity = 0.0

    return similarity


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
