from .utils import Atom, Residue, ActiveSite
import numpy as np

def compute_similarity(site_a, site_b, redo_mean_dev=False):
    """
    Compute the similarity between two given ActiveSite instances.

    Input: two ActiveSite instances
        -NOTE, optional "redo_mean_dev". Default: FALSE = normalize according to 
            the factors derived from the original dataset. Otherwise (e.g. if True), 
            normalize according to data currently in data folder
            
    Output: 
        -the similarity between them (a floating point number)
    
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
        -get vectors normalized to whole population (or saved file)
        -Find the Euclidean distance between the vectors
        -Transform distance to similarity on range [0,1]
    """
    
    #fetch normalized metrics
    metrics = [site_a.get_norm_metrics(redo_mean_dev), site_b.get_norm_metrics(redo_mean_dev)]
    
    #take the Euclidean distance between them, then normalize to [0,1] range
    distance = np.linalg.norm(metrics[0] - metrics[1])
    similarity = 1.0/(1.0+distance)
    
    return similarity


def find_active_sites(activeSites, targetName):
    """
    Find an active site in a list of AS's given the name
    
    Input: list of active sites, name of target
    Output: list of indeces of matched active sites
    """
    
    return [i for i, j in enumerate(activeSites) if j.name==targetName]
