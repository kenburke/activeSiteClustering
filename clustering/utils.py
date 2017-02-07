# Some utility classes to represent a PDB structure
import io
import numpy as np
from itertools import combinations as combo

class Atom:
    """
    A simple class for an amino acid residue
    """

    def __init__(self, type):
        self.type = type
        self.coords = (0.0, 0.0, 0.0)

    # Overload the __repr__ operator to make printing simpler.
    def __repr__(self):
        return self.type

class Residue:
    """
    A simple class for an amino acid residue
    """

    def __init__(self, type, number):
        self.type = type
        self.number = number
        self.atoms = []

    # Overload the __repr__ operator to make printing simpler.
    def __repr__(self):
        return "{0} {1}".format(self.type, self.number)

class ActiveSite:
    """
    A simple class for an active site
    """

    def __init__(self, name):
        self.name = name
        self.residues = []
        self.metrics = [None,None,None]

    # Overload the __repr__ operator to make printing simpler.
    def __repr__(self):
        return self.name
        
    def get_raw_metrics(self):
        if self.metrics[0] is None:
            self._individual_metrics()
        return self.metrics[0] 
    
    def get_norm_metrics(self,redo_mean_dev=False):
        if self.metrics[2] is None or redo_mean_dev:
            self._individual_metrics()
            self._flatten_metrics()
            self._normalize_metrics(redo_mean_dev)
        return self.metrics[2]
    
    def _individual_metrics(self,printMe=None):
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
        numRes = len(self.residues)
    
        for residue in self.residues:
            res_type[chem_properties[residue.type.strip()]] += 1
            residue_charge = np.append(residue_charge, charge[chem_properties[residue.type.strip()]])
            residue_polarity = np.append(residue_polarity, polarity[chem_properties[residue.type.strip()]])
    
        for key in res_type.keys():
            res_type[key] /= numRes
        
        self._update_dict(metrics,
            ('meanChar','varChar','meanPol','varPol','nonpolarFrac','polarFrac','acidicFrac','basicFrac'),
            (np.mean(residue_charge),np.var(residue_charge),
                np.mean(residue_polarity),np.var(residue_polarity),
                res_type['nonpolar'],res_type['polar'],res_type['acidic'],res_type['basic'])
            )
      
        #now calculate the fraction of atoms of carbon/nitrogen/oxygen in the active zone
        numAtoms = 0
        
        for residue in self.residues:
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
    
        self._update_dict(metrics,
            ('elemFracC','elemFracN','elemFracO','numAtoms','numRes'),
            (atom_type['C'],atom_type['N'],atom_type['O'],numAtoms,numRes)
            )
    
        #now calculate the "center of mass" for each residue.
        #then calculate the mean/var of the pairwise Euclidean distance between residue COMs.  

        for residue in self.residues:
            numAtoms = 0
            centerOfMass = np.zeros(3)
        
            for atom in residue.atoms:
                numAtoms += 1
                centerOfMass += np.asarray(atom.coords)
        
            if not numAtoms==0:
                centerOfMass /= numAtoms
        
            residue.com = tuple(centerOfMass)
    
        residue_dist = []
    
        for pair in combo(self.residues,2):
            distance = np.linalg.norm(np.asarray(pair[0].com)-np.asarray(pair[1].com))
            residue_dist.append(distance)  
       
        self._update_dict(metrics,
            ('meanResDist','varResDist'),
            (np.mean(residue_dist),np.var(residue_dist))
            )
    
        if printMe:    
            for metric in metrics.keys():
                print('{0}\n\t\t{1}\n'.format(metric,metrics[metric]))
        
        #store
        self.metrics[0] = metrics
        
        return self.metrics[0]
    
    def _flatten_metrics(self):
        """
        Flattens metric dictionary into a numpy array
        """
    
        #NOTE: When converting into array, follows the metric_names list
        #so that you can consistently compare numpy arrays element-wise

        metric_names = ['meanChar','varChar','meanPol','varPol',
            'elemFracC','elemFracN','elemFracO','numAtoms','numRes','meanResDist','varResDist', 
            'nonpolarFrac','polarFrac','acidicFrac','basicFrac']
    
        metric_arr = np.ndarray(0)
    
        for metric in metric_names:
            metric_arr = np.append(metric_arr,self.metrics[0][metric])
            
        #store
        self.metrics[1] = metrics_arr
        
        return self.metrics[1]
    
    def _normalize_metrics(self,redo_mean_dev=False):
        """
        Normalizes metric arrays to saved (or custom) dataset
        """
    
        #if you want to recalculate mean/devs
        if redo_mean_dev:
            io.gen_mean_dev_normalizations()
    
        means = io.load_obj('means_arr')
        devs = io.load_obj('devs_arr')
    
        #normalize
        
        self.metrics[2] = (self.metrics[1] - means)/devs

        return self.metrics[2]
    
    def _update_dict(self,dictionary,keys,values):
        """
        Update dictionary given lists of keys and values
    
        Input: dictionary, keys, values
        Output: None (dictionary updated automatically)
        """

        dictionary.update(dict(zip(keys,values)))
        