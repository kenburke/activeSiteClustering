from clustering import cluster, io, similarity
import os
from random import random, seed
from math import floor
import pytest


@pytest.mark.parametrize("sites", [
    (["276","1806","70919"]), #first two small and basic, last is big and acidic
    (["3733","3458","71389"]) 
])
def test_similarity(sites):
    """
    Test similarity function for expected identities and behavior with test data
    
    Input: a list of 3 sites, the first two being more similar than the third.
    """
    
    activeSites = io.read_active_sites('data/')
    
    tolerance = 0.001 # for rounding errors
    
    # get random sites from activeSites    
    seed()
    rand1 = floor(random()*len(activeSites))
    rand2 = floor(random()*len(activeSites))
    
    while rand2 == rand1:
        rand2 = floor(random()*len(activeSites))

    #Range
    sim_y = similarity.compute_similarity(activeSites[rand1], activeSites[rand2])
    assert (0-tolerance) <= sim_y <= (1+tolerance)
    
    # Reflexive
    sim_ref = similarity.compute_similarity(activeSites[rand1],activeSites[rand1])
    assert sim_ref > (1-tolerance) 
    
    # Symmetric
    sim_x = similarity.compute_similarity(activeSites[rand1], activeSites[rand2])
    sim_y = similarity.compute_similarity(activeSites[rand2], activeSites[rand1])
    assert abs(sim_x - sim_y) < (0+tolerance)
    
    
    # get known sites from activeSites
    ind_a = similarity.find_active_sites(activeSites,sites[0])[0]
    ind_b = similarity.find_active_sites(activeSites,sites[1])[0]
    ind_c = similarity.find_active_sites(activeSites,sites[2])[0]

    sim_ab = similarity.compute_similarity(activeSites[ind_a],activeSites[ind_b])
    sim_bc = similarity.compute_similarity(activeSites[ind_b],activeSites[ind_c])
    sim_ac = similarity.compute_similarity(activeSites[ind_a],activeSites[ind_c])
    
    # Similar ones are closer than dissimilar ones
    assert sim_ab - sim_ac > tolerance
    assert sim_ab - sim_bc > tolerance
    
    # Triangle inequality
    assert (1-sim_ab)-(1-sim_ac+1-sim_bc) < tolerance 
    assert (1-sim_ac)-(1-sim_ab+1-sim_bc) < tolerance
    assert (1-sim_bc)-(1-sim_ab+1-sim_ac) < tolerance

def test_partition_clustering():
    # tractable subset
    pdb_ids = [276, 4629, 10701]

    active_sites = []
    for id in pdb_ids:
        filepath = os.path.join("data", "%i.pdb"%id)
        active_sites.append(io.read_active_site(filepath))

    # update this assertion
    assert cluster.cluster_by_partitioning(active_sites) == []

def test_hierarchical_clustering():
    # tractable subset
    pdb_ids = [276, 4629, 10701]

    active_sites = []
    for id in pdb_ids:
        filepath = os.path.join("data", "%i.pdb"%id)
        active_sites.append(io.read_active_site(filepath))

    # update this assertion
    assert cluster.cluster_hierarchically(active_sites) == []
