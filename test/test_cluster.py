from clustering import cluster
from clustering import io
import os
from random import random, seed
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
    
    tolerance = 0.001
    
    # get random sites from activeSites    
    seed()
    rand1 = round(random()*len(activeSites))
    rand2 = round(random()*len(activeSites))
    
    while rand2 == rand1:
        rand2 = round(random()*len(activeSites))

    #Range
    sim_y, metrics_y = cluster.compute_similarity(activeSites[rand1], activeSites[rand2])
    assert (0-tolerance) <= sim_y <= (1+tolerance) #tolerance for rounding errors
    
    # Reflexive
    sim_ref, metrics_ref = cluster.compute_similarity(activeSites[rand1],activeSites[rand1])
    assert sim_ref > (1-tolerance) # tolerance for rounding errors
    
    # Symmetric
    sim_x, metrics_x = cluster.compute_similarity(activeSites[rand1], activeSites[rand2])
    sim_y, metrics_y = cluster.compute_similarity(activeSites[rand2], activeSites[rand1])
    assert abs(sim_x - sim_y) < (0+tolerance) # tolerance for rounding errors
    
    
    # get known sites from activeSites
    ind_a = cluster.find_active_sites(activeSites,sites[0])[0]
    ind_b = cluster.find_active_sites(activeSites,sites[1])[0]
    ind_c = cluster.find_active_sites(activeSites,sites[2])[0]

    sim_ab, met_ab = cluster.compute_similarity(activeSites[ind_a],activeSites[ind_b])
    sim_bc, met_bc = cluster.compute_similarity(activeSites[ind_b],activeSites[ind_c])
    sim_ac, met_ac = cluster.compute_similarity(activeSites[ind_a],activeSites[ind_c])
    
    # Similar ones are closer than dissimilar ones
    assert sim_ab - sim_ac > tolerance
    assert sim_ab - sim_bc > tolerance
    
    # Triangle inequality
    assert (1-sim_ab)-(1-sim_ac+1-sim_bc) < tolerance # small tolerance for rounding errors
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
