import sys
from .io import read_active_sites, write_clustering, gen_mean_dev_normalizations
from .cluster import cluster_by_partitioning, cluster_hierarchically

# Some quick stuff to make sure the program is called correctly
if len(sys.argv) < 4:
    print("Usage: python -m clustering [-P| -H] <pdb directory> <output file>")
    sys.exit(0)

# initializes normalization files for files in data/
gen_mean_dev_normalizations()

active_sites = read_active_sites(sys.argv[2])

# Choose clustering algorithm
if sys.argv[1][0:2] == '-P':
    print("Clustering using Partitioning method")
    clustering = cluster_by_partitioning(active_sites)
    write_clustering(sys.argv[3], clustering)

if sys.argv[1][0:2] == '-H':
    print("Clustering using hierarchical method")
    clusterings = cluster_hierarchically(active_sites)
    write_clustering(sys.argv[3], clusterings)
