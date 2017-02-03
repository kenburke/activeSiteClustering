# Clustering Homework

[![Build
Status](https://travis-ci.org/kenburke/clustering.svg?branch=master)](https://travis-ci.org/kenburke/clustering)

Clustering homework assignment, with continuous integration testing.

## assignment

1. Implement a similarity metric
2. Implement a clustering method based on a partitioning algorithm
3. Implement a clustering method based on a hierarchical algorithm
4. Implement a function to measure the quality of the clustering
5. Implement a function to compare the two clusterings to each other


## structure

The main file that you will need to modify is `cluster.py` and the corresponding `test_cluster.py`. `utils.py` contains helpful classes that you can use to represent Active Sites. `io.py` contains some reading and writing files for interacting with PDB files and writing out cluster info.

```
.
├── README.md
├── data
│   ...
├── clustering
│   ├── __init__.py
│   ├── __main__.py
│   ├── cluster.py
│   ├── io.py
│   └── utils.py
└── test
    ├── test_cluster.py
    └── test_io.py
```
## usage

To use the package, first create a conda environment

```
conda env create
```

to install all dependencies outlined in `environment.yml`. Then activate the env

```
source activate activeSiteClustering
```

Then the package's main function (located in `clustering/__main__.py`) 
can be run as follows

```
python -m clustering [-P| -H] data test.txt
```

where ``-P`` and ``-H`` use partition or hierarchical clustering, respectively.

## testing

Testing is as simple as running

```
python -m pytest
```

from the root directory of this project.


## contributors

Original design by Scott Pegg. Refactored and updated by Tamas Nagy.
