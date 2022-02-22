#!/bin/bash
#SBATCH -p fat -c 18 -J count_size --mail-user=markus.rainer.schmidt@gmail.com --mail-type END --time=240:00:00 -o slurm_count_size-%j.out

python3 count_size.py


# two options for the datalayout
# 1) 'naive'
#    btree where each point is encoded via the coordinates as the key and a pointer to the text or nullptr & the 
#    continuous sum.
# 2) 'complex'
#    Split individual coordinates of the dimensions up and store them seperatly in a 'nested' vector.
#    since i don't want the nested vector in multiple files on disk i will concatenate all nested vectors to one.
#    so child nodes will not have pointers to the sub vector but the indices of the sub vector in the concatenated one.
#    the final layer's indices will then refer to the index of the points that is in this bin and the continuous sum.
#    the points will be in a seperate vector. if there are multiple points in the exact same location the last layer 
#    will have them one after the other.
#    The keys will be the positions in the given dimensions (order of dimensions is predefined so no need to store it)
#    optional: to make binarry search more efficient via Eytzinger Binary Search.
