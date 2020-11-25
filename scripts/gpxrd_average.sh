#!/bin/bash

# This script replaces the average functionality of the gpxrd.py script.
# Instead of handling a whole trajectory in memory, it splits the trajectory in snapshots (prepare).
# Then, the pattern for each of snapshots is calculated (sp).
# Finally, the average pattern is calculated (post).
#
# This script creates two folders, cifs/ and frames/,
# containing the cif files of the snapshots and the PXRD patterns of the snapshots respectively.
#
# This script also output a log file containing the heuristic values for the averaged pattern.
# The average pattern can be found in avg.dat, and favg.dat,
# which differ only at low 2theta values as favg.dat contains the full 2theta range,
# whereas avg.dat only contains the same 2theta values as the reference pattern.

# Example usage: $ gpxrd.sh traj_md.h5 50
# You can adapt the code below to suit your needs.

#   Change plot=True to plot=name to save the plot to a file instead of a pop-up
#   Adapt exp.tsv to the name of your reference pattern (should have two columns, 2theta and intensity or counts)
#   Adapt runuptime based on how equilibrated your initial structure is.


traj=$1
snapshots=$2

gpxrd.py prepare ${traj} --runuptime=1000 --snapshots=${snapshots}

for i in $(seq 0 $(( ${snapshots}-1 ))); do
    gpxrd.py sp $i exp.tsv > log.log
done

gpxrd.py post exp.tsv --plot=True
