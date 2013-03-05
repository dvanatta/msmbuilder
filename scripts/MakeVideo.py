#!/home/dvanatta/epd-7.1-1-rh5-x86_64/bin/python
# This file is part of MSMBuilder.
#
# Copyright 2011 Stanford University
#
# MSMBuilder is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

import sys
import os
import MSMLib 
import Trajectory 
import scipy 
import numpy
from msmbuilder import Project
from msmbuilder import Serializer
from msmbuilder import arglib
from MSMLib import Sample

def run(project, assignments, start, stop, output, format, matrix):
    #arglib.die_if_path_exists(output)
    num_states = max(assignments.flatten()) + 1
    print "sampling, stop=", stop 
    T=scipy.io.mmread(matrix)
    sample=MSMLib.Sample(T,start, 50000)
    R1 = Trajectory.Trajectory.LoadFromPDB("./active.pdb")
    sample2=[]
    for i in range(len(sample)):
        sample2.append(sample[i])
        if sample[i]==stop:
	    print "SAMPLE=STOP AT step", i
            break
    print sample, len(sample) 
    print sample2, len(sample2)
    ind=0
    for i in range(len(sample2)): 
	#print sample2[i] 
	conf = project.GetRandomConfsFromState(assignments, sample2[i], 1)
    #    print "conf=",conf 
        if i==0:  #prints first frame as pdb to feed into vmd or w/e:
            R1["XYZList"] = conf 
            outputname=output+str(len(sample2))+".pdb"
	    ind=ind+1
	    R1.SaveToPDB(outputname)
        if format == 'xtc' and i>0:
            R1['XYZList']=numpy.vstack((R1['XYZList'],conf))
     #       print "lenR1=", len(R1['XYZList'])
        if format == 'pdb':
            R1["XYZList"] = conf 
            outputname=output+str(ind)+".pdb"
	    ind=ind+1
	#    print outputname
	    R1.SaveToPDB(outputname)
    #elif format == 'lh5':
    #    states.SaveToLHDF(output)
    R1["XYZList"]=R1["XYZList"].astype(numpy.float32)
    R1.SaveToXTC(output+str(len(sample2))+".xtc")
    #else:
    #    raise ValueError('Unrecognized format')
    print "Wrote output to:", output

    return


if __name__ == "__main__":
    parser = arglib.ArgumentParser(description="""
Pulls a certain number of random conformations from each cluster. Returns these
as an HDF5/PDB/XTC file that contains one long chain of these conformations that looks
like a Trajectory. If you selected to sample N conformations from each cluster,
the first N conformations are from cluster 0, the next N from cluster 1, etc.

Output default: XRandomConfs.lh5, where X=Number of Conformations.""")
    parser.add_argument('project')
    parser.add_argument('assignments', default='Data/Assignments.Fixed.h5')
    parser.add_argument('matrix', default='Data/tProb.mtx')
    parser.add_argument('output', description="""The name of the RandomConfs
        trajectory (.lh5) to write. XRandomConfs.lh5, where X=Number of
        Conformations.""", default='movie')
    parser.add_argument('stop', description='''Index of stop state''',type=int)
    parser.add_argument('start', description='''Index of start state
        ''', type=int)
    parser.add_argument('format', description='''Format to output the data in. Note
        that the PDB format is uncompressed and not efficient. For XTC, you can view
        the trajectory using your project's topology file''', default='pdb',
        choices=['pdb', 'xtc', 'lh5'])    
    args = parser.parse_args()
    
#    if args.output == 'XRandomConfs':
#            args.output = '%dRandomConfs.%s' % (args.conformations_per_state, args.format)

    run(args.project, args.assignments['Data'], args.start, args.stop,
        args.output, args.format, args.matrix)
