#!/usr/bin/python
# -*- coding: utf-8 -*-
# Kimmo Sääskilahti, 2015

from __future__ import division
import numpy as np
from SHCPostProcFv_kys import SHCPostProcFv

def main(fileprefix):

    # fileprefix='140415a'
    # Datafolder 
    dataFolder='DATA'
    outputFolder=dataFolder+'/'+fileprefix+'_tar'
    # Post-processor searches/saves file "KijFilePrefix.Kij.npy"
    
    # Create the data folder
    from subprocess import call
    command=["mkdir","-p",outputFolder]
    print " ".join(command)
    call(command)
    
    if 1:
        dt_md=0.002 # Timestep used in MD, affects the frequency grid
        widthWin=0.5 # Width of the Daniell smoothing window
    else:
        dt_md=1e-15 # Timestep used in MD, affects the frequency grid
        widthWin=0.25e12
    # The force-velocity dump file from LAMMPS
    fileFV=fileprefix+'.Fv.dat' 
    # The compactly formatted velocity file, produced using a C++ script if not found
    fileCompactFv=fileprefix+'.Fv.dat.compact'

    # Correct the units
    scaleFactor=1 # 

    # Prepare the post-processor
    pP=SHCPostProcFv(fileCompactFv,
                   dt_md=dt_md,scaleFactor=scaleFactor,
                   LAMMPSDumpFile=fileFV,widthWin=widthWin,
                   NChunks=200,chunkSize=25000,
                   backupPrefix=fileprefix,
                   reCalcFv=False)
    # Post-process
    pP.postProcess() # All variables will be contained in the object pP
    
    # Various output options
    
    # Pickling the post-processing instance
    import cPickle as pickle
    with open(outputFolder+'/'+fileprefix+'_PP_kys.pckl','w') as f:
        pickle.dump(pP,f)
    
    # Saving into numpy files
    np.save(outputFolder+'/'+fileprefix+'_oms.npy',pP.oms_fft)
    np.save(outputFolder+'/'+fileprefix+'_SHC.npy',pP.SHC_smooth)

    # Saving to file
    print "Saving to file "+outputFolder+'/'+fileprefix+'_SHC.txt'
    np.savetxt(outputFolder+'/'+fileprefix+'_SHC.txt',np.column_stack((pP.oms_fft,pP.SHC_smooth)))

    # Tar the output folder without the DATA folder 
    command=["tar","-czvf",fileprefix+'_tar.tgz','--directory='+dataFolder,outputFolder.strip(dataFolder+'/')]
    print " ".join(command)
    call(command)

    # Plotting if available
    # import matplotlib.pylab as plt
    # plt.plot(pP.oms_fft/(2*np.pi*1.0e12),pP.SHC_smooth)
    # plt.xlabel('Frequency (THz)')
    # plt.ylabel('Spectral current')
    # plt.savefig(fileprefix+'_SHC.eps')

if __name__=='__main__':

    import argparse

    parser=argparse.ArgumentParser()
    parser.add_argument("filePrefix",help="The used file prefix.")
    # parser.add_argument("KijFilePrefix",nargs='?',help="The prefix for the used force constant file, script searches for *.Kij.npy where *=argument (optional, default=fileprefix)")
    args=parser.parse_args()
    filePrefix=args.filePrefix
    
    main(fileprefix=filePrefix)
