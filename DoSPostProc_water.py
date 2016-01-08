# -*- coding: utf-8 -*-
# Kimmo Sääskilahti, 2015
from __future__ import division
import numpy as np

class DoSPostProc_water(object):
    '''
    Post-process the data produced using LAMMPS Molecular Dynamics simulation to calculate the spectral heat current.

    The forces and velocities are read from the "compact" file produced with the C++-code compactify_vels.cpp from a LAMMPS dump file. If the file does not exist, it is produced by calling the binary "compactify_vels", which must be found in the environment's $PATH.

    Minimal usage in Python:
      pP=DoSPostProc(compactVelsFile) # See the documentation for arguments below
      pP.postProcess() # Calculate the heat current spectrum

    Public attributes:
      DoS_smooth (numpy float array): The chunk-averaged, smoothened density of states
      DoS_average (numpy float array): The chunk-averaged density of states without smoothing
      oms_fft (numpy float array): The angular frequency grid (in the units of Hz if dt_md is given in the units of seconds in the initialization)
    '''

    def __init__(self,compactVelsFile,reCalcVels=False,**args):
        '''
        Positional arguments:
           compactVelsFile (str): The file where the velocities are read. Produced using the binary compactify_vels if the file does not exist. In this case, you must also supply the keyword argument LAMMPSDumpFile containing the velocities produced using LAMMPS. 

        Keyword arguments:

           dt_md (float): Timestep used in the NEMD simulation (seconds), used for inferring the sampling timestep and the frequency grid (default 1.0)
           scaleFactor (float): Multiply the spectral heat current by this factor to convert to correct units (default 1.0)
           LAMMPSDumpFile (str): Use this velocity dump file produced by LAMMPS for post-processing, needed only if the compact velocity file cannot be found (default None)
           widthWin (float): Use this width for the smoothing window (Hz) (default 1.0)
           chunkSize (int): Used chunk size for reading the velocities, affects the frequency grid (default 50000). Performing FFT is faster if chunkSize is a power of 2.
           NChunks (int): The number of chunks to be read, this should be set to a sufficiently large value if the whole velocity file should be read (default 20)
           backupPrefix (str): Prefix for pickling the post-processing object to a file after the read of each chunk (default None)
           reCalcVels (boolean): If set to True, the re-creation of the compact force-velocity file is enforced (default False)
        '''

        # Attributes set by positional arguments
        self.compactVelsFile=compactVelsFile
        
        # Attributes set by keyword parameters below
        self.dt_md=1.0 # Default
        self.scaleFactor=1.0 # Default
        self.LAMMPSDumpFile=None
        self.widthWin=1.0 # Default
        self.chunkSize=50000
        self.NChunks=20
        self.backupPrefix=None
        
        self.reCalcVels=reCalcVels

        self.NL=None
        self.NR=None
        
        self.oms_fft=None
        
        self.DoS_smooth=None
        self.DoS_average=None

        self.DoS_O_smooth=None
        self.DoS_O_average=None
        
        self.DoS_H_smooth=None
        self.DoS_H_average=None

        for key,value in args.items():
            if not hasattr(self,key):
                raise ValueError, "Invalid argument " + key + " to PostProc!"
            print "Using the value "+key+"="+str(value)+"."
            setattr(self,key,value)

        import os
        if self.reCalcVels or not os.path.isfile(self.compactVelsFile): # Check if the force-velocity file exists
            # Check that the LAMMPS Dump file exists
            if self.LAMMPSDumpFile is None or not os.path.isfile(self.LAMMPSDumpFile):
                raise ValueError, "You must give the LAMMPS velocity dump file as an argument to create the file "+self.compactVelsFile+"!"
            #print self.compactVelocityFile + " does not exist, creating by reading from file " + self.LAMMPSDumpFile
            # Run the C++ script
            self._compactVels(self.LAMMPSDumpFile,self.compactVelsFile)
        else:
            print self.compactVelsFile + " exists, using the file for post-processing."

    def __enter__(self):
        return self

    def __exit__(self,t1,t2,t3):
        return False

    def _compactVels(self,fileVels,finalFileVels):
        from subprocess import call
        command=["compactify_vels_notype",fileVels,finalFileVels]
        print "Running "+" ".join(command)
        call(command)

    def _smoothen(self,df,func,widthWin):
        Nwindow=np.ceil(widthWin/df)
        daniellWindow=np.ones(Nwindow)/Nwindow
        # daniellWindow/=np.sum(daniellWindow)
        # Smooth the value           
        smooth=np.convolve(func,daniellWindow,'same')
        return smooth

    def postProcess(self):
        
        print "Reading the compact velocity file "+self.compactVelsFile+"."
        f=open(self.compactVelsFile,'r')
        s=f.readline().split()
        self.NAtoms=int(s[1])
        
        print "NAtoms=%d." % (self.NAtoms)
        
        s=f.readline()
        print s
        s=s.split()
        self.sampleTimestep=int(s[1])*self.dt_md

        s=f.readline() # Atom ids:
        print s
        # Read the atom ids
        indArray=np.fromfile(f,dtype=int,count=self.NAtoms,sep=" ")
        # print indArray
        
        s=f.readline() # ------
        print s
  
        # Total number of degrees of freedom
        NDOF=3*(self.NAtoms)

        self.oms_fft=np.fft.rfftfreq(self.chunkSize,d=self.sampleTimestep)*2*np.pi
        Nfreqs=np.size(self.oms_fft)
        # Initialize the arrays
        self.DoS_smooth=np.zeros(Nfreqs)
        self.DoS_average=np.zeros(Nfreqs)

        self.DoS_O_smooth=np.zeros(Nfreqs)
        self.DoS_O_average=np.zeros(Nfreqs)

        self.DoS_H_smooth=np.zeros(Nfreqs)
        self.DoS_H_average=np.zeros(Nfreqs)

        exitFlag=False

        for k in np.arange(self.NChunks): # Start the iteration over chunks
#        for k in range(0,2): # Start the iteration over chunks
            print "Chunk %d/%d." % (k+1,self.NChunks)
            # Read a chunk of velocitites
            FV_Array=np.fromfile(f,dtype=np.dtype('f8'),count=1*self.chunkSize*NDOF,sep='\n')
            # print np.size(FV_Array)

            if np.size(FV_Array)==0:
                print "Nothing more to read, finishing."
                self.NChunks=k-1
                break

            velArray=FV_Array
            
            print "Size of velocity array: %d" % (np.size(velArray))
            # print "Size of force array: %d" % (np.size(forceArray))
            
            # Prepare for exit if the read size does not match the chunk size
            
            if np.size(velArray)!=self.chunkSize*NDOF:
                # Reaching the end of file  
                print "Reaching the end of file!"
                self.chunkSize=int(np.size(velArray)/NDOF)             
                if k>0: # Not the first chunk
                    self.NChunks=k-1
                    break
                else:
                    exitFlag=True
                    self.oms_fft=np.fft.rfftfreq(self.chunkSize,d=self.sampleTimestep)*2*np.pi
                    Nfreqs=np.size(self.oms_fft)
                    print "Changing chunk size to "+str(int(np.size(velArray)/NDOF))+"!"
                    
 
            # Reshape the array so that each row corresponds to different degree of freedom (e.g. particle 1, direction x etc.)
            velArray=np.reshape(velArray,(NDOF,self.chunkSize),order='F')
            
            # FFT with respect to the second axis (NOTE THE USE OF RFFT)
            velFFT=np.fft.rfft(velArray,axis=1)
            velFFT*=self.sampleTimestep

            #forceFFT=np.fft.rfft(forceArray,axis=1)
            #forceFFT*=self.sampleTimestep
   
            # Spectral heat current for the specific chunk
            # SHC=np.zeros(Nfreqs)

            #for ki in range(1,Nfreqs): # Skip the first one with zero frequency
            #    SHC[ki]=2*np.real(np.dot(forceFFT[:,ki],np.conj(velFFT[:,ki])))
                # SHC[ki]=2*np.real(np.dot(forceFFT[1::3,ki],np.conj(velFFT[1::3,ki])))
            # SHC=2*np.real(np.sum(forceFFT*np.conj(velFFT),axis=0))
            DoS=np.real(np.sum(velFFT*np.conj(velFFT),axis=0))
            DoS_O=16*np.real(np.sum(velFFT[0::9,:]*np.conj(velFFT[0::9,:]),axis=0))
            DoS_O+=16*np.real(np.sum(velFFT[1::9,:]*np.conj(velFFT[1::9,:]),axis=0))
            DoS_O+=16*np.real(np.sum(velFFT[2::9,:]*np.conj(velFFT[2::9,:]),axis=0))
            
            DoS_H=np.real(np.sum(velFFT[3::9,:]*np.conj(velFFT[3::9,:]),axis=0))
            for m in range(4,9):
                DoS_H+=np.real(np.sum(velFFT[m::9,:]*np.conj(velFFT[m::9,:]),axis=0))
            # Normalize correctly
            # SHC/=(self.chunkSize*self.sampleTimestep)
            DoS/=(self.chunkSize*self.sampleTimestep)
            DoS_O/=(self.chunkSize*self.sampleTimestep)
            DoS_H/=(self.chunkSize*self.sampleTimestep)
            

            DoS_orig=DoS.copy()
            DoS_O_orig=DoS_O.copy()
            DoS_H_orig=DoS_H.copy()
            # Smooth the value           
            # SHC=np.convolve(SHC,daniellWindow,'same')
            df=(self.oms_fft[1]-self.oms_fft[0])/(2*np.pi)
            # SHC=self._smoothen(df,SHC,self.widthWin)
            DoS=self._smoothen(df,DoS,self.widthWin)
            DoS_O=self._smoothen(df,DoS_O,self.widthWin)
            DoS_H=self._smoothen(df,DoS_H,self.widthWin)
            # SHC_k0=self._smoothen(df,SHC_k0,self.widthWin)
            # SHC_k0_x=self._smoothen(df,SHC_k0_x,self.widthWin)
            

            if not exitFlag: # If Nfreqs has changed, the running averaging cannot be performed

                self.DoS_smooth=(k*self.DoS_smooth+DoS)/(k+1.0)
                self.DoS_average=(k*self.DoS_average+DoS_orig)/(k+1.0)

                self.DoS_O_smooth=(k*self.DoS_O_smooth+DoS_O)/(k+1.0)
                self.DoS_O_average=(k*self.DoS_O_average+DoS_O_orig)/(k+1.0)

                self.DoS_H_smooth=(k*self.DoS_H_smooth+DoS_H)/(k+1.0)
                self.DoS_H_average=(k*self.DoS_H_average+DoS_H_orig)/(k+1.0)
                
                if self.backupPrefix is not None:
                    np.save(self.backupPrefix+'_backup_oms.npy',self.oms_fft)
                    np.save(self.backupPrefix+'_backup_DoS.npy',self.DoS_smooth)
                    import cPickle as pickle
                    with open(self.backupPrefix+'_run_PP.pckl','w') as pf:
                        pickle.dump(self,pf)
            elif exitFlag and k==0: # First chunk and new chunk size, needs re-initializing the vectors as Nfreqs may have changed               
                self.DoS_smooth=DoS
                self.DoS_average=DoS_orig 
                self.DoS_O_smooth=DoS_O
                self.DoS_O_average=DoS_O_orig 
                self.DoS_H_smooth=DoS_H
                self.DoS_H_average=DoS_H_orig 
  
                self.NChunks=1
                break
            else: # This should never be reached
                assert False, "DoSPostProc should not reach here (exitFlag=True and k>0)."
                break

        print "Finished post-processing."
