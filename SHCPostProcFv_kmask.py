# -*- coding: utf-8 -*-
# Kimmo Sääskilahti, 2015
from __future__ import division
import numpy as np

class SHCPostProcFv(object):
    '''
    Post-process the data produced using LAMMPS Molecular Dynamics simulation to calculate the spectral heat current.

    The forces and velocities are read from the "compact" file produced with the C++-code compactify_Fv.cpp from a LAMMPS dump file. If the file does not exist, it is produced by calling the binary "compactify_Fv", which must be found in the environment's $PATH.

    Minimal usage in Python:
      pP=SHCPostProcFv(compactFvFile) # See the documentation for arguments below
      pP.postProcess() # Calculate the heat current spectrum

    Public attributes:
      SHC_smooth (numpy float array): The chunk-averaged, smoothened spectral heat current
      SHC_smooth2 (numpy float array): Square of the chunk-averaged, smoothened spectral heat current, used for estimating the error from the between-chunk variance
      SHC_average (numpy float array): The chunk-averaged spectral heat current without smoothing
      SHC_error (numpy float array): The estimated error from the between-chunk variance, None if only one chunk evaluated
      oms_fft (numpy float array): The angular frequency grid (in the units of Hz if dt_md is given in the units of seconds in the initialization)
    '''

    def __init__(self,compactFvFile,reCalcFv=False,**args):
        '''
        Positional arguments:
           compactVelocityFile (str): The file where the velocities are read. Produced using the binary compactify_vels if the file does not exist. In this case, you must also supply the keyword argument LAMMPSDumpFile containing the velocities produced using LAMMPS. 

        Keyword arguments:

           dt_md (float): Timestep used in the NEMD simulation (seconds), used for inferring the sampling timestep and the frequency grid (default 1.0)
           scaleFactor (float): Multiply the spectral heat current by this factor to convert to correct units (default 1.0)
           LAMMPSDumpFile (str): Use this velocity dump file produced by LAMMPS for post-processing, needed only if the compact velocity file cannot be found (default None)
           widthWin (float): Use this width for the smoothing window (Hz) (default 1.0)
           chunkSize (int): Used chunk size for reading the velocities, affects the frequency grid (default 50000). Performing FFT is faster if chunkSize is a power of 2.
           NChunks (int): The number of chunks to be read, this should be set to a sufficiently large value if the whole velocity file should be read (default 20)
           backupPrefix (str): Prefix for pickling the post-processing object to a file after the read of each chunk (default None)
           hstep (float): The displacement used in calculating the force constants by the finite-displacement method (default 0.001)
           reCalcFv (boolean): If set to True, the re-creation of the compact force-velocity file is enforced (default False)
        '''

        # Attributes set by positional arguments
        self.compactFvFile=compactFvFile
        
        # Attributes set by keyword parameters below
        self.dt_md=1.0 # Default
        self.scaleFactor=1.0 # Default
        self.LAMMPSDumpFile=None
        self.widthWin=1.0 # Default
        self.chunkSize=50000
        self.NChunks=20
        self.backupPrefix=None
        self.hstep=0.001
        self.reCalcFv=reCalcFv

        self.NL=None
        self.NR=None
        self.SHC_smooth=None
        self.SHC_smooth2=None
        self.SHC_average=None
        self.SHC_error=None
        self.oms_fft=None
        self.SHC_k0_smooth=None
        self.SHC_k0_smooth2=None
        self.SHC_k0_average=None
        self.SHC_k0_error=None

        for key,value in args.items():
            if not hasattr(self,key):
                raise ValueError, "Invalid argument " + key + " to PostProc!"
            print "Using the value "+key+"="+str(value)+"."
            setattr(self,key,value)

        import os
        if self.reCalcFv or not os.path.isfile(self.compactFvFile): # Check if the force-velocity file exists
            # Check that the LAMMPS Dump file exists
            if self.LAMMPSDumpFile is None or not os.path.isfile(self.LAMMPSDumpFile):
                raise ValueError, "You must give the LAMMPS velocity dump file as an argument to create the file "+self.compactFvFile+"!"
            #print self.compactVelocityFile + " does not exist, creating by reading from file " + self.LAMMPSDumpFile
            # Run the C++ script
            self._compactFv(self.LAMMPSDumpFile,self.compactFvFile)
        else:
            print self.compactFvFile + " exists, using the file for post-processing."

    def __enter__(self):
        return self

    def __exit__(self,t1,t2,t3):
        return False

    def _compactFv(self,fileFv,finalFileFv):
        from subprocess import call
        command=["compactify_Fv",fileFv,finalFileFv]
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
        
        print "Reading the compact force-velocity file "+self.compactFvFile+"."
        f=open(self.compactFvFile,'r')
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
        # Initialize the spectral heat current arrays
        self.SHC_smooth=np.zeros(Nfreqs)
        self.SHC_smooth2=np.zeros(Nfreqs)
        self.SHC_average=np.zeros(Nfreqs)

        self.SHC_k0_smooth=np.zeros(Nfreqs)
        self.SHC_k0_smooth2=np.zeros(Nfreqs)
        self.SHC_k0_average=np.zeros(Nfreqs)

        exitFlag=False

        for k in np.arange(self.NChunks): # Start the iteration over chunks
#        for k in range(0,2): # Start the iteration over chunks
            print "Chunk %d/%d." % (k+1,self.NChunks)
            # Read a chunk of velocitites
            FV_Array=np.fromfile(f,dtype=np.dtype('f8'),count=2*self.chunkSize*NDOF,sep=" ")
            # print np.size(FV_Array)

            if np.size(FV_Array)==0:
                print "Nothing more to read, finishing."
                self.NChunks=k-1
                break

            forceArray=FV_Array[0::2]
            velArray=FV_Array[1::2]
            
            print "Size of velocity array: %d" % (np.size(velArray))
            print "Size of force array: %d" % (np.size(forceArray))
            
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
            forceArray=np.reshape(forceArray,(NDOF,self.chunkSize),order='F')
            
            # FFT with respect to the second axis (NOTE THE USE OF RFFT)
            velFFT=np.fft.rfft(velArray,axis=1)
            velFFT*=self.sampleTimestep

            forceFFT=np.fft.rfft(forceArray,axis=1)
            forceFFT*=self.sampleTimestep
   
            # Spectral heat current for the specific chunk
            SHC=np.zeros(Nfreqs)

            #for ki in range(1,Nfreqs): # Skip the first one with zero frequency
            #    SHC[ki]=2*np.real(np.dot(forceFFT[:,ki],np.conj(velFFT[:,ki])))
                # SHC[ki]=2*np.real(np.dot(forceFFT[1::3,ki],np.conj(velFFT[1::3,ki])))
            SHC=2*np.real(np.sum(forceFFT*np.conj(velFFT)),axis=0)

            # Normalize correctly
            SHC/=(self.chunkSize*self.sampleTimestep)

            # Change units             
            SHC*=self.scaleFactor
      
            
            # print "Found file "+self.kmask_file+", proceeding."
            vels_k0=np.zeros((12,Nfreqs),dtype=complex)
            Fs_k0=np.zeros((12,Nfreqs),dtype=complex)
            SHC_k0=np.zeros(Nfreqs)
            if NDOF%12!=0:
                print "Number of atoms at the interface not divisible by 12, cannot calculate heat current for transverse wavevector zero."
            else:
                print "Proceeding to calculate contribution of ky=kz=0, this only makes sense if only a single unit cell is included at the interface."
                for kj in range(0,12):
                    # np.exp(1j*(ky*kmask_kys[kj::12]+kz*kmask_kzs[kj::12])
                    # Wavevector decomposition for transverse wavevector zero
                    vels_k0[kj,:]=np.sum(velFFT[kj::12,:],axis=0)
                    Fs_k0[kj,:]=np.sum(forceFFT[kj::12,:],axis=0)

                
                for ki in range(1,Nfreqs): # Skip the first one with zero frequency
                    SHC_k0[ki]=2*np.real(np.dot(Fs_k0[0::1,ki],np.conj(vels_k0[0::1,ki])))
                    # SHC_k0[ki]=2*np.real(np.dot(Fs_k0[0::3,ki],np.conj(vels_k0[0::3,ki])))
                # Normalize correctly
                SHC_k0/=(self.chunkSize*self.sampleTimestep)

                # Change units             
                SHC_k0*=self.scaleFactor

            # daniellWindow=np.ones(np.ceil(self.widthWin*2*np.pi/(self.oms_fft[1]-self.oms_fft[0])))
            # daniellWindow/=np.sum(daniellWindow)

            SHC_orig=SHC.copy()
            SHC_k0_orig=SHC_k0.copy()
            # Smooth the value           
            # SHC=np.convolve(SHC,daniellWindow,'same')
            df=(self.oms_fft[1]-self.oms_fft[0])/(2*np.pi)
            SHC=self._smoothen(df,SHC,self.widthWin)
            SHC_k0=self._smoothen(df,SHC_k0,self.widthWin)
            

            if not exitFlag: # If Nfreqs has changed, the running averaging cannot be performed
                self.SHC_smooth=(k*self.SHC_smooth+SHC)/(k+1.0)
                # The square
                self.SHC_smooth2=(k*self.SHC_smooth2+SHC**2)/(k+1.0)
                # The non-smoothened average
                self.SHC_average=(k*self.SHC_average+SHC_orig)/(k+1.0)
                
                self.SHC_k0_smooth=(k*self.SHC_k0_smooth+SHC_k0)/(k+1.0)
                # The square
                self.SHC_k0_smooth2=(k*self.SHC_k0_smooth2+SHC_k0**2)/(k+1.0)
                # The non-smoothened average
                self.SHC_k0_average=(k*self.SHC_k0_average+SHC_k0_orig)/(k+1.0)

                if self.backupPrefix is not None:
                    np.save(self.backupPrefix+'_backup_oms.npy',self.oms_fft)
                    np.save(self.backupPrefix+'_backup_SHC.npy',self.SHC_smooth)
                    import cPickle as pickle
                    with open(self.backupPrefix+'_run_PP.pckl','w') as pf:
                        pickle.dump(self,pf)
            elif exitFlag and k==0: # First chunk and new chunk size, needs re-initializing the vectors as Nfreqs may have changed               
                self.SHC_smooth=SHC
                self.SHC_smooth2=SHC**2
                self.SHC_average=SHC_orig
                self.SHC_k0_smooth=SHC_k0
                self.SHC_k0_smooth2=SHC_k0**2
                self.SHC_k0_average=SHC_k0_orig
                self.NChunks=1
                break
            else: # This should never be reached
                assert False, "SHCPostProc should not reach here (exitFlag=True and k>0)."
                break

        # Calculate the error estimate at each frequency from the between-chunk variances
        if self.NChunks>1:
            print "Calculating error estimates."
            samplevar=(self.NChunks/(self.NChunks-1.0))*(self.SHC_smooth2-self.SHC_smooth**2)
            self.SHC_error=np.sqrt(samplevar)/np.sqrt(self.NChunks)
            samplevar=(self.NChunks/(self.NChunks-1.0))*(self.SHC_k0_smooth2-self.SHC_k0_smooth**2)
            self.SHC_k0_error=np.sqrt(samplevar)/np.sqrt(self.NChunks)
            
        else:
            self.SHC_error=None
            self.SHC_k0_error=None

        print "Finished post-processing."
