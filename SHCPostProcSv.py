# -*- coding: utf-8 -*-
# Kimmo Sääskilahti, 2015
from __future__ import division
import numpy as np

class SHCPostProcSv(object):
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
        
        self.SHC_x_smooth=None
        self.SHC_x_smooth2=None

        self.DoS_smooth=None
        self.DoS_average=None

        for key,value in args.items():
            if not hasattr(self,key):
                raise ValueError, "Invalid argument " + key + " to PostProc!"
            print "Using the value "+key+"="+str(value)+"."
            setattr(self,key,value)
        # Smoothing width for ky-kz calculations
        self.widthWin2=self.widthWin/10.0
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
        smooth[0:round(Nwindow/2.0)]=func[0:round(Nwindow/2.0)]
        return smooth

    def initialize_arrays(self,Nfreqs):
        # Initialize the spectral heat current arrays
        self.SHC_smooth=np.zeros(Nfreqs)
        self.SHC_smooth2=np.zeros(Nfreqs)
        self.SHC_average=np.zeros(Nfreqs)

        self.SHC_x_smooth=np.zeros(Nfreqs)
        self.SHC_x_smooth2=np.zeros(Nfreqs)

        self.DoS_smooth=np.zeros(Nfreqs)
        self.DoS_average=np.zeros(Nfreqs)  


    def update_averages(self,k,SHC,SHC_x,DoS):
        SHC2=self._smoothen(self.df,SHC,self.widthWin)
        SHC2_x=self._smoothen(self.df,SHC_x,self.widthWin)
        DoS2=self._smoothen(self.df,DoS,self.widthWin)
        self.SHC_smooth=(k*self.SHC_smooth+SHC2)/(k+1.0)

        self.SHC_x_smooth=(k*self.SHC_x_smooth+SHC2_x)/(k+1.0)
        # The square
        self.SHC_smooth2=(k*self.SHC_smooth2+SHC2**2)/(k+1.0)
        self.SHC_x_smooth2=(k*self.SHC_x_smooth2+SHC2_x**2)/(k+1.0)
        # The non-smoothened average
        self.SHC_average=(k*self.SHC_average+SHC)/(k+1.0)

        self.DoS_smooth=(k*self.DoS_smooth+DoS2)/(k+1.0)
        self.DoS_average=(k*self.DoS_average+DoS)/(k+1.0)


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
  
        self.initialize_arrays(Nfreqs)
        exitFlag=False

        for k in np.arange(self.NChunks): # Start the iteration over chunks
#        for k in range(0,2): # Start the iteration over chunks
            print "Chunk %d/%d." % (k+1,self.NChunks)
            # Read a chunk of velocitites
            print "Reading forces and velocities..."
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
            print "Fourier transforming the arrays..."
            velFFT=np.fft.rfft(velArray,axis=1)
            velFFT*=self.sampleTimestep

            forceFFT=np.fft.rfft(forceArray,axis=1)
            forceFFT*=self.sampleTimestep
   
            print "Calculating the spectral heat current..."
            SHC=2*np.real(np.sum(forceFFT*np.conj(velFFT),axis=0))
            SHC_x=2*np.real(np.sum(forceFFT[0::3,:]*np.conj(velFFT[0::3,:]),axis=0))
            DoS=np.real(np.sum(velFFT*np.conj(velFFT),axis=0))

            # Normalize correctly
            SHC/=(self.chunkSize*self.sampleTimestep)
            SHC_x/=(self.chunkSize*self.sampleTimestep)
            DoS/=(self.chunkSize*self.sampleTimestep)

            # Change units             
            SHC*=self.scaleFactor
            
                      
            # For smoothing used in the update           
            self.df=(self.oms_fft[1]-self.oms_fft[0])/(2*np.pi)
            
            if not exitFlag: # If Nfreqs has changed, the running averaging cannot be performed
                self.update_averages(k,SHC,SHC_x,DoS)                

                if self.backupPrefix is not None:
                    np.save(self.backupPrefix+'_backup_oms.npy',self.oms_fft)
                    np.save(self.backupPrefix+'_backup_SHC.npy',self.SHC_smooth)
                    import cPickle as pickle
                    with open(self.backupPrefix+'_run_PP.pckl','w') as pf:
                        pickle.dump(self,pf)
            elif exitFlag and k==0: # First chunk and new chunk size, needs re-initializing the vectors as Nfreqs may have changed               
                self.SHC_smooth=SHC
                self.SHC_smooth2=SHC**2
                self.SHC_x_smooth=SHC_x
                self.SHC_x_smooth2=SHC_x**2
                self.SHC_average=SHC_orig
                
                
                self.DoS_smooth=DoS
                self.DoS_average=DoS_orig
                
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
            samplevar=(self.NChunks/(self.NChunks-1.0))*(self.SHC_x_smooth2-self.SHC_x_smooth**2)
            self.SHC_x_error=np.sqrt(samplevar)/np.sqrt(self.NChunks)
        else:
            self.SHC_error=None
            self.SHC_x_error=None
            
        print "Finished post-processing."
