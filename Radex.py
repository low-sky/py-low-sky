import sys
import numpy as np
import os
import random
import subprocess
import shutil
import re

def WriteRadex(InputFileName,OutputFileName,Molecule = 'hco+', 
               UpperFrequencyGHz = 500, LowerFrequencyGHz = 50, 
               KineticTemperature = 10.0, NumberDensity = 1e4, 
               Colliders = 'H2', BackgroundTemperature = 2.73,
               ColumnDensity = 1e13, LineWidth = 1.0):
    InputFile = open(InputFileName,'w')
    InputFile.write(Molecule+'.dat\n')
    InputFile.write(OutputFileName+'\n')
    InputFile.write(str(LowerFrequencyGHz)+' '+str(UpperFrequencyGHz)+'\n')
    InputFile.write(str(KineticTemperature)+'\n')
    InputFile.write('1\n') #Number of Collider Species
    InputFile.write(Colliders+'\n')
    InputFile.write(str(NumberDensity)+'\n')
    InputFile.write(str(BackgroundTemperature)+'\n')
    InputFile.write(str(ColumnDensity)+'\n')
    InputFile.write(str(LineWidth)+'\n')
    InputFile.write('0') # End Calculation
    InputFile.close()
    return

def ReadRadex(OutputFile):
    template = np.zeros(1,dtype={'names':['Transition','Frequency',
                                           'ExcitationTemperature',
                                           'OpticalDepth',
                                           'RadiationTemperature',
                                           'Flux'],
                                  'formats':['S10','f4','f4','f4','f4','f4']})
    for line in OutputFile:
        if re.match('[0-9]',line[0]):
            RadexData = np.copy(template)
            parts = line.split()
            RadexData['Transition'] = parts[0]+parts[1]+parts[2]
            EnergyUp = float(parts[3])
            RadexData['Frequency'] = float(parts[4])
            RadexData['ExcitationTemperature'] = float(parts[6])
            RadexData['OpticalDepth'] = float(parts[7])
            RadexData['RadiationTemperature'] = float(parts[8])
            RadexData['Flux'] = float(parts[11])
            try:
                output = np.concatenate((output,RadexData))
            except NameError:
                output = RadexData
    return output

def RunRadex(Molecule = 'hco+', 
             UpperFrequencyGHz = 500, LowerFrequencyGHz = 50, 
             KineticTemperature = 10.0, NumberDensity = 1e4, 
             Colliders = 'H2', BackgroundTemperature = 2.73,
             ColumnDensity = 1e13, LineWidth = 1.0):
    
# Set location of radex binary
#    RadexBinary = '/Users/snc/astro/m33/jcmt/radex/bin/radex'
    RadexBinary = os.getenv("HOME")+'/astro/radex/bin/radex'
# TempName = os.tmpnam()
    dirname = '/var/tmp/RadexRun'+str(random.getrandbits(32))
    os.mkdir(dirname)
    InputFileName = dirname+'/radex.inp'
    OutputFileName = dirname+'/radex.out'
    WriteRadex(InputFileName,OutputFileName,Molecule = Molecule, 
               UpperFrequencyGHz = UpperFrequencyGHz, 
               LowerFrequencyGHz = LowerFrequencyGHz,
               KineticTemperature = KineticTemperature, 
               NumberDensity = NumberDensity,
               Colliders = Colliders,
               BackgroundTemperature = BackgroundTemperature,
               ColumnDensity = ColumnDensity, 
               LineWidth = LineWidth)
    InputFile = open(InputFileName,'r')
    OutputLog = open(dirname+'/radex.log','w')
    ExitCode = subprocess.call(RadexBinary,stdin = InputFile,stdout = OutputLog)
    InputFile.close()
    OutputLog.close()

    try:
        OutputFile=open(OutputFileName,'r')
    except IOError:
        print("BadRun")
        shutil.rmtree(dirname)        
        return(None)

    RadexDict = ReadRadex(OutputFile)
    shutil.rmtree(dirname)        
    return(RadexDict)
