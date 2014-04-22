#!/usr/bin/env python

import os
import os.path
import sys
import subprocess
from numpy import *
import dicom

from operator import itemgetter
import pandas as pd


'''
Compilation of functions

"usage:" 
import processDicoms
then on script use:
SUPPORT FUNCTIONS:
- processDicoms.FileCheck(argvs)
- processDicoms.is_number(argvs)
- processDicoms.get_immediate_subdirectories(argvs)
- processDicoms.get_only_linksindirectory(argvs)
- processDicoms.get_immediate_subdirectories(argvs)
- processDicoms.find(argvs)
- processDicoms.get_display_series(argvs)
- processDicoms.get_display(argvs)
- processDicoms.get_series(argvs)

PARSING DICOM SERIES:
- processDicoms.get_slices_at_all_locs(argvs)
- processDicoms.get_slices_for_volumes(argvs)

% Copyright (C) Cristina Gallego, University of Toronto, 2012, 2013-2014
----------------------------------------------------------------------
'''

# Function that checks whether file DICOMDIR.txt exits 
def FileCheck(filename):       
       try:
           fn=open(filename,"r") 
           fn.close()
           return True
       except IOError: 
           print "Error: DICOMDIR.txt doesn't exit, loading Series using vtkDICOMImageReader by setDirectoryName."
       return False

# Checks whether a number is float
def is_number(s):
    try:
        float(s)
        return True
    except ValueError:
        return False
        
# Gets immediate_subdirectories of folder mydir
def get_immediate_subdirectories(mydir):
    return [name for name in os.listdir(mydir) 
            if os.path.isdir(os.path.join(mydir, name))]

# Gets only_links in directory of folder mydir
def get_only_linksindirectory(mydir):
    return [name for name in os.listdir(mydir) 
            if os.path.islink(os.path.join(mydir, name))]

# Gets only_files in directory of folder mydir, excluding subdirectories folders
def get_only_filesindirectory(mydir):
     return [name for name in os.listdir(mydir) 
            if os.path.isfile(os.path.join(mydir, name))]

# Finds a substring within a string of chars
def find(strng, ch):
    index = 0
    while index < len(strng):
        if strng[index] == ch:
            return index
        index += 1                  
    return -1
    
# lists and displays the summary of DICOM Series available for display    
def get_display_series(abspath_SeriesID):
    arranged_folders = get_immediate_subdirectories(abspath_SeriesID);

    # Initialize series count
    s=0;
    print "Total number of series: %d" % len(arranged_folders)
    print " "
    print "%s     %s        %s     %s " % ('n', 'Series#', '#Images', 'SeriesDescription')
            
    # Iterate for each series in ExamID
    for arrangedF in arranged_folders:
        path_arrangedF_ID = abspath_SeriesID+os.sep+arrangedF
        #print path_SeriesID
        
        # Get total number of files
        listSeries_files = get_only_filesindirectory(path_arrangedF_ID)
                    
        # Use only the first slice one file and get DICOM DICTIONARY
        if(listSeries_files != []):
            path_filenameID = abspath_SeriesID+os.sep+arrangedF+os.sep+listSeries_files[0]
            
            NumberOfVolumes = 1 # default
            
            # iterate number of slices and get full volume (#slices, each slice loc)
            slices = []
            num_images=0;
            for filename in listSeries_files:
                num_images = num_images+1
                    
            # Print series info                        
            print "%d    %d        %s" % (s, num_images, arrangedF) 
            # increment series number
            s=s+1;    
        else:
            print "%d    %s        %d        %s" % (s, "NONE", 0, "NULL") 
            # increment series number
            s=s+1;    
        
    # Go back to rootfolder
    os.chdir(path_rootFolder)
    
    return arranged_folders

# Find all dicoms Series in sudyfolder abspath_SeriesID and Iterates
def get_display(abspath_SeriesID):
    arranged_folders = get_immediate_subdirectories(abspath_SeriesID);
    
    # Initialize series count
    s=0;
    print " "
    print "%s \t %s \t\t %s " % ('n', 'Series#', '#Images')
    
    # Find all dicoms Series in sudyfolder (will go into including subfolders)
    # Iterate
    for arrangedF in arranged_folders:
        #print "arrangedFolder: %s" % arrangedF
        path_arrangedF_ID = abspath_SeriesID+os.sep+arrangedF
        #print path_arrangedF_ID
        
        # Enter Studyfolder to process all Dicom Series of ExamID
        path_arrangedF_images = get_only_filesindirectory(path_arrangedF_ID);
        #print path_arrangedF_images

        print "%s \t %s \t\t %s " % (s, arrangedF, len(path_arrangedF_images))
        
        # increment series number
        s=s+1;    
        
    # Go back to rootfolder
    os.chdir(path_rootFolder)
    
    return arranged_folders            

# obtain all Series in subdirs in the StudyID directory 
def get_series(StudyID,img_folder):
    # obtain subdires in the StudyID directory, correspond to ExamsID 
    # check for one in subdirectory tree, (e.g Mass or NonMass) 
    global abspath_SeriesID
    
    path_studyID = img_folder+StudyID
    studyFolder = os.path.abspath(path_studyID)
    print studyFolder
    
    ExamsID = get_immediate_subdirectories(path_studyID);
    #print ExamsID
    c = 0
    if(len(ExamsID)>1):
        
        print "%s     %s    " % ('n', 'Series#')
        for iexam in ExamsID:
            print "%d    %s " % (c, str(iexam)) 
            c=c+1
                    
        choseSerie = raw_input('Enter n Series to load (0-n), or x to exit: ')
        if(choseSerie != 'x'):
            c = 0
            for iexam in ExamsID:
                if(int(choseSerie) == c):
                    eID = iexam
                c=c+1    
        else:
            return '', 0, '', studyFolder
            
        
        print "ExamID: %s" % eID
        path_ExamID = img_folder+StudyID+os.sep+eID
        abspath_ExamID = os.path.abspath(path_ExamID)
        print abspath_ExamID
        
        # Enter Studyfolder to process all Dicom Series of ExamID
        SeriesID = get_immediate_subdirectories(abspath_ExamID);
        #print SeriesID
        
        # Initialize series count
        s=0;
        print "Total number of series: %d" % len(SeriesID)
        print " "
        print "%s     %s        %s     %s " % ('n', 'Series#', '#Images', 'SeriesDescription')
                
        # Iterate for each series in ExamID
        for sID in SeriesID:
            if sID != 'DynPhases':
                path_SeriesID = img_folder+StudyID+os.sep+eID+os.sep+sID
                                
                abspath_SeriesID = os.path.abspath(path_SeriesID)
                
                # Get total number of files
                listSeries_files = get_only_filesindirectory(abspath_SeriesID)
                            
                # Use only the first slice one file and get DICOM DICTIONARY
                if(listSeries_files != []):
                    path_filenameID = img_folder+StudyID+os.sep+eID+os.sep+sID+os.sep+listSeries_files[0]
                
                    try:
                        dicomInfo = dicom.read_file(os.path.abspath(path_filenameID))
                    except ValueError:
                        dicomInfo = []
                        
                    # Get the main dataset (they are in fact two separate datasets in the DICOM standard). 
                    # That dicom dataset is now stored in the file_meta attribute of the dataset
                                                    
                    # Get structure of study (all files in directory consistent with studyID and patientID
                    studyTree = []
                    FileNames=listSeries_files;
                    if("PatientID" in dicomInfo):
                        PatientID = dicomInfo.PatientID#
                    else:    PatientID=''
                    if("SeriesNumber" in dicomInfo):
                        SeriesNumber = dicomInfo.SeriesNumber#
                    else:    SeriesNumber=''
                    if("SeriesDescription" in dicomInfo):
                        SeriesDescription = dicomInfo.SeriesDescription; #
                    else:    SeriesDescription=''
                    if("SliceLocation" in dicomInfo):
                        MinSliceLocation = dicomInfo.SliceLocation; # Infos to identify number of slices/volumes
                    else:    MinSliceLocation=''
                    if('ImageOrientationPatient' in dicomInfo):
                        ImageOrientationPatient = dicomInfo.ImageOrientationPatient; # Infos to identify number of slices/volumes
                    else:    ImageOrientationPatient=''
                    
                    NumberOfVolumes = 1 # default
                    
                # iterate number of slices and get full volume (#slices, each slice loc)
                slices = []
                num_images=0;
                for filename in listSeries_files:
                    num_images = num_images+1
                        
                # Print series info                        
                print "%d    %s        %d        %s" % (s, SeriesNumber, num_images, SeriesDescription) 
                # increment series number
                s=s+1;    
                
    else:    
        studyFolder = os.path.abspath(path_studyID)
        print studyFolder
        
        # Find all dicoms Series in sudyfolder (will go into including subfolders)
        # Iterate
        for eID in ExamsID:
            print "ExamID: %s" % eID
            path_ExamID = img_folder+StudyID+os.sep+eID
            abspath_ExamID = os.path.abspath(path_ExamID)
            print abspath_ExamID
            
            # Enter Studyfolder to process all Dicom Series of ExamID
            SeriesID = get_immediate_subdirectories(path_ExamID);
            #print SeriesID
            
            # Initialize series count
            s=0;
            print "Total number of series: %d" % len(SeriesID)
            print " "
            print "%s     %s        %s     %s " % ('n', 'Series#', '#Images', 'SeriesDescription')
                    
            # Iterate for each series in ExamID
            for sID in SeriesID:
    
                path_SeriesID = img_folder+StudyID+os.sep+eID+os.sep+sID
                
                abspath_SeriesID = os.path.abspath(path_SeriesID)
                #print abspath_SeriesID
                
                # Get total number of files
                listSeries_files = get_only_filesindirectory(abspath_SeriesID)
            
                # Use only the first slice one file and get DICOM DICTIONARY
                if(listSeries_files != []):
                    path_filenameID = img_folder+StudyID+os.sep+eID+os.sep+sID+os.sep+listSeries_files[0]
                    try:
                        dicomInfo = dicom.read_file(os.path.abspath(path_filenameID)) 
                    except ValueError:
                        dicomInfo = []
                                        
                    # Get the main dataset (they are in fact two separate datasets in the DICOM standard). 
                    # That dicom dataset is now stored in the file_meta attribute of the dataset
                        
                    # Get structure of study (all files in directory consistent with studyID and patientID
                    studyTree = []
                    FileNames=listSeries_files;
                    
                    if("PatientID" in dicomInfo):
                        PatientID = dicomInfo.PatientID#
                    else:    PatientID=''
                    if("SeriesNumber" in dicomInfo):
                        SeriesNumber = dicomInfo.SeriesNumber#
                    else:    SeriesNumber=''
                    if("SeriesNumber" in dicomInfo):
                        SeriesNumber = dicomInfo.SeriesNumber
                    else:    SeriesNumber=''
                    if("SeriesDescription" in dicomInfo):
                        SeriesDescription = dicomInfo.SeriesDescription;
                    else:    SeriesDescription=''
                    if("SliceLocation" in dicomInfo):
                        MinSliceLocation = dicomInfo.SliceLocation; # Infos to identify number of slices/volumes
                    else:    MinSliceLocation=''
                    if('ImageOrientationPatient' in dicomInfo):
                        ImageOrientationPatient = dicomInfo.ImageOrientationPatient; # Infos to identify number of slices/volumes
                    else:    ImageOrientationPatient=''
                    
                    NumberOfVolumes = 1 # default
                    
                    # iterate number of slices and get full volume (#slices, each slice loc)
                    slices = []
                    num_images=0;
                    for filename in listSeries_files:
                        num_images = num_images+1
                            
                    # Print series info                        
                    print "%d    %s        %d        %s" % (s, SeriesNumber, num_images, SeriesDescription) 
                    # increment series number
                    s=s+1;    
                else:
                    print "%d    %s        %d        %s" % (s, "NONE", 0, "NULL") 
                    # increment series number
                    s=s+1;    
                
        # Go back to rootfolder
        #os.chdir(path_rootFolder)    
                
    return abspath_ExamID, eID, SeriesID, studyFolder, dicomInfo

def get_slices_at_all_locs(img_folder, abspath_SeriesID, len_listSeries_files, StudyID, eID, SeriesID, choseSerie, listSeries_files, abspath_ExamID):
    # enter folder of Series selection /examID/S_XXX 
    print "\n------------------------------------------------------------------"   
    print "\n entering folder abspath_SeriesID"
                    
    # Obtain all datasets in current series
    # iterate number of slices and get full volume (#slices, each slice loc)
    slices = []
    FileNms_slices =  []
                
    for n in range(len_listSeries_files):
        # Use all DICOM slices on series
        ''' EXTRACT DICOM SLICE LOCATION '''
        absp_fsID = img_folder+StudyID+'/'+eID+'/'+SeriesID[int(choseSerie)]+'/'+listSeries_files[n]
                            
        dInfo = dicom.read_file(absp_fsID)
        slices.append(dInfo.SliceLocation)
        FileNms_slices.append(listSeries_files[n])
        FileNms_slices.append(dInfo.SliceLocation)
        
    print '''\nPROCESS STACKS BY SLICE LOCATIONS '''
    FileNms_slices_stack = reshape(FileNms_slices, [len_listSeries_files,2])
    
    FileNms_slices_sorted = sorted(FileNms_slices_stack, key=itemgetter(1))
    #print FileNms_slices_sorted
    #time.sleep(5) #will sleep for 5 seconds
    
    FileNms_slices_sorted_stack = reshape(FileNms_slices_sorted, [len_listSeries_files,2])
    #print FileNms_slices_sorted_stack
    #time.sleep(5) #will sleep for 5 seconds
    
    current_slice = FileNms_slices_sorted_stack[0,1]
    print current_slice
    stack_byLocation = []
    name_byLocation = []
    scount = 0
        
    for sliceS in FileNms_slices_sorted_stack:
        # Get the num_series = 5 name_byLocations for a given Location
        if( current_slice == sliceS[1]):
            print "Location: %s" % sliceS[1]
            #print "Slice_loc: %s" % sliceS[1]
            stack_byLocation.append(sliceS[1])
            name_byLocation.append(sliceS[0])
            scount = scount+1
    
        # Finish getting all series for a given Location
        
    # finalize the last location move
    print '''\nFINISH get_slices_at_all_locs \n'''   
    print "Number of pre+post-contrast scans: %d" % scount
    
    return scount

def get_slices_for_volumes(scount, choseSerie, img_folder, abspath_SeriesID, len_listSeries_files, StudyID, eID, SeriesID, listSeries_files, abspath_ExamID ):
    global num_locs_per_Vol
    print "\n entering folder abspath_SeriesID"
                
    # Obtain all datasets in current series
    # iterate number of slices and get full volume (#slices, each slice loc)
    slices = []
    FileNms_slices =  []
        
    for n in range(len_listSeries_files):
        # Use all DICOM slices on series
        ''' EXTRACT DICOM SLICE LOCATION '''
        absp_fsID = img_folder+StudyID+'/'+eID+'/'+SeriesID[int(choseSerie)]+'/'+listSeries_files[n]
            
        dInfo = dicom.read_file(listSeries_files[n])
        slices.append(dInfo.SliceLocation)
        FileNms_slices.append(listSeries_files[n])
        FileNms_slices.append(dInfo.SliceLocation)
    
    print '''\nPROCESS STACKS BY SLICE LOCATIONS '''
    print "Total number of locations found: %d" % scount
    FileNms_slices_stack = reshape(FileNms_slices, [len_listSeries_files,2])
    #print FileNms_slices_stack
    FileNms_slices_sorted = sorted(FileNms_slices_stack, key=itemgetter(0))
    FileNms_slices_sorted_stack = reshape(FileNms_slices_sorted, [len_listSeries_files,2])
    print   FileNms_slices_sorted_stack

    stack_byLocation = []
    name_byLocation = []
    
    # Extracting the number of slices per bilateral 3D volumes based on num_series for 280/5 = 56
    num_locs_per_Vol = int(float(len_listSeries_files)/float(scount))
    scount = int(scount)
    
    stack_byBilatVol = []
        
    k=0
    # Get the folder names based on num_series
    for numlocs in range(scount):
        if ( numlocs == 0):    
            stack_byBilatVol.append('pre-Contrast')
        else:# Now link slices at location to folder
            stack_byBilatVol.append('post_Contrast-'+str(k))
        k=k+1
            
    # Initialized BilatVol
    print stack_byBilatVol
    slice_i = 0
    
    for Vol_i in range(scount):    
        print'''-----\t NOW HAVE ALL SLICES IN A bilateral vol '''
        # Get the new Location folder
        print "bit vol %d" % Vol_i
        current_vol = stack_byBilatVol[Vol_i]
                                
        # Makedir of current loc
        # if current loc folder doesn't exist create it        
        if not os.path.exists(str(current_vol)):
            os.makedirs(str(current_vol))
        
        # Get inside location directory
        os.chdir(str(current_vol))
        print os.getcwd()
        
        #FileNms_slices_sorted_stack[slice_i,0]
        # Now link slices at location to folder
        filename = str(FileNms_slices_sorted_stack[slice_i,0])
        filename = filename[0:-10]
        file_ending = '.MR.dcm'
        
        # Save the file list to read as series later
        filename_series = 'DIRCONTENTS.txt'
        file_series = open(filename_series, 'a')
        
        for j in range(num_locs_per_Vol):
            link_to = '../'+FileNms_slices_sorted_stack[slice_i,0]
            if ( j < 9):
                name4link_to = filename+'00'+str(j+1)+file_ending
            else:
                name4link_to = filename+'0'+str(j+1)+file_ending
            print "linking file: %s to: %s" % (link_to, name4link_to)
            ln_subp = subprocess.Popen(['ln', link_to, name4link_to], stdout=subprocess.PIPE)
            ln_subp.wait()
            file_series.write(str(name4link_to)+'\n')
            slice_i = slice_i+1
            
        file_series.close()
            
        # Get back inside the  Series directory
        os.chdir(abspath_SeriesID)
                        
    print "\n------------------------------------------------------------------"                        
    print '''FINISH get_slices_for_volumes \n'''    
    print "------\tNumber of Dynamic Volumes (series time points: including pre-contrast): %d" % scount 
    print "------\tNumber of Locations per Bilateral Volume: %d " % num_locs_per_Vol
    os.chdir(abspath_ExamID)
    
    return stack_byBilatVol

def ReadDicomfiles(abspath_PhaseID):
    slices = []
    FileNms_slices =  []
    
    listSeries_files = sorted(list(get_only_filesindirectory(str(abspath_PhaseID))))
    len_listSeries_files = len(listSeries_files)
                
    for n in range(len_listSeries_files):
        # Use all DICOM slices on series
        ''' EXTRACT DICOM SLICE LOCATION '''
        absp_fsID = str(abspath_PhaseID)+os.sep+listSeries_files[n]
                
        if(absp_fsID!=str(abspath_PhaseID)+os.sep+'DIRCONTENTS.txt'):
            dInfo = dicom.read_file(absp_fsID)
            slices.append(dInfo.SliceLocation)
            FileNms_slices.append(listSeries_files[n])
        else:
            cmd = 'rm '+str(abspath_PhaseID)+os.sep+'DIRCONTENTS.txt'
            print "cmd -> " + cmd
            p1 = subprocess.Popen(cmd, shell=True, stdin=subprocess.PIPE)
            p1.wait()
            len_listSeries_files = len_listSeries_files-1

    print "Total images in series: %d " % len_listSeries_files

    '''\nPROCESS STACKS BY SLICE LOCATIONS '''
    FileNms_slices_stack = pd.DataFrame({'slices': FileNms_slices,
                                         'location': slices})
            
    return len_listSeries_files, FileNms_slices_stack
    
## END

