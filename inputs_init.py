# -*- coding: utf-8 -*-
"""
 USAGE:
=============
from inputs_init import *
input = Inputs_init()
input.readVolumes(StudyID, DicomExamNumber, SeriesID, Lesions_id)
input.loadSegmentation(lesionID_path)

Class Methods:
=============
readVolumes(StudyID, DicomExamNumber, SeriesID, Lesions_id)
loadSegmentation(lesionID_path)

Class Instance Attributes:
===============
'slice_thickn': 3.0, 
'image_ori_pat': ['-0', '1', '0', '-0', '-0', '-1'], 
'spacing': (0.44920000433921814, 0.44920000433921814, 3.0), 
'dims': (512, 512, 96), 
'image_pos_pat': ['145.059', '-167.095', '69.3364'], 
'DICOMImages': list[(vtkImageData)
'lesion3D_mesh': (vtkPolyData)
    
Created on Wed Apr 02 13:42:50 2014

@ author (C) Cristina Gallego, University of Toronto, 2013
----------------------------------------------------------------------
"""

import os, os.path
import sys
import string
import dicom
import vtk
import pandas as pd
import subprocess
from glob import glob

#!/usr/bin/env python
class Inputs_init:
    """
    USAGE:
    =============
    input = Inputs_init()
    input.readVolumes(StudyID, DicomExamNumber, SeriesID, Lesions_id)
    
    """
    def __init__(self):
        self.DICOMImages = []          
            
    def __call__(self):       
        """ Turn Class into a callable object """
        Inputs_init()
    
    # Gets only_files in directory of folder mydir, excluding subdirectories folders
    def get_only_filesindirectory(self, mydir):
        return [name for name in os.listdir(mydir) 
            if os.path.isfile(os.path.join(mydir, name))]
       
    def ReadDicomfiles(self, abspath_PhaseID):
        slices = []
        FileNms_slices =  []
        
        listSeries_files = sorted(list(self.get_only_filesindirectory(str(abspath_PhaseID))))
        len_listSeries_files = len(listSeries_files)
                    
        for n in range(len_listSeries_files):
            # Use all DICOM slices on series
            ''' EXTRACT DICOM SLICE LOCATION '''
            absp_fsID = str(abspath_PhaseID)+os.sep+listSeries_files[n]
            dInfo = dicom.read_file(absp_fsID)
            slices.append(dInfo.SliceLocation)
            FileNms_slices.append(listSeries_files[n])
    
        print "Total images in series: %d " % len_listSeries_files
    
        '''\nPROCESS STACKS BY SLICE LOCATIONS '''
        slices_stack = pd.DataFrame({'slices': FileNms_slices,
                                             'location': slices})
        # sort
        FileNms_slices_stack = slices_stack.sort(['location'], ascending=1)
            
        return len_listSeries_files, FileNms_slices_stack
    
        
    def readVolumes(self, path_rootFolder, StudyID, DicomExamNumber, SeriesID, Lesions_id):
        """
        ARGUMENTS:
        =============
        StudyID (str)           Study CAD patient number      
        DicomExamNumber (str)   Imaging Study Number    
        SeriesID  (str)         Dynamic Series number (e.g S600, S3)     
        Lesions_id (str)        Older database LesionID identifier
        
        OUTPUTS:
        =============
        series_path (str)       Path to series location      
        phases_series (list)      list of pre and post contrast volume names      
        lesionID_path  (str)      path to lesion segmentation      
        DICOMImages list(vtkImageData)      List of vtkImageData objects corresponding to dynamic series
        
        """
        series_path = path_rootFolder+os.sep+str(StudyID)+os.sep+str(DicomExamNumber)
        print series_path
        # test
        ###############################################################
        os.chdir(series_path)   
        if os.path.exists('DynPhases'):
            print '''DynPhases'''
            
            # Get series to load
            series_path = path_rootFolder+os.sep+str(StudyID)+os.sep+str(DicomExamNumber)
            lesionID_path = series_path+os.sep+'DynPhases'+os.sep+'VOIlesions_id'+str(Lesions_id)
            
            phases_series=[]
            testSID = str(SeriesID)
            if 'S' in str(testSID):
                #print testSID[1:]
                chosen_phase = int(testSID[1:])
            else:
                chosen_phase = int(testSID)
            
            if(testSID[0] == 'S'):
                phases_series.append('S'+str(chosen_phase))
                                
                for chSer in [chosen_phase+1, chosen_phase+2, chosen_phase+3, chosen_phase+4]:
                    phases_series.append( 'S'+str(chSer) )    
            else:
                phases_series.append(str(chosen_phase))
                                
                for chSer in [chosen_phase+1, chosen_phase+2, chosen_phase+3, chosen_phase+4]:
                    phases_series.append( str(chSer) )

        if not os.path.exists('DynPhases'):
            print '''SeriesPhases'''
            
            # Get series to load
            series_path = path_rootFolder+os.sep+str(StudyID)+os.sep+str(DicomExamNumber)+os.sep+str(SeriesID)
            lesionID_path = series_path+os.sep+'VOIlesions_id'+str(Lesions_id)           
                                                        
            # process all Volumes when in stacks of Dyn Volumes
            if os.path.exists(series_path+os.sep+'pre-Contrast'):
                phases_series = []
                phases_series.append('pre-Contrast')
                        
                #"Arranging series scans"
                for i in range(1,5):
                    phases_series.append('post_Contrast-'+str(i))
        
        # Get total number of files and some DICOM tags needed fro DICOM coords
        pre_abspath_PhaseID = series_path+os.sep+phases_series[0]
        [len_listSeries_files, FileNms_slices_sorted_stack] = self.ReadDicomfiles(pre_abspath_PhaseID)
        mostleft_slice = FileNms_slices_sorted_stack.slices[0]
        
        # Get dicom header, retrieve: image_pos_pat and image_ori_pat
        dicomInfo_series = dicom.read_file(pre_abspath_PhaseID+os.sep+str(mostleft_slice)) 
        self.slice_thickn = float(dicomInfo_series[0x0018,0x0050].value)  
        
        # Image Position (0020,0032): specifies the x, y, and z coordinates of the upper left hand corner of the image. This tag specifies the coordinates 
        # of the the first voxel transmitted.
        self.image_pos_pat = list(dicomInfo_series[0x0020,0x0032].value)
        # Image Orientation (0020,0037): specifies the direction cosines 
        # of the first row and the first column with respect to the patient. 
        # The direction of the axes are defined by the patients orientation 
        # to ensure LPS system ( x-axis increasing to the left hand side of the patient, 
        # y-axis increasing to the posterior side of the patient and z-axis increasing toward
        # the head of the patient )
        self.image_ori_pat = list(dicomInfo_series[0x0020,0x0037].value)
        
        for i in range(0,len(phases_series)):
            abspath_PhaseID = series_path+os.sep+phases_series[i] 
            print abspath_PhaseID
                                     
            # Done reading volume filenames now load them to DICOMimageReader                                         
            os.chdir(abspath_PhaseID)               
            dicomReader  = vtk.vtkDICOMImageReader()
            dicomReader.SetDirectoryName( abspath_PhaseID )
            dicomReader.Update()
            
            im = vtk.vtkImageData()
            im =  dicomReader.GetOutput()
            self.dims = im.GetDimensions()
            self.spacing = im.GetSpacing()
            
            print "VTK Dimensions im.GetDimensions(): %d %d %d" % self.dims
            print "VTK Spacing im.GetSpacing(): %f %f %f\n" % self.spacing
            
            # Append to objects image            
            self.DICOMImages.append( im )
            
        
        return(series_path, phases_series, lesionID_path)        
        
        
    def readT2(self, path_T2Series):
        """
        ARGUMENTS:
        =============
        StudyID (str)           Study CAD patient number      
        DicomExamNumber (str)   Imaging Study Number    
        SeriesID  (str)         T2 Series number (e.g S600, S3)     
        Lesions_id (str)        Older database LesionID identifier
        
        OUTPUTS:
        =============
        DICOMImages list(vtkImageData)      List of vtkImageData objects corresponding to dynamic series
        
        """
        # Get total number of files and some DICOM tags needed fro DICOM coords
        [len_listSeries_files, FileNms_slices_sorted_stack] = self.ReadDicomfiles(path_T2Series)
        mostleft_slice = FileNms_slices_sorted_stack.slices[0]
        
        # Get dicom header, retrieve: image_pos_pat and image_ori_pat
        dicomInfoT2 = dicom.read_file(path_T2Series+os.sep+str(mostleft_slice)) 
        self.slice_thicknT2 = float(dicomInfoT2[0x0018,0x0050].value)  
        
        # Image Position (0020,0032): specifies the x, y, and z coordinates of the upper left hand corner of the image. This tag specifies the coordinates 
        self.T2image_pos_pat = list(dicomInfoT2[0x0020,0x0032].value)
        self.T2image_ori_pat = list(dicomInfoT2[0x0020,0x0037].value)
        self.T2fatsat = dicomInfoT2[0x0019,0x10a4].value
       
        os.chdir(path_T2Series)               
        dicomReader  = vtk.vtkDICOMImageReader()
        dicomReader.SetDirectoryName( path_T2Series )
        dicomReader.Update()
        
        im = vtk.vtkImageData()
        im =  dicomReader.GetOutput()
        self.T2dims = im.GetDimensions()
        self.T2spacing = im.GetSpacing()
        
        print "VTK Dimensions im.GetDimensions(): %d %d %d" % self.T2dims
        print "VTK Spacing im.GetSpacing(): %f %f %f\n" % self.T2spacing
        
        # to objects image            
        self.T2Images = im 
        
        return self.T2Images
        
        
    def loadSegmentation(self, lesionID_path, lesionname):
        """
        ARGUMENTS:
        =============
        lesionID_path: (str)        path to lesion segmentation

        OUTPUT:
        =============
        lesion3D (vtkPolyData)      3D lesion segmentation as a vtkPolyData object       
        """ 
        # need to locate and read Lesion Seg
        if(lesionname):
            VOIlesion=' '
        else:            
            VOIlesion = os.sep+'VOIlesion_selected.vtk'
            
        lesion3D_reader = vtk.vtkPolyDataReader()
        lesion3D_reader.SetFileName( lesionID_path+VOIlesion )
        lesion3D_reader.Update()
        
        # Extract the polydata from the PolyDataReader        
        self.lesion3D_mesh = lesion3D_reader.GetOutput()
      
        return self.lesion3D_mesh
        
       
    
    
    
    
    
    
    
    
    
    
    
    
        
        
    