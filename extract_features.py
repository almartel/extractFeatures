# -*- coding: utf-8 -*-
"""
Created on Tue Apr 01 10:47:40 2014

This script will load volumes, Load a lesion Segmentation (VOI), Visualize volumes
and then extract Dynamic, Morphology and Texture features from the VOI.

Arguments:
============
sys.argv[1] = input text file with one case per line in the following format:
StudyNumber    DicomExamNumber    MRN    chosen_lesions_id    StudyDate    SeriesID

@ author (C) Cristina Gallego, University of Toronto, 2013
----------------------------------------------------------------------
"""

import os, os.path
import sys
import string

from inputs_init import *
from display import *
from features_dynamic import *
from features_morphology import *
from features_texture import *
     

# Open filename list
file_ids = open(sys.argv[1],"r")
for fileline in file_ids:
    # Get the line: StudyNumber    DicomExamNumber    MRN    chosen_lesions_id    StudyDate    SeriesID    image_pos_pat    image_ori_pat
    fileline = fileline.split()
    StudyID = fileline[0] 
    AccessionN = fileline[1]
    Lesions_id = fileline[3]
    SeriesID = fileline[5] # corresponds to dynamic sequence;
        
    ###### Loading 
    print "Start by loading volumes..."
    # Get Root folder ( the directory of the script being run)
    path_rootFolder = 'Z:\Cristina\MassNonmass\mass'
    load = Inputs_init()
    [series_path, phases_series, lesionID_path, lesionseg_name] = load.readVolumes(path_rootFolder, StudyID, AccessionN, SeriesID, Lesions_id)
    print "Path to series location: %s" % series_path 
    print "List of pre and post contrast volume names: %s" % phases_series
    print "Path to lesion segmentation: %s" % lesionID_path
    
    print "\n Load Segmentation..."
    lesion3D = load.loadSegmentation(lesionID_path, lesionseg_name)
    print "Data Structure: %s" % lesion3D.GetClassName()
    print "Number of points: %d" % int(lesion3D.GetNumberOfPoints())
    print "Number of cells: %d" % int(lesion3D.GetNumberOfCells())
    
    print "\n Visualize volumes..."
    loadDisplay = Display()
    lesion3D_mesh = loadDisplay.addSegment(lesion3D, (0,1,0), interact=False)
    loadDisplay.visualize(load.DICOMImages, load.image_pos_pat, load.image_ori_pat, sub=True, postS=3, interact=True)

    ###### Extract Dynamic features
    print "\n Extract Dynamic features..."
    loadDynamic = Dynamic()
    dynamicfeatures_contour = loadDynamic.extractfeatures_contour(load.DICOMImages, load.image_pos_pat, load.image_ori_pat, series_path, phases_series, lesion3D)
    print dynamicfeatures_contour
    dynamicfeatures_inside = loadDynamic.extractfeatures_inside(load.DICOMImages, load.image_pos_pat, load.image_ori_pat, series_path, phases_series, lesion3D)
    print dynamicfeatures_inside
    
    ###### Extract Morphology features
    print "\n Extract Morphology features..."
    loadMorphology = Morphology()
    morphofeatures = loadMorphology.extractfeatures(load.DICOMImages, load.image_pos_pat, load.image_ori_pat, series_path, phases_series, lesion3D)
    print morphofeatures
    
    ###### Extract Texture features
    print "\n Extract Texture features..."
    loadTexture = Texture()
    texturefeatures = loadTexture.extractfeatures(load.DICOMImages, load.image_pos_pat, load.image_ori_pat, series_path, phases_series, lesion3D, loadMorphology.VOI_efect_diameter, loadMorphology.lesion_centroid)
    print texturefeatures
    
    
    
    
    
    
    
    
    
    
    
    
    
    
        
        
    