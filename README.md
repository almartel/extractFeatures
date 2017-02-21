Extract_Features Module: 

Description: Module that loads up data and segmented VOI to then extract dynamic morphological and texture features.

# Modules:
* load:  inputs_init.py
* display: display.py
* DICOM reader: processDicoms.py

# Requires (installed python packages):
* vtk
* pydicom or dicom
* pandas

# Usage: Main 
extract_features.py

This script will load all modules, MRI volumes, Load a lesion Segmentation (VOI), Visualize volumes
and then extract Dynamic, Morphology and Texture features from the VOI.

e.g 
###### Loading DICOMS
import os, os.path
import sys
import string

from inputs_init import *
from display import *
from features_dynamic import *
from features_morphology import *

load = Inputs_init()
# Get Data folder ( the directory of the DICOMS)
path_rootFolder = 'Z:\Breast\DICOMS'
Lesions_id = 12
StudyID = '0114'	
AccessionN = 6896014	
SeriesID = 600

[series_path, phases_series, lesionID_path, lesionseg_name] = load.readVolumes(path_rootFolder, StudyID, AccessionN, SeriesID, Lesions_id)
print "Path to series location: %s" % series_path 
print "List of pre and post contrast volume names: %s" % phases_series
print "Path to lesion segmentation: %s" % lesionID_path

######  Load Segmentation (specified in lesionID_path)
lesion3D = load.loadSegmentation(lesionID_path, lesionseg_name)
print "Data Structure: %s" % lesion3D.GetClassName()
print "Number of points: %d" % int(lesion3D.GetNumberOfPoints())
print "Number of cells: %d" % int(lesion3D.GetNumberOfCells())

###### Visualize volumes 
loadDisplay = Display()
lesion3D_mesh = loadDisplay.addSegment(lesion3D, (0,1,0))
loadDisplay.visualize(load.DICOMImages, load.image_pos_pat, load.image_ori_pat, sub=True, postS=4, interact=False)

###### Extract Dynamic features
print "\n Extract Dynamic features..."
loadDynamic = Dynamic()
dynamicfeatures_contour = loadDynamic.extractfeatures_contour(load.DICOMImages, load.image_pos_pat, load.image_ori_pat, series_path, phases_series, lesion3D)
print dynamicfeatures_contour
dynamicfeatures_inside = loadDynamic.extractfeatures_inside(load.DICOMImages, load.image_pos_pat, load.image_ori_pat, series_path, phases_series, lesion3D)
print dynamicfeatures_inside

###### Extract Morphology features
loadMorphology = Morphology()
morphofeatures = loadMorphology.extractfeatures(load.DICOMImages, load.image_pos_pat, load.image_ori_pat, series_path, phases_series, lesion3D)
print morphofeatures

# e.g To extract 2D texture:
from features_texture import *
loadTexture = Texture()
texturefeatures = loadTexture.extractfeatures(load.DICOMImages, load.image_pos_pat, load.image_ori_pat, series_path, phases_series, lesion3D, loadMorphology.VOI_efect_diameter, loadMorphology.lesion_centroid)

# e.g To extract 3D texture:
from features_3Dtexture import *
loadTexture3D = Texture3D()
texture3Dfeatures = loadTexture3D.extractfeatures(load.DICOMImages, load.image_pos_pat, load.image_ori_pat, series_path, phases_series, lesion3D, loadMorphology.VOI_efect_diameter, loadMorphology.lesion_centroid)
