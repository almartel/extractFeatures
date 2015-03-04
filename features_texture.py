# -*- coding: utf-8 -*-
"""
Created on Wed Apr 02 09:40:11 2014
@ author (C) Cristina Gallego, University of Toronto
----------------------------------------------------------------------
"""


from numpy import *
import vtk
from vtk.util.numpy_support import vtk_to_numpy
from pandas import DataFrame
from pylab import *
import matplotlib.pyplot as plt
from scipy import stats
from skimage.feature import greycomatrix, greycoprops

from display import *

#!/usr/bin/env python
class Texture(object):
    """
    USAGE:
    =============
    texture_features = Texture()
    """
    def __init__(self): 
        """ initialize Texture """
        self.texture_features = []
    
        
    def __call__(self):       
        """ Turn Class into a callable object """
        Texture()
        
    def histeq(self, im, nbr_bins=256):
        #get image histogram  
        imhist,bins = histogram(im.flatten(),nbr_bins,normed=True)
        cdf = imhist.cumsum() #cumulative distribution function
        cdf = 255 * cdf / cdf[-1] #normalize
        
        #use linear interpolation of cdf to find new pixel values
        im2 = interp(im.flatten(),bins[:-1],cdf)
        
        return im2.reshape(im.shape), cdf 

   
    def extractfeatures(self, DICOMImages, image_pos_pat, image_ori_pat, series_path, phases_series, VOI_mesh, VOI_efect_diameter, lesion_centroid):
        ################################### 
        # Haralick et al. defined 10 grey-level co-occurrence matrix (GLCM) enhancement features (energy, maximum probability, contrast, homogeneity, entropy, correlation, sum average, sum variance, difference average and difference variance) to describe texture
        # For both mass and non-mass lesions
        averg_texture_contrast=array([0,0,0,0]).reshape(1,4)    
        averg_texture_homogeneity=array([0,0,0,0]).reshape(1,4)
        averg_texture_dissimilarity=array([0,0,0,0]).reshape(1,4)
        averg_texture_correlation=array([0,0,0,0]).reshape(1,4)
        averg_texture_ASM=array([0,0,0,0]).reshape(1,4)
        averg_texture_energy=array([0,0,0,0]).reshape(1,4)
        
        # N is the number of distinct gray levels in the histogram equalized image;
        for i in range(1,len(DICOMImages)):
            # obtain vols of interest
            subtractedImage = Display().subImage(DICOMImages, i)
            [sub_pre_transformed_image, transform_cube] = Display().dicomTransform(subtractedImage, image_pos_pat, image_ori_pat) 
            
            print "\nLoading VOI_mesh MASK for Series %s... " % str(i)
            VOIPnt = [0,0,0]; pixVals_margin = [];  pixVals = []
            
            #################### HERE GET INTERNAL PIXELS IT AND MASK IT OUT
            VOI_scalars = sub_pre_transformed_image.GetPointData().GetScalars()
            np_VOI_imagedata = vtk_to_numpy(VOI_scalars)     
            
            dims = sub_pre_transformed_image.GetDimensions()
            spacing = sub_pre_transformed_image.GetSpacing()
            np_VOI_imagedata = np_VOI_imagedata.reshape(dims[2], dims[1], dims[0]) 
            np_VOI_imagedata = array(np_VOI_imagedata.transpose(2,1,0))
            
            # Prepare lesion localization and PATCH size for Texture analysis    
            # lesion_centroid        
            lesion_radius = VOI_efect_diameter/(spacing[0])
            print "VOI_efect_diameter %s... " % str(lesion_radius)
            
            lesion_location = lesion_centroid
            print "lesion_location %s... " % str(lesion_location)
            lesion_patches = []
            
            ######### translate lesion location to ijk coordinates
            # sub_pre_transformed_image.FindPoint(lesion_location
            pixId = sub_pre_transformed_image.FindPoint(lesion_location[0], lesion_location[1], lesion_location[2])
            sub_pt = [0,0,0]            
            sub_pre_transformed_image.GetPoint(pixId, sub_pt)
            
            ijk = [0,0,0]
            pco = [0,0,0]
                    
            inorout = sub_pre_transformed_image.ComputeStructuredCoordinates( sub_pt, ijk, pco)
            print "coresponding ijk_vtk indices:"
            print ijk
            ijk_vtk = ijk            

            # compute texture cdf        
            eq_numpy_pre_VOI_imagedata, cdf = self.histeq(np_VOI_imagedata[:,:,int(ijk_vtk[2])])    
            plt.figure()
            plt.subplot(231)
            n, bins, patches = plt.hist(array(np_VOI_imagedata[:,:,int(ijk_vtk[2])].flatten()), 50, normed=1, facecolor='green', alpha=0.75)
            
            plt.subplot(233)
            n, bins, patches = plt.hist(array(eq_numpy_pre_VOI_imagedata.flatten()), 50, normed=1, facecolor='green', alpha=0.75)
             
            plt.subplot(234)
            plt.imshow(np_VOI_imagedata[:,:,int(ijk_vtk[2])])
            plt.gray()
            
            plt.subplot(236)
            plt.imshow(eq_numpy_pre_VOI_imagedata)
            plt.gray()
            
            # FInally display
            plt.show()            

            #****************************"
            print "Normalizing data..."
            np_VOI_imagedataflat = np_VOI_imagedata.flatten().astype(float)
            np_VOI_imagedata_num = (np_VOI_imagedataflat-min(np_VOI_imagedataflat))
            np_VOI_imagedata_den = (max(np_VOI_imagedataflat)-min(np_VOI_imagedataflat))
            np_VOI_imagedata_flatten = 255*(np_VOI_imagedata_num/np_VOI_imagedata_den)
            np_VOI_imagedata = np_VOI_imagedata_flatten.reshape(dims[0], dims[1], dims[2])
            print np_VOI_imagedata
            
            # Do histogram cast
            eq_numpy_pre_VOI_imagedata, cdf = self.histeq(np_VOI_imagedata[:,:,int(ijk_vtk[2])])
            
            
            # Get the shape of the numpy image data to confirm dimensions
            eq_numpy_pre_VOI_imagedata = eq_numpy_pre_VOI_imagedata.reshape(dims[0], dims[1], 1)
            
            # Perform texture classification using grey level co-occurrence matrices (GLCMs).
            # A GLCM is a histogram of co-occurring greyscale values at a given offset over an image.
            # compute some GLCM properties each patch
            # p(i,j) is the (i,j)th entry in a normalized spatial gray-level dependence matrix; 
            lesion_patches = []
            lesion_patches = np_VOI_imagedata[
                    int(ijk_vtk[0] - lesion_radius):int(ijk_vtk[0] + lesion_radius),
                    int(ijk_vtk[1] - lesion_radius):int(ijk_vtk[1] + lesion_radius),
                    int(ijk_vtk[2])]
                                    
            print '\n Lesion_patches:'
            print lesion_patches
    
            #skimage.feature.greycomatrix(image, distances, angles, levels=256, symmetric=False, normed=False)
            glcm = greycomatrix(lesion_patches, [3], [0, 45*pi/180, 90*pi/180, 135*pi/180], 256, symmetric=True, normed=True)
            texture_contrast = greycoprops(glcm, 'contrast')
            texture_homogeneity = greycoprops(glcm, 'homogeneity')
            texture_dissimilarity = greycoprops(glcm, 'dissimilarity')
            texture_correlation = greycoprops(glcm, 'correlation')
            texture_ASM = greycoprops(glcm, 'ASM')
            texture_energy = greycoprops(glcm, 'energy')
                        
            """TEXTURE FEATURES accumulated"""
            averg_texture_contrast = averg_texture_contrast+texture_contrast
            averg_texture_homogeneity = averg_texture_homogeneity+ texture_homogeneity
            averg_texture_dissimilarity = averg_texture_dissimilarity+texture_dissimilarity
            averg_texture_correlation = averg_texture_correlation+texture_correlation
            averg_texture_ASM = averg_texture_ASM+texture_ASM
            averg_texture_energy = averg_texture_energy+texture_energy
                        
            # display the image patches
            plt.subplot(3, len(lesion_patches), len(lesion_patches) * 1 + i + 1)
            plt.imshow(lesion_patches, cmap=plt.cm.gray, interpolation='nearest',
                   vmin=0, vmax=255)
            plt.xlabel('lesion_patches %d' % (i + 1))
            
            # display original image with locations of patches
            plt.subplot(3, 2, 1)
            plt.imshow(np_VOI_imagedata[:,:,int(ijk_vtk[2])] , cmap=plt.cm.gray, interpolation='nearest',
                   vmin=0, vmax=255)
            # Plot
            # create the figure
            plt.plot(ijk_vtk[0] - lesion_radius/2, ijk_vtk[1] - lesion_radius/2, 'gs')
            plt.xlabel('Original Image')
            plt.xticks([])
            plt.yticks([])
            plt.axis('image')
            
            # FInally display
            plt.show()
            
            #plt.close()     
            
        # writing to file from row_lesionID Drow_PathRepID
        print "\n Average texture features for each orientation"
        contrast = averg_texture_contrast/4
        print contrast
        [self.contrast_zero, self.contrast_quarterRad, self.contrast_halfRad, self.contrast_halfRad] =  contrast[0,0], contrast[0,1], contrast[0,2], contrast[0,3]
        
        homogeneity = averg_texture_homogeneity/4
        print homogeneity                
        [self.homogeneity_zero, self.homogeneity_quarterRad, self.homogeneity_halfRad, self.homogeneity_threeQuaRad] = homogeneity[0,0], homogeneity[0,1], homogeneity[0,2], homogeneity[0,3]

        dissimilarity = averg_texture_dissimilarity/4
        print dissimilarity                
        [self.dissimilarity_zero, self.dissimilarity_quarterRad, self.dissimilarity_halfRad, self.dissimilarity_threeQuaRad] = dissimilarity[0,0], dissimilarity[0,1], dissimilarity[0,2], dissimilarity[0,3]
        
        correlation = averg_texture_correlation/4
        print correlation  
        [self.correlation_zero, self.correlation_quarterRad, self.correlation_halfRad, self.correlation_threeQuaRad] = correlation[0,0], correlation[0,1], correlation[0,2], correlation[0,3]
        
        ASM = averg_texture_ASM/4
        print ASM  
        [self.ASM_zero, self.ASM_quarterRad, self.ASM_halfRad, self.ASM_threeQuaRad] = ASM[0,0], ASM[0,1], ASM[0,2], ASM[0,3]
        
        energy = averg_texture_energy/4
        print energy  
        [self.energy_zero, self.energy_quarterRad, self.energy_halfRad, self.energy_threeQuaRad] = energy[0,0], energy[0,1], energy[0,2], energy[0,3]

        
      
        ##################################################
        # orgamize into dataframe
        self.texture_features = DataFrame( data=array([[ self.contrast_zero, self.contrast_quarterRad, self.contrast_halfRad, self.contrast_halfRad, 
                                                   self.homogeneity_zero, self.homogeneity_quarterRad, self.homogeneity_halfRad, self.homogeneity_threeQuaRad,
                                                   self.dissimilarity_zero, self.dissimilarity_quarterRad, self.dissimilarity_halfRad, self.dissimilarity_threeQuaRad,
                                                   self.correlation_zero, self.correlation_quarterRad, self.correlation_halfRad, self.correlation_threeQuaRad,
                                                   self.ASM_zero, self.ASM_quarterRad, self.ASM_halfRad, self.ASM_threeQuaRad,
                                                   self.energy_zero, self.energy_quarterRad, self.energy_halfRad, self.energy_threeQuaRad ]]), 
        columns=['texture_contrast_zero', 'texture_contrast_quarterRad', 'texture_contrast_halfRad', 'texture_contrast_threeQuaRad', 'texture_homogeneity_zero', 'texture_homogeneity_quarterRad', 'texture_homogeneity_halfRad', 'texture_homogeneity_threeQuaRad', 'texture_dissimilarity_zero', 'texture_dissimilarity_quarterRad', 'texture_dissimilarity_halfRad', 'texture_dissimilarity_threeQuaRad', 'texture_correlation_zero', 'texture_correlation_quarterRad', 'texture_correlation_halfRad', 'texture_correlation_threeQuaRad', 'texture_ASM_zero', 'texture_ASM_quarterRad', 'texture_ASM_halfRad', 'texture_ASM_threeQuaRad', 'texture_energy_zero', 'texture_energy_quarterRad', 'texture_energy_halfRad', 'texture_energy_threeQuaRad'])


        return self.texture_features