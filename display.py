# -*- coding: utf-8 -*-
"""
Create visualization with standard vtk actors, renders, windowsn, interactors

USAGE:
=============
from display import *
loadDisplay = Display()  
loadDisplay.dicomTransform(image, image_pos_pat, image_ori_pat)
loadDisplay.addSegment(lesion3D)
loadDisplay.subImage(Images2Sub, timep)                  
loadDisplay.visualize(images, image_pos_pat, image_ori_pat, sub, postS, interact)

Class Methods:
=============
dicomTransform(image, image_pos_pat, image_ori_pat)
addSegment(lesion3D)
subImage(Images2Sub, timep)                  
visualize(images, image_pos_pat, image_ori_pat, sub, postS, interact)

Class Instance Attributes:
===============
'origin': (-167.0, -69.0, -145.0)
'spacing': (0.44920000433921814, 0.44920000433921814, 3.0)
'dims': (512, 512, 96), 

VTK Instance objects:
=============
'xImagePlaneWidget': (vtkImagePlaneWidget)
'yImagePlaneWidget': (vtkImagePlaneWidget)
'zImagePlaneWidget': (vtkImagePlaneWidget)
'picker': (vtkCellPicker)
'iren1': (vtkWin32RenderWindowInteractor)
'camera': (vtkOpenGLCamera)
'mapper_mesh': (vtkPainterPolyDataMapper)
'actor_mesh': (vtkOpenGLActor)
'renWin1': (vtkWin32OpenGLRenderWindow)
'renderer1': (vtkOpenGLRenderer)


Created on Tue Apr 01 10:18:34 2014
@ author (C) Cristina Gallego, University of Toronto, 2013
----------------------------------------------------------------------
"""

import os, os.path
import sys
import string
from sys import argv, stderr, exit
import vtk
import numpy as np

#!/usr/bin/env python
class Display(object):
    """
    USAGE:
    =============
    loadDisplay = Display()
    """
    def __init__(self): 
        """ initialize visualization with standard vtk actors, renders, windowsn, interactors """           
        # use cell picker for interacting with the image orthogonal views.
        self.picker = vtk.vtkCellPicker()
        self.picker.SetTolerance(0.005) 
        
        # Create 3 orthogonal view using the ImagePlaneWidget
        self.xImagePlaneWidget = vtk.vtkImagePlaneWidget()
        self.yImagePlaneWidget = vtk.vtkImagePlaneWidget()
        self.zImagePlaneWidget = vtk.vtkImagePlaneWidget()
        
        #  The 3 image plane widgets
        self.xImagePlaneWidget.DisplayTextOn();
        self.xImagePlaneWidget.SetPicker(self.picker);
        self.xImagePlaneWidget.RestrictPlaneToVolumeOn();
        self.xImagePlaneWidget.SetKeyPressActivationValue('x');
        self.xImagePlaneWidget.GetPlaneProperty().SetColor(1, 0, 0);
        self.xImagePlaneWidget.SetResliceInterpolateToNearestNeighbour();
        
        self.yImagePlaneWidget.DisplayTextOn();
        self.yImagePlaneWidget.SetPicker(self.picker);
        self.yImagePlaneWidget.RestrictPlaneToVolumeOn();
        self.yImagePlaneWidget.SetKeyPressActivationValue('y');
        self.yImagePlaneWidget.GetPlaneProperty().SetColor(0, 1, 0);
        self.yImagePlaneWidget.SetLookupTable(self.xImagePlaneWidget.GetLookupTable());
        
        self.zImagePlaneWidget.DisplayTextOn();
        self.zImagePlaneWidget.SetPicker(self.picker);
        self.zImagePlaneWidget.SetKeyPressActivationValue('z');
        self.zImagePlaneWidget.GetPlaneProperty().SetColor(0, 0, 1);
        self.zImagePlaneWidget.SetLookupTable(self.xImagePlaneWidget.GetLookupTable());
        self.zImagePlaneWidget.SetRightButtonAutoModifier(1);
        
        # Create a renderer, render window, and render window interactor to
        # display the results.
        self.renderer1 = vtk.vtkRenderer()
        self.renWin1 = vtk.vtkRenderWindow()
        self.iren1 = vtk.vtkRenderWindowInteractor()
        
        self.renWin1.SetSize(1000, 800);
        self.renWin1.AddRenderer(self.renderer1);
        self.iren1.SetRenderWindow(self.renWin1);
        
        self.xImagePlaneWidget.SetInteractor( self.iren1 )
        self.yImagePlaneWidget.SetInteractor( self.iren1 )
        self.zImagePlaneWidget.SetInteractor( self.iren1 )
        
        # Set Up Camera view
        self.camera = self.renderer1.GetActiveCamera()
        self.renderer1.SetBackground(0.0, 0.0, 0.0)
        self.iren1.SetPicker(self.picker)
                
        self.origin=[]
        
    def __call__(self):       
        """ Turn Class into a callable object """
        Display()
        
        
    def dicomTransform(self, image, image_pos_pat, image_ori_pat):
        """ dicomTransform: transforms an image to a DICOM coordinate frame
        
        INPUTS:
        =======        
        image: (vtkImageData)    Input image to Transform
        image_pos_pat: (list(dicomInfo[0x0020,0x0032].value)))  Image position patient Dicom Tag
        image_ori_pat: (list(dicomInfo[0x0020,0x0037].value))   Image oreintation patient Dicom Tag
        
        OUTPUTS:
        =======
        transformed_image (vtkImageData)    Transformed imaged mapped to dicom coords frame
        transform (vtkTransform)            Transform used
        
        """ 
        # If one considers the localizer plane as a "viewport" onto the DICOM 3D coordinate space, then that viewport is described by its origin, its row unit vector, column unit vector and a normal unit vector (derived from the row and column vectors by taking the cross product). Now if one moves the origin to 0,0,0 and rotates this viewing plane such that the row vector is in the +X direction, the column vector the +Y direction, and the normal in the +Z direction, then one has a situation where the X coordinate now represents a column offset in mm from the localizer's top left hand corner, and the Y coordinate now represents a row offset in mm from the localizer's top left hand corner, and the Z coordinate can be ignored. One can then convert the X and Y mm offsets into pixel offsets using the pixel spacing of the localizer imag
        # Initialize Image orienation
        "Image Orientation Patient Matrix"
        IO = np.matrix(    [[0, 0,-1, 1],
                [1, 0, 0, 1],
                [0,-1, 0, 1],
                [0, 0, 0, 1]])
        # Assign the 6-Image orientation patient coordinates (from Dicomtags)
        IO[0,0] = image_ori_pat[0]; IO[0,1] = image_ori_pat[1]; IO[0,2] = image_ori_pat[2]; 
        IO[1,0] = image_ori_pat[3]; IO[1,1] = image_ori_pat[4]; IO[1,2] = image_ori_pat[5]; 
        
        # obtain thrid column as the cross product of column 1 y 2
        IO_col1 = [image_ori_pat[0], image_ori_pat[1], image_ori_pat[2]]
        IO_col2 = [image_ori_pat[3], image_ori_pat[4], image_ori_pat[5]]
        IO_col3 = np.cross(IO_col1, IO_col2)
        
        # assign column 3    
        IO[2,0] = IO_col3[0]; IO[2,1] = IO_col3[1]; IO[2,2] = IO_col3[2]; 
        
        IP =  np.array([0, 0, 0, 1]) # Initialization Image Position
        IP[0] = image_pos_pat[0]; IP[1] = image_pos_pat[1]; IP[2] = image_pos_pat[2];  
        IO[0,3] = -image_pos_pat[0]; IO[1,3] = -image_pos_pat[1]; IO[2,3] = -image_pos_pat[2]
               
        "Compute Volume Origin"
        origin = IP*IO.I
        
        # Create matrix 4x4
        DICOM_mat = vtk.vtkMatrix4x4();
        DICOM_mat.SetElement(0, 0, IO[0,0])
        DICOM_mat.SetElement(0, 1, IO[0,1])
        DICOM_mat.SetElement(0, 2, IO[0,2])
        DICOM_mat.SetElement(0, 3, IO[0,3])
        
        DICOM_mat.SetElement(1, 0, IO[1,0])
        DICOM_mat.SetElement(1, 1, IO[1,1])
        DICOM_mat.SetElement(1, 2, IO[1,2])
        DICOM_mat.SetElement(1, 3, IO[1,3])
        
        DICOM_mat.SetElement(2, 0, IO[2,0])
        DICOM_mat.SetElement(2, 1, IO[2,1])
        DICOM_mat.SetElement(2, 2, IO[2,2])
        DICOM_mat.SetElement(2, 3, IO[2,3])
        
        DICOM_mat.SetElement(3, 0, IO[3,0])
        DICOM_mat.SetElement(3, 1, IO[3,1])
        DICOM_mat.SetElement(3, 2, IO[3,2])
        DICOM_mat.SetElement(3, 3, IO[3,3])
        #DICOM_mat.Invert()
        
        # Set up the axes    
        transform = vtk.vtkTransform()
        transform.Concatenate(DICOM_mat)
        transform.Update()
        
        # Set up the cube (set up the translation back to zero    
        DICOM_mat_cube = vtk.vtkMatrix4x4();
        DICOM_mat_cube.DeepCopy(DICOM_mat)
        DICOM_mat_cube.SetElement(0, 3, 0)
        DICOM_mat_cube.SetElement(1, 3, 0)
        DICOM_mat_cube.SetElement(2, 3, 0)
            
        transform_cube = vtk.vtkTransform()
        transform_cube.Concatenate(DICOM_mat_cube)
        transform_cube.Update()
         
        # Change info
        # Flip along Y-Z-axis: VTK uses computer graphics convention where the first pixel in memory is shown 
        # in the lower left of the displayed image.
        flipZ_image = vtk.vtkImageFlip()
        flipZ_image.SetInput(image)
        flipZ_image.SetFilteredAxis(2)
        flipZ_image.Update() 
        
        flipY_image = vtk.vtkImageFlip()
        flipY_image.SetInput(flipZ_image.GetOutput())
        flipY_image.SetFilteredAxis(1)
        flipY_image.Update() 
          
        # Change info origin
        flipY_origin_image = vtk.vtkImageChangeInformation()
        flipY_origin_image.SetInput( flipY_image.GetOutput() );
        flipY_origin_image.SetOutputOrigin(origin[0,0], origin[0,1], origin[0,2])
        flipY_origin_image.Update()
        
        return flipY_origin_image.GetOutput(), transform_cube   
            
    def addSegment(self, lesion3D):        
       
        # Set the planes based on seg bounds
        bounds = lesion3D.GetBounds()
        print "\n Mesh DICOM bounds: "
        print "xmin, xmax= [%d, %d]" % (bounds[0], bounds[1])
        print "yin, ymax= [%d, %d]" %  (bounds[2], bounds[3]) 
        print "zmin, zmax= [%d, %d]" % (bounds[4], bounds[5])
        
        # Add ICPinit_mesh.vtk to the render
        self.mapper_mesh = vtk.vtkPolyDataMapper()
        self.mapper_mesh.SetInput( lesion3D )
        self.mapper_mesh.ScalarVisibilityOff()
        
        self.actor_mesh = vtk.vtkActor()
        self.actor_mesh.SetMapper(self.mapper_mesh)
        self.actor_mesh.GetProperty().SetColor(0, 1, 0)    #R,G,B
        self.actor_mesh.GetProperty().SetOpacity(0.6)
        self.actor_mesh.GetProperty().SetPointSize(5.0)
        self.actor_mesh.GetProperty().SetRepresentationToWireframe()
        
        self.xImagePlaneWidget.SetSliceIndex(0)
        self.yImagePlaneWidget.SetSliceIndex(0)
        self.zImagePlaneWidget.SetSliceIndex( 0 )
        
        self.renderer1.AddActor(self.actor_mesh)
        self.renWin1.Render()
                
        return
        
    def subImage(self, Images2Sub, timep):
        #subtract volumes based on indicated postS            
        sub_preMat = vtk.vtkImageMathematics()
        sub_preMat.SetOperationToSubtract()
        sub_preMat.SetInput1(Images2Sub[timep])
        sub_preMat.SetInput2(Images2Sub[0])
        sub_preMat.Update()
                    
        sub_pre = vtk.vtkImageData()
        sub_pre =  sub_preMat.GetOutput()
        # define image based on subtraction of postS -preS
        subtractedImage = sub_pre
        
        return subtractedImage
               
    def visualize(self, images, image_pos_pat, image_ori_pat, sub, postS, interact):
        
        if(sub):
            #subtract volumes based on indicated postS            
            # define image based on subtraction of postS -preS
            image = self.subImage(images, postS)
        else:
            image = images[postS]            
             
        # Proceed to build reference frame for display objects based on DICOM coords   
        [transformed_image, transform_cube] = self.dicomTransform(image, image_pos_pat, image_ori_pat)
        
        # get info from image before visualization
        transformed_image.UpdateInformation()
        self.dims = transformed_image.GetDimensions()
        print "Image Dimensions"
        print self.dims
        (xMin, xMax, yMin, yMax, zMin, zMax) = transformed_image.GetWholeExtent()
        print "Image Extension"
        print xMin, xMax, yMin, yMax, zMin, zMax
        self.spacing = transformed_image.GetSpacing()
        print "Image Spacing"
        print self.spacing
        self.origin = transformed_image.GetOrigin()
        print "Image Origin"
        print self.origin
            
        # Set up ortogonal planes
        self.xImagePlaneWidget.SetInput( transformed_image )
        self.xImagePlaneWidget.SetPlaneOrientationToXAxes()
        self.xImagePlaneWidget.SetSliceIndex(0)
        self.yImagePlaneWidget.SetInput( transformed_image )
        self.yImagePlaneWidget.SetPlaneOrientationToYAxes()
        self.yImagePlaneWidget.SetSliceIndex(0)
        self.zImagePlaneWidget.SetInput( transformed_image )
        self.zImagePlaneWidget.SetPlaneOrientationToZAxes()
        self.zImagePlaneWidget.SetSliceIndex(0)
            
        self.xImagePlaneWidget.On()
        self.yImagePlaneWidget.On()
        self.zImagePlaneWidget.On()
        
        # set up cube actor with Orientation(A-P, S-I, L-R) using transform_cube
        # Set up to ALS (+X=A, +Y=S, +Z=L) source:
        cube = vtk.vtkAnnotatedCubeActor()
        cube.SetXPlusFaceText( "L" );
        cube.SetXMinusFaceText( "R" );
        cube.SetYPlusFaceText( "A" );
        cube.SetYMinusFaceText( "P" );
        cube.SetZPlusFaceText( "S" );
        cube.SetZMinusFaceText( "I" );
        cube.SetFaceTextScale( 0.5 );
        cube.GetAssembly().SetUserTransform( transform_cube );
            
        # Set UP the axes
        axes2 = vtk.vtkAxesActor()
        axes2.SetShaftTypeToCylinder();
        #axes2.SetUserTransform( transform_cube );         
        axes2.SetTotalLength( 1.5, 1.5, 1.5 );
        axes2.SetCylinderRadius( 0.500 * axes2.GetCylinderRadius() );
        axes2.SetConeRadius( 1.025 * axes2.GetConeRadius() );
        axes2.SetSphereRadius( 1.500 * axes2.GetSphereRadius() );
    
        tprop2 = axes2.GetXAxisCaptionActor2D()
        tprop2.GetCaptionTextProperty();
    
        assembly = vtk.vtkPropAssembly();
        assembly.AddPart( axes2 );
        assembly.AddPart( cube );
        
        widget = vtk.vtkOrientationMarkerWidget();
        widget.SetOutlineColor( 0.9300, 0.5700, 0.1300 );
        widget.SetOrientationMarker( assembly );
        widget.SetInteractor( self.iren1 );
        widget.SetViewport( 0.0, 0.0, 0.4, 0.4 );
        widget.SetEnabled( 1 );
        widget.InteractiveOff();
                    
        # Create a text property for both cube axes
        tprop = vtk.vtkTextProperty()
        tprop.SetColor(1, 1, 1)
        tprop.ShadowOff()
        
        # Create a vtkCubeAxesActor2D.  Use the outer edges of the bounding box to
        # draw the axes.  Add the actor to the renderer.
        axes = vtk.vtkCubeAxesActor2D()
        axes.SetInput(transformed_image)
        axes.SetCamera(self.renderer1.GetActiveCamera())
        axes.SetLabelFormat("%6.4g")
        axes.SetFlyModeToOuterEdges()
        axes.SetFontFactor(1.2)
        axes.SetAxisTitleTextProperty(tprop)
        axes.SetAxisLabelTextProperty(tprop)      
        self.renderer1.AddViewProp(axes)
        
        ############
        # bounds and initialize camera
        bounds = transformed_image.GetBounds()
        self.renderer1.ResetCamera(bounds)    
        self.renderer1.ResetCameraClippingRange()
        self.camera.SetViewUp(0.0,-1.0,0.0)
        self.camera.Azimuth(315)
        
        # Initizalize
        self.renWin1.Render()
        self.renderer1.Render()
        
        if(interact==True):
            self.iren1.Start()  
                            
        return

