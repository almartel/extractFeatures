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
from numpy import *
import re
import collections
import math

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
        self.renWin1.AddRenderer(self.renderer1)
        self.iren1.SetRenderWindow(self.renWin1)
        
        self.xImagePlaneWidget.SetInteractor( self.iren1 )
        self.yImagePlaneWidget.SetInteractor( self.iren1 )
        self.zImagePlaneWidget.SetInteractor( self.iren1 )
        
        # Set Up Camera view
        self.camera = self.renderer1.GetActiveCamera()
        self.renderer1.SetBackground(0.0, 0.0, 0.0)
        self.iren1.SetPicker(self.picker)
                
        self.T1origin = [0,0,0]
        self.T2origin = [0,0,0]
        self.T2extent = [0,0,0,0,0,0]
        self.T1extent = [0,0,0,0,0,0]
        self.T1spacing = [0,0,0]
        
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
        IO = matrix(    [[0, 0, 0, 0],
                [0, 0, 0, 0],
                [0, 0, 0, 0],
                [0, 0, 0, 1]])
        # Assign the 6-Image orientation patient coordinates (from Dicomtags)
        IO[0,0] = image_ori_pat[0]; IO[0,1] = image_ori_pat[1]; IO[0,2] = image_ori_pat[2]; 
        IO[1,0] = image_ori_pat[3]; IO[1,1] = image_ori_pat[4]; IO[1,2] = image_ori_pat[5]; 
        
        # obtain thrid column as the cross product of column 1 y 2
        IO_col1 = [image_ori_pat[0], image_ori_pat[1], image_ori_pat[2]]
        IO_col2 = [image_ori_pat[3], image_ori_pat[4], image_ori_pat[5]]
        IO_col3 = cross(IO_col1, IO_col2)
        
        # assign column 3    
        IO[2,0] = IO_col3[0]; IO[2,1] = IO_col3[1]; IO[2,2] = IO_col3[2]; 
        
        IP =  array([0, 0, 0, 1]) # Initialization Image Position
        IP[0] = image_pos_pat[0]; IP[1] = image_pos_pat[1]; IP[2] = image_pos_pat[2];  
        IO[0,3] = image_pos_pat[0]; IO[1,3] = image_pos_pat[1]; IO[2,3] = image_pos_pat[2]
               
        print "Compute: Position, Orientation, matrix, & Volume Origin"
        print image_pos_pat
        print image_ori_pat
        print IO
        origin = IP*IO.I
        print origin[0,0], origin[0,1], origin[0,2]
        
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
        
        transformed_image = flipY_origin_image.GetOutput()
        
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
        
        return transformed_image, transform_cube   
        
    def mhaTransform(self, image, image_pos_pat, image_ori_pat):
        """ mhaTransform: transforms an image from mha coordinate frame"""
        "Image Orientation Patient Matrix"
        IO = matrix(    [[0, 0, 0, 0],
                [0, 0, 0, 0],
                [0, 0, 0, 0],
                [0, 0, 0, 1]])
        # Assign the 6-Image orientation patient coordinates (from Dicomtags)
        IO[0,0] = image_ori_pat[0]; IO[0,1] = image_ori_pat[1]; IO[0,2] = image_ori_pat[2]; 
        IO[1,0] = image_ori_pat[3]; IO[1,1] = image_ori_pat[4]; IO[1,2] = image_ori_pat[5]; 
        
        # obtain thrid column as the cross product of column 1 y 2
        IO_col1 = [image_ori_pat[0], image_ori_pat[1], image_ori_pat[2]]
        IO_col2 = [image_ori_pat[3], image_ori_pat[4], image_ori_pat[5]]
        IO_col3 = cross(IO_col1, IO_col2)
        
        # assign column 3    
        IO[2,0] = IO_col3[0]; IO[2,1] = IO_col3[1]; IO[2,2] = IO_col3[2]; 
        
        IP =  array([0, 0, 0, 1]) # Initialization Image Position
        IP[0] = image_pos_pat[0]; IP[1] = image_pos_pat[1]; IP[2] = image_pos_pat[2];  
        IO[0,3] = image_pos_pat[0]; IO[1,3] = image_pos_pat[1]; IO[2,3] = image_pos_pat[2]
               
        print "Compute: Position, Orientation, matrix, & Volume Origin"
        print image_pos_pat
        print image_ori_pat
        print IO
        origin = IP*IO.I
        print origin[0,0], origin[0,1], origin[0,2]
        
        # Change info
        # No need to Flip along Y-Z-axis like form DICOM          
        # Change info origin
        origin_image = vtk.vtkImageChangeInformation()
        origin_image.SetInput( image )
        origin_image.SetOutputOrigin(origin[0,0], origin[0,1], origin[0,2])
        origin_image.Update()
        
        transformed_image = origin_image.GetOutput()
        
        transformed_image.UpdateInformation()
        self.T2dims = transformed_image.GetDimensions()
        print "Image Dimensions"
        print self.T2dims
        (xMin, xMax, yMin, yMax, zMin, zMax) = transformed_image.GetWholeExtent()
        print "Image Extension"
        print xMin, xMax, yMin, yMax, zMin, zMax
        self.T2spacing = transformed_image.GetSpacing()
        print "Image Spacing"
        print self.T2spacing
        self.T2origin = transformed_image.GetOrigin()
        print "Image Origin"
        print self.T2origin
        
        return transformed_image 
        
        
          
    def addSegment(self, lesion3D, color, interact):        
        '''Add segmentation to current display'''
        # Set the planes based on seg bounds
        self.lesion_bounds = lesion3D.GetBounds()
        print "\n Mesh DICOM bounds: "
        print "xmin, xmax= [%d, %d]" % (self.lesion_bounds[0], self.lesion_bounds[1])
        print "yin, ymax= [%d, %d]" %  (self.lesion_bounds[2], self.lesion_bounds[3]) 
        print "zmin, zmax= [%d, %d]" % (self.lesion_bounds[4], self.lesion_bounds[5])
        
        ### GEt semgnetation information
        self.no_pts_segm = lesion3D.GetNumberOfPoints()
        print "no pts %d" % self.no_pts_segm
        
        # get VOI volume
        VOI_massProperty = vtk.vtkMassProperties()
        VOI_massProperty.SetInput(lesion3D)
        VOI_massProperty.Update()
               
        # VTK is unitless. The units you get out are the units you put in.
        # If your input polydata has points defined in terms of millimetres, then
        # the volume will be in cubic millimetres. 
        self.VOI_vol = VOI_massProperty.GetVolume() # mm3
        self.VOI_surface = VOI_massProperty.GetSurfaceArea() # mm2
    
        # just print the results
        print "\nVolume lesion = ", self.VOI_vol
        print "Surface lesion  = ", self.VOI_surface
        
        # Calculate the effective diameter of the surface D=2(sqrt3(3V/(4pi))) 
        diam_root = (3*self.VOI_vol)/(4*pi)
        self.VOI_efect_diameter = 2*pow(diam_root,1.0/3) 
        print "VOI_efect_diameter = ", self.VOI_efect_diameter
            
        centerOfMassFilter = vtk.vtkCenterOfMass()
        centerOfMassFilter.SetInput( lesion3D )
        centerOfMassFilter.SetUseScalarsAsWeights(False)
        centerOfMassFilter.Update()
        
        # centroid of lesion 
        self.lesion_centroid = [0,0,0]
        self.lesion_centroid = centerOfMassFilter.GetCenter()
        print "lesion_centroid = ", self.lesion_centroid
                
        # Add ICPinit_mesh.vtk to the render
        self.mapper_mesh = vtk.vtkPolyDataMapper()
        self.mapper_mesh.SetInput( lesion3D )
        self.mapper_mesh.ScalarVisibilityOff()
        
        self.actor_mesh = vtk.vtkActor()
        self.actor_mesh.SetMapper(self.mapper_mesh)
        self.actor_mesh.GetProperty().SetColor(color)    #R,G,B
        self.actor_mesh.GetProperty().SetOpacity(0.3)
        self.actor_mesh.GetProperty().SetPointSize(5.0)
        self.actor_mesh.GetProperty().SetRepresentationToWireframe()
        
        self.xImagePlaneWidget.SetSliceIndex(0)
        self.yImagePlaneWidget.SetSliceIndex(0)
        self.zImagePlaneWidget.SetSliceIndex( 0 )
        
        self.renderer1.AddActor(self.actor_mesh)
        
        # Initizalize
        self.renderer1.Modified()
        self.renWin1.Render()
        self.renderer1.Render()
        
        if(interact==True):
            self.iren1.Start()
            
                
        return 
        
        
    def subImage(self, Images2Sub, timep):
        '''subtract volumes based on indicated postS'''
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
        

    def display_pick(self, images, image_pos_pat, image_ori_pat, postS, LesionZslice):
        '''Display a z-slice and use picker to pick coordinates with a mouse right-click'''
        #subtract volumes based on indicated postS            
        # define image based on subtraction of postS -preS
        image = self.subImage(images, postS)    

        # Proceed to build reference frame for display objects based on DICOM coords   
        [transformed_image, transform_cube] = self.dicomTransform(image, image_pos_pat, image_ori_pat)
                
        # Calculate the center of the volume
        transformed_image.UpdateInformation() 
    
        # Set up ortogonal planes
        self.xImagePlaneWidget.SetInput( transformed_image )
        self.yImagePlaneWidget.SetInput( transformed_image )
        self.zImagePlaneWidget.SetInput( transformed_image )
        
        self.zImagePlaneWidget.SetSliceIndex( LesionZslice )
        self.xImagePlaneWidget.On()
        self.yImagePlaneWidget.On()
        self.zImagePlaneWidget.On()
        
        ############
        self.textMapper = vtk.vtkTextMapper()
        tprop = self.textMapper.GetTextProperty()
        tprop.SetFontFamilyToArial()
        tprop.SetFontSize(10)
        tprop.BoldOn()
        tprop.ShadowOn()
        tprop.SetColor(1, 0, 0)
           
        # initialize 
        self.seeds = vtk.vtkPoints()  
        self.textActor = vtk.vtkActor2D()
        self.textActor.VisibilityOff() 
        self.textActor.SetMapper(self.textMapper)

        # Initizalize
        self.iren1.SetPicker(self.picker)
        self.picker.AddObserver("EndPickEvent", self.annotatePick)
        self.renWin1.Render()
        self.renderer1.Render()
        self.iren1.Start()
                
        return self.seeds
        
    
    def annotatePick(self, object, event):
        '''Auxiliary function of picker, print coords chosen'''        
        if(self.picker.GetCellId() < 0):
            self.textActor.VisibilityOff()     
        else:
            selPt = self.picker.GetSelectionPoint()
            pickPos = self.picker.GetPickPosition()
            self.seeds.InsertNextPoint(pickPos[0], pickPos[1], pickPos[2] )
            print pickPos
        
            self.textMapper.SetInput("(%.6f, %.6f, %.6f)"%pickPos)
            self.textActor.SetPosition(selPt[:2])
            self.textActor.VisibilityOn()
        
        return 
      
      
    def visualize(self, images, image_pos_pat, image_ori_pat, sub, postS, interact):
        '''Display and render volumes, reference frames, actors and widgets'''
        if(sub):
            #subtract volumes based on indicated postS            
            # define image based on subtraction of postS -preS
            image = self.subImage(images, postS)
        else:
            image = images[postS]            
             
        # Proceed to build reference frame for display objects based on DICOM coords   
        [self.transformed_image, transform_cube] = self.dicomTransform(image, image_pos_pat, image_ori_pat)
        
        # get info from image before visualization
        self.transformed_image.UpdateInformation()
        self.dims = self.transformed_image.GetDimensions()
        print "Image Dimensions"
        print self.dims
        self.T1spacing = self.transformed_image.GetSpacing()
        print "Image Spacing"
        print self.T1spacing
        self.T1origin = self.transformed_image.GetOrigin()
        print "Image Origin"
        print self.T1origin
        self.T1extent = list(self.transformed_image.GetWholeExtent())
        print "Image Extent"
        print self.T1extent
            
        # Set up ortogonal planes
        self.xImagePlaneWidget.SetInput( self.transformed_image )
        self.xImagePlaneWidget.SetPlaneOrientationToXAxes()
        self.xImagePlaneWidget.SetSliceIndex(0)
        self.yImagePlaneWidget.SetInput( self.transformed_image )
        self.yImagePlaneWidget.SetPlaneOrientationToYAxes()
        self.yImagePlaneWidget.SetSliceIndex(0)
        self.zImagePlaneWidget.SetInput( self.transformed_image )
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
        axes.SetInput(self.transformed_image)
        axes.SetCamera(self.renderer1.GetActiveCamera())
        axes.SetLabelFormat("%6.4g")
        axes.SetFlyModeToOuterEdges()
        axes.SetFontFactor(1.2)
        axes.SetAxisTitleTextProperty(tprop)
        axes.SetAxisLabelTextProperty(tprop)      
        self.renderer1.AddViewProp(axes)
        
        ############
        # bounds and initialize camera
        bounds = self.transformed_image.GetBounds()
        self.renderer1.ResetCamera(bounds)    
        self.renderer1.ResetCameraClippingRange()
        self.camera.SetViewUp(0.0,-1.0,0.0)
        self.camera.Azimuth(315)
        
        # Initizalize
        self.renWin1.Modified()
        self.renWin1.Render()
        self.renderer1.Render()
        
        if(interact==True):
            interactor = self.renWin1.GetInteractor()
            interactor.Start()
                            
        return
        
    def addT2visualize(self, T2images, image_pos_pat, image_ori_pat, T2dims, T2spacing, interact):
        '''Added to build second reference frame and display T2 overlayed into T1 reference frame'''
        # Proceed to build reference frame for display objects based on DICOM coords   
        [self.transformed_T2image, transform_cube] = self.dicomTransform(T2images, image_pos_pat, image_ori_pat)
        
        self.T2origin = list(self.transformed_T2image.GetOrigin())
        print "T2 Extent"
        self.T2extent = list(self.transformed_T2image.GetExtent())
        print self.T2extent        
        
        # Set up ortogonal planes
        self.xImagePlaneWidget.SetInput( self.transformed_T2image )
        self.xImagePlaneWidget.SetSliceIndex(0)
        self.yImagePlaneWidget.SetInput( self.transformed_T2image )
        self.yImagePlaneWidget.SetSliceIndex(0)
        self.zImagePlaneWidget.SetInput( self.transformed_T2image )
        self.zImagePlaneWidget.SetSliceIndex(0)
                    
        # Create a text property for both cube axes
        tprop = vtk.vtkTextProperty()
        tprop.SetColor(0, 1, 1)
        tprop.ShadowOff()
        
        # Update the reneder window to receive new image !Important*****
        self.renderer1.Modified()
        self.renWin1.Modified()
        
        # Create a vtkCubeAxesActor2D.  Use the outer edges of the bounding box to
        # draw the axes.  Add the actor to the renderer.
        axesT2 = vtk.vtkCubeAxesActor2D()
        axesT2.SetInput(self.transformed_T2image)
        axesT2.SetCamera(self.renderer1.GetActiveCamera())
        axesT2.SetLabelFormat("%6.4g")
        axesT2.SetFlyModeToOuterEdges()
        axesT2.SetFontFactor(1.2)
        axesT2.SetAxisTitleTextProperty(tprop)
        axesT2.SetAxisLabelTextProperty(tprop)      
        self.renderer1.AddViewProp(axesT2)

        ############                
        if(interact==True):
            interactor = self.renWin1.GetInteractor()
            interactor.Start()
            
        return
        
    def addT2transvisualize(self, T2images, image_pos_pat, image_ori_pat, T2dims, T2spacing, interact):
        '''Added to build second reference frame and display T2 overlayed into T1 reference frame'''
        # Proceed to build reference frame for display objects based on DICOM coords   
        [transformed_T2image, transform_cube] = self.dicomTransform(T2images, image_pos_pat, image_ori_pat)
        
        self.T2origin = list(transformed_T2image.GetOrigin())
        print "T2 Extent"
        self.T2extent = list(transformed_T2image.GetExtent())
        print self.T2extent        
        
        alignR = int(raw_input('\nAlign right? Yes:1 or align with T1w: !=1 : '))
        if alignR:
            zf1 = self.T1spacing[2]*self.T1extent[5] + self.T1origin[2]
            self.T2origin[2] = zf1 - T2spacing[2]*self.T2extent[5] # this is z-span
        else:
            self.T2origin[2] = self.T1origin[2]
                
        # Change info origin
        transformedInfo_T2image = vtk.vtkImageChangeInformation()
        transformedInfo_T2image.SetInput( transformed_T2image )
        transformedInfo_T2image.SetOutputOrigin(self.T2origin)
        transformedInfo_T2image.Update()
        
        self.transformed_T2image = transformedInfo_T2image.GetOutput()
        
        # Set up ortogonal planes
        self.xImagePlaneWidget.SetInput( self.transformed_T2image )
        self.xImagePlaneWidget.SetSliceIndex(0)
        self.yImagePlaneWidget.SetInput( self.transformed_T2image )
        self.yImagePlaneWidget.SetSliceIndex(0)
        self.zImagePlaneWidget.SetInput( self.transformed_T2image )
        self.zImagePlaneWidget.SetSliceIndex(0)
                    
        # Create a text property for both cube axes
        tprop = vtk.vtkTextProperty()
        tprop.SetColor(0.5, 0.5, 0)
        tprop.ShadowOff()
        
        # Update the reneder window to receive new image !Important*****
        self.renderer1.Modified()
        self.renWin1.Modified()
        
        # Create a vtkCubeAxesActor2D.  Use the outer edges of the bounding box to
        # draw the axes.  Add the actor to the renderer.
        axesT2 = vtk.vtkCubeAxesActor2D()
        axesT2.SetInput(self.transformed_T2image)
        axesT2.SetCamera(self.renderer1.GetActiveCamera())
        axesT2.SetLabelFormat("%6.4g")
        axesT2.SetFlyModeToOuterEdges()
        axesT2.SetFontFactor(1.2)
        axesT2.SetAxisTitleTextProperty(tprop)
        axesT2.SetAxisLabelTextProperty(tprop)      
        self.renderer1.AddViewProp(axesT2)
        
        ### Update T2Images
        t_T2images = vtk.vtkImageChangeInformation()
        t_T2images.SetInput( T2images )
        t_T2images.SetOutputOrigin(self.T2origin)
        t_T2images.Update()        

        ############                
        if(interact==True):
            interactor = self.renWin1.GetInteractor()
            interactor.Start()
            
        return 

        
    def visualizemha(self, img_path, image_pos_pat, image_ori_pat, interact):
        '''Display and render volumes, reference frames, actors and widgets'''  
             
        # Proceed to build reference frame for display objects based on DICOM coords   
        mhareader = vtk.vtkMetaImageReader()
        mhareader.SetFileName(img_path)
        mhareader.Update()
        
        self.warpT2_mha = self.mhaTransform(mhareader.GetOutput(), image_pos_pat, image_ori_pat)
        
        self.T2origin = list(self.warpT2_mha.GetOrigin())
        print "T2 Extent"
        self.T2extent = list(self.warpT2_mha.GetExtent())
        print self.T2extent        

        # Initizalize
        self.xImagePlaneWidget.Modified()
        self.yImagePlaneWidget.Modified()
        self.zImagePlaneWidget.Modified()
                
        # Set up ortogonal planes
        self.xImagePlaneWidget.SetInput(self.warpT2_mha )
        self.xImagePlaneWidget.SetSliceIndex(0)
        self.yImagePlaneWidget.SetInput(self.warpT2_mha )
        self.yImagePlaneWidget.SetSliceIndex(0)
        self.zImagePlaneWidget.SetInput( self.warpT2_mha )
        self.zImagePlaneWidget.SetSliceIndex(0)
                    
        # Create a text property for both cube axes
        tprop = vtk.vtkTextProperty()
        tprop.SetColor(0.5, 0.5, 0)
        tprop.ShadowOff()
    
        # Create a vtkCubeAxesActor2D.  Use the outer edges of the bounding box to
        # draw the axes.  Add the actor to the renderer.
        axesT2 = vtk.vtkCubeAxesActor2D()
        axesT2.SetInput(self.warpT2_mha)
        axesT2.SetCamera(self.renderer1.GetActiveCamera())
        axesT2.SetLabelFormat("%6.4g")
        axesT2.SetFlyModeToOuterEdges()
        axesT2.SetFontFactor(1.2)
        axesT2.SetAxisTitleTextProperty(tprop)
        axesT2.SetAxisLabelTextProperty(tprop)      
        self.renderer1.AddViewProp(axesT2)     
        
        if(interact==True):
            interactor = self.renWin1.GetInteractor()
            interactor.Start()     
                            
        return

        
    def extract_annot(self, list_annots):
        '''Parse list of annotations, put markers according to notes and color code according to sequence order'''
        
        annots_dict_list=[]
        a_count = 1
        for one_annot in list_annots:
            annots_dict = {}
            print one_annot
            
            # iterate throuhg attributes of annotation
            annots_dict['AccessionNumber'] = one_annot[one_annot.find("':")+4:one_annot.find("',")] 
            one_annot = one_annot[one_annot.find("',")+2:]
            
            annots_dict['SeriesDate'] = one_annot[one_annot.find("':")+4:one_annot.find("',")] 
            one_annot = one_annot[one_annot.find("',")+2:]
            
            annots_dict['SeriesNumber'] = one_annot[one_annot.find("':")+4:one_annot.find("',")] 
            one_annot = one_annot[one_annot.find("',")+2:]
            
            annots_dict['SliceLocation'] = one_annot[one_annot.find("':")+4:one_annot.find("',")] 
            one_annot = one_annot[one_annot.find("',")+2:]
            
            annots_dict['SeriesDescription'] = one_annot[one_annot.find("':")+4:one_annot.find("',")] 
            one_annot = one_annot[one_annot.find("',")+2:]
            
            annots_dict['PatientID'] = one_annot[one_annot.find("':")+4:one_annot.find("',")] 
            one_annot = one_annot[one_annot.find("',")+2:]
            
            annots_dict['StudyID'] = one_annot[one_annot.find("':")+4:one_annot.find("',")] 
            one_annot = one_annot[one_annot.find("',")+2:]
                            
            # get the type of annotation: e.g CALLIPER, ELLIPSE, ARROW
            annots_dict['note'] = one_annot[one_annot.find("':")+4:one_annot.find("\\")]
            
            # extract annotation coordinate location
            coords_str = one_annot[one_annot.find("\\"):one_annot.find("',")]
            non_dec = re.compile(r'[^\d.]+')
            coords = non_dec.sub(',', coords_str).split(',')
            print coords
            if coords != ['']:
                annots_dict['xi']=float(coords[1])
                annots_dict['yi']=float(coords[2])
                annots_dict['xf']=float(coords[3])
                annots_dict['yf']=float(coords[4])
            
            # finish last attribute of annotations
            one_annot = one_annot[one_annot.find("',")+2:]
            annots_dict['SeriesInstanceUID'] = one_annot[one_annot.find("':")+4:]
            annots_dict_list.append(annots_dict)
            a_count+=1
            
        return annots_dict_list
        
    
    def display_annot(self, images,  image_pos_pat, image_ori_pat, annots_dict_list, interact):
        '''Display list of annotations, put markers according to notes and color code according to sequence order'''
        # define image based on subtraction of postS -preS          
        image = images[4]            
             
        # Proceed to build reference frame for display objects based on DICOM coords   
        [self.transformed_image, transform_cube] = self.dicomTransform(image, image_pos_pat, image_ori_pat)
        
        # supports 56 annotations
        color_list = [ [0,0,0], [1,0,0], [0,0,1], [0,1,1], [1,1,0], [1,0,1], [1,1,1], [0.5,0.5,0.5], [0.5,0.5,0],    [1,0.2,0], [1,1,0], [0,1,1],[0,1,0.6], [0,1,1], [1,0.4,0], [0.6,0,0.2], [1,1,0], [0,1,1],[0,1,0], [0,1,1], [1,0,0.1], [1,0.4,0],[0,1,1],[0,1,0], [0,1,1], [1,0,0.1], [1,0.4,0], 
                        [0,0,0], [1,0,0], [0,0,1], [0,1,1], [1,1,0], [1,0,1], [1,1,1], [0.5,0.5,0.5], [0.5,0.5,0],    [1,0.2,0], [1,1,0], [0,1,1],[0,1,0.6], [0,1,1], [1,0.4,0], [0.6,0,0.2], [1,1,0], [0,1,1],[0,1,0], [0,1,1], [1,0,0.1], [1,0.4,0],[0,1,1],[0,1,0], [0,1,1], [1,0,0.1], [1,0.4,0] ]
        a_count = 1
        annot_pts_lbl = vtk.vtkPoints()
        for annots_dict in annots_dict_list:
            try:
                float(annots_dict['SliceLocation'])
                print '\n=========#'+str(a_count)
                print annots_dict            
                ######################
                ## Display in graphics
                ######################
                im_pt = [0,0,0]
                ijk = [0,0,0]
                pco = [0,0,0]
                pi_2display=[0,0,0]
                pf_2display=[0,0,0]
                
                # extract Slice locaton
                pixId_sliceloc = self.transformed_image.FindPoint(self.origin[0], self.origin[1], float(annots_dict['SliceLocation']))
                self.transformed_image.GetPoint(pixId_sliceloc, im_pt) 
                io = self.transformed_image.ComputeStructuredCoordinates( im_pt, ijk, pco)
                
                # mark initial
                print "Point init"
                ijk[0] = int(annots_dict['xi'])
                ijk[1] = int(annots_dict['yi'])
                print ijk
                annots_dict['pi_ijk']=ijk
                pixId = self.transformed_image.ComputePointId(ijk)
                pi_2display = self.transformed_image.GetPoint(pixId)
                annots_dict['pi_2display']=pi_2display
                print pi_2display
                
                # mark final
                print "Point final"
                ijk[0] = int(annots_dict['xf'])
                ijk[1] = int(annots_dict['yf'])
                print ijk
                annots_dict['pf_ijk']=ijk
                pixId = self.transformed_image.ComputePointId(ijk)
                pf_2display = self.transformed_image.GetPoint(pixId)
                annots_dict['pf_2display']=pf_2display
                print pf_2display
                
                # Create a graphial line between the two points
                annot_pts = vtk.vtkPoints()
                annot_pts.InsertNextPoint(pi_2display)
                annot_pts.InsertNextPoint(pf_2display)
      
                annot_ln = vtk.vtkLine()
                annot_ln.GetPointIds().SetId(0,0)
                annot_ln.GetPointIds().SetId(1,1)
                note_lines = vtk.vtkCellArray()
                note_lines.InsertNextCell(annot_ln)
                
                annot_poly = vtk.vtkPolyData()
                annot_poly.SetPoints(annot_pts)
                annot_poly.SetLines(note_lines)
                annot_poly.Update()
    
                # Create mappers and actors
                annot_mapper_mesh = vtk.vtkPolyDataMapper()
                annot_mapper_mesh.SetInput( annot_poly )
                            
                self.annot_actor = vtk.vtkActor()
                self.annot_actor.SetMapper(annot_mapper_mesh)
                self.annot_actor.GetProperty().SetColor(color_list[a_count])
                self.annot_actor.GetProperty().SetLineWidth(3)
                self.annot_actor.GetProperty().SetOpacity(0.6)
                self.annot_actor.GetProperty().SetPointSize(7.0)
                self.annot_actor.GetProperty().SetRepresentationToWireframe()
                
                ############
                # Generate data arrays containing label ids
                annot_pts_lbl.InsertPoint(a_count, pi_2display)
                                   
                # add annotation to scene
                print annots_dict
                self.renderer1.AddActor(self.annot_actor)
                        
                # Initizalize
                self.renWin1.Render()
                self.renderer1.Render()
                a_count +=1
            
            except ValueError:
                a_count +=1
                pass
            
        ############
        print annot_pts_lbl.GetNumberOfPoints()
        annot_lbl_poly = vtk.vtkPolyData()
        annot_lbl_poly.SetPoints(annot_pts_lbl)
        annot_lbl_poly.Update()
                
        # Generate data arrays containing label ids                
        ids = vtk.vtkIdFilter()
        ids.SetInput(annot_lbl_poly)
        ids.PointIdsOn()
        ids.CellIdsOff()
        
        # Create labels for points
        visPts = vtk.vtkSelectVisiblePoints()
        visPts.SetInput(ids.GetOutput())
        visPts.SetRenderer(self.renderer1)
        visPts.SelectionWindowOff()
        
        # Create the mapper to display the point ids.  Specify the format to
        # use for the labels.  Also create the associated actor.
        ldm = vtk.vtkLabeledDataMapper()
        ldm.SetInput(visPts.GetOutput())
        ldm.SetLabelModeToLabelFieldData()
        pointLabels = vtk.vtkActor2D()
        pointLabels.SetMapper(ldm)
                        
        # initialize 
        self.renderer1.AddActor2D(pointLabels)
        
        print "\n====== Color codes:\n "
        print '\033[1;31m 1) Red '
        print '\033[1;34m 2) Blue '
        print '\033[1;36m 3) Cyan '
        print '\033[1;33m 4) Yellow '
        print '\033[1;35m 5) Fushia '
        print '\033[1;37m 6) White '
        print '\033[1;30m 7) Gray '
        print '\033[1;0m'
        
        ############                
        if(interact==True):
            interactor = self.renWin1.GetInteractor()
            interactor.Start()
            
        return
        
        
    def extract_segment_dims(self, lesion3D):
        '''Extract mayor dimensions of automatic lesion segmentation for validation'''
        # define image based on subtraction of postS -preS          
        axis_lenghts = array( [0,0,0] ).astype(float)
        
        l_bounds = lesion3D.GetBounds()
        axis_lenghts[0] = (l_bounds[1]-l_bounds[0])**2     
        axis_lenghts[1] = (l_bounds[3]-l_bounds[2])**2
        axis_lenghts[2] = (l_bounds[5]-l_bounds[4])**2 
        
        return axis_lenghts
        