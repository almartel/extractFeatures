# -*- coding: utf-8 -*-
"""
================
glcm3d
================

-standard 2D directions:

========== ===============
[0 0 1]    0 degrees
[0 1 -1]   45 degrees
[0 1 0]    90 degrees
[0 1 1]    135 degrees
========== ===============

- additional 9 directions 

=========== ============= =============           
            horizontal    vertical
[1 0 -1]    0 degrees     45 degrees
[1 0 0]     straight up
[1 0 1]     0 degree      135 degrees
[1 1 0]     90 degrees    45 degrees
[1 -1 0]     90 degrees    135 degrees
[1 1 -1]    45 degrees    45 degrees
[1 -1 1]    45 degrees    135 degrees
[1 1 1]     135 degrees   45 degrees
[1 -1 -1]   135 degrees   135 degrees
=========== ============= =============   

C++ Implementation
========================

.. code:: cpp module

    #include <stdio.h>
    #include <string.h>
    #include <math.h> 
    
    void glcm3d_c(int levels, int sizez, int sizex, int sizey, int *arr, int h, int *offset, int h1, int w1, int *result){
        
        int i,j,k,aux,aux2,rows,cols,slices;
        for(k=0;k<sizez;k++){
            for(i=0;i<sizex;i++){
                for(j=0;j<sizey;j++){
                    aux = arr[k*sizex*sizey + i*sizey+j];
                    //printf("%d (%d,%d,%d)\n",aux,k,i,j);
                    slices = k+offset[0];
                    rows = i+offset[1];
                    cols = j+offset[2];
                    //printf("(%d,%d,%d \n)",slices,rows,cols);
                    if(slices >= 0 && slices < sizez){
                        if(rows >= 0 && rows < sizex){
                            if(cols >= 0 && cols < sizey){
                                aux2 = arr[slices*sizex*sizey + rows*sizey+cols];
                                //printf("%d (%d,%d,%d)\n",aux2,slices,rows,cols);
                                if(aux >= 0 && aux < levels){
                                    if(aux2 >= 0 && aux2 < levels){
                                        result[aux*levels+aux2]+=1;
                                    }
                                 }
                             }
                         }  
                     } 
                 }
            }
       }
   }

    
.. code:: cpp header
    :show_code: yes
    :show_output: yes
    
    
    void glcm3d_c(int levels, int DIM1, int DIM2,int DIM3, int *IN_ARRAY3, int DIM1, int *IN_ARRAY1,int DIM1, int DIM2, int *INPLACE_ARRAY2);
    


.. code:: python module

   import glcm3d as bib
   import numpy as np
   
   def glcm3(levels,offsets,img):
      result = np.zeros((levels,levels), dtype = np.int32)    
      bib.glcm3d_c(levels, img,offsets, result)
      return result
      
      
.. code:: python 

    import numpy as np
    from glcm3d import glcm3
    
    arr = zeros((3,5,6), dtype=np.uint8)
    #arr = np.array([[[1,2,2,0,0,1],[0,0,1,2,2,1],[1,1,0,0,0,2],[1,1,1,2,2,2],[1,1,2,2,0,0]],
    #       [[1,2,2,0,0,1],[0,0,1,2,2,1],[1,1,0,0,0,2],[1,1,1,2,2,2],[1,1,2,2,0,0]]], dtype=np.uint8)
    arr = np.array([[[1,2,2,0,0,1],[0,0,1,2,2,1],[1,1,0,0,0,2],[1,1,1,2,2,2],[1,1,2,2,0,0]],
           [[1,2,2,0,0,1],[0,0,1,2,2,1],[1,1,0,0,0,2],[1,1,1,2,2,2],[1,1,2,2,0,0]],
           [[1,2,2,0,0,1],[0,0,1,2,2,1],[1,1,0,0,0,2],[1,1,1,2,2,2],[1,1,2,2,0,0]]], dtype=np.uint8)
    print arr.shape
    levels = int(arr.max()+1)

Standard directions
======================

.. code:: python 
  
    offset = np.array([0,0,1], dtype = np.int32) # orientation 0 degrees (example same slices: equal to Nslices*0degree 2D case)
    result = glcm3(levels, offset, arr)
    print result
    print 
    offset = np.array([0,1,-1], dtype = np.int32) # orientation 45 degrees (example same slices: equal to Nslices*45degree 2D case)
    result = glcm3(levels, offset, arr)
    print result
    print 
    offset = np.array([0,1,0], dtype = np.int32) # orientation 90 degrees (example same slices: equal to Nslices*90degree 2D case)
    result = glcm3(levels, offset, arr)
    print result
    print 
    offset = np.array([0,1,1], dtype = np.int32) # orientation 135 degrees (example same slices: equal to Nslices*135degree 2D case)
    result = glcm3(levels, offset, arr)
    print result
    print

Addiction 9 directions
========================

.. code:: python 
  
    offset = np.array([1,0,-1], dtype = np.int32) # 0 degrees/45 degrees (example same slices: equal to (Nslices-1)*0degree 2D case)
    result = glcm3(levels, offset, arr)
    print result
    print 
    offset = np.array([1,0,0], dtype = np.int32) # straight up (example same slices: equal to np.unique())
    result = glcm3(levels, offset, arr)
    print result
    print 
    offset = np.array([1,0,1], dtype = np.int32) # 0 degree/135 degrees (example same slices: equal to (Nslices-1)*transpose of 0degree 2D case)
    result = glcm3(levels, offset, arr)
    print result
    print 
    offset = np.array([1,1,0], dtype = np.int32) # 90 degrees/45 degrees (example same slices: equal to (Nslices-1)*90 degree 2D case)
    result = glcm3(levels, offset, arr)
    print result
    print
    
    offset = np.array([1,-1,0], dtype = np.int32) # 90 degrees/135 degrees (example same slices: equal to (Nslices-1)*transpose of 90 degree 2D case)
    result = glcm3(levels, offset, arr)
    print result
    print 
    offset = np.array([1,1,-1], dtype = np.int32) # 45 degrees/45 degrees (example same slices: equal to (Nslices-1)*45 degree 2D case)
    result = glcm3(levels, offset, arr)
    print result
    print 
    offset = np.array([1,-1,1], dtype = np.int32) # 45 degree/135 degrees (example same slices: equal to (Nslices-1)*transpose of 45 degree 2D case)
    result = glcm3(levels, offset, arr)
    print result
    print 
    offset = np.array([1,1,1], dtype = np.int32) # 135 degrees/45 degrees (example same slices: equal to (Nslices-1)*135 degree 2D case)
    result = glcm3(levels, offset, arr)
    print result
    print
    offset = np.array([1,-1,-1], dtype = np.int32) # 135 degrees/135 degrees (example same slices: equal to (Nslices-1)*transpose of 135 degree 2D case)
    result = glcm3(levels, offset, arr)
    print result
    print
    
    
    
Created on Wed Feb 24 17:25:49 2016

@author: windows
"""


def glcm3d(levels, img, offsets): 
    ''' test:
    import numpy as np
    levels=3
    img = np.array([ [[1,2,2,0,1],[0,0,1,2,1],[1,1,0,0,2],[1,1,1,2,2]],
                     [[2,2,0,0,1],[0,1,2,2,1],[1,0,0,0,2],[1,1,2,2,0]]], dtype=np.uint8)
                    
    # img.shape = 2 x 4 x 5
    offsets = np.array([0,0,1], dtype = np.int32)
    result = np.zeros((levels,levels), dtype = np.int32)
    
    '''
    import numpy as np
    sizex=img.shape[0]
    sizey=img.shape[1]
    sizez=img.shape[2]
    arr = img.flatten()
    result = np.zeros((levels,levels), dtype = np.int32)
    arrresult = result.flatten()
    
    for k in range(sizez):
        for i in range(sizex):
            for j in range(sizey):
                aux = arr[k*sizex*sizey + i*sizey+j]
                slices = k+offsets[0]
                rows = i+offsets[1]
                cols = j+offsets[2]
                if(slices >= 0 and slices < sizez):
                    if(rows >= 0 and rows < sizex):
                        if(cols >= 0 and cols < sizey):
                            aux2 = arr[slices*sizex*sizey + rows*sizey+cols]
                            if(aux >= 0 and aux < levels):
                                if(aux2 >= 0 and aux2 < levels):
                                    arrresult[aux*levels+aux2]+=1
        
    result = arrresult.reshape(result.shape)
    return result
    

def glcmdesc(f,offset=[],mask=[]):
    '''
    import numpy as np
    f = np.array([ [[1,2,2,0,1],[0,0,1,2,1],[1,1,0,0,2],[1,1,1,2,2]],
                     [[2,2,0,0,1],[0,1,2,2,1],[1,0,0,0,2],[1,1,2,2,0]]], dtype=np.uint8)
                    
    # img.shape = 2 x 4 x 5
    offset = np.array([0,0,1], dtype = np.int32)
    '''

    import numpy as np
    from graycomatrix3D import glcm3d
    
    lev = int(f.max()+1) # levels
    g = glcm3d(lev,f,offset) 
    
    ### glcm normalization ###
    if g.sum() != 0:
        g = g.astype(float)/g.sum()
    
    ### compute auxiliary variables ###
    (num_level, num_level2) = g.shape
    I, J = np.ogrid[0:num_level, 0:num_level]
    I = 1+ np.array(range(num_level)).reshape((num_level, 1))
    J = 1+ np.array(range(num_level)).reshape((1, num_level))
    diff_i = I - np.apply_over_axes(np.sum, (I * g), axes=(0, 1))[0, 0]
    diff_j = J - np.apply_over_axes(np.sum, (J * g), axes=(0, 1))[0, 0]
    std_i = np.sqrt(np.apply_over_axes(np.sum, (g * (diff_i) ** 2),axes=(0, 1))[0, 0])
    std_j = np.sqrt(np.apply_over_axes(np.sum, (g * (diff_j) ** 2),axes=(0, 1))[0, 0])
    cov = np.apply_over_axes(np.sum, (g * (diff_i * diff_j)),axes=(0, 1))[0, 0]
    
    gxy = np.zeros(2*g.shape[0]-1)   ### g x+y
    gx_y = np.zeros(g.shape[0])  ### g x-y       
    for i in xrange(g.shape[0]):
        for j in xrange(g.shape[0]):
            gxy[i+j] += g[i,j]
            gx_y[np.abs(i-j)] += g[i,j]  
    mx_y = (gx_y*np.arange(len(gx_y))).sum()
    v = np.zeros(11)
    i,j = np.indices(g.shape)+1
    ii = np.arange(len(gxy))+2
    ii_ = np.arange(len(gx_y))
    

    ### compute descriptors ###
    v[0] = np.apply_over_axes(np.sum, (g ** 2), axes=(0, 1))[0, 0] # Angular second moment
    v[1] = np.apply_over_axes(np.sum, (g * ((I - J) ** 2)), axes=(0, 1))[0, 0] # Contrast
    if std_i>1e-15 and std_j>1e-15: # handle the special case of standard deviations near zero
        v[2] = cov/(std_i*std_j)#v[2] = greycoprops(g,'correlation') # Correlation
    else:
        v[2] = 1
    v[3] = np.apply_over_axes(np.sum, (g* (diff_i) ** 2),axes=(0, 1))[0, 0]# Sum of squares
    v[4] = np.sum(g * (1. / (1. + (I - J) ** 2))) # Inverse difference moment
    v[5] = (gxy*ii).sum() # Sum average
    v[6] = ((ii-v[5])*(ii-v[5])*gxy).sum() # Sum variance
    v[7] = -1*(gxy*np.log10(gxy+ np.spacing(1))).sum() # Sum entropy
    v[8] = -1*(g*np.log10(g+np.spacing(1))).sum() # Entropy
    v[9] = ((ii_-mx_y)*(ii_-mx_y)*gx_y).sum() # Difference variance
    v[10] = -1*(gx_y*np.log10(gx_y++np.spacing(1))).sum() # Difference entropy
    
    return g,v
 


