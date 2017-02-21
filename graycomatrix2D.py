# -*- coding: utf-8 -*-
"""
================
graycomatrix
================

C functions
==============

-`glcm2d glcm2d`


2D offsets
=============

- offset = [0 D] -> zero degrees, distance D
- offset = [D 0] -> 90 degrees, distance D
- offset = [D D] -> 135 degrees, distance D
- offset = [D -D] -> 45 degrees, distance D 


Function code
================

.. code:: python module

   from iatexture.glcm2d import glcm2d_c
   import numpy as np
   
   def glcm(levels,offsets,img):
      result = np.zeros((levels,levels), dtype = np.int32)
      if len(img.shape) ==2:
          if offsets == []:
              offsets = np.array([0,1], dtype = np.int32)
          glcm2d_c(levels, img, offsets, result)
      elif len(img.shape) ==3:
          if offsets == []:
              offsets = np.array([0,0,1], dtype = np.int32)
          glcm3d_c(levels, img, offsets, result)
      else:
          print 'Invalid input. Please use another array'
      return result
      
    #include <stdio.h>
    #include <string.h>
    #include <math.h> 
    
    void glcm2d_c(int levels, int sizex, int sizey, int *arr, int h, int *offset, int h1, int w1, int *result){
       
        int i,j,k,aux,aux2,rows,cols;
    
        for(i=0;i<sizex;i++){
            k = i*sizey;
            for(j=0;j<sizey;j++){
                aux = arr[k+j];
                rows = i+offset[0];
                cols = j+offset[1];
                if(rows >= 0 && rows < sizex){
                    if(cols >= 0 && cols < sizey){
                        aux2 = arr[rows*sizey+cols];
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
    
    ##########
    def glcm2d(levels, img, offsets, result):
        ''' test:
        import numpy as np
        levels=3
        img = np.array([[1,2,2,0,1],[0,0,1,2,1],[1,1,0,0,2],[1,1,1,2,2]], dtype = np.uint8)
        # img.shape = 4 x 5
        offsets = np.array([0,1], dtype = np.int32)
        result = np.zeros((levels,levels), dtype = np.int32)
        
        '''
        sizex=img.shape[0]
        sizey=img.shape[1]
        arr = img.flatten()
        arrresult = result.flatten()
        
        for i in range(sizex):
            k = i*sizey
            for j in range(sizey):
                aux = arr[k+j]
                print(aux)
                rows = i+offsets[0]
                print(rows)
                cols = j+offsets[1]
                print(cols)
                if(rows >= 0 and rows < sizex):
                    if(cols >= 0 and cols < sizey):
                        aux2 = arr[rows*sizey+cols]
                        if(aux >= 0 and aux < levels):
                            if(aux2 >= 0 and aux2 < levels):
                                arrresult[aux*levels+aux2]+=1
            
        result = arrresult.reshape(result.shape)
        
    return result


Examples
============

2D examples
--------------
      
.. code:: python 

    import numpy as np
    from iatexture.graycomatrix import glcm
    
    #levels = 3
    arr = np.array([[1,2,2,0,0,1],[0,0,1,2,2,1],[1,1,0,0,0,2],[1,1,1,2,2,2],[1,1,2,2,0,0]], dtype = np.uint8)
    levels = int(arr.max()+1)
    
.. code:: python 
  
    offset = np.array([0,1], dtype = np.int32) # orientation 0 degrees
    result = glcm(levels, offset, arr)
    print result

    offset = np.array([1,-1], dtype = np.int32) # orientation 45 degrees
    result = glcm(levels, offset, arr)
    print result

    offset = np.array([1,0], dtype = np.int32) # orientation 90 degrees
    result = glcm(levels, offset, arr)
    print result

    offset = np.array([1,1], dtype = np.int32) # orientation 135 degrees
    result = glcm(levels, offset, arr)
    print result



Equation
========

.. equation:: latex
    :idpi: 100

The co-occurence matrix becomes the estimate of the joint probability, :eq:`p_{d\theta}`, of 
two pixels, a distance :eq:`d` apart along a given direction :eq:`\theta` having particular 
(co-occurring) values :eq:`i` and :eq:`j`.

Let :eq:`\mu_{x}`, :eq:`\mu_{y}` e :eq:`\rho_{x}`, :eq:`\rho_{y}` denote the mean and standard 
deviations of the row and column sums of the co-occurrence matrix, respectively [related to the 
marginal distributions :eq:`p_{x}(i)` and :eq:`p_{y}(j)`], and :eq:`N_{g}` the number of 
discrete intensity levels.

Angular second moment
-----------------------

.. equation:: latex

    AngScMom = \sum_{i=1}^{N_{g}} \sum_{i=1}^{N_{g}}p(i,j)^{2} 

Contrast
---------

.. equation:: latex

    Contrast = \sum_{n=0}^{N_{g}-1} n^{2} \sum_{i=1}^{N_{g}} \sum_{i=j}^{N_{g}}p(i,j) \\ where \ n = \left | i-j \right |  

Correlation
------------

.. equation:: latex
    
    Correlation = \frac{\sum_{i=1}^{N_{g}} \sum_{i=1}^{N_{g}}ijp(i,j) -\mu_{x}\mu_{y}}{\rho_{x}\rho_{y}}

Sum of squares
---------------

.. equation:: latex

    SumOfSqs = \sum_{i=1}^{N_{g}} \sum_{i=1}^{N_{g}}(i -\mu_{x})^2 p(i,j)

Inverse difference moment
----------------------------

.. equation:: latex

    InvDfMom = \sum_{i=1}^{N_{g}} \sum_{i=1}^{N_{g}}\frac{1}{1+(i-j)^{2}}p(i,j)

Sum average
--------------

.. equation:: latex
    
    SumAverg = \sum_{i=1}^{2N_{g}}ip_{x+y}(i)) \\ where \ p_{x+y}(k) = \sum_{i=1}^{N_{g}}\sum_{j=1}^{N_{g}}p(i,j) \\ k = i+j

Sum variance
---------------

.. equation:: latex

    SumVarnc = \sum_{i=1}^{2N_{g}} (i-SumAverg)^{2}p_{x+y}(i)

Sum entropy
------------
    
.. equation:: latex
    
    SumEntrp = -\sum_{i=1}^{2N_{g}} p_{x+y}(i)log_{10}(p_{x+y}(i))

Entropy
---------
    
.. equation:: latex

    Entropy = -\sum_{i=1}^{N_{g}}\sum_{j=1}^{N_{g}} p(i,j)log_{10}(p(i,j))

Difference variance
----------------------
    
.. equation:: latex

    DifVarnc = \sum_{i=1}^{N_{g}-1}(i-\mu_{x-y})^2p_{x-y}(i)) \\ where \ p_{x-y}(k) = \sum_{i=1}^{N_{g}}\sum_{j=1}^{N_{g}}p(i,j)\\ k = \left | i-j \right |

Difference entropy
-------------------

.. equation:: latex    

    DifEntrp = -\sum_{i=1}^{N_{g}}p_{x-y}(i)log_{10}(p_{x-y}(i)) 
        
References
==========

- `http://www.mathworks.com/help/images/ref/graycomatrix.html Matlab GLCM`
- `http://www.fp.ucalgary.ca/mhallbey/tutorial.htm Tutorial about GLCM`

See Also
========

- `iatexture:graycomatrix Calculation of the gray co-occurence matrix`
- `http://scikit-image.org/docs/dev/api/skimage.feature.texture.html#skimage.feature.texture.greycomatrix scikits image`

Adapted from Contributions by:
==============
- Mariana Bento, october, 23th 2013: initial function.
- Mariana Leite, september, 9th 2014: C implementation and inclusion of 3D images.  
    
    
Created on Wed Feb 24 17:25:49 2016

@author: windows
"""
   
def glcm2d(levels, img, offsets):
    ''' test:
    import numpy as np
    levels=3
    img = np.array([[1,2,2,0,1],[0,0,1,2,1],[1,1,0,0,2],[1,1,1,2,2]], dtype = np.uint8)
    # img.shape = 4 x 5
    offsets = np.array([0,1], dtype = np.int32)
    result = np.zeros((levels,levels), dtype = np.int32)
    
    '''
    import numpy as np
    sizex=img.shape[0]
    sizey=img.shape[1]
    arr = img.flatten()
    result = np.zeros((levels,levels), dtype = np.int32)
    arrresult = result.flatten()
    
    for i in range(sizex):
        k = i*sizey
        for j in range(sizey):
            aux = arr[k+j]
            print(aux)
            rows = i+offsets[0]
            print(rows)
            cols = j+offsets[1]
            print(cols)
            if(rows >= 0 and rows < sizex):
                if(cols >= 0 and cols < sizey):
                    aux2 = arr[rows*sizey+cols]
                    if(aux >= 0 and aux < levels):
                        if(aux2 >= 0 and aux2 < levels):
                            arrresult[aux*levels+aux2]+=1
        
    result = arrresult.reshape(result.shape)
    return result


def glcm(levels, offsets, img):
    import numpy as np
    from graycomatrix2D import glcm2d
    # glcm will have dimensions based on number of greylevel levels
    if len(img.shape) == 2:
        if offsets == []:
            offsets = np.array([0,1], dtype = np.int32)
        # call main
        result = glcm2d(levels, img, offsets)
        print(result)
    
    return result



def glcmdesc(f,offsets=[],mask=[]):
    '''
    ## test: 2D GLCM = [[2,2,1],[1,3,3],[1,1,2]]
    import numpy as np
    f = np.array([[1,2,2,0,1],[0,0,1,2,1],[1,1,0,0,2],[1,1,1,2,2]], dtype = np.uint8)
    # f.shape = 4 x 5
    offsets = np.array([0,1], dtype = np.int32)
    mask=[]
    '''
    import numpy as np
    from graycomatrix2D import glcm
    
    lev = int(f.max()+1) # levels
    if mask != []: # if there is a mask, it is necessary to set a invalid value to all the pixels that are not in the mask (this invalid value can be any value that it is not present inside the mask)
        aux = f.copy()
        hist,_ = np.histogram(np.ravel(f[mask]),bins=lev,range=(0,lev-1)) # look for an empty level inside the mask
        if np.argwhere(hist==0).shape[0]!=0: # if there is an empty levels
            new_v = np.argwhere(hist==0).max() # it is selected the highest empty gray level 
        else: # # if there is not an empty levels
            new_v = lev # the chosen invalid value is f.max()+1
            lev+=1 # in this case, the number of gray levels increases
        
        aux[~mask] = new_v # all pixels that are not inside the mask receive this invalid value
        g = glcm(lev,offset,aux) # the glcm matrix is computed
        g[new_v,:] = 0 # the rows and columns that corresponds to the invalid value (not inside the mask) are nulled 
        g[:,new_v] = 0 # it is important to notice that null rows or columns do not change the descriptors
    else: # if there is no mask
        g = glcm(lev, offsets, f) 
    
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
 


