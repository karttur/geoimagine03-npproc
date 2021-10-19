'''
Created on 28 Mar 2021

@author: thomasgumbricht
'''

import os

from sys import exit

import numpy as np

import gc

#from numpy import pi, log, tan, empty, float32, arctan, rad2deg, gradient

from scipy.ndimage import gaussian_gradient_magnitude

from math import sqrt

from geoimagine.gis import RasterOpenGetFirstLayer

class Filter:
    '''
    '''
    
    def __init__(self):
        '''
        '''
        pass
    
    def _MovingWindow(self, radius, smooth=True):
            
        # ----------  create the moving window  ------------
        #r = self.pp.process.parameters.pixelradius
        
        if smooth:
            
            win = np.ones((2* radius +1, 2* radius +1))
            
        else:
            
            win = []
            
            for col in range (-radius,radius+1):
                
                colL = []
                
                for row in range  (-radius,radius+1):
                    
                    if abs(col) == radius and abs(row) == radius:
                        
                        colL.append(1)
                        
                    elif abs(col) == radius and row == 0:
                        
                        colL.append(1)
                        
                    elif abs(row) == radius and col == 0:
                        
                        colL.append(1)
                        
                    else:
                        
                        colL.append(0)
                    
                win.append(colL)
                
            win = np.asarray(win)
        
        '''
        kernelWin = []
        for col in range (-radius,radius+1):
            colL = []
            for row in range  (-radius,radius+1):
                if abs(col) == radius and abs(row) == radius:
                    
                colL.append([abs(col),abs(row)])
            kernelWin.append(colL)
            
            
        print ('kernelWin', kernelWin)   
        
        kernelWin = np.asanyarray(kernelWin)
        '''
                     
        distWin = []
        for col in range (-radius,radius+1):
            colL = []
            for row in range  (-radius,radius+1):
                colL.append(sqrt(col**2+row**2))
            distWin.append(colL)
            
        distWin = np.asanyarray(distWin)
        
        # ----------   or, copy paste your window matrix -------------
        # win = np.array( [    [0, 1, 1, 1, 0]
        #                      [1, 1, 1, 1, 1],
        #                      [1, 1, 0, 1, 1],
        #                      [1, 1, 1, 1, 1],
        #                      [0, 1, 1, 1, 0]  ])
        

        if self.pp.process.parameters.kernelform == 'square':
             
            pass # = default
        
        elif self.pp.process.parameters.kernelform == 'round':
            
            # Only keep valid values where distWin <= radius
            SNULLE
            # If not smooth = extract
            # Only retain the central smoothing cell for each kernel
    
        # window radius is needed for the function,
        # deduce from window size (can be different for height and width…)
        r_y, r_x  = win.shape[0]//2, win.shape[1]//2
        
        if self.pp.process.parameters.mode in ['TPI','TRI'] and not smooth:
            
            win[r_y, r_x  ]=0  # set the central cell to weight = 0 
            
            
        print('    Kernel:', np.array2string(win, prefix='    Kernel: '))
                
        return win, r_x, r_y
        
    def _View (self, offset_y, offset_x, shape, step=1):
        """
        Function returning two matching numpy views for moving window routines.
        - 'offset_y' and 'offset_x' refer to the shift in relation to the analysed (central) cell 
        - 'shape' are 2 dimensions of the data matrix (not of the window!)
        - 'view_in' is the shifted view and 'view_out' is the position of central cells
        (see on LandscapeArchaeology.org/2018/numpy-loops/)
        """
        size_y, size_x = shape
        x, y = abs(offset_x), abs(offset_y)
        
        x_in = slice(x , size_x, step) 
        x_out = slice(0, size_x - x, step)
    
        y_in = slice(y, size_y, step)
        y_out = slice(0, size_y - y, step)
     
        # the swapping trick    
        if offset_x < 0: x_in, x_out = x_out, x_in                                 
        if offset_y < 0: y_in, y_out = y_out, y_in
     
        # return window view (in) and main view (out)
        return np.s_[y_in, x_in], np.s_[y_out, x_out]

class ProcessNumpy(Filter):
    '''
    '''
    
    def __init__(self, pp, session):
        '''
        '''
        
        self.session = session
                
        self.pp = pp  
        
        self.verbose = self.pp.process.verbose 
        
        self.session._SetVerbosity(self.verbose) 
        
        Filter.__init__(self)      
        
        # Direct to subprocess
        if self.pp.process.processid == 'NumpyDemRegion':
            
            self._DEMRegion()
            
        elif self.pp.process.processid == 'NumpyGeomorphology':
            
            self._Geomorphology()
            
        elif self.pp.process.processid == 'NumpyDemTiles':
            
            self._DEMTiles()
            
        else:
        
            print (self.pp.process.processid)

            exit('Numpy process not recognized in processNumpy')
                   
    def _SetFileName(self,mode,dstLayerFPN,radius,radiusL):
        
        FP, FN = os.path.split(dstLayerFPN)
        
        baseN, ext = os.path.splitext(FN)
          
        if mode in ['TPI','TRI','roughness']:
                    
            radiusAddon = 'np-%sx%sx%s' %(radius,3,3)
                
        elif mode in ['slope','elev']:
            
            radiusAddon = 'np-%sx%s' %(radius*2+1,radius*2+1)
            
        elif mode[0:11] == 'landformTPI':
            
            if self.pp.process.parameters.standardize:
        
                radiusAddon = 'np-stnd-%s+%s' %(radiusL[0],radiusL[1])
                
            else:
                
                radiusAddon = 'np-raw-%s+%s' %(radiusL[0],radiusL[1])
                
        else:
            
            exitstr = "EXITING, unrecognized mode in NumpyProcess: %s" %(mode)
            
            exit(exitstr)
            
        radiusFN = '%s-%s%s' %(baseN,radiusAddon,ext)

        return os.path.join(FP,radiusFN)
            
    def _DstFPN(self,dstLayer,radiusL):
        '''
        '''
        dstDict = {}
        
        for radius in radiusL:
            
            radiusFPN = self._SetFileName(self.pp.process.parameters.mode,dstLayer.FPN,radius,radiusL)
            
            if os.path.exists(radiusFPN) and not self.pp.process.overwrite:
                
                if self.verbose:
                    
                    infostr = '    Destination layer already exists:\n        %s' %(radiusFPN)
                    
                    print (infostr)
                
                continue
            
            dstDict[radius] = radiusFPN
        
        return dstDict

    def _DEMTiles(self):
        '''
        '''
        
        for srcLocus in self.pp.srcLayerD:
            
            if self.verbose:
                
                infostr = '    Locus: %s' %(srcLocus)
                
                print (infostr)
            
            radiusL = self.pp.process.parameters.radiuscsv.split(',')
                    
            radiusL = [int(r) for r in radiusL]
            
            dstLayer = self._LoopSingleStaticDstLayer(srcLocus)
                
            dstDict = self._DstFPN(dstLayer,radiusL)
            
            if len(dstDict) == 0:
                
                continue
                        
            srcLayer = self._LoopSingleStaticSrcLayer(srcLocus)
            
            if srcLayer is None:
                
                continue
                
            srcLayer.RasterOpenGetFirstLayer()
                
            srcLayer.layer.ReadBand()
            
            for radius in dstDict:
                
                if self.pp.process.parameters.mode == 'elev':
                    
                    dstArr = self._SmoothBand(srcLayer.layer.NPBAND,radius)
                            
                elif self.pp.process.parameters.mode == 'TPI':
                    
                    dstArr = self._TPIfromDEM(srcLayer.layer.NPBAND,radius)
                                                        
                elif self.pp.process.parameters.mode == 'slope':
                    
                    dstArr = self._SlopefromDEM(srcLayer.layer.NPBAND,radius)
                    
                elif self.pp.process.parameters.mode == 'landformTPI':
                    
                    dstArr = self._IniTPILFfromDEM(srcLayer, dstLayer.FPN, dstDict)
                     
                else:
                    
                    exitstr = 'EXITING - unrecoognized mode in NumpyProcess._DEMregion: %s' %(self.pp.process.parameters.mode)
                    
                    exit(exitstr)
                
                if dstArr is None:
                    
                    continue
                
                #Create the dst layer
                #dstLayer.layer = lambda:None
                dstLayer.EmptyLayer()
                    
                #copy the geoformat from the src layer
                self._SetDstMeta(srcLayer, dstLayer)
                           
                #Set the np array as the band
                dstLayer.layer.NPBAND = dstArr
                
                #Reset the dstlayer name
                dstLayer.FPN = dstDict[radius]
                                
                dstLayer.RasterCreateWithFirstLayer()
                                    
                if self.verbose:
                    
                    infostr = '    Destination layer created:\n        %s' %(dstLayer.FPN)
                        
                    print (infostr)
                    
                if self.pp.process.parameters.mode == 'LandformTPI':
                    
                    break
                
                srcLayer.layer.NPBAND = None
                
                dstLayer.layer = None
                
                dstLayer = None
            
                              
            srcLayer.layer.NPBAND = None
                
            srcLayer.layer = None
                
            srcLayer = None
            
            # Collect the garbage
            gc.collect()   
    
    def _LoopSingleStaticSrcLayer(self,srcLocus):
        '''
        '''
                               
        for datum in self.pp.srcLayerD[srcLocus]:
                      
            for comp in self.pp.srcLayerD[srcLocus][datum]:
                                    
                if not os.path.exists(self.pp.srcLayerD[srcLocus][datum][comp].FPN):
                
                    infostr = '    SKIPPING - region layer for tiling missing\n        %s' %(self.pp.srcLayerD[srcLocus][datum][comp].FPN)
        
                    print (infostr)
                    
                    return None
                    
                return self.pp.srcLayerD[srcLocus][datum][comp]
            
    def _LoopSingleStaticDstLayer(self,dstLocus):
        '''
        '''
                               
        for datum in self.pp.dstLayerD[dstLocus]:
                      
            for comp in self.pp.dstLayerD[dstLocus][datum]:
                                                        
                return self.pp.dstLayerD[dstLocus][datum][comp]
            
    def _LoopMultiStaticSrcLayers(self,srcLocus):
        '''
        '''
        self.srcLayerD  = {}  
                       
        for datum in self.pp.srcLayerD[srcLocus]:
                      
            for comp in self.pp.srcLayerD[srcLocus][datum]:
                     
                self.srcLayerD[comp] =  self.pp.srcLayerD[srcLocus][datum][comp]                                  
                            
    def _SlopefromDEM(self, BAND,radius):
        '''
        '''

        infostr = '    Calculating slope with numpy/scipy, radius = %s' %(radius)
        
        print (infostr)
        
        slope = gaussian_gradient_magnitude(BAND, sigma=radius, mode='nearest')
        
        if self.pp.process.parameters.trigonometric:
            
            slope = np.rad2deg(np.arctan(slope / 100))

        return slope
                     
    def _TPIfromDEM(self, BAND, radius):
        '''
        '''
        
        # http://bepls.com/october2014bepls/6.pdf
        
        # see https://landscapearchaeology.org/2019/tpi/
        # Open a DEM and let a kernel slide over to det 
        #- slope 
        #- TPI at high resol
        #- TPI at average resol
        #- Derive landform directly
        
        if radius > 1:
            
            if (radius-1) % 2 != 0:
                
                exitstr = 'EXITING: radius = %s can not form a half size kernel' %(radius)
                
                exit(exitstr)
            
            smoothRadius = int((radius-1)/2)
            
            BAND = self._SmoothBand(BAND,smoothRadius)
            
        if self.verbose:
            
            if radius > 1:
            
                infostr = '    Calculating TPI with smoothed radius: %s (kernel: %sx%s)' %(radius, radius*2+1, radius*2+1)
            
            else:
                
                infostr = '    Calculating TPI with radius: %s (kernel: %sx%s)' %(radius, radius*2+1,radius*2+1)
   
            print (infostr)
        
        win, r_x, r_y = self._MovingWindow(radius, False)
                
        mx_temp = np.zeros(BAND.shape)
        
        mx_count = np.zeros(BAND.shape)
          
        # loop through window and accumulate values
        for (y,x), weight in np.ndenumerate(win):
            
            if weight == 0 : continue  #skip zero values !
            # determine views to extract data 
            
            view_in, view_out = self._View(y - r_y, x - r_x, BAND.shape)
            #print ('view_in',view_in)
            #print ('view_out',view_out)
            # using window weights (eg. for a Gaussian function)
            mx_temp[view_out] += BAND[view_in]  * weight
            
            # track the number of neighbours 
            # (this is used for weighted mean : Σ weights*val / Σ weights)
            mx_count[view_out] += weight
            
        # this is TPI (spot height – average neighbourhood height)
        TPI = BAND - mx_temp / mx_count
        
        mx_temp = None
        
        mx_count = None
        
        return (TPI)
   
    def _SmoothBand(self, BAND, radius):
        ''' Smooth DEM by averaging while keeping original resolution
        '''        
        
        if self.verbose:
            
            infostr = '    Smoothing DEM with radius %s' %(radius)
            
            print (infostr)
        
        win, r_x, r_y = self._MovingWindow(radius, True)
                
        mx_temp = np.zeros(BAND.shape)
        
        mx_count = np.zeros(BAND.shape)
           
        # loop through window and accumulate values
        for (y,x), weight in np.ndenumerate(win):
            
            if weight == 0 : continue  #skip zero values !
            
            # determine views to extract data 
            view_in, view_out = self._View(y - r_y, x - r_x, BAND.shape)

            # using window weights (eg. for a Gaussian function)
            mx_temp[view_out] += BAND[view_in]  * weight
            
            # track the number of neighbours 
            # (this is used for weighted mean : Σ weights*val / Σ weights)
            mx_count[view_out] += weight

        # return the smoothed average
        
        mx_arr = mx_temp / mx_count
        
        # Delete and garbage collect
        
        del mx_temp
        
        del mx_count
        
        gc.collect()
         
        return mx_arr

    def _IniTPILFfromDEM(self, srcLayer, dstLayerFPN, dstDict):
        ''' TPI landform from DEM 
        '''
        
        '''
        # TPI landform from DEM requires 3 input layers  
        # - detailed TPI
        # - blunt TPI
        # - slope
        '''
        
        inputD = {'dTPI':{'FPN':False},'bTPI':{'FPN':False},'slope':{'FPN':False}}
        
        # detailed TPI
        radiusL = self.pp.process.parameters.radiuscsv.split(',')
                    
        radiusL = [int(r) for r in radiusL]
            
        # Reconstruct the name of the input TPI layer  
        FP,FN = os.path.split(dstLayerFPN)
        
        FNparts = FN.split('_')
        
        TPIFN = 'tpi_%s_%s_%s_%s' %(FNparts[1],FNparts[2],FNparts[3],FNparts[4])
        
        TPIFPN = os.path.join(FP,TPIFN)
        
        detailedTPIFPN = self._SetFileName('TPI',TPIFPN,radiusL[0],radiusL)
                
        inputD['dTPI']['radius'] = radiusL[0]
        
        if os.path.exists(detailedTPIFPN):
            
            inputD['dTPI']['FPN'] = detailedTPIFPN
        
        bluntTPIFPN = self._SetFileName('TPI',TPIFPN,radiusL[1],radiusL)
        
        inputD['bTPI']['radius'] = radiusL[1]
        
        if os.path.exists(bluntTPIFPN):
            
            inputD['bTPI']['FPN'] = bluntTPIFPN
        
        base,ext = os.path.splitext(FNparts[4])
        
        slopeFN = 'slope_%s_%s_%s_%s-%s%s' %(FNparts[1],FNparts[2],FNparts[3], base, self.pp.process.parameters.slope, ext)
        
        slopeFPN = os.path.join(FP,slopeFN)
        
        inputD['slope']['radius'] = int(self.pp.process.parameters.sloperadius)
        
        if os.path.exists(slopeFPN):
            
            inputD['slope']['FPN'] = slopeFPN
        
        # Loop to see which layers must be created on the fly
        
        for comp in inputD:
            
            if not inputD[comp]['FPN']:
                
                inputD[comp]['layer'] = lambda:None
                
                if comp == 'slope':
                    
                    inputD[comp]['cellnull'] = srcLayer.layer.cellnull
                    
                    inputD[comp]['layer'].NPBAND = self._SlopefromDEM(srcLayer.layer.NPBAND,inputD[comp]['radius'])
                    
                else:
                    
                    inputD[comp]['cellnull'] = srcLayer.layer.cellnull
                    
                    inputD[comp]['layer'].NPBAND = self._TPIfromDEM(srcLayer.layer.NPBAND,inputD[comp]['radius'])
        
        # Read the preexisting bands
        for comp in inputD:
            
            if inputD[comp]['FPN']:
                
                inputD[comp]['DS'], inputD[comp]['layer'] = RasterOpenGetFirstLayer(inputD[comp]['FPN'],{'mode':'read'})
                
                inputD[comp]['layer'].ReadBand()
                                
                inputD[comp]['cellnull'] = inputD[comp]['layer'].cellnull
                
        #Create the dst array
        dstArr = np.empty_like(inputD[comp]['layer'].NPBAND, dtype=np.uint8)
          
        if self.pp.process.parameters.standardize:
            
            #Standardize the TPI prior to Landform classification
            if self.pp.process.parameters.dtpistd and self.pp.process.parameters.btpistd:
                
                self._LandformGlobalStandardized(inputD['dTPI']['layer'].NPBAND,
                                inputD['bTPI']['layer'].NPBAND,
                                inputD['slope']['layer'].NPBAND,
                                dstArr)
            
            elif self.pp.process.parameters.dtpistd and not self.pp.process.parameters.btpistd:
            
                print ('    Relative standardized')

                self._LandformRelativeStandardized(inputD['dTPI']['layer'].NPBAND,
                                inputD['bTPI']['layer'].NPBAND,
                                inputD['slope']['layer'].NPBAND,
                                dstArr)
            
            else:
                print ('    Local standardized')

                #Calcualte local std for TPI (large and small)
                self._LandformStandardized(inputD['dTPI']['layer'].NPBAND,
                                inputD['bTPI']['layer'].NPBAND,
                                inputD['slope']['layer'].NPBAND,
                                dstArr)
      
        else:
            
            #Landform classification against original TPI data
            
            print ('    Raw TPI classification')
            
            self._LandformTPI(inputD['dTPI']['layer'].NPBAND,
                                inputD['bTPI']['layer'].NPBAND,
                                inputD['slope']['layer'].NPBAND,
                                dstArr)
                


        dstArr[ np.where( inputD['dTPI']['layer'].NPBAND == inputD['dTPI']['cellnull'] ) ] = 255
        dstArr[np.where( inputD['bTPI']['layer'].NPBAND == inputD['bTPI']['cellnull'] ) ] = 255
        dstArr[np.where( inputD['slope']['layer'].NPBAND == inputD['slope']['cellnull'] ) ] = 255

        return dstArr
    
    
    def _LandformTPI(self,dTPI,bTPI,SLOPE,dstArr):
        '''
        '''   
        # http://bepls.com/october2014bepls/6.pdf
        
        # see https://landscapearchaeology.org/2019/tpi/
           
        tpiTH = self.pp.process.parameters.tpithreshold
        
        slopeTH = self.pp.process.parameters.slopethreshold
        
        #Plain = 50 (alternative = 10)
        dstArr[np.where( (dTPI > -tpiTH) & (dTPI < tpiTH) & (bTPI > -tpiTH) & (bTPI < tpiTH) & (SLOPE <= slopeTH))] = 50
        
        #Open slope = 6 (alternative = 11)
        dstArr[np.where( (dTPI > -tpiTH) & (dTPI < tpiTH) & (bTPI > -tpiTH) & (bTPI < tpiTH) & (SLOPE > slopeTH))] = 60


        # Mesa or flat ridge = 7 (alternative = 12)
        #upper slope edge = 71
        dstArr[np.where( (dTPI > -tpiTH) & (dTPI < tpiTH) & (bTPI >= tpiTH) & (SLOPE > slopeTH))] = 71
        #mesa or flat ridge top = 72
        dstArr[np.where( (dTPI > -tpiTH) & (dTPI < tpiTH) & (bTPI >= tpiTH) & (SLOPE <= slopeTH))] = 72
        
        
        #U-shaped valley = 4  (alternative = 9)      
        #valley floor edge = 41
        dstArr[np.where( (dTPI > -tpiTH) & (dTPI < tpiTH) & (bTPI <= -tpiTH) & (SLOPE > slopeTH))] = 41
        #central valley floor  = 42
        dstArr[np.where( (dTPI > -tpiTH) & (dTPI < tpiTH) & (bTPI <= -tpiTH) & (SLOPE <= slopeTH))] = 42
        
        # shallow valley / midslope drainage = 2 (alternative = 6)
        #shallow valley edge = 21
        dstArr[np.where( (dTPI <= -tpiTH) & (bTPI > -tpiTH) & (bTPI < tpiTH) & (SLOPE > slopeTH))] = 21
        #central shallow valley = 21
        dstArr[np.where( (dTPI <= -tpiTH) & (bTPI > -tpiTH) & (bTPI < tpiTH) & (SLOPE <= slopeTH))] = 22
        
        #midslope ridge / hills in valleys = 9 (alterantive = 14)
        #midslope ridge / hills in valleys = 91
        dstArr[np.where( (dTPI >= tpiTH) & (bTPI > -tpiTH) & (bTPI < tpiTH) & (SLOPE > slopeTH))] = 91
        #midslope ridge / hills in valleys = 92
        dstArr[np.where( (dTPI >= tpiTH) & (bTPI > -tpiTH) & (bTPI < tpiTH) & (SLOPE <= slopeTH))] = 92

        #upland drainage = 3 (alternative = 8)
        #upland drainage slope = 31
        dstArr[np.where( (dTPI <= -tpiTH) & (bTPI >= tpiTH) & (SLOPE > slopeTH))] = 31
        #upland drainage flat = 32
        dstArr[np.where( (dTPI <= -tpiTH) & (bTPI >= tpiTH) & (SLOPE <= slopeTH))] = 32
        
        #canyon / incised stream = 1 (alternative = 5)
        #incised canyon - rapid = 11
        dstArr[np.where( (dTPI <= -tpiTH) & (bTPI <= -tpiTH) & (SLOPE > slopeTH))] = 11
        #incised canyon - calm = 12
        dstArr[np.where( (dTPI <= -tpiTH) & (bTPI <= -tpiTH) & (SLOPE <= slopeTH))] = 12
      
        #peak / mt top = 10 (alterantive = 16)
        #peak / mt top = 101
        dstArr[np.where( (dTPI >= tpiTH) & (bTPI >= tpiTH) & (SLOPE > slopeTH))] = 101
        #peak / mt top = 102
        dstArr[np.where( (dTPI >= tpiTH) & (bTPI >= tpiTH) & (SLOPE <= slopeTH))] = 102
        
        #local ridge = 8 (alterantive = 13)
        #local ridge = 81
        dstArr[np.where( (dTPI >= tpiTH) & (bTPI < -tpiTH) & (SLOPE > slopeTH))] = 81
        #local ridge = 82
        dstArr[np.where( (dTPI >= tpiTH) & (bTPI < -tpiTH) & (SLOPE <= slopeTH))] = 82 
    
    def _Geomorphology(self):
        '''
        '''
        
        for srcLocus in self.pp.srcLayerD:
            
            dstLayer = self._LoopSingleStaticDstLayer(srcLocus)
            
            if dstLayer._Exists() and not self.pp.process.overwrite:
                
                if self.verbose:
                    
                    infostr = '    Destination layer already exists:\n        %s' %(dstLayer.FPN)
                    
                    print (infostr)
                
                continue
                        
            self._LoopMultiStaticSrcLayers(srcLocus)
            
            for layer in self.srcLayerD:
                
                #print (self.srcLayerD[layer].FPN)
                
                self.srcLayerD[layer].RasterOpenGetFirstLayer()
                
                self.srcLayerD[layer].layer.ReadBand()
            
            #Create the dst array
            dstArr = np.empty_like(self.srcLayerD[layer].layer.NPBAND, dtype=np.uint8)
            
            if self.pp.process.parameters.mode == 'weiss':
 
                self._WeissbTPIF(dstArr)
                
            else:
                
                exitstr = 'EXITING - Geomorphology model %s does not exist (npproc)' %(self.pp.process.parameters.mode)
                
                exit(exitstr)
                
            #Create the dst layer            
            dstLayer.EmptyLayer()
              
            # copy the geoformat from any src layer
            # MOVE TO RASTERLAYER FUNCTION
            self._SetDstMeta(self.srcLayerD[layer], dstLayer)
            
            # Close the srclayers

            for layer in self.srcLayerD:
            
                self.srcLayerD[layer].layer.NPBAND = None
                
                self.srcLayerD[layer].layer = None
                
                self.srcLayerD[layer] = None
                
            # collect the garbage
            gc.collect()
          
            #Set the np array as the band
            dstLayer.layer.NPBAND = dstArr
            
            # For some reason I must send the layer along
            dstLayer.RasterCreateWithFirstLayer()
                                            
            if self.verbose:
                
                infostr = '    Destination layer created:\n        %s' %(dstLayer.FPN)
                    
                print (infostr)
                
    def _LandformGlobalStandardized(self,dTPI,bTPI,SLOPE,dstArr):
        '''
        '''
        dtpistd = self.pp.process.parameters.dtpistd
        
        btpistd = self.pp.process.parameters.btpistd
        
        #dTPIf = dTPI.astype(np.float32)

        dTPIstandard = (dTPI/dtpistd)*100+0.5
        
        #bTPIf = bTPI.astype(np.float32)

        bTPIstandard = (bTPI/btpistd)*100+0.5
        
        self._LandformTPI(dTPIstandard,bTPIstandard,SLOPE,dstArr)
                   
                        
    def _SetDstMeta(self, srcLayer, dstLayer):
        '''
        '''
        #Transfer the geoformat for this locus
        
        l = srcLayer.layer
        
        geoFormatD = {'lins':l.lins,'cols':l.cols,'projection':l.projection,'geotrans':l.geotrans,'cellsize':l.cellsize}
                                      
        for item in geoFormatD:
                    
            setattr(dstLayer.layer, item, geoFormatD[item]) 