'''
Created on 11 Jan 2019

@author: thomasgumbricht
'''


import os
from sys import exit
#import geoimagine.gis.mj_gis_v80 as mj_gis 
#from geoimagine.gdalutilities import GDALstuff
import numpy as np
import gc

GDALpath = '/Library/Frameworks/GDAL.framework/Versions/2.1/Programs'

class ProcessDEM:
    '''class for DEM processing'''   
    def __init__(self, process, session, verbose):
        self.session = session
        self.verbose = verbose
        self.process = process
        #DEM processing is restricted to tiles only
        print ('        ProcessDEM',self.process.proc.processid) 
        if self.process.proc.processid[0:3].lower() in ['tpi','tri','rou']:
            self._LoopOneToMany() 
        elif self.process.proc.processid[0:8].lower() in ['landform']:
            self._LoopManyToOne()
        else:
            self._LoopOneToOne()
            
            
    def _LoopManyToOne(self):
        for locus in self.process.srcLayerD:
            for datum in self.process.srcLayerD[locus]:
                srcCompD = {}
                dstCompL = []
                #Loop the targets for this layer

                for dstComp in self.process.dstLayerD[locus][datum]:

                    if not self.process.dstLayerD[locus][datum][dstComp]._Exists() or self.process.overwrite:
                        dstCompL.append(dstComp)
                        print ('    Processing',self.process.dstLayerD[locus][datum][dstComp].FPN)
                    elif self.process.dstLayerD[locus][datum][dstComp]._Exists():
                        print ('    Done',self.process.dstLayerD[locus][datum][dstComp].FPN)
                        self.session._InsertLayer(self.process.dstLayerD[locus][datum][dstComp], self.process.overwrite, self.process.delete)

                if len(dstCompL) == 0:
                    print ('nothing to create')
                    continue

                dstComp = dstCompL[0]

                for srcComp in self.process.srcLayerD[locus][datum]:
                    #Check if the file is there
                    if os.path.exists(self.process.srcLayerD[locus][datum][srcComp].FPN):
                        compid = self.process.srcLayerD[locus][datum][srcComp].comp.id
                        srcCompD[compid] = srcComp
                        srcCompD[compid] = self.process.srcLayerD[locus][datum][srcComp]
                        srcCompD[compid].ReadRasterLayer()
                    else:
                        print ('missing',self.process.srcLayerD[locus][datum][srcComp].FPN)
                        SNULLEBULLE
                        break
                        
                if self.process.proc.processid[0:8].lower() == 'landform':
                    self._LandformTPIini(locus,datum,srcCompD,dstComp)
                
                else:
                    print (self.process.proc.processid)
                    SNULLE

                #Register the layer
                self.session._InsertLayer(self.process.dstLayerD[locus][datum][dstComp], self.process.overwrite, self.process.delete)
                self.process.dstLayerD[locus][datum][dstComp].layer.NPBAND = None
                self.process.dstLayerD[locus][datum][dstComp] = None
                for comp in srcCompD:
                    srcCompD[comp].layer.NPBAND = None
                    srcCompD[comp] = None
        
                #grabage collect
                gc.collect()
            
    def _LoopOneToMany(self):
        for locus in self.process.srcLayerD:
            for datum in self.process.srcLayerD[locus]:
                dstCompL = []
                for srcComp in self.process.srcLayerD[locus][datum]:
                    #Check if the file is there
                    if os.path.exists(self.process.srcLayerD[locus][datum][srcComp].FPN):
                        #Loop the targets for this layer
                        for dstComp in self.process.dstLayerD[locus][datum]:

                            if not self.process.dstLayerD[locus][datum][dstComp]._Exists() or self.process.overwrite:
                                dstCompL.append(dstComp)
                                print ('    Processing',self.process.dstLayerD[locus][datum][dstComp].FPN)
                            elif self.process.dstLayerD[locus][datum][dstComp]._Exists():
                                self.session._InsertLayer(self.process.dstLayerD[locus][datum][dstComp], self.process.overwrite, self.process.delete)

                        if len(dstCompL) == 0:
                            continue
                        #Get the raster metadata for the source layer
                        self.process.srcLayerD[locus][datum][srcComp].GetRastermetadata()
                        if self.process.proc.processid[0:3].lower() == 'tpi':
                            self._TPITRI(locus,datum,srcComp,dstCompL,'tpi')
                        elif self.process.proc.processid[0:3].lower() == 'tri':
                            self._TPITRI(locus,datum,srcComp,dstCompL,'tri')
                        elif self.process.proc.processid[0:9].lower() == 'roughness':
                            self._TPITRI(locus,datum,srcComp,dstCompL,'rn')
                        elif self.process.proc.processid[0:9].lower() == 'hillshade':
                            self._HillShade(locus,datum,srcComp,dstCompL)
                        else:
                            print (self.process.proc.processid)
                            SNULLE
                        for dstcomp in dstCompL:
                            self.session._InsertLayer(self.process.dstLayerD[locus][datum][dstcomp], self.process.overwrite, self.process.delete)

                            
    def _LoopOneToOne(self):
        for locus in self.process.srcLayerD:
            for datum in self.process.srcLayerD[locus]:
                dstCompL = []
                for srcComp in self.process.srcLayerD[locus][datum]:
                    #Check if src file is there
                    if os.path.exists(self.process.srcLayerD[locus][datum][srcComp].FPN):
                        #Loop the targets for this layer
                        for dstComp in self.process.dstLayerD[locus][datum]:
                            if not self.process.dstLayerD[locus][datum][dstComp]._Exists() or self.process.overwrite:
                                dstCompL.append(dstComp)
                                print ('    Processing',self.process.dstLayerD[locus][datum][dstComp].FPN)
                            elif self.process.dstLayerD[locus][datum][dstComp]._Exists():
                                self.session._InsertLayer(self.process.dstLayerD[locus][datum][dstComp], self.process.overwrite, self.process.delete)

                        if len(dstCompL) == 0:
                            continue
                        #Get the raster metadata for the source layer
                        self.process.srcLayerD[locus][datum][srcComp].GetRastermetadata()
                        if self.process.proc.processid[0:9].lower() == 'hillshade':
                            self._HillShade(locus,datum,srcComp,dstCompL)
                        elif self.process.proc.processid[0:5].lower() == 'slope':
                            self._Slope(locus,datum,srcComp,dstCompL)
                        elif self.process.proc.processid[0:6].lower() == 'aspect':
                            self._Aspect(locus,datum,srcComp,dstCompL)
                        else:
                            print (self.process.proc.processid)
                            SNULLE
                        for dstcomp in dstCompL:
                            self.session._InsertLayer(self.process.dstLayerD[locus][datum][dstcomp], self.process.overwrite, self.process.delete)


                            
    def _TPITRI(self,locus,datum,srcComp,dstCompL,txi):
        lins = self.process.srcLayerD[locus][datum][srcComp].comp.metadata.lins
        cols = self.process.srcLayerD[locus][datum][srcComp].comp.metadata.cols
        if self.process.params.mosaic:
            srcFPN = self._MosaicAdjacentTiles()
        else:
            srcFPN = self.process.srcLayerD[locus][datum][srcComp].FPN
        #The analysis as such is run as GDAL utility

        for dstComp in dstCompL:
            ot =  self.process.dstLayerD[locus][datum][dstComp].comp.celltype
            #The principal steps include 1. changing resolution, 2. Running analysis, 
            #3 resampling back to original resolutuon, and 4 cut out tile
            if self.process.proc.tpiD[dstComp]['resolfac'] == 1:
                scaledSrcFPN = srcFPN

                dsttempFP = os.path.split(self.process.dstLayerD[locus][datum][dstComp].FPN)[0]
                scaledDstFPN = os.path.join(dsttempFP,'scaledTPI.tif')
            else:
                xsize = int(cols/self.process.proc.tpiD[dstComp]['resolfac'])
                ysize = int(lins/self.process.proc.tpiD[dstComp]['resolfac'])
                srctempFP = os.path.split(srcFPN)[0]
                scaledSrcFPN = os.path.join(srctempFP,'scaledDEM.tif')
                dsttempFP = os.path.split(self.process.dstLayerD[locus][datum][dstComp].FPN)[0]
                scaledDstFPN = os.path.join(dsttempFP,'scaledTPI.tif')
                gdaltranslate = GDALstuff(srcFPN, scaledSrcFPN, self.process.params)
                gdaltranslate.TransformOutSize(xsize,ysize,'average',ot)
                
            gdaltxi = GDALstuff(scaledSrcFPN, scaledDstFPN, self.process.params)
            if txi == 'tpi':
                #Run the GDAL TPI utility     
                gdaltxi.TPI()
            elif txi == 'tri':
                #Run the GDAL TPI utility   
                gdaltxi.TRI()
            elif txi == 'rn':
                #Run the GDAL TPI utility   
                gdaltxi.Roughness()
            else:
                exit('unknown terrain index')
            #Copy dst layer and clean up
            if self.process.proc.tpiD[dstComp]['resolfac'] != 1:
                #Resample back to the original size and delete intermediate layers
                xsize = cols
                ysize = lins
                gdaltranslate = GDALstuff(scaledDstFPN, self.process.dstLayerD[locus][datum][dstComp].FPN, self.process.params)
                gdaltranslate.TransformOutSize(xsize,ysize,'nearest',ot)
                os.remove(scaledSrcFPN)
            else:
                gdaltranslate = GDALstuff(scaledDstFPN, self.process.dstLayerD[locus][datum][dstComp].FPN, self.process.params)
                gdaltranslate.TransformOT(ot)
            os.remove(scaledDstFPN)
                
    def _HillShade(self,locus,datum,srcComp,dstCompL):
        if self.process.params.mosaic:
            srcFPN = self._MosaicAdjacentTiles()
        else:
            srcFPN = self.process.srcLayerD[locus][datum][srcComp].FPN
        #The analysis as such is run as GDAL utility
        for dstComp in dstCompL: 
            gdalHillShade = GDALstuff(srcFPN, self.process.dstLayerD[locus][datum][dstComp].FPN, self.process.params)
            gdalHillShade.HillShade()
      
    def _Slope(self,locus,datum,srcComp,dstCompL):
        if self.process.params.mosaic:
            srcFPN = self._MosaicAdjacentTiles()
        else:
            srcFPN = self.process.srcLayerD[locus][datum][srcComp].FPN
        #The analysis as such is run as GDAL utility
        for dstComp in dstCompL: 
            gdalSlope = GDALstuff(srcFPN, self.process.dstLayerD[locus][datum][dstComp].FPN, self.process.params)
            gdalSlope.Slope()
            
    def _Aspect(self,locus,datum,srcComp,dstCompL):
        if self.process.params.mosaic:
            srcFPN = self._MosaicAdjacentTiles()
        else:
            srcFPN = self.process.srcLayerD[locus][datum][srcComp].FPN
        #The analysis as such is run as GDAL utility
        for dstComp in dstCompL: 
            gdalAspect = GDALstuff(srcFPN, self.process.dstLayerD[locus][datum][dstComp].FPN, self.process.params)
            gdalAspect.Aspect()
                   
    def _LandformTPIini(self,locus,datum,srcCompD,dstComp):
        '''
        '''
        standardize = self.process.params.standardize
        expanded = self.process.params.expanded
        tpiTH = self.process.params.tpithreshold
        slopeTH = self.process.params.slopethreshold
        tpisstd = self.process.params.tpisstd
        tpilstd = self.process.params.tpilstd
        for comp in srcCompD:
            if comp == 'tpis':
                TPIS = srcCompD[comp].layer.NPBAND
                tpisnull = srcCompD[comp].layer.cellnull
            elif comp =='tpim':
                TPIL = srcCompD[comp].layer.NPBAND
                tpilnull = srcCompD[comp].layer.cellnull
            elif comp =='slope':
                SLOPE = srcCompD[comp].layer.NPBAND
                slopenull = srcCompD[comp].layer.cellnull
            else:
                exit('unknown input band in landform')
        #Create the dst array
        dstArr = np.empty_like(TPIS, dtype=np.uint8)
        
        if standardize:
            #Standardize the TPI prior to Landform classification
            if tpisstd and tpilstd:
                print ('    Global standardized')

                self._LandformGlobalStandardized(TPIS,TPIL,SLOPE,dstArr,tpiTH,slopeTH,tpisnull,tpilnull,expanded)
            elif tpisstd and not tpilstd:
                print ('    Relative standardized')

                self._LandformRelativeStandardized(TPIS,TPIL,SLOPE,dstArr,tpiTH,slopeTH,tpisnull,tpilnull,expanded)
            else:
                print ('    Local standardized')

                #Calcualte local std for TPI (large and small)
                self._LandformStandardized(TPIS,TPIL,SLOPE,dstArr,tpiTH,slopeTH,tpisnull,tpilnull,expanded)
      
        else:
            print ('    Raw TPI classification')

            #Landform classification against original TPI data
            if expanded:
                self._LandformTPI12(TPIS,TPIL,SLOPE,dstArr,tpiTH,slopeTH)
            else:
                #Landform classification against original TPI data
                self._LandformTPI10(TPIS,TPIL,SLOPE,dstArr,tpiTH,slopeTH)

        dstArr[np.where( TPIS == tpisnull )] = 255
        dstArr[np.where( TPIL == tpilnull )] = 255
        dstArr[np.where( SLOPE == slopenull )] = 255

        #Create the dst layer
        self.process.dstLayerD[locus][datum][dstComp].layer = lambda:None
        #Set the np array as the band
        self.process.dstLayerD[locus][datum][dstComp].layer.NPBAND = dstArr
        
        #copy the geoformat from the src layer
        self.process.dstLayerD[locus][datum][dstComp].CopyGeoformatFromSrcLayer(srcCompD[comp].layer)
        #write the results
        self.process.dstLayerD[locus][datum][dstComp].CreateDSWriteRasterArray()
        
    def _LandformStandardized(self,TPIS,TPIL,SLOPE,dstArr,tpiTH,slopeTH,tpisnull,tpilnull,expanded):
        '''
        '''
        TPISf = TPIS.astype(np.float32)
        TPISf[TPISf==tpisnull] = np.nan
        
        tpiSmean = np.nanmean(TPISf)
        tpiSstd = np.nanstd(TPISf)
        
        print ('    small TPI mean std',tpiSmean, tpiSstd)
        TPISstandard = ((TPISf-tpiSmean)/tpiSstd)*50+0.5
        
        TPILf = TPIL.astype(np.float32)
        TPILf[TPILf==tpilnull] = np.nan
        
        tpiLmean = np.nanmean(TPILf)
        tpiLstd = np.nanstd(TPILf)
        print ('    large TPI mean std',tpiLmean, tpiLstd)

        TPILstandard = ((TPILf-tpiLmean)/tpiLstd)*50+0.5
        if expanded:
            self._LandformTPI12(TPISstandard,TPILstandard,SLOPE,dstArr,tpiTH,slopeTH)
        else:
            #Landform classification against original TPI data
            self._LandformTPI10(TPISstandard,TPILstandard,SLOPE,dstArr,tpiTH,slopeTH)
            
            
    def _LandformGlobalStandardized(self,TPIS,TPIL,SLOPE,dstArr,tpiTH,slopeTH,tpisnull,tpilnull,expanded):
        '''
        '''
        tpisstd = self.process.params.tpisstd
        tpilstd = self.process.params.tpilstd 
        TPISf = TPIS.astype(np.float32)

        TPISstandard = (TPISf/tpisstd)*50+0.5
        
        TPILf = TPIL.astype(np.float32)

        TPILstandard = (TPILf/tpilstd)*50+0.5
        if expanded:
            self._LandformTPI12(TPISstandard,TPILstandard,SLOPE,dstArr,tpiTH,slopeTH)
        else:
            #Landform classification against original TPI data
            self._LandformTPI10(TPISstandard,TPILstandard,SLOPE,dstArr,tpiTH,slopeTH)
            
    def _LandformRelativeStandardized(self,TPIS,TPIL,SLOPE,dstArr,tpiTH,slopeTH,tpisnull,tpilnull,expanded):
        '''
        '''
        tpisstd = self.process.params.tpisstd
        
        TPISf = TPIS.astype(np.float32)
        TPISf[TPISf==tpisnull] = np.nan

        TPISstandard = (TPISf/tpisstd)*50+0.5
        
        TPILf = TPIL.astype(np.float32)
        TPILf[TPILf==tpilnull] = np.nan
        
        
        Sstd = np.nanstd(TPISf)
        
        Lstd = np.nanstd(TPILf)
        relativeSD =  Lstd*tpisstd/Sstd
        
        TPILstandard = (TPILf/relativeSD)*50+0.5
                
        if expanded:
            self._LandformTPI12(TPISstandard,TPILstandard,SLOPE,dstArr,tpiTH,slopeTH)
        else:
            #Landform classification against original TPI data
            self._LandformTPI10(TPISstandard,TPILstandard,SLOPE,dstArr,tpiTH,slopeTH)
            
    def _LandformTPIWeissOriginal(self,TPIS,TPIL,SLOPE,dstArr,tpiTH,slopeTH):
                
        #Plain = 5
        dstArr[np.where( (TPIS > -tpiTH) & (TPIS < tpiTH) & (TPIL > -tpiTH) & (TPIL < tpiTH) & (SLOPE <= slopeTH))] = 50

        #Open slope = 6
        dstArr[np.where( (TPIS > -tpiTH) & (TPIS < tpiTH) & (TPIL > -tpiTH) & (TPIL < tpiTH) & (SLOPE > slopeTH))] = 60
        
        # Mesa or flat ridge = 7
        #upper slope edge = 71
        dstArr[np.where( (TPIS > -tpiTH) & (TPIS < tpiTH) & (TPIL >= tpiTH) & (SLOPE > slopeTH))] = 71
        #mesa or flat ridge top = 72
        dstArr[np.where( (TPIS > -tpiTH) & (TPIS < tpiTH) & (TPIL >= tpiTH) & (SLOPE <= slopeTH))] = 72
        
        #U-shaped valley = 4        
        #valley floor edge = 41
        dstArr[np.where( (TPIS > -tpiTH) & (TPIS < tpiTH) & (TPIL <= -tpiTH) & (SLOPE > slopeTH))] = 41
        #central valley floor  = 42
        dstArr[np.where( (TPIS > -tpiTH) & (TPIS < tpiTH) & (TPIL <= -tpiTH) & (SLOPE <= slopeTH))] = 42
        
        # shallow valley / midslope drainage = 2
        #shallow valley edge = 21
        dstArr[np.where( (TPIS <= -tpiTH) & (TPIL > -tpiTH) & (TPIL < tpiTH) & (SLOPE > slopeTH))] = 21
        #central shallow valley = 21
        dstArr[np.where( (TPIS <= -tpiTH) & (TPIL > -tpiTH) & (TPIL < tpiTH) & (SLOPE <= slopeTH))] = 22
        
        #midslope ridge / hills in valleys = 9
        #midslope ridge / hills in valleys = 91
        dstArr[np.where( (TPIS >= tpiTH) & (TPIL > -tpiTH) & (TPIL < tpiTH) & (SLOPE > slopeTH))] = 91
        #midslope ridge / hills in valleys = 92
        dstArr[np.where( (TPIS >= tpiTH) & (TPIL > -tpiTH) & (TPIL < tpiTH) & (SLOPE <= slopeTH))] = 92
        
        #upland drainage = 3
        #upland drainage slope = 31
        dstArr[np.where( (TPIS <= -tpiTH) & (TPIL >= tpiTH) & (SLOPE > slopeTH))] = 31
        #upland drainage flat = 32
        dstArr[np.where( (TPIS <= -tpiTH) & (TPIL >= tpiTH) & (SLOPE <= slopeTH))] = 32
        
        #canyon / incised stream = 1
        #incised canyon - rapid = 11
        dstArr[np.where( (TPIS <= -tpiTH) & (TPIL <= -tpiTH) & (SLOPE > slopeTH))] = 11
        #incised canyon - calm = 12
        dstArr[np.where( (TPIS <= -tpiTH) & (TPIL <= -tpiTH) & (SLOPE <= slopeTH))] = 12
        
        #peak / mt top = 10
        #peak / mt top = 101
        dstArr[np.where( (TPIS >= tpiTH) & (TPIL >= tpiTH) & (SLOPE > slopeTH))] = 101
        #peak / mt top = 102
        dstArr[np.where( (TPIS >= tpiTH) & (TPIL >= tpiTH) & (SLOPE <= slopeTH))] = 102
        
        #local ridge = 8
        #local ridge = 81
        dstArr[np.where( (TPIS >= tpiTH) & (TPIL < -tpiTH) & (SLOPE > slopeTH))] = 81
        #local ridge = 82
        dstArr[np.where( (TPIS >= tpiTH) & (TPIL < -tpiTH) & (SLOPE <= slopeTH))] = 82  
        
    def _LandformTPI10(self,TPIS,TPIL,SLOPE,dstArr,tpiTH,slopeTH):
        '''
        '''
        #Plain 
        dstArr[np.where( (TPIS > -tpiTH) & (TPIS < tpiTH) & (TPIL > -tpiTH) & (TPIL < tpiTH) & (SLOPE <= slopeTH) )] = 10
        #Open slope
        dstArr[np.where( (TPIS > -tpiTH) & (TPIS < tpiTH) & (TPIL > -tpiTH) & (TPIL < tpiTH) & (SLOPE > slopeTH) )] = 11
        # Mesa or flat ridge
        dstArr[np.where( (TPIS > -tpiTH) & (TPIS < tpiTH) & (TPIL >= tpiTH) )] = 12
        #U-shaped valley       
        dstArr[np.where( (TPIS > -tpiTH) & (TPIS < tpiTH) & (TPIL <= -tpiTH) )] = 9
        # shallow valley / midslope drainage
        dstArr[np.where( (TPIS <= -tpiTH) & (TPIL > -tpiTH) & (TPIL < tpiTH) )] = 6
        #midslope ridge / hills in valleys
        dstArr[np.where( (TPIS >= tpiTH) & (TPIL > -tpiTH) & (TPIL < tpiTH) )] = 14
        #upland drainage
        dstArr[np.where( (TPIS <= -tpiTH) & (TPIL >= tpiTH) )] = 8
        #canyon / incised stream
        dstArr[np.where( (TPIS <= -tpiTH) & (TPIL <= -tpiTH) )] = 5
        #peak / mt top
        dstArr[np.where( (TPIS >= tpiTH) & (TPIL >= tpiTH) )] = 16
        #local ridge
        dstArr[np.where( (TPIS >= tpiTH) & (TPIL <= -tpiTH) )] = 13
        
    def _LandformTPI12(self, TPIS, TPIL, SLOPE, dstArr,tpiTH,slopeTH):
        '''
        '''
        #Incised stream /canyon
        dstArr[np.where( (TPIS <= -tpiTH) & (TPIL <= -tpiTH) )] = 5
        #channel / shallow valley
        dstArr[np.where( (TPIS <= -tpiTH) & (TPIL > -tpiTH) & (TPIL < tpiTH) & (SLOPE < slopeTH))] = 6
        #Midslope drainage
        dstArr[np.where( (TPIS <= -tpiTH) & (TPIL > -tpiTH) & (TPIL < tpiTH) & (SLOPE >= slopeTH))] = 7
        #Upland drainage, Headwater 
        dstArr[np.where( (TPIS <= -tpiTH) & (TPIL >= tpiTH))] = 8
        #U-shaped valley 
        dstArr[np.where( (TPIS > -tpiTH) & (TPIS < tpiTH) & (TPIL <= -tpiTH) )] = 9
        #Plain 
        dstArr[np.where( (TPIS > -tpiTH) & (TPIS < tpiTH) & (TPIL > -tpiTH) & (TPIL < tpiTH) & (SLOPE <= slopeTH))] = 10
        #Open slope = 6, changed to 7
        dstArr[np.where( (TPIS > -tpiTH) & (TPIS < tpiTH) & (TPIL > -tpiTH) & (TPIL < tpiTH) & (SLOPE > slopeTH))] = 11
        #mesa or flat ridge top = 8 (should include peat dome edges)
        dstArr[np.where( (TPIS > -tpiTH) & (TPIS < tpiTH) & (TPIL >= tpiTH))] = 12
        #local ridge
        dstArr[np.where( (TPIS >= tpiTH) & (TPIL <= -tpiTH) )] = 13
        #hills in valleys 
        dstArr[np.where( (TPIS >= tpiTH) & (TPIL > -tpiTH) & (TPIL < tpiTH) & (SLOPE < slopeTH))] = 14
        #midslope ridge = 92, new = 11
        dstArr[np.where( (TPIS >= tpiTH) & (TPIL > -tpiTH) & (TPIL < tpiTH) & (SLOPE <= slopeTH))] = 15
        #Mountain top
        dstArr[np.where( (TPIS >= tpiTH) & (TPIL >= tpiTH) )] = 16
        
        