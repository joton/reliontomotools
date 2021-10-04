#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import re
from glob import glob
import numpy as np
import pandas as pd
from docopt import docopt
from scipy.interpolate import interpn
import mrcfile
from tqdm import tqdm
from ._version import __version__
from .utils import rotX, rotY, rotZ, getShiftMatrix
from .fileIO import readStarFile, writeStarFile, WarpXMLHandler

__all__ = ['WarpTomo2Relion', 'warpTomo2RelionProgram']
# % Warp To Relion Converter


class WarpTomo2Relion():

    def __init__(self, origStackFn, xmlFname, thickness, tomoName=None,
                 angpix=-1,
                 flipZ=False, flipYZ=False, flipAngles=False, hand=-1,
                 aliDims=None, offset3D=[0, 0, 0]):

        self.warp = WarpXMLHandler(xmlFname)

        self.origStackFn = origStackFn

        self._tomoName = tomoName

        self.angpix = angpix if angpix > 0 else self.warp.pixelSize

        self.flipAngles = np.logical_xor(self.warp.AreAnglesInverted,
                                         flipAngles)

        # properties flipz and flipYZ are regarding to IMOD, which are True by
        # default. Input parameters are regarding to Warp
        self.flipZ = not flipZ
        self.flipYZ = not flipYZ
        self.hand = hand

        mrc = mrcfile.open(origStackFn)
        self.w_ts = mrc.header['nx']
        self.h_ts = mrc.header['ny']
        self.fc_ts = mrc.header['nz']
        self.origCenter = np.array([self.w_ts - 1, self.h_ts - 1])/2.

        self.aliDims = aliDims

        self.thickness = thickness
        self.offset3D = np.array(offset3D)

        self.processRelionVolToProjTransforms()

        self.processWarpAnglesCorrectionTransforms()

    def getImodVolToProjTransforms(self):

        angles = self.warp.AxisAngle
        # xyOffsets in Warp are already transformed
        xOff = self.warp.AxisOffsetX/self.angpix
        yOff = self.warp.AxisOffsetY/self.angpix

        tiltAngles = self.warp.Angles
        if self.flipAngles:
            tiltAngles *= -1

        fc = len(angles)

        if fc != self.fc_ts:
            raise Exception(f'getImodVolToProjTransforms: Tilt series size in '
                            'Warp {fc} and stack size {self.fc_ts} mismatch.')

        origOffProjAdd = getShiftMatrix([*self.origCenter, 0])

        self.flipXYDim = False

        if self.aliDims is not None:
            self.w_ali = self.aliDims[0]
            self.h_ali = self.aliDims[1]
        else:
            if abs(angles[0]) % 90 > 45:
                self.flipXYDim = True
                self.w_ali = self.h_ts
                self.h_ali = self.w_ts
            else:
                self.w_ali = self.w_ts
                self.h_ali = self.h_ts

        self.aliCenter = np.array([self.w_ali - 1, self.h_ali - 1])/2.

        origOffVolSub = np.asmatrix(np.identity(4))
        origOffVolSub[[0, 2], 3] = -self.aliCenter
        origOffVolSub[1, 3] = -(self.thickness)/2.0

        rot90X = rotX(-90)

        xforms = list()
        for k in range(fc):
            Aaxis = rotZ(angles[k])
            # Ai = np.linalg.inv(A)
            S = getShiftMatrix((xOff[k], yOff[k], 0))

            Atilt = rotY(tiltAngles[k])

            xforms.append(origOffProjAdd@S@Aaxis@Atilt@rot90X@origOffVolSub)

        return xforms

    def processRelionVolToProjTransforms(self):

        toImodOrigin3D = getShiftMatrix((-1, 0, -1))

        imodTransforms = self.getImodVolToProjTransforms()

        if (self.flipYZ):
            if (self.flipZ):
                AflipYZ = np.matrix([[1, 0, 0, 0],
                                    [0, 0, -1, self.thickness-1],
                                    [0, 1, 0, 0],
                                    [0, 0, 0, 1]])
            else:
                AflipYZ = np.matrix([[1, 0, 0, 0],
                                    [0, 0, 1, 0],
                                    [0, 1, 0, 0],
                                    [0, 0, 0, 1]])

            self.h_out = self.h_ali
            self.d_out = self.thickness
        else:
            AflipYZ = np.asmatrix(np.identity(4))
            self.h_out = self.thickness
            self.d_out = self.h_ali

        self.w_out = self.w_ali

        Aoffset3D = getShiftMatrix(self.offset3D)

        Apre = toImodOrigin3D@AflipYZ@Aoffset3D

        relionTransforms = list()

        for k in range(self.fc_ts):
            relionTransforms.append(imodTransforms[k]@Apre)

        self.relionTransforms = relionTransforms

    def processWarpAnglesCorrectionTransforms(self):

        if 'GridAngleX' not in self.warp.params:
            self.warpAngleTransforms = np.identity(4)
            return

        # Warp angle reference system is inverted
        angleX = -np.squeeze(self.warp.GridAngleX)
        angleY = -np.squeeze(self.warp.GridAngleY)
        angleZ = -np.squeeze(self.warp.GridAngleZ)

        angleXndims = angleX.ndim
        # If only one value, we use it for all tilt angles
        if angleXndims == 0:
            angleX = angleX * np.ones(self.fc_ts)
            angleY = angleY * np.ones(self.fc_ts)
            angleZ = angleZ * np.ones(self.fc_ts)
        elif angleXndims > 1:
            raise Exception('Local Warp angles correction is not currently '
                            'supported.')

        origOffProjAdd = getShiftMatrix([*self.origCenter, 0])
        origOffProjSub = getShiftMatrix([*-self.origCenter, 0])

        warpAngleTransforms = list()

        for k in range(self.fc_ts):

            Acor = origOffProjAdd@\
                rotZ(angleZ[k])@rotY(angleY[k])@rotX(angleX[k])@origOffProjSub

            warpAngleTransforms.append(Acor)

        self.warpAngleTransforms = warpAngleTransforms

    def getTiltsCTF(self):

        defocus = np.squeeze(self.warp.GridCTF*1.e4)
        defDeltah = np.squeeze(self.warp.GridCTFDefocusDelta/2*1.e4)
        defAngle = np.squeeze(self.warp.GridCTFDefocusAngle)
        dose = np.squeeze(self.warp.Dose)

        if len(defocus.shape) > 1:
            raise Exception('Local defocus estimation is not currently '
                            'supported.')

        defMax = defocus + defDeltah
        defMin = defocus - defDeltah

        return np.column_stack((defMax, defMin, defAngle, dose))

    def getInterpolated(self, param, gridCoords):

        if param not in self.warp.params:
            raise Exception(f'{param} not present in warp data.')

        gridData = self.warp[param]
        # If only one single value instead of a grid, return it
        sqData = np.squeeze(gridData)
        if sqData.ndim == 0:
            return sqData

        iSize = gridData.shape
        XXref = list()
        for size in iSize:
            XXref.append(np.linspace(0, 1, size))

        return interpn(XXref, gridData, gridCoords)

    def getWarpLocalCorrection(self, dataPart, applyGlobWarp=True,
                               applyLocalWarpVol=True, applyLocalWarpImg=True):

        nPart = len(dataPart)

        coordDiff = np.zeros((nPart, self.fc_ts, 3))

        coords = dataPart[['_rlnCoordinateX',
                           '_rlnCoordinateY',
                           '_rlnCoordinateZ']].values.astype(np.float)

        coordsTmp = coords

        if applyLocalWarpVol:
            dose = self.warp.Dose
            doseMax = dose.max()
            doseMin = dose.min()
            doseStep = 1./(doseMax - doseMin)

            gridCoords4Norm = np.empty((nPart, 4))
            gridCoords4Norm[:, 0] = coords[:, 0]/self.w_out
            gridCoords4Norm[:, 1] = coords[:, 1]/self.h_out
            gridCoords4Norm[:, 2] = coords[:, 2]/self.d_out

        if applyLocalWarpImg:
            gridStep = 1./(self.fc_ts - 1)

        for kt in range(self.fc_ts):

            if applyGlobWarp:
                rotProj2Vol = (self.warpAngleTransforms[kt]@
                               self.relionTransforms[kt]).I
            else:
                rotProj2Vol = self.relionTransforms[kt].I

            rotVol2Proj = self.relionTransforms[kt]

            if applyLocalWarpVol:
                gridCoords4Norm[:, 3] = (dose[kt] - doseMin)*doseStep

                gridVolWarp = np.empty_like(coords)
                gridVolWarp[:, 0] = self.getInterpolated('GridVolumeWarpX',
                                                         gridCoords4Norm)
                gridVolWarp[:, 1] = self.getInterpolated('GridVolumeWarpY',
                                                         gridCoords4Norm)
                gridVolWarp[:, 2] = self.getInterpolated('GridVolumeWarpZ',
                                                         gridCoords4Norm)
                coordsTmp = coords + gridVolWarp/self.angpix

            for kp in range(nPart):
                coordVolOrig = coords[kp, :]
                coordVol = coordsTmp[kp, :]

                coordProj = rotVol2Proj@getShiftMatrix(coordVol)

                if applyLocalWarpImg:
                    gridCoordProjNorm = np.array([
                        coordProj[0, 3]/self.w_ts,
                        coordProj[1, 3]/self.h_ts,
                        kt * gridStep])

                    # if np.any(gridCoordProjNorm > 1):
                    #     print(f'kt = {kt} - kp = {kp}')
                    #     print(f'coordVolOrig={coordVolOrig}')
                    #     print(f'coordVol={coordVol}')
                    #     print(f'coordProj={coordProj}')

                    gridMovWarpX = self.getInterpolated('GridMovementX',
                                                        gridCoordProjNorm)
                    gridMovWarpY = self.getInterpolated('GridMovementY',
                                                        gridCoordProjNorm)
                    coordProj[0, 3] -= gridMovWarpX/self.angpix
                    coordProj[1, 3] -= gridMovWarpY/self.angpix

                newCoordVol = (rotProj2Vol@coordProj)[0:3, 3].T

                coordDiff[kp, kt, :] = newCoordVol - coordVolOrig

        coordDiff *= self.angpix

        return coordDiff

    def getRelionTomoTables(self, applyGlobWarp=True):

        tomoName = self._tomoName if self._tomoName is not None else \
            self.warp.tomoName

        labelsGlobal = ['_rlnTomoName',
                        '_rlnTomoTiltSeriesName',
                        '_rlnTomoFrameCount',
                        '_rlnTomoSizeX',
                        '_rlnTomoSizeY',
                        '_rlnTomoSizeZ',
                        '_rlnTomoHand',
                        '_rlnOpticsGroupName',
                        '_rlnTomoTiltSeriesPixelSize',
                        '_rlnVoltage',
                        '_rlnSphericalAberration',
                        '_rlnAmplitudeContrast',
                        '_rlnTomoImportFractionalDose']

        dataGlobal = pd.DataFrame(index=range(1), columns=labelsGlobal)

        dataGlobal['_rlnTomoName'] = tomoName
        dataGlobal['_rlnTomoTiltSeriesName'] = self.origStackFn
        dataGlobal['_rlnTomoFrameCount'] = str(self.fc_ts)
        dataGlobal['_rlnTomoSizeX'] = str(self.w_out)
        dataGlobal['_rlnTomoSizeY'] = str(self.h_out)
        dataGlobal['_rlnTomoSizeZ'] = str(self.d_out)
        dataGlobal['_rlnTomoHand'] = str(self.hand)
        dataGlobal['_rlnOpticsGroupName'] = 'opticsGroup1'
        dataGlobal['_rlnTomoTiltSeriesPixelSize'] = '{:.6g}'.format(self.angpix)
        dataGlobal['_rlnVoltage'] = '{:.6g}'.format(self.warp.Voltage)
        dataGlobal['_rlnSphericalAberration'] = '{:.6g}'.format(self.warp.Cs)
        dataGlobal['_rlnAmplitudeContrast'] = '{:.6g}'.format(self.warp.Amplitude)
        dataGlobal['_rlnTomoImportFractionalDose'] = '{:.6g}'.\
            format(self.warp.fracDose)

        labelsTomo = ['_rlnTomoProjX',
                      '_rlnTomoProjY',
                      '_rlnTomoProjZ',
                      '_rlnTomoProjW',
                      '_rlnDefocusU',
                      '_rlnDefocusV',
                      '_rlnDefocusAngle',
                      '_rlnCtfScalefactor',
                      '_rlnMicrographPreExposure']

        dataTomo = pd.DataFrame(index=range(self.fc_ts), columns=labelsTomo)

        ctfData = self.getTiltsCTF()

        for k in range(self.fc_ts):
            dataTilt = dataTomo.iloc[k]

            if applyGlobWarp:
                rotMat = self.warpAngleTransforms[k]@self.relionTransforms[k]
            else:
                rotMat = self.relionTransforms[k]

            dataTilt['_rlnTomoProjX'] = '[{:.13g},{:.13g},{:.13g},{:.13g}]'.\
                format(*rotMat.A[0, :])
            dataTilt['_rlnTomoProjY'] = '[{:.13g},{:.13g},{:.13g},{:.13g}]'.\
                format(*rotMat.A[1, :])
            dataTilt['_rlnTomoProjZ'] = '[{:.13g},{:.13g},{:.13g},{:.13g}]'.\
                format(*rotMat.A[2, :])
            dataTilt['_rlnTomoProjW'] = '[{:.13g},{:.13g},{:.13g},{:.13g}]'.\
                format(*rotMat.A[3, :])
            dataTilt['_rlnDefocusU'] = '{:.6g}'.format(ctfData[k, 0])
            dataTilt['_rlnDefocusV'] = '{:.6g}'.format(ctfData[k, 1])
            dataTilt['_rlnDefocusAngle'] = '{:.6g}'.format(ctfData[k, 2])
            dataTilt['_rlnCtfScalefactor'] = '1.0'
            dataTilt['_rlnMicrographPreExposure'] = '{:.6g}'.format(ctfData[k, 3])

        tomoTables = dict()
        tomoTables['global'] = dataGlobal
        tomoTables[tomoName] = dataTomo

        return tomoTables

    def getRelionMotionTable(self, dataPart, applyGlobWarp=True,
                             applyLocalWarpVol=True, applyLocalWarpImg=True):

        if '_rlnTomoParticleName' not in dataPart.columns:
            raise Exception('ERROR: getRelionMotionTable requires '
                            '_rlnTomoParticleName label in particles file.')

        labelsMotion = ['_rlnOriginXAngst',
                        '_rlnOriginYAngst',
                        '_rlnOriginZAngst']

        motionTable = dict()
        nPart = len(dataPart)

        motion = self.getWarpLocalCorrection(dataPart, applyGlobWarp,
                                             applyLocalWarpVol,
                                             applyLocalWarpImg).astype(str)

        for k in range(nPart):
            partName = dataPart.loc[k, '_rlnTomoParticleName']

            motionTable[partName] = pd.DataFrame(motion[k],
                                                 columns=labelsMotion)

        return motionTable

    def writeTomogramStarFile(self, tomoOutFname, particlesFn=None,
                              motionOutFn=None, applyGlobWarp=True,
                              applyLocalWarpVol=True, applyLocalWarpImg=True):

        if particlesFn is not None:
            dataPart = readStarFile(particlesFn, 'particles')
            dataPart = dataPart.loc[dataPart._rlnTomoName == self._tomoName]
            dataPart = dataPart.reset_index(drop=True)
            _applyGlobWarp = applyGlobWarp
        else:
            dataPart = None
            _applyGlobWarp = False

        tomoTables = self.getRelionTomoTables(_applyGlobWarp)

        writeStarFile(tomoOutFname, tomoTables)

        if dataPart is not None:
            motion = self.getRelionMotionTable(dataPart, _applyGlobWarp,
                                               applyLocalWarpVol,
                                               applyLocalWarpImg)

            nPart = len(dataPart)
            motion = dict(general=pd.DataFrame(str(nPart),
                                               index=['_rlnParticleNumber'],
                                               columns=['value']),
                          **motion)

            writeStarFile(motionOutFn, motion)


def warpTomo2RelionProgram(args=None):

    doc = """warptomo2relion: converts warp tomo metadata to Relion.

    Usage:
      warptomo2relion -i <xml_template> -s <ts_template> -d <thickness> -o <outDir> [options]

    Arguments:
         -i <xml_template>  Single Warp XML filename or file list/template with
                            wild cards within quotes.
          -s <ts_template>  Single tilt series MRC filename or file list/template
                            with wild cards within quotes.
    -d --depth <thickness>  Tomogram thickness used for all tomograms, in pixels
                            at bin1.
               -o <outDir>  Directory where output files are created.
    Options:
           --tn <tomoName>  Tomo name template to match XML and tilt series files.
                            It is used as tomoName label. [Default: TS_01]
         -p <particles_fn>  Relion particle set to convert Warp/M local
                            improvements. This option creates a trajectory file.
      --ignore_global_warp  Do not apply global warp angles correction.
         --ignore_vol_warp  Do not apply local volume warp offsets correction.
         --ignore_img_warp  Do not apply local image warp offsets correction.
                   --flipZ  Change the sign of Z coordinate.
                 --flipAng  Change the sign of tilt angles.
                --hand <h>  Handedness of the tilt geometry. [Default: -1]

                 -h --help  Show this screen.
              -v --version  Show version.

    """

    arguments = docopt(doc, version=__version__)
    # print(arguments)

    xmlTmpl = arguments['-i']
    tsTmpl = arguments['-s']
    outRoot = arguments['-o']
    thickness = int(arguments['--depth'])
    tomoName = arguments['--tn']
    particlesFn = arguments['-p']
    doTraject = particlesFn is not None

    flipZ = arguments['--flipZ']
    flipAngles = arguments['--flipAng']
    hand = arguments['--hand']

    os.makedirs(outRoot, exist_ok=True)
    tomoOutFname = os.path.join(outRoot, 'tomograms.star')
    opSetOutFname = os.path.join(outRoot, 'optimisation_set.star')

    if doTraject:
        motionOutFn = os.path.join(outRoot, 'motion.star')
        applyGlobalWarp = not arguments['--ignore_global_warp']
        applyLocalWarpVol = not arguments['--ignore_vol_warp']
        applyLocalWarpImg = not arguments['--ignore_img_warp']

        if os.path.exists(motionOutFn):
            os.remove(motionOutFn)
    else:
        motionOutFn = None
        applyGlobalWarp = False
        applyLocalWarpVol = False
        applyLocalWarpImg = False

    xmlList = glob(xmlTmpl)
    tsList = glob(tsTmpl)
    nTomos = len(xmlList)

    # Get Tomo numbers and labels from xml/ts filenames
    p = 0
    for c in tomoName:
        if c.isdigit():
            break
        p += 1
    tnTmpl = tomoName[0:p]

    xmlNums = np.array([int(re.findall(f'{tnTmpl}\d+', name)[0][p:])
                        for name in xmlList])
    xmlLabels = np.array([re.findall(f'{tnTmpl}\d+', name)[0]
                         for name in xmlList])
    tsLabels = np.array([re.findall(f'{tnTmpl}\d+', name)[0]
                         for name in tsList])

    if nTomos == 1:

        tomoLabel = str(xmlLabels[0])
        warpTomo = WarpTomo2Relion(tsList[0], xmlList[0], thickness,
                                   tomoName=tomoLabel, flipZ=flipZ,
                                   flipAngles=flipAngles, hand=hand)

        warpTomo.writeTomogramStarFile(tomoOutFname, particlesFn,
                                       motionOutFn, applyGlobalWarp,
                                       applyLocalWarpVol,
                                       applyLocalWarpImg)
    else:

        xmlOrder = np.argsort(xmlNums)

        if doTraject:
            dataPart = readStarFile(particlesFn, 'particles')
            motion = dict()
        else:
            dataPartTomo = None

        tomoStarTables = dict()
        tomoStarTables['global'] = None

        globalList = list()

        pbar = tqdm(total=nTomos, unit="Tomos", ncols=80)

        for kx in xmlOrder:

            tomoLabel = str(xmlLabels[kx])
            pbar.write(f'Tomogram {tomoLabel}')

            kt = np.nonzero(tomoLabel == tsLabels)[0][0]

            warpTomo = WarpTomo2Relion(tsList[kt], xmlList[kx], thickness,
                                       tomoName=tomoLabel, flipZ=flipZ,
                                       flipAngles=flipAngles, hand=hand)

            tomoData = warpTomo.getRelionTomoTables(applyGlobalWarp)

            globalList.append(tomoData['global'])
            tomoStarTables[tomoLabel] = tomoData[tomoLabel]

            if doTraject:
                dataPartTomo = dataPart.loc[dataPart._rlnTomoName == tomoLabel]
                dataPartTomo = dataPartTomo.reset_index(drop=True)

                motionData = warpTomo.getRelionMotionTable(dataPartTomo,
                                                           applyGlobalWarp,
                                                           applyLocalWarpVol,
                                                           applyLocalWarpImg)

                motion = dict(motion, **motionData)

            del warpTomo
            pbar.update()

        tomoStarTables['global'] = pd.concat(globalList)
        writeStarFile(tomoOutFname, tomoStarTables)

        if doTraject:
            nMotionPart = len(motion)
            motion = dict(general=pd.DataFrame(str(nMotionPart),
                                               index=['_rlnParticleNumber'],
                                               columns=['value']),
                          **motion)
            writeStarFile(motionOutFn, motion)

    opSetData = pd.DataFrame(tomoOutFname, index=['_rlnTomoTomogramsFile'],
                             columns=['value'])
    if doTraject:
        opSetData = opSetData.append(pd.DataFrame([particlesFn, motionOutFn],
                                      index=['_rlnTomoParticlesFile',
                                             '_rlnTomoTrajectoriesFile'],
                                      columns=['value']))

    writeStarFile(opSetOutFname, opSetData)
