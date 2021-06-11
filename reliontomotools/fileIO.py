#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import pandas as pd
import numpy as np
import io
import xmltodict
from glob import glob

from ._version import __version__


__all__ = ['readStarFile', 'writeStarFile', 'WarpXMLHandler']


def getTable(lines):

    nLines = len(lines)
    lId = 0
    labels = list()

    if 'loop_' in lines[lId]:
        lId += 1
        while lines[lId][0] == '_':
            labels.append(lines[lId].split()[0])
            lId += 1

        datalist = [lines[k].split() for k in range(lId, nLines)]
        data = pd.DataFrame(datalist, columns=labels)
    else:
        values = list()
        while lId < nLines and lines[lId][0] == '_':
            label, value = lines[lId].split()
            labels.append(label)
            values.append(value)
            lId += 1
        data = pd.DataFrame(1, columns=['value'], index=labels)
        data['value'] = values

    return data


def readStarFile(fName, mytable=''):

    with open(fName, 'r') as file:
        lines = file.read().splitlines()

    # Remove empty lines
    lines = [line for line in lines if line]
    # Remove single space lines
    lines = [line for line in lines if line not in ' ']
    # Remove comment lines
    lines = [line for line in lines if line[0] not in '#']
    nLines = len(lines)

    # Get tables
    tables = list()
    tLines = list()  # lines where tables start
    for k, line in enumerate(lines):
        if 'data_' in line[:10]:
            tables.append(line[5:])
            tLines.append(k)

    nTables = len(tables)
    tRanges = [None]*nTables

    for k in range(nTables-1):
        tRanges[k] = (tLines[k] + 1, tLines[k+1])
    tRanges[-1] = (tLines[-1] + 1, nLines)

    data = dict()

    for k, table in enumerate(tables):
        if np.diff(tRanges[k]) > 1:
            data[table] = getTable(lines[slice(*tRanges[k])])
        else:
            data[table] = ''

    if mytable:
        return data[mytable]
    else:
        return data


def writeStarFile(fName, data, tableName=''):

    if isinstance(data, dict):
        mydata = data
    elif isinstance(data, pd.core.frame.DataFrame):
        mydata = {}
        mydata[tableName] = data

    file = open(fName, 'w')

    for tName, table in mydata.items():

        file.write(f'\n\n# reliontomotools {__version__}\n\n')
        file.write(f'data_{tName}\n')

        ncol = len(table.columns)
        if ncol == 1 and not isinstance(table.index[0], int):

            file.write('\n')

            for label, value in table.itertuples():
                file.write(f'{label}\t{value}\n'.expandtabs(60))

            file.write('\n')

        else:  # regular table

            labels = table.columns

            file.write('\nloop_\n')

            for k, label in enumerate(labels):
                file.write(f'{label} #{k+1}\n')

            for k, line in table.iterrows():
                linesep = "\t".join(line)
                file.write(f'{linesep}\n')

    file.close()


class WarpXMLHandler():

    def __init__(self, xmlFname):

        self._xmlFname = xmlFname
        file = io.open(xmlFname, 'r', encoding='utf-16-le')
        self.data = xmltodict.parse(file.read())['TiltSeries']

    @property
    def params(self):
        return list(self.data)

    def getSimpleAttrib(self, param):

        if param in ['@AreAnglesInverted', '@UnselectFilter']:
            return self.getBoolAttrib(param)
        elif param in ['@PlaneNormal']:
            return self.getParamArray(param, ',')
        else:
            return float(self.data[param])

    def getBoolAttrib(self, param):

        return self.data[param] == 'True'

    def getParamArray(self, param, sep='\n'):

        items = self.data[param].split(sep)

        if items[0] in ['True', 'False']:
            return np.array(items) == 'True'
        elif items[0][0] in [str(x) for x in set('-0123456789.')]:
            return np.array(items, dtype=np.float)
        else:
            return items

    def getParam(self, param):

        if param in self.data:
            if param[0] == '@':
                return self.getSimpleAttrib(param)
            else:
                return self.getParamArray(param)
        else:
            if param[0] != '@':
                return self.getParam('@'+param)
            else:
                raise Exception(f'Parameter {param} not present in '
                                '{self._xmlFname}.')

    def getParamGrid(self, param):

        gridParam = self.data[param]
        w = int(gridParam['@Width'])
        h = int(gridParam['@Height'])
        d = int(gridParam['@Depth'])

        if '@Duration' in gridParam:
            t = int(gridParam['@Duration'])
            dims = 4
            shape = (w, h, d, t)
        else:
            dims = 3
            shape = (w, h, d)

        if np.prod(shape) == 1:
            gridParamList = [gridParam['Node']]
        else:
            gridParamList = gridParam['Node']

        gridData = np.empty(shape)

        for item in gridParamList:
            x = int(item['@X'])
            y = int(item['@Y'])
            z = int(item['@Z'])
            v = float(item['@Value'])
            pos = (x, y, z)

            if dims > 3:
                tw = int(item['@W'])
                pos = pos + (tw,)

            gridData[pos] = v

        return gridData

    def getParamCTF(self, param):

        paramsCTF = self.data['CTF']['Param']

        for item in paramsCTF:
            if item['@Name'] == param:
                return float(item['@Value'])

        raise Exception(f'CTF param {param} not present in '
                        'file {self._xmlFname}')

    @property
    def pixelSize(self):
        return self.getParamCTF('PixelSize')

    @property
    def Amplitude(self):
        return self.getParamCTF('Amplitude')

    @property
    def Cs(self):
        return self.getParamCTF('Cs')

    @property
    def Voltage(self):
        return self.getParamCTF('Voltage')

    @property
    def fracDose(self):

        fracDoseV = np.diff(np.sort(self.Dose))
        p = np.flatnonzero(np.isclose(np.diff(fracDoseV), 0))[0]

        return fracDoseV[p]

    @property
    def tomoName(self):
        tiltFname = self.MoviePath[0]
        p = tiltFname.rfind('_')

        return tiltFname[:p]

    def __getattr__(self, name):

        return self.__getitem__(name)

    def __getitem__(self, name):

        # try:
        if 'Grid' in name:
            return self.getParamGrid(name)
        else:
            return self.getParam(name)
        # except Exception as ins:
            # print(ins)


def cleanDir(path):

    files = glob(os.path.join(path, '*'))

    for f in files:
        try:
            os.remove(f)
        except OSError as e:
            print("Error: %s : %s" % (f, e.strerror))

