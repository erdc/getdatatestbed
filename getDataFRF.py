# -*- coding: utf-8 -*-
"""This is a class definition designed to get data from the FRF thredds server.

@author: Spicer Bak, PhD
@contact: spicer.bak@usace.army.mil
@organization: USACE CHL FRF
"""
import collections
import datetime as DT
import os
import pickle as pickle
import time
import warnings
from posixpath import join as urljoin
import socket
import netCDF4 as nc
import numpy as np
import pandas as pd

from testbedutils import geoprocess as gp, sblib as sb

def gettime(allEpoch, epochStart, epochEnd):
    """This function opens the netcdf file, and retrieves time.
    
    It pulls the dates of interest from the THREDDS (data loc) server based on d1,d2, and data location it returns
     the indices in the NCML file of the dates d1>=time>d2

    It was modified to check if there are duplicate times, and only produces indices with unique
    times


    Args:
        allEpoch (list, float): a list of floats that has epoch times in it
        epochStart (float): start time in epoch
        epochEnd (float): end time in epoch

    Returns:
        index  of dates between

    """
    try:
        mask = (allEpoch >= epochStart) & (allEpoch < epochEnd)
        idx = np.argwhere(mask).squeeze()
        if np.size(idx) == 0:
            idx = None
    except TypeError:  # when None's are handed for allEpoch
        idx = None
    finally:
        return idx

def getnc(dataLoc, callingClass, dtRound=60, **kwargs):
    """Function grabs the netCDF file interested.
    
    Responsible for drilling down to specific monthly file if applicable to speed things up.

    Args:
        dataLoc (str):
        THREDDS (str): a key associated with the server location
        callingClass (str): which class calls this
        dtRound(int): rounding the times returned from the server (Default=60 (s))

    Keyword Args:
        start: if given, will parse out to monthly netCDF file (if query is in same month)
        end: if given, will parse out to monthly netCDF file (if query is in same month)

    Returns:
        object:

    TODO: could use thredds crawler to more efficiently pick files to pull from.  This would save query time
    """
    # toggle my data location
    start = kwargs.get('start', None)
    end = kwargs.get('end', None)
    FRFdataloc = u'http://134.164.129.55/thredds/dodsC/'
    chlDataLoc = u'https://chldata.erdc.dren.mil/thredds/dodsC/'
    # a list of data sets (just the ncml) that shouldn't drill down to monthly file
    doNotDrillList = ['survey']

    # chose which server to select based on IP
    ipAddress = socket.gethostbyname(socket.gethostname())
    if ipAddress.startswith('134.164.129'):  # FRF subdomain
        THREDDSloc = FRFdataloc
        pName = u'FRF'
    else:
        THREDDSloc = chlDataLoc
        pName = u'frf'
        
    if callingClass == 'getDataTestBed':  # overwrite pName if calling for model data
        pName = u'cmtb'
    
    # now set URL for netCDF file call,
    if start is None and end is None:
        ncfileURL = urljoin(THREDDSloc, pName, dataLoc)
    elif isinstance(start, float) and isinstance(end, float):  # then we assume epoch
        raise NotImplementedError(
            'check conversion for floats (epoch time), currently needs to be datetime object')
        # ncfileURL = urljoin(THREDDSloc, pName, monthlyPath)

    elif isinstance(start, DT.datetime) and isinstance(end, DT.datetime) \
        and (start.year == end.year and start.month == end.month) \
        and ~np.in1d(doNotDrillList, dataLoc.split('/')).any():
        # this section dives to the specific month's datafile if it's within the same month
        dataLocSplit = os.path.split(dataLoc)
        fileparts = dataLocSplit[0].split('/')
        if fileparts[0] == 'oceanography':
            field = 'ocean'
        else:
            field = fileparts[0]
        try:  # this will work for get Obs
            fname = "{}-{}_{}_{}_{}{:02d}.nc".format(pName.upper(), field, fileparts[1],
                                                     fileparts[2], start.year,
                                                     start.month)
        except IndexError:  # works for getDataTestBed class
            fname = u"{}-{}_{}_{}{:02d}.nc".format(pName.upper(), field, fileparts[1], start.year,
                                                   start.month)
        
        ncfileURL = urljoin(THREDDSloc, pName, dataLocSplit[0], str(start.year), fname)
    else:  # function couldn't be more efficient, default to old way
        ncfileURL = urljoin(THREDDSloc, pName, dataLoc)
    
    # ___________________ go now to open file   ___________________________________________
    finished, n, maxTries = False, 0, 3  # initializing variables to iterate over
    ncFile, allEpoch = None, None  # will return None's when URL doesn't exist
    while not finished and n < maxTries:
        try:
            ncFile = nc.Dataset(ncfileURL)  # get the netCDF file
            allEpoch = ncFile['time'][:]
            finished = True
        except IOError:
            print('Error reading {}, trying again {}/{}'.format(ncfileURL, n + 1, maxTries))
            time.sleep(5)  # time in seconds to wait
            n += 1  # iteration number
    
    return ncFile, sb.baseRound(allEpoch, base=dtRound) # round to nearest minute

def removeDuplicatesFromDictionary(inputDict):
    """This function checks through the data and will remove duplicates from key 'epochtime's.
    
    A place holder to check, and remove duplicate times from this whole class. It needs to be though through still,
    but the code below is used to do it from an exterior script and would be a good place to start.

    Args:
        inputDict (dict): to check this if its duplicate

    Returns:
        inputdict (dict): same dictionary with-out duplicates in time

    References:
        https://www.peterbe.com/plog/fastest-way-to-uniquify-a-list-in-python-3.6
        https://stackoverflow.com/questions/480214/how-do-you-remove-duplicates-from-a-list
        -whilst-preserving-order

    """
    from collections.abc import Iterable  # this is not available in Python 2.7?
    if inputDict is not None:
        if 'epochtime' in inputDict:
            key = 'epochtime'
        elif 'time' in inputDict:
            key = 'time'
            warnings.warn(
                'Removing duplicates is faster using numeric time, failed looking for "epochtime" '
                'key')
        else:
            raise NotImplementedError('Requires keys "time" or "epochtime"')
        
        if isinstance(inputDict[key], Iterable) and np.size(set(np.array(inputDict[key]))) != np.size(inputDict[key]):
            # we have duplicate times in dictionary
            print(' Removing Duplicates from {}'.format(inputDict['name']))  # find the duplicates
            _, idxUnique = np.unique(inputDict[key], return_index=True)
            inputDict = sb.reduceDict(inputDict, idxUnique)
            # try:  # python 3.6 + only
            #    ans2 = list(dict.fromkeys(inputDict[key]))
            #    # nonzero(np.in1d(inputDict[key], ans2)) #<===================================
            #    this leaves duplicates
            # except: # a slower way
            #    seen = set()
            #    ans3b = [x for x in inputDict[key] if x not in seen and not seen.add(x)]
            #    idxObs = np.nonzero(np.in1d(inputDict[key], ans3b))
            # original Way  --- super slow
            # dupes = np.array([x for n, x in enumerate(inputDict[key]) if x in inputDict[key][
            # :n]]).squeeze()
            # idxObs = np.delete(np.arange(len(inputDict[key])),
            #                    np.argwhere(np.in1d(inputDict[key], dupes).squeeze())[
            #                    ::2].squeeze())  # delete
            #                    every other duplicate record
            # inputDict = sb.reduceDict(inputDict, idxObs)
    return inputDict


class getObs:
    """Class focused on retrieving observational data."""
    def __init__(self, d1, d2):
        """Data are returned in self.dataindex are inclusive at start, exclusive at end."""
        # this is active wave gauge list for looping through as needed
        self.waveGaugeList = ['waverider-26m', 'waverider-17m', 'awac-11m', '8m-array',
                              'awac-6m', 'awac-4.5m', 'adop-3.5m', 'xp200m', 'xp150m', 'xp125m']
        
        self.directionalWaveGaugeList = ['waverider-26m', 'waverider-17m', 'awac-11m', '8m-array',
                                         'awac-6m', 'awac-4.5m', 'adop-3.5m']
        
        self.currentsGaugeList = ['awac-11m', 'awac-6m', 'awac-4.5m', 'adop-3.5m']
        #self.rawdataloc_wave = []
        #self.outputdir = []  # location for outputfiles
        self.d1 = d1  # start date for data grab
        self.d2 = d2  # end data for data grab
        self.timeunits = 'seconds since 1970-01-01 00:00:00'
        self.epochd1 = nc.date2num(self.d1, self.timeunits)
        self.epochd2 = nc.date2num(self.d2, self.timeunits)
        self.callingClass = 'getObs'
        self.FRFdataloc = 'http://134.164.129.55/thredds/dodsC/FRF/'
        self.crunchDataLoc = 'http://134.164.129.55/thredds/dodsC/cmtb/'
        self.chlDataLoc = 'https://chlthredds.erdc.dren.mil/thredds/dodsC/frf/'  #

        self._comp_time()
        assert type(self.d2) == DT.datetime, 'd1 need to be in python "Datetime" data types'
        assert type(self.d1) == DT.datetime, 'd2 need to be in python "Datetime" data types'
    
    def _comp_time(self):
        """Test if times are backwards."""
        assert self.d2 >= self.d1, 'finish time: end needs to be after start time: start'
    
    def _roundtime(self, dt=None, roundto=60):
        """Round a datetime object to any time laps in seconds.
        
        Author: Thierry Husson 2012 - Use it as you want but don't blame me.

        Args:
          dt: datetime.datetime object, default now.
          roundto: Closest number of SECONDS to round to, default 1 minute

        Returns:
            datetime object that is rounded
            
        """
        if dt is None:
            dt = DT.datetime.now()
        seconds = (dt - dt.min).seconds
        # // is a floor division, not a comment on following line:
        rounding = (seconds + roundto / 2) // roundto * roundto
        return dt + DT.timedelta(0, rounding - seconds, -dt.microsecond)
    
    def getWaveSpec(self, gaugenumber=0, roundto=30, removeBadDataFlag=4, **kwargs):
        """This function pulls down the data from the thredds server and puts the data into a dictionary.
        
        TODO: Set optional date input from function arguments to change self.start self.end

        Args:
          gaugenumber: wave gauge numbers pulled from self.waveGaugeURLlookup
               see help on self.waveGaugeURLlookup for possible gauge names (Default value = 0)
          roundto: this is duration in minutes which data are expected.  times are rounded to nearest
             30 minute increment (data on server are not even times) (Default value = 30)
          removeBadDataFlag (int): this will remove data with a directional flag of 3/4 signaling questionable or
             failed directional spectra (default = 4, remove failed (directional) spectral data time periods)
             valid values: [3, 4, False]  False will not remove any data

        Keyword Args:
            "returnAB" (bool): if this is True function will return a's and b's for time period (default=False)
            "specOnly" (bool); if this is True function will not return bulk statistics (default=False)

        Returns:
          dictionary with following keys for all gauges
            'time' (array): time in datetime objects

            'epochtime' (array): time in epoch time

            'name' (str): gauge name

            'wavefreqbin' (array): wave frequencys associated with 2D spectra

            'wavedirbin' (array): wave direction bin associated with 2D spectra

            'xFRF' (float): x location in FRF coordinates

            'yFRF' (float): y location in FRF coordinates

            'lat' (float): latitude

            'lon' (float): longitude

            'depth' (float): nominal water dept

            'Hs' (array): wave height

            'peakf' (array): wave peak frequency

        """
        returnAB = kwargs.get('returnAB', False)
        specOnly = kwargs.get('specOnly', False)
        # Making gauges flexible
        self._waveGaugeURLlookup(gaugenumber)
        # parsing out data of interest in time
        self.ncfile, self.allEpoch = getnc(dataLoc=self.dataloc, callingClass=self.callingClass,
                                           dtRound=roundto * 60, start=self.d1, end=self.d2)
        try:
            self.wavedataindex = gettime(allEpoch=self.allEpoch, epochStart=self.epochd1,
                                         epochEnd=self.epochd2)
            assert np.array(self.wavedataindex).all() is not None, 'there''s no data in your time period'
            
            if np.size(self.wavedataindex) >= 1:
                # consistant for all wave gauges
                if np.size(self.wavedataindex) == 1:
                    self.wavedataindex = np.expand_dims(self.wavedataindex, axis=0)
                self.snaptime = nc.num2date(self.allEpoch[self.wavedataindex],
                                            self.ncfile['time'].units,
                                            only_use_cftime_datetimes=False)
                try:
                    depth = self.ncfile['nominalDepth'][:]  # this should always go
                except IndexError:
                    try:
                        depth = self.ncfile['gaugeDepth'][:]  # non directionalWaveGaugeList gauges
                    except IndexError:
                        depth = -999  # fill value
                try:
                    wave_coords = gp.FRFcoord(self.ncfile['longitude'][:],
                                              self.ncfile['latitude'][:])
                except IndexError:
                    wave_coords = gp.FRFcoord(self.ncfile['lon'][:], self.ncfile['lat'][:])
                #######################################################################################################
                # now that wave data index is resolved, go get data
                self.snaptime = nc.num2date(self.allEpoch[self.wavedataindex],
                                            self.ncfile['time'].units,
                                            only_use_cftime_datetimes=False)
                wavespec = {'time':        self.snaptime,  # note this is new variable names??
                            'epochtime':   self.allEpoch[self.wavedataindex],
                            'name':        str(self.ncfile.title),
                            'wavefreqbin': self.ncfile['waveFrequency'][:],
                            'xFRF':        wave_coords['xFRF'],
                            'yFRF':        wave_coords['yFRF'],
                            'lat':         self.ncfile['latitude'][:],
                            'lon':         self.ncfile['longitude'][:],
                            'depth':       depth,
                            'Hs':          self.ncfile['waveHs'][self.wavedataindex],
                            'peakf':       1 / self.ncfile['waveTp'][self.wavedataindex],
                }
                
                # now do directionalWaveGaugeList gauge try
                try:  # pull time specific data based on self.wavedataindex
                    wavespec['depth'] = self.ncfile['nominalDepth'][:]  # this should always go with directional gauges
                    wavespec['wavedirbin'] = self.ncfile['waveDirectionBins'][:]
                    wavespec['fspec'] = self.ncfile['waveEnergyDensity'][self.wavedataindex, :]
                    wavespec['qcFlagE'] = self.ncfile['qcFlagE'][self.wavedataindex]
                    wavespec['qcFlagD'] = self.ncfile['qcFlagD'][self.wavedataindex]
                    wavespec['dWED'] = self.ncfile['directionalWaveEnergyDensity'][self.wavedataindex, :, :]


                    if wavespec['dWED'].ndim < 3:
                        wavespec['dWED'] = np.expand_dims(wavespec['dWED'], axis=0)
                        wavespec['fspec'] = np.expand_dims(wavespec['fspec'], axis=0)

                    if specOnly is True:
                        return removeDuplicatesFromDictionary(wavespec)  # pull out here if specOnly is true (saves time)

                    wavespec['waveDp'] = self.ncfile['wavePeakDirectionPeakFrequency'][self.wavedataindex]
                    wavespec['waveDm'] = self.ncfile['waveMeanDirection'][self.wavedataindex]
                    wavespec['Tm'] = self.ncfile['waveTm'][self.wavedataindex]
                    
                    if returnAB is True:
                        wavespec['a1'] = self.ncfile['waveA1Value'][self.wavedataindex, :]
                        wavespec['a2'] = self.ncfile['waveA2Value'][self.wavedataindex, :]
                        wavespec['b1'] = self.ncfile['waveB1Value'][self.wavedataindex, :]
                        wavespec['b2'] = self.ncfile['waveB2Value'][self.wavedataindex, :]
                # this should throw when gauge is non directionalWaveGaugeList
                except IndexError:  # if error its non-directional gauge
                    # lidar guages don't have this variable.
                    if 'nominalDepth' in self.ncfile.variables.keys():
                        wavespec['depth'] = self.ncfile['nominalDepth'][:]  # non directional gauges
                    else:
                        # leave it blank if lidar wave gauge.
                        wavespec['depth'] = np.nan
                        
                    wavespec['wavedirbin'] = np.arange(0, 360, 90)  # 90 degree bins
                    wavespec['waveDp'] = np.ones_like(self.wavedataindex) * -999
                    try:
                        wavespec['fspec'] = self.ncfile['waveEnergyDensity'][self.wavedataindex, :]
                    except(RuntimeError):  # handle n-1 index error with Thredds
                        wavespec['fspec'] = self.ncfile['waveEnergyDensity'][self.wavedataindex[:-1], :]
                        wavespec['fspec'] = np.append(wavespec['fspec'],
                                                      self.ncfile['waveEnergyDensity'][
                                                      self.wavedataindex[-1], :][
                                                      np.newaxis, :], axis=0)
                    if wavespec['fspec'].ndim < 2:
                        wavespec['fspec'] = np.expand_dims(wavespec['fspec'], axis=0)
                    
                    # multiply the freq spectra for all directions
                    wavespec['dWED'] = np.ones([np.size(self.wavedataindex), np.size(wavespec['wavefreqbin']),
                                                np.size(wavespec['wavedirbin'])])
                    wavespec['dWED'] = wavespec['dWED']*wavespec['fspec'][:, :, np.newaxis]/len(wavespec['wavedirbin'])
                    if 'qcFlagE' in self.ncfile.variables.keys():
                        # lidar wave gauges don't have this variable.
                        wavespec['qcFlagE'] = self.ncfile['qcFlagE'][self.wavedataindex]
                    else:
                        # lidar wave gauges have waterLevelQCFlag and spectralQCFlag
                        wavespec['qcFlagE'] = self.ncfile['waterLevelQCFlag'][self.wavedataindex]
                if removeBadDataFlag is not False:
                    # Energy should not be needed
                    try:
                        # find data that are below bad data threshold
                        idx = np.argwhere(wavespec['qcFlagD'] < removeBadDataFlag).squeeze()
                        if np.size(idx) > 0:
                            wavespec = sb.reduceDict(wavespec,
                                                     idx)  # if there are values, keep good ones
                        idx = np.argwhere(wavespec['qcFlagE'] < removeBadDataFlag).squeeze()
                        if np.size(idx) > 0:
                            wavespec = sb.reduceDict(wavespec, idx)
                    except(KeyError):
                        pass  # non -directional gauge
                wavespec = removeDuplicatesFromDictionary(wavespec)
        
        except (RuntimeError, AssertionError):
            print('     ---- Problem Retrieving wave data from %s\n    - in this time period start: %s  End: %s' % (
                    gaugenumber, self.d1, self.d2))
            
            try:
                wavespec = {'lat':  self.ncfile['latitude'][:],
                            'lon':  self.ncfile['longitude'][:],
                            'name': str(self.ncfile.title), }
            except TypeError:  # when self.ncfile is None
                wavespec = None
            except KeyError:
                wavespec = {'lat':  self.ncfile['lat'][:],
                            'lon':  self.ncfile['lon'][:],
                            'name': str(self.ncfile.title), }
        
        return wavespec
    
    def getCurrents(self, gaugenumber=5, roundto=1):
        """This function pulls down the currents data from the Thredds Server.

        Args:
          gaugenumber: a string or number to get ocean currents from look up table

            gaugenumber = [2, 'awac-11m']

            gaugenumber = [3, 'awac-8m']

            gaugenumber = [4, 'awac-6m']

            gaugenumber = [5, 'awac-4.5m']

            gaugenumber = [6, 'adop-3.5m'] (Default value = 5)

          roundto: the time over which the wind record exists, ie data is collected in 10 minute
          increments
            data is rounded to the nearst [roundto] (default 1 min)

        Returns:
            dict, None if error is encoutered
                'name' (str): gauge name

                'time' (obj): date time objects time stamp

                'epochtime' (float): unix epoch time

                'aveU' (array): average cross-shore current from collection

                'aveV' (array): average along-shore current from collection

                'speed' (array): average speed  [m/s]

                'dir' (array):  current direction (TN)

                'lat' (float): latitude of gauge

                'lon' (float): longitude of gauge

                'xFRF' (float): cross-shore coordinate of gauge

                'yFRF' (float): along-shore coordinate of gauge

                'depth' (float): gauge nominal depth Depth is calculated by:
                            depth = -xducerD + blank + (binSize/2) + (numBins * binSize)

                'meanP' (array): mean pressure

        """
        assert gaugenumber.lower() in [2, 3, 4, 5, 6, 'awac-11m', 'awac-8m', 'awac-6m', 'awac-4.5m',
                                       'adop-3.5m'], 'Input string/number is not a valid gage ' \
                                                     'name/number'
        
        if gaugenumber in [2, 'awac-11m']:
            # gname = 'AWAC04 - 11m'
            self.dataloc = 'oceanography/currents/awac-11m/awac-11m.ncml'
        elif gaugenumber in [3, 'awac-8m']:
            # gname = 'AWAC 8m'
            self.dataloc = 'oceanography/currents/awac-8m/awac-8m.ncml'
        elif gaugenumber in [4, 'awac-6m']:
            # gname = 'AWAC 6m'
            self.dataloc = 'oceanography/currents/awac-6m/awac-6m.ncml'
        elif gaugenumber in [5, 'awac-4.5m']:
            # gname = 'AWAC 4.5m'
            self.dataloc = 'oceanography/currents/awac-4.5m/awac-4.5m.ncml'
        elif gaugenumber in [6, 'adop-3.5m']:
            # gname = 'Aquadopp 3.5m'
            self.dataloc = 'oceanography/currents/adop-3.5m/adop-3.5m.ncml'
        else:
            raise NameError('Check gauge name')
        
        self.ncfile, self.allEpoch = getnc(dataLoc=self.dataloc, callingClass=self.callingClass,
                                           dtRound=roundto * 60)  # start=self.d1, end=self.d2) <
        # -- needs to be tested
        currdataindex = gettime(allEpoch=self.allEpoch, epochStart=self.epochd1,
                                epochEnd=self.epochd2)
        
        # _______________________________________
        # get the actual current data
        if np.size(currdataindex) > 1:
            curr_aveU = self.ncfile['aveU'][
                currdataindex]  # pulling depth averaged Eastward current
            curr_aveV = self.ncfile['aveV'][
                currdataindex]  # pulling depth averaged Northward current
            curr_spd = self.ncfile['currentSpeed'][currdataindex]  # currents speed [m/s]
            curr_dir = self.ncfile['currentDirection'][
                currdataindex]  # current from direction [deg]
            self.curr_time = nc.num2date(self.allEpoch[currdataindex], self.ncfile['time'].units,
                                         self.ncfile['time'].calendar,
                                         only_use_cftime_datetimes=False)
            # for num in range(0, len(self.curr_time)):
            #     self.curr_time[num] = self.roundtime(self.curr_time[num], roundto=roundto * 60)
            
            curr_coords = gp.FRFcoord(self.ncfile['lon'][0], self.ncfile['lat'][0])
            
            self.curpacket = {
                'name':      str(self.ncfile.title),
                'time':      self.curr_time,
                'epochtime': self.allEpoch[currdataindex],
                'aveU':      curr_aveU,
                'aveV':      curr_aveV,
                'speed':     curr_spd,
                'dir':       curr_dir,
                'lat':       self.ncfile['lat'][0],
                'lon':       self.ncfile['lon'][0],
                'xFRF':      curr_coords['xFRF'],
                'yFRF':      curr_coords['yFRF'],
                'depth':     self.ncfile['depth'][:],
                # Depth is calculated by: depth = -xducerD + blank + (binSize/2) + (numBins *
                # binSize)
                'meanP':     self.ncfile['meanPressure'][currdataindex]}
            
            return self.curpacket
        
        else:
            print('ERROR: There is no current data for this time period!!!')
            self.curpacket = None
            return self.curpacket
    
    def getWind(self, gaugenumber=0, collectionlength=10):
        """This function retrieves the wind data.
        
        Collection length is the time over which the wind record exists ie data is collected in 10 minute increments
        data is rounded to the nearst [collectionlength] (default 10 min).

        Args:
          collectionlength: Default value = 10)
          gaugenumber: (Default value = 0)

            gauge number in ['derived', 'Derived', 0]

            '932 wind gauge' in [1]

            '832 wind gauge' in [2]

            '732 wind gauge' in [3]

        Returns:
            dict, will return None if an error is encountered
                'name' (str):  station name

                'time' (obj): datetime object time stamp

                'vecspeed' (array):  Vector Averaged Wind Speed

                'windspeed' (array): Mean Wind Speed

                'windspeed_corrected' (array): corrected 10m windspeed

                'winddir' (array):  Wind direction from true north

                'windgust' (array):  5 second largest mean wind speed

                'qcflagS' (array): QC flag for speed

                'qcflagD' (array): qcflag for directions

                'stdspeed' (array): std dev of 10 min wind record

                'minspeed' (array):  min speed in 10 min avg

                'maxspeed' (array):  max speed in 10 min avg

                'sustspeed' (array): 1 min largest mean wind speed

                'lat' (float):  latitude

                'lon' (float):  longitde

                'gaugeht' (float): gauge height for uncorrected wind measurements

        """
        # Making gauges flexible
        # different Gauges
        if gaugenumber in ['derived', 'Derived', 0]:
            self.dataloc = 'meteorology/wind/derived/derived.ncml'  # 932 wind gauge
            gname = 'Derived wind gauge '
        elif gaugenumber == 1:
            self.dataloc = 'meteorology/wind/D932/D932.ncml'  # 932 wind gauge
            gname = '932 wind gauge'
        elif gaugenumber == 2:
            gname = '832 wind gauge'
            self.dataloc = 'meteorology/wind/D832/D832.ncml'
        elif gaugenumber == 3:
            gname = '732 wind gauge'
            self.dataloc = 'meteorology/wind/D732/D732.ncml'
        else:
            raise NameError('Specifiy proper Gauge number')
        
        self.ncfile, self.allEpoch = getnc(dataLoc=self.dataloc, callingClass=self.callingClass,
                                           dtRound=collectionlength * 60)
        
        self.winddataindex = gettime(allEpoch=self.allEpoch, epochStart=self.epochd1,
                                     epochEnd=self.epochd2)
        # remove nan's that shouldn't be there
        # ______________________________________
        if np.size(self.winddataindex) > 0 and self.winddataindex is not None:
            # TODO: why are we removing nan's here. this should be resolved down stream if they're causing problems
            # do they even come up as nans? i thought they returned as masked arrays
            self.winddataindex = self.winddataindex[~np.isnan(self.ncfile['windDirection'][self.winddataindex])]
            if np.size(self.winddataindex) == 0:
                # return None is he wind direction is associated with the wind is no good!
                return None
            # MPG: moved inside if statement b/c call to gettime possibly returns None.
            self.winddataindex = self.winddataindex[~np.isnan(self.ncfile['windDirection'][self.winddataindex])]
            
            windvecspd = self.ncfile['vectorSpeed'][self.winddataindex]
            windspeed = self.ncfile['windSpeed'][self.winddataindex]  # wind speed
            winddir = self.ncfile['windDirection'][self.winddataindex]  # wind direction
            windgust = self.ncfile['windGust'][self.winddataindex]  # 5 sec largest mean speed
            stdspeed = self.ncfile['stdWindSpeed'][self.winddataindex]  # std dev of 10 min avg
            qcflagS = self.ncfile['qcFlagS'][self.winddataindex]  # qc flag
            qcflagD = self.ncfile['qcFlagD'][self.winddataindex]
            minspeed = self.ncfile['minWindSpeed'][self.winddataindex]  # min wind speed in 10 min avg
            maxspeed = self.ncfile['maxWindSpeed'][self.winddataindex]  # max wind speed in 10 min avg
            sustspeed = self.ncfile['sustWindSpeed'][self.winddataindex]  # 1 minute largest mean wind speed
            gaugeht = self.ncfile.geospatial_vertical_max
            
            self.windtime = nc.num2date(self.allEpoch[self.winddataindex],
                                        self.ncfile['time'].units,
                                        only_use_cftime_datetimes=False)
            
            # correcting for wind elevations from Johnson (1999) - Simple Expressions for
            # correcting wind speed data
            # for elevation
            if gaugeht <= 20:
                windspeed_corrected = windspeed * (10 / gaugeht) ** (1 / 7)
            else:
                windspeed_corrected = 'No Corrections done for gauges over 20m, please read: ' \
                                      '\nJohnson (1999) - ' \
                                      'Simple Expressions for correcting wind speed data for ' \
                                      'elevation'
            windpacket = {
                'name':                str(self.ncfile.title),  # station name
                'time':                self.windtime,  # time
                'epochtime':           self.allEpoch[self.winddataindex],
                'vecspeed':            windvecspd,  # Vector Averaged Wind Speed
                'windspeed':           windspeed,  # Mean Wind Speed
                'windspeed_corrected': windspeed_corrected,  # corrected windspeed
                'winddir':             winddir,  # Wind direction from true nort
                'windgust':            windgust,  # 5 second largest mean wind speed
                'qcflagS':             qcflagS,  # QC flag
                'qcflagD':             qcflagD,
                'stdspeed':            stdspeed,  # std dev of 10 min wind record
                'minspeed':            minspeed,  # min speed in 10 min avg
                'maxspeed':            maxspeed,  # max speed in 10 min avg
                'sustspeed':           sustspeed,  # 1 min largest mean wind speed
                'lat':                 self.ncfile['latitude'][:],  # latitude
                'lon':                 self.ncfile['longitude'][:],  # longitde
                'gaugeht':             gaugeht,
                }
            if (windpacket['qcflagD'] == 3).all() or (windpacket['qcflagS'] == 3).all():
                print("Wind querey returned all bad data for speed or direction")
                windpacket = None
            return windpacket
        else:
            print('     ---- ERROR: Problem finding wind !!!')
            return None
    
    def getWL(self, collectionlength=6):
        """This function retrieves the water level data from the server.
        
        WL data on server is NAVD88
        
        collection length is the time over which the wind record exists
            ie data is collected in 10 minute increments
            data is rounded to the nearst [collectionlength] (default 6 min)

        Args:
          collectionlength (int): dictates what value to round time to (Default value = 6)

        Returns:
          dictionary with keys
            'name': gauge name - taken from title

            'WL': measured water level (NAVD88) [m]

            'time': datetime object

            'epochtime': time in seconds since 1970-01-01 (float)

            'lat': latitude

            'lon':  longitude

            'residual': water level residual

            'predictedWL': predicted tide

        """
        # this is the back end of the url for waterlevel
        self.dataloc = 'oceanography/waterlevel/eopNoaaTide/eopNoaaTide.ncml'

        self.ncfile, self.allEpoch = getnc(dataLoc=self.dataloc, callingClass=self.callingClass,
                                           dtRound=collectionlength * 60, start=self.d1,
                                           end=self.d2)
        self.WLdataindex = gettime(allEpoch=self.allEpoch, epochStart=self.epochd1,
                                   epochEnd=self.epochd2)
        
        if np.size(self.WLdataindex) > 1:
            self.WLtime = nc.num2date(self.allEpoch[self.WLdataindex], self.ncfile['time'].units,
                                      only_use_cftime_datetimes=False)
            self.WLpacket = {
                'name':        str(self.ncfile.title),
                'WL':          self.ncfile['waterLevel'][self.WLdataindex],
                # why does this call take so long for even 10 data points?
                'time':        self.WLtime,
                'epochtime':   self.allEpoch[self.WLdataindex],
                'lat':         self.ncfile['latitude'][:],
                'lon':         self.ncfile['longitude'][:],
                'predictedWL': self.ncfile['predictedWaterLevel'][self.WLdataindex],
                    }
            # this is faster to calculate myself, than pull from server
            self.WLpacket['residual'] = self.WLpacket['WL'] - self.WLpacket['predictedWL']
        elif self.WLdataindex is not None and np.size(self.WLdataindex) == 1:
            raise BaseException('you have 1 WL point, can the above be a >= logic or does 1 cause problems')
        else:
            print('ERROR: there is no WATER level Data for this time period!!!')
            self.WLpacket = None
        return self.WLpacket
    
    def getGaugeWL(self, gaugenumber=5, roundto=1):
        """This function pulls down the water level data at a particular gauge from the Server.
        
        Args:
            gaugenumber (int/str) describing the location (default=5 End of pier)
            roundto: the time over which the wind record exists ie data is collected in 10 minute
            increments
                        data is rounded to the nearst [roundto] (default 1 min)

        Returns
            wlpacket (dict) with keys below
                'name': gagename

                'time': datetime of the measurements

                'epochtime': epochtime of the measurements

                'wl': water level at the gage (NAVD88?)

                'lat': latitude of the gage

                'lon': longitude of the gage

                'xFRF': xFRF position of the gage

                'yFRF': yFRF position of the gage

        """
        # Making gauges flexible
        self._wlGageURLlookup(gaugenumber)
        # parsing out data of interest in time
        
        self.ncfile, self.allEpoch = getnc(dataLoc=self.dataloc, callingClass=self.callingClass,
                                           dtRound=roundto * 60)
        
        try:
            self.wldataindex = gettime(allEpoch=self.allEpoch, epochStart=self.epochd1,
                                       epochEnd=self.epochd2)
            assert np.array(self.wldataindex).all() != None, 'there''s no data in your time period'
            if np.size(self.wldataindex) >= 1:
                # consistant for all wl gauges
                
                # do we need this?
                # it is causing some of the stuff further down to crash
                # if you have only one data point in the date range?
                # if np.size(self.wldataindex) == 1:
                # self.wldataindex = np.expand_dims(self.wldataindex, axis=0)
                
                self.snaptime = nc.num2date(self.allEpoch[self.wldataindex],
                                            self.ncfile['time'].units,
                                            only_use_cftime_datetimes=False)
                try:
                    wl_coords = gp.FRFcoord(self.ncfile['longitude'][:], self.ncfile['latitude'][:])
                except IndexError:
                    wl_coords = gp.FRFcoord(self.ncfile['lon'][:], self.ncfile['lat'][:])
                
                wlpacket = {'time':      self.snaptime,  # note this is new variable names??
                            'epochtime': self.allEpoch[self.wldataindex],
                            'name':      str(self.ncfile.title),
                            'xFRF':      wl_coords['xFRF'],
                            'yFRF':      wl_coords['yFRF'],
                            'lat':       self.ncfile['latitude'][:],
                            'lon':       self.ncfile['longitude'][:],
                            'wl':        self.ncfile['waterLevel'][self.wldataindex]}
                return wlpacket
        
        except (RuntimeError, AssertionError):
            print(
                '     ---- Problem Retrieving water level data from %s\n    - in this time period '
                'start: %s  End: %s'
                % (gaugenumber, self.d1, self.d2))
            try:
                wlpacket = {'lat':  self.ncfile['latitude'][:],
                            'lon':  self.ncfile['longitude'][:],
                            'name': str(self.ncfile.title), }
            except:
                wlpacket = {'lat':  self.ncfile['lat'][:],
                            'lon':  self.ncfile['lon'][:],
                            'name': str(self.ncfile.title), }
            return wlpacket
    
    def getBathyTransectFromNC(self, profilenumbers=None, method=1, forceReturnAll=False):
        """This function gets the bathymetric data from the server.

        Args:
          profilenumbers: Default value = None)
          method: bathymetry selection method (Default value = 1)
               method == 1  - > 'Bathymetry is taken as closest in HISTORY - operational'

               method == 0  - > 'Bathymetry is taken as closest in TIME - NON-operational'
          forceReturnAll (bool): (Default Value = False)
                This will force the survey to take and return all indices between start and end,
                not the single
        Returns:
          dictionary with keys, will return None if call fails
            'xFRF': x coordinate in frf

            'yFRF': y coordiante in Frf

            'elevation': bathy elevation

            'time': time in date time object

            'lat': lat,

            'lon': lon,

            'northing': NC northing

            'easting': NC easting

            'profileNumber': FRF profile number

            'surveyNumber': FRF survey Number

            'Ellipsoid': which ellipsoid is used

        """
        # do check here on profile numbers
        # acceptableProfileNumbers = [None, ]
        self.dataloc = 'geomorphology/elevationTransects/survey/surveyTransects.ncml'  # location
        # of the gridded surveys
        self.ncfile, self.allEpoch = getnc(dataLoc=self.dataloc, callingClass=self.callingClass,
                                           dtRound=1 * 60)
        try:
            self.bathydataindex = gettime(allEpoch=self.allEpoch, epochStart=self.epochd1,
                                          epochEnd=self.epochd2)
        except IOError:  # when data are not on CHL thredds
            self.bathydataindex = None
        # returning None object is convention and must be followed/handled down the line
        # if self.bathydataindex is None:
        #     self.bathydataindex = []
        
        # logic to handle no transects in date range
        if forceReturnAll == True:
            idx = self.bathydataindex
        elif np.size(self.bathydataindex) == 1 and self.bathydataindex is not None:
            idx = self.bathydataindex
        elif (np.size(self.bathydataindex) < 1 & method == 1) or (
            self.bathydataindex is None and method == 1):
            # there's no exact bathy match so find the max negative number where the negative
            # numbers are historical and the max would be the closest historical
            val = (max([n for n in (self.ncfile['time'][:] - self.epochd1) if n < 0]))
            idx = np.where((self.ncfile['time'][:] - self.epochd1) == val)[0][0]
            # print 'Bathymetry is taken as closest in HISTORY - operational'
        elif (np.size(self.bathydataindex) < 1 and method == 0) or (
            self.bathydataindex is None and method == 1):
            idx = np.argmin(np.abs(self.ncfile['time'][:] - self.d1))  # closest in time
            # print 'Bathymetry is taken as closest in TIME - NON-operational'
        elif np.size(self.bathydataindex) > 1:  # if dates fall into d1,d2 bounds,
            idx = self.bathydataindex[
                0]  # return a single index. this means there was a survey between d1,d2
        
        if forceReturnAll is not True:
            # find the whole survey (via surveyNumber) and assign idx to return the whole survey
            idxSingle = idx
            idx = np.argwhere(
                self.ncfile['surveyNumber'][:] == self.ncfile['surveyNumber'][idxSingle]).squeeze()
            if np.size(idx) == 0:
                print('The closest in history to your start date is %s\n' % nc.num2date(
                    self.gridTime[idx],
                    self.ncfile['time'].units))
                raise NotImplementedError('Please End new simulation with the date above')
                idx = self.bathydataindex
        
        # else:
        #     # Now that indices of interest are sectioned off, find the survey number that
        #     matches them and return
        #     whole survey
        #     idxSingle = idx
        #     # idx = np.argwhere(self.ncfile['surveyNumber'][:] == self.ncfile['surveyNumber'][
        #     idxSingle]).squeeze()
        # isolate specific profile numbers if necessicary
        if profilenumbers != None:
            assert pd.Series(profilenumbers).isin(np.unique(self.ncfile['profileNumber'][
                                                                idx])).all(), 'given profiles ' \
                                                                              'don''t Match ' \
                                                                              'profiles ' \
                                                                              'in database'  # if
            # all of the profile
            # numbers match
            idx2mask = np.in1d(self.ncfile['profileNumber'][idx],
                               profilenumbers)  # boolean true/false of time and profile number
            idx = idx[idx2mask]
        # elif pd.Series(profileNumbers).isin(np.unique(self.cshore_ncfile['profileNumber'][
        # :])).any(): #if only some
        # of the profile numbers match
        #     print 'One or more input profile numbers do not match those in the FRF transects!
        #     Fetching data for
        #     those that do.'
        #     mask = (self.alltime >= self.start) & (self.alltime < self.end) & np.in1d(
        #     self.cshore_ncfile[
        #     'profileNumber'][:],profileNumbers)  # boolean true/false of time and profile number
        
        # now retrieve data with idx
        if np.size(idx) > 0 and idx is not None:
            elevation_points = self.ncfile['elevation'][idx]
            xCoord = self.ncfile['xFRF'][idx]
            yCoord = self.ncfile['yFRF'][idx]
            lat = self.ncfile['lat'][idx]
            lon = self.ncfile['lon'][idx]
            northing = self.ncfile['northing'][idx]
            easting = self.ncfile['easting'][idx]
            profileNum = self.ncfile['profileNumber'][idx]
            surveyNum = self.ncfile['surveyNumber'][idx]
            Ellipsoid = self.ncfile['Ellipsoid'][idx]
            time = nc.num2date(self.ncfile['time'][idx], self.ncfile['time'].units,
                               only_use_cftime_datetimes=False)
            
            profileDict = {'xFRF':          xCoord,
                           'yFRF':          yCoord,
                           'elevation':     elevation_points,
                           'epochtime':     self.allEpoch[idx],
                           'time':          time,
                           'lat':           lat,
                           'lon':           lon,
                           'northing':      northing,
                           'easting':       easting,
                           'profileNumber': profileNum,
                           'surveyNumber':  surveyNum,
                           'Ellipsoid':     Ellipsoid, }
        
        else:
            profileDict = None
        
        return profileDict
    
    def getBathyTransectProfNum(self, method=1):
        """This function gets the bathymetric data from the server, currently designed for the bathy duck experiment.
        
        Just gets profile numbers only.

        Args:
            method (int): approach to select which method of how to select bathymetry
                method == 1  - > 'Bathymetry is taken as closest in HISTORY - operational'

                method == 0  - > 'Bathymetry is taken as closest in TIME - NON-operational' (
                Default value = 1)

        Returns:
            prof_nums (array): an array of profile numbers over which a single survey was taken

        """
        # do check here on profile numbers
        # acceptableProfileNumbers = [None, ]
        self.dataloc = 'geomorphology/elevationTransects/survey/surveyTransects.ncml'  # location
        # of the gridded surveys
        self.ncfile, self.allEpoch = getnc(dataLoc=self.dataloc, callingClass=self.callingClass,
                                           dtRound=1 * 60)
        
        try:
            self.bathydataindex = gettime(allEpoch=self.allEpoch, epochStart=self.epochd1,
                                          epochEnd=self.epochd2)
        except IOError:  # when data are not on CHL thredds
            self.bathydataindex = []
        
        if self.bathydataindex is None:
            self.bathydataindex = []
        else:
            pass
        
        # logic to handle no transects in date range
        if len(self.bathydataindex) == 1:
            idx = self.bathydataindex
        elif len(self.bathydataindex) < 1 & method == 1:
            # there's no exact bathy match so find the max negative number where the negative
            # numbers are historical and the max would be the closest historical
            val = (max([n for n in (self.ncfile['time'][:] - self.epochd1) if n < 0]))
            idx = np.where((self.ncfile['time'][:] - self.epochd1) == val)[0][0]
            print('Bathymetry is taken as closest in HISTORY - operational')
        
        elif len(self.bathydataindex) < 1 and method == 0:
            idx = np.argmin(np.abs(self.ncfile['time'][:] - self.d1))  # closest in time
            print('Bathymetry is taken as closest in TIME - NON-operational')
        
        elif len(self.bathydataindex) > 1:
            try:
                # switch back to the FRF cshore_ncfile?
                self.ncfile = nc.Dataset(self.FRFdataloc + self.dataloc)
            except:
                pass
            raise NotImplementedError('DLY NOTE')
            
            # DLY Note - this section of the script does NOT work
            # (bb.e., if you DO have a survey during your date range!!!)
            timeunits = 'seconds since 1970-01-01 00:00:00'
            d1Epoch = nc.date2num(self.d1, timeunits)
            val = (max([n for n in (self.ncfile['time'][:] - d1Epoch) if n < 0]))
            idx = np.where((self.ncfile['time'][:] - d1Epoch) == val)[0][0]
        
        # returning whole survey
        idxSingle = idx
        idx = np.argwhere(
            self.ncfile['surveyNumber'][:] == self.ncfile['surveyNumber'][idxSingle]).squeeze()
        
        # what profile numbers are in this survey?
        prof_nums = np.unique(self.ncfile['profileNumber'][idx])
        
        return prof_nums
    
    def getBathyGridFromNC(self, method, removeMask=True):
        """This function gets the frf krigged grid product, it will currently break with the present link.
        
        Bathymetric data from the server.

        Args:
          method: defines which choice method to use
             method == 1  - > 'Bathymetry is taken as closest in HISTORY - operational'

             method == 0  - > 'Bathymetry is taken as closest in TIME - NON-operational'

          removeMask (bool): remove data that are masked (Default value = True)

        Returns:
            gridDict (dict): colleciton of variables with keys below, will return None if an
            error occurs
               'xFRF': xCoord,

               'yFRF': yCoord,

               'elevation': elevation_points,

               'time': time,

               'lat': lat,

               'lon': lon,

               'northing': northing,

               'easting': easting
               
        """
        self.dataloc = 'survey/gridded/gridded.ncml'  # location of the gridded surveys
        self.ncfile, self.allEpoch = getnc(dataLoc=self.dataloc, callingClass=self.callingClass,
                                           dtRound=1 * 60)
        try:
            self.bathydataindex = gettime(allEpoch=self.allEpoch, epochStart=self.epochd1,
                                          epochEnd=self.epochd2)  # getting the index of the grid
        except IOError:
            self.bathydataindex = []
        if self.bathydataindex != None and len(self.bathydataindex) == 1:
            idx = self.bathydataindex
        elif (self.bathydataindex == None or len(self.bathydataindex) < 1) & method == 1:
            # there's no exact bathy match so find the max negative number where the negitive
            # numbers are historical and the max would be the closest historical
            val = (max([n for n in (self.ncfile['time'][:] - self.epochd1) if n < 0]))
            idx = np.where((self.ncfile['time'][:] - self.epochd1) == val)[0][0]
            print('Bathymetry is taken as closest in HISTORY - operational')
        elif (self.bathydataindex == None or len(self.bathydataindex) < 1) and method == 0:
            idx = np.argmin(np.abs(self.ncfile['time'][:] - self.d1))  # closest in time
            print('Bathymetry is taken as closest in TIME - NON-operational')
        elif self.bathydataindex != None and len(self.bathydataindex) > 1:
            val = (max([n for n in (self.ncfile['time'][:] - self.d1) if n < 0]))
            idx = np.where((self.ncfile['time'] - self.d1) == val)[0][0]
            
            print('The closest in history to your start date is %s\n' % nc.num2date(
                self.gridTime[idx],
                self.ncfile['time'].units))
            print('Please End new simulation with the date above')
            raise Exception
        # the below line was in place, it should be masking nan's but there is not supposed to be
        # nan's
        # in the data, should only be fill values (-999)
        # elevation_points = np.ma.array(cshore_ncfile['elevation'][idx,:,:], mask=np.isnan(
        # cshore_ncfile[
        # 'elevation'][idx,:,:]))
        # remove -999's
        elevation_points = self.ncfile['elevation'][idx, :, :]
        if type(elevation_points) != np.ma.core.MaskedArray:
            maskedElev = (elevation_points == self.ncfile['elevation']._FillValue)
            elevation_points = np.ma.array(elevation_points, mask=maskedElev)
        if elevation_points.ndim == 3:
            elevation_points = elevation_points[0]
        xCoord = self.ncfile['xFRF'][:]
        yCoord = self.ncfile['yFRF'][:]
        lat = self.ncfile['latitude'][:]
        lon = self.ncfile['longitude'][:]
        northing = self.ncfile['northing'][:]
        easting = self.ncfile['easting'][:]
        if removeMask == True:
            xCoord = xCoord[~np.all(elevation_points.mask, axis=0)]
            yCoord = yCoord[~np.all(elevation_points.mask, axis=1)]
            lon = lon[~np.all(elevation_points.mask, axis=0), :]
            lat = lat[:, ~np.all(elevation_points.mask, axis=1)]
            northing = northing[~np.all(elevation_points.mask, axis=0), :]
            easting = easting[:, ~np.all(elevation_points.mask, axis=1)]
            elevation_points = elevation_points[~np.all(elevation_points.mask, axis=1), :]  #
            elevation_points = elevation_points[:, ~np.all(elevation_points.mask, axis=0)]
        if elevation_points.ndim == 2:
            elevation_points = np.ma.expand_dims(elevation_points, axis=0)
        
        time = (self.ncfile['time'][idx], self.ncfile['time'].units)
        print('Sim start: %s\nSim End: %s\nSim bathy chosen: %s' % (self.d1, self.d2,
                                                                    nc.num2date(
                                                                        self.ncfile['time'][idx],
                                                                        self.ncfile['time'].units)))
        print('Bathy is %s old' % (
            self.d2 - nc.num2date(self.ncfile['time'][idx], self.ncfile['time'].units,
                                  only_use_cftime_datetimes=False)))
        
        gridDict = {'xFRF':      xCoord,
                    'yFRF':      yCoord,
                    'elevation': elevation_points,
                    'time':      time,
                    'lat':       lat,
                    'lon':       lon,
                    'northing':  northing,
                    'easting':   easting
                    }
        return gridDict
    
    def _waveGaugeURLlookup(self, gaugenumber):
        """A lookup table function that sets the URL backend for get wave spec and get wave gauge locations.

        Args:
          gaugenumber: a string or number that refers to a specific gauge and will set a url
                Available values inclue:

                26m waverider    can be [0, 'waverider-26m', 'Waverider-26m', '26m']

                17m waverider    can be [1, 'Waverider-17m', 'waverider-17m']

                11m AWAC         can be [2, 'AWAC-11m', 'awac-11m', 'Awac-11m']

                8m AWAC          can be [3, 'awac-8m', 'AWAC-8m']

                6m AWAC          can be [4, 'awac-6m', 'AWAC-6m']

                4.5m AWAC        can be [5, 'awac-4.5m', 'Awac-4.5m']

                3.5m aquadopp    can be [6, 'adop-3.5m', 'aquadopp 3.5m']

                200m pressure    can be [8, 'xp200m', 'xp200']

                150m pressure    can be [9, 'xp150m', 'xp150']

                125m pressure    can be [10, 'xp125m', 'xp125']

                100m pressure    can be [11, 'xp100m']

                8m array         can be [8, '8m-Array', '8m Array', '8m array', '8m-array']

                oregon inlet WR  can be ['oregonInlet', 'OI', 'oi']

        Returns:
          Nothing, this just sets the self.dataloc data member

        """
        if str(gaugenumber).lower() in ['0', 'waverider-26m', '26m']:
            # 26 m wave rider
            self.dataloc = 'oceanography/waves/waverider-26m/waverider-26m.ncml'  #
            # 'oceanography/waves/waverider430/waverider430.ncml'  # 26m buoy
        elif str(gaugenumber).lower() in ['1', 'waverider-17m', '17m']:
            # 2D 17m waverider
            self.dataloc = 'oceanography/waves/waverider-17m/waverider-17m.ncml'  # 17 m buoy
        elif str(gaugenumber).lower() in ['2', 'awac-11m', '11m']:
            self.dataloc = 'oceanography/waves/awac-11m/awac-11m.ncml'
        elif str(gaugenumber).lower() in ['3', 'awac-8m']:
            self.dataloc = 'oceanography/waves/awac-8m/awac-8m.ncml'
        elif str(gaugenumber).lower() in ['4', 'awac-6m', 'awac 6m']:
            self.dataloc = 'oceanography/waves/awac-6m/awac-6m.ncml'
        elif str(gaugenumber).lower() in ['5', 'awac-4.5m', 'awac_4.5m']:
            self.dataloc = 'oceanography/waves/awac-4.5m/awac-4.5m.ncml'
        elif str(gaugenumber).lower() in ['6', 'adop-3.5m', 'aquadopp 3.5m']:
            self.dataloc = 'oceanography/waves/adop-3.5m/adop-3.5m.ncml'
        elif str(gaugenumber).lower() in ['7', 'adop-2m']:
            self.dataloc = 'oceanography/waves/adop01/adop01.ncml'
        elif str(gaugenumber).lower() in ['8', 'xp200m', 'xp200']:
            self.dataloc = 'oceanography/waves/xp200m/xp200m.ncml'
        elif str(gaugenumber).lower() in ['9', 'xp150m', 'xp150']:
            self.dataloc = 'oceanography/waves/xp150m/xp150m.ncml'
        elif str(gaugenumber).lower() in ['10', 'xp125m', 'xp125']:
            self.dataloc = 'oceanography/waves/xp125m/xp125m.ncml'
        elif str(gaugenumber).lower() in ['11', 'xp100m']:
            self.dataloc = 'oceanography/waves/xp100m/xp100m.ncml'
        elif str(gaugenumber).lower() in ['12', '8m array', '8m-array']:
            self.dataloc = 'oceanography/waves/8m-array/8m-array.ncml'
        # lidar wave gauges - 140 m
        elif str(gaugenumber).lower() in ['lidarwavegauge140', 'lidargauge140',
                                          'lidarwavegauge140m', 'lidargauge140m']:
            self.dataloc = 'oceanography/waves/lidarWaveGauge140/lidarWaveGauge140.ncml'
        # lidar wave gauges - 110 m
        elif str(gaugenumber).lower() in ['lidarwavegauge110', 'lidargauge110',
                                          'lidarwavegauge110m', 'lidargauge110m']:
            self.dataloc = 'oceanography/waves/lidarWaveGauge110/lidarWaveGauge110.ncml'
        # lidar wave gauges - 100 m
        elif str(gaugenumber).lower() in ['lidarwavegauge100', 'lidargauge100',
                                          'lidarwavegauge100m', 'lidargauge100m']:
            self.dataloc = 'oceanography/waves/lidarWaveGauge100/lidarWaveGauge100.ncml'
        # lidar wave gauges - 90 m
        elif str(gaugenumber).lower() in ['lidarwavegauge90', 'lidargauge90', 'lidarwavegauge90m',
                                          'lidargauge90m']:
            self.dataloc = 'oceanography/waves/lidarWaveGauge90/lidarWaveGauge90.ncml'
        # lidar wave gauges - 80 m
        elif str(gaugenumber).lower() in ['lidarwavegauge80', 'lidargauge80', 'lidarwavegauge80m',
                                          'lidargauge80m']:
            self.dataloc = 'oceanography/waves/lidarWaveGauge80/lidarWaveGauge80.ncml'
        elif str(gaugenumber).lower() in ['oregoninlet', 'oi']:
            self.dataloc = 'oceanography/waves/waverider-oregon-inlet-nc/waverider-oregon-inlet' \
                           '-nc.ncml'
        elif str(gaugenumber).lower() in ['lidarwavegauge080', 'lidarwavegauge80']:
            self.dataloc = "oceanography/waves/lidarWaveGauge080/lidarWaveGauge080.ncml"
        elif str(gaugenumber).lower() in ['lidarwavegauge090', 'lidarwavegauge90']:
            self.dataloc = "oceanography/waves/lidarWaveGauge090/lidarWaveGauge090.ncml"
        elif str(gaugenumber).lower() in ['lidarwavegauge100']:
            self.dataloc = "oceanography/waves/lidarWaveGauge100/lidarWaveGauge100.ncml"
        elif str(gaugenumber).lower() in ['lidarwavegauge110']:
            self.dataloc = "oceanography/waves/lidarWaveGauge110/lidarWaveGauge110.ncml"
        elif str(gaugenumber).lower() in ['lidarwavegauge140']:
            self.dataloc = "oceanography/waves/lidarWaveGauge140/lidarWaveGauge140.ncml"
        else:
            self.gname = 'There Are no Gauge numbers here'
            raise NameError('Bad Gauge name, specify proper gauge name/number, or add capability')
    
    def _wlGageURLlookup(self, gaugenumber):
        """A lookup table function that sets the URL backend for getGageWL.
        
        Args:
            gaugenumber: a string or number that refers to a specific gauge and will set a url
               Available values include:
                   11m AWAC         can be [2, 'AWAC-11m', 'awac-11m', 'Awac-11m']
                   8m AWAC          can be [3, 'awac-8m', 'AWAC-8m']
                   6m AWAC          can be [4, 'awac-6m', 'AWAC-6m']
                   4.5m AWAC        can be [5, 'awac-4.5m', 'Awac-4.5m']
                   3.5m aquadopp    can be [6, 'adop-3.5m', 'aquadopp 3.5m']
                   200m pressure    can be [8, 'xp200m', 'xp200']
                   150m pressure    can be [9, 'xp150m', 'xp150']
                   125m pressure    can be [10, 'xp125m', 'xp125']
                   100m pressure    can be [11, 'xp100m']
                   8m array         can be [8, '8m-Array', '8m Array', '8m array', '8m-array']
        
        Returns:
             Nothing, this just sets the self.dataloc data member

        """
        if gaugenumber in [2, 'AWAC-11m', 'awac-11m', 'Awac-11m']:
            gname = 'AWAC 11m'
            self.dataloc = 'oceanography/waves/awac-11m/awac-11m.ncml'
        elif gaugenumber in [3, 'awac-8m', 'AWAC-8m']:
            self.gname = 'AWAC 8m'
            self.dataloc = 'oceanography/waves/awac-8m/awac-8m.ncml'
        elif gaugenumber in [4, 'awac-6m', 'AWAC-6m']:
            self.gname = 'AWAC 6m'
            self.dataloc = 'oceanography/waves/awac-6m/awac-6m.ncml'
        elif gaugenumber in [5, 'awac-4.5m', 'Awac-4.5m', 'awac_4.5m']:
            self.gname = 'AWAC 4.5m'
            self.dataloc = 'oceanography/waves/awac-4.5m/awac-4.5m.ncml'
        elif gaugenumber in [6, 'adop-3.5m', 'aquadopp 3.5m']:
            self.gname = 'Aquadopp 3.5m'
            self.dataloc = 'oceanography/waves/adop-3.5m/adop-3.5m.ncml'
        elif gaugenumber in [7, 'adop-2m']:
            self.gname = 'Aquadopp01 - 2m'
            self.dataloc = 'oceanography/waves/adop01/adop01.ncml'
        elif gaugenumber in [8, 'xp200m', 'xp200']:
            self.gname = 'Paros xp200m'
            self.dataloc = 'oceanography/waves/xp200m/xp200m.ncml'
        elif gaugenumber in [9, 'xp150m', 'xp150']:
            self.gname = 'Paros xp150m'
            self.dataloc = 'oceanography/waves/xp150m/xp150m.ncml'
        elif gaugenumber in [10, 'xp125m', 'xp125']:
            self.gname = 'Paros xp125m'
            self.dataloc = 'oceanography/waves/xp125m/xp125m.ncml'
        elif gaugenumber in [11, 'xp100m']:
            self.gname = 'Paros xp100m'
            self.dataloc = 'oceanography/waves/xp100m/xp100m.ncml'
        elif gaugenumber in [12, '8m-Array', '8m Array', '8m array', '8m-array']:
            self.gname = "8m array"
            self.dataloc = 'oceanography/waves/8m-array/8m-array.ncml'
        else:
            self.gname = 'There Are no Gauge numbers here'
            raise NameError('Bad Gauge name, specify proper gauge name/number')
    
    def getBathyDuckLoc(self, gaugenumber):
        """This function pulls the stateplane location (if desired) from the survey.
        
        FRF coords from deployed ADV's, These are data owned by WHOI and kept on private local server

        Args:
          gaugenumber: a gauge number with associated data from bathyduck experiment

        Returns:
          location dictionary with keys
            'xFRF' (int): cross-shore coordinate of gauge

            'yFRF' (int): alongshore coordinate of gauge

        Notes:
            these data are only available by request

        """
        if type(gaugenumber) != str:
            gaugenumber = str(gaugenumber)
        assert gaugenumber in ['11', '12', '13', '14',  # middle transect
                               '21', '22', '23', '24',  # south transect
                               '83', '84'], 'gauge number not recognized'  # near pier
        try:
            loc = str(
                self.FRFdataloc + "projects/bathyduck/data/BathyDuck-ocean_waves_p%s_201510.nc" %
                gaugenumber)
            ncfile = nc.Dataset(loc)
            xloc = ncfile['xFRF'][:]
            yloc = ncfile['yFRF'][:]
        except:
            loc = str(
                self.chlDataLoc + 'projects/bathyduck/data/BathyDuck-ocean_waves_p%s_201510.nc' %
                gaugenumber)
            ncfile = nc.Dataset(loc)
            xloc = ncfile['xloc'][
                   :]  # these are hard coded in these files [do not change w/o recreating the file]
            yloc = ncfile['yloc'][:]
        assert len(np.unique(xloc)) == 1, "there are different locations in the netCDFfile"
        assert len(np.unique(yloc)) == 1, "There are different Y locations in the NetCDF file"
        locDict = gp.FRFcoord(xloc[0], yloc[0])
        
        return locDict
    
    def getWaveGaugeLoc(self, gaugenumber):
        """This function gets gauge location data quickly, faster than getwavespec.

        Args:
          gaugenumber (str, int): wave gauge numbers
                pulled from self.waveGaugeURLlookup
        
        Returns:
          dictionary with keys
            lat: latitude
            lon: longitude

        Notes:
            see help on self.waveGaugeURLlookup for gauge keys
        """
        self._waveGaugeURLlookup(gaugenumber)
        try:
            ncfile = nc.Dataset(self.FRFdataloc + self.dataloc)
        except IOError:
            ncfile = nc.Dataset(self.chlDataLoc + self.dataloc)
        out = {'Lat': ncfile['latitude'][:],
               'Lon': ncfile['longitude'][:]}
        return out
    
    def get_sensor_locations_from_thredds(self):
        """Retrieves lat/lon coordinates for each gauge in gauge_list.
        
        Function converts to state plane and frf coordinates and creates a dictionary containing
        all three coordinates types with gaugenumbers as keys.

        Args:
            None

        Returns:
            loc_dict (dict): Dictionary containing lat/lon, state plane, and frf coordinates
                for each available gaugenumber with gaugenumbers as keys.

              'lat': latitude

              'lon': longitude

              'spE': North Carolina StatePlane easting

              'spN': North Carolina StatePlane Northing

              'xFRF': FRF local coordinate system - cross-shore

              'yFRF': FRF local coordinate system - alongshore

        """
        loc_dict = {}
        for g in self.waveGaugeList:
            loc_dict[g] = {}
            data = loc_dict[g]
            
            # Get latlon from Thredds server.
            try:
                if g in ['11', '12', '13', '14', '21', '22', '23', '24']:
                    latlon = self.getBathyDuckLoc(gaugenumber=g)
                else:
                    latlon = self.getWaveGaugeLoc(g)
            except IOError:
                continue
            # lat and lon values currently stored as 1 element arrays.
            # Cast to float for consistency.
            lat = float(latlon['Lat'])
            lon = float(latlon['Lon'])
            
            # Covert latlon to stateplane.
            coords = gp.LatLon2ncsp(lon, lat)
            spE = coords['StateplaneE']
            spN = coords['StateplaneN']
            
            # Convert stateplane to frf coords.
            frfcoords = gp.ncsp2FRF(spE, spN)
            xfrf = frfcoords['xFRF']
            yfrf = frfcoords['yFRF']
            
            data['lat'] = lat
            data['lon'] = lon
            data['spE'] = spE
            data['spN'] = spN
            data['xFRF'] = xfrf
            data['yFRF'] = yfrf
        
        return loc_dict
    
    def get_sensor_locations(self, datafile='frf_sensor_locations.pkl', window_days=14):
        """Retrieve sensor coordinate dictionary from file if there is an entry.
        
        Look within window_days of the specified timestampstr. Otherwise query the
        Thredds server for location information and update archived data
        accordingly.

        Args:
          datafile (str): Name of file containing archived sensor location data.
               Updates datafile when new information is retrieved.
               (Default value = 'frf_sensor_locations.pkl')
          window_days (int): Maximum interval between desired timestamp and closest timestamp
               in datafile to use archived data. If this interval is larger than
               window_days days, query the Thredds server. (Default value = 14)

        Returns:
          sensor_locations (dict):  Coordinates in lat/lon, stateplane, and frf for each available
               gaugenumber (gauges 0 to 12).
          
        Notes:
            Updates datafile when new information is retrieved.

        """
        try:
            with open(datafile, 'rb') as fid:  # this will close a file when done loading it
                loc_dict = pickle.load(fid)
        except (IOError, UnicodeDecodeError):
            # now create pickle if one's not around
            loc_dict = {}
            # this should be date of searching( or date of gauge placement more accurately)
            timestamp = self.d1  # DT.strptime(timestampstr, '%Y%m%d_%H%M')
            loc_dict[timestamp] = self.get_sensor_locations_from_thredds()  # timestamp)
            with open(datafile, 'wb') as fid:
                pickle.dump(loc_dict, fid)
        
        available_timestamps = np.array(list(loc_dict.keys()))
        idx = np.abs(self.d1 - available_timestamps).argmin()
        nearest_timestamp = available_timestamps[idx]
        if abs(self.d1 - nearest_timestamp).days < window_days:
            # if there is data, and its within the window
            archived_sensor_locations = loc_dict[nearest_timestamp]
            # MPG: only use locations specified in self.waveGaugeList (for the case
            # that there are archived locations that should not be used).
            sensor_locations = collections.OrderedDict()
            
            for g in self.waveGaugeList:
                if g in archived_sensor_locations:
                    sensor_locations[g] = archived_sensor_locations[g]
                elif g not in archived_sensor_locations:
                    sensor_locations = collections.OrderedDict()
                    sensor_locations = self.get_sensor_locations_from_thredds()
                    return sensor_locations
                else:
                    # MPG: use empty dict as a placeholder to indicate that no
                    # data is available.
                    sensor_locations[g] = {}
        else:
            sensor_locations = self.get_sensor_locations_from_thredds()
            loc_dict[self.d1] = sensor_locations
            with open(datafile, 'wb') as fid:
                pickle.dump(loc_dict, fid)
        
        return sensor_locations
    
    def getLidarRunup(self, removeMasked=True):
        """This function will get the wave runup measurements from the lidar mounted in the dune.

        Args:
          removeMasked: if data come back as masked, remove from the arrays removeMasked will
            toggle the removing of data points from the tsTime series based on the flag
            status (Default value = True)

        Returns:
          dictionary with collected data.  keys listed below (for more info see the netCDF file
          metadata)
            'name': gauge name

            'lat': latitude for points

            'lon': longitude for points

            'lidarX': the x coordinate of the lidar (making the measurement)

            'lidarY': the y coordinate of the lidar (making the measurement)

            'time': time stamp for data in datetime object

            'totalWaterLevel': total elevation of water dimensioned by time

            'elevation': 1D array of swash elevation at times listed in tsTime.

            'samplingTime': time stamp for time series data

            'xFRF': x location of the data points

            'yFRF': y location of the data points

            'totalWaterLevelQCflag': qc flag following quartod standards for total water level

            'percentMissing': percent of missing data can be used as a confidence factor for
            measurement

        """
        self.dataloc = 'oceanography/waves/lidarWaveRunup/lidarWaveRunup.ncml'
        self.ncfile, self.allEpoch = getnc(dataLoc=self.dataloc, callingClass=self.callingClass,
                                           dtRound=1 * 60)
        self.lidarIndex = gettime(allEpoch=self.allEpoch, epochStart=self.epochd1,
                                  epochEnd=self.epochd2)
        
        if np.size(self.lidarIndex) > 0 and self.lidarIndex is not None:
            out = {'name':                  nc.chartostring(self.ncfile['station_name'][:]),
                   'lat':                   self.ncfile['lidarLatitude'][:],  # Coordintes
                   'lon':                   self.ncfile['lidarLongitude'][:],
                   'lidarX':                self.ncfile['lidarX'][:],
                   'lidarY':                self.ncfile['lidarY'][:],
                   'time':                  nc.num2date(self.allEpoch[self.lidarIndex], self.ncfile['time'].units,
                                                        self.ncfile['time'].calendar, only_use_cftime_datetimes=False),
                   'epochtime':             self.allEpoch[self.lidarIndex],
                   'totalWaterLevel':       self.ncfile['totalWaterLevel'][self.lidarIndex],
                   'elevation':             self.ncfile['elevation'][self.lidarIndex, :],
                   'xFRF':                  self.ncfile['xFRF'][self.lidarIndex, :],
                   'yFRF':                  self.ncfile['yFRF'][self.lidarIndex, :],
                   'samplingTime':          self.ncfile['tsTime'][:],
                   'totalWaterLevelQCflag': self.ncfile['totalWaterLevelQCFlag'][self.lidarIndex],
                   'percentMissing':        self.ncfile['percentTimeSeriesMissing'][self.lidarIndex],
                   }
            
            if removeMasked:
                if isinstance(out['elevation'], np.ma.MaskedArray):
                    # out['elevation'] = np.array(out['elevation'][~out['elevation'].mask])
                    out['elevation'] = out['elevation'].filled(np.nan)
                else:
                    pass
                if isinstance(out['totalWaterLevel'], np.ma.MaskedArray):
                    out['totalWaterLevel'] = np.array(
                        out['totalWaterLevel'][~out['totalWaterLevel'].mask])
                else:
                    pass
                if isinstance(out['xFRF'], np.ma.MaskedArray):
                    out['xFRF'] = np.array(out['xFRF'][~out['xFRF'].mask])
                else:
                    pass
                if isinstance(out['yFRF'], np.ma.MaskedArray):
                    out['yFRF'] = np.array(out['yFRF'][~out['yFRF'].mask])
                else:
                    pass
                # if isinstance(out['runupDownLine'], np.ma.MaskedArray):
                #     out['runupDownLine'] = np.array(out['runupDownLine'][~out[
                #     'runupDownLine'].mask])
                # else:
                #     pass
                if isinstance(out['samplingTime'], np.ma.MaskedArray):
                    out['samplingTime'] = np.array(out['samplingTime'][~out['samplingTime'].mask])
                else:
                    pass
            else:
                pass
        
        else:
            print('There is no LIDAR data during this time period')
            out = None
        return out
    
    def getCTD(self):
        """THIS FUNCTION IS CURRENTLY BROKEN.
        
        THE PROBLEM IS THAT self.cshore_ncfile does not have any keys?
        TODO fix this function
        This function gets the CTD data from the thredds server
        
        Args:  None

        Returns:
            dict: output dictionary with keys listed below, will return None if error happens
                'depth': depth of sample

                'temp': water temperature

                'time' (obj): datetime object time stamp

                'lat': latitude

                'lon': longitude

                'salin': salinity

                'soundSpeed': speed of sound

                'sigmaT':

        """
        # do check here on profile numbers
        # acceptableProfileNumbers = [None, ]
        self.dataloc = 'oceanography/ctd/eop-ctd/eop-ctd.ncml'  # location of the gridded surveys
        
        try:
            self.ncfile = self.FRFdataloc + self.dataloc
            val = (max([n for n in (self.ncfile['time'][:] - self.epochd1) if n < 0]))
            idx = np.where((self.ncfile['time'][:] - self.epochd1) == val)[0][0]
            print('CTD data is closest in HISTORY - operational')
        
        except (RuntimeError, NameError, AssertionError,
                TypeError):  # if theres any error try to get good data from next location
            try:
                self.ncfile = self.chlDataLoc + self.dataloc
                val = (max([n for n in (self.ncfile['time'][:] - self.epochd1) if n < 0]))
                idx = np.where((self.ncfile['time'][:] - self.epochd1) == val)[0][0]
                print('CTD data closest in HISTORY - operational')
            except (RuntimeError, NameError, AssertionError,
                    TypeError):  # if there are still errors, give up
                idx = []
        
        if np.size(idx) > 0:
            # now retrieve data with idx
            depth = self.ncfile['depth'][idx]
            lat = self.ncfile['lat'][idx]
            lon = self.ncfile['lon'][idx]
            time = nc.num2date(self.ncfile['time'][idx], self.ncfile['time'].units, only_use_cftime_datetimes=False)
            temp = self.ncfile['waterTemperature'][idx]
            salin = self.ncfile['salinity'][idx]
            soundSpeed = self.ncfile['soundSpeed'][idx]
            sigmaT = self.ncfile['sigmaT'][idx]
            
            ctd_Dict = {'depth':      depth,
                        'temp':       temp,
                        'time':       time,
                        'lat':        lat,
                        'lon':        lon,
                        'salin':      salin,
                        'soundSpeed': soundSpeed,
                        'sigmaT':     sigmaT}
        else:
            ctd_Dict = None
        
        return ctd_Dict
    
    def getALT(self, gaugeName=None, removeMasked=True):
        """This function gets the Altimeter data from the thredds server.

        Args:
          gaugeName (str):  This is just the name of the altimeter we want to use
             available gauge names listed below (Default value = None)
             'Alt03', 'Alt04', 'Alt05', 'Alt769-150', 'Alt769-200', 'Alt769-250','Alt769-300',
             'Alt769-350', 'Alt861-150', 'Alt861-200', 'Alt861-250', 'Alt861-300', 'Alt861-350'

          removeMasked (bool): remove the data that are masked (Default value = True)

        Returns:
          a dictionary with below keys with selected data, for more info see netCDF files on server
            'name': file title

            'time': date time objects

            'epochtime': time in epoch, seconds since 1970

            'lat': latitude of location of data

            'PKF': Kalman Filtered Error covariance estimate

            'lon': longitude location of data

            'xFRF': x location of data

            'yFRF': y location of data

            'stationName: station name variable

            'timeStart': start time of the sample

            'timeEnd': end time of the sample

            'bottomElev': Kalman filtered elevation

        """
        # location of the data
        gauge_list = ['Alt03', 'Alt04', 'Alt05',
                      'Alt769-150', 'Alt769-200', 'Alt769-250', 'Alt769-300', 'Alt769-350',
                      'Alt861-150', 'Alt861-200', 'Alt861-250', 'Alt861-300', 'Alt861-350']
        if gaugeName not in gauge_list:
            raise NotImplementedError('Input string is not a valid gage name, please check')
        if gaugeName in ['Alt05']:
            self.dataloc = 'geomorphology/altimeter/Alt05-altimeter/Alt05-altimeter.ncml'
        elif gaugeName in ['Alt04']:
            self.dataloc = 'geomorphology/altimeter/Alt04-altimeter/Alt04-altimeter.ncml'
        elif gaugeName in ['Alt03']:
            self.dataloc = 'geomorphology/altimeter/Alt03-altimeter/Alt03-altimeter.ncml'
        elif gaugeName in ['Alt769-150']:
            self.dataloc = 'geomorphology/altimeter/Alt769-150-altimeter/Alt769-150-altimeter.ncml'
        elif gaugeName in ['Alt769-200']:
            self.dataloc = 'geomorphology/altimeter/Alt769-200-altimeter/Alt769-200-altimeter.ncml'
        elif gaugeName in ['Alt769-250']:
            self.dataloc = 'geomorphology/altimeter/Alt769-250-altimeter/Alt769-250-altimeter.ncml'
        elif gaugeName in ['Alt769-300']:
            self.dataloc = 'geomorphology/altimeter/Alt769-300-altimeter/Alt769-300-altimeter.ncml'
        elif gaugeName in ['Alt769-350']:
            self.dataloc = 'geomorphology/altimeter/Alt769-350-altimeter/Alt769-350-altimeter.ncml'
        elif gaugeName in ['Alt861-150']:
            self.dataloc = 'geomorphology/altimeter/Alt861-150-altimeter/Alt861-150-altimeter.ncml'
        elif gaugeName in ['Alt861-200']:
            self.dataloc = 'geomorphology/altimeter/Alt861-200-altimeter/Alt861-200-altimeter.ncml'
        elif gaugeName in ['Alt861-250']:
            self.dataloc = 'geomorphology/altimeter/Alt769-250-altimeter/Alt769-250-altimeter.ncml'
        elif gaugeName in ['Alt861-300']:
            self.dataloc = 'geomorphology/altimeter/Alt769-300-altimeter/Alt769-300-altimeter.ncml'
        elif gaugeName in ['Alt861-350']:
            self.dataloc = 'geomorphology/altimeter/Alt769-350-altimeter/Alt769-350-altimeter.ncml'
        else:
            raise NotImplementedError('Please use one of the following keys\n'.format(gauge_list))
        
        self.ncfile, self.allEpoch = getnc(dataLoc=self.dataloc, callingClass=self.callingClass,
                                           dtRound=1 * 60)
        altdataindex = gettime(allEpoch=self.allEpoch, epochStart=self.epochd1,
                               epochEnd=self.epochd2)
        
        # get the actual current data
        if np.size(altdataindex) > 1:
            alt_lat = self.ncfile['Latitude'][0]  # pulling latitude
            alt_lon = self.ncfile['Longitude'][0]  # pulling longitude
            alt_be = self.ncfile['bottomElevation'][altdataindex]  # pulling bottom elevation
            alt_pkf = self.ncfile['PKF'][altdataindex]  # ...
            alt_stationname = nc.chartostring(self.ncfile['station_name'][:])  # name of the station
            self.alt_timestart = nc.num2date(self.ncfile['timestart'][altdataindex],
                                             self.ncfile['timestart'].units,
                                             self.ncfile['time'].calendar, only_use_cftime_datetimes=False)
            self.alt_timeend = nc.num2date(self.ncfile['timeend'][altdataindex],
                                           self.ncfile['timeend'].units,
                                           self.ncfile['time'].calendar, only_use_cftime_datetimes=False)
            self.alt_time = nc.num2date(self.ncfile['time'][altdataindex],
                                        self.ncfile['time'].units,
                                        self.ncfile['time'].calendar, only_use_cftime_datetimes=False)
            for num in range(0, len(self.alt_time)):
                self.alt_time[num] = self._roundtime(self.alt_time[num], roundto=1 * 60)
                self.alt_timestart[num] = self._roundtime(self.alt_timestart[num], roundto=1 * 60)
                self.alt_timeend[num] = self._roundtime(self.alt_timeend[num], roundto=1 * 60)
            
            alt_coords = gp.FRFcoord(alt_lon, alt_lat)
            
            if removeMasked:
                altpacket = {'name':        str(self.ncfile.title),
                             'time':        np.array(self.alt_time[~alt_be.mask]),
                             'epochtime':   np.array(self.allEpoch[altdataindex][~alt_be.mask]),
                             'lat':         alt_lat,
                             'PKF':         np.array(alt_pkf[~alt_be.mask]),
                             'lon':         alt_lon,
                             'xFRF':        alt_coords['xFRF'],
                             'yFRF':        alt_coords['yFRF'],
                             'stationName': alt_stationname,
                             'timeStart':   np.array(self.alt_timestart[~alt_be.mask]),
                             'timeEnd':     np.array(self.alt_timeend[~alt_be.mask]),
                             'bottomElev':  np.array(alt_be[~alt_be.mask])}
            else:
                altpacket = {'name':        str(self.ncfile.title),
                             'time':        self.alt_time,
                             'lat':         alt_lat,
                             'PKF':         alt_pkf,
                             'lon':         alt_lon,
                             'xFRF':        alt_coords['xFRF'],
                             'yFRF':        alt_coords['yFRF'],
                             'stationName': alt_stationname,
                             'timeStart':   self.alt_timestart,
                             'timeEnd':     self.alt_timeend,
                             'bottomElev':  alt_be}
            
            return altpacket
        else:
            print('No %s data found for this period' % (gaugeName))
            self.altpacket = None
            return self.altpacket
    
    def getLidarWaveProf(self, removeMasked=True):
        """Grabs wave profile data from Lidar gauge.

        Args:
          removeMasked: Removes masked values from data. Default value = (True)

        Returns:
          dictionary with keys below, will None if there's no data or an error
            'name': station name

            'lat': latitude

            'lon': longitude

            'lidarX': station frf X coordinate

            'lidarY': station frf y coordinate

            'xFRF' (array): cross shore coordinate

            'yFRF' (array): along-shore coordinate

            'runupDownLine':  wave runDown (need more details here)

            'waveFrequency' (array): wave frequencies for specturm

            'time' (obj): date time object time stamp

            'hydroQCflag': QC flag for hydro

            'waterLevel': water level [time, xFRF]

            'waveHs': wave height in the incident (gravity) frequency band

            'waveHsIG': wave height inf the infragravity band

            'waveHsTotal': total frequency band wave height

            'waveSkewness':  need more detail here

            'waveAsymmetry':  Need more detail here

            'waveEnergyDensity' (array): frequency spectra [time, wave frequency]

            'percentMissing': quality check - tells user how much of the time series was missing
            while these stats
            were calculated

        """
        self.dataloc = 'oceanography/waves/lidarHydrodynamics/lidarHydrodynamics.ncml'
        self.ncfile, self.allEpoch = getnc(dataLoc=self.dataloc, callingClass=self.callingClass,
                                           dtRound=1 * 60)
        self.lidarIndex = gettime(allEpoch=self.allEpoch, epochStart=self.epochd1,
                                  epochEnd=self.epochd2)
        if np.size(self.lidarIndex) > 0 and self.lidarIndex is not None:
            out = {'name':              nc.chartostring(self.ncfile['station_name'][:]),
                   'time':              nc.num2date(self.ncfile['time'][self.lidarIndex],
                                                    self.ncfile['time'].units,
                                                    self.ncfile['time'].calendar, only_use_cftime_datetimes=False),
                   'lat':               self.ncfile['lidarLatitude'][:],  # Coordinates
                   'lon':               self.ncfile['lidarLongitude'][:],
                   'lidarX':            self.ncfile['lidarX'][:],
                   'lidarY':            self.ncfile['lidarY'][:],
                   'xFRF':              self.ncfile['xFRF'][:],
                   'yFRF':              self.ncfile['yFRF'][:],
                   'waveFrequency':     self.ncfile['waveFrequency'][:],
                   'hydroQCflag':       self.ncfile['hydrodynamicsFlag'][self.lidarIndex],
                   'waterLevel':        self.ncfile['waterLevel'][self.lidarIndex, :],
                   'waveHs':            self.ncfile['waveHs'][self.lidarIndex, :],
                   'waveHsIG':          self.ncfile['waveHsIG'][self.lidarIndex, :],
                   'waveHsTotal':       self.ncfile['waveHsTotal'][self.lidarIndex, :],
                   'waveSkewness':      self.ncfile['waveSkewness'][self.lidarIndex, :],
                   'waveAsymmetry':     self.ncfile['waveAsymmetry'][self.lidarIndex, :],
                   'waveEnergyDensity': self.ncfile['waveEnergyDensity'][self.lidarIndex, :, :],
                   'percentMissing':    self.ncfile['percentTimeSeriesMissing'][self.lidarIndex, :],
                   }
            
            if removeMasked:
                # TODO lets put this into a loop
                if isinstance(out['waterLevel'], np.ma.MaskedArray):
                    out['waterLevel'] = np.array(out['waterLevel'][~out['waterLevel'].mask])
                
                if isinstance(out['waveHs'], np.ma.MaskedArray):
                    out['waveHs'] = np.array(out['waveHs'][~out['waveHs'].mask])
                
                if isinstance(out['waveHsIG'], np.ma.MaskedArray):
                    out['waveHsIG'] = np.array(out['waveHsIG'][~out['waveHsIG'].mask])
                
                if isinstance(out['waveHsTotal'], np.ma.MaskedArray):
                    # DLY note 01092019 - this bit of codes turns the 2d array out['waveHsTotal']
                    # of time by distance
                    # into a 1-d array?
                    # i dont think we can have this be an option for 2d data?
                    out['waveHsTotal'] = np.array(out['waveHsTotal'][~out['waveHsTotal'].mask])
                
                if isinstance(out['waveSkewness'], np.ma.MaskedArray):
                    out['waveSkewness'] = np.array(out['waveSkewness'][~out['waveSkewness'].mask])
                
                if isinstance(out['waveAsymmetry'], np.ma.MaskedArray):
                    out['waveAsymmetry'] = np.array(
                        out['waveAsymmetry'][~out['waveAsymmetry'].mask])
                
                if isinstance(out['waveEnergyDensity'], np.ma.MaskedArray):
                    out['waveEnergyDensity'] = np.array(
                        out['waveEnergyDensity'][~out['waveEnergyDensity'].mask])
        
        else:
            print('There is no LIDAR data during this time period')
            out = None
        return out
    
    def getLidarTopo(self, **kwargs):
        """This function will get the lidar DEM data, beach topography data.

        This function is not finished being developed

        Args: None

        Keyword Args:
           'xbounds': frf cross-shore bounds
           'ybounds': frf alongshore bounds
            lidarLoc':
        Returns:
          dictionary with lidar beach topography
            keys to be determined

        """
        lidarLoc = kwargs.get('lidarLoc', 'dune')
        if lidarLoc == 'pier':
            self.dataloc = 'geomorphology/DEMs/pierLidarDEM/pierLidarDEM.ncml'
        elif lidarLoc == 'dune':
            self.dataloc = 'geomorphology/DEMs/duneLidarDEM/duneLidarDEM.ncml'
        elif lidarLoc == 'claris':
            self.dataloc = 'geomorphology/DEMs/clarisDEM/clarisDEM.ncml'
        else:
            raise NotImplementedError('valid lidar locs are "dune" and "pier"')
        
        self.ncfile, self.allEpoch = getnc(dataLoc=self.dataloc, callingClass=self.callingClass,
                                           dtRound=1 * 60)
        self.idxDEM = gettime(allEpoch=self.allEpoch, epochStart=self.epochd1,
                              epochEnd=self.epochd2)
        self.DEMtime = nc.num2date(self.allEpoch[self.idxDEM], 'seconds since 1970-01-01',
                                   only_use_cftime_datetimes=False)
        
        if 'xbounds' in kwargs and np.array(kwargs['xbounds']).size == 2:
            if kwargs['xbounds'][0] > kwargs['xbounds'][1]:
                kwargs['xbounds'] = np.flip(kwargs['xbounds'], axis=0)
            # first min of x
            if (kwargs['xbounds'][0] < self.ncfile['xm'][:]).all():
                # then set xmin to 0
                removeMinX = 0
            else:  # <= used here to handle inclusive initial index inherant in python
                removeMinX = np.argwhere(
                    self.ncfile['xFRF'][:] <= kwargs['xbounds'][0]).squeeze().max()
            # now max of x
            if (kwargs['xbounds'][1] > self.ncfile['xFRF'][:]).all():
                removeMaxX = None
            else:
                removeMaxX = np.argwhere(
                    self.ncfile['xFRF'][:] >= kwargs['xbounds'][
                        1]).squeeze().min() + 1  # python indexing
            xs = slice(removeMinX, removeMaxX)
        else:
            xs = slice(None)
        
        if 'ybounds' in kwargs and np.array(kwargs['ybounds']).size == 2:
            if kwargs['ybounds'][0] > kwargs['ybounds'][1]:
                kwargs['ybounds'] = np.flip(kwargs['ybounds'], axis=0)
            # first min of y
            if (kwargs['ybounds'][0] < self.ncfile['yFRF'][:]).all():
                # then set the ymin to first index [0]
                removeMinY = 0  # ie get all data
            else:
                removeMinY = np.argwhere(
                    self.ncfile['yFRF'][:] <= kwargs['ybounds'][0]).squeeze().max()
            ## now max of y
            if (kwargs['ybounds'][1] > self.ncfile['yFRF'][:]).all():
                removeMaxY = None
            else:
                removeMaxY = np.argwhere(
                    self.ncfile['yFRF'][:] >= kwargs['ybounds'][
                        1]).squeeze().min() + 1  # python indexing
            ys = slice(removeMinY, removeMaxY)
        else:
            ys = slice(None)
        
        DEMdata = {
                'xFRF':      self.ncfile['xFRF'][self.idxDEM, xs],
                'yFRF':      self.ncfile['yFRF'][self.idxDEM, ys],
                'elevation': self.ncfile['elevation'][self.idxDEM, ys, xs],
                'time':      self.DEMtime,
                'epochtime': self.allEpoch[self.idxDEM],
        
                'lat':       self.ncfile['yFRF'][self.idxDEM, ys],
                'lon':       lon
                }
        return DEMdata
    
    def getBathyRegionalDEM(self, utmEmin, utmEmax, utmNmin, utmNmax):
        """Grabs bathymery from the regional background grid.

        Args:
          utmEmin: left side of DEM bounding box in UTM
          utmEmax: right side of DEM bounding box in UTM
          utmNmin: bottom of DEM bounding box in UTM
          utmNmax: top of DEM bounding box in UTM

        Returns:
          dictionary comprising a smaller rectangular piece of the DEM data, bounded by inputs above
          'utmEasting': UTM Easting [meters]

          'utmNorthing': UTM  Northing [meters]

          'latitude': self explanatory

          'longitude': self explanatory

          'bottomElevation': elevation [m] (NAVD88)


        """
        self.dataloc = 'integratedBathyProduct/RegionalBackgroundDEM/backgroundDEM.nc'
        self.ncfile = nc.Dataset(self.crunchDataLoc + self.dataloc)
        
        # get a 1D ARRAY of the utmE and utmN of the rectangular grid (NOT the full grid!!!)
        utmE_all = self.ncfile['utmEasting'][0, :]
        utmN_all = self.ncfile['utmNorthing'][:, 0]
        
        # find indices I need to pull...
        ni_min = np.where(utmE_all >= utmEmin)[0][0]
        ni_max = np.where(utmE_all <= utmEmax)[0][-1]
        nj_min = np.where(utmN_all <= utmNmax)[0][0]
        nj_max = np.where(utmN_all >= utmNmin)[0][-1]
        
        assert (np.size(ni_min) >= 1) & (np.size(ni_max) >= 1) & (np.size(nj_min) >= 1) & (
            np.size(
                nj_max) >= 1), 'getBathyDEM Error: bounding box is too close to edge of DEM domain'
        
        out = {}
        out['utmEasting'] = self.ncfile['utmEasting'][nj_min:nj_max + 1, ni_min:ni_max + 1]
        out['utmNorthing'] = self.ncfile['utmNorthing'][nj_min:nj_max + 1, ni_min:ni_max + 1]
        out['latitude'] = self.ncfile['latitude'][nj_min:nj_max + 1, ni_min:ni_max + 1]
        out['longitude'] = self.ncfile['longitude'][nj_min:nj_max + 1, ni_min:ni_max + 1]
        out['bottomElevation'] = self.ncfile['bottomElevation'][nj_min:nj_max + 1,
                                 ni_min:ni_max + 1]
        
        return out
    
    def getBathyGridcBathy(self, **kwargs):
        """This function gets the cbathy data from the below address, assumes fill value of -999.

        Keyword Args:
            xbound: = [xmin, xmax]  which will truncate the cbathy domain to xmin, xmax (frf coord)

            ybound: = [ymin, ymax]  which will truncate the cbathy domain to ymin, ymax (frf coord)

        Returns:
            dictionary with keys below, will return None if error is found
                'time': time

                'xm':  frf xoordinate x's

                'ym': frf ycoordinates

                'depth': raw cbathy depths

                'depthfC: same as depth

                'depthKF':  kalman filtered hourly depth

                'depthKFError': errors associated with the kalman filter

                'fB':  ?

                'k':  ??

        """
        fillValue = -999  # assumed fill value from the rest of the files taken as less than or
        # equal to
        self.dataloc = 'projects/bathyduck/data/cbathy_old/cbathy.ncml'
        self.ncfile, self.allEpoch = getnc(dataLoc=self.dataloc, callingClass=self.callingClass,
                                           dtRound=30 * 60)
        self.cbidx = gettime(allEpoch=self.allEpoch, epochStart=self.epochd1, epochEnd=self.epochd2)
        
        self.cbtime = nc.num2date(self.allEpoch[self.cbidx], 'seconds since 1970-01-01',
                                  only_use_cftime_datetimes=False)
        # mask = (time > start) & (time < end)
        # assert (emask == mask).all(), 'epoch time is not working'
        # idx = np.where(emask)[0] # this leaves a list that keeps the data iteratable with a
        # size 1.... DON'T CHANGE
        if np.size(self.cbidx) == 1 and self.cbidx == None:
            cbdata = None  # throw a kick out if there's no data avaiable
            return cbdata
        # truncating data from experimental parameters to
        if 'xbounds' in kwargs and np.array(kwargs['xbounds']).size == 2:
            if kwargs['xbounds'][0] > kwargs['xbounds'][1]:
                kwargs['xbounds'] = np.flip(kwargs['xbounds'], axis=0)
            # first min of x
            if (kwargs['xbounds'][0] < self.ncfile['xm'][:]).all():
                # then set xmin to 0
                removeMinX = 0
            else:  # <= used here to handle inclusive initial index inherant in python
                removeMinX = np.argwhere(
                    self.ncfile['xm'][:] <= kwargs['xbounds'][0]).squeeze().max()
            # now max of x
            if (kwargs['xbounds'][1] > self.ncfile['xm'][:]).all():
                removeMaxX = None
            else:
                removeMaxX = np.argwhere(
                    self.ncfile['xm'][:] >= kwargs['xbounds'][
                        1]).squeeze().min() + 1  # python indexing
            xs = slice(removeMinX, removeMaxX)
        else:
            xs = slice(None)
        
        if 'ybounds' in kwargs and np.array(kwargs['ybounds']).size == 2:
            if kwargs['ybounds'][0] > kwargs['ybounds'][1]:
                kwargs['ybounds'] = np.flip(kwargs['ybounds'], axis=0)
            # first min of y
            if (kwargs['ybounds'][0] < self.ncfile['ym'][:]).all():
                # then set the ymin to first index [0]
                removeMinY = 0  # ie get all data
            else:
                removeMinY = np.argwhere(
                    self.ncfile['ym'][:] <= kwargs['ybounds'][0]).squeeze().max()
            ## now max of y
            if (kwargs['ybounds'][1] > self.ncfile['ym'][:]).all():
                removeMaxY = None
            else:
                removeMaxY = np.argwhere(
                    self.ncfile['ym'][:] >= kwargs['ybounds'][
                        1]).squeeze().min() + 1  # python indexing
            ys = slice(removeMinY, removeMaxY)
        else:
            ys = slice(None)
        
        try:
            cbdata = {'time':         self.cbtime,  # round the time to the nearest 30 minutes
                      'epochtime':    self.allEpoch[self.cbidx],
                      'xm':           self.ncfile['xm'][xs],
                      'ym':           self.ncfile['ym'][ys],
                      'depthKF':      np.ma.array(self.ncfile['depthKF'][self.cbidx, ys, xs],
                                                  mask=(self.ncfile['depthKF'][
                                                            self.cbidx, ys, xs] <= fillValue),
                                                  fill_value=np.nan),
                      'depthKFError': np.ma.array(self.ncfile['depthKF'][self.cbidx, ys, xs],
                                                  mask=(self.ncfile['depthKF'][
                                                            self.cbidx, ys, xs] <= fillValue),
                                                  fill_value=np.nan),
                      'depthfC':      np.ma.array(self.ncfile['depthfC'][self.cbidx, ys, xs],
                                                  mask=(self.ncfile['depthfC'][
                                                            self.cbidx, ys, xs] <= fillValue),
                                                  fill_value=np.nan),
                      'depthfCError': np.ma.array(self.ncfile['depthErrorfC'][self.cbidx, ys, xs],
                                                  mask=(self.ncfile['depthErrorfC'][
                                                            self.cbidx, ys, xs] <= fillValue),
                                                  fill_value=np.nan),
                      'fB':           np.ma.array(self.ncfile['fB'][self.cbidx, ys, xs, :],
                                                  mask=(self.ncfile['fB'][self.cbidx, ys, xs,
                                                        :] <= fillValue),
                                                  fill_value=np.nan),
                      'k':            np.ma.array(self.ncfile['k'][self.cbidx, ys, xs, :],
                                                  mask=(self.ncfile['k'][self.cbidx, ys, xs, :] <= fillValue),
                                                  fill_value=np.nan),
                      'P':            np.ma.array(self.ncfile['PKF'][self.cbidx, ys, xs],
                                                  mask=(self.ncfile['PKF'][self.cbidx, ys, xs] <= fillValue),
                                                  fill_value=np.nan)}  # may need to be masked
            
            assert ~cbdata[
                'depthKF'].mask.all(), 'all Cbathy kalman filtered data retrieved are masked '
            print('Grabbed cBathy Data, successfully')
        
        except (IndexError, AssertionError):  # there's no data in the Cbathy
            cbdata = None
        
        return cbdata
    
    def getArgus(self, type, **kwargs):
        """Grabs argus data from the bathyDuck time period, particularly staple products.
        
        Currently this is only retrieves variance and timex images.

        Args:
            type (str): this is a string that describes the video product eg var, timex

        Keyword Args:
            xbound: = [xmin, xmax]  which will truncate the cbathy domain to xmin, xmax (frf coord)

            ybound: = [ymin, ymax]  which will truncate the cbathy domain to ymin, ymax (frf coord)

        Returns:
            None

        """
        from skimage import color
        if type.lower() not in ['var', 'timex']:
            raise NotImplementedError(
                "These data are not currently available through this function")
        elif type.lower() in ['var', 'variance']:
            self.dataloc = "projects/bathyduck/data/argus/variance/variance.ncml"
        elif type.lower() in ['timex']:
            self.dataloc = "projects/bathyduck/data/argus/timex/timex.ncml"
        
        ################ go get data index
        self.ncfile, self.allEpoch = getnc(dataLoc=self.dataloc, callingClass=self.callingClass,
                                           dtRound=1 * 60)
        self.idxArgus = gettime(allEpoch=self.allEpoch, epochStart=self.epochd1,
                                epochEnd=self.epochd2)
        
        ###### sub divide bounds in kwargs
        if 'xbounds' in kwargs and np.array(kwargs['xbounds']).size == 2:
            if kwargs['xbounds'][0] > kwargs['xbounds'][1]:
                kwargs['xbounds'] = np.flip(kwargs['xbounds'], axis=0)
            # first min of x
            if (kwargs['xbounds'][0] < self.ncfile['x'][:]).all():
                # then set xmin to 0
                removeMinX = 0
            else:  # <= used here to handle inclusive initial index inherant in python
                removeMinX = np.argwhere(
                    self.ncfile['x'][:] <= kwargs['xbounds'][0]).squeeze().max()
            # now max of x
            if (kwargs['xbounds'][1] > self.ncfile['x'][:]).all():
                removeMaxX = None
            else:
                removeMaxX = np.argwhere(
                    self.ncfile['x'][:] >= kwargs['xbounds'][
                        1]).squeeze().min() + 1  # python indexing
            xs = slice(removeMinX, removeMaxX)
        else:
            xs = slice(None)
        
        if 'ybounds' in kwargs and np.array(kwargs['ybounds']).size == 2:
            if kwargs['ybounds'][0] > kwargs['ybounds'][1]:
                kwargs['ybounds'] = np.flip(kwargs['ybounds'], axis=0)
            # first min of y
            if (kwargs['ybounds'][0] < self.ncfile['y'][:]).all():
                # then set the ymin to first index [0]
                removeMinY = 0  # ie get all data
            else:
                removeMinY = np.argwhere(
                    self.ncfile['y'][:] <= kwargs['ybounds'][0]).squeeze().max()
            ## now max of y
            if (kwargs['ybounds'][1] > self.ncfile['y'][:]).all():
                removeMaxY = None
            else:
                removeMaxY = np.argwhere(
                    self.ncfile['y'][:] >= kwargs['ybounds'][
                        1]).squeeze().min() + 1  # python indexing
            ys = slice(removeMinY, removeMaxY)
        else:
            ys = slice(None)
        try:
            timeArgus = nc.num2date(self.allEpoch[self.idxArgus], 'seconds since 1970-01-01',
                                    only_use_cftime_datetimes=False)
            Ip = self.ncfile['Ip'][self.idxArgus, xs, ys]
            out = {'time':      timeArgus,
                   'epochtime': self.allEpoch[self.idxArgus],
                   'rgb':       Ip,
                   'bw':        color.rgb2gray(Ip),
                   'xFRF':      self.ncfile['x'][xs],
                   'yFRF':      self.ncfile['y'][ys], }
        
        except(IndexError, AssertionError):
            out = None
        return out


class getDataTestBed:
    """Retrieves model data."""
    def __init__(self, d1, d2):
        """Initialization description here.
        
        Data are returned in self.datainex are inclusive at d1,d2
        """
        self.rawdataloc_wave = []
        self.outputdir = []  # location for outputfiles
        self.start = d1  # start date for data grab
        self.end = d2  # end data for data grab
        self.timeunits = 'seconds since 1970-01-01 00:00:00'
        self.epochd1 = nc.date2num(self.start, self.timeunits)
        self.epochd2 = nc.date2num(self.end, self.timeunits)
        self.comp_time()
        # if THREDDS is None:
        #     ipAddress = socket.gethostbyname(socket.gethostname())
        #     if ipAddress.startswith('134.164.129'):  # FRF subdomain
        #         self.THREDDS = 'FRF'
        #     else:
        #         self.THREDDS = 'CHL'
        #
        self.callingClass = 'getDataTestBed'
        self.FRFdataloc = u'http://134.164.129.55/thredds/dodsC/FRF/'
        self.crunchDataLoc = u'http://134.164.129.55/thredds/dodsC/cmtb/'
        self.chlDataLoc = u'https://chlthredds.erdc.dren.mil/thredds/dodsC/frf/'
        
        assert type(self.end) == DT.datetime, 'end dates need to be in python "Datetime" data types'
        assert type(
            self.start) == DT.datetime, 'start dates need to be in python "Datetime" data types'
    
    def comp_time(self):
        """Test if times are backwards."""
        assert self.end >= self.start, 'finish time: end needs to be after start time: start'
    
    def gettime(self, dtRound=60):
        """This function opens the netcdf file, pulls down all of the time, then pulls the dates of interest.
        
        from the THREDDS (data loc) server based on start,end, and data location
        it returns the indicies in the NCML file of the dates start>=time>end

        Args:
            dtRound: the time delta of the data out of interest, default minute (60 second)

        Returns:
            index (bytearray): indicies for time of interest

        """
        raise NotImplementedError('please use master get time that is not a member of this class ')
        # TODO find a way to pull only hourly data or regular interval of desired time
        # todo this use date2index and create a list of dates see help(nc.date2index)
        try:
            self.ncfile = nc.Dataset(
                self.crunchDataLoc + self.dataloc)  # loads all of the netCDF file
            #            try:
            self.allEpoch = sb.baseRound(self.ncfile['time'][:],
                                         base=dtRound)  # round to nearest minute
            
            # now find the boolean!
            mask = (self.allEpoch >= self.epochd1) & (self.allEpoch < self.epochd2)
            idx = np.argwhere(mask).squeeze()
            # old slow way of doing time!
            # self.alltime = nc.num2date(self.cshore_ncfile['time'][:], self.cshore_ncfile[
            # 'time'].units,
            #                            self.cshore_ncfile['time'].calendar) # converts all
            #                            epoch time to datetime
            #                            objects
            # for bb, date in enumerate(self.alltime):  # rounds time to nearest
            #     self.alltime[bb] = self.roundtime(dt=date, roundto=dtRound)
            #
            # mask = (self.alltime >= self.start) & (self.alltime < self.end)  # boolean
            # true/false of time
            # if (np.argwhere(mask).squeeze() == idx).all():
            #     print '.... old Times match New Times' % np.argwhere(mask).squeeze()
            assert np.size(idx) > 0, 'no data locally, check CHLthredds'
            print("Data Gathered From Local Thredds Server")
        
        except (IOError, RuntimeError, NameError,
                AssertionError):  # if theres any error try to get good data from next location
            try:
                # MPG: Use self.chlDataLoc with 'frf/' removed from string for correct url.
                self.ncfile = nc.Dataset(self.chlDataLoc.replace('frf/', 'cmtb/') + self.dataloc)
                self.allEpoch = sb.baseRound(self.ncfile['time'][:], base=dtRound)  # round to nearest minute
                # now find the boolean !
                emask = (self.allEpoch >= self.epochd1) & (self.allEpoch < self.epochd2)
                idx = np.argwhere(emask).squeeze()
                
                # self.alltime = nc.num2date(self.cshore_ncfile['time'][:], self.cshore_ncfile[
                # 'time'].units,
                #                            self.cshore_ncfile['time'].calendar)
                # for bb, date in enumerate(self.alltime):
                #     self.alltime[bb] = self.roundtime(dt=date, roundto=dtRound)
                # # mask = (sb.roundtime(self.cshore_ncfile['time'][:]) >= self.epochd1) & (
                # sb.roundtime(
                # self.cshore_ncfile['time'][:]) < self.epochd2)\
                #
                # mask = (self.alltime >= self.start) & (self.alltime < self.end)  # boolean
                # true/false of time
                #
                # idx = np.argwhere(mask).squeeze()
                
                try:
                    assert np.size(
                        idx) > 0, ' There are no data within the search parameters for this gauge'
                    print("Data Gathered from CHL thredds Server")
                except AssertionError:
                    idx = None
            except IOError:  # this occors when thredds is down
                print(' Trouble Connecteing to data on CHL Thredds')
                idx = None
        
        self.ncfile = nc.Dataset(self.crunchDataLoc + self.dataloc)
        # switch us back to the local THREDDS if it moved us to CHL
        
        return idx
    
    def getGridCMS(self, method):
        """This function will grab data from the CMS grid folder on the server.
        
        This Function is depricated

        Args:
          method: can be [1, historical, history]  for historical
        can be [0, 'time'] for non oporational consideration

        Returns:
          key 'xCoord': x in FRF
          key 'yCoord': y in FRF
          key 'elevation': elevation NAVD 88
          key 'time': time in date time object
          key 'lat': latitude
          key 'lon': longitude
          key 'northing': NC stateplane Northing
          key 'easting': NC stateplane Easting
          key 'x0': origin in x (stateplane easting)
          key 'azimuth': grid orientation
          key 'y0': origin in y (stateplane northing)

        """
        self.dataloc = 'grids/CMSwave_v1/CMSwave_v1.ncml'
        self.ncfile, self.allEpoch = getnc(dataLoc=self.dataloc, callingClass=self.callingClass,
                                           dtRound=1 * 60)
        try:
            self.bathydataindex = gettime(allEpoch=self.allEpoch, epochStart=self.epochd1,
                                          epochEnd=self.epochd2)  # getting the index of the grid
        except IOError:
            self.bathydataindex = None
        if self.bathydataindex != None and np.size(self.bathydataindex) == 1:
            idx = self.bathydataindex.squeeze()
        elif (self.bathydataindex == None or len(self.bathydataindex) < 1) and method in [1,
                                                                                          'historical',
                                                                                          'history']:
            # there's no exact bathy match so find the max negative number where the negitive
            # numbers are historical and the max would be the closest historical
            val = (max([n for n in (self.ncfile['time'][:] - self.epochd1) if n < 0]))
            idx = np.where((self.ncfile['time'][:] - self.epochd1) == val)[0][0]
            print('Bathymetry is taken as closest in HISTORY - operational')
        elif (self.bathydataindex == None or np.size(self.bathydataindex) < 1) and method == 0:
            idx = np.argmin(np.abs(self.ncfile['time'][:] - self.start))  # closest in time
            print('Bathymetry is taken as closest in TIME - NON-operational')
        elif self.bathydataindex != None and len(self.bathydataindex) > 1:
            val = (max([n for n in (self.ncfile['time'][:] - self.start) if n < 0]))
            idx = np.where((self.ncfile['time'] - self.start) == val)[0][0]
            
            #
            # if self.bathydataindex is not None and  len(self.bathydataindex) == 1:
            #     idx = self.bathydataindex
            # elif self.bathydataindex is not None and len(self.bathydataindex) < 1 and method in
            # [1, 'historical',
            # 'history']:
            #     # there's no exact bathy match so find the max negative number where the negitive
            #     # numbers are historical and the max would be the closest historical
            #     val = (max([nHs for nHs in (self.cshore_ncfile['time'][:] - self.epochd1) if
            #     nHs < 0]))
            #     idx = np.where((self.cshore_ncfile['time'][:] - self.epochd1) == val)[0][0]
            #     print 'Bathymetry is taken as closest in HISTORY - operational'
            # elif self.bathydataindex is not None and len(self.bathydataindex) < 1 and method == 0:
            #     idx = np.argmin(np.abs(self.cshore_ncfile['time'][:] - self.start))  # closest
            #     in time
            #     print 'Bathymetry is taken as closest in TIME - NON-operational'
            # elif self.bathydataindex is not None and len(self.bathydataindex) > 1:
            #     val = (max([nHs for nHs in (self.cshore_ncfile['time'][:] - self.start) if nHs
            #     < 0]))
            #     idx = np.where((self.cshore_ncfile['time'] - self.start) == val)[0][0]
            
            print('The closest in history to your start date is %s\n' % nc.num2date(
                self.gridTime[idx],
                self.ncfile['time'].units))
            print('Please End new simulation with the date above')
            
            raise Exception
        if np.size(idx) > 0 and idx is not None:
            # now retrieve data with idx
            elevation_points = self.ncfile['elevation'][idx]
            xCoord = self.ncfile['xFRF'][:]
            yCoord = self.ncfile['yFRF'][:]
            lat = self.ncfile['latitude'][:]
            lon = self.ncfile['longitude'][:]
            northing = self.ncfile['northing'][:]
            easting = self.ncfile['easting'][:]
            
            time = nc.num2date(self.ncfile['time'][idx], self.ncfile['time'].units, only_use_cftime_datetimes=False)
            
            gridDict = {'xCoord':    xCoord,
                        'yCoord':    yCoord,
                        'elevation': elevation_points,
                        'time':      time,
                        'lat':       lat,
                        'lon':       lon,
                        'northing':  northing,
                        'easting':   easting,
                        'x0':        self.ncfile['x0'][:],
                        'azimuth':   self.ncfile['azimuth'][:],
                        'y0':        self.ncfile['y0'][:],
                        }
            return gridDict
    
    def getBathyIntegratedTransect(self, method=1, ForcedSurveyDate=None, **kwargs):
        """This function gets the integraated bathy, using the plant (2009) method.

        Args:
            method (int): a key which determines which method to find bathymetry with (Default
            value = 1)
                method == 1  - > 'Bathymetry is taken as closest in HISTORY - operational'

                method == 0  - > 'Bathymetry is taken as closest in TIME - NON-operational'

            ForcedSurveyDate (str): This is to force a date of survey gathering (Default value =
            None) if set to 'all'
                it will pull multiple bathys, this can choke up the server, and is recommended to
                set xbounds and
                ybounds

        Keyword Args:
           'cBKF': if true will get cBathy original Kalman Filter

           'cBKF_T': if true will get wave height thresholded Kalman filter

            'xbound': = [xmin, xmax]  which will truncate the cbathy domain to xmin, xmax (frf
            coord)

            'ybound': = [ymin, ymax]  which will truncate the cbathy domain to ymin, ymax (frf
            coord)

            'forceReturnAll' (bool): return all survey grids, not just nearest in time by method
            above

            'forceReturnAllPlusOne' (bool): return all survey grids, but also the one just before
            your inital survey
            date

        Returns:
          dictionary with keys
            'xFRF': x coordinate in FRF

            'yFRF': y coorindate in FRF

            'elevation': bathymetry

            'time': time in Datetime objects

            'lat': latitude

            'lon': longitude

            'northing': NC stateplane northing

            'easting': NC stateplane Easting

            'surveyNumber': FRF survey number (metadata)

        """
        forceReturnAll = kwargs.get('forceReturnAll', False)  # returns all survey
        forceReturnAllPlusOne = kwargs.get('forceReturnAllPlusOne', False)  # returns all surveys
        verbose = kwargs.get('verbose', False)
        
        if ForcedSurveyDate != None:
            # start is used in the gettime function,
            # to force a selection of survey date self.start/end is changed to the forced
            # survey date and then changed back using logged start/stop
            # a check is in place to ensure that the retieved time is == to the forced time
            oldD1 = self.start
            oldD2 = self.end
            oldD1epoch = self.epochd1
            oldD2epoch = self.epochd2
            self.start = ForcedSurveyDate  # change time one
            self.end = ForcedSurveyDate + DT.timedelta(0, 1)  # and time 2
            self.epochd1 = nc.date2num(self.start, 'seconds since 1970-01-01')
            self.epochd2 = nc.date2num(self.end, 'seconds since 1970-01-01')
            print('!!!Forced bathy date %s' % ForcedSurveyDate)
        
        ####################################################################
        #  Set URL based on Keyword, Default to surveyed bathymetry        #
        ####################################################################
        if 'cBKF_T' in kwargs and kwargs['cBKF_T'] == True:
            self.dataloc = 'integratedBathyProduct/cBKF-T/cBKF-T.ncml'
        elif 'cBKF' in kwargs and kwargs['cBKF'] == True:
            self.dataloc = 'integratedBathyProduct/cBKF/cBKF.ncml'
        else:
            self.dataloc = 'integratedBathyProduct/survey/survey.ncml'
        ####################################################################
        #   go get the index and return based on method chosen             #
        ####################################################################
        # go ahead and assign the ncfile first....
        self.ncfile, self.allEpoch = getnc(dataLoc=self.dataloc, callingClass=self.callingClass,
                                           dtRound=1 * 60, start=self.start, end=self.end)
        
        try:
            self.bathydataindex = gettime(allEpoch=self.allEpoch, epochStart=self.epochd1,
                                          epochEnd=self.epochd2)  # getting the index of the grid
        except IOError:
            self.bathydataindex = []  # when a server is not available
        
        if (forceReturnAll == True and self.bathydataindex is not None) or (
            np.size(self.bathydataindex) == 1 and self.bathydataindex != None):
            idx = self.bathydataindex.squeeze()

        elif forceReturnAllPlusOne == True and self.bathydataindex is not None:
            idx = np.append(self.bathydataindex.squeeze().min() - 1, self.bathydataindex.squeeze())

        elif np.size(self.bathydataindex) > 1:
            val = (max([n for n in (self.ncfile['time'][:] - self.epochd1) if n < 0]))
            # idx = np.where((self.ncfile['time'][:] - self.epochd1) == val)[0][0]
            idx = np.argmin(np.abs(self.ncfile['time'][:] - self.epochd1))  # closest in time
            warnings.warn('Pulled multiple bathymetries')
            print('   The nearest bathy to your simulation start date is %s' % nc.num2date(
                self.allEpoch[idx],
                self.ncfile['time'].units))
            print(
                '   Please End new simulation with the date above, so it does not pull multiple '
                'bathymetries')
            raise NotImplementedError
        
        elif (self.bathydataindex == None or len(self.bathydataindex) < 1) & method == 1:
            # there's no exact bathy match so find the max negative number where the negitive
            # numbers are historical and the max would be the closest historical
            val = (max([n for n in (self.ncfile['time'][:] - self.epochd1) if n < 0]))
            idx = np.where((self.ncfile['time'][:] - self.epochd1) == val)[0][0]
            print('Bathymetry is taken as closest in HISTORY - operational')
        elif (self.bathydataindex == None or np.size(self.bathydataindex) < 1) and method == 0:
            idx = np.argmin(np.abs(self.ncfile['time'][:] - self.epochd1))  # closest in time
            print('Bathymetry is taken as closest in TIME - NON-operational')
        elif self.bathydataindex != None and len(self.bathydataindex) > 1:
            val = (max([n for n in (self.ncfile['time'][:] - self.epochd1) if n < 0]))
            idx = np.where((self.ncfile['time'][:] - self.epochd1) == val)[0][0]
            
            print('The closest in history to your start date is %s\n' % nc.num2date(
                self.ncfile['time'][idx],
                self.ncfile['time'].units))
            print('Please End new simulation with the date above')
            raise Exception

        # after all that logic, if it's still None, just return it
        if idx is None:
            return None
        ###############################################################################################################
        # bound it if requested
        ###############################################################################################################
        if 'xbounds' in kwargs and np.array(kwargs['xbounds']).size == 2:
            if kwargs['xbounds'][0] > kwargs['xbounds'][1]:
                kwargs['xbounds'] = np.flip(kwargs['xbounds'], axis=0)
            # first min of x
            if (kwargs['xbounds'][0] < self.ncfile['xFRF'][:]).all():
                # then set xmin to 0
                removeMinX = 0
            else:  # <= used here to handle inclusive initial index inherant in python
                removeMinX = np.argwhere(
                    self.ncfile['xFRF'][:] <= kwargs['xbounds'][0]).squeeze().max()
            # now max of x
            if (kwargs['xbounds'][1] > self.ncfile['xFRF'][:]).all():
                removeMaxX = None
            else:
                removeMaxX = np.argwhere(
                    self.ncfile['xFRF'][:] >= kwargs['xbounds'][
                        1]).squeeze().min() + 1  # python indexing
            xs = slice(removeMinX, removeMaxX)
        else:
            xs = slice(None)
        
        if 'ybounds' in kwargs and np.array(kwargs['ybounds']).size == 2:
            if kwargs['ybounds'][0] > kwargs['ybounds'][1]:
                kwargs['ybounds'] = np.flip(kwargs['ybounds'], axis=0)
            # first min of y
            if (kwargs['ybounds'][0] < self.ncfile['yFRF'][:]).all():
                # then set the ymin to first index [0]
                removeMinY = 0  # ie get all data
            else:
                removeMinY = np.argwhere(
                    self.ncfile['yFRF'][:] <= kwargs['ybounds'][0]).squeeze().max()
            ## now max of y
            if (kwargs['ybounds'][1] > self.ncfile['yFRF'][:]).all():
                removeMaxY = None
            else:
                removeMaxY = np.argwhere(
                    self.ncfile['yFRF'][:] >= kwargs['ybounds'][
                        1]).squeeze().min() + 1  # python indexing
            ys = slice(removeMinY, removeMaxY)
        else:
            ys = slice(None)
        ###############################################################################################################
        ###############################################################################################################
        # the below line was in place, it should be masking nan's but there is not supposed to be
        # nan's
        # in the data, should only be fill values (-999)
        # elevation_points = np.ma.array(cshore_ncfile['elevation'][idx,:,:], mask=np.isnan(
        # cshore_ncfile[
        # 'elevation'][idx,:,:]))
        # remove -999's

        elevation_points = self.ncfile['elevation'][idx, ys, xs]
        updateTime = self.ncfile['updateTime'][idx, ys, xs]
        xCoord = self.ncfile['xFRF'][xs]
        yCoord = self.ncfile['yFRF'][ys]
        lat = self.ncfile['latitude'][ys, xs]
        lon = self.ncfile['longitude'][ys, xs]
        
        # putting dates and times back for all the other instances that use get time
        if ForcedSurveyDate != None:
            self.start = oldD1
            self.end = oldD2
            self.epochd2 = oldD2epoch
            self.epochd1 = oldD1epoch
        
        bathyT = nc.num2date(self.allEpoch[idx],
                             'seconds since 1970-01-01', only_use_cftime_datetimes=False)  # This one is rounded
        # appropraitely
        # this comes directly from file (useful if server is acting funny)
        # bathyT = nc.num2date(self.ncfile['time'][idx], 'seconds since 1970-01-01')
        
        if verbose: print('  Measured Bathy is %s old' % (self.end - bathyT))
        
        gridDict = {'xFRF':      xCoord,
                    'yFRF':      yCoord,
                    'elevation': elevation_points,
                    'time':      bathyT,
                    'lat':       lat,
                    'lon':       lon, }
        if ('cBKF_T' not in kwargs) and (f'cBKF' not in kwargs):  # then its a survey, get the survey number
            gridDict['surveyNumber'] = self.ncfile['surveyNumber'][idx]
        
        return gridDict
    
    def getStwaveField(self, var, prefix, local=True, ijLoc=None, model='STWAVE'):
        """Depricated."""
        warnings.warn('Using depricated function name: getStwaveField')
        return self.getModelField(var, prefix, local, ijLoc, model)
    
    def getModelField(self, var, prefix, local=True, ijLoc=None, model='STWAVE', **kwargs):
        """Retrives data from spatial data CMSWave and STWAVE model.

        Args:
            local (bool): defines whether the data is from the nested simulation or the regional
            simulation (Default
            value = True)

            ijLoc (tuple, int): x or y or (x,y) tuple of location of interest in FRF Coordinates
                if None, will grab everything (expensive) (Default value = None)

            var (str): which variable to get from the spatial data ncml file

            prefix (str): this dictates which model run data are retrieved from

            local (bool): pull from the nested grid or the regional grid (Default value = True)

            model (str): one of: STWAVE, CMS (other models can be added)
        Keyword Args:
            xbound: = [xmin, xmax]  which will truncate the cbathy domain to xmin, xmax (frf coord)

            ybound: = [ymin, ymax]  which will truncate the cbathy domain to ymin, ymax (frf coord)

        Returns:
            a dictionary with keys below, see netCDF file for more metadata

            'time' (obj):  date time object

            'epochtime' (float): epoch time

            'var' (str): variable of interest as put in to the function

            'xFRF' (int): x location of data

            'yFRF' (int): y location of data

        """
        if local == True:
            grid = 'Local'
        elif local == False:
            grid = 'Regional'
        ############## setting up the ncfile ############################
        if prefix == 'CBHPStatic' and local == True:  # this is needed because projects are
            # stored in weird place
            fname = 'http://bones/thredds/dodsC/CMTB/projects/bathyDuck_SingleBathy_CBHP' \
                    '/Local_Field/Local_Field.ncml'
        elif prefix == 'CBHPStatic' and local == False:
            fname = 'http://bones/thredds/dodsC/CMTB/projects/bathyDuck_SingleBathy_CBHP' \
                    '/Regional_Field' \
                    '/Regional_Field.ncml'
        elif model == 'STWAVE':  # this is standard operational model url Structure
            fname = self.crunchDataLoc + u'waveModels/%s/%s/%s-Field/%s-Field.ncml' % (
                model, prefix, grid, grid)
        elif model == 'CMS':  # this is standard operational model url Structure
            fname = self.crunchDataLoc + u'waveModels/%s/%s/Field/Field.ncml' % (model, prefix)
        finished = False
        n = 0
        while not finished and n < 15:
            try:
                ncfile = nc.Dataset(fname)
                finished = True
            except IOError:
                print('Error reading {}, trying again'.format(fname))
                time.sleep(10)
                n += 1
        if not finished:
            raise (RuntimeError, 'Data not accessible right now')
        
        assert var in ncfile.variables.keys(), 'variable called is not in file please use\n%s' % \
                                               ncfile.variables.keys()
        
        mask = (ncfile['time'][:] >= nc.date2num(self.start, ncfile['time'].units)) & (
            ncfile['time'][:] <= nc.date2num(self.end, ncfile['time'].units))
        idx = np.where(mask)[0]
        assert np.size(idx > 0), " there's no data"
        if model == 'STWAVE':
            print('getting %s %s  %s %s Data' % (prefix, model, grid, var))
        elif model == 'CMS':
            print('getting %s %s  %s Data' % (prefix, model, var))
        ##############################################################################################################
        # now creating tool to remove single data point
        if ijLoc != None:
            assert len(
                ijLoc) == 2, 'if giving a postion, must be a tuple of i, j location (of length 2)'
            if type(ijLoc[0]) == int:
                x = ncfile[var].shape[1] - ijLoc[
                    0]  # the data are stored with inverse indicies to grid node locations
                y = ijLoc[1]  # use location given by function call
            else:  # ijLoc[0] == slice:
                x = ijLoc[0]
                y = np.argmin(np.abs(ncfile['yFRF'][:] - ijLoc[1]))
        
        else:
            x = slice(None)  # take entire data
            y = slice(None)  # take entire data
        ############################################ xbounds/ybounds
        # #################################################
        ###### sub divide bounds by x and y to get subdomain
        if 'xbounds' in kwargs and np.array(kwargs['xbounds']).size == 2:
            if kwargs['xbounds'][0] > kwargs['xbounds'][1]:
                kwargs['xbounds'] = np.flip(kwargs['xbounds'], axis=0)
            # first min of x
            if (kwargs['xbounds'][0] < ncfile['xFRF'][:]).all():
                # then set xmin to 0
                removeMinX = 0
            else:  # <= used here to handle inclusive initial index inherant in python
                removeMinX = np.argwhere(ncfile['xFRF'][:] <= kwargs['xbounds'][0]).squeeze().max()
            # now max of x
            if (kwargs['xbounds'][1] > ncfile['xFRF'][:]).all():
                removeMaxX = None
            else:
                removeMaxX = np.argwhere(
                    ncfile['x'][:] >= kwargs['xbounds'][1]).squeeze().min() + 1  # python indexing
            x = slice(removeMinX, removeMaxX)
        else:
            x = slice(None)
        
        if 'ybounds' in kwargs and np.array(kwargs['ybounds']).size == 2:
            if kwargs['ybounds'][0] > kwargs['ybounds'][1]:
                kwargs['ybounds'] = np.flip(kwargs['ybounds'], axis=0)
            # first min of y
            if (kwargs['ybounds'][0] < ncfile['yFRF'][:]).all():
                # then set the ymin to first index [0]
                removeMinY = 0  # ie get all data
            else:
                removeMinY = np.argwhere(ncfile['yFRF'][:] <= kwargs['ybounds'][0]).squeeze().max()
            ## now max of y
            if (kwargs['ybounds'][1] > ncfile['yFRF'][:]).all():
                removeMaxY = None
            else:
                removeMaxY = np.argwhere(
                    ncfile['yFRF'][:] >= kwargs['ybounds'][
                        1]).squeeze().min() + 1  # python indexing
            y = slice(removeMinY, removeMaxY)
        else:
            y = slice(None)
        ################################################################################################################
        if ncfile[var].ndim > 2 and idx.shape[0] > 100:  # looping through ... if necessary
            list = np.round(
                np.linspace(idx.min(), idx.max(), idx.max() - idx.min(), endpoint=True, dtype=int))
            # if idx.max() not in list:
            #     list = np.append(list, idx.max())
            xFRF = ncfile['xFRF'][x]
            yFRF = ncfile['yFRF'][y]
            if len(list) < 100:
                dataVar = np.array(ncfile[var][np.squeeze(list)])
                timeVar = nc.num2date(ncfile['time'][np.squeeze(list)], ncfile['time'].units,
                                      only_use_cftime_datetimes=False)
            else:
                for num, minidx in enumerate(list):
                    if len(list) < 100:
                        dataVar = np.array(ncfile[var][np.squeeze(list)])
                        timeVar = nc.num2date(ncfile['time'][np.squeeze(list)],
                                              ncfile['time'].units,
                                              only_use_cftime_datetimes=False)
                    elif num == 0:
                        dataVar = np.array(ncfile[var][range(minidx, list[num + 1]), y, x])
                        timeVar = nc.num2date(
                            np.array(ncfile['time'][range(minidx, list[num + 1])]),
                            ncfile['time'].units,
                            only_use_cftime_datetimes=False)
                    elif minidx == list[-1]:
                        lastIdx = (idx - minidx)[(idx - minidx) >= 0] + minidx
                        dataVar = np.append(dataVar, ncfile[var][lastIdx, y, x], axis=0)
                        timeVar = np.append(timeVar, nc.num2date(ncfile['time'][lastIdx],
                                                                 ncfile['time'].units,
                                                                 only_use_cftime_datetimes=False), axis=0)
                    else:
                        dataVar = np.append(dataVar,
                                            ncfile[var][range(minidx, list[num + 1]), y, x], axis=0)
                        timeVar = np.append(timeVar, nc.num2date(
                            ncfile['time'][range(minidx, list[num + 1])],
                            ncfile['time'].units,
                            only_use_cftime_datetimes=False), axis=0)
        
        elif ncfile[var].ndim > 2:
            dataVar = ncfile[var][idx, y, x]
            xFRF = ncfile['xFRF'][x]
            yFRF = ncfile['yFRF'][y]
            timeVar = nc.num2date(ncfile['time'][np.squeeze(idx)], ncfile['time'].units,
                                  only_use_cftime_datetimes=False)
        else:  # probably bathymetry date variable
            dataVar = ncfile[var][idx]
            xFRF = ncfile['xFRF'][x]
            yFRF = ncfile['yFRF'][y]
            timeVar = nc.num2date(ncfile['time'][np.squeeze(idx)], ncfile['time'].units,
                                  only_use_cftime_datetimes=False)
        # package for output
        field = {'time':      timeVar,
                 'epochtime': ncfile['time'][idx],  # pulling down epoch time of interest
                 var:         dataVar,
                 'xFRF':      xFRF,
                 'yFRF':      yFRF,
                 }
        try:
            field['bathymetryDate'] = ncfile['bathymetryDate'][idx]
        except IndexError:
            field['bathymetryDate'] = np.ones_like(field['time'])
        
        assert field[var].shape[0] == len(
            field['time']), " the indexing is wrong for pulling down the spatial output"
        field = removeDuplicatesFromDictionary(field)
        return field
    
    def getWaveSpecSTWAVE(self, prefix, gaugenumber, local=True, model='STWAVE'):
        """Depricated."""
        warnings.warn('Using depricated function name')
        return self.getWaveSpecModel(prefix, gaugenumber, model)
    
    def getWaveSpecModel(self, prefix, gaugenumber, model='STWAVE', removeBadWLFlag=True):
        """This function pulls down the data from the thredds server and puts the data into proper places.
        
        To be read for STwave Scripts this will return the wavespec with dir/freq bin and directionalWaveGaugeList wave
            energy

        Args:
            prefix (str): a 'key' to select which version of the simulations to pull data from
                available values are listed in the table below
                ['CB', 'HP', 'CBHP', 'FP', 'S' %any date string, 'CBThresh_0']

            local (bool): this denotes whether the waves are pulled from the nested domian or (
            Default value = True)

            gaugenumber (int, str): keys associated with data
                26m waverider can be [0, 'waverider-26m', 'Waverider-26m', '26m']
                17m waverider can be [1, 'Waverider-17m', 'waverider-17m']
                11m AWAC      can be [2, 'AWAC-11m', 'awac-11m', 'Awac-11m']
                8m array      can be [3, '8m-Array', '8m Array', '8m array', '8m-array']
                6m AWAC       can be [4, 'awac-6m', 'AWAC-6m']
                4.5m AWAC     can be [5, 'awac-4.5m', 'Awac-4.5m']
                3.5m aquadopp can be [6, 'adop-3.5m', 'aquadopp 3.5m']
                200m pressure can be [8, 'xp200m', 'xp200']
                150m pressure can be [9, 'xp150m', 'xp150']
                125m pressure can be [10, 'xp125m', 'xp125']

            model (str): one of: STWAVE, CMS
            
            removeBadWLFlag(bool): run logic to remove bad data from returned dictionary (default=True)

        Returns: return dictionary with packaged data following keys
            'epochtime': time in epoch ('second since 1970-01-01

            'time': time in date time object

            'name': name

            'wavefreqbin': wave frequency bins

            'Hs': wave Height

            'peakf': peak frequencies or 1/waveTp

            'wavedirbin': wave direction bins for dwed

            'dWED': directionalWaveGaugeList wave energy density - 2d spectra [t, freq, dir]

            'waveDm': wave mean direction

            'waveTm': wave mean period

            'waveTp': wave peak period

            'WL': water level (NAVD 88)

            'Umag': wind speed [m/s]

            'Udir': wind direction (True north)

            'fspec': frequency spectra [t, nfreq]

            'qcFlag': qc flags


        """
        # Making gauges flexible
        if prefix in ['CB', 'HP', 'CBHP', 'FP', 'CBThresh']:
            urlFront = 'waveModels/%s/%s' % (model, prefix)
        elif prefix.startswith('S') and prefix[1].isdigit():  # this is static bathy
            urlFront = 'projects/%s/CBHP/SingleBathy_%s' % (model, prefix[1:])
        elif prefix in ['CBThresh_0']:
            urlFront = 'projects/%s/%s' % (model, prefix)
        elif prefix.lower() in ['cbthresh_0_oversmoothed']:
            urlFront = 'projects/%s/CBThresh_0_oversmoothed' % model
        elif prefix.lower() in ['cbthresh_0_papersubmittedv1']:
            urlFront = 'projects/%s/CBThresh_0_paperSubmittedV1' % model
        
        ############### now identify file name #################
        if gaugenumber in [0, 'waverider-26m', 'Waverider-26m', '26m']:
            # 26 m wave rider
            fname = 'waverider-26m/waverider-26m.ncml'
            gname = '26m Waverider Buoy'
        elif gaugenumber in [1, 'Waverider-17m', 'waverider-17m']:
            # 2D 17m waverider
            fname = 'waverider-17m/waverider-17m.ncml'
            gname = '17m Waverider Buoy'
        elif gaugenumber in [2, 'AWAC-11m', 'awac-11m', 'Awac-11m']:
            gname = 'AWAC04 - 11m'
            fname = 'awac-11m/awac-11m.ncml'
        elif gaugenumber in [3, 'awac-8m', 'AWAC-8m', 'Awac-8m', 'awac 8m',
                             '8m-Array', '8m Array', '8m array', '8m-array']:
            gname = 'AWAC 8m'
            fname = '8m-array/8m-array.ncml'
        elif gaugenumber in [4, 'awac-6m', 'AWAC-6m']:
            gname = 'AWAC 6m'
            fname = 'awac-6m/awac-6m.ncml'
        elif gaugenumber in [5, 'awac-4.5m', 'Awac-4.5m']:
            gname = 'AWAC 4.5m'
            fname = 'awac-4.5m/awac-4.5m.ncml'
        elif gaugenumber in [6, 'adop-3.5m', 'aquadopp 3.5m']:
            gname = 'Aquadopp 3.5m'
            fname = 'adop-3.5m/adop-3.5m.ncml'
        elif gaugenumber in [8, 'xp200m', 'xp200']:
            gname = 'Paros xp200m'
            fname = 'xp200m/xp200m.ncml'
        elif gaugenumber in [9, 'xp150m', 'xp150']:
            gname = 'Paros xp150m'
            fname = 'xp150m/xp150m.ncml'
        elif gaugenumber == 10 or gaugenumber == 'xp125m':
            gname = 'Paros xp125m'
            fname = 'xp125m/xp125m.ncml'
        elif gaugenumber == 11 or gaugenumber == 'xp100m':
            gname = 'Paros xp100m'
            fname = 'xp100m/xp100m.ncml'
        else:
            gname = 'There Are no Gauge numbers here'
            raise NameError('Bad Gauge name, specify proper gauge name/number')
        # parsing out data of interest in time
        self.dataloc = urlFront + '/' + fname
        self.ncfile, self.allEpoch = getnc(dataLoc=self.dataloc, callingClass=self.callingClass, dtRound=1 * 60)
        try:
            # go get indices of interest
            self.wavedataindex = gettime(allEpoch=self.allEpoch, epochStart=self.epochd1,
                                         epochEnd=self.epochd2)
            assert np.array(
                self.wavedataindex).all() != None, 'there''s no data in your time period'
            if np.size(self.wavedataindex) >= 1:
                wavespec = {'epochtime':   self.ncfile['time'][self.wavedataindex],
                            'time':        nc.num2date(self.allEpoch[self.wavedataindex],
                                                       self.ncfile['time'].units,
                                                       only_use_cftime_datetimes=False),
                            'name':        nc.chartostring(self.ncfile['station_name'][:]),
                            'wavefreqbin': self.ncfile['waveFrequency'][:],
                            # 'lat': self.ncfile['lat'][:],
                            # 'lon': self.ncfile['lon'][:],
                            'Hs':          self.ncfile['waveHs'][self.wavedataindex],
                            'peakf':       self.ncfile['waveTp'][self.wavedataindex],
                            'wavedirbin':  self.ncfile['waveDirectionBins'][:],
                            'dWED':        self.ncfile['directionalWaveEnergyDensity'][self.wavedataindex,
                                           :, :],
                            'waveDm':      self.ncfile['waveDm'][self.wavedataindex],
                            'waveTm':      self.ncfile['waveTm'][self.wavedataindex],
                            'waveTp':      self.ncfile['waveTp'][self.wavedataindex],
                            'WL':          self.ncfile['waterLevel'][self.wavedataindex],
                            'qcFlagWL':    self.ncfile['qcFlag'][self.wavedataindex, 2],
                            'qcFlagWind':  self.ncfile['qcFlag'][self.wavedataindex, 1]}
                # make frequency spectra from directional energy spectrum
                wavespec['fspec'] = wavespec['dWED'].sum(axis=2) * np.median(
                    np.diff(np.array(wavespec['wavedirbin'])))
                if model == 'STWAVE':
                    wavespec['Umag'] = self.ncfile['Umag'][self.wavedataindex]
                    wavespec['Udir'] = self.ncfile['Udir'][self.wavedataindex]
                wavespec['dWED'][wavespec['dWED'] == 0] = 1e-6
                wavespec['fspec'][wavespec['fspec'] == 0] = 1e-6
            qcFlags = self.ncfile['qcFlag'][self.wavedataindex]
            if removeBadWLFlag is not False:
                idxGood = np.argwhere(qcFlags[:, 2] <= 5).squeeze()
                wavespec = sb.reduceDict(wavespec, idxGood)
            return removeDuplicatesFromDictionary(wavespec)
        
        except (RuntimeError, AssertionError) as err:
            print(err)
            print('ERROR Retrieving data from %s\n in this time period start: %s  End: %s' % (
                gname, self.start, self.end))
            return None
    
    def getCSHOREOutput(self, prefix):
        """Retrieves data from spatial data CSHORE model.
        
        Args:
            prefix (str): a 'key' to select which version of the simulations to pull data from
                available value is only 'MOBILE_RESET' for now but more could be
                added in the future

        Returns:
            dictionary with packaged data following keys:
                'epochtime' (float):  epoch time
    
                'time' (obj): datetime of model output
    
                'xFRF' (float): x location of data
    
                'Hs' (float): significant wave height
    
                'zb' (float): bed elevation
    
                'WL' (float): water level
    
                'bathyTime' (ojb): datetime of bathymetric survey used
    
                'setup' (float): wave induced setup height
    
                'aveN' (float): average northward current
    
                'stdN' (float): standard deviation of northward current
    
                'runupMean' (float): mean runup elevation
    
                'runup2perc' (float): 2% runup elevation

        """
        dataLoc = 'morphModels/CSHORE/{0}/{0}.ncml'.format(prefix)
        ncfile, allEpoch = getnc(dataLoc, self.THREDDS, self.callingClass)
        dataIndex = gettime(allEpoch, epochStart=self.epochd1, epochEnd=self.epochd2)
        if dataIndex is None:
            print(('There\'s no data in time period ' + self.start.strftime('%Y-%m-%dT%H%M%SZ') +
                   ' to ' + self.end.strftime('%Y-%m-%dT%H%M%SZ')))
            return {}
        if isinstance(ncfile['bottomElevation'][dataIndex, :], np.ma.masked_array):
            dataIndex = dataIndex[~ncfile['bottomElevation'][dataIndex, :].mask.any(1)]
        
        if len(dataIndex) == 0:
            print(('There\'s no data in time period ' + self.start.strftime('%Y-%m-%dT%H%M%SZ') +
                   ' to ' + self.end.strftime('%Y-%m-%dT%H%M%SZ')))
            return {}
        mod = {'epochtime':  ncfile['time'][dataIndex],
               'time':       nc.num2date(ncfile['time'][dataIndex], ncfile['time'].units,
                                         only_use_cftime_datetimes=False),
               'xFRF':       ncfile['xFRF'][:],
               'Hs':         ncfile['waveHs'][dataIndex, :],
               'zb':         ncfile['bottomElevation'][dataIndex, :].data,
               'WL':         ncfile['waterLevel'][dataIndex, :],
               'bathyTime':  nc.num2date(ncfile['bathymetryDate'][dataIndex],
                                         ncfile['bathymetryDate'].units,
                                         only_use_cftime_datetimes=False),
               'setup':      ncfile['setup'][dataIndex, :],
               'aveN':       ncfile['aveN'][dataIndex, :],
               'stdN':       ncfile['stdN'][dataIndex, :],
               'runupMean':  ncfile['runupMean'][dataIndex],
               'runup2perc': ncfile['runup2perc'][dataIndex]}
        return mod
