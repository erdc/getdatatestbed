# -*- coding: utf-8 -*-
"""
Created on Tue Jul 21 13:38:05 2015
This is a class definition designed to get data from the FRF thredds server 

@author: Spicer Bak, PhD
@contact: spicer.bak@usace.army.mil
@organization: USACE CHL FRF


"""
import datetime as DT
import sys
from subprocess import check_output
import collections
import netCDF4 as nc
import numpy as np
import pandas as pd
from sblib import sblib as sb
from sblib import geoprocess as gp

# MPG: import cPickle for file i/o. 
import cPickle as pickle

class getObs:
    """
    Note d1 and d2 have to be in date-time formats
    are all data set times in UTC?
    need to write error handling, what to do if there's no data ??

    """

    def __init__(self, d1, d2):
        """
        Initialization description here
        Data are returned in self.datainex are inclusive at d1,
         exclusive at d2
        """
        # this is active wave gauge list for doing wave rider
        self.gaugelist = [
            'waverider-17m', 
            'awac-11m', 
            '8m-array', 
            'awac-6m', 
            'awac-4.5m', 
            'adop-3.5m', 
            'xp200m', 
            'xp150m', 
            'xp125m', 
            'waverider-26m'
            ]
        self.directional = ['waverider-26m', 'waverider-17m', 'awac-11m', '8m-array', 'awac-6m', 'awac-4.5m',
                            'adop-3.5m']
        self.rawdataloc_wave = []
        self.outputdir = []  # location for outputfiles
        self.d1 = d1  # start date for data grab
        self.d2 = d2  # end data for data grab
        self.timeunits = 'seconds since 1970-01-01 00:00:00'
        self.epochd1 = nc.date2num(self.d1, self.timeunits)
        self.epochd2 = nc.date2num(self.d2, self.timeunits)
        self.comp_time()
        self.FRFdataloc = u'http://134.164.129.55/thredds/dodsC/FRF/'
        self.crunchDataLoc = u'http://134.164.129.62:8080/thredds/dodsC/CMTB'
        self.chlDataLoc = u'https://chlthredds.erdc.dren.mil/thredds/dodsC/frf/' #'http://10.200.23.50/thredds/dodsC/frf/'
        assert type(self.d2) == DT.datetime, 'd1 need to be in python "Datetime" data types'
        assert type(self.d1) == DT.datetime, 'd2 need to be in python "Datetime" data types'

    def comp_time(self):
        """
        Test if times are backwards
        """
        assert self.d2 >= self.d1, 'finish time: d2 needs to be after start time: d1'

    def roundtime(self, dt=None, roundto=60):
        """
        Round a datetime object to any time laps in seconds
        Author: Thierry Husson 2012 - Use it as you want but don't blame me.
        :rtype: object

        :param dt:
            datetime.datetime object, default now.
        :param roundto:
            Closest number of SECONDS to round to, default 1 minute
        """
        if dt is None:
            dt = DT.datetime.now()
        seconds = (dt - dt.min).seconds
        # // is a floor division, not a comment on following line:
        rounding = (seconds + roundto / 2) // roundto * roundto
        return dt + DT.timedelta(0, rounding - seconds, -dt.microsecond)

    def gettime(self, dtRound=60):
        """
        this function opens the netcdf file, pulls down all of the time, then pulls the dates of interest
        from the THREDDS (data loc) server based on d1,d2, and data location
        it returns the indicies in the NCML file of the dates d1>=time>d2
        INPUTS:

             :param dtRound: the time delta of the data out of interest IN SECONDS, default minute (60 second)

        """
        # TODO find a way to pull only hourly data or regular interval of desired time
        # todo this use date2index and create a list of dates see help(nc.date2index)
        try:

            self.ncfile = nc.Dataset(self.FRFdataloc + self.dataloc) #loads all of the netCDF file
            #            try:
            self.allEpoch = sb.myround(self.ncfile['time'][:], base=dtRound) # round to nearest minute
            # now find the boolean!
            mask = (self.allEpoch >= self.epochd1) & (self.allEpoch < self.epochd2)
            idx = np.argwhere(mask).squeeze()
            # old slow way of doing time!
            # self.alltime = nc.num2date(self.ncfile['time'][:], self.ncfile['time'].units,
            #                            self.ncfile['time'].calendar) # converts all epoch time to datetime objects
            # for i, date in enumerate(self.alltime):  # rounds time to nearest
            #     self.alltime[i] = self.roundtime(dt=date, roundto=dtRound)
            #
            # mask = (self.alltime >= self.d1) & (self.alltime < self.d2)  # boolean true/false of time
            # if (np.argwhere(mask).squeeze() == idx).all():
            #     print '.... old Times match New Times' % np.argwhere(mask).squeeze()
            assert np.size(idx) > 0, 'no data locally, check CHLthredds'
            print "Data Gathered From Local Thredds Server"

        except (IOError, RuntimeError, NameError, AssertionError):  # if theres any error try to get good data from next location
            try:
                self.ncfile = nc.Dataset(self.chlDataLoc + self.dataloc)
                self.allEpoch = sb.myround(self.ncfile['time'][:], base=dtRound) # round to nearest minute
                # now find the boolean !
                emask = (self.allEpoch >= self.epochd1) & (self.allEpoch < self.epochd2)
                idx = np.argwhere(emask).squeeze()

                # self.alltime = nc.num2date(self.ncfile['time'][:], self.ncfile['time'].units,
                #                            self.ncfile['time'].calendar)
                # for i, date in enumerate(self.alltime):
                #     self.alltime[i] = self.roundtime(dt=date, roundto=dtRound)
                # # mask = (sb.roundtime(self.ncfile['time'][:]) >= self.epochd1) & (sb.roundtime(self.ncfile['time'][:]) < self.epochd2)\
                #
                # mask = (self.alltime >= self.d1) & (self.alltime < self.d2)  # boolean true/false of time
                #
                # idx = np.argwhere(mask).squeeze()


                try:
                    assert np.size(idx) > 0, ' There are no data within the search parameters for this gauge'
                    print "Data Gathered from CHL thredds Server"
                except AssertionError:
                    idx = None
            except IOError:  # this occors when thredds is down
                print ' Trouble Connecteing to data on CHL Thredds'
                idx = None

        return idx

    def getWaveSpec(self, gaugenumber=0, roundto=30):
        """
        This function pulls down the data from the thredds server and puts the data into proper places
        to be read for STwave Scripts
        this will return the wavespec with dir/freq bin and directional wave energy

        TO DO:
        Set optional date input from function arguments

        :param gaugenumber:
            gaugenumber = 0, 26m wave rider
            gaugenumber = 1, 17m waverider
            gaugenumber = 2, 'awac-11m'
            gaugenumber = 3, awac3 - 8m
            gaugenumber = 4, awac2 - 6m
            gaugenumber = 5, awac1 - 5m
            gaugenumber = 6, adopp2 - 3m
            gaugenumber = 7, adopp1 - 2m
            gaugenumber = 8,  Paros xp200m
            gaugenumber = 9,  Paros xp150m
            gaugenumber = 10, Paros xp125m
            gaugenumber = 11, Paros xp100m
            gaugenumber = 12, 8 m array
        :param collectionlength:
            s the time over which the wind record exists
            ie data is collected in 10 minute increments time is rounded to nearest 10min increment
            data is rounded to the nearst [collectionlength] (default 30 min)
        """
        # Making gauges flexible

        if gaugenumber in [0, 'waverider-26m', 'Waverider-26m', '26m']:
            # 26 m wave rider
            self.dataloc = 'oceanography/waves/waverider-26m/waverider-26m.ncml'  # 'oceanography/waves/waverider430/waverider430.ncml'  # 26m buoy
            gname = '26m Waverider Buoy'
        elif gaugenumber == 1 or gaugenumber == 'waverider-17m':
            # 2D 17m waverider
            self.dataloc = 'oceanography/waves/waverider-17m/waverider-17m.ncml'  # 17 m buoy
            gname = '17m Waverider Buoy'
        elif gaugenumber == 2 or gaugenumber == 'awac-11m':
            gname = 'AWAC 11m'
            self.dataloc = 'oceanography/waves/awac-11m/awac-11m.ncml'
        elif gaugenumber == 3 or gaugenumber == 'awac-8m':
            gname = 'AWAC 8m'
            self.dataloc = 'oceanography/waves/awac-8m/awac-8m.ncml'
        elif gaugenumber == 4 or gaugenumber == 'awac-6m':
            gname = 'AWAC 6m'
            self.dataloc = 'oceanography/waves/awac-6m/awac-6m.ncml'
        elif gaugenumber  in [5, 'awac-4.5m', 'awac_4.5m']:
            gname = 'AWAC 4.5m'
            self.dataloc = 'oceanography/waves/awac-4.5m/awac-4.5m.ncml'
        elif gaugenumber == 6 or gaugenumber == 'adop-3.5m':
            gname = 'Aquadopp 3.5m'
            self.dataloc = 'oceanography/waves/adop-3.5m/adop-3.5m.ncml'
        elif gaugenumber == 7 or gaugenumber == 'adop-2m':
            gname = 'Aquadopp01 - 2m'
            self.dataloc = 'oceanography/waves/adop01/adop01.ncml'
        elif gaugenumber == 8 or gaugenumber == 'xp200m':
            gname = 'Paros xp200m'
            self.dataloc = 'oceanography/waves/xp200m/xp200m.ncml'
        elif gaugenumber == 9 or gaugenumber == 'xp150m':
            gname = 'Paros xp150m'
            self.dataloc = 'oceanography/waves/xp150m/xp150m.ncml'
        elif gaugenumber == 10 or gaugenumber == 'xp125m':
            gname = 'Paros xp125m'
            self.dataloc = 'oceanography/waves/xp125m/xp125m.ncml'
        elif gaugenumber == 11 or gaugenumber == 'xp100m':
            gname = 'Paros xp100m'
            self.dataloc = 'oceanography/waves/xp100m/xp100m.ncml'
        elif gaugenumber == 12 or gaugenumber == '8m-array':
            gname = "8m array"
            self.dataloc = 'oceanography/waves/8m-array/8m-array.ncml'
        elif gaugenumber in ['oregonInlet', 'OI', 'oi']:
            gname = 'Oregon Inlet'
            self.dataloc = 'oceanography/waves/waverider-oregon-inlet-nc/waverider-oregon-inlet-nc.ncml'
        else:
            gname = 'There Are no Gauge numbers here'
            raise NameError('Bad Gauge name, specify proper gauge name/number')
        # parsing out data of interest in time
        try:
            self.wavedataindex = self.gettime(dtRound=roundto * 60)
            assert np.array(self.wavedataindex).all() != None, 'there''s no data in your time period'
            if np.size(self.wavedataindex) >= 1:
                # consistant for all wave gauges
                if np.size(self.wavedataindex) == 1:
                    self.wavedataindex = np.expand_dims(self.wavedataindex, axis=0)
                self.snaptime = nc.num2date(self.allEpoch[self.wavedataindex], self.ncfile['time'].units)
                try:
                    depth = self.ncfile['nominalDepth'][:]  # this should always go
                except IndexError:
                    depth = self.ncfile['gaugeDepth'][:]  # non directional gauges
                try:
                    wave_coords = gp.FRFcoord(self.ncfile['longitude'][:], self.ncfile['latitude'][:])
                except IndexError:
                    wave_coords = gp.FRFcoord(self.ncfile['lon'][:], self.ncfile['lat'][:])
                # try:   # try new variable names
                wavespec = {'time': self.snaptime,                # note this is new variable names??
                            'epochtime': nc.date2num(self.snaptime, self.ncfile['time'].units),
                            'name': str(self.ncfile.title),
                            'wavefreqbin': self.ncfile['waveFrequency'][:],
                            'xFRF': wave_coords['xFRF'],
                            'yFRF': wave_coords['yFRF'],
                            'lat': self.ncfile['latitude'][:],
                            'lon': self.ncfile['longitude'][:],
                            'depth': depth,
                            'Hs': self.ncfile['waveHs'][self.wavedataindex],}
                try:
                    wavespec['peakf'] =  1/self.ncfile['waveTp'][self.wavedataindex]
                except:
                    wavespec['peakf'] = 1/self.ncfile['waveTpPeak'][self.wavedataindex]
                    # except IndexError:
                #     wavespec = {'time': self.snaptime,                # note this is old Variable names remove soon
                #         'epochtime': nc.date2num(self.snaptime, self.ncfile['time'].units),
                #         'name': str(self.ncfile.title),
                #         'wavefreqbin': self.ncfile['waveFrequency'][:],
                #         'xFRF': wave_coords['xFRF'],
                #         'yFRF': wave_coords['yFRF'],
                #         'lat': self.ncfile['lat'][:],
                #         'lon': self.ncfile['lon'][:],
                #         'depth': depth,
                #         'Hs': self.ncfile['waveHs'][self.wavedataindex],
                #         'peakf': self.ncfile['wavePeakFrequency'][self.wavedataindex]}
                try:  # pull time specific data based on self.wavedataindex
                    wavespec['wavedirbin'] = self.ncfile['waveDirectionBins'][:]
                    wavespec['waveDp'] = self.ncfile['wavePeakDirectionPeakFrequency'][self.wavedataindex]
                    wavespec['fspec'] = self.ncfile['waveEnergyDensity'][self.wavedataindex, :]
                    wavespec['waveDm'] = self.ncfile['waveMeanDirection'][self.wavedataindex]
                    wavespec['qcFlagE'] = self.ncfile['qcFlagE'][self.wavedataindex]
                    wavespec['qcFlagD'] = self.ncfile['qcFlagD'][self.wavedataindex]
                    wavespec['a1'] = self.ncfile['waveA1Value'][self.wavedataindex, :]
                    wavespec['a2'] = self.ncfile['waveA2Value'][self.wavedataindex, :]
                    wavespec['b1'] = self.ncfile['waveB1Value'][self.wavedataindex, :]
                    wavespec['b2'] = self.ncfile['waveB2Value'][self.wavedataindex, :]
                    wavespec['dWED'] = self.ncfile['directionalWaveEnergyDensity'][self.wavedataindex, :, :]
                    if wavespec['dWED'].ndim < 3:
                        wavespec['dWED'] = np.expand_dims(wavespec['dWED'], axis=0)
                        wavespec['fspec'] = np.expand_dims(wavespec['fspec'], axis=0)

                except IndexError:
                    # this should throw when gauge is non directional
                    wavespec['wavedirbin'] = np.arange(0, 360, 90)  # 90 degree bins
                    wavespec['waveDp'] = np.zeros(np.size(self.wavedataindex)) * -999
                    wavespec['fspec'] = self.ncfile['waveEnergyDensity'][self.wavedataindex, :]
                    if wavespec['fspec'].ndim < 2:
                        wavespec['fspec'] = np.expand_dims(wavespec['fspec'], axis=0)
                    # multiply the freq spectra for all directions
                    wavespec['dWED'] = np.ones(
                         [np.size(self.wavedataindex), np.size(wavespec['wavefreqbin']), np.size(wavespec['wavedirbin'])])  # *
                    wavespec['dWED'] = wavespec['dWED'] * wavespec['fspec'][:, :, np.newaxis]/len(wavespec['wavedirbin'])
                    wavespec['qcFlagE'] = self.ncfile['qcFlagE'][self.wavedataindex]

                return wavespec

        except (RuntimeError, AssertionError):
            print '     ---- Problem Retrieving wave data from %s\n    - in this time period start: %s  End: %s' % (
            gname, self.d1, self.d2)
            try:
                wavespec = {'lat': self.ncfile['latitude'][:],
                            'lon': self.ncfile['longitude'][:],
                            'name': str(self.ncfile.title),}
            except:
                wavespec = {'lat': self.ncfile['lat'][:],
                            'lon': self.ncfile['lon'][:],
                            'name': str(self.ncfile.title),}
            return wavespec

    def getCurrents(self, gaugenumber=5, roundto=1):
        """
        This function pulls down the currents data from the Thredds Server

            :param gaugenumber:
            gaugenumber = 2, 'awac-11m'
            gaugenumber = 3, awac3 - 8m
            gaugenumber = 4, awac2 - 6m
            gaugenumber = 5, awac1 - 4.5m
            gaugenumber = 6, adopp2 - 3m
            gaugenumber = 7, adopp1 - 2m
            gaugenumber = 13, awac - 5m
            
            :param roundto:
                the time over which the wind record exists
                ie data is collected in 10 minute increments
                data is rounded to the nearst [roundto] (default 1 min)
        """
        assert gaugenumber in [2, 3, 4, 5, 6, 'awac-11m', 'awac-8m', 'awac-6m', 'awac-4.5m', 'adop-3.5m'], 'Input string/number is not a valid gage name/number'

        if gaugenumber == 2 or gaugenumber == 'awac-11m':
            gname = 'AWAC04 - 11m'
            self.dataloc = 'oceanography/currents/awac-11m/awac-11m.ncml'
        elif gaugenumber == 3 or gaugenumber == 'awac-8m':
            gname = 'AWAC 8m'
            self.dataloc = 'oceanography/currents/awac-8m/awac-8m.ncml'
        elif gaugenumber == 4 or gaugenumber == 'awac-6m':
            gname = 'AWAC 6m'
            self.dataloc = 'oceanography/currents/awac-6m/awac-6m.ncml'
        elif gaugenumber == 5 or gaugenumber == 'awac-4.5m':
            gname = 'AWAC 4.5m'
            self.dataloc = 'oceanography/currents/awac-4.5m/awac-4.5m.ncml'
        elif gaugenumber == 6 or gaugenumber == 'adop-3.5m':
            gname = 'Aquadopp 3.5m'
            self.dataloc = 'oceanography/currents/adop-3.5m/adop-3.5m.ncml'

        currdataindex = self.gettime(dtRound=roundto * 60)
        # _______________________________________
        # get the actual current data
        if np.size(currdataindex) > 1:
            curr_aveU = self.ncfile['aveU'][currdataindex]  # pulling depth averaged Eastward current
            curr_aveV = self.ncfile['aveV'][currdataindex]  # pulling depth averaged Northward current
            curr_spd = self.ncfile['currentSpeed'][currdataindex]  # currents speed [m/s]
            curr_dir = self.ncfile['currentDirection'][currdataindex]  # current from direction [deg]
            self.curr_time = nc.num2date(self.ncfile['time'][currdataindex], self.ncfile['time'].units,
                                    self.ncfile['time'].calendar)
            for num in range(0, len(self.curr_time)):
                self.curr_time[num] = self.roundtime(self.curr_time[num], roundto=roundto * 60)

            curr_coords = gp.FRFcoord(self.ncfile['lon'][0], self.ncfile['lat'][0])

            self.curpacket = {
                'name': str(self.ncfile.title),
                'time': self.curr_time,
                'aveU': curr_aveU,
                'aveV': curr_aveV,
                'speed': curr_spd,
                'dir': curr_dir,
                'lat': self.ncfile['lat'][0],
                'lon': self.ncfile['lon'][0],
                'FRF_X': curr_coords['xFRF'],
                'FRF_Y': curr_coords['yFRF'],
                'depth': self.ncfile['depth'][:],
                # Depth is calculated by: depth = -xducerD + blank + (binSize/2) + (numBins * binSize)
                'meanP': self.ncfile['meanPressure'][currdataindex]}

            return self.curpacket

        else:

            print 'ERROR: There is no current data for this time period!!!'
            self.curpacket = None
            return self.curpacket

    def getWind(self, collectionlength=10, gaugenumber=0):
        """
        this function retrieves the wind data from the FDIF server
        collection length is the time over which the wind record exists
            ie data is collected in 10 minute increments
            data is rounded to the nearst [collectionlength] (default 10 min)

            gauge 0 = 932
        """
        # Making gauges flexible
        # different Gauges
        if gaugenumber in ['derived', 'Derived', 0]:
            self.dataloc = u'meteorology/wind/derived/derived.ncml'  # 932 wind gauge
            gname = 'Derived wind gauge '
        elif gaugenumber == 1:
            self.dataloc = u'meteorology/wind/D932/D932.ncml'  # 932 wind gauge
            gname = '932 wind gauge'
        elif gaugenumber == 2:
            gname = '832 wind gauge'
            self.dataloc = u'meteorology/wind/D832/D832.ncml'
        elif gaugenumber == 3:
            gname = '632 wind gauge'
            self.dataloc = u'meteorology/wind/D732/D732.ncml'
        else:
            self.winddataindex = []
            print '<EE>ERROR Specifiy proper Gauge number'

        self.winddataindex = self.gettime(dtRound=collectionlength * 60)
        # remove nan's that shouldn't be there
        # ______________________________________
        if np.size(self.winddataindex) > 0 and self.winddataindex is not None:
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

            self.windtime = nc.num2date(self.allEpoch[self.winddataindex], self.ncfile['time'].units)

            # correcting for wind elevations from Johnson (1999) - Simple Expressions for correcting wind speed data for elevation
            if gaugeht <= 20:
                windspeed_corrected = windspeed * (10 / gaugeht) ** (1 / 7)
            else:
                windspeed_corrected = 'No Corrections done for gauges over 20m, please read: \nJohnson (1999) - Simple Expressions for correcting wind speed data for elevation'
            windpacket = {
                'name': str(self.ncfile.title),  # station name
                'time': self.windtime,  # time
                'vecspeed': windvecspd,  # Vector Averaged Wind Speed
                'windspeed': windspeed,  # Mean Wind Speed
                'windspeed_corrected': windspeed_corrected,  # corrected windspeed
                'winddir': winddir,  # Wind direction from true nort
                'windgust': windgust,  # 5 second largest mean wind speed
                'qcflagS': qcflagS,  # QC flag
                'qcflagD': qcflagD,
                'stdspeed': stdspeed,  # std dev of 10 min wind record
                'minspeed': minspeed,  # min speed in 10 min avg
                'maxspeed': maxspeed,  # max speed in 10 min avg
                'sustspeed': sustspeed,  # 1 min largest mean wind speed
                'lat': self.ncfile['latitude'][:],  # latitude
                'lon': self.ncfile['longitude'][:],  # longitde
                'gaugeht': gaugeht,
            }
            if (windpacket['qcflagD'] == 3).all() or (windpacket['qcflagS'] == 3).all():
                print "Wind querey returned all bad data for speed or direction"
                windpacket = None
            return windpacket
        else:
            print  '     ---- Problem finding wind !!!'
            windpacket = None
            return windpacket

    def getWL(self, collectionlength=6):
        """
        This function retrieves the water level data from the FDIF server
        WL data on server is NAVD88

        collection length is the time over which the wind record exists
            ie data is collected in 10 minute increments
            data is rounded to the nearst [collectionlength] (default 6 min)
        """
        self.dataloc = 'oceanography/waterlevel/eopNoaaTide/eopNoaaTide.ncml'  # this is the back end of the url for waterlevel
        self.WLdataindex = self.gettime(dtRound=collectionlength * 60)

        if np.size(self.WLdataindex) > 1:
            # self.WL = self.ncfile['waterLevel'][self.WLdataindex]
            # self.WLtime = nc.num2date(self.ncfile['time'][self.WLdataindex], self.ncfile['time'].units,
            #                           self.ncfile['time'].calendar)
            # for num in range(0, len(self.WLtime)):
            #     self.WLtime[num] = self.roundtime(self.WLtime[num], roundto=collectionlength * 60)
            self.WLtime = nc.num2date(self.allEpoch[self.WLdataindex], self.ncfile['time'].units)
            self.WLpacket = {
                'name': str(self.ncfile.title),
                'WL': self.ncfile['waterLevel'][self.WLdataindex],
                'time': self.WLtime,
                'lat': self.ncfile['latitude'][:],
                'lon': self.ncfile['longitude'][:],
                'residual': self.ncfile['residualWaterLevel'][self.WLdataindex],
                'predictedWL': self.ncfile['predictedWaterLevel'][self.WLdataindex],
                'gapNum': self.ncfile['gapGauge'][self.WLdataindex],
            }

        else:
            print 'ERROR: there is no WATER level Data for this time period!!!'
            self.WLpacket = None
        return self.WLpacket

    def getBathyFromArcServer(self, output_location, grid_data, method=1):
        """
        This function is designed to pull the raw gridded text file from the Mobile, AL geospatial data server between
        the times of interest (d1, d2) or the most recent file there in
        method = 0 uses the nearest in time to d1
        method = 1 uses the most recent historical survey but not future to d1
        grid_data = must be true or false, true returns gridded data file, false returns transect data
        """
        from getdatatestbed import download_grid_data as DGD
        # url for raw grid data setup on geospatial database
        if grid_data == True:
            service_url = u'http://gis.sam.usace.army.mil/server/rest/services/FRF/FRF/FeatureServer/4'
        elif grid_data == False:
            service_url = u'http://gis.sam.usace.army.mil/server/rest/services/FRF/FRF_DEV2/FeatureServer/4'
        else:
            print 'grid data must be True (returns gridded data) or False (returns transect data)'
        # query the survey to get the file name and the ID of the file name back for the most recent survey on location
        gridID_list, grid_fname_list, grid_date_list = DGD.query_survey_data(service_url, grid_data=grid_data)
        #
        # do logic here for which survey to pull
        #
        mask = (grid_date_list >= self.epochd1) & (grid_date_list < self.epochd2)  # boolean true/false of time
        maskids = np.where(mask)[0]  # where the true values are
        if len(maskids) == 1:  # there is 1 record found between the dates of interest
            print "One bathymetry surveys found between %s and %s" % (self.d1, self.d2)
            gridID = gridID_list[maskids[0]]
            grid_fname = grid_fname_list[maskids[0]]
        elif len(maskids) < 1:
            print "No bathymetry surveys found between %s and %s" % (self.d1, self.d2)
            print "Latest survey found is %s" % sorted(grid_fname_list)[-1]
            if method == 0:
                idx = np.argmin(np.abs(grid_date_list - self.epochd1))  # closest in time
                print 'Bathymetry is taken as closest in TIME - NON-operational'
            # or
            elif method == 1:
                val = (max([n for n in (grid_date_list - self.epochd1) if n < 0]))
                idx = np.where((grid_date_list - self.epochd1) == val)[0]
                if len(idx) > 1:
                    if grid_fname_list[idx[0]] == grid_fname_list[idx[-1]]:
                        idx = idx[0]
                    else:
                        print 'Multiple grids are returned on the Bathy Server, they are not the same, this will cause an error'
                print 'Bathymetry is taken as closest in HISTORY - operational'

            grid_fname = grid_fname_list[int(idx)]
            gridID = gridID_list[int(idx)]
            gridtime = nc.num2date(grid_date_list[int(idx)], 'seconds since 1970-01-01')
            if grid_data == True:
                print "Downloading Bathymetry GRID from %s" % gridtime
            elif grid_data == False:
                print "Downloading SURVEY TRANSECT from %s" % gridtime
            print "This survey is %s  old" % ((self.d1 - gridtime))
        else:
            print ' There Are Multiple Surveys between %s and %s\nPlease Break Simulation up into Multiple Parts.' % (
            self.d1, self.d2)
            print 'The latest survey is %s' % grid_fname_list[maskids[0]]
            raise
        #
        # download the file name and the ID
        #
        DGD.download_survey(gridID, grid_fname, output_location)  # , grid_data)
        return grid_fname  # file name returned w/o prefix simply the name

    def getBathyTransectFromNC(self, profilenumbers=None, method=1, timewindow=None):
        """
        This function gets the bathymetric data from the thredds server, currently designed for the bathy duck experiment
        method == 1  - > 'Bathymetry is taken as closest in HISTORY - operational'
        method == 0  - > 'Bathymetry is taken as closest in TIME - NON-operational'
        :param
        :return:

        """
        # do check here on profile numbers
        # acceptableProfileNumbers = [None, ]
        self.dataloc = u'geomorphology/elevationTransects/survey/surveyTransects.ncml'  # location of the gridded surveys

        try:
            self.bathydataindex = self.gettime()
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

            try:
                # switch back to the FRF ncfile?
                self.ncfile = nc.Dataset(self.FRFdataloc + self.dataloc)
            except:
                pass

            # there's no exact bathy match so find the max negative number where the negitive
            # numbers are historical and the max would be the closest historical
            val = (max([n for n in (self.ncfile['time'][:] - self.epochd1) if n < 0]))
            idx = np.where((self.ncfile['time'][:] - self.epochd1) == val)[0][0]
            print 'Bathymetry is taken as closest in HISTORY - operational'
        elif len(self.bathydataindex) < 1 and method == 0:

            try:
                # switch back to the FRF ncfile?
                self.ncfile = nc.Dataset(self.FRFdataloc + self.dataloc)
            except:
                pass

            idx = np.argmin(np.abs(self.ncfile['time'][:] - self.d1))  # closest in time
            print 'Bathymetry is taken as closest in TIME - NON-operational'
        elif len(self.bathydataindex) > 1:

            try:
                # switch back to the FRF ncfile?
                self.ncfile = nc.Dataset(self.FRFdataloc + self.dataloc)
            except:
                pass

            # DLY Note - this section of the script does NOT work
            # (i.e., if you DO have a survey during your date range!!!)
            timeunits = 'seconds since 1970-01-01 00:00:00'
            d1Epoch = nc.date2num(self.d1, timeunits)
            val = (max([n for n in (self.ncfile['time'][:] - d1Epoch) if n < 0]))
            idx = np.where((self.ncfile['time'][:] - d1Epoch) == val)[0][0]

        # try:
        #     assert profilenumbers in acceptableProfileNumbers, 'Ch3eck numbers should be in %s' % acceptableProfileNumbers
        #     self.bathydataindex = self.gettime(profilenumbers)  # getting the index of the grid
        # except IOError:
                #     self.bathydataindex = []

        # returning whole survey
        idxSingle = idx
        idx = np.argwhere(self.ncfile['surveyNumber'][:] == self.ncfile['surveyNumber'][idxSingle]).squeeze()

        if profilenumbers != None:
            assert pd.Series(profilenumbers).isin(np.unique(self.ncfile['profileNumber'][idx])).all(), 'given profiles don''t Match profiles in database'  # if all of the profile numbers match
            idx2mask = np.in1d(self.ncfile['profileNumber'][idx], profilenumbers)  # boolean true/false of time and profile number
            idx = idx[idx2mask]
        # elif pd.Series(profileNumbers).isin(np.unique(self.ncfile['profileNumber'][:])).any(): #if only some of the profile numbers match
        #     print 'One or more input profile numbers do not match those in the FRF transects!  Fetching data for those that do.'
        #     mask = (self.alltime >= self.d1) & (self.alltime < self.d2) & np.in1d(self.ncfile['profileNumber'][:],profileNumbers)  # boolean true/false of time and profile number


        if np.size(idx) == 0:

            print 'The closest in history to your start date is %s\n' % nc.num2date(self.gridTime[idx],self.ncfile['time'].units)
            print 'Please End new simulation with the date above'
            raise Exception
            idx = self.bathydataindex

        if len(idx) > 0 and idx is not None:

            # now retrieve data with idx
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
            time = nc.num2date(self.ncfile['time'][idx], self.ncfile['time'].units)

            profileDict = {'xFRF': xCoord,
                           'yFRF': yCoord,
                           'elevation': elevation_points,
                           'time': time,
                           'lat': lat,
                           'lon': lon,
                           'northing': northing,
                           'easting': easting,
                           'profileNumber': profileNum,
                           'surveyNumber': surveyNum,
                           'Ellipsoid': Ellipsoid,
                            }
        else:
            profileDict = None

        return profileDict

    def getBathyTransectProfNum(self, method=1):
        """
        This function gets the bathymetric data from the thredds server, currently designed for the bathy duck experiment
        method == 1  - > 'Bathymetry is taken as closest in HISTORY - operational'
        method == 0  - > 'Bathymetry is taken as closest in TIME - NON-operational'
        :param
        :return:

        """
        # do check here on profile numbers
        # acceptableProfileNumbers = [None, ]
        self.dataloc = u'geomorphology/elevationTransects/survey/surveyTransects.ncml'  # location of the gridded surveys

        try:
            self.bathydataindex = self.gettime()
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

            try:
                # switch back to the FRF ncfile?
                self.ncfile = nc.Dataset(self.FRFdataloc + self.dataloc)
            except:
                pass

            # there's no exact bathy match so find the max negative number where the negative
            # numbers are historical and the max would be the closest historical
            val = (max([n for n in (self.ncfile['time'][:] - self.epochd1) if n < 0]))
            idx = np.where((self.ncfile['time'][:] - self.epochd1) == val)[0][0]
            print 'Bathymetry is taken as closest in HISTORY - operational'

        elif len(self.bathydataindex) < 1 and method == 0:

            try:
                # switch back to the FRF ncfile?
                self.ncfile = nc.Dataset(self.FRFdataloc + self.dataloc)
            except:
                pass

            idx = np.argmin(np.abs(self.ncfile['time'][:] - self.d1))  # closest in time
            print 'Bathymetry is taken as closest in TIME - NON-operational'

        elif len(self.bathydataindex) > 1:

            try:
                # switch back to the FRF ncfile?
                self.ncfile = nc.Dataset(self.FRFdataloc + self.dataloc)
            except:
                pass

            # DLY Note - this section of the script does NOT work
            # (i.e., if you DO have a survey during your date range!!!)
            timeunits = 'seconds since 1970-01-01 00:00:00'
            d1Epoch = nc.date2num(self.d1, timeunits)
            val = (max([n for n in (self.ncfile['time'][:] - d1Epoch) if n < 0]))
            idx = np.where((self.ncfile['time'][:] - d1Epoch) == val)[0][0]

        # returning whole survey
        idxSingle = idx
        idx = np.argwhere(self.ncfile['surveyNumber'][:] == self.ncfile['surveyNumber'][idxSingle]).squeeze()

        # what profile numbers are in this survey?
        prof_nums = np.unique(self.ncfile['profileNumber'][idx])

        return prof_nums

    def getBathyGridFromNC(self, method, removeMask=True):
        """
        This function gets the frf krigged grid product, it will currently break with the present link
        bathymetric data from the thredds server, currently designed for the bathy duck experiment
        method == 1  - > 'Bathymetry is taken as closest in HISTORY - operational'
        method == 0  - > 'Bathymetry is taken as closest in TIME - NON-operational'
        :param
        :return:

        """
        self.dataloc = u'survey/gridded/gridded.ncml'  # location of the gridded surveys
        try:
            self.bathydataindex = self.gettime()  # getting the index of the grid
        except IOError:
            self.bathydataindex = []
        if self.bathydataindex != None and len(self.bathydataindex) == 1:
            idx = self.bathydataindex
        elif (self.bathydataindex == None or len(self.bathydataindex) < 1) & method == 1:
            # there's no exact bathy match so find the max negative number where the negitive
            # numbers are historical and the max would be the closest historical
            val = (max([n for n in (self.ncfile['time'][:] - self.epochd1) if n < 0]))
            idx = np.where((self.ncfile['time'][:] - self.epochd1) == val)[0][0]
            print 'Bathymetry is taken as closest in HISTORY - operational'
        elif (self.bathydataindex == None or len(self.bathydataindex) < 1) and method == 0:
            idx = np.argmin(np.abs(self.ncfile['time'][:] - self.d1))  # closest in time
            print 'Bathymetry is taken as closest in TIME - NON-operational'
        elif self.bathydataindex != None and len(self.bathydataindex) > 1:
            val = (max([n for n in (self.ncfile['time'][:] - self.d1) if n < 0]))
            idx = np.where((self.ncfile['time'] - self.d1) == val)[0][0]

            print 'The closest in history to your start date is %s\n' % nc.num2date(self.gridTime[idx], self.ncfile['time'].units)
            print 'Please End new simulation with the date above'
            raise Exception
        # the below line was in place, it should be masking nan's but there is not supposed to be nan's
        # in the data, should only be fill values (-999)
        # elevation_points = np.ma.array(ncfile['elevation'][idx,:,:], mask=np.isnan(ncfile['elevation'][idx,:,:]))
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
        print 'Sim start: %s\nSim End: %s\nSim bathy chosen: %s' % (self.d1, self.d2,
                                                                    nc.num2date(self.ncfile['time'][idx],
                                                                                self.ncfile['time'].units))
        print 'Bathy is %s old' % (self.d2 - nc.num2date(self.ncfile['time'][idx], self.ncfile['time'].units))

        gridDict = {'xCoord': xCoord,
                    'yCoord': yCoord,
                    'elevation': elevation_points,
                    'time': time,
                    'lat': lat,
                    'lon': lon,
                    'northing': northing,
                    'easting': easting
                    }
        return gridDict

    def getBathyDuckLoc(self, gaugenumber):
        """
        this function pulls the stateplane location (if desired) from the surveyed
        FRF coords
        :param gaugenumber:
        :return:
        """
        if type(gaugenumber) != str:
            gaugenumber = str(gaugenumber)
        assert gaugenumber in ['11', '12', '13', '14',  # middle transect
                               '21', '22', '23', '24',  # south transect
                               '83', '84'], 'gauge number not recognized'  # near pier
        try:
            loc = str(self.FRFdataloc + u"projects/bathyduck/data/BathyDuck-ocean_waves_p%s_201510.nc" % gaugenumber)
            ncfile = nc.Dataset(loc)
            xloc = ncfile['xloc'][:]
            yloc = ncfile['yloc'][:]
        except:
            loc = str(self.chlDataLoc + u'projects/bathyduck/data/BathyDuck-ocean_waves_p%s_201510.nc' % gaugenumber)
            ncfile = nc.Dataset(loc)
            xloc = ncfile['xloc'][:]
            yloc = ncfile['yloc'][:]
        assert len(np.unique(xloc)) == 1, "there are different locations in the netCDFfile"
        assert len(np.unique(yloc)) == 1, "There are different Y locations in the NetCDF file"
        locDict = gp.FRFcoord(xloc[0], yloc[0])

        return locDict

    def getWaveGaugeLoc(self, gaugenumber):
        """
        This function gets gauge location data quickly, faster than get wave data

        :param gaugenumber:
        :return:
        """
        if gaugenumber in [0, 'waverider-26m', 'Waverider-26m', '26m']:
            # 26 m wave rider
            self.dataloc = 'oceanography/waves/waverider-26m/waverider-26m.ncml'  # 'oceanography/waves/waverider430/waverider430.ncml'  # 26m buoy
            gname = '26m Waverider Buoy'
        elif gaugenumber == 1 or gaugenumber == 'waverider-17m':
            # 2D 17m waverider
            self.dataloc = 'oceanography/waves/waverider-17m/waverider-17m.ncml'  # 17 m buoy
            gname = '17m Waverider Buoy'
        elif gaugenumber == 2 or gaugenumber == 'awac-11m':
            gname = 'AWAC 11m'
            self.dataloc = 'oceanography/waves/awac-11m/awac-11m.ncml'
        elif gaugenumber == 3 or gaugenumber == 'awac-8m':
            gname = 'AWAC 8m'
            self.dataloc = 'oceanography/waves/awac-8m/awac-8m.ncml'
        elif gaugenumber == 4 or gaugenumber == 'awac-6m':
            gname = 'AWAC 6m'
            self.dataloc = 'oceanography/waves/awac-6m/awac-6m.ncml'
        elif gaugenumber  in [5, 'awac-4.5m', 'awac_4.5m']:
            gname = 'AWAC 4.5m'
            self.dataloc = 'oceanography/waves/awac-4.5m/awac-4.5m.ncml'
        elif gaugenumber == 6 or gaugenumber == 'adop-3.5m':
            gname = 'Aquadopp 3.5m'
            self.dataloc = 'oceanography/waves/adop-3.5m/adop-3.5m.ncml'
        elif gaugenumber == 7 or gaugenumber == 'adop-2m':
            gname = 'Aquadopp01 - 2m'
            self.dataloc = 'oceanography/waves/adop01/adop01.ncml'
        elif gaugenumber == 8 or gaugenumber == 'xp200m':
            gname = 'Paros xp200m'
            self.dataloc = 'oceanography/waves/xp200m/xp200m.ncml'
        elif gaugenumber == 9 or gaugenumber == 'xp150m':
            gname = 'Paros xp150m'
            self.dataloc = 'oceanography/waves/xp150m/xp150m.ncml'
        elif gaugenumber == 10 or gaugenumber == 'xp125m':
            gname = 'Paros xp125m'
            self.dataloc = 'oceanography/waves/xp125m/xp125m.ncml'
        elif gaugenumber == 11 or gaugenumber == 'xp100m':
            gname = 'Paros xp100m'
            self.dataloc = 'oceanography/waves/xp100m/xp100m.ncml'
        elif gaugenumber == 12 or gaugenumber == '8m-array':
            gname = "8m array"
            self.dataloc = 'oceanography/waves/8m-array/8m-array.ncml'
        elif gaugenumber in ['oregonInlet', 'OI', 'oi']:
            gname = 'Oregon Inlet'
            self.dataloc = 'oceanography/waves/waverider-oregon-inlet-nc/waverider-oregon-inlet-nc.ncml'
        else:
            gname = 'There Are no Gauge numbers here'
            raise NameError('Bad Gauge name, specify proper gauge name/number')
        try:
            ncfile = nc.Dataset(self.FRFdataloc + self.dataloc)
        except IOError:
            ncfile = nc.Dataset(self.chlDataLoc + self.dataloc)
        out = {'lat': ncfile['latitude'][:],
               'lon': ncfile['longitude'][:]}
        return out

    def get_sensor_locations_from_thredds(self):

        """ 
        Retrieves lat/lon coordinates for each gauge in gauge_list, converts 
        to state plane and frf coordinates and creates a dictionary containing 
        all three coordinates types with gaugenumbers as keys.

        Returns
        -------
        loc_dict : dict
            Dictionary containing lat/lon, state plane, and frf coordinates
            for each available gaugenumber with gaugenumbers as keys.
        """

        loc_dict = {}

        for g in self.gaugelist:
            loc_dict[g] = {}
            data = loc_dict[g]

            # Get latlon from Thredds server.
            try:
                latlon = self.getWaveGaugeLoc(g)
            except IOError:
                continue

            # lat and lon values currently stored as 1 element arrays. 
            # Cast to float for consistency.
            lat = float(latlon['lat'])
            lon = float(latlon['lon'])

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
        """
        Retrieve sensor coordinate dictionary from file if there is an entry
        within window_days of the specified timestampstr. Otherwise query the 
        Thredds server for location information and update archived data 
        accordingly.

        Parameters
        ----------
        timestamp : datetime.datetime
            timestamp for which coordinates are desired.
        datafile : str
            Name of file containing archived sensor location data.
        window_days : int
            Maximum interval between desired timestamp and closest timestamp
            in datafile to use archived data. If this interval is larger than 
            window_days days, query the Thredds server.

        Returns
        -------
        sensor_locations : dict
            Coordinates in lat/lon, stateplane, and frf for each available 
            gaugenumber (gauges 0 to 12).

        Updates datafile when new information is retrieved.
        """
        try:
            with open(datafile, 'rb') as fid:  # this will close a file when done loading it
                loc_dict = pickle.load(fid)
        except IOError:
            loc_dict = {}
            # now create pickle if one's not around
            # this should be date of searching( or date of gauge placement more acccureately)
            timestamp = self.d1 # DT.strptime(timestampstr, '%Y%m%d_%H%M')
            loc_dict[timestamp] = [self.get_sensor_locations_from_thredds()] #timestamp)
            with open(datafile, 'wb') as fid:
                pickle.dump(loc_dict, fid)

        # loc_dict = pickle.load(open(datafile, 'rb'))
        available_timestamps = np.array(loc_dict.keys()) 
        idx = np.abs(self.d1 - available_timestamps).argmin()
        nearest_timestamp = available_timestamps[idx]
        if abs(self.d1 - nearest_timestamp).days < window_days:
            archived_sensor_locations = loc_dict[nearest_timestamp]
            # MPG: only use locations specified in self.gaugelist (for the case 
            # that there are archived locations that should not be used).
            sensor_locations = collections.OrderedDict()
            for g in self.gaugelist:
                if g in archived_sensor_locations:
                    sensor_locations[g] = archived_sensor_locations[g]
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

    def getBathyGridcBathy(self, **kwargs):
        """
        this functin gets the cbathy data from the below address, assumes fill value of -999

        This function accepts kwargs
            xbound = [xmin, xmax]  which will truncate the cbathy domain to xmin, xmax (frf coord)
            ybound = [ymin, ymax]  which will truncate the cbathy domain to ymin, ymax (frf coord)

        :return:  dictionary with keys:
                'time': time
                'xm':  frf xoordinate x's
                'ym': frf ycoordinates
                'depth': raw cbathy depths
                'depthKF':  kalman filtered hourly depth
                'depthKFError': errors associated with the kalman filter
                'fB':  ?
                'k':  ??
        """
        fillValue = -999  # assumed fill value from the rest of the files taken as less than or equal to

        self.dataloc = u'projects/bathyduck/data/cbathy_old/cbathy.ncml'
        self.cbidx = self.gettime(dtRound=30*60)

        self.cbtime = nc.num2date(self.allEpoch[self.cbidx], 'seconds since 1970-01-01')
        # mask = (time > d1) & (time < d2)
        # assert (emask == mask).all(), 'epoch time is not working'
        # idx = np.where(emask)[0] # this leaves a list that keeps the data iteratable with a size 1.... DON'T CHANGE
        if np.size(self.cbidx) == 1 and self.cbidx == None :
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
            else:# <= used here to handle inclusive initial index inherant in python
                removeMinX = np.argwhere(self.ncfile['xm'][:] <= kwargs['xbounds'][0]).squeeze().max()
            # now max of x
            if (kwargs['xbounds'][1] > self.ncfile['xm'][:]).all():
                removeMaxX = None
            else:
                removeMaxX = np.argwhere(self.ncfile['xm'][:] >= kwargs['xbounds'][1]).squeeze().min() + 1 # python indexing
            xs = slice(removeMinX, removeMaxX)
        else:
            xs = slice(None)
            # cbdata['xm'] = cbdata['xm'][removeMinX:removeMaxX]  # sectioning off data from min to max
            # cbdata['depthKF'] = cbdata['depthKF'][:, :, removeMinX:removeMaxX]
            # cbdata['depthKFError'] = cbdata['depthKFError'][:, :, removeMinX:removeMaxX]
            # cbdata['depth'] = cbdata['depth'][:, :, removeMinX:removeMaxX]
            # cbdata['k'] = cbdata['k'][:, :, removeMinX:removeMaxX, :]
            # cbdata['fB'] = cbdata['fB'][:, :, removeMinX:removeMaxX, :]

        if 'ybounds' in kwargs and np.array(kwargs['ybounds']).size == 2:
            if kwargs['ybounds'][0] > kwargs['ybounds'][1]:
                kwargs['ybounds'] = np.flip(kwargs['ybounds'],axis=0)
            # first min of y
            if (kwargs['ybounds'][0] < self.ncfile['ym'][:]).all():
                # then set the ymin to first index [0]
                removeMinY = 0  # ie get all data
            else:
                removeMinY = np.argwhere(self.ncfile['ym'][:] <= kwargs['ybounds'][0]).squeeze().max()
            ## now max of y
            if (kwargs['ybounds'][1] > self.ncfile['ym'][:]).all():
                removeMaxY = None
            else:
                removeMaxY = np.argwhere(self.ncfile['ym'][:] >= kwargs['ybounds'][1]).squeeze().min()+1  # python indexing
            ys = slice(removeMinY, removeMaxY)
        else:
            ys = slice(None)

            # # <= used here to handle inclusive initial index inherant in python
            # removeMinY = np.argwhere(cbdata['ym'] <= kwargs['ybounds'][0]).squeeze().max()
            # # < used here to handle exclusive ending indexing inherant in python
            # removeMaxY = np.argwhere(cbdata['ym'] > kwargs['ybounds'][1]).squeeze().min()
            # cbdata['ym'] = cbdata['ym'][removeMinY:removeMaxY]
            # cbdata['depthKF'] = cbdata['depthKF'][:, removeMinY:removeMaxY, :]
            # cbdata['depthKFError'] = cbdata['depthKFError'][:, removeMinY:removeMaxY, :]
            # cbdata['depth'] = cbdata['depth'][:, removeMinY:removeMaxY, :]
            # cbdata['k'] = cbdata['k'][:, removeMinY:removeMaxY, :, :]
            # cbdata['fB']= cbdata['fB'][:, removeMinY:removeMaxY, :, :]



        try:
            cbdata = {'time': self.cbtime,  # round the time to the nearest 30 minutes
                      'epochtime': self.allEpoch[self.cbidx],
                      'xm': self.ncfile['xm'][xs],
                      'ym': self.ncfile['ym'][ys],
                      'depth': np.ma.array(self.ncfile['depthfC'][self.cbidx, ys, xs], mask=(self.ncfile['depthfC'][self.cbidx, ys, xs] <= fillValue)), # has different fill value
                      'depthKF': np.ma.array(self.ncfile['depthKF'][self.cbidx, ys, xs], mask=(self.ncfile['depthKF'][self.cbidx, ys, xs] <= fillValue)),
                      'depthKFError': np.ma.array(self.ncfile['depthKF'][self.cbidx, ys, xs], mask=(self.ncfile['depthKF'][self.cbidx, ys, xs] <= fillValue)),
                      'depthfC': np.ma.array(self.ncfile['depthfC'][self.cbidx, ys, xs], mask=(self.ncfile['depthfC'][self.cbidx, ys, xs] <= fillValue)),
                      'depthfCError': np.ma.array(self.ncfile['depthErrorfC'][self.cbidx, ys, xs], mask=(self.ncfile['depthErrorfC'][self.cbidx, ys, xs] <= fillValue)),
                      'fB': np.ma.array(self.ncfile['fB'][self.cbidx, ys, xs, :],  mask=(self.ncfile['fB'][self.cbidx, ys, xs, :] <= fillValue)),
                      'k': np.ma.array(self.ncfile['k'][self.cbidx, ys, xs, :], mask=(self.ncfile['k'][self.cbidx, ys, xs, :] <= fillValue)),
                      'P': np.ma.array(self.ncfile['PKF'][self.cbidx, ys, xs], mask=(self.ncfile['PKF'][self.cbidx, ys, xs] <= fillValue))}  # may need to be masked
            print 'Grabbed cBathy Data, successfully'
        except IndexError:  # there's no data in the Cbathy
            cbdata = None

        return cbdata

    def getLidarRunup(self, removeMasked=True):
        """
        :param: removeMasked will toggle the removing of data points from the tsTime series based on the flag status
        :return:
        """

        self.dataloc = 'oceanography/waves/lidarWaveRunup/lidarWaveRunup.ncml'
        self.lidarIndex = self.gettime(dtRound=60)

        if np.size(self.lidarIndex) > 0 and self.lidarIndex is not None:

            out = {'name': nc.chartostring(self.ncfile[u'station_name'][:]),
                   'lat': self.ncfile[u'lidarLatitude'][:],  # Coordintes
                   'lon': self.ncfile[u'lidarLongitude'][:],
                   'lidarX': self.ncfile[u'lidarX'][:],
                   'lidarY': self.ncfile[u'lidarY'][:],
                   'time': nc.num2date(self.ncfile['time'][self.lidarIndex], self.ncfile['time'].units, self.ncfile['time'].calendar),
                   'totalWaterLevel': self.ncfile['totalWaterLevel'][self.lidarIndex],
                   'elevation': self.ncfile['elevation'][self.lidarIndex, :],
                   'xFRF': self.ncfile[u'xFRF'][self.lidarIndex, :],
                   'yFRF': self.ncfile[u'yFRF'][self.lidarIndex, :],
                   'samplingTime': self.ncfile['tsTime'][:],
                   'runupDownLine': self.ncfile['downLineDistance'][self.lidarIndex, :],
                   'totalWaterLevelQCflag': self.ncfile['totalWaterLevelQCFlag'][self.lidarIndex],
                   'percentMissing': self.ncfile['percentTimeSeriesMissing'][self.lidarIndex],
                   }

            if removeMasked:

                if isinstance(out['elevation'], np.ma.MaskedArray):
                    out['elevation'] = np.array(out['elevation'][~out['elevation'].mask])
                else:
                    pass
                if isinstance(out['totalWaterLevel'], np.ma.MaskedArray):
                    out['totalWaterLevel'] = np.array(out['totalWaterLevel'][~out['totalWaterLevel'].mask])
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
                if isinstance(out['runupDownLine'], np.ma.MaskedArray):
                    out['runupDownLine'] = np.array(out['runupDownLine'][~out['runupDownLine'].mask])
                else:
                    pass
                if isinstance(out['samplingTime'], np.ma.MaskedArray):
                    out['samplingTime'] = np.array(out['samplingTime'][~out['samplingTime'].mask])
                else:
                    pass
            else:
                pass

        else:
            print 'There is no LIDAR data during this time period'
            out = None
        return out

    def getCTD(self):
        
        #THIS FUNCTION IS CURRENTLY BROKEN - THE PROBLEM IS THAT self.ncfile does not have any keys?

        """
        This function gets the CTD data from the thredds server
        :param
        :return:
        """

        # do check here on profile numbers
        # acceptableProfileNumbers = [None, ]
        self.dataloc = u'oceanography/ctd/eop-ctd/eop-ctd.ncml'  # location of the gridded surveys

        try:
            self.ncfile = self.FRFdataloc + self.dataloc
            val = (max([n for n in (self.ncfile['time'][:] - self.epochd1) if n < 0]))
            idx = np.where((self.ncfile['time'][:] - self.epochd1) == val)[0][0]
            print 'CTD data is closest in HISTORY - operational'

        except (RuntimeError, NameError, AssertionError, TypeError):  # if theres any error try to get good data from next location
            try:
                self.ncfile = self.chlDataLoc + self.dataloc
                val = (max([n for n in (self.ncfile['time'][:] - self.epochd1) if n < 0]))
                idx = np.where((self.ncfile['time'][:] - self.epochd1) == val)[0][0]
                print 'CTD data closest in HISTORY - operational'
            except (RuntimeError, NameError, AssertionError, TypeError): # if there are still errors, give up
                idx = []

        if np.size(idx) > 0:
            # now retrieve data with idx
            depth = self.ncfile['depth'][idx]
            lat = self.ncfile['lat'][idx]
            lon = self.ncfile['lon'][idx]
            time = nc.num2date(self.ncfile['time'][idx], self.ncfile['time'].units)
            temp = self.ncfile['waterTemperature'][idx]
            salin = self.ncfile['salinity'][idx]
            soundSpeed = self.ncfile['soundSpeed'][idx]
            sigmaT = self.ncfile['sigmaT'][idx]

            ctd_Dict = {'depth': depth,
                        'temp': temp,
                        'time': time,
                        'lat': lat,
                        'lon': lon,
                        'salin': salin,
                        'soundSpeed': soundSpeed,
                        'sigmaT': sigmaT}
        else:
            ctd_Dict = None

        return ctd_Dict

    def getALT(self, gagename='Alt05', removeMasked=True):

        """
        This function gets the Altimeter data from the thredds server
        :param:
        gagename - 'Alt03, Alt04, Alt05'  This is just the name of the altimeter we want to use
        :return:
        """
        # location of the data
        gage_list = ['Alt03', 'Alt04', 'Alt05']
        assert gagename in gage_list, 'Input string is not a valid gage name'
        if gagename == 'Alt05':
            a = 1
            b = 0
            self.dataloc = u'geomorphology/altimeter/Alt05-altimeter/Alt05-altimeter.ncml'

        elif gagename == 'Alt04':
            self.dataloc = u'geomorphology/altimeter/Alt04-altimeter/Alt04-altimeter.ncml'
        elif gagename == 'Alt03':
            self.dataloc = u'geomorphology/altimeter/Alt03-altimeter/Alt03-altimeter.ncml'

        altdataindex = self.gettime(dtRound=1 * 60)

        # get the actual current data
        if np.size(altdataindex) > 1:

            alt_lat = self.ncfile['Latitude'][0]  # pulling latitude
            alt_lon = self.ncfile['Longitude'][0]  # pulling longitude
            alt_be = self.ncfile['bottomElevation'][altdataindex]  # pulling bottom elevation
            alt_pkf = self.ncfile['PKF'][altdataindex]  # i have no idea what this stands for...
            alt_stationname = self.ncfile['station_name'][0]  # name of the station
            self.alt_timestart = nc.num2date(self.ncfile['timestart'][altdataindex], self.ncfile['timestart'].units, self.ncfile['time'].calendar)
            self.alt_timeend = nc.num2date(self.ncfile['timeend'][altdataindex], self.ncfile['timeend'].units, self.ncfile['time'].calendar)
            self.alt_time = nc.num2date(self.ncfile['time'][altdataindex], self.ncfile['time'].units, self.ncfile['time'].calendar)
            for num in range(0, len(self.alt_time)):
                self.alt_time[num] = self.roundtime(self.alt_time[num], roundto=1 * 60)
                self.alt_timestart[num] = self.roundtime(self.alt_timestart[num], roundto=1 * 60)
                self.alt_timeend[num] = self.roundtime(self.alt_timeend[num], roundto=1 * 60)

            alt_coords = gp.FRFcoord(alt_lon, alt_lat)

            if removeMasked:
                self.altpacket = {
                    'name': str(self.ncfile.title),
                    'time': np.array(self.alt_time[~alt_be.mask]),
                    'lat': alt_lat,
                    'PKF': np.array(alt_pkf[~alt_be.mask]),
                    'lon': alt_lon,
                    'xFRF': alt_coords['xFRF'],
                    'yFRF': alt_coords['yFRF'],
                    'stationName': alt_stationname,
                    'gageName': gagename,
                    'timeStart': np.array(self.alt_timestart[~alt_be.mask]),
                    'timeEnd': np.array(self.alt_timeend[~alt_be.mask]),
                    'bottomElev': np.array(alt_be[~alt_be.mask])}
            else:
                self.altpacket = {
                    'name': str(self.ncfile.title),
                    'time': self.alt_time,
                    'lat': alt_lat,
                    'PKF': alt_pkf,
                    'lon': alt_lon,
                    'xFRF': alt_coords['xFRF'],
                    'yFRF': alt_coords['yFRF'],
                    'stationName': alt_stationname,
                    'gageName': gagename,
                    'timeStart': self.alt_timestart,
                    'timeEnd': self.alt_timeend,
                    'bottomElev': alt_be}

            return self.altpacket
        else:
            print 'No %s data found for this period' %(gagename)
            self.altpacket = None
            return self.altpacket

    def getLidarWaveProf(self, removeMasked=True):

        """
        :param: removeMasked will toggle the removing of data points from the tsTime series based on the flag status
        :return:
        """
        self.dataloc = 'oceanography/waves/lidarHydrodynamics/lidarHydrodynamics.ncml'
        self.lidarIndex = self.gettime(dtRound=60)
        if np.size(self.lidarIndex) > 0 and self.lidarIndex is not None:

            out = {'name': nc.chartostring(self.ncfile[u'station_name'][:]),
                   'lat': self.ncfile[u'lidarLatitude'][:],  # Coordinates
                   'lon': self.ncfile[u'lidarLongitude'][:],
                   'lidarX': self.ncfile[u'lidarX'][:],
                   'lidarY': self.ncfile[u'lidarY'][:],
                   'xFRF': self.ncfile[u'xFRF'][:],
                   'yFRF': self.ncfile[u'yFRF'][:],
                   'runupDownLine': self.ncfile['downLineDistance'][:],
                   'waveFrequency': self.ncfile['waveFrequency'][:],
                   'time': nc.num2date(self.ncfile['time'][self.lidarIndex], self.ncfile['time'].units, self.ncfile['time'].calendar),
                   'hydroQCflag': self.ncfile['hydrodynamicsFlag'][self.lidarIndex],
                   'waterLevel': self.ncfile['waterLevel'][self.lidarIndex, :],
                   'waveHs': self.ncfile['waveHs'][self.lidarIndex, :],
                   'waveHsIG': self.ncfile['waveHsIG'][self.lidarIndex, :],
                   'waveHsTotal': self.ncfile['waveHsTotal'][self.lidarIndex, :],
                   'waveSkewness': self.ncfile['waveSkewness'][self.lidarIndex, :],
                   'waveAsymmetry': self.ncfile['waveAsymmetry'][self.lidarIndex, :],
                   'waveEnergyDensity': self.ncfile['waveEnergyDensity'][self.lidarIndex, :, :],
                   'percentMissing': self.ncfile['percentTimeSeriesMissing'][self.lidarIndex, :],
                   }


            if removeMasked:

                if isinstance(out['waterLevel'], np.ma.MaskedArray):
                    out['waterLevel'] = np.array(out['waterLevel'][~out['waterLevel'].mask])
                else:
                    pass
                if isinstance(out['waveHs'], np.ma.MaskedArray):
                    out['waveHs'] = np.array(out['waveHs'][~out['waveHs'].mask])
                else:
                    pass
                if isinstance(out['waveHsIG'], np.ma.MaskedArray):
                    out['waveHsIG'] = np.array(out['waveHsIG'][~out['waveHsIG'].mask])
                else:
                    pass
                if isinstance(out['waveHsTotal'], np.ma.MaskedArray):
                    out['waveHsTotal'] = np.array(out['waveHsTotal'][~out['waveHsTotal'].mask])
                else:
                    pass
                if isinstance(out['waveSkewness'], np.ma.MaskedArray):
                    out['waveSkewness'] = np.array(out['waveSkewness'][~out['waveSkewness'].mask])
                else:
                    pass
                if isinstance(out['waveAsymmetry'], np.ma.MaskedArray):
                    out['waveAsymmetry'] = np.array(out['waveAsymmetry'][~out['waveAsymmetry'].mask])
                else:
                    pass
                if isinstance(out['waveEnergyDensity'], np.ma.MaskedArray):
                    out['waveEnergyDensity'] = np.array(out['waveEnergyDensity'][~out['waveEnergyDensity'].mask])
                else:
                    pass

            else:
                pass
        else:
            print 'There is no LIDAR data during this time period'
            out = None
        return out

    def getBathyRegionalDEM(self, utmEmin, utmEmax, utmNmin, utmNmax):

        """
        :param utmEmin: left side of DEM bounding box in UTM
        :param utmEmax: right side of DEM bounding box in UTM
        :param utmNmin: bottom of DEM bounding box in UTM
        :param utmNmax: top of DEM bounding box in UTM
        
        :return:
          dictionary comprising a smaller rectangular piece of the DEM data, bounded by inputs above
        """

        self.dataloc = u'grids/RegionalBackgroundDEM/backgroundDEM.nc'
        self.ncfile = nc.Dataset(self.crunchDataLoc + self.dataloc)

        # get a 1D ARRAY of the utmE and utmN of the rectangular grid (NOT the full grid!!!)
        utmE_all = self.ncfile['utmEasting'][0, :]
        utmN_all = self.ncfile['utmNorthing'][:, 0]

        # find indices I need to pull...
        ni_min = np.where(utmE_all >= utmEmin)[0][0]
        ni_max = np.where(utmE_all <= utmEmax)[0][-1]
        nj_min = np.where(utmN_all <= utmNmax)[0][0]
        nj_max = np.where(utmN_all >= utmNmin)[0][-1]

        assert (np.size(ni_min) >= 1) & (np.size(ni_max) >= 1) & (np.size(nj_min) >= 1) & (np.size(nj_max) >= 1), 'getBathyDEM Error: bounding box is too close to edge of DEM domain'

        out = {}
        out['utmEasting'] = self.ncfile['utmEasting'][nj_min:nj_max + 1, ni_min:ni_max+1]
        out['utmNorthing'] = self.ncfile['utmNorthing'][nj_min:nj_max + 1, ni_min:ni_max+1]
        out['latitude'] = self.ncfile['latitude'][nj_min:nj_max + 1, ni_min:ni_max + 1]
        out['longitude'] = self.ncfile['longitude'][nj_min:nj_max + 1, ni_min:ni_max + 1]
        out['bottomElevation'] = self.ncfile['bottomElevation'][nj_min:nj_max + 1, ni_min:ni_max + 1]

        return out

class getDataTestBed:

    def __init__(self, d1, d2):
        """
        Initialization description here
        Data are returned in self.datainex are inclusive at d1,d2
        Data comes from waverider 632 (26m?)
        """

        self.rawdataloc_wave = []
        self.outputdir = []  # location for outputfiles
        self.d1 = d1  # start date for data grab
        self.d2 = d2  # end data for data grab
        self.timeunits = 'seconds since 1970-01-01 00:00:00'
        self.epochd1 = nc.date2num(self.d1, self.timeunits)
        self.epochd2 = nc.date2num(self.d2, self.timeunits)
        self.comp_time()
        self.FRFdataloc = u'http://134.164.129.55/thredds/dodsC/FRF/'
        self.crunchDataLoc = u'http://134.164.129.55/thredds/dodsC/cmtb/'
        self.chlDataLoc = u'https://chlthredds.erdc.dren.mil/thredds/dodsC/cmtb/' #'http://10.200.23.50/thredds/dodsC/frf/'
        assert type(self.d2) == DT.datetime, 'd1 need to be in python "Datetime" data types'
        assert type(self.d1) == DT.datetime, 'd2 need to be in python "Datetime" data types'

    def comp_time(self):
        """
        Test if times are backwards
        """
        assert self.d2 >= self.d1, 'finish time: d2 needs to be after start time: d1'

    def roundtime(self, dt=None, roundto=60):
        """
        Round a datetime object to any time laps in seconds
        Author: Thierry Husson 2012 - Use it as you want but don't blame me.
        :rtype: object

        :param dt:
            datetime.datetime object, default now.
        :param roundto:
            Closest number of SECONDS to round to, default 1 minute
        """
        if dt is None:
            dt = DT.datetime.now()
        seconds = (dt - dt.min).seconds
        # // is a floor division, not a comment on following line:
        rounding = (seconds + roundto / 2) // roundto * roundto
        return dt + DT.timedelta(0, rounding - seconds, -dt.microsecond)

    def gettime(self, dtRound=60):
        """
        this function opens the netcdf file, pulls down all of the time, then pulls the dates of interest
        from the THREDDS (data loc) server based on d1,d2, and data location
        it returns the indicies in the NCML file of the dates d1>=time>d2
        INPUTS:

             :param dtRound: the time delta of the data out of interest, default minute (60 second)

        """
        # TODO find a way to pull only hourly data or regular interval of desired time
        # todo this use date2index and create a list of dates see help(nc.date2index)
        try:

            self.ncfile = nc.Dataset(self.crunchDataLoc + self.dataloc) #loads all of the netCDF file
            #            try:

            self.allEpoch = sb.myround(self.ncfile['time'][:], base=dtRound) # round to nearest minute
            # now find the boolean!
            mask = (self.allEpoch >= self.epochd1) & (self.allEpoch < self.epochd2)
            idx = np.argwhere(mask).squeeze()
            # old slow way of doing time!
            # self.alltime = nc.num2date(self.ncfile['time'][:], self.ncfile['time'].units,
            #                            self.ncfile['time'].calendar) # converts all epoch time to datetime objects
            # for i, date in enumerate(self.alltime):  # rounds time to nearest
            #     self.alltime[i] = self.roundtime(dt=date, roundto=dtRound)
            #
            # mask = (self.alltime >= self.d1) & (self.alltime < self.d2)  # boolean true/false of time
            # if (np.argwhere(mask).squeeze() == idx).all():
            #     print '.... old Times match New Times' % np.argwhere(mask).squeeze()
            assert np.size(idx) > 0, 'no data locally, check CHLthredds'
            print "Data Gathered From Local Thredds Server"

        except (IOError, RuntimeError, NameError, AssertionError):  # if theres any error try to get good data from next location
            try:
                # MPG: Use self.chlDataLoc with 'frf/' removed from string for correct url.
                self.ncfile = nc.Dataset(self.chlDataLoc.replace('frf/', 'cmtb/') + self.dataloc)
                self.allEpoch = sb.myround(self.ncfile['time'][:], base=dtRound) # round to nearest minute
                # now find the boolean !
                emask = (self.allEpoch >= self.epochd1) & (self.allEpoch < self.epochd2)
                idx = np.argwhere(emask).squeeze()

                # self.alltime = nc.num2date(self.ncfile['time'][:], self.ncfile['time'].units,
                #                            self.ncfile['time'].calendar)
                # for i, date in enumerate(self.alltime):
                #     self.alltime[i] = self.roundtime(dt=date, roundto=dtRound)
                # # mask = (sb.roundtime(self.ncfile['time'][:]) >= self.epochd1) & (sb.roundtime(self.ncfile['time'][:]) < self.epochd2)\
                #
                # mask = (self.alltime >= self.d1) & (self.alltime < self.d2)  # boolean true/false of time
                #
                # idx = np.argwhere(mask).squeeze()


                try:
                    assert np.size(idx) > 0, ' There are no data within the search parameters for this gauge'
                    print "Data Gathered from CHL thredds Server"
                except AssertionError:
                    idx = None
            except IOError:  # this occors when thredds is down
                print ' Trouble Connecteing to data on CHL Thredds'
                idx = None

        self.ncfile = nc.Dataset(self.crunchDataLoc + self.dataloc)
        # switch us back to the local THREDDS if it moved us to CHL

        return idx

    def getGridCMS(self, method):
        """

        :param method: can be [1, historical, history]  for historical
                     can be [0, 'time'] for non oporational consideration
        :return:
        """
        self.dataloc = 'grids/CMSwave_v1/CMSwave_v1.ncml'
        try:
            self.bathydataindex = self.gettime()  # getting the index of the grid
        except IOError:
            self.bathydataindex = None
        if self.bathydataindex != None and np.size(self.bathydataindex) == 1:
            idx = self.bathydataindex.squeeze()
        elif (self.bathydataindex == None or len(self.bathydataindex) < 1) and method in [1, 'historical', 'history']:
            # there's no exact bathy match so find the max negative number where the negitive
            # numbers are historical and the max would be the closest historical
            val = (max([n for n in (self.ncfile['time'][:] - self.epochd1) if n < 0]))
            idx = np.where((self.ncfile['time'][:] - self.epochd1) == val)[0][0]
            print 'Bathymetry is taken as closest in HISTORY - operational'
        elif (self.bathydataindex == None or np.size(self.bathydataindex) < 1) and method == 0:
            idx = np.argmin(np.abs(self.ncfile['time'][:] - self.d1))  # closest in time
            print 'Bathymetry is taken as closest in TIME - NON-operational'
        elif self.bathydataindex != None and len(self.bathydataindex) > 1:
            val = (max([n for n in (self.ncfile['time'][:] - self.d1) if n < 0]))
            idx = np.where((self.ncfile['time'] - self.d1) == val)[0][0]

        #
        # if self.bathydataindex is not None and  len(self.bathydataindex) == 1:
        #     idx = self.bathydataindex
        # elif self.bathydataindex is not None and len(self.bathydataindex) < 1 and method in [1, 'historical', 'history']:
        #     # there's no exact bathy match so find the max negative number where the negitive
        #     # numbers are historical and the max would be the closest historical
        #     val = (max([n for n in (self.ncfile['time'][:] - self.epochd1) if n < 0]))
        #     idx = np.where((self.ncfile['time'][:] - self.epochd1) == val)[0][0]
        #     print 'Bathymetry is taken as closest in HISTORY - operational'
        # elif self.bathydataindex is not None and len(self.bathydataindex) < 1 and method == 0:
        #     idx = np.argmin(np.abs(self.ncfile['time'][:] - self.d1))  # closest in time
        #     print 'Bathymetry is taken as closest in TIME - NON-operational'
        # elif self.bathydataindex is not None and len(self.bathydataindex) > 1:
        #     val = (max([n for n in (self.ncfile['time'][:] - self.d1) if n < 0]))
        #     idx = np.where((self.ncfile['time'] - self.d1) == val)[0][0]


            print 'The closest in history to your start date is %s\n' % nc.num2date(self.gridTime[idx], self.ncfile['time'].units)
            print 'Please End new simulation with the date above'

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

            time = nc.num2date(self.ncfile['time'][idx], self.ncfile['time'].units)

            gridDict = {'xCoord': xCoord,
                        'yCoord': yCoord,
                        'elevation': elevation_points,
                        'time': time,
                        'lat': lat,
                        'lon': lon,
                        'northing': northing,
                        'easting': easting,
                        'x0': self.ncfile['x0'][:],
                        'azimuth': self.ncfile['azimuth'][:],
                        'y0': self.ncfile['y0'][:],
                        }
            return gridDict

    def getBathyIntegratedTransect(self, method=1, ForcedSurveyDate=None):
        """
        This function gets the integraated bathy, useing the plant (2009) method.
        :param method: method == 1  - > 'Bathymetry is taken as closest in HISTORY - operational'
                       method == 0  -  > 'Bathymetry is taken as closest in TIME - NON-operational'
        :param ForcedSurveyDate:  This is to force a date of survey gathering


        :return:
        """
        if ForcedSurveyDate != None:
            # d1 is used in the gettime function,
            # to force a selection of survey date self.d1/d2 is changed to the forced
            # survey date and then changed back using logged start/stop
            # a check is in place to ensure that the retieved time is == to the forced time
            oldD1 = self.d1
            oldD2 = self.d2
            oldD1epoch = self.epochd1
            oldD2epoch = self.epochd2
            self.d1 = ForcedSurveyDate  # change time one
            self.d2 = ForcedSurveyDate + DT.timedelta(0,1)  # and time 2
            self.epochd1 = nc.date2num(self.d1, 'seconds since 1970-01-01')
            self.epochd2 = nc.date2num(self.d2, 'seconds since 1970-01-01')
            print '!!!Forced bathy date %s' % ForcedSurveyDate

        self.dataloc = 'integratedBathyProduct/survey/survey.ncml'
        try:
            self.bathydataindex = self.gettime()  # getting the index of the grid
        except IOError:
            self.bathydataindex = []  # when a server is not available

        if self.bathydataindex != None and np.size(self.bathydataindex) == 1:
            idx = self.bathydataindex.squeeze()
        elif (self.bathydataindex == None or len(self.bathydataindex) < 1) & method == 1:
            # there's no exact bathy match so find the max negative number where the negitive
            # numbers are historical and the max would be the closest historical
            val = (max([n for n in (self.ncfile['time'][:] - self.epochd1) if n < 0]))
            idx = np.where((self.ncfile['time'][:] - self.epochd1) == val)[0][0]
            print 'Bathymetry is taken as closest in HISTORY - operational'
        elif (self.bathydataindex == None or np.size(self.bathydataindex) < 1) and method == 0:
            idx = np.argmin(np.abs(self.ncfile['time'][:] - self.epochd1))  # closest in time
            print 'Bathymetry is taken as closest in TIME - NON-operational'
        elif self.bathydataindex != None and len(self.bathydataindex) > 1:
            val = (max([n for n in (self.ncfile['time'][:] - self.epochd1) if n < 0]))
            idx = np.where((self.ncfile['time'][:] - self.epochd1) == val)[0][0]

            print 'The closest in history to your start date is %s\n' % nc.num2date(self.ncfile['time'][idx],
                                                                                    self.ncfile['time'].units)
            print 'Please End new simulation with the date above'
            raise Exception

        # the below line was in place, it should be masking nan's but there is not supposed to be nan's
        # in the data, should only be fill values (-999)
        # elevation_points = np.ma.array(ncfile['elevation'][idx,:,:], mask=np.isnan(ncfile['elevation'][idx,:,:]))
        # remove -999's
        elevation_points = self.ncfile['elevation'][idx, :, :]
        xCoord = self.ncfile['xFRF'][:]
        yCoord = self.ncfile['yFRF'][:]
        lat = self.ncfile['latitude'][:]
        lon = self.ncfile['longitude'][:]
        try:
            northing = self.ncfile['northing'][:]
            easting = self.ncfile['easting'][:]
        except IndexError:
            northing = None
            easting = None


        # putting dates and times back for all the other instances that use get time
        if ForcedSurveyDate != None:
            self.d1 = oldD1
            self.d2 = oldD2
            self.epochd2 = oldD2epoch
            self.epochd1 = oldD1epoch

        # this commented out section will work once the times on the CHL THREDDS are fixed.
        # Until then, it will error because self.allEpoch was obtained
        # from the CHL THREDDS times, and they are screwed all up!
        # bathyT = nc.num2date(self.allEpoch[idx], 'seconds since 1970-01-01')
        bathyT = nc.num2date(self.ncfile['time'][idx], 'seconds since 1970-01-01')

        print '  Measured Bathy is %s old' % (self.d2 - bathyT)

        gridDict = {'xFRF': xCoord,
                    'yFRF': yCoord,
                    'elevation': elevation_points,
                    'time': bathyT,
                    'lat': lat,
                    'lon': lon,
                    'northing': northing,
                    'easting': easting,
                    'surveyNumber': self.ncfile['surveyNumber'][idx]
                    }
        return gridDict

    def getStwaveField(self, var, prefix, local=True, ijLoc=None):
        """

        :param Local: defines whether the data is from the nested simulation or the regional simulation
        :param ijLoc:  x or y or (x,y) tuple of location of interest in FRF Coordinates
        :return:
        """
        if local == True:
            grid = 'Local'
        elif local == False:
            grid = 'Regional'
        ############## setting up the ncfile ############################
        if prefix == 'CBHPStatic' and local == True:  # this is needed because projects are stored in weird place
            ncfile = nc.Dataset('http://crunchy:8080/thredds/dodsC/CMTB/projects/bathyDuck_SingleBathy_CBHP/Local_Field/Local_Field.ncml')
        elif prefix =='CBHPStatic' and local == False:
            ncfile = nc.Dataset('http://crunchy:8080/thredds/dodsC/CMTB/projects/bathyDuck_SingleBathy_CBHP/Regional_Field/Regional_Field.ncml')
        else:  # this is standard operational model url Structure
            ncfile = nc.Dataset(self.crunchDataLoc + u'waveModels/STWAVE/%s/%s-Field/%s-Field.ncml' % (prefix, grid, grid))
        assert var in ncfile.variables.keys(), 'variable called is not in file please use\n%s' % ncfile.variables.keys()
        mask = (ncfile['time'][:] >= nc.date2num(self.d1, ncfile['time'].units)) & (
            ncfile['time'][:] <= nc.date2num(self.d2, ncfile['time'].units))
        idx = np.where(mask)[0]
        assert np.size(idx > 0), " there's no data"
        print 'getting %s STWAVE  %s %s Data' % (prefix, grid, var)
        # now creating tool to remove single data point
        if ijLoc != None:
            assert len(ijLoc) == 2, 'if giving a postion, must be a tuple of i, j location (of length 2)'
            if type(ijLoc[0]) == int:
                x = ncfile[var].shape[1] - ijLoc[0]  # the data are stored with inverse indicies to grid node locations
                y = ijLoc[1]  # use location given by function call
            else: # ijLoc[0] == slice:
                x = ijLoc[0]
                y = np.argmin(np.abs(ncfile['yFRF'][:] - ijLoc[1]))

        else:
            x = slice(None)  # take entire data
            y = slice(None)  # take entire data

        if ncfile[var][idx].ndim > 2 and ncfile[var][idx].shape[0] > 100: # looping through ... if necessicary
            list = np.round(np.linspace(idx.min(), idx.max(), idx.max()-idx.min(), endpoint=True, dtype=int))
            # if idx.max() not in list:
            #     list = np.append(list, idx.max())
            xFRF = ncfile['xFRF'][x]
            yFRF = ncfile['yFRF'][y]
            if len(list) < 100:
                    bathy = np.array(ncfile[var][np.squeeze(list)])
                    time = nc.num2date(ncfile['time'][np.squeeze(list)], ncfile['time'].units)
            else:
                for num, minidx in enumerate(list):
                    if len(list) < 100:
                        bathy = np.array(ncfile[var][np.squeeze(list)])
                        time = nc.num2date(ncfile['time'][np.squeeze(list)], ncfile['time'].units)
                    elif num == 0:
                        bathy = np.array(ncfile[var][range(minidx, list[num + 1]), y, x])
                        time = nc.num2date(np.array(ncfile['time'][range(minidx, list[num + 1])]), ncfile['time'].units)
                    elif minidx == list[-1]:
                        lastIdx = (idx - minidx)[(idx - minidx) >= 0] + minidx
                        bathy = np.append(bathy, ncfile[var][lastIdx, y, x], axis=0)
                        time = np.append(time, nc.num2date(ncfile['time'][lastIdx], ncfile['time'].units), axis=0)
                    else:
                        bathy = np.append(bathy, ncfile[var][range(minidx, list[num + 1]), y, x], axis=0)
                        time = np.append(time,
                                         nc.num2date(ncfile['time'][range(minidx, list[num + 1])], ncfile['time'].units),
                                         axis=0)
        else:
            bathy = ncfile[var][idx, y, x]
            xFRF = ncfile['xFRF'][x]
            yFRF = ncfile['yFRF'][y]
            time = nc.num2date(ncfile['time'][np.squeeze(idx)], ncfile['time'].units)
        # package for output
        field = {'time': time,
                 'epochtime': ncfile['time'][idx], # pulling down epoch time of interest
                 var: bathy,
                 'xFRF': xFRF,
                 'yFRF': yFRF,
                 }
        try:
            field['bathymetryDate'] = ncfile['bathymetryDate'][idx]
        except IndexError:
            field['bathymetryDate'] = np.ones_like(field['time'])

        assert field[var].shape[0] == len(field['time']), " the indexing is wrong for pulling down bathy"
        return field

    def getWaveSpecSTWAVE(self, prefix, gaugenumber, local=True):
            """
            This function pulls down the data from the thredds server and puts the data into proper places
            to be read for STwave Scripts
            this will return the wavespec with dir/freq bin and directional wave energy

            :param gaugenumber:
                gaugenumber = 0, 26m wave rider
                gaugenumber = 1, 17m waverider
                gaugenumber = 2, awac4 - 11m
                gaugenumber = 3, awac3 - 8m
                gaugenumber = 4, awac2 - 6m
                gaugenumber = 5, awac1 - 5m
                gaugenumber = 6, adopp2 - 3m
                gaugenumber = 7, adopp1 - 2m
                gaugenumber = 8,  Paros xp200m
                gaugenumber = 9,  Paros xp150m
                gaugenumber = 10, Paros xp125m
                gaugenumber = 11, Paros xp100m
                gaugenumber = 12, 8 m array
            :param collectionlength:
                s the time over which the wind record exists
                ie data is collected in 10 minute increments time is rounded
                to nearest 10min increment
                data is rounded to the nearst [collectionlength] (default 30 min)
            """
            # Making gauges flexible
            if prefix in ['CB', 'HP', 'CBHP', 'FP']:
                model = 'STWAVE'
                urlFront = 'waveModels/%s/%s' %(model, prefix)
            elif prefix.startswith('S') and prefix[1].isdigit():  # this is static bathy
                model = 'STWAVE'
                urlFront = 'projects/%s/CBHP/SingleBathy_%s' %(model, prefix[1:])
            elif prefix in ['CBThresh_0']:
                model = 'STWAVE'
                urlFront = 'projects/STWAVE/CBThresh_0'
            ############### now identify file name #################
            if gaugenumber in [0, 'waverider-26m', 'Waverider-26m', '26m']:
                # 26 m wave rider
                fname = 'waverider-26m/waverider-26m.ncml'
                gname = '26m Waverider Buoy'
            elif gaugenumber in [1, 'Waverider-17m', 'waverider-17m']:
                # 2D 17m waverider
                fname = 'waverider-17/waverider-17m.ncml'
                gname = '17m Waverider Buoy'
            elif gaugenumber in [2, 'AWAC-11m', 'awac-11m', 'Awac-11m']:
                gname = 'AWAC04 - 11m'
                fname = 'awac11m/awac11m.ncml'
            elif gaugenumber in [3, 'awac-8m', 'AWAC-8m', 'Awac-8m', 'awac 8m',
                                 '8m-Array', '8m Array', '8m array', '8m-array']:
                gname = 'AWAC 8m'
                fname = 'wac-8m/awac-8m.ncml'
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
            self.dataloc = urlFront +'/'+ fname
            try:
                self.wavedataindex = self.gettime()
                assert np.array(self.wavedataindex).all() != None, 'there''s no data in your time period'

                if np.size(self.wavedataindex) >= 1:
                    wavespec = {'epochtime': self.ncfile['time'][self.wavedataindex],
                                'time': nc.num2date(self.ncfile['time'][self.wavedataindex], self.ncfile['time'].units),
                                'name': nc.chartostring(self.ncfile['station_name'][:]),
                                'wavefreqbin': self.ncfile['waveFrequency'][:],
                                #'lat': self.ncfile['lat'][:],
                                #'lon': self.ncfile['lon'][:],
                                'Hs': self.ncfile['waveHs'][self.wavedataindex],
                                'peakf': self.ncfile['waveTp'][self.wavedataindex],
                                'wavedirbin': self.ncfile['waveDirectionBins'][:],
                                'dWED': self.ncfile['directionalWaveEnergyDensity'][self.wavedataindex, :, :],
                                # 'waveDp': self.ncfile['wavePeakDirectionPeakFrequency'][self.wavedataindex],  # 'waveDp'][self.wavedataindex]
                                'waveDm': self.ncfile['waveDm'][self.wavedataindex],
                                'waveTm': self.ncfile['waveTm'][self.wavedataindex],
                                'waveTp': self.ncfile['waveTp'][self.wavedataindex],
                                'WL': self.ncfile['waterLevel'][self.wavedataindex],
                                'Umag': self.ncfile['Umag'][self.wavedataindex],
                                'Udir': self.ncfile['Udir'][self.wavedataindex],
                                'fspec': self.ncfile['directionalWaveEnergyDensity'][self.wavedataindex, :,:].sum(axis=2),
                                'qcFlag': self.ncfile['qcFlag'][self.wavedataindex]}

            except (RuntimeError, AssertionError):
                print '<<ERROR>> Retrieving data from %s\n in this time period start: %s  End: %s' % (
                    gname, self.d1, self.d2)
                wavespec = None
            return wavespec

    def getLidarRunup(self, removeMasked=True):
        """

        :return:
        """
        self.dataloc = u'projects/tucker/runup/test.ncml'
        self.lidarIndex = self.gettime(dtRound=60)
        if np.size(self.lidarIndex) > 0 and self.lidarIndex is not None:

            out = {'name': nc.chartostring(self.ncfile[u'station_name'][:]),
                   'lat': self.ncfile[u'lidarLatitude'][:],  # Coordintes
                   'lon': self.ncfile[u'lidarLongitude'][:],
                   'lidarX': self.ncfile[u'lidarX'][:],
                   'lidarY': self.ncfile[u'lidarY'][:],
                   'time': self.ncfile[u'time'][self.lidarIndex],
                   'totalWaterLevel': self.ncfile['totalWaterLevel'][self.lidarIndex],
                   'elevation': self.ncfile['elevation'][self.lidarIndex],

                   'samplingTime': self.ncfile['tsTime'][self.lidarIndex, :], # this will need to be changed once Tucker uploads the new ncml file, will not be dimensioned in time!
                   # 'samplingTime': self.ncfile['tsTime'][:],

                   'frfX': self.ncfile[u'xFRF'][self.lidarIndex],
                   'frfY': self.ncfile[u'yFRF'][self.lidarIndex],
                   'runupDownLine': self.ncfile['downLineDistance'][self.lidarIndex],
                   'totalWaterLevelQCflag': self.ncfile['totalWaterLevelQCFlag'][self.lidarIndex],
                   'percentMissing': self.ncfile['percentTimeSeriesMissing'][self.lidarIndex],
                   }


            if removeMasked:
                # copy mask over to the sampling time!
                mask = np.ma.getmask(out['elevation'])
                out['samplingTime'] = np.ma.masked_array(out['samplingTime'], mask)

                out['elevation'] = np.ma.compress_cols(out['elevation'])
                out['samplingTime'] = np.ma.compress_cols(out['samplingTime'])
                out['frfX'] = np.ma.compress_cols(out['frfX'])
                out['frfY'] = np.ma.compress_cols(out['frfY'])
                out['runupDownLine'] = np.ma.compress_cols(out['runupDownLine'])


            else:
                pass

        else:
            print 'There is no LIDAR data during this time period'
            out = None
        return out

    def getLidarWaveProf(self, removeMasked=True):

        """
        :return:
        """

        self.dataloc = u'projects/tucker/waveprofile/test.ncml'
        self.lidarIndex = self.gettime(dtRound=60)

        if np.size(self.lidarIndex) > 0 and self.lidarIndex is not None:

            out = {'name': nc.chartostring(self.ncfile[u'station_name'][:]),
                   'lat': self.ncfile[u'lidarLatitude'][:],
                   'lon': self.ncfile[u'lidarLongitude'][:],
                   'lidarX': self.ncfile[u'lidarX'][:],
                   'lidarY': self.ncfile[u'lidarY'][:],
                   'frfX': self.ncfile[u'xFRF'][:],
                   'frfY': self.ncfile[u'yFRF'][:],
                   'runupDownLine': self.ncfile['downLineDistance'][:],
                   'waveFreq': self.ncfile['waveFrequency'][:],
                   'time': self.ncfile[u'time'][self.lidarIndex],
                   'WaterLevel': self.ncfile['waterLevel'][self.lidarIndex],
                   'waveHs': self.ncfile['waveHs'][self.lidarIndex],
                   'waveHsIG': self.ncfile['waveHsIG'][self.lidarIndex],
                   'waveHsTot': self.ncfile['waveHsTotal'][self.lidarIndex],
                   'waveSkewness': self.ncfile['waveSkewness'][self.lidarIndex],
                   'waveAsymmetry': self.ncfile['waveAsymmetry'][self.lidarIndex],
                   'waveEnergyDens': self.ncfile['waveEnergyDensity'][self.lidarIndex],
                   'hydroFlag': self.ncfile['hydrodynamicsFlag'][self.lidarIndex],
                   'percentMissing': self.ncfile['percentTimeSeriesMissing'][self.lidarIndex],
                   }


            if removeMasked:
                # copy mask over to the sampling time!
                mask = np.ma.getmask(out['elevation'])
                out['samplingTime'] = np.ma.masked_array(out['samplingTime'], mask)

                out['elevation'] = np.ma.compress_cols(out['elevation'])
                out['samplingTime'] = np.ma.compress_cols(out['samplingTime'])
                out['frfX'] = np.ma.compress_cols(out['frfX'])
                out['frfY'] = np.ma.compress_cols(out['frfY'])
                out['runupDownLine'] = np.ma.compress_cols(out['runupDownLine'])


            else:
                pass

        else:
            print 'There is no LIDAR data during this time period'
            out = None
        return out
