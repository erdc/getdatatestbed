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
import netCDF4 as nc
import numpy as np
import pandas as pd

try:
    import sblib as sb
except ImportError:
    whoami = check_output('whoami', shell=True)
    sys.path.append('C:\Users\spike\Documents\Code_Repositories\sblib')
    sys.path.append('/home/%s/repos/sblib' % whoami[:whoami.index('\n')])
    import sblib as sb

class getObs:
    """
    Note d1 and d2 have to be in date-time formats
    are all data set times in UTC?
    need to write error handling, what to do if there's no data ??

    """

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

             :param dtRound: the time delta of the data out of interest, default minute (60 second)

        """
        # TODO find a way to pull only hourly data or regular interval of desired time
        # todo this use date2index and create a list of dates see help(nc.date2index)
        try:

            self.ncfile = nc.Dataset(self.FRFdataloc + self.dataloc) #loads all of the netCDF file
            #            try:
            self.alltime = nc.num2date(self.ncfile['time'][:], self.ncfile['time'].units,
                                       self.ncfile['time'].calendar) # converts all epoch time to datetime objects
            for i, date in enumerate(self.alltime):  # rounds time to nearest
                self.alltime[i] = self.roundtime(dt=date, roundto=dtRound)

            mask = (self.alltime >= self.d1) & (self.alltime < self.d2)  # boolean true/false of time
            idx = np.where(mask)[0]

            assert len(idx) > 0, 'no data locally, check CHLthredds'
            print "Data Gathered From Local Thredds Server"

        except (RuntimeError, NameError, AssertionError):  # if theres any error try to get good data from next location
            self.ncfile = nc.Dataset(self.chlDataLoc + self.dataloc)
            self.alltime = nc.num2date(self.ncfile[ 'time'][:], self.ncfile['time'].units,
                                       self.ncfile['time'].calendar)
            for i, date in enumerate(self.alltime):
                self.alltime[i] = self.roundtime(dt=date, roundto=dtRound)
            # mask = (sb.roundtime(self.ncfile['time'][:]) >= self.epochd1) & (sb.roundtime(self.ncfile['time'][:]) < self.epochd2)\

            mask = (self.alltime >= self.d1) & (self.alltime < self.d2)  # boolean true/false of time
            idx = np.where(mask)[0]

            try:
                assert len(idx) > 0, ' There are no data within the search parameters for this gauge'
                print "Data Gathered from CHL thredds Server"
            except AssertionError:
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
            gname = 'AWAC04 - 11m'
            self.dataloc = 'oceanography/waves/awac-11m/awac-11m.ncml'
        elif gaugenumber == 3 or gaugenumber == 'awac-8m':
            gname = 'AWAC 8m'
            self.dataloc = 'oceanography/waves/awac-8m/awac-8m.ncml'
        elif gaugenumber == 4 or gaugenumber == 'awac-6m':
            gname = 'AWAC 6m'
            self.dataloc = 'oceanography/waves/awac-6m/awac-6m.ncml'
        elif gaugenumber == 5 or gaugenumber == 'awac-4.5m':
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
        else:
            gname = 'There Are no Gauge numbers here'
            raise NameError('Bad Gauge name, specify proper gauge name/number')
        # parsing out data of interest in time
        try:
            self.wavedataindex = self.gettime(dtRound=roundto * 60)
            assert np.array(self.wavedataindex).all() != None, 'there''s no data in your time period'

            if np.size(self.wavedataindex) >= 1:
                # consistant for all wave gauges
                self.snaptime = nc.num2date(self.ncfile['time'][self.wavedataindex], self.ncfile['time'].units,
                                            self.ncfile['time'].calendar)
                for num in range(0, len(self.snaptime)):
                    self.snaptime[num] = self.roundtime(self.snaptime[num], roundto=roundto * 60)
                #                if roundto != 30:
                #                    self.wavedataindex=self.cliprecords(self.snaptime)
                try:
                    depth = self.ncfile['depth'][:]
                except IndexError:
                    depth = self.ncfile['depthP'][-1]

                wavespec = {'time': self.snaptime,
                            'epochtime': self.ncfile['time'][self.wavedataindex],
                            'name': str(self.ncfile.title),
                            'wavefreqbin': self.ncfile['waveFrequency'][:],
                            'lat': self.ncfile['lat'][:],
                            'lon': self.ncfile['lon'][:],
                            'depth': depth,
                            'Hs': self.ncfile['waveHs'][self.wavedataindex],
                            'peakf': self.ncfile['wavePeakFrequency'][self.wavedataindex]
                            }
                try:
                    wavespec['wavedirbin'] = self.ncfile['waveDirectionBins'][:]
                    wavespec['dWED'] = self.ncfile['directionalWaveEnergyDensity'][self.wavedataindex, :, :]
                    wavespec['waveDp'] = self.ncfile['wavePeakDirectionPeakFrequency'][self.wavedataindex]
                    wavespec['fspec'] = self.ncfile['waveEnergyDensity'][self.wavedataindex, :]
                    wavespec['waveDm'] = self.ncfile['waveMeanDirection'][self.wavedataindex]
                    wavespec['qcFlagE'] = self.ncfile['qcFlagE'][self.wavedataindex]
                    wavespec['qcFlagD'] = self.ncfile['qcFlagD'][self.wavedataindex]
                    wavespec['a1'] = self.ncfile['waveA1Value'][self.wavedataindex, :]
                    wavespec['a2'] = self.ncfile['waveA2Value'][self.wavedataindex, :]
                    wavespec['b1'] = self.ncfile['waveB1Value'][self.wavedataindex, :]
                    wavespec['b2'] = self.ncfile['waveB2Value'][self.wavedataindex, :]
                except IndexError:
                    # this should throw when gauge is non directional
                    # wavespec['peakf'] = self.ncfile['waveFp'][self.wavedataindex],
                    wavespec['wavedirbin'] = np.arange(0, 360, 90)  # 90 degree bins
                    wavespec['dWED'] = np.ones(
                        [len(self.wavedataindex), len(wavespec['wavefreqbin']), len(wavespec['wavedirbin'])])  # * 1e-8
                    wavespec['waveDp'] = np.zeros(len(self.wavedataindex))
                    wavespec['fspec'] = self.ncfile['waveEnergyDensity'][self.wavedataindex, :]
                    wavespec['depthp'] = self.ncfile['depthP'][self.wavedataindex]
                    # wavespec['qcFlagE'] = self.ncfile['qcFlagE'][self.wavedataindex]
                    # multiply the freq spectra for all directions
                    wavespec['dWED'] = wavespec['dWED'] * wavespec['fspec'][:, :, np.newaxis]/len(wavespec['wavedirbin'])
                    wavespec['qcFlagE'] = self.ncfile['qcFlagE'][self.wavedataindex]

                return wavespec

        except (RuntimeError, AssertionError):
            print '<<ERROR>> Retrieving data from %s\n in this time period start: %s  End: %s' % (
            gname, self.d1, self.d2)
            wavespec = None
            return wavespec

    def getCurrents(self, roundto=1):
        """
        This function pulls down the currents data from the Thredds Server


            :param roundto:
                the time over which the wind record exists
                ie data is collected in 10 minute increments
                data is rounded to the nearst [roundto] (default 1 min)
        """
        self.dataloc = 'oceanography/currents/awac04/awac04.ncml'

        currdataindex = self.gettime(dtRound=roundto * 60)
        # _______________________________________
        # get the actual current data
        if np.size(currdataindex) > 1:
            curr_aveU = self.ncfile['aveU'][currdataindex]  # pulling depth averaged Eastward current
            curr_aveV = self.ncfile['aveV'][currdataindex]  # pulling depth averaged Northward current
            curr_spd = self.ncfile['currentSpeed'][currdataindex]  # currents speed [m/s]
            curr_dir = self.ncfile['currentDirection'][currdataindex]  # current from direction [deg]
            curr_time = nc.num2date(self.ncfile['time'][currdataindex], self.ncfile['time'].units,
                                    self.ncfile['time'].calendar)
            for num in range(0, len(self.curr_time)):
                self.curr_time[num] = self.roundtime(self.curr_time[num], roundto=roundto * 60)
            self.curpacket = {
                'name': str(self.ncfile.title),
                'time': curr_time,
                'aveU': curr_aveU,
                'aveV': curr_aveV,
                'speed': curr_spd,
                'dir': curr_dir,
                'lat': self.ncfile['lat'][:],
                'lon': self.ncfile['lon'][:],
                'depth': self.ncfile['depth'][:],
                # Depth is calculated by: depth = -xducerD + blank + (binSize/2) + (numBins * binSize)
                'meanP': self.ncfile['meanPressure'][currdataindex],

            }
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

        # ______________________________________
        if np.size(self.winddataindex) > 0 and self.winddataindex is not None:
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

            self.windtime = nc.num2date(self.ncfile['time'][self.winddataindex], self.ncfile['time'].units,
                                        self.ncfile['time'].calendar)  # wind time
            for num in range(0, len(self.windtime)):
                self.windtime[num] = self.roundtime(self.windtime[num], roundto=collectionlength * 60)
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
                'lat': self.ncfile['lat'][:],  # latitude
                'lon': self.ncfile['lon'][:],  # longitde
                'gaugeht': gaugeht,
            }
            if (windpacket['qcflagD'] == 3).all() or (windpacket['qcflagS'] == 3).all():
                print "Wind querey returned all bad data for speed or direction"
                windpacket = None
            return windpacket
        else:
            print 'ERROR: There is no Wind Data for this time period !!!'
            windpacket = None
            return windpacket

    def getWL(self, collectionlength=6):
        """
        This function retrieves the water level data from the FDIF server
        WL data on server is NAVD

        collection length is the time over which the wind record exists
            ie data is collected in 10 minute increments
            data is rounded to the nearst [collectionlength] (default 6 min)
        """
        self.dataloc = 'oceanography/waterlevel/11/11.ncml'  # this is the back end of the url for waterlevel
        self.WLdataindex = self.gettime(dtRound=collectionlength * 60)

        if np.size(self.WLdataindex) > 1:
            self.WL = self.ncfile['waterLevelHeight'][self.WLdataindex]
            self.WLtime = nc.num2date(self.ncfile['time'][self.WLdataindex], self.ncfile['time'].units,
                                      self.ncfile['time'].calendar)
            for num in range(0, len(self.WLtime)):
                self.WLtime[num] = self.roundtime(self.WLtime[num], roundto=collectionlength * 60)

            self.WLpacket = {
                'name': str(self.ncfile.title),
                'WL': self.ncfile['waterLevelHeight'][self.WLdataindex],
                'time': self.WLtime,
                'lat': self.ncfile['lat'][:],
                'lon': self.ncfile['lon'][:],
                # 'surge': self.ncfile['surge'][self.WLtime],
                'predictedWL': self.ncfile['predictedWaterLevelHeight'][self.WLdataindex]
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
        from GetDataTestBed import download_grid_data as DGD
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
        self.dataloc = u'survey/transect/transect.ncml'  # location of the gridded surveys
        try:
            self.bathydataindex = self.gettime()
        except IOError:  # when data are not on CHL thredds
            self.bathydataindex = []
        # logic to handle no transects in date range
        if len(self.bathydataindex) == 1:
            idx = self.bathydataindex
        elif len(self.bathydataindex) < 1 & method == 1:
            # there's no exact bathy match so find the max negative number where the negitive
            # numbers are historical and the max would be the closest historical
            val = (max([n for n in (self.ncfile['time'][:] - self.epochd1) if n < 0]))
            idx = np.where((self.ncfile['time'][:] - self.epochd1) == val)[0][0]
            print 'Bathymetry is taken as closest in HISTORY - operational'
        elif len(self.bathydataindex) < 1 and method == 0:
            idx = np.argmin(np.abs(self.ncfile['time'][:] - self.d1))  # closest in time
            print 'Bathymetry is taken as closest in TIME - NON-operational'
        elif len(self.bathydataindex) > 1:
            val = (max([n for n in (self.ncfile['time'][:] - self.d1) if n < 0]))
            idx = np.where((self.ncfile['time'] - self.d1) == val)[0][0]
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


        if np.size(idx) > 0:
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

    def getBathyGridFromNC(self, method, removeMask=True):
        """
        This function gets the bathymetric data from the thredds server, currently designed for the bathy duck experiment
        method == 1  - > 'Bathymetry is taken as closest in HISTORY - operational'
        method == 0  - > 'Bathymetry is taken as closest in TIME - NON-operational'
        :param
        :return:

        """
        self.dataloc = u'survey/transect/transect.ncml'  # location of the gridded surveys
        try:
            self.bathydataindex = self.gettime()  # getting the index of the grid
        except IOError:
            self.bathydataindex = []
        if len(self.bathydataindex) == 1:
            idx = self.bathydataindex
        elif len(self.bathydataindex) < 1 & method == 1:
            # there's no exact bathy match so find the max negative number where the negitive
            # numbers are historical and the max would be the closest historical
            val = (max([n for n in (self.ncfile['time'][:] - self.epochd1) if n < 0]))
            idx = np.where((self.ncfile['time'][:] - self.epochd1) == val)[0][0]
            print 'Bathymetry is taken as closest in HISTORY - operational'
        elif len(self.bathydataindex) < 1 and method == 0:
            idx = np.argmin(np.abs(self.ncfile['time'][:] - self.d1))  # closest in time
            print 'Bathymetry is taken as closest in TIME - NON-operational'
        elif len(self.bathydataindex) > 1:
            val = (max([n for n in (self.ncfile['time'][:] - self.d1) if n < 0]))
            idx = np.where((self.ncfile['time'] - self.d1) == val)[0][0]

            print 'The closest in history to your start date is %s\n' % nc.num2date(self.gridTime[idx],
                                                                                    self.ncfile['time'].units)
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
        xCoord = self.ncfile['FRF_Xshore'][:]
        yCoord = self.ncfile['FRF_Yshore'][:]
        lat = self.ncfile['lat'][:]
        lon = self.ncfile['lon'][:]
        northing = self.ncfile['northing'][:]
        easting = self.ncfile['easting'][:]
        if removeMask == True:
            xCoord = xCoord[~np.all(elevation_points.mask, axis=0)]
            yCoord = yCoord[~np.all(elevation_points.mask, axis=1)]
            lon = lon[~np.all(elevation_points.mask, axis=0)]
            lat = lat[~np.all(elevation_points.mask, axis=1)]
            northing = northing[~np.all(elevation_points.mask, axis=0)]
            easting = easting[~np.all(elevation_points.mask, axis=1)]
            elevation_points = elevation_points[~np.all(elevation_points.mask, axis=1), :]  #
            elevation_points = elevation_points[:, ~np.all(elevation_points.mask, axis=0)]
        if elevation_points.ndim == 2:
            elevation_points = np.ma.expand_dims(elevation_points, axis=0)

        time = (self.ncfile['time'][idx], self.ncfile['time'].units)
        print 'Sim start: %s\nSim End: %s\nSim bathy chosen: %s' % (self.d1, self.d2,
                                                                    nc.num2date(self.ncfile['time'][idx],
                                                                                self.ncfile['time'].units))
        print 'Bathy is %s old' % (self.d2 - nc.num2date(self.ncfile['time'][idx], self.ncfile['time'].units))[0]

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
        locDict = sb.FRFcoord(xloc[0], yloc[0])
        return locDict

    def getBathyGridcBathy(self):
        """
        this functin gets the cbathy data from the below address, assumes fill value of -999
        :return:
        """
        fillValue = -999.99  # assumed fill value from
        try:
            cbathyloc2 = self.FRFdataloc + u'projects/bathyduck/BathyDuck-ocean_bathy_argus_201510.nc'
            cbfile = nc.Dataset(cbathyloc2)
        except:
            cbathyloc2 = self.chlDataLoc + u'projects/bathyduck/BathyDuck-ocean_bathy_argus_201510.nc'
            cbfile = nc.Dataset(cbathyloc2)
        ed1 = nc.date2num(self.d1, 'seconds since 1970-01-01')
        ed2 = nc.date2num(self.d2, 'seconds since 1970-01-01')
        emask = (cbfile['time'][:] >= ed1) & (cbfile['time'][:] < ed2)
        # mask = (time > d1) & (time < d2)
        # assert (emask == mask).all(), 'epoch time is not working'
        idx = np.where(emask)[0]
        try:
            maskedElev = (cbfile['depthKF'][idx, :, :] == fillValue)
            depthKF = np.ma.array(cbfile['depthKF'][idx, :, :], mask=maskedElev)
            time = (cbfile['time'][idx], 'seconds since 1970-01-01')
            cbdata = {'time': time,
                      'xm': cbfile['xm'][:],
                      'ym': cbfile['ym'][:],
                      'depth': cbfile['depth'][idx, :, :],
                      'depthKF': depthKF,
                      'depthKFError': cbfile['depthKFError'][idx, :, :],
                      'fB': cbfile['fB'][idx, :, :, :],
                      'k': cbfile['k'][idx, :, :, :]}
            print 'Grabbed BathyDuck, cBathy Data'
        except IndexError:  # there's no data in the Cbathy
            cbdata = None

        return cbdata

    def getLidarRunup(self):
        """

        :return:
        """
        self.dataloc = 'oceanography/waves/lidarWaveRunup/lidarWaveRunup.ncml'
        self.lidarIndex = self.gettime(dtRound=60)
        if np.size(self.lidarIndex) > 0 and self.lidarIndex is not None:
            out = {'name': nc.chartostring(self.ncfile[u'station_name'][:]),
                   'lat': self.ncfile[u'lat'][:],  # Coordintes
                   'lon': self.ncfile[u'lon'][:],
                   'frfX': self.ncfile[u'frfX'][:],
                   'frfY': self.ncfile[u'frfY'][:],
                   'time': self.ncfile[u'time'][self.lidarIndex],
                   'totalWaterLevel': self.ncfile['totalWaterLevel'][self.lidarIndex],
                   'QCtotalWaterLevel': self.ncfile['QCTotalWaterLevel'][self.lidarIndex],
                   'percentMissing': self.ncfile['percentTimeSeriesMissing'][self.lidarIndex],
                   'samplingTime': self.ncfile['samplingTime'][self.lidarIndex, :],
                   'runupX': self.ncfile['runupX'][self.lidarIndex, :],
                   'runupY': self.ncfile['runupY'][self.lidarIndex, :],
                   'runupZ': self.ncfile['runupZ'][self.lidarIndex, :],
                   'runupDownLine': self.ncfile['runupDownLineDistance'][self.lidarIndex, :],
                   }
        else:
            print 'There is no Data during this time period'
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
        # location of the gridded surveys
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

            alt_coords = sb.sblib.FRFcoord(alt_lon, alt_lat)

            if removeMasked:
                self.altpacket = {
                    'name': str(self.ncfile.title),
                    'time': np.array(self.alt_time[~alt_be.mask]),
                    'lat': alt_lat,
                    'PKF': np.array(alt_pkf[~alt_be.mask]),
                    'lon': alt_lon,
                    'FRF_X': alt_coords['FRF_X'],
                    'FRF_Y': alt_coords['FRF_Y'],
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
                    'FRF_X': alt_coords['FRF_X'],
                    'FRF_Y': alt_coords['FRF_Y'],
                    'stationName': alt_stationname,
                    'gageName': gagename,
                    'timeStart': self.alt_timestart,
                    'timeEnd': self.alt_timeend,
                    'bottomElev': alt_be}

            return self.altpacket
        else:
            print 'No %s data found for this period' %(gagename)


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
        self.crunchDataLoc = u'http://134.164.129.62:8080/thredds/dodsC/CMTB/'
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

    def gettime(self, dtRound=30):
        """
        this function opens the netcdf file, pulls down all of the time, then pulls the dates of interest
        from the THREDDS (data loc) server based on d1,d2, and data location
        it returns the indicies in the NCML file of the dates d1>=time>d2
        INPUTS:

            d1: start time - pulled from self
            d2: end time  - pulled from self
            dataloc: location of the data to search through
            :param dtRound: the time delta of the data out of interest

        """
        # TODO find a way to pull only hourly data or regular interval of desired time
        # todo this use date2index and create a list of dates see help(nc.date2index)
        try:

            self.ncfile = nc.Dataset(self.FRFdataloc + self.dataloc)
            #            try:
            self.alltime = nc.num2date(self.ncfile['time'][:], self.ncfile['time'].units,
                                       self.ncfile['time'].calendar)
            for i, date in enumerate(self.alltime):
                self.alltime[i] = self.roundtime(dt=date, roundto=dtRound)

            mask = (self.alltime >= self.d1) & (self.alltime < self.d2)  # boolean true/false of time
            # mask = (sb.roundtime(self.ncfile['time'][:]) >= self.epochd1) & (sb.roundtime(self.ncfile['time'][:]) < self.epochd2)\

            idx = np.where(mask)[0]
            assert len(idx) > 0, 'no data locally, checking CHLthredds'
            print "Data Gathered From Local Thredds Server"
        except (RuntimeError, NameError, AssertionError):  # if theres any error try to get good data from next location
            self.ncfile = nc.Dataset(self.chlDataLoc + self.dataloc)
            self.alltime = nc.num2date(self.ncfile['time'][:], self.ncfile['time'].units,
                                       self.ncfile['time'].calendar)
            for i, date in enumerate(self.alltime):
                self.alltime[i] = self.roundtime(dt=date, roundto=dtRound)
            # mask = (sb.roundtime(self.ncfile['time'][:]) >= self.epochd1) & (sb.roundtime(self.ncfile['time'][:]) < self.epochd2)\

            mask = (self.alltime >= self.d1) & (self.alltime < self.d2)  # boolean true/false of time
            idx = np.where(mask)[0]

            try:
                assert len(idx) > 0, ' There are no data within the search parameters for this gauge'
                print "Data Gathered from CHL thredds Server"
            except AssertionError:
                idx = None

        return idx

    def getStwaveField(self, var, prefix, local=True, ijLoc=None):
        """

        :param Local: defines whether the data is from the nested simulation or the regional simulation
        : param ijLoc: x or y or (x,y) tuple of location of interest
        :return:
        """
        if local == True:
            grid = 'Local'
        elif local == False:
            grid = 'Regional'
        if prefix == 'CBHPStatic' and local == True:  # this is needed because projects are stored in weird place
            ncfile = nc.Dataset('http://crunchy:8080/thredds/dodsC/CMTB/projects/bathyDuck_SingleBathy_CBHP/Local_Field/Local_Field.ncml')
        elif prefix =='CBHPStatic' and local == False:
            ncfile = nc.Dataset('http://crunchy:8080/thredds/dodsC/CMTB/projects/bathyDuck_SingleBathy_CBHP/Regional_Field/Regional_Field.ncml')
        else:
            ncfile = nc.Dataset(self.crunchDataLoc + u'/%s_STWAVE_data/%s_Field/%s_Field.ncml' % (prefix, grid, grid))
        assert var in ncfile.variables.keys(), 'variable called is not in file please use\n%s' % ncfile.variables.keys()
        mask = (ncfile['time'][:] >= nc.date2num(self.d1, ncfile['time'].units)) & (
            ncfile['time'][:] <= nc.date2num(self.d2, ncfile['time'].units))
        idx = np.where(mask)[0]
        print 'getting %s STWAVE  %s %s Data' % (prefix, grid, var)
        # now creating tool to remove single data point
        if ijLoc != None:
            assert len(ijLoc) == 2, 'if giving a postion, must be a tuple of i, j location (of length 2)'
            if type(ijLoc[0]) == int:
                x = ncfile[var].shape[1] - ijLoc[0]  # the data are stored with inverse indicies to grid node locations
                y = ijLoc[1]  # use location given by function call
            else: # ijLoc[0] == slice:
                x = ijLoc[0]
                y = ijLoc[1]

        else:
            x = slice(None)  # take entire data
            y = slice(None)  # take entire data

        if ncfile[var].shape[0] > 100: # looping through ... if necessicary
            list = np.round(np.linspace(idx.min(), idx.max(), idx.max()-idx.min(), endpoint=True, dtype=int))
            # if idx.max() not in list:
            #     list = np.append(list, idx.max())
            if len(list) < 100:
                    bathy = np.array(ncfile[var][np.squeeze(list)])
                    time = nc.num2date(ncfile['time'][np.squeeze(list)], ncfile['time'].units)
            else:
                for num, minidx in enumerate(list):
                    if len(list) < 100:
                        bathy = np.array(ncfile[var][np.squeeze(list)])
                        time = nc.num2date(ncfile['time'][np.squeeze(list)], ncfile['time'].units)
                    elif num == 0:
                        bathy = np.array(ncfile[var][range(minidx, list[num + 1]), x, y])
                        time = nc.num2date(np.array(ncfile['time'][range(minidx, list[num + 1])]), ncfile['time'].units)
                    elif minidx == list[-1]:
                        lastIdx = (idx - minidx)[(idx - minidx) >= 0] + minidx
                        bathy = np.append(bathy, ncfile[var][lastIdx, x, y], axis=0)
                        time = np.append(time, nc.num2date(ncfile['time'][lastIdx], ncfile['time'].units), axis=0)
                    else:
                        bathy = np.append(bathy, ncfile[var][range(minidx, list[num + 1]), x, y], axis=0)
                        time = np.append(time,
                                         nc.num2date(ncfile['time'][range(minidx, list[num + 1])], ncfile['time'].units),
                                         axis=0)
        else:
            bathy = ncfile[var][:, x, y]

        field = {'time': time,
                 var: bathy,
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

            if gaugenumber in [0, 'waverider-26m', 'Waverider-26m', '26m']:
                # 26 m wave rider
                self.dataloc = '%s_data/FRF_CMTB_%s_waverider-26m.ncml'  %(prefix, prefix)
                gname = '26m Waverider Buoy'
            elif gaugenumber in [1, 'Waverider-17m', 'waverider-17m']:
                # 2D 17m waverider
                self.dataloc = '%s_data/FRF_CMTB_%s_waverider-17m.ncml'  %(prefix, prefix)
                gname = '17m Waverider Buoy'
            elif gaugenumber in [2, 'AWAC-11m', 'awac-11m', 'Awac-11m']:
                gname = 'AWAC04 - 11m'
                self.dataloc = '%s_data/FRF_CMTB_%s_awac11m.ncml'  %(prefix, prefix)
            elif gaugenumber in [3, 'awac-8m', 'AWAC-8m', 'Awac-8m', 'awac 8m',
                                 '8m-Array', '8m Array', '8m array', '8m-array']:
                gname = 'AWAC 8m'
                self.dataloc = '%s_data/FRF_CMTB_%s_awac-8m.ncml' %(prefix, prefix)
            elif gaugenumber in [4, 'awac-6m', 'AWAC-6m']:
                gname = 'AWAC 6m'
                self.dataloc = '%s_data/FRF_CMTB_%s_awac-6m.ncml' %(prefix, prefix)
            elif gaugenumber in [5, 'awac-4.5m', 'Awac-4.5m']:
                gname = 'AWAC 4.5m'
                self.dataloc = '%s_data/FRF_CMTB_%s_awac-4.5m.ncml' %(prefix, prefix)
            elif gaugenumber [6, 'adop-3.5m', 'aquadopp 3.5m']:
                gname = 'Aquadopp 3.5m'
                self.dataloc = '%s_data/FRF_CMTB_%s_adop-3.5m.ncml' %(prefix, prefix)
            elif gaugenumber in [8, 'xp200m', 'xp200']:
                gname = 'Paros xp200m'
                self.dataloc = '%s_data/FRF_CMTB_%s_xp200m.ncml' %(prefix, prefix)
            elif gaugenumber in [9, 'xp150m', 'xp150']:
                gname = 'Paros xp150m'
                self.dataloc = '%s_data/FRF_CMTB_%s_xp150m.ncml' %(prefix, prefix)
            elif gaugenumber == 10 or gaugenumber == 'xp125m':
                gname = 'Paros xp125m'
                self.dataloc = '%s_data/FRF_CMTB_%s_xp125m.ncml' %(prefix, prefix)
            elif gaugenumber == 11 or gaugenumber == 'xp100m':
                gname = 'Paros xp100m'
                self.dataloc = '%s_data/FRF_CMTB_%s_xp100m.ncml' %(prefix, prefix)
            else:
                gname = 'There Are no Gauge numbers here'
                raise NameError('Bad Gauge name, specify proper gauge name/number')
            # parsing out data of interest in time
            try:
                self.wavedataindex = self.gettime()
                assert np.array(self.wavedataindex).all() != None, 'there''s no data in your time period'

                if np.size(self.wavedataindex) >= 1:
                    wavespec = {'time': nc.num2date(self.ncfile['time'], self.ncfile['time'].units),
                                'name': nc.chartostring(self.ncfile['station_name'][:]),
                                'wavefreqbin': self.ncfile['waveFrequency'][:],
                                'lat': self.ncfile['lat'][:],
                                'lon': self.ncfile['lon'][:],
                                'Hs': self.ncfile['waveHs'][self.wavedataindex],
                                'peakf': self.ncfile['wavePeakFrequency'][self.wavedataindex],
                                'wavedirbin': self.ncfile['waveDirectionBins'][:],
                                'dWED': self.ncfile['directionalWaveEnergyDensity'][self.wavedataindex, :, :],
                                'waveDp': self.ncfile['wavePeakDirectionPeakFrequency'][self.wavedataindex],  # 'waveDp'][self.wavedataindex]
                                'waveDm': self.ncfile['waveMeanDirection'][self.wavedataindex],
                                'WED': self.ncfile['waveEnergyDensity'][self.wavedataindex, :],}

            except (RuntimeError, AssertionError):
                print '<<ERROR>> Retrieving data from %s\n in this time period start: %s  End: %s' % (
                    gname, self.d1, self.d2)
                wavespec = None
            return wavespec