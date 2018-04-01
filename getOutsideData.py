import datetime as DT
import netCDF4 as nc
import os
import numpy as np
import sys

class forecastData:
    def __init__(self, d1):
        """
        Initialization description here
        Data are returned in self.datainex are inclusive at start,end
        """
        self.rawdataloc_wave = []
        self.outputdir = []  # location for outputfiles
        self.d1 = d1  # start date for data grab
        self.timeunits = 'seconds since 1970-01-01 00:00:00'
        self.epochd1 = nc.date2num(self.d1, self.timeunits)
        self.dataLocFRF = u'http://134.164.129.55/thredds/dodsC/FRF/'
        self.dataLocTB = u'http://134.164.129.62:8080/thredds/dodsC/CMTB'
        self.dataLocCHL = u'https://chlthredds.erdc.dren.mil/thredds/dodsC/frf/' #'http://10.200.23.50/thredds/dodsC/frf/'
        self.dataLocNCEP = u'http://nomads.ncep.noaa.gov/pub/data/nccf/com/wave/prod/'#ftpprd.ncep.noaa.gov/pub/data/nccf/com/wave/prod/multi_1.'
        self.dataLocECWMF = u'ftp://data-portal.ecmwf.int/20170808120000/'  # ECMWF forecasts
        assert type(self.d1) == DT.datetime, 'end need to be in python "Datetime" data types'

    def getWW3(self, forecastHour, buoyNumber=44100):
        """
        This funcion will get spectral forecasts from the NCEP nomads server and 
        parse it out to geographic coordinate system. Currently, the data is 
        transformed from oceanographic to meteorological coordinates and from 
        units of m^2 s rad^-1 to m^2 s deg^-1 to maintain FRF gauge data 
        conventions. Spectra are also sorted in ascending order by frequency and
        direction. The functionality associated with transforming the data may be 
        more appropriately located in cmtb/PrepData.
        
        :param forecastHour:
        :param buoyNumber:
        :return: A dictionary with wave directions, frequencies, directional wave
            :key 'wavedirbin':
            :key 'wavefreqbin':
            :key 'dWED': 2dimensional wave spectra [t, freq, dir]
            :key 'lat': latitude
            :key 'lon': longitude
            :key 'time': date time
            spectra, and the timestamps for each spectrum.
        """
        import urllib
        assert type(forecastHour) is str, 'Forecast hour variable must be a string'
        forecastHour = forecastHour.zfill(2)
        urlBack = '/bulls.t%sz/' %forecastHour +'multi_1.%d.spec' %buoyNumber
        ftpURL = self.dataLocNCEP + 'multi_1.' + self.d1.strftime('%Y%m%d') + urlBack
        ftpstream = urllib.urlopen(ftpURL)  # open url
        lines = ftpstream.readlines()  # read the lines of the url into an array of lines
        ftpstream.close()  # close connection with the server
        # # # # # # # # # # # # now the forecast spectra are in lines # # # # # # # # # # # #
        frequencies, directions, forcastDates, forcastDateLines = [], [], [], []

        for ii, line in enumerate(lines):  # read through each line
            split = line.split()   # split the current line up
            if split[0].strip("'")  == 'WAVEWATCH':  # is it the header of the file?
                nFreq = int(split[3])  # number of Frequencies
                nDir = int(split[4])
            elif len(split) == 8 or nFreq - len(frequencies) == len(split) and len(frequencies) != nFreq: # this is frequencies
                frequencies.extend(split)
            elif (len(split) == 7 or nDir - len(directions) == len(split)) and len(directions) != nDir:  # this is directions
                directions.extend(split)
            elif len(split[0]) == 8  and len(split) == 2: ## this is the date line for the beggining of a spectra
                # MPG: Include time component (second entry in date line).
                timestampstr = split[0] + split[1] 
                timestamp = DT.datetime.strptime(timestampstr, '%Y%m%d%H%M%S')
                forcastDates.append(timestamp)
                forcastDateLines.append(ii)
        
        # MPG: convert directions and frequencies from list(string) to 
        # np.array(float).
        directions = np.array(directions).astype('float')
        frequencies = np.array(frequencies).astype('float')

        # MPG: convert directions from radians to degrees.
        directions = np.rad2deg(directions)

        # MPG: convert directions from oceanographic to meteorological 
        # convention to be consistent w/ FRF wave gauge data.
        small_angle = directions < 180.0
        directions[small_angle] = 180.0 + directions[small_angle]
        directions[~small_angle] = directions[~small_angle] - 180.0

        # MPG: sort directions and frequencies.
        didx = directions.argsort()
        fidx = frequencies.argsort()
        directions = directions[didx]
        frequencies = frequencies[fidx]
        
        ## now go back through 'lines' and parse spectra
        spectra = np.ones((len(forcastDateLines), nFreq, nDir), dtype=float) * 1e-8
        buoyNum, lon, lat, Depth, Hm0, Dp, b, c = [], [], [], [], [], [], [], []
        for ll in forcastDateLines:
            numLinesPerSpec = np.ceil(nFreq*float(nDir)/len(lines[ll+2].split())).astype(int)
            buoyStats = lines[ll+1].split()
            if ll == forcastDateLines[0]:  # if its not going to change only grab it once
                buoyNum = int(buoyStats[0].strip("'"))
                lon = float(buoyStats[2])
                lat = float(buoyStats[3])
                Depth = float(buoyStats[4])
            # Hm0.append(float(buoyStats[5])) # these need to be converted to meters ... is this actually wind?
            # Dp.append(float(buoyStats[6])) # these are not exported   ... are these wind?
            # b.append(float(buoyStats[7]))  # not sure what this field is .... wind speed?
            # c.append(float(buoyStats[8]))  # not sure what this field is  ... wind dir?

            tt = np.floor(float(ll) / (numLinesPerSpec - 1)).astype(int)  # time index
            linear = []
            for ss in range(numLinesPerSpec):
                # data =
                linear.extend(lines[ss + ll + 2].split())
            spectra[tt] = np.array(linear, dtype=float).reshape(nDir, nFreq).T
            spectra[tt] = spectra[tt][fidx][:,didx]

        # MPG: convert dWED from rad^-1 to deg^-1 to be consistent w/
        # FRF wave gauge data.
        spectra = spectra*2*np.pi / 180.0
        
        out = {'wavedirbin': directions,
               'wavefreqbin': frequencies,
               'buoyNum': buoyNum,
               'dWED': spectra,
               'lat': lat,
               'lon': lon,
               'Depth': Depth,
               'time': np.array(forcastDates)}

        return out

    def get_CbathyFromFTP(self, dlist, path, timex=True):
        """
        this function downloads argus cbathy bathy data from the argus ftp server
        times must be on the hour or half hour, it will return dates from a list
        provided as dlist.  dlist can be a single point (not in list) in time or
        a list of datetimes
        # written by Ty Hesser
        # modified by Spicer Bak

        dlist: a list of  datetime dataList for cbathy data to be collected
        path: directory to put the cbathy file(s)
        """

        curdir = os.getcwd()  # remembering where i am now
        if not os.path.exists(path):
            os.mkdir(path)
        os.chdir(path)  # changing locations to where data should be downloaded to
        # defining month string to month numbers
        mon = {1: 'Jan', 2: 'Feb', 3: 'Mar', 4: 'Apr', 5: 'May', 6: 'Jun', 7: 'Jul',
               8: 'Aug', 9: 'Sep', 10: 'Oct', 11: 'Nov', 12: 'Dec'}
        # defining days of the week to day number
        dow = {0: 'Mon', 1: 'Tue', 2: 'Wed', 3: 'Thu', 4: 'Fri', 5: 'Sat', 6: 'Sun'}
        # quick data check
        if type(dlist) == DT.datetime:
            dlist = [dlist]  # making it into a list if its a single value
        assert type(dlist[0]) == DT.datetime, 'This function requires datetime dataList'
        # begin looping through data, acquiring cbathy data
        oflist = []
        for ii in range(0, len(dlist)):
            # assert dlist[ii].minute == 0 or dlist[ii].minute == 30, 'the minutes on your datetime object are not 0 or 30'
            if timex == True:
                din_m = DT.timedelta(0, seconds=1) + dlist[
                    ii]  # changing to data processed on the mintue and 31 of the hour
            else:
                din_m = dlist[ii] - DT.timedelta(0, 60)

            # creating month/day hours of timestamp that is being looked for
            yearc = din_m.strftime('%Y')  # string year
            # monthc = din_m.strftime('%m') # string month
            monthc = mon[din_m.month]  # making a month string
            tt = din_m.timetuple()  # making time tuple to make more strings
            dayc = din_m.strftime('%d')  # making a day string
            hourc = din_m.strftime('%H')
            mmc = str(tt.tm_min)  # making a minute string
            if len(mmc) == 1:
                mmc = '0' + mmc
            ssc = str(tt.tm_sec)
            if len(ssc) == 1:
                ssc = '0' + ssc
            # creating epoch time
            eptm = str(int(nc.date2num(din_m, 'seconds since 1970-01-01')))

            # creating the url to download the cbathy data from
            # frfserver = "'\\134.164.129.42\cil\argus02b\'"                       # at the frf
            OSUserver = "ftp://cil-ftp.coas.oregonstate.edu/pub/argus02b/"  # the oregon state server
            svr = OSUserver + yearc + "/cx/"  # server base
            daynum = str(tt.tm_yday)  # day number in a year
            fldr = daynum + "_" + monthc + "." + dayc  # defining the folder (date) structure to be used
            fname = "/" + eptm + '.' + dow[tt.tm_wday] + '.' + monthc + '.' + dayc + \
                    "_" + hourc + '_' + mmc + \
                    '_' + ssc + '.GMT.' + yearc + ".argus02b.cx.cBathy.mat"
            # fname = '/1445709540.Sat.Oct.24_17_59_00.GMT.2015.argus02b.cx.cBathy.mat'  # copied and pasted
            if timex == True:
                fname = '/*timex.merge.mat'
            addr = svr + fldr + fname
            print "checking " + fldr + fname
            if os.path.isfile(fname[1:]):
                print "already downloaded: %s" % fname
            elif not os.path.isfile(fname[1:]):
                # try:
                # dlfname = wget.download(addr)
                os.system('wget %s' % addr)
                # oflist.append(dlfname)
                # print 'Retrieved %s' %dlfname
                # except IOError:
                #    print "There is no file on the server, It's probably dark outside"

        os.chdir(curdir)
        #    import urllib
        #    urllib.urlretrieve(fn4)
        return oflist
