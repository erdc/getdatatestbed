import datetime as DT
import netCDF4 as nc
import os

class forecastData:
    def __init__(self, d1):
        """
        Initialization description here
        Data are returned in self.datainex are inclusive at d1,d2
        """
        self.rawdataloc_wave = []
        self.outputdir = []  # location for outputfiles
        self.d1 = d1  # start date for data grab
        self.d2 = d2  # end data for data grab
        self.timeunits = 'seconds since 1970-01-01 00:00:00'
        self.epochd1 = nc.date2num(self.d1, self.timeunits)
        self.epochd2 = nc.date2num(self.d2, self.timeunits)
        self.comp_time()
        self.dataLocFRF = u'http://134.164.129.55/thredds/dodsC/FRF/'
        self.dataLocTB = u'http://134.164.129.62:8080/thredds/dodsC/CMTB'
        self.dataLocCHL = u'https://chlthredds.erdc.dren.mil/thredds/dodsC/frf/' #'http://10.200.23.50/thredds/dodsC/frf/'
        self.dataLocNCEP = u'ftp://ftpprd.ncep.noaa.gov/pub/data/nccf/com/wave/prod/multi_1.'
        assert type(self.d2) == DT.datetime, 'd1 need to be in python "Datetime" data types'
        assert type(self.d1) == DT.datetime, 'd2 need to be in python "Datetime" data types'

    def getWWIII(self):
        import urllib, tarfile

        forecastHour = '00'
        ftpURL = self.dataLocNCEP + d1.strftime('%Y%m%d') \
                 + '/multi_1.t%sz.spec_tar.gz' % forecastHour
        spectralFile = './multi_1.44100.spec'
        ftpstream = urllib.urlopen(ftpURL)
        tar = tarfile.open(fileobj=ftpstream, mode='r|bz2')
        ftpstream.close()
        f = tar.extractfile(spectralFile)
        content = f.readlines()
        tar.close()

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
