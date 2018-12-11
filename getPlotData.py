from .getDataFRF import getObs
from testbedutils.anglesLib import vectorRotation
import datetime as DT
import numpy as np
from testbedutils.gridTools import findNearestUnstructNode
import testbedutils.sblib as sb
import netCDF4 as nc

def alt_PlotData(name, mod_time, mod_times, THREDDS='FRF'):

    """
    This function is just to remove clutter in my plot functions
    all it does is pull out altimeter data and put it into the appropriate dictionary keys.
    If None, it will return masked arrays.

    :param name: name of the altimeter you want (Alt03, Alt04, Alt05)
    :param mod_time: start time of the model
    :param time: array of model datetimes
    :param mod: array of model observations at that instrument location corresponding to variable "comp_time"

    :return: Altimeter data dictionary with keys:
            'zb' - elevation
            'name' - gage name
            'time' - timestamps of the data
            'xFRF' - position of the gage
            'plot_ind' - this just tells it which data point it should plot for the snapshots
    """
    t1 = mod_times[0] - DT.timedelta(days=0, hours=0, minutes=3)
    t2 = mod_times[-1] + DT.timedelta(days=0, hours=0, minutes=3)
    frf_Data = getObs(t1, t2, THREDDS)

    try:
        dict = {}
        alt_data = frf_Data.getALT(name)
        dict['zb'] = alt_data['bottomElev']
        dict['time'] = alt_data['time']
        dict['name'] = alt_data['gageName']
        dict['xFRF'] = round(alt_data['xFRF'])
        plot_ind = np.where(abs(dict['time'] - mod_time) == min(abs(dict['time'] - mod_time)), 1, 0)
        dict['plot_ind'] = plot_ind
        dict['TS_toggle'] = True


    except:
        # just make it a masked array
        dict = {}
        comp_time = [mod_times[ii] + (mod_times[ii + 1] - mod_times[ii]) / 2 for ii in range(0, len(mod_times) - 1)]
        dict['time'] = comp_time
        fill_x = np.ones(np.shape(dict['time']))
        dict['zb'] = np.ma.array(fill_x, mask=np.ones(np.shape(dict['time'])))
        dict['name'] = name
        dict['xFRF'] = np.ma.array(np.ones(1), mask=np.ones(1))
        fill_ind = np.zeros(np.shape(dict['time']))
        fill_ind[0] = 1
        dict['plot_ind'] = fill_ind
        dict['TS_toggle'] = False


    return dict

def wave_PlotData(name, mod_time, time, THREDDS='FRF'):
    """
    This function is just to remove clutter in my plotting scripts
    all it does is pull out altimeter data and put it into the appropriate dictionary keys.
    If None, it will return masked arrays.

    :param t1: start time you want to pull (a datetime, not a string)
    :param t2: end time you want to pull (a datetime, not a string)
    :param name: name of the wave gage you want
    :param mod_time: start time of the model
    :param time: array of model datetimes

    :return: Altimeter data dictionary with keys:
            'Hs' - significant wave height
            'name' - gage name
            'wave_time' - timestamps of the data
            'xFRF' - position of the gage
            'plot_ind' - this just tells it which data point it should plot for the snapshots
    """

    t1 = time[0] - DT.timedelta(days=0, hours=0, minutes=3)
    t2 = time[-1] + DT.timedelta(days=0, hours=0, minutes=3)

    frf_Data = getObs(t1, t2, THREDDS)

    try:

        dict = {}
        wave_data = frf_Data.getWaveSpec(gaugenumber=name)
        cur_data = frf_Data.getCurrents(name)

        dict['name'] = name
        dict['wave_time'] = wave_data['time']
        dict['cur_time'] = cur_data['time']
        dict['Hs'] = wave_data['Hs']
        dict['xFRF'] = wave_data['xFRF']
        dict['plot_ind'] = np.where(abs(dict['wave_time'] - mod_time) == min(abs(dict['wave_time'] - mod_time)), 1, 0)
        dict['plot_ind_V'] = np.where(abs(dict['cur_time'] - mod_time) == min(abs(dict['cur_time'] - mod_time)), 1, 0)
        # rotate my velocities!!!
        test_fun = lambda x: vectorRotation(x, theta=360 - (71.8 + (90 - 71.8) + 71.8))
        newV = [test_fun(x) for x in zip(cur_data['aveU'], cur_data['aveV'])]
        dict['U'] = zip(*newV)[0]
        dict['V'] = zip(*newV)[1]
        dict['TS_toggle'] = True

    except:
        print(('No data at %s!  Will return masked array.') %name)
        # just make it a masked array
        dict = {}
        dict['wave_time'] = time
        dict['cur_time'] = time
        fill_x = np.ones(np.shape(dict['wave_time']))
        dict['Hs'] = np.ma.array(fill_x, mask=np.ones(np.shape(dict['wave_time'])))
        dict['name'] = 'AWAC8m'
        dict['xFRF'] = np.ma.array(np.ones(1), mask=np.ones(1))
        fill_ind = np.zeros(np.shape(dict['wave_time']))
        fill_ind[0] = 1
        dict['plot_ind'] = fill_ind
        dict['plot_ind_V'] = fill_ind
        dict['U'] = np.ma.array(fill_x, mask=np.ones(np.shape(dict['wave_time'])))
        dict['V'] = np.ma.array(fill_x, mask=np.ones(np.shape(dict['wave_time'])))
        dict['TS_toggle'] = False

    return dict

def lidar_PlotData(time, THREDDS='FRF'):

    t1 = time[0] - DT.timedelta(days=0, hours=0, minutes=3)
    t2 = time[-1] + DT.timedelta(days=0, hours=0, minutes=3)

    frf_Data = getObs(t1, t2, THREDDS)

    try:
        dict = {}
        lidar_data_RU = frf_Data.getLidarRunup()
        dict['runup2perc'] = lidar_data_RU['totalWaterLevel']
        dict['runupTime'] = lidar_data_RU['time']
        dict['runupMean'] = np.nanmean(lidar_data_RU['elevation'], axis=1)

        lidar_data_WP = frf_Data.getLidarWaveProf()
        dict['waveTime'] = lidar_data_WP['time']
        dict['xFRF'] = lidar_data_WP['xFRF']
        dict['yFRF'] = lidar_data_WP['yFRF']
        dict['Hs'] = lidar_data_WP['waveHsTotal']
        dict['WL'] = lidar_data_WP['waterLevel']

    except:
        # just make it a masked array
        dict['runupTime'] = np.zeros(20)
        fill_x = np.ones(np.shape(dict['runupTime']))
        dict['runup2perc'] = np.ma.array(fill_x, mask=np.ones(np.shape(dict['runupTime'])))
        dict['runupMean'] = np.ma.array(fill_x, mask=np.ones(np.shape(dict['runupTime'])))
        dict['waveTime'] = np.zeros(20)
        dict['xFRF'] = np.ones(np.shape(dict['waveTime']))
        dict['yFRF'] = np.ones(np.shape(dict['waveTime']))
        fill_x = np.ones((np.shape(dict['waveTime'])[0], np.shape(dict['xFRF'])[0]))
        dict['Hs'] = np.ma.array(fill_x, mask=np.ones(np.shape(fill_x)))
        dict['WL'] = np.ma.array(fill_x, mask=np.ones(np.shape(fill_x)))
        dict['TS_toggle'] = False

    return dict

def CMSF_velData(cmsfDict, station, dThresh=None):

    """
    this is a little function I wrote that will do the heavy lifting of pulling the current data from a particular gage,
    finds the closest model node to that gage, time matches the data, and returns the variables that need to be
    handed to obsVmod_TS to make pretty plots.
    :param cmsfDict: keys (that are used...):
                    'time' - this needs to be in epochtime
                    'aveE' - average eastward velocity
                    'aveN' - average northward velocity
    :param station: this is the stationname that will get handed to getCurrents, a gagenumber would (should?) also work
    :return: dictionary with keys:
             'time' - epochtimes of the matched data
             'aveEobs' - time-matched observed eastward velocity
             'aveNobs' - time-matched observed northward velocity
             'aveEmod' - time-matched model eastward velocity
             'aveNmod' - time-matched model northward velocity
    """
    # initialize this class
    timeunits = 'seconds since 1970-01-01 00:00:00'
    modTime = nc.num2date(cmsfDict['time'], timeunits)
    go = getObs(modTime[0] - DT.timedelta(minutes=3), modTime[-1] + DT.timedelta(minutes=3))

    # get my obs_dict
    obsV = go.getCurrents(station)
    aveUobs = obsV['aveU']
    aveVobs = obsV['aveV']
    obsTime = obsV['time']
    xFRFobs = obsV['xFRF']
    yFRFobs = obsV['yFRF']

    # find the closest node and pull that data
    ind, dist = findNearestUnstructNode(xFRFobs, yFRFobs, cmsfDict)
    if dThresh is None:
        pass
    else:
        assert dist <= dThresh, 'Error: this grid has no nodes within %s of gage %s.' %(dThresh, station)

    modTime = cmsfDict['time']
    aveUmod = cmsfDict['aveE'][:][ind]
    aveVmod = cmsfDict['aveN'][:][ind]

    # run the time matching.
    out = {}
    out['time'], out['aveEobs'], out['aveEmod'] = sb.timeMatch(obsTime, aveUobs, modTime, aveUmod)
    time, out['aveNobs'], out['aveNmod'] = sb.timeMatch(obsTime, aveVobs, modTime, aveVmod)

    # make it output a datetime
    modTime = [nc.num2date(ii, timeunits) for ii in out['time']]
    del out['time']
    out['time'] = modTime

    return out

def CMSF_wlData(cmsfDict, station, dThresh=None):
    """
    this is a little function I wrote that will do the heavy lifting of pulling the wl data from a particular gage,
    finds the closest model node to that gage, time matches the data, and returns the variables that need to be
    handed to obsVmod_TS to make pretty plots.
    :param cmsfDict: keys (that are used...):
                    'time' - this needs to be in epochtime
                    'waterLevel' - water level from the CSMF model
    :param station: this is the stationname that will get handed to getGageWL, a gagenumber would (should?) also work
    :return: dictionary with keys:
             'time' - epochtimes of the matched data
             'obsWL' - time-matched observed eastward velocity
             'modWL' - time-matched model northward velocity
    """
    # initialize this class
    timeunits = 'seconds since 1970-01-01 00:00:00'
    modTime = nc.num2date(cmsfDict['time'], timeunits)
    go = getObs(modTime[0] - DT.timedelta(minutes=3), modTime[-1] + DT.timedelta(minutes=3))

    # get my obs_dict
    obsWLdict = go.getGageWL(station)
    obsWL = obsWLdict['wl']
    obsTime = obsWLdict['time']
    xFRFobs = obsWLdict['xFRF']
    yFRFobs = obsWLdict['yFRF']

    # find the closest node and pull that data
    ind, dist = findNearestUnstructNode(xFRFobs, yFRFobs, cmsfDict)
    if dThresh is None:
        pass
    else:
        assert dist <= dThresh, 'Error: this grid has no nodes within %s of gage %s.' %(dThresh, station)

    modTime = cmsfDict['time']
    modWL = cmsfDict['waterLevel'][:][ind]

    # run the time matching.
    out = {}
    out['time'], out['obsWL'], out['modWL'] = sb.timeMatch(obsTime, obsWL, modTime, modWL)

    # make it output a datetime
    modTime = [nc.num2date(ii, timeunits) for ii in out['time']]
    del out['time']
    out['time'] = modTime

    return out

