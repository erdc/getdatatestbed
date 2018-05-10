from getDataFRF import getObs
from sblib.anglesLib import vectorRotation
import datetime as DT
import numpy as np
from sblib.gridTools import findNearestUnstructNode
import sblib.sblib as sb

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
        print('No data at %s!  Will return masked array.') %name
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
    write stuff here
    :param cmsfDict:
    :param station:
    :return:
    """
    # get my obs_dict
    obsV = getObs.getCurrents(station)
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

    return out


