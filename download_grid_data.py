# -------------------------------------------------------------------------------
# Name:        download_grid_data.py
# Purpose:     Downloads the most recent FRF GRID Survey data in TXT format.
# Author:      k5opjakk
# Created:     12/08/2015
# Edited:      Spicer Bak 12/25/15
#
# -------------------------------------------------------------------------------
import urllib
import urllib2
import json
import time
import gzip
import os
import sys
import numpy as np
import zipfile

__author__ = 'k5opjakk'


def download_survey(objectid, survey_filename, save_path):
    """
    Downloads the survey data extracts file to specified folder location
    :param objectid: objectid of survey to download
    :param survey_filename: filename of survey
    :param save_path: folder path to save file
    :return: path to downloaded survey file
    """
    task_url = "http://gis.sam.usace.army.mil/server/rest/services/FRF/DownloadSourceData/" \
               "GPServer/Download%20Source%20Data"
    submit_url = task_url + "/submitJob"
    try:
        urllib2.urlopen(submit_url)
    except:
        sys.exit('Service is down, please try again later')
    data = {'objectid': objectid, 'f': 'pjson'}
    submit_response = urllib.urlopen(submit_url, urllib.urlencode(data))
    submit_json = json.loads(submit_response.read())
    if 'jobId' in submit_json:
        job_id = submit_json['jobId']
        status = submit_json['jobStatus']
        job_url = task_url + "/jobs/" + job_id

        while status == "esriJobSubmitted" or status == "esriJobExecuting":
            print("checking to see if job is completed...")
            time.sleep(4)

            job_response = urllib.urlopen(job_url, "f=json")
            job_json = json.loads(job_response.read())

            if 'jobStatus' in job_json:
                status = job_json['jobStatus']

                if status == "esriJobSucceeded":
                        if 'results' in job_json:
                            results_url = job_url + "/results/"
                            results_json = job_json['results']
                            for param_name in results_json.keys():
                                result_url = results_url + param_name
                                result_response = urllib.urlopen(result_url, "f=json")
                                result_json = json.loads(result_response.read())
                                # print result_json
                                urllib.urlretrieve(result_json['value']['url'],
                                                   os.path.join(save_path, survey_filename + '.zip'))
                            z = zipfile.ZipFile(os.path.join(save_path, survey_filename + '.zip'), 'r')
                            gzip_file = z.namelist()[0]
                            z.extractall(save_path)
                            z.close()
                            in_file = gzip.GzipFile(os.path.join(save_path, gzip_file), 'rb')
                            s = in_file.read()
                            in_file.close()
                            out_file = file(os.path.join(save_path, survey_filename), 'wb')
                            out_file.write(s)
                            out_file.close()
                            os.remove(os.path.join(save_path, survey_filename + '.zip'))
                            os.remove(os.path.join(save_path, gzip_file))
                            print("File located here: {}".format(os.path.join(save_path, survey_filename)))
                            return os.path.join(save_path, survey_filename)

                if status == "esriJobFailed":
                        if 'messages' in job_json:
                            print("Job Failed")
                            raise SystemError('Esri Server Down, Contact Mobile Alabama 1-800-USACE-Mobile')
                            #  print job_json['messages']
    else:
        print "no jobId found in the response"


def query_survey_data(service_url, grid_data=True):
    """
    Returns the objectid, survey_filename of the most recent survey GRID
    :param service_url: url to the service to query
    :return: type: number, objectid; type: string, filename
    """
    try:
        urllib2.urlopen(service_url)
    except:
        sys.exit('Service is down, please try again later')

    rest_operation = {
        'query': service_url + '/query?',
        'identify': service_url + '/identify?'
    }

    parameters = urllib.urlencode({
        'where': """SURVEYTYPE = 'GRID'""",
        'geometry': '',
        'geometryType': '',
        'inSR': '',
        'spatialRel': '',
        'outFields': 'OBJECTID, SURVEYDATE, TEXTFILENAME',
        'returnIdsOnly': 'false',
        'returnGeometry': 'false',
        'orderByFields': 'SURVEYDATE DESC',
        'f': 'json'
    })
    query_url = urllib2.urlopen(rest_operation['query'] + parameters) # this was comma'd
    data = json.load(query_url)
    survey_filename, objectid, surveydate = [], [], []
    for i in range(0, len(data['features'])):
        survey_filename.append( data['features'][i]['attributes']['TEXTFILENAME'])
        objectid.append(data['features'][i]['attributes']['OBJECTID'])
        surveydate.append(data['features'][i]['attributes']['SURVEYDATE']/1000)
    # turn surveydate into numpy array
    surveydate = np.array(surveydate)

    return objectid, survey_filename, surveydate

if __name__ == '__main__':
    objectid, survey_filename, surveydate = query_survey_data(
        service_url='http://gis.sam.usace.army.mil/server/rest/services/FRF/FRF/FeatureServer/4')
    #output_file = download_survey(objectid=objectid, survey_filename=survey_filename, save_path=r'D:\temp')

    for file in survey_filename:
        print file + '\n'