import os
import io
import pandas as pd
from datetime import datetime
import pytz
import numpy as np

import wind_ae.McAstro.utils.constants as const
import wind_ae.McAstro.utils.timedate as timedate
from wind_ae.McAstro.utils.data_download import request_content

_lisird_directory = (os.path.dirname(os.path.abspath(__file__))
                     +'/lisird_spectrums/')

# LISIRD website and LATIS API
_lisird_url = 'http://lasp.colorado.edu/lisird/'
_latis_url = _lisird_url+'latis/dap/'

# Dailies
_maven_url = 'mvn_euv_l3_daily'
_see_url = 'timed_see_ssi_l3'
_eve_url = 'sdo_eve_ssi_1nm_l3'
_fism_url = 'fism_p_ssi_earth'
_fism2_url = 'fism_daily_hr'

_maven_date_func = timedate.UT_to_str
_see_date_func = timedate.YYYYDDD_to_str
_eve_date_func = timedate.YYYYDDD_to_str
_fism_date_func = timedate.UT_to_str
_fism2_date_func = timedate.YYYYDDD_to_str

_col_type1 = ['date', 'wl', 'F_wl', 'model_unc', 'total_unc']
_col_type2 = ['date', 'wl', 'F_wl', 'unc']
_maven_cols = _col_type1
_see_cols = _col_type2
_eve_cols = _col_type2
_fism_cols = _col_type1
_fism2_cols = _col_type2

_missions = {'maven':(_maven_url, _maven_cols, _maven_date_func),
             'see':(_see_url, _see_cols, _see_date_func),
             'eve':(_eve_url, _eve_cols, _eve_date_func),
             'fism':(_fism_url, _fism_cols, _fism_date_func),
             'fism2':(_fism2_url, _fism2_cols, _fism2_date_func)}

class lisird_spectrum:
    def __init__(self, mission='fism2', date='2009-07-19', sort_values='nu',
                 update=False, single_date=True):
        """
        Description:
            A class relating to LISIRD spectrums. The class will fetch,
            save, and load the lisird spectrums. Preps the LISIRD data into
            a standardized pandas DataFrame.
            
        Keyword arguments:
            mission: Which mission's data to fetch
            date: Date of observations (string format: %Y-%m-%d)
            sort_values: Which value dataframe is sorted by
            update: If data should be refetched if already exist (boolean)
            single_date: Only grab a single day of data (boolean)
        """
        self.mission = mission
        self.date = date
        if self.mission not in _missions.keys():
            print('Unknown mission, requires adding backend.')
            return
        date = list(map(int, date.split('-')))
        mission_filename = _lisird_directory+mission+'/'+mission
        if single_date:
            filters = [f'&time>={date[0]:04d}-{date[1]:02d}-{date[2]:02d}T',
                       f'&time<{date[0]:04d}-{date[1]:02d}-{(date[2]+1):02d}T']
            mission_filename += f'_{date[0]:04d}-{date[1]:02d}-{date[2]:02d}'
        elif mission == 'fism2': # Too big to be gotten in one go
            decade = date[0]-date[0]%10
            filters = ['&time>={}-01-01T'.format(decade),
                       '&time<{}-01-01T'.format(decade+10)]
            mission_filename += '_'+str(decade)
        else:
            filters = None
        mission_filename += '.parquet'
        mission_url, mission_cols, date_func = _missions[mission]
        if not os.path.exists(_lisird_directory+mission+'/'):
            os.makedirs(_lisird_directory+mission+'/', exist_ok=True)
        if update or not os.path.isfile(mission_filename):
            url = _latis_url+mission_url+'.csv'
            if filters is not None:
                url += '?'
                for f in filters:
                    url += f
            mission_csv = request_content(url)
            self.data = pd.read_csv(io.BytesIO(mission_csv))
            # Fix LISIRD column names
            self.data.columns = mission_cols
            # change nm to cm
            self.data['wl'] = self.data['wl']*1e-7
            self.data['nu'] = const.c/self.data['wl']
            # change seconds since unix epoch to julian date
            self.data['date'] = (
                date_func(self.data['date'].values, fmt='%Y-%m-%d')
            )
            self.data['date'] = pd.to_datetime(self.data['date'],
                                               format='%Y-%m-%d')
            # change W/m^2/nm to erg/s/cm^2/cm
            #b/c Lisird spectra are in 1 angstrom bins 
            self.data['F_wl'] = self.data['F_wl']*1e10 
            self.data['F_nu'] = self.data['F_wl']*self.data['wl']**2/const.c
            self.data.to_parquet(mission_filename, index=False)
        else:
            self.data = pd.read_parquet(mission_filename)
        self.date_min = self.data.iloc[0]['date']
        self.date_max = self.data.iloc[-1]['date']
        date = datetime(*date)
        date = date.replace(tzinfo=pytz.UTC)
        if date < self.date_min or date > self.date_max:#.tz_localize(None):
            print(f'Date given, {date.strftime("%Y-%m-%d")}, outside of '
                  f'{self.mission} date range.\n'
                  f'  Data date range: [{self.date_min.strftime("%Y-%m-%d")}, '
                  f'{self.date_max.strftime("%Y-%m-%d")}].')
            return
#         print(date,self.data['date'])
#         print(np.shape(self.data))
        self.data = self.data.loc[self.data['date']==date]
#         print(np.shape(self.data))
        self.data.drop('date', axis=1, inplace=True)
        self.data.reset_index(drop=True, inplace=True)
        # Drop bad data points
        self.data = self.data.loc[self.data['F_wl'] > 0]
        self.data = self.data.sort_values('wl')
        try:
            self.wl_min = self.data.iloc[0]['wl']
#             print("wl_min in lisard_spectrum",self.wl_min)
            self.wl_max = self.data.iloc[-1]['wl']
        except: # Empty DataFrame
            self.wl_min = 0
            self.wl_max = -1
        if sort_values != 'wl':
            self.data = self.data.sort_values(sort_values)