from pandas import pd
from datetime import datetime


def fetch_core_data():
    _fetch_atmoic_masses()
    # Currently ftp://gradj.pa.uky.edu//dima//photo//photo.dat
    # from http://www.pa.uky.edu/~verner/photo.html
    # is down so cannot automate fetching Verner tables

    
def _fetch_atmoic_masses():
    atomic_masses = pd.read_html('https://www.ciaaw.org/atomic-masses.htm',
                                 encoding='ISO-8859-1')[0]
    atomic_masses.columns = ['Z', 'Symbol', 'Element', 'A', 'mass']
    atomic_masses.drop(index=len(atomic_masses)-1, inplace=True)
    atomic_masses.drop(index=len(atomic_masses)-1, inplace=True)
    atomic_masses = ((atomic_masses['mass'].str.split('(', expand=True)
                     .rename(columns={0:'mass',1:'unc'}))
                     .combine_first(atomic_masses))
    atomic_masses['mass'] = atomic_masses['mass'].str.replace(u'\xa0', u'')
    atomic_masses.mass = atomic_masses.mass.astype(float)
    atomic_masses['unc'] = atomic_masses['unc'].str.replace(')','')
    atomic_masses['decays'] = atomic_masses.apply(lambda row: '*' in row['A'],
                                                  axis=1)
    atomic_masses['A'] = atomic_masses.apply(lambda row: row['A'].replace('*',''),
                                             axis=1)
    atomic_masses = atomic_masses[['Symbol', 'Element', 'decays',
                                   'A', 'Z', 'mass', 'unc']]
    # Write atomic_masses to file
    with open('atomic_masses.csv', 'w') as file_:
        file_.write('# mass in units of Daltons\n')
        file_.write('# Lasted pulled from '
                    'https://www.ciaaw.org/atomic-masses.htm on '
                    f'{datetime.today().strftime("%Y-%m-%d")}\n')
        atomic_masses.to_csv(file_, index=False)