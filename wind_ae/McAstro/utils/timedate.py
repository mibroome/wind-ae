import locale
from pytz import timezone as pytz
from datetime import datetime, timedelta

locale.setlocale(locale.LC_ALL, '')
_unix_epoch = datetime(1970, 1, 1, tzinfo=pytz('UTC'))

def UT_to_str(time, timezone='UTC', fmt='%Y-%m-%d %H:%M:%S %Z',
              epoch=_unix_epoch):
    """
    Arguments:
        time: time in seconds since epoch
    
    Keyword arguments:
        timezone: one of the pytz timezones (default: UTC)
        fmt: how time string is formatted (default: %Y-%m-%d %H:%M:%S %Z)
        epoch: when t=0 (default: unix, 1970-01-01 00:00:00 UTC)
        
    Returns:
        String of date
        
    Example:
        UT_to_str(345459600, timezone='America/New_York')
    """
    if type(time) is not float:
        date = []
        for t in time:
            unix_time = epoch + timedelta(seconds=float(t))
            date.append(unix_time.astimezone(tz=pytz(timezone)).strftime(fmt))
        return date
    unix_time = epoch + timedelta(seconds=time)
    return unix_time.astimezone(tz=pytz(timezone)).strftime(fmt)

def str_to_UT(time, timezone='UTC', fmt='%Y-%m-%d %H:%M:%S',
              epoch=_unix_epoch):
    """
    Arguments:
        time: string to be convert to epoch time
        
    Keyword argumens:
        timezone: one of the pytz timezones (default: UTC)
        fmt: how time string is formatted (default: %Y-%m-%d %H:%M:%S)
        epoch: when t=0 in UTC (default: 1970-01-01)
        
    Returns:
        Integer number of seconds since epoch
        
    Example:
        str_to_UT('2020-01-01 00:00:00', timezone='US/Pacific')
    """
    tz = pytz(timezone)
    time = tz.localize(datetime.strptime(time, fmt))
    return int((time-epoch).total_seconds())


def YYYYDDD_to_str(time, fmt='%Y-%m-%d'):
    if type(time) is not float:
        date = []
        for t in time:
            days = t%1000
            year = (t-days)//1000
            ordinal = datetime(year, 1, 1).toordinal()+days-1
            date.append(datetime.fromordinal(ordinal).strftime(fmt))
        return date
    days = t%1000
    year = (t-days)//1000
    ordinal = datetime(year, 1, 1).toordinal()+days-1
    return datetime.fromordinal(ordinal).strftime(fmt)