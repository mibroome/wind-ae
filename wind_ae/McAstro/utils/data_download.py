import time
import random
import requests

def request_content(url, **kwargs):
    """
    Description:
        Given a url, we will perform a request using the requests library.
        This handles error catching and sleeps to avoid spamming, in hopes
        of avoiding being blocked from future request.
        
    Arguement:
        url: url to perform a request on
        
    Keyword argument:
        kwargs: keyword arguments passed to request.get()
        
    Notes:
        ^Sleep [0.5, 1) seconds after a request to avoid spamming server
        
    Returns:
        A response's content
    """
    try:
        response = requests.get(url, **kwargs)
        response.raise_for_status()
    except requests.exceptions.HTTPError as error:
        raise SystemExit(error)
    except requests.exceptions.RequestException as error:
        raise SystemExit(error)
    print('[HTTP STATUS:{}] {} (elapsed={})'
          .format(response.status_code, url, response.elapsed))
    time.sleep(0.5*(1+random.random()))
    return response.content