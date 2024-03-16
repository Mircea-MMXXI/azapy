from .MkTreader import MkTreader

def readMkT(symbol=[], sdate="2012-01-01", edate='today', calendar=None,
            output_format='frame', source=None, force=False, save=True,
            file_dir="outDir", file_format='csv', api_key=None, param=None,  
            verbose=True):
    """
    Retrieves market data for a set of stock symbols.\n
    It is a wrapper for `MkTreader` class returning directly the requested
    historical time series. The function call variables are the same as for 
    'MkTreader' member function 'get'.
    
    Parameters
    ----------
    symbol : `str` or `list` of `str`, optional
        Stock symbols to be uploaded.
        The default is `[]`.
    sdate : date like, optional
        The start date of historical time series.
        The default is `"2012-01-01"`.
    edate : date like, optional
        The end date of historical time series (must: `sdate` >= `edate`)
        The default is `'today'`.
    calendar : `numpy.busdaycalendar`, optional
        Exchange business day calendar. If set to `None` it will default to 
        the NY stock exchange business calendar (provided by the azapy 
        function NYSEgen).
        The default is `None`.
    output_format : `str`, optional
        The function output format. It can be:
            - `'frame'` - `pandas.DataFrame`
            - `'dict'` - `dict` of `pandaws.DataFrame`. 
              The symbols are the keys.
              
        The default is `'frame'`
    source : `str` or `dict`, optional
        If it is a `str`, then it represents the market data provider for 
        all historical prices request. Possible values are: `'yahoo'`, 
        `'alphavantage'`, `'alphavantage_yahoo'`, `'eodhistoricaldata'`,
        `'eodhistoricaldata_yahoo'` and `'marketstack'`. If set to `None` 
        it will default to `'yahoo'`.\n
        It can be set to a `dict` containing specific instructions for 
        each stock symbol. 
        The `dict` keys are the symbols and the values are 'dict'
        instructions specific to each symbol. 
        Valid keys for the instructions `dict` are the names of
        this function call variables except `'sdate'`, `'edate'`, 
        `'calendar'` and `'output_format'`. 
        The actual set of stock symbols is given by the union 
        of variable `'symbol'` and the keys of the dict `'source'`. Missing  
        values in the symbol instruction dict's will be filled with the 
        values of the function call variables. 
        The values of the function call variables act as 
        generic values to be used in absence of specific instructions 
        in the `'source'` dict. The default is `None`.\n
        *Example* of dict `'source'`: \n
        source = {'AAPL': {'source': 'eodhistoricaldata, 'verbose': `True`},
        'SPY': {'source': 'yahoo', 'force': `True`}}\n
        In this case there are 2 symbols that will be added (union) to 
        the set of symbols defined by 'symbol' variable. For symbol 'AAPL' 
        the provider source is eodhistoricaldata and the 'verbose' 
        instruction 
        is set to `True`. The rest of the instructions: 'force', 'save',
        'file_dir', 'file_format', 'api_key' and 'param' are set 
        to the values of the corresponding function call variables.
        Similar for symbol 'SPY'. The instructions for the rest of the 
        symbols that may be specified in the 'symbol' variable will be
        set according to the values of the function call variables.
    force : `Boolean`, optional
            - `True`: will try to collect historical prices exclusive from
              the market data providers.
            - `False`: first it will try to load the historical 
              prices from a local saved file. If such a file does not exist
              the market data provider will be accessed.  
              
        If the file exists but the saved historical 
        data is too short then it will try to collect the missing values 
        only from the market data provider.
        The default is `False`.
    save : `Boolean`, optional
            - `True`: It will try to save the historical price collected from 
              the providers to a local file.
            - `False`: No attempt to save the data is made.
            
        The default is `True`.
    file_dir : `str`, optional
        Directory with (to save) historical market data. If it does not 
        exists then it will be created.
        The default is "outDir".
    file_format : `str`, optional
        The saved file format for the historical prices. The following 
        files formats are supported: csv, json and feather
        The default is 'csv'.
    api_key : `str`, optional
        Provider API key (where is required). If set to `None`  
        then the API key is set to the value of global environment variables\n 
            - `APLPHAVANTAGE_API_KEY` for alphavantage,
            - `EODHISTORICALDATA_API_KEY` for eodhistoricaldata,
            - `MARKETSTACK_API_KEY` for marketstack.
        
        The default is `None`.
    param : `dict`, optional
        Set of additional information to access the market data provider.  
        At this point in time only accessing alphavantage provider requires 
        an additional parameter specifying the maximum number of API 
        (symbols) requested per minute. 
        It varies with the level of access 
        corresponding to the API key. The minimum value is 5 for a free key 
        and starts at 75 for premium keys. This value is stored in
        `max_req_per_min` variable.\n
        *Example*: param = {'max_req_per_min': 5}\n
        This is also the default vale for alphavantage, if `param` is set to 
        `None`.
        The default is `None`.   
    verbose : `Boolean`, optional
        If set to `True`, then additional information will be printed  
        during the loading of historical prices.
        The default is `True`.


    Returns
    -------
    `pandas.DataFrame` or 'dict' `pandas.DataFrame` : Historical market data. 
        The output format is designated by the value of the input parameter 
        `output_format`.
    """
    return MkTreader().get(symbol=symbol, sdate=sdate, edate=edate, 
                           calendar=calendar, output_format=output_format, 
                           source=source, force=force, save=save, 
                           file_dir=file_dir, file_format=file_format, 
                           api_key=api_key, param=param,  verbose=verbose)