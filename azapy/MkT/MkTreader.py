import pandas as pd
import os
import yfinance as yf
import requests
import numpy as np
import pathlib
import time
import logging
from copy import deepcopy

from azapy.MkT.MkTcalendar import NYSEgen


class MkTreader:
    """ 
    Collects historical market prices from market data providers such as
    'yahoo', 'eodhistoricaldata', 'alphavantage' and 'marketstack'.
    
    **Attributs**
        - `dsource` : dict of request instructions per symbol
        - `delta_time` : execution time of the request in seconds
        - `rout` : pandas.DataFrame containing historical prices for all
          symbols. It is created during the call of `get` function.
        - `rout_status` : request status information. It is created during 
          the call of `get_request_status` function or during the 
          call of function `get` with option `verbose=True`.
        - `error_log` : contains lists of missing historical observation 
          dates. It is created together with `rout_status`.
    """
    def __init__(self, verbose=True):
        '''
        Constructor
        
        Parameters
        ----------
        verbose : `Boolean`, optional
            If set to `True`, additional information will be printed  during 
            the loading of historical prices. The default value is `True`.

        Returns
        -------
        The MkTreder object
        '''
        self._verbose = False if verbose is None else verbose
        self._set_verbose(self._verbose)
            
        self.dsource = None
        self.delta_time = None
        self.rout = None
        self.rout_status = None
        self.error_log = None
        
        self.sdate = None
        self.edate = None
        
        self._bday = None
        
        self._col = ['open', 'high', 'low', 'close', 'volume', 'adjusted', 
                     'divd', 'split']
        self._out_col = ['symbol'] + self._col + ['source', 'recordDate']
        self._alphavantage_max_req_per_min = 5
 
        
    def get(self, symbol=[], sdate="2012-01-01", edate='today', calendar=None,
            output_format='frame', source=None, force=False, save=True,
            file_dir="outDir", file_format='csv', api_key=None, param=None,  
            verbose=None):
        """
        Retrieves market data for a set of stock symbols.
        
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
                - `'dict'` - `dict` of `pandas.DataFrame`. 
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
            If set `True`, the additional information will be printed  
            during the loading of historical prices. If `None` it is 
            ignored, otherwise it overwrites the value set by the constructor.
            The default value is `None`.

        Returns
        -------
        `pandas.DataFrame` or 'dict' `pandas.DataFrame` : Historical market data. 
            The output format is designated by the value of the input parameter 
            `output_format`.
        """
        # Process the inputs
        self._set_verbose(verbose)
            
        if (calendar is None) or (calendar == 'NYSE'):
            calendar = NYSEgen()
        if isinstance(calendar, np.busdaycalendar):
            self._bday = pd.offsets.CustomBusinessDay(calendar=calendar)
        elif isinstance(calendar, pd.offsets.CustomBusinessDay):
            self._bday = calendar
        else:
            raise ValueError(f"Wrong calendar type {type(calendar)} "
                             "must be numpy.busdaycalendar")
            
        self.sdate, self.edate = self._mkt_sedate(sdate, edate)
        if self.sdate > self.edate:
            raise ValueError("Wrong rage of dates -"
                             f" start date {self.sdate} > end date {self.edate}")
 
        outputFormats = ('frame', 'dict')
        if output_format not in outputFormats:
            raise ValueError( f"Wrong output format {output_format} -"
                             f" must be one of {outputFormats}")
        #self._output_format = output_format
        
        fileFormats = ('csv', 'json', 'feather')
        if file_format not in fileFormats:
            raise ValueError(f"Wrong file format {file_format} -"
                             f" must be one of {fileFormats}")
        #self._file_format = file_format

        if isinstance(symbol, str):
            symbol = [symbol]
        elif not isinstance(symbol, list):
            raise ValueError(f"Wrong symbol type: {type(symbol)} "
                             "must be str or a list of str")
        elif not all([isinstance(sy, str) for sy in symbol]):
            raise ValueError("Wrong element type in the symbol list"
                             " - all must be str")
        
        allSources = ('yahoo', 'alphavantage', 'eodhistoricaldata', 'marketstack',
                      'alphavantage_yahoo', 'eodhistoricaldata_yahoo')
        if (source is None) or isinstance(source, dict):
            lsource = 'yahoo'
        elif source in allSources:
            lsource = source
        else:
            raise ValueError(f"Wrong source type {type(source)}: "
                             "must be None, str, or dict")

        sy_inf = {'source': lsource,
                  'force': force,
                  'save': save,
                  'file_dir': file_dir,
                  'file_format': file_format,
                  'api_key': api_key,
                  'param': param}
        self.dsource = {sy: deepcopy(sy_inf) for sy in symbol}
        if isinstance(source, dict):
            self.dsource = source.update(self.dsource)

        dft = pd.DataFrame.from_dict(self.dsource, orient='index')
        dft.index.name = 'symbol'
        dft.reset_index(inplace=True)

        # prep computation
        def _srr2(nn):
            if (dft.loc[nn,'source'] == 'alphavantage') or \
               (dft.loc[nn,'source'] == 'alphavantage_yahoo'):
                return (dft.loc[nn, 'api_key'], 'alpha')
        
            return (dft.loc[nn, 'api_key'], 'regular')
 
        dfg = dft.groupby(_srr2)
        
        # main computation loop
        tic = time.perf_counter()
        rout = []
        for nm, gr in dfg: 
            match nm[1]:
                case 'alpha':
                    rout.append(self._alphavantage_process(gr))
                case 'regular':
                    rout.append(self._regular_process(gr))
                case _:
                    # you should not be here - ever
                    raise ValueError(f"Unknown category {nm}\n {gr}")

        if len(rout) == 0:
            if self._verbose:
                print("Warning: no mkt data was fund!")
            return {} if output_format == 'dict' else pd.DataFrame()
        
        self.rout = pd.concat(rout)
        toc = time.perf_counter()
        self.delta_time = toc - tic

        # output
        if self._verbose:
            self.get_request_status()
            with pd.option_context("display.max_columns", None):
                print("\nRequest between "
                      f"{self.sdate.date()} : {self.edate.date()}\n"
                      f"{self.rout_status}\n"
                      f"extraction time {round(self.delta_time, 3)} s")
            
        if output_format == 'dict':
            if self.rout.empty: return {}
            return dict(tuple(self.rout.groupby('symbol')))
    
        return self.rout
    
    
    def _mkt_sedate(self, sdate, edate):
        dn = pd.Timestamp.now(tz='America/New_York')
        refdate = dn.replace(hour=16, minute=0, second=0, microsecond=0)
        if dn > refdate:
            refdate = self._bday.rollback(refdate).normalize().tz_localize(None)
        else:
            refdate = self._bday.rollback(refdate - pd.Timedelta(1, 'day')).normalize().tz_localize(None)

        edate = self._bday.rollback(min(pd.to_datetime(edate).normalize().tz_localize(None), refdate))

        if sdate == '': 
            sdate = edate
        elif isinstance(sdate, int) and (sdate <= 0):
            sdate = edate + sdate * self._bday
        else:
            sdate = self._bday.rollforward(sdate).normalize().tz_localize(None)

        return sdate, edate


    def _set_verbose(self, verbose):
        if verbose is None:
            return
        self._verbose = verbose
        if self._verbose:
            logging.getLogger('yfinance').setLevel(logging.INFO)
        else:
            logging.getLogger('yfinance').setLevel(logging.CRITICAL)

    
    def get_request_status(self, verbose=None):
        '''
        Reports abbreviated information about request status.
        
        verbose : `Boolean`, optional
            If set to `True`, additional information will be printed  
            during the function execution. In set to 'None', it will be 
            ignored, otherwise it will overwrite the value set by the 
            constructor.
            The default value is `None`.

        Returns
        -------
        `pandas.DataFrame` : The status report.
        The column names are the symbols for which the data was requested.
        The rows contain the actual input parameters per symbol as well
        as:   
            - `'nrow'` : the length of historical time series.
            - `'sdate'` : first date in the time series.
            - `'edate'` : end date of the time series.
            - `'error'` : if there are missing data. If its value is `'Yes'`
              then the actual list of missing date per symbol can
              be obtained by calling `get_error_log`.
        '''
        self._set_verbose(verbose)
            
        if self.rout is None or self.rout.empty:
            if self._verbose:
                print("Warning: request was returned empty")
            self.rout_status = pd.DataFrame(self.dsource)
            return self.rout_status
        
        self.rout_status = deepcopy(self.dsource)
        for kk in self.rout_status.keys():
            sblock = self.rout.loc[self.rout['symbol'] == kk]
            if sblock.empty: continue
            if not self.rout_status[kk]['param'] is None: 
                for pkk in self.rout_status[kk]['param'].keys():
                    self.rout_status[kk]['param:' + pkk] = \
                         self.rout_status[kk]['param'][pkk]
            del self.rout_status[kk]['param']
            self.rout_status[kk]['nrow']  = sblock.shape[0]
            self.rout_status[kk]['sdate'] = sblock.index[0].date()
            self.rout_status[kk]['edate'] = sblock.index[-1].date()
            self._ts_analyzer()
            
        self.rout_status = pd.DataFrame(self.rout_status)
        return self.rout_status
      
    
    def get_error_log(self):
        '''
        Returns lists of missing historical observation dates per
        symbol

        Returns
        -------
        `dict` : The error-log.
        
        If it is an empty `dict` then there are no missing dates in the
        collected historical time series.
        Otherwise, the keys of the `dict` are the symbols that have missing 
        dates. The values for these keys are also `dict` with the following 
        fields: 
            - `'back'`: a list of missing date at the tail of the time series
            - `'front'` : a list of missing data at the head of the time series
            - `'mid'` : a list of missing data in the middle of the time series
        
        Fields with empty list of dates are omitted.
        '''
        return self.error_log
      
    
    def _ts_analyzer(self):
        self.error_log = {}
        lbday = pd.bdate_range(self.sdate, self.edate, freq='C', 
                               holidays=self._bday.holidays)
        for symbol in self.dsource.keys():
            symbData = self.rout.loc[self.rout['symbol'] == symbol]
            errdif = lbday.difference(symbData.index)
            if errdif.empty: 
                self.rout_status[symbol]['error'] = 'No'
            elif symbData.empty:
                self.rout_status[symbol]['error'] = 'Wrong symbol'
                self.error_log[symbol] = {}
                self.error_log[symbol]['mid'] = errdif
            else:
                self.rout_status[symbol]['error'] = 'Yes'
                self.error_log[symbol] = {}
                eerr = errdif[errdif < symbData.index[0]]
                if not eerr.empty:
                    self.error_log[symbol]['back'] = eerr
                eerr = errdif[errdif > symbData.index[-1]]
                if not eerr.empty:
                    self.error_log[symbol]['front'] = eerr
                eerr = errdif[(errdif >= symbData.index[0]) & 
                              (errdif <= symbData.index[-1])]
                if not eerr.empty:
                    self.error_log[symbol]['mid'] = eerr
    
    
    def _alphavantage_process(self, data):
        counter = 0
        timer = 0.
        #min req per min is 5
        max_rqs = max([self._alphavantage_max_req_per_min] + 
                [x['max_req_per_min'] for x in data['param'] if not x is None])
        
        rout = []
        for row in data.index:
            tic = time.perf_counter()
            py = data.iloc[row]
            kod, pri = self._reader_symb(py.symbol, py.source, 
                                         py.force, py.save, 
                                         py.file_dir, py.file_format,
                                         py.api_key, py.verbose)

            if kod:
                timer += (time.perf_counter() - tic)
                counter += 1
                
            if pri.empty: continue
 
            pri = pri.iloc[[self._bday.is_on_offset(x) for x in pri.index]]
            rout.append(pri)
            
            if counter >= max_rqs:
                if timer < 60.:
                    if self._verbose: 
                        print(f"time to sleep {60. - timer} s")
                        
                    time.sleep(60. - timer)
                timer = 0.
                counter = 0
                
        if len(rout) == 0:
            return pd.DataFrame()
        else:   
            return pd.concat(rout)


    def _regular_process(self, data):
        rout = []
        for row in data.index:
            py = data.loc[row]
            kod, pri = self._reader_symb(py.symbol, py.source,
                                         py.force, py.save, 
                                         py.file_dir, py.file_format,
                                         py.api_key)
            if pri.empty: continue
            pri = pri.iloc[[self._bday.is_on_offset(x) for x in pri.index]]
            rout.append(pri)

        if len(rout) == 0:
            return pd.DataFrame()
        else:   
            return pd.concat(rout) 
    
    
    def _reader_symb(self, symbol, source, 
                     flag_web_only=True, flag_save=False,
                     file_dir="./", file_format='csv',
                     api_key=None):
        kod_web = False

        if flag_web_only:
            kod_web = True
            if self._verbose:
                print(f"get {symbol} from {source} only")
                
            rout_web = self._reader_web(source, symbol, api_key=api_key)
            if rout_web.empty:
                return kod_web, pd.DataFrame()
            
            if flag_save:
                if self._verbose: 
                    print(f"save {symbol} data to file")
                rout_web['adjusetd'] = self._adjustDividend(rout_web)
                self._writer_disk(rout_web, symbol, file_dir, file_format)

            rout = rout_web.loc[self.sdate : self.edate].copy()
            if rout.empty:
                if self._verbose:
                    print("no data in the range "
                          f"{self.sdate.date()} : {self.edate.date()} "
                          f"for {symbol} from source {source}")
                return kod_web, pd.DataFrame()
            
            rout['adjusted'] = self._adjustDividend(rout)
            return kod_web, rout[['symbol'] + self._col]
        else:
            rout_disk = self._reader_disk(symbol, file_dir, file_format)
            if rout_disk.empty:
                if self._verbose:
                    print(f"no saved data for {symbol}")
                rout = pd.DataFrame()
                sdate_web = pd.to_datetime("1900-01-01")
            else:
                if self._verbose:
                    print(f"read {symbol} data from file")
                    
                rout = rout_disk.loc[:self.edate].copy()
                rout['adjusted'] = self._adjustDividend(rout)   
                sdate_web = rout_disk.index[-1] + self._bday
                
            if sdate_web <= self.edate:
                kod_web = True
                if self._verbose:
                    print(f"get {symbol} updates from {source}")
                    
                rout_web = self._reader_web(source, symbol, 
                                            [sdate_web, self.edate], api_key)
                if rout_web.empty:
                    if self._verbose:
                        print(f"{symbol}:{source} no data in range "
                              f"{sdate_web.date()}:{self.edate.date()}") 
                    return False, pd.DataFrame()
                else:
                    rout = pd.concat([rout, rout_web])
                    
                    rout['adjusted'] = self._adjustDividend(rout)    
                    if flag_save:
                        if self._verbose:
                            print(f"save {symbol} updated data to file")
                            
                        self._writer_disk(rout[self._out_col], symbol, 
                                          file_dir, file_format)
                
            return kod_web, rout.loc[self.sdate: , ['symbol'] + self._col]
      
        
    def _reader_web(self, source, symbol, drange=None, api_key=None):
        if source == 'yahoo':
            price = self._yahoo_finance(symbol)
        elif source == 'eodhistoricaldata':
            price = self._eodhistoricaldata(symbol, api_key)
        elif source == 'eodhistoricaldata_yahoo':
            price = self._eodhistoricaldata_yahoo(symbol, api_key)
        elif source == 'alphavantage':
            price = self._alphavantage(symbol, api_key)
        elif source == 'alphavantage_yahoo':
            price = self._alphavantage_yahoo(symbol, api_key)
        elif source == 'marketstack':
            price = self._marketstack(symbol, api_key)
        else:
            raise ValueError(f"Unknown source {source}")
        
        if price.empty:
            return pd.DataFrame()
        
        price['symbol'] = symbol
        price['source'] = source
        price['recordDate'] = pd.to_datetime('today').strftime("%Y-%m-%d %X")
        price.index = pd.to_datetime(price.index.date, utc=False)
        price.index.name = 'date'

        if drange is None:
            return price[self._out_col]
        else:
            return price.loc[drange[0]:drange[1], self._out_col]
    

    def _reader_disk(self, symbol, file_dir, file_format):
        file_name = file_dir + "/" + symbol + "." + file_format
        if os.path.isfile(file_name):
            if file_format == 'csv':
                price = pd.read_csv(file_name, index_col='date')
                price.index = pd.to_datetime(price.index)
                return price
            elif file_format == 'feather':
                price = pd.read_feather(file_name).set_index('date')
                return price
            elif file_format == 'json':
               price = pd.read_json(file_name).set_index('date')
               return price
            else:
               raise ValueError(f"Unknown file format {file_format} "
                                "supported are: csv, ft and json")
        
        return pd.DataFrame()


    def _writer_disk(self, data, symbol, file_dir, file_format):
        pathlib.Path(file_dir).mkdir(parents=True, exist_ok=True) 
        
        file_name = file_dir + '/' + symbol + '.' + file_format
        if file_format == 'csv':
            data[self._out_col].to_csv(file_name, index=True)
        elif file_format == 'feather':
            data.reset_index()[['date'] + self._out_col].to_feather(file_name)
        elif file_format == 'json':
            data.reset_index()[['date'] + self._out_col].to_json(file_name)
        else:
            raise ValueError(f"unknown output file format {file_format} "
                             "supported are: csv, ft and json")
    
            
    def _yahoo_finance(self, symbol):
        '''returns maximum period adjusted for split only'''
        
        ysymb = yf.Ticker(symbol)
        yrprice = ysymb.history(period="max", 
                                auto_adjust=False, rounding=False)
        if yrprice.empty:
            return pd.DataFrame()
    
        yrprice.rename(columns={'Open': 'open', 'High': 'high',
                                'Low': 'low', 'Close': 'close',
                                'Adj Close': 'adjusted', 'Volume': 'volume',
                                'Dividends': 'divd', 'Stock Splits': 'split'},
                       inplace=True)
        yrprice.index.name = 'date'
        yrprice.loc[yrprice['split'] == 0., 'split'] = 1.
    
        return yrprice[self._col]
    
    
    def _eodhistoricaldata(self, symbol, api_key=None):
        '''
        returns maximum period adjusted for splits only,
        test key OeAFFmMliFG5orCUuwAKQ8l4WWFQ67YX for MCA only,
        needs premium key for anything else.
        '''
        if api_key is None:
            api_key = os.getenv('EODHISTORICALDATA_API_KEY') 
            if api_key is None:
                raise ValueError("Worning: no API key set as "
                                 "global environment variable")
                return pd.DataFrame()
                
        eodh = ('https://eodhistoricaldata.com/api/eod/' 
                + symbol + '.US?api_token=' + api_key
                + '&order=a')
        try:
            eeprice = pd.read_csv(eodh)[:-1]
        except:
            if self._verbose:
                print(f"Warning: {symbol} "
                      "no price data in eodhistoricaldata")
            return pd.DataFrame()
        
        eeprice.rename(columns={'Date': 'date', 
                                'Open': 'open', 'High': 'high',
                                'Low': 'low', 'Close': 'close',
                                'Adjusted_close': 'adjusted', 
                                'Volume': 'volume'},
                       inplace=True)
        eeprice.set_index('date', inplace=True)
        eeprice.index = pd.to_datetime(eeprice.index)
        
        eodd = ('https://eodhistoricaldata.com/api/div/' 
                + symbol + '.US?api_token=' + api_key)
        try:
            edprice = pd.read_csv(eodd)[:-1]
        except:
            if self._verbose:
                print(f"Warning: {symbol} "
                      "no dividend data in eodhistoricaldata")
            return pd.DataFrame()
    
        if edprice.empty:
            eeprice['divd'] = 0.
        else:
            edprice.rename(columns={'Date': 'date', 'Dividends': 'divd'},
                           inplace=True)
            edprice.set_index('date', inplace=True)
            edprice.index = pd.to_datetime(edprice.index)
            eeprice = pd.merge(eeprice, edprice, how='left', 
                               left_index=True, right_index=True)
            eeprice['divd'].fillna(0., inplace=True)
        
        eods = ('https://eodhistoricaldata.com/api/splits/' 
                + symbol + '.US?api_token=' + api_key) 
        try:
            esprice = pd.read_csv(eods)[:-1]
        except:
            if self._verbose:
                print(f"Warning: {symbol} "
                      "no split data in eodhistoricaldata")
            return pd.DataFrame()
    
        if esprice.empty:
            eeprice['split'] = 1.
        else:
            esprice.rename(columns={'Date': 'date', 'Stock Splits': 'split'},
                           inplace=True)
            esprice.set_index('date', inplace=True)
            esprice = esprice.apply(lambda x: [eval(y) for y in x])
            esprice.index = pd.to_datetime(esprice.index)
            eeprice = pd.merge(eeprice, esprice, how='left', 
                               left_index=True, right_index=True)
            eeprice['split'].fillna(1., inplace=True)
            
        self._adjustSplit(eeprice, volume=False, inplace=True)
            
        return eeprice[self._col]
    
    
    def _eodhistoricaldata_yahoo(self, symbol, api_key=None):
        '''
        returns maximum period adjusted for splits only,
        with prices from edohistoricaldata and dividends + splits from yahoo
        test key OeAFFmMliFG5orCUuwAKQ8l4WWFQ67YX for MCA only,
        can work with a free key but only for 1y!!
        '''
        if api_key is None:
            api_key = os.getenv('EODHISTORICALDATA_API_KEY') 
            if api_key is None:
                raise ValueError("Worning: no API key set as "
                                 "global environment variable")
                
        eodh = ('https://eodhistoricaldata.com/api/eod/' 
                + symbol + '.US?api_token=' + api_key
                + '&order=a')
        try:
            eeprice = pd.read_csv(eodh)[:-1]
        except:
            if self._verbose:
                print(f"Warning: {symbol} "
                      "no price data in eodhistoricaldata")
            return pd.DataFrame()
        
        eeprice.rename(columns={'Date': 'date',
                                'Open': 'open', 'High': 'high',
                                'Low': 'low', 'Close': 'close',
                                'Adjusted_close': 'adjusted', 
                                'Volume': 'volume'},
                       inplace=True)
        eeprice.set_index('date', inplace=True)
        eeprice.index = pd.to_datetime(eeprice.index)
     
        ysymb = yf.Ticker(symbol)
        
        edprice = pd.DataFrame({'divd': ysymb.dividends})
        if edprice.empty:
            eeprice['divd'] = 0.
        else:
            edprice.index = pd.to_datetime(edprice.index.date, utc=False)
            edprice.index.name = 'date'
            eeprice = pd.merge(eeprice, edprice, how='left', 
                               left_index=True, right_index=True)
            eeprice['divd'].fillna(0., inplace=True)
            
        esprice = pd.DataFrame({'split': ysymb.splits})
        if esprice.empty:
            eeprice['split'] = 1.
        else:
            esprice.index = pd.to_datetime(esprice.index.date, utc=False)
            esprice.index.name = 'date'
            eeprice = pd.merge(eeprice, esprice, how='left', 
                               left_index=True, right_index=True)
            eeprice['split'].fillna(1., inplace=True)
            
        self._adjustSplit(eeprice, volume=False, inplace=True)
            
        return eeprice[self._col]
    
    
    def _alphavantage(self, symbol, api_key=None):
        ''' 
        returns maximum period adjusted for splits only,
        with data from alphavantage
        key "demo" works only for"IBM"
        for anything else a premium key is requierd
        '''
        if api_key is None:
            api_key = os.getenv('ALPHAVANTAGE_API_KEY')
            if api_key is None:
                raise ValueError("Warning: no api_key is set as a "
                                 "global environment variable")    
                return pd.DataFrame()
                
        url = ('https://www.alphavantage.co/query?'
               + 'function=TIME_SERIES_DAILY_ADJUSTED'
               + '&symbol=' + symbol 
               + '&outputsize=full' 
               + '&apikey=' + api_key)
    
        req = requests.get(url)
        if req.status_code != 200:
            if self._verbose:
                print(f"request status code = {req.status_code} no data"
                      f" collected for {symbol} from alphavantage")
            return pd.DataFrame()
        
        rout = req.json()
        if len(rout.keys()) == 1:
            if self._verbose:
                print(f"Warning: {symbol} no data in alphavantage")
            return pd.DataFrame()
    
        avpp = pd.DataFrame.from_dict(rout['Time Series (Daily)'], 
                                      orient='index', dtype='float')[::-1]
        
        avpp.rename(columns={'1. open': 'open', '2. high': 'high',
                             '3. low': 'low', '4. close': 'close',
                             '5. adjusted close': 'adjusted', 
                             '6. volume': 'volume',
                             '7. dividend amount': 'divd', 
                             '8. split coefficient': 'split'}, inplace=True)
        avpp.index = pd.to_datetime(avpp.index)
        
        self._adjustSplit(avpp, volume=False, inplace=True)
    
        return avpp[self._col]
    
    
    def _alphavantage_yahoo(self, symbol, api_key=None):
        ''' 
        returns maximum period adjusted for splits only,
        prices are from alphavantage and dividends + splits are from yahoo
        works with a free key
        '''
        if api_key is None:
            api_key = os.getenv('ALPHAVANTAGE_API_KEY')
            if api_key is None:
                raise ValueError("Worning: no API key set as "
                                 "global environment variable")
                return pd.DataFrame()
            
        url = ('https://www.alphavantage.co/query?function=TIME_SERIES_DAILY'
               + '&symbol=' + symbol 
               + '&outputsize=full&datatype=csv'
               + '&apikey=' + api_key)
    
        aprice = pd.read_csv(url)
        if aprice.shape[1] < 5:
            print(aprice.loc[0][0])
            return pd.DataFrame()
    
        aprice = aprice[::-1]
        aprice.rename(columns={'timestamp': 'date'}, inplace=True)
        aprice.set_index('date', inplace=True)
        aprice.index = pd.to_datetime(aprice.index)
    
        ysymb = yf.Ticker(symbol)
        
        edprice = pd.DataFrame({'divd': ysymb.dividends})
        if edprice.empty:
            aprice['divd'] = 0.
        else:
            edprice.index = pd.to_datetime(edprice.index.date, utc=False)
            edprice.index.name = 'date'
            aprice = pd.merge(aprice, edprice, how='left', 
                              left_index=True, right_index=True)
            aprice['divd'].fillna(0., inplace=True)
            
        esprice = pd.DataFrame({'split': ysymb.splits})
        if esprice.empty:
            aprice['split'] = 1.
        else:
            esprice.index.name = 'date'
            esprice.index = pd.to_datetime(esprice.index.date, utc=False)
            aprice = pd.merge(aprice, esprice, how='left', 
                              left_index=True, right_index=True)
            aprice['split'].fillna(1., inplace=True)
        
        aprice['adjusted'] = aprice['close']
        self._adjustSplit(aprice, inplace=True)
        
        return aprice[self._col]
    
    
    def _marketstack(self, symbol, api_key=None):
        ''' 
        returns maximum period adjusted for splits only,
        free key works only for 1y of data
        '''
        if api_key is None:
            api_key = os.getenv('MARKETSTACK_API_KEY')
            if api_key is None:
                raise ValueError("Worning: no API key set as "
                                 "global environment variable")
                return pd.DataFrame()
                
        params = {'access_key': api_key, 'limit': 1000, 'offset': 0}
    
        req = requests.get('http://api.marketstack.com/v1/tickers/'
                            + str(symbol) + '/eod', params)
        if req.status_code != 200:
            if self._verbose:
                print(f"request status code = {req.status_code} no data"
                      f" collected for {symbol} from alphavantage")
            return pd.DataFrame()
        
        rout = req.json()
        if 'error' in rout.keys():
            if self._verbose:
                print(f"Worning: {rout['error']['message']}"
                      f" with error code: {rout['error']['code']}")
            return pd.DataFrame()
        
        pag = rout['pagination']
        nr = pag['total'] // pag['limit']
        res = [pd.DataFrame(rout['data']['eod'])]
        
        for k in range(nr):
            params['offset'] += params['limit']
            api_result = requests.get('http://api.marketstack.com/v1/tickers/'
                                      + str(symbol) + '/eod', params)
            api_response = api_result.json()
            res.append(pd.DataFrame(api_response['data']['eod']))
            
        xpp = pd.concat(res)[::-1]
        xpp['date'] = pd.to_datetime([x[0:10] for x in xpp['date']])
        xpp.set_index('date', inplace=True)
        xpp.rename(columns={'split_factor': 'split', 'dividend': 'divd',
                            'adj_close': 'adjusted'},
                   inplace=True)
        xpp = self._adjustSplit(xpp)
        
        return xpp[self._col]    
    
    
    def _adjustSplit(self, data, volume=True, unsplit=False, inplace=False):
        if inplace: 
            ldata = data
        else: 
            ldata = data.copy()
            
        spl = ldata.split.shift(-1, fill_value=1.)[::-1].cumprod()[::-1]
        if unsplit: spl = 1. / spl
        
        ldata.open /= spl
        ldata.close /= spl
        ldata.high /= spl
        ldata.low /= spl
        if volume: ldata.volume *= spl
        ldata.divd /= spl
            
        return ldata


    def _adjustDividend(self, data, field='close', divd='divd'):
        # A(i-1) = [ Z(i-1) - D(i)] * A(i) / Z(i)
        return (1 - data[divd].shift(-1, fill_value=0.) 
                / data[field])[::-1].cumprod()[::-1] * data[field]