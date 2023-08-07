from azapy.Generators.Port_Generator import Port_Generator


class _Port_Generator(Port_Generator):
    
    def __init__(self, mktdata, symb=None, sdate=None, edate=None, 
                 col_price='close', col_divd='divd', col_ref='adjusted',
                 col_calib='adjusted',
                 pname='Port', pcolname=None, capital=100000, 
                 schedule=None,
                 freq='Q', noffset=-3, fixoffset=-1, histoffset=3.25, 
                 calendar=None, multithreading=True):
        """
        Constructor
    
        Parameters
        ----------
        mktdata : `pandas.DataFrame`
            MkT data in the format "symbol", "date", "open", "high", "low",
            "close", "volume", "adjusted", "divd", "split" (e.g., as returned
            by `azapy.readMkT` function).
        symb : `list`, optional
            List of symbols for the basket components. All symbols MkT data
            should be included in mktdata. If set to `None` the `symb` will be
            set to include all the symbols from `mktdata`. The default
            is `None`.
        sdate : date like, optional
            Start date for historical data. If set to `None` the `sdate` will
            be set to the earliest date in mktdata. The default is `None`.
        edate : date like, optional
            End date for historical dates and so the simulation. Must be
            greater than  `sdate`. If it is `None` then `edate` will be set
            to the latest date in mktdata. The default is `None`.
        col_price : `str`, optional
            Column name in the mktdata DataFrame that will be considered
            for portfolio aggregation. The default is `'close'`.
        col_divd :  `str`, optional
            Column name in the mktdata DataFrame that holds the dividend
            information. The default is `'dvid'`.
        col_ref : `str`, optional
            Column name in the mktdata DataFrame that will be used as a price
            reference for portfolio components. The default is `'adjusted'`.
        col_calib : `str`, optional
            Column name used for historical weights calibrations. 
            The default is `'adjusted'`.
        pname : `str`, optional
            The name of the portfolio. The default is `'Port'`.
        pcolname : `str`, optional
            Name of the portfolio price column. If it set to `None` then
            `pcolname=pname`. The default is `None`.
        capital : `float`, optional
            Initial portfolio Capital in dollars. The default is `100000`.
        schedule : `pandas.DataFrame`, optional
            Rebalancing schedule, with columns for `'Droll'` rolling date and
            `'Dfix'` fixing date. If it is `None` than the schedule will be set
            using the `freq`, `noffset`, `fixoffset` and `calendar`
            information. The default is `None`.
        freq : `str`, optional
            Rebalancing frequency. It can be `'Q'` for quarterly or `'M'` for
            monthly rebalancing, respectively. It is relevant only if the 
            schedule is `None`. The default is `'Q'`.
        noffset : `int`, optional
            Number of business days offset for rebalancing date `'Droll'`
            relative to the end of the period (quart or month). A positive
            value add business days beyond the calendar end of the period while
            a negative value subtracts business days. It is relevant only if
            the schedule is `None`. The default is `-3`.
        fixoffset : `int`, optional
            Number of business day offset of fixing date `'Dfix'` relative to
            the rebalancing date `'Droll'`. It can be 0 or negative. It is
            relevant only if the schedule is `None`. The default is `-1`.
        calendar : `numpy.busdaycalendar`, optional
            Business calendar. If it is `None` then it will be set to NYSE
            business calendar. The default
            value is `None`.
        multitreading : Boolean, optional
            If it is `True`, then the  rebalancing weights will 
            be computed concurrent. The default is `True`.
    
        Returns
        -------
        The object.
        """
        
        super().__init__(mktdata, symb=symb, sdate=sdate, edate=edate, 
                     col_price=col_price, col_divd=col_divd, col_ref=col_ref,
                     pname=pname, pcolname=pcolname, capital=capital, 
                     schedule=schedule,
                     freq=freq, noffset=noffset, fixoffset=fixoffset, 
                     histoffset=histoffset, calendar=calendar,
                     multithreading=multithreading)
        self.col_calib = col_calib
        
        
    def get_weights(self, fancy=False):
        """
        Returns the portfolio weights at each rebalancing period.
        
        Parameters
        ----------
        fancy : Boolean, optional
            - `False`: reports the weights in algebraic format.
            - `True`: reports the weights in percentage rounded to 2 decimals.  
            
            The default is `False`.

        Returns
        -------
        `pandas.DataFrame` : portfolio weights per symbol.
        """
        return super().get_weights(fancy=fancy).drop('_CASH_', axis=1)
