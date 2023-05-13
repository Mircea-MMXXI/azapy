import azapy as az

sdate = "2012-01-01"
edate = 'today'
symb = ['GLD', 'TLT', 'XLV', 'VGT', 'PSJ']

mktdir = "../../MkTdata"

# simple calls (most often used) 
# returns a pd.DataFrame
mktdata = az.readMkT(symb, sdate=sdate, edate=edate, file_dir=mktdir)

# returns a dict of pd.DataFrame
mktdata_dict = az.readMkT(symb, sdate=sdate, edate=edate, file_dir=mktdir,
                          output_format='dict')

# complex call (extremely customized)
# Note: 'alhavantage' and 'eofhistoricaldata' may requier valid API keys
source = {'GLD': {'force': True,
                  'save': True,
                  'file_dir': '../../MkTdata_yahoo',
                  'foramt_foramt' : 'feather'
                 },
          'TLT': {'source' : 'alhavantage',
                  'force': False,
                  'save': True,
                  'file_dir': '../../MkTdata_av',
                  'foramt_foramt': 'json',
                  'param': {'max_req_per_min': 75}
                 },
          'XLV': {'source': 'eofhistoricaldata'}
         }

symb = ['VGT', 'PSJ']

mktdata = az.readMkT(symb, sdate=sdate, edate=edate, source=source,
                     file_dir=mktdir, file_format='csv')