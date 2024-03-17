import os
from azapy.MkT.MkTreader import MkTreader


def update_MkTdata(mktdir, source=None, api_key=None, param=None,
                   except_file=[], verbose=True):
    """
    Updates all mkt data saved in a directory.

    Parameters
    ----------
    mktdir : `str`
        Mkt data directory.
    source : `str`, optional
        Mkt data provider. 
        For more details see the `azapy.MkTreader.get` function doc. 
        The default is `None`.
    api_key : `str`, optional
        Mkt data provider API key.
        For more details see the `azapy.MkTreader.get` function doc. 
        The default is `None`
    param : `dict`, optional
        Additional parameters required by mkt data provider.
        For more details see the `azapy.MkTreader.get` function doc. 
        The default is `None`.
    except_file : `list`, optional
        List of symbols to be omitted from the update. The default is [].
    verbose : `Boolean`, optional
            - `True` will print a progress report,
            - `False` suppress any printing to the terminal.
        
        The default is `True`.

    Returns
    -------
    int : Error code
        - 200 : successful, everything updated
        - 201 : some (or all) were not completely updated 
        - 101 : the `mktdir` does not exists
        - 102 : unsupported mkt data sources
            
    Notes
    -----
    Files with unsupported extensions (see `azapy.MkTreader.get` function)
    are silently omitted from the update.
    """
    if not os.path.isdir(mktdir):
        if verbose:
            print(f"Warning: doesn\'t exist dir: {mktdir}")
        return 101
    if source is None:
        source = 'yahoo'
    elif not source in ['yahoo', 
                        'eodhistoricaldata', 'eodhistoricaldata_yahoo',
                        'alphavantage', 'alphavantage_yahoo',
                        'marketstack']:
        raise ValueError(f"Unsupported mkt data source {source}")

    if except_file is None:
        except_file = []
    
    sce = {}
    for file in os.listdir(mktdir):
        if verbose:
            print(file)
        fname, fext = os.path.splitext(file)
        if fext[1:] not in {'csv', 'json', 'feather'}:
            if verbose:
                print(f"file {file} omitted from update")
            continue
        if fname in except_file:
            if verbose:
                print(f"file {file} excepted from update")
            continue
        ld = {}
        ld['source'] = source
        ld['file_format'] = fext[1:]
        sce[fname] = ld

    mkt = MkTreader()
    mkt.get(source=sce, file_dir=mktdir, api_key=api_key, 
            param=param, verbose=verbose)
    
    elo = mkt.get_error_log()
    if verbose:
        sta = mkt.get_request_status()
        print("\nThe following symbols were updated:")
        for symb in sta.columns:
            print(f"symb {symb} : errors {sta.loc['error', symb]}")
        if len(elo) != 0:
            print(f"errors: {elo}")
          
    if len(elo) == 0: 
        return 200
    return 201