## CONF
"""
    specific configuration for pyAlceli.py
    NOTE: pyAlceli.py reads physical parameters from SIMULINAC's setutil.py now
"""
CONF = {
    'twiss_filename'          : 'twiss.dat',
    'result_filename'         : 'result.dat',
    'bunchOut_filename'       : 'bunchf.dat'   ,
    'bunchIn_filename'        : 'bunchi.dat',
    'title'                   : 'ALCELI',

    'dumpBunchIN'             : True,
    'dumpBunchOUT'            : True,
    'twissPlot'               : False,
    
    # display limits
    'limx'                    : 101,
    'limxp'                   : 101,
    'limy'                    : 101,
    'limyp'                   : 101,
    'limz'                    : 101,
    'limzp'                   : 1010,
    }
