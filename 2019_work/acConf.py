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
    'title'                   : 'pyALCELI',

    'dumpBunchIN'             : True,
    'dumpBunchOUT'            : True,
#todo: twissplot
    'twissPlot'               : False,  

    # display limits
    'ingnore_limits'          : False,
    'limx'                    : 201,
    'limxp'                   : 201,
    'limy'                    : 201,
    'limyp'                   : 201,
    'limz'                    : 161,
    'limzp'                   : 121,
    }
