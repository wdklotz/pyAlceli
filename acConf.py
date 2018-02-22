## CONF
CONF = {      ## CONFIG constants, names and units....
    'lichtgeschwindigkeit'    : 299792458.,          # [m/s] const
    'elementarladung'         : 1.602176565e-19,     # [coulomb] const
    'proton_mass'             : 938.272,             # [MeV/c**2] const
    'electron_mass'           : 0.5109989,           # [MeV/c**2] const
    'frequenz'                : 816.6e6,             # [Hz] frequency
    'injection_energy'        : 70.,                 # [MeV] eKin-in
    'emitx_i'                 : 2.0,                 # [mm*mrad] emittance-x @ entrance
    'emity_i'                 : 2.0,                 # [mm*mrad] emittance-y @ entrance
    'emitz_i'                 : 0.2,                 # [mm*mrad] emittance-z @ entrance
    # 'betax_i'                 : 0.7084,              # [m] twiss betax @ entrance
    # 'betay_i'                 : 0.7084,              # [m] twiss betax @ entrance
    'betax_i'                 : 2.8000,              # [m] twiss betax @ entrance
    'betay_i'                 : 0.2000,              # [m] twiss betax @ entrance
    'betaz_i'                 : 0.01,                # [m] twiss betaz @ entrance
    'alfax_i'                 : 0.0,                 # twiss alphax @ entrance
    'alfay_i'                 : 0.0,                 # twiss alphaxy @ entrance
    'alfaz_i'                 : 0.0,                 # twiss alphaxz @ entrance

     # 'quad_gradient'           : 38.8,                # [T/m] DB/dr
    # 'E0TL'                  : -0.0500e-3,          # [GeV] E0TL effective gap voltage
    # 'E0L'                   : -0.0500e-3,          # [GeV] E0L average gap voltage

    'twiss_filename'          : 'twiss.dat',
    'result_filename'         : 'result.dat',
    'bunchOut_filename'       : 'bunchf.dat'   ,
    'bunchIn_filename'        : 'bunchi.dat',
    'title'                   : 'ALCELI',

    'dumpBunchIN'             : True,
    'dumpBunchOUT'            : True,
    'twissPlot'               : False,
    
    # display limits
    'limx'                    : 10,
    'limxp'                   : 10,
    'limy'                    : 10,
    'limyp'                   : 10,
    'limz'                    : 5000,
    'limzp'                   : 5000,
    }

#     'dWf': False,                # acceleration off/on flag
#     'E0z_feld': 1.14,            # [MV/m] E(r=0,z)-field strength
#     'spalt_laenge': 0.044,       # [m] rf gap length
#     'cav_laenge': 0.044,         # [m] cavity length
#     'quad_laenge': 0.1,          # [m] full quadrupole length
#     'drfs_laenge': 0.08,         # [m] drift before and after rf section
#     'soll_phase': -90.0,         # [deg] sync phase
