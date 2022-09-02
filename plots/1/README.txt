
16 table entries with NaN and 1270/3000 with coordinates off-limits

--------------------------------------------------------------------------------------
Content of acConf.py:
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



-----------------------------------------------------------------------------

Content of orbit.yml:
# Copyright 2015 Wolf-Dieter Klotz <wdklotz@gmail.com>
# This file is part of the SIMULINAC code

# SIMULINAC is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# any later version.

# SIMULINAC is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with SIMULINAC.  If not, see <http://www.gnu.org/licenses/>.
#
# Input file for FODO linac simulator
# Input file follows YAML syntax (http:\\http://yaml.org/)
#
# Note:
#   Repeated nodes are initially denoted by an ampersand ( & )
#   and thereafter referenced with an asterisk ( * ).
#..........................................................................................
## RUN configuration

flags:
    # - accON:        False           # {True}  acceleration on/off flag
    # - egf:          True            # {False} emittance growth flag
    # - sigma:        False           # {True}  beam sizes by sigma-matrix
    # - KVprint:      True            # {False} print a dictionary of Key-Value pairs, no display
    # - periodic:     True            # {False} treat lattice as periodic cell sequence(else as transfer line)
    - express:      False           # {True}  use express version of thin quads
    # - verbose:      2               # {0}     print flag (0 = minimal print), try 0,1,2,3
#...........................................................................................
## SECTION definitions
sections:
- [&LE LE,  &HE HE]
#...........................................................................................
## INPUT parameter definitions

parameters:
    - aperture:             0.01        # [m] aperture = bore radius
    - Tkin:                 70.0        # [MeV] energy @ entrance (injection)
    - frequency:    &p01    816.e+6     # [Hz] frequency
    - ql0:          &p02    0.10        # [m] quad-length
    - ql:           &p03    0.05        # [m] 1/2 quad-length
    - emity_i:      &emy    2.e-6       # [m*rad] transverse emittance @ entrance x
    - emitx_i:              *emy        # [m*rad] transverse emittance @ entrance y
    # - betax_i:      &bex    7.084       # [m] twiss beta @ entrance x
    - betax_i:      &bex    2.800       # [m] twiss beta @ entrance x
    # - betay_i:              0.04        # [m] twiss beta @ entrance y (zero RF)
    - betay_i:              0.202       # [m] twiss beta @ entrance y (zero RF)
    - alfax_i:              0.0         # [1] twiss alpha x @ entrance
    - alfay_i:              0.0         # [1] twiss alpha y @ entrance
    - sigmaz_i:             1.e-3       # [m] longitudinal displacement z @ entrance
    # - sigmaz_i:             8.1795e-3   # separatrix? TTFG
    # - sigmaz_i:             8.1330e-3   # separatrix? RFG-RB
    - dp2p_i:               0.1         # [%] longitidinal dp/p spread    @ entrance
    - phi_sync:     &p06   -20.         # [deg] synchronous phase
    - windings:             30          # [1] quad-coil windings
    - GAP:          &p15    0.048       # [m] RF gap
#...........................................................................................
## ELEMENT definitions

elements:
# HE
    - D3:        &D3             # ID: &alias
        - type:     D            # D class
        # - type:     SIXD       # SIXD class
        - length:   0.08         # [m]
        - sec:      *HE          # is part of section
    - D5:    &D5                 # ID: &alias
        - type:     D            # D class
        # - type:     SIXD       # SIXD class
        - length:   0.022        # [m] half cavity
        - sec:      *HE          # is part of section
    - QFH:   &QFH                # ID: &alias
        - type:     QF           # QF class
        - length:   *p03         # [m]
        - B':       &Bgrad   40. # [T/m] quadrupole gradient
        # - B':       &Bgrad   30. # [T/m] quadrupole gradient
        # - B':       &Bgrad   38.8 # [T/m] quadrupole gradient  (zero RF)
        # - B':       &Bgrad  34.5 # [T/m] quadrupole gradient  (zero RF)
        - slices:   0            # if slices == 1: one thick element else: many thin elements
        - sec:      *HE          # is part of section
    - QDH:   &QDH                # ID: &alias
        - type:     QD           # QD class
        - length:   *p02         # [m] length
        - B':       *Bgrad       # [T/m] quadrupole gradient
        - slices:   0            # if slices == 1: one thick element else: many thin elements
        - sec:      *HE          # is part of section
    # - RFGH:    &RFGH
    #     - type:     D
    #     - length:   0.0
    #     - sec:      *HE
    - RFGH:  &RFGH               # ID: &alias
        - type:     RFG          # RFG class
        - Ez:       2.31         # [MV/m] (Ez)av electric field
        - PhiSync:  *p06         # [deg] synchronous phase
        - fRF:      *p01         # [Hz] frequency
        - gap:      *p15         # [m] length
        - SFdata:   SF_WDK2g44.TBL # superfish tbl-data file
        - Ezpeak:   4.00         # [MV/m] # corresdponds to (Ez)av = 2.3[MV/m] electric field
        - mapping:   t3d         # Trace 3D linear map model
        # - mapping:   simple      # Shishlo/Holmes linear map model
        - mapping:   base        # Shishlo/Holmes base map model
        # - mapping:   ttf         # Shishlo/Holmes three point TTF RF gap-model
        # - mapping:   dyn         # Tanke/Valero RF gap-model
        - sec:      *HE          # is part of section
    - RFCH:  &RFCH               # ID: &alias
        - type:     RFC          # RFC class
        - Ez:       2.31         # [MV/m] (Ez)av electric field
        - PhiSync:  *p06         # [deg] synchronous phase
        - fRF:      *p01         # [Hz] frequency
        - gap:      *p15         # [m] gap length
        - length:   0.044        # [m] cavity length
        - sec:      *HE          # is part of section
    - M:  &M                     # ID: &alias
        - type:    MRK           # MRK class
        - sec:     *HE
        - actions:
            - sigma_x
            - sigma_y
            - Tkin
    - F:  &Folie                 # ID: &alias
        - type:    MRK           # MRK class
        - sec:     *HE
#...........................................................................................
## SEGMENT definitions

segments:
# LE            # empty section
# HE            # high energy section
    - SEG1H:
        - *QFH
        - *D3
    - SEG2H:
        - *D3
        - *QDH
        - *D3
    - SEG3H:
        - *D3
        - *QFH
    - RFGH:
        - *D5
        - *RFGH
        - *D5
        #
        - *D5
        - *RFGH
        - *D5
        #
        - *D5
        - *RFGH
        - *D5
        #
        - *D5
        - *RFGH
        - *D5
        #
        - *D5
        - *RFGH
        - *D5
        #
        - *D5
        - *RFGH
        - *D5
        #
        - *D5
        - *RFGH
        - *D5
        #
        - *D5
        - *RFGH
        - *D5
        #
        - *D5
        - *RFGH
        - *D5
        #
        # - *D5     # 10th cavity makes it unstable!!
        # - *RFGH
        # - *D5
        #
    - RFCAVH:
        - *RFCH
        - *RFCH
        - *RFCH
        - *RFCH
        - *RFCH
        - *RFCH
        - *RFCH
        - *RFCH
        - *RFCH
    - MARK:
        - *Folie
#...........................................................................................
## LATTICE definition

lattice:
    - title:   1.0.0        # description -  DON't remove or reposition!
    # - [80,    SEG1H, RFCAVH, SEG2H, RFCAVH, SEG3H]     # cavities
    - [80,    SEG1H, RFGH, SEG2H, RFGH, SEG3H, MARK]   # gaps
    # - [80,    SEG1H, RFGH, SEG2H, RFGH, SEG3H]   # gaps, w/o MARKs

