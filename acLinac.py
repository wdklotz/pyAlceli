#! /usr/bin/env python

"""
This script will track the bunch through the ALCELI Linac.

At the beginning the lattice can be modified by replacing
the BaseRF_Gap nodes with AxisFieldRF_Gap nodes for
the selected sequences. These nodes will use the
RF fields at the axis of the RF gap to track the bunch.
The usual BaseRF_Gap nodes have a zero length.

The apertures are added to the lattice.
"""

import os
import sys
import math
import random
import time
import json

# from linac & bunch import the C++ RF gap classes
from linac import BaseRfGap, MatrixRfGap, RfGapTTF
from bunch import Bunch
from bunch import BunchTwissAnalysis
# from orbit import the python modules
from orbit.bunch_generators import TwissContainer
from orbit.bunch_generators import WaterBagDist3D, GaussDist3D, KVDist3D
from orbit.lattice import AccLattice, AccNode, AccActionsContainer
from orbit.py_linac.lattice_modifications import Add_quad_apertures_to_lattice
from orbit.py_linac.lattice_modifications import Add_rfgap_apertures_to_lattice
from orbit.py_linac.lattice_modifications import AddMEBTChopperPlatesAperturesToSNS_Lattice
from orbit.py_linac.lattice_modifications import AddScrapersAperturesToLattice
# Option: BaseRF_Gap to  AxisFieldRF_Gap replacement
from orbit.py_linac.lattice_modifications import Replace_BaseRF_Gap_to_AxisField_Nodes
from orbit.py_linac.lattice.LinacAccNodes import TiltElement,FringeField
import orbit
# import python modules customized for ALCELI
from acBunchGenerator import AcLinacBunchGenerator
from acLatticeFactory import AcLinacLatticeFactory
from acConf  import CONF
# import from SIMULINAC
from setutil import PARAMS,WConverter

# DEBUG
from acDebugHelpers import caller_name, lineno, DEBUG_ON, DEBUG_OFF, DEXIT
DEBUG_MAIN = DEBUG_OFF

# root dir of SIMULINAC
simulinacRoot = os.getenv('SIMULINAC_ROOT')

def tblprnt(headr,records):
    """
    Custom helper to print a nice table to memory (clever!?)
    IN:
        headr   = table header [...]
        records = table rows [[...]]
    OUT:
        s = the table as a string
    """
    rows = []; s=''
    rows.append(headr)
    for record in records:
        row = record
        rows.append(row)
    widths = [max(map(len, col)) for col in zip(*rows)]
    for row in rows:
            s+=" | ".join((val.ljust(width) for val,width in zip(row, widths)))+'\n'
    return s

def dumpBunch(bunch,fileName):
    placeholders = [
    '% placeholer --> PARTICLE_ATTRIBUTES_CONTROLLERS_NAMES',
    '% placeholer --> BUNCH_ATTRIBUTE_DOUBLE charge',
    '% placeholer --> BUNCH_ATTRIBUTE_DOUBLE classical_radius',
    '% placeholer --> BUNCH_ATTRIBUTE_DOUBLE macro_size',
    '% placeholer --> BUNCH_ATTRIBUTE_DOUBLE m0c2',
    '% placeholer --> SYNC_PART_COORDS x, y, z positions in [m]',
    '% placeholer --> SYNC_PART_MOMENTUM px, py, pz momentum component in GeV/c',
    '% placeholer --> SYNC_PART_X_AXIS x-axis ort coordinates',
    '% placeholer --> info only: energy of the synchronous particle [GeV]',
    '% placeholer --> info only: momentum of the synchronous particle [GeV/c]',
    '% placeholer --> info only: beta=v/c of the synchronous particle',
    '% placeholer --> info only: gamma=1/sqrt(1-(v/c)**2) of the synchronous particle',
    '% placeholer --> SYNC_PART_TIME time in [sec]',
    '% x[m] px[rad] y[m] py[rad] z[m]  (pz or dE [GeV]) '
    ]
    nParticles = bunch.getSize()
    with open(fileName,'w') as file:
        for placeholder in placeholders:
            file.write(placeholder+'\n')
        for i in range(nParticles):
            x  = bunch.x(i)
            xp = bunch.px(i)
            y  = bunch.y(i)
            yp = bunch.py(i)
            z  = bunch.z(i)
            zp = bunch.pz(i)
            file.write('{} {} {} {} {} {}'.format(x, xp, y, yp, z, zp)+ ' \n')
    print '-> bunch with {} partcles dumped to {}'.format(nParticles, fileName)

def action_exit(paramsDict):
    bunch = paramsDict["bunch"]
    gamma = bunch.getSyncParticle().gamma()
    beta  = bunch.getSyncParticle().beta()
    node  = paramsDict["node"]
    m0c2  = paramsDict['m0c2']
    Tkfin = m0c2*(gamma-1.)   # m0c2 in [MeV]
    if isinstance(node, FringeField):
        DEBUG_MAIN(__file__,lineno(),'exit action at node: {} --> usage: {}'.format(node.getName(),node.getUsage()))
    elif isinstance(node,TiltElement):
        DEBUG_MAIN(__file__,lineno(),'exit action at node: {} --> tilt angle: {}'.format(node.getName(),node.getTiltAngle()))
    else:
        DEBUG_MAIN(__file__,lineno(),'exit action at node: {} --> tkin[MeV] {}'.format(node.getName(),Tkfin))
    

#todo: use WConverter
#todo: strukturieren - zu viel sphargetti code!
#todo: use AxisField models
#todo: make twiss plots
#todo: read parameter from simu.py instead from xml-input
def main():
    random.seed(100)

    # section list
    names = ["S25to200"]

    #---- the XML input file name with the linac structure
    xml_file_name = simulinacRoot+"/lattice.xml"

    #---- create the FACTORY instance
    linac_factory = AcLinacLatticeFactory()
    linac_factory.setMaxDriftLength(0.01)
    
    #---- call FACTORY
    (accLattice,acc_da) = linac_factory.getLinacAccLattice(names,xml_file_name)
    DEBUG_MAIN(__file__,lineno(),accLattice)
    print "Linac lattice is ready. L=",accLattice.getLength()

    #----set up RF Gap Model -------------
    #---- There are three available models at this moment
    #---- MatrixRfGap uses a matrix approach like envelope codes
    #---- BaseRfGap  uses only E0TL*cos(phi)*J0(kr) with E0TL = const
    #---- RfGapTTF uses Transit Time Factors (TTF) like PARMILA
    cppGapModel = MatrixRfGap()
    cppGapModel = BaseRfGap()
    # cppGapModel = RfGapTTF
    rf_gaps = accLattice.getRF_Gaps()
    for rf_gap in rf_gaps:
        # DEBUG_MAIN(__file__,lineno(),rf_gap)
        rf_gap.setCppGapModel(cppGapModel)
    quads = accLattice.getQuads()
    for cnt,quad in enumerate(quads):
        quad.setUsageFringeFieldOUT(usage = False)
        quad.setUsageFringeFieldIN(usage  = False)
        if cnt == 1:
            # DEBUG_MAIN(__file__,lineno(),quad.__dict__)
            # DEBUG_MAIN(__file__,lineno(),quad.getNodeFringeFieldIN().__dict__)
            # DEBUG_MAIN(__file__,lineno(),quad.getNodeFringeFieldOUT().__dict__)
            # DEBUG_MAIN(__file__,lineno(),quad.getNodeTiltIN().__dict__)
            # DEBUG_MAIN(__file__,lineno(),quad.getNodeTiltOUT().__dict__)
            pass

    # get PARAMS from xml-lattice
    [params_da] = acc_da.childAdaptors(name='PARAMS')
    DEBUG_MAIN(__file__,lineno(),params_da.getAttributes())

    # twiss parameters at the entrance
    tkin      = params_da.doubleValue('injection_energy')        # in [MeV]
    m0c2      = params_da.doubleValue('proton_mass')             # in [MeV]
    frequency = params_da.doubleValue('frequenz')                # in [Hz]
    clight    = params_da.doubleValue('clight')                  # in [m/sec]
    lamb      = clight/frequency                                 # in [m]
    gamma     = (m0c2 + tkin)/m0c2
    pi        = math.pi
    beta      = math.sqrt(gamma**2 - 1.0)/gamma
    betax_i   = params_da.doubleValue('betax_i')    # [m]
    betay_i   = params_da.doubleValue('betay_i')    # [m]
    betaz_i   = params_da.doubleValue('betaz_i')    # [m/rad]
    alfax_i   = params_da.doubleValue('alfax_i')    # []
    alfay_i   = params_da.doubleValue('alfay_i')    # []
    alfaz_i   = params_da.doubleValue('alfaz_i')    # []
    emitx_i   = params_da.doubleValue('emitx_i')    # [m*rad]
    emity_i   = params_da.doubleValue('emity_i')    # [m*rad]
    emitz_i   = params_da.doubleValue('emitz_i')    # [m*rad]
    emitw_i   = params_da.doubleValue('emitw_i')    # [rad]
    
    print "At injection: T= {}[GeV], gamma= {}, beta= {}".format(tkin, gamma, beta)

    print " ========= Twiss parameters at injection ==========="
    print " aplha beta emitt[mm*mrad] X= %6.3g %6.3g %6.3g "%(alfax_i,betax_i,emitx_i*1.0e+6)
    print " aplha beta emitt[mm*mrad] Y= %6.3g %6.3g %6.3g "%(alfay_i,betay_i,emity_i*1.0e+6)
    print " aplha beta emitt[mm*rad]  Z= %6.3g %6.3g %6.3g "%(alfaz_i,betaz_i,emitz_i*1.0e+3)

    #-----TWISS Parameters at the entrance of MEBT ---------------
    #-----transverse emittances are unnormalized and in [pi*mm*mrad]
    #-----longitudinal emittance is in [pi*m*GeV]

    #---- transform to pyORBIT (apparently {z-DW} phase space)
    emitzW  = m0c2*gamma*beta**2*emitz_i*1.e-3        # [m*GeV]
    betazW  = 1./(m0c2*gamma*beta**2)*betaz_i*1.e+3   # [m/GeV]
    
    print " ========= PyORBIT parameters at injection ==========="
    print " aplha beta[mm/mrad] emitt[mm*mrad] X= %6.3g %6.3g %6.3g "%(alfax_i,betax_i,emitx_i*1.0e+6)
    print " aplha beta[mm/mrad] emitt[mm*mrad] Y= %6.3g %6.3g %6.3g "%(alfay_i,betay_i,emity_i*1.0e+6)
    print " aplha beta[m/Gev]   emitt[m*GeV]   Z= %6.3g %6.3g %6.3g "%(alfaz_i,betazW,emitzW)

    #-----longitudinal emittance is in [pi*m*GeV]
    twissX = TwissContainer(alfax_i,betax_i,emitx_i)
    twissY = TwissContainer(alfay_i,betay_i,emity_i)
    twissZ = TwissContainer(alfaz_i,betazW,emitzW)

    # BUNCH generation
    print "-> Start Bunch Generation"
    bunch_gen  = AcLinacBunchGenerator(twissX,twissY,twissZ,frequency=frequency)
    paramsDict = {'BunchGenerator': bunch_gen}
    #----------------------------------------
    # set the initial kinetic energy in [GeV]
    bunch_gen.setKinEnergy(tkin*1.e-3)
    #----------------------------------------
    #set the beam peak current in mA
    # bunch_gen.setBeamCurrent(PARAMS['elementarladung']*PARAMS['frequenz']*1.e3)   # 1 e-charge per bunch
    bunch_gen.setBeamCurrent(10.)
    bunch = bunch_gen.getBunch(nParticles = 5000, distributorClass = GaussDist3D)
    # print '\npossible particle attributes names:\n'+''.join(['\t"{}"\n'.format(i) for i in bunch.getPossiblePartAttrNames()])

    # DUMP bunch at lattice entrance
    if CONF['dumpBunchIN']:
        bunch.dumpBunch(CONF['bunchIn_filename'])
    print "-> Bunch Generation finished"

    # DESIGN tracking
    print "-> Design tracking started"
    accLattice.trackDesignBunch(bunch)
    print "-> Design tracking finished "

    # BUNCH tracking preparation
    accLattice.setLinacTracker(switch=False)    # use TeapotBase (TPB) tracking
    paramsDict = {"old_pos":-1.,"count":0,"pos_step":0.1,'m0c2':m0c2}
    last_node_index = len(accLattice.getNodes())-1
    nodes           = accLattice.getNodes()[:last_node_index-1]
    last_node       = accLattice.getNodes()[last_node_index]
    # DEBUG_MAIN(__file__,lineno(),nodes)
    DEBUG_MAIN(__file__,lineno(),'last node: {}'.format(last_node.getName()))

    # BUNCH tracking
    print "-> Bunch tracking started "
    time_start = time.clock()
    # all but last node
    for node in nodes:
        node.trackBunch(bunch, paramsDict=paramsDict)
    # last node action
    actionsContainer = AccActionsContainer("Bunch Tracking")
    actionsContainer.addAction(action_exit, AccActionsContainer.EXIT)    
    last_node.trackBunch(bunch, paramsDict=paramsDict, actionContainer=actionsContainer)
    time_exec = time.clock() - time_start
    print "-> Bunch tracking finished in {:4.2f} [sec], T-final[MeV] {}".format(time_exec,bunch.getSyncParticle().kinEnergy()*1.e3)

    # DUMP bunch at lattice end
    if CONF['dumpBunchOUT']:
        # DEBUG_MAIN(__file__,lineno(),'bunch.getSize(): {}'.format(bunch.getSize()))
        bunch.dumpBunch(CONF['bunchOut_filename'])
        # dumpBunch(bunch,CONF['bunchOut_filename'])

if __name__ == '__main__':
    main()