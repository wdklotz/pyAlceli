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
from setutil import PARAMS
# DEBUG
from debugHelpers import caller_name, lineno, DEBUG_ON, DEBUG_OFF, DEXIT
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
    '% placeholer --> BUNCH_ATTRIBUTE_DOUBLE mass',
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
    Tkfin = m0c2*(gamma-1.)*1.e3
    if isinstance(node, FringeField):
        DEBUG_MAIN(__file__,lineno(),'exit action at node: {} --> usage: {}'.format(node.getName(),node.getUsage()))
    elif isinstance(node,TiltElement):
        DEBUG_MAIN(__file__,lineno(),'exit action at node: {} --> tilt angle: {}'.format(node.getName(),node.getTiltAngle()))
    else:
        DEBUG_MAIN(__file__,lineno(),'exit action at node: {} --> Tkin[MeV] {}'.format(node.getName(),Tkfin))
    

def main():
    random.seed(100)

    # section list
    names = ["HE"]

    #---- the XML input file name with the linac structure
    xml_file_name = simulinacRoot+"/lattice.xml"

    #---- create the FOCTORY instance
    linac_factory = AcLinacLatticeFactory()
    linac_factory.setMaxDriftLength(0.01)
    
    #---- call FACTORY
    accLattice = linac_factory.getLinacAccLattice(names,xml_file_name)
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
        
    # twiss parameters at the entrance
    # transverse emittances are unnormalized and in pi*mm*mrad
    # longitudinal emittance is in pi*eV*sec
    Tkin      = PARAMS['injection_energy']*1.e-3  # in [GeV]
    mass      = PARAMS['proton_mass']*1.e-3       # in [GeV]
    frequency = PARAMS['frequenz']                # in [Hz]
    clight    = PARAMS['lichtgeschwindigkeit']    # in [m/sec]
    gamma     = (mass + Tkin)/mass
    beta      = math.sqrt(gamma**2 - 1.0)/gamma
    betax_i   =PARAMS['betax_i']    # [m]
    betay_i   =PARAMS['betay_i']    # [m]
    betaz_i   =PARAMS['betaz_i']    # [m]
    alfax_i   =PARAMS['alfax_i']    # []
    alfay_i   =PARAMS['alfay_i']    # []
    alfaz_i   =PARAMS['alfaz_i']    # []
    emitx_i   =PARAMS['emitx_i']*(gamma*beta)     # [m*rad]
    emity_i   =PARAMS['emity_i']*(gamma*beta)     # [m*rad]
    emitz_i   =PARAMS['emitz_i']*(gamma**3*beta)  # [m*rad]
    print "At injection: T= {}[GeV], gamma= {}, beta= {}".format(Tkin, gamma, beta)

    #------ emittances normalized - transverse by gamma*beta and long. by gamma**3*beta
    (alphaX,betaX,emittX) = (alfax_i, betax_i, emitx_i)
    (alphaY,betaY,emittY) = (alfay_i, betay_i, emity_i)
    (alphaZ,betaZ,emittZ) = (alfaz_i, betaz_i, emitz_i)

    alphaZ = - alphaZ
    
    #---make emittances un-normalized XAL units [m*rad]
    emittX = emittX/(gamma*beta)
    emittY = emittY/(gamma*beta)
    emittZ = emittZ/(gamma**3*beta)
    print " ========= Twiss parameters at injection ==========="
    print " aplha beta emitt[mm*mrad] X= %6.4f %6.4f %6.4f "%(alphaX,betaX,emittX*1.0e+6)
    print " aplha beta emitt[mm*mrad] Y= %6.4f %6.4f %6.4f "%(alphaY,betaY,emittY*1.0e+6)
    print " aplha beta emitt[mm*mrad] Z= %6.4f %6.4f %6.4f "%(alphaZ,betaZ,emittZ*1.0e+6)

    #---- long. size in [mm]
    sizeZ = math.sqrt(emittZ*betaZ)*1.0e+3
    #---- transform to pyORBIT emittance[GeV*m]
    emittZ = emittZ*gamma**3*beta**2*mass
    betaZ  = betaZ/(gamma**3*beta**2*mass)

    print " ========= PyORBIT parameters at injection ==========="
    print " aplha beta emitt[mm*mrad] X= %6.4f %6.4f %6.4f "%(alphaX,betaX,emittX*1.0e+6)
    print " aplha beta emitt[mm*mrad] Y= %6.4f %6.4f %6.4f "%(alphaY,betaY,emittY*1.0e+6)
    print " aplha beta emitt[mm*MeV]  Z= %6.4f %6.4f %6.4f "%(alphaZ,betaZ,emittZ*1.0e+6)

    twissX = TwissContainer(alphaX,betaX,emittX)
    twissY = TwissContainer(alphaY,betaY,emittY)
    twissZ = TwissContainer(alphaZ,betaZ,emittZ)

    # BUNCH generation
    print "-> Start Bunch Generation"
    bunch_gen  = AcLinacBunchGenerator(twissX,twissY,twissZ,frequency=PARAMS['frequenz'])
    paramsDict = {'BunchGenerator': bunch_gen}
    #set the initial kinetic energy in GeV
    bunch_gen.setKinEnergy(Tkin)
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
    paramsDict = {"old_pos":-1.,"count":0,"pos_step":0.1,'m0c2':mass}
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
    # last node
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