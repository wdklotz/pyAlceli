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

# import python modules customized for ALCELI
from acBunchGenerator import AcLinacBunchGenerator
from acLatticeFactory import AcLinacLatticeFactory
from acConf import CONF

# DEBUG
from orbit.utils.debugHelpers import caller_name, lineno, DEBUG_ON, DEBUG_OFF, DEXIT

DEBUG_MAIN = DEBUG_ON

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

def action_exit(paramsDict):
    bunch = paramsDict["bunch"]
    gamma = bunch.getSyncParticle().gamma()
    beta  = bunch.getSyncParticle().beta()
    node  = paramsDict["node"]
    DEBUG_ON(__file__,lineno(),'exit action at node: {}'.format(node.getName()))
    

def main():
    random.seed(100)

    # section list
    names = ["HE"]

    #---- the XML input file name with the linac structure
    xml_file_name = "../lattice.xml"

    #---- create the factory instance
    linac_factory = AcLinacLatticeFactory()
    linac_factory.setMaxDriftLength(0.01)
    
    #---- call FACTORY
    accLattice = linac_factory.getLinacAccLattice(names,xml_file_name)
    DEBUG_MAIN(__file__,lineno(),accLattice)
    print "Linac lattice is ready. L=",accLattice.getLength()

    #----set up RF Gap Model -------------
    #---- There are three available models at this moment
    #---- BaseRfGap  uses only E0TL*cos(phi)*J0(kr) with E0TL = const
    #---- MatrixRfGap uses a matrix approach like envelope codes
    #---- RfGapTTF uses Transit Time Factors (TTF) like PARMILA
    # cppGapModel = BaseRfGap()
    cppGapModel = MatrixRfGap()
    # cppGapModel = RfGapTTF
    rf_gaps = accLattice.getRF_Gaps()
    for rf_gap in rf_gaps:
        # DEBUG_MAIN(__file__,lineno(),rf_gap)
        rf_gap.setCppGapModel(cppGapModel)
    quads = accLattice.getQuads()
    for cnt,quad in enumerate(quads):
        quad.setUsageFringeFieldOUT(usage = False)
        quad.setUsageFringeFieldIN(usage = False)
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
    Tkin      = 70.e-3                          # in [GeV]
    mass      = CONF['proton_mass']*1.e-3       # in [GeV]
    frequency = 816.e+6                         # in [Hz]
    Clight    = CONF['lichtgeschwindigkeit']    # in [m/sec]
    gamma     = (mass + Tkin)/mass
    beta      = math.sqrt(gamma*gamma - 1.0)/gamma
    print "At injection: T= {}[GeV], gamma= {}, beta= {}".format(Tkin, gamma, beta)

    #------ emittances are normalized - transverse by gamma*beta and long. by gamma**3*beta
    (alphaX,betaX,emittX) = (-1.9620, 0.1831, 0.21)
    (alphaY,betaY,emittY) = ( 1.7681, 0.1620, 0.21)
    (alphaZ,betaZ,emittZ) = ( 0.0196, 0.5844, 0.24153)

    #---make emittances un-normalized XAL units [m*rad]
    emittX = 1.0e-6*emittX/(gamma*beta)
    emittY = 1.0e-6*emittY/(gamma*beta)
    emittZ = 1.0e-6*emittZ/(gamma**3*beta)
    print " ========= Twiss parameters at injection ==========="
    print " aplha beta emitt[mm*mrad] X= %6.4f %6.4f %6.4f "%(alphaX,betaX,emittX*1.0e+6)
    print " aplha beta emitt[mm*mrad] Y= %6.4f %6.4f %6.4f "%(alphaY,betaY,emittY*1.0e+6)
    print " aplha beta emitt[mm*mrad] Z= %6.4f %6.4f %6.4f "%(alphaZ,betaZ,emittZ*1.0e+6)

    #---- long. size in mm
    sizeZ = math.sqrt(emittZ*betaZ)*1.0e+3

    #---- transform to pyORBIT emittance[GeV*m]
    emittZ = emittZ*gamma**3*beta**2*mass
    betaZ = betaZ/(gamma**3*beta**2*mass)

    print " ========= PyORBIT parameters at injection ==========="
    print " aplha beta emitt[mm*mrad] X= %6.4f %6.4f %6.4f "%(alphaX,betaX,emittX*1.0e+6)
    print " aplha beta emitt[mm*mrad] Y= %6.4f %6.4f %6.4f "%(alphaY,betaY,emittY*1.0e+6)
    print " aplha beta emitt[mm*MeV] Z= %6.4f %6.4f %6.4f "%(alphaZ,betaZ,emittZ*1.0e+6)

    twissX = TwissContainer(alphaX,betaX,emittX)
    twissY = TwissContainer(alphaY,betaY,emittY)
    twissZ = TwissContainer(alphaZ,betaZ,emittZ)

    # BUNCH generation
    print "Start Bunch Generation",
    bunch_gen  = AcLinacBunchGenerator(twissX,twissY,twissZ)
    paramsDict = {'BunchGenerator': bunch_gen}
    #set the initial kinetic energy in GeV
    bunch_gen.setKinEnergy(Tkin)
    #set the beam peak current in mA
    bunch_gen.setBeamCurrent(1.0)
    bunch = bunch_gen.getBunch(nParticles = 100000, distributorClass = GaussDist3D)

    # DUMP bunch at lattice entrance
    if CONF['dumpBunchIN']:
        bunch.dumpBunch(CONF['bunchIn_filename'])
    print " - finished"

    # DESIGN tracking
    print "Design tracking started",
    accLattice.trackDesignBunch(bunch)
    print " - finished"

    # BUNCH tracking preparation
    bunch                = bunch_gen.getBunch(nParticles = 100000, distributorClass = GaussDist3D)
    paramsDict ['bunch'] = bunch
    actionContainer      = AccActionsContainer("Bunch Tracking")
    nodes                = accLattice.getNodes()
    # DEBUG_MAIN(__file__,lineno(),nodes)

    # BUNCH tracking
    print "Bunch tracking started",
    time_start = time.clock()
    # all but last node
    for node in nodes[:-2]:
        node.trackActions(actionContainer, paramsDict=paramsDict)
    # last node
    end_node = nodes[-1]
    actionContainer.addAction(action_exit, AccActionsContainer.EXIT)    
    end_node.trackActions(actionContainer, paramsDict=paramsDict)
    time_exec = time.clock() - time_start
    print " - finished in {:4.2f} [sec]".format(time_exec)
    # DEBUG_ON(__file__,lineno(),'last node: {}'.format(end_node.getName()))

    # DUMP bunch at lattice end
    if CONF['dumpBunchOUT']:
        # DEBUG_MAIN(__file__,lineno(),'bunch.getSize(): {}'.format(bunch.getSize()))
        bunch.dumpBunch(CONF['bunchOut_filename'])

if __name__ == '__main__':
    main()