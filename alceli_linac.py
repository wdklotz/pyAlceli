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
# import python customized for ALCELI
from alceli_linac_bunch_generator import ALCELI_Linac_BunchGenerator
from alceli_linac_lattice_factory import ALCELI_LinacLatticeFactory
from alceli_conf import CONF

# DEBUG
import pprint
import inspect
PP = pprint.PrettyPrinter(indent=4).pprint
def lineno():
    return inspect.currentframe().f_back.f_lineno
def DEBUG(line,arg):
    file = os.path.basename(__file__)
    if isinstance(arg,str):
        print'DEBUG[{}-{}]: '.format(line,file)+arg
    elif isinstance(arg,(tuple,list,dict)):
        print'DEBUG[{}-{}]: '.format(line,file),
        for i in arg:
            PP(i)
    else:
        print'DEBUG[{}-{}]: '.format(line,file)+repr(arg)
def DEBUG_ON(*args):
    DEBUG(*args)
def DEBUG_OFF(*args):
    pass

DEBUG_MAIN = DEBUG_OFF

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
    node      = paramsDict["node"]
    bunch     = paramsDict["bunch"]
    lattice   = paramsDict["lattice"]
    bunch_gen = paramsDict["BunchGenerator"]
    pos_start = paramsDict['pos_start']
    end_node = lattice.getNodes()[-1]
    if CONF['dumpBunchOUT'] and end_node == node:
        DEBUG_MAIN('bunch.getSize(): ',bunch.getSize())
        bunch.dumpBunch(CONF['bunchOut_filename'])
    pos = paramsDict["path_length"]
    DEBUG_MAIN(lineno(),'path_length {}'.format(pos))
    if(paramsDict["old_pos"] == pos): return
    if(paramsDict["old_pos"] + paramsDict["pos_step"] > pos): return
    paramsDict["old_pos"] = pos
    paramsDict["count"]  += 1
    gamma = bunch.getSyncParticle().gamma()
    beta  = bunch.getSyncParticle().beta()
    twiss_analysis = BunchTwissAnalysis()
    twiss_analysis.analyzeBunch(bunch)
    x_rms = math.sqrt(twiss_analysis.getTwiss(0)[1]*twiss_analysis.getTwiss(0)[3])*1000.
    y_rms = math.sqrt(twiss_analysis.getTwiss(1)[1]*twiss_analysis.getTwiss(1)[3])*1000.
    z_rms = math.sqrt(twiss_analysis.getTwiss(2)[1]*twiss_analysis.getTwiss(2)[3])*1000.
    z_to_phase_coeff = bunch_gen.getZtoPhaseCoeff(bunch)
    z_rms_deg = z_to_phase_coeff*z_rms/1000.0
    nParts = bunch.getSizeGlobal()
    (alphaX,betaX,emittX) = (twiss_analysis.getTwiss(0)[0],twiss_analysis.getTwiss(0)[1],twiss_analysis.getTwiss(0)[3]*1.0e+6)
    (alphaY,betaY,emittY) = (twiss_analysis.getTwiss(1)[0],twiss_analysis.getTwiss(1)[1],twiss_analysis.getTwiss(1)[3]*1.0e+6)
    (alphaZ,betaZ,emittZ) = (twiss_analysis.getTwiss(2)[0],twiss_analysis.getTwiss(2)[1],twiss_analysis.getTwiss(2)[3]*1.0e+6)
    norm_emittX           = emittX*gamma*beta
    norm_emittY           = emittY*gamma*beta
    #---- phi_de_emittZ will be in [pi*deg*MeV]
    phi_de_emittZ = z_to_phase_coeff*emittZ
    eKin          = bunch.getSyncParticle().kinEnergy()*1.0e+3

    # dump into result file
    s  = " %35s %4.5f "%(node.getName(),pos+pos_start)
    s += "%6.4f %6.4f %6.4f %6.4f "%(alphaX,betaX,emittX,norm_emittX)
    s += "%6.4f %6.4f %6.4f %6.4f "%(alphaY,betaY,emittY,norm_emittY)
    s += "%6.4f %6.4f %6.4f %6.4f "%(alphaZ,betaZ,emittZ,phi_de_emittZ)
    s += "%5.3f %5.3f %5.3f "%(x_rms,y_rms,z_rms_deg)
    s += "%10.6f %8d"%(eKin,nParts)
    CONF['result_file'].write(s +"\n")
    CONF['result_file'].flush()

    # result to console
    console_buffer_row = []
    console_buffer_row.append(  '%5d'%paramsDict["count"])
    console_buffer_row.append(  '%25s'%node.getName())
    console_buffer_row.append( '%4.5f'%(pos+pos_start))
    console_buffer_row.append( '%5.3f'%x_rms)
    console_buffer_row.append( '%5.3f'%y_rms)
    console_buffer_row.append( '%5.3f'%z_rms_deg)
    console_buffer_row.append('%10.6f'%eKin)
    console_buffer_row.append(   '%8d'%nParts)
    CONF['console_buffer_rows'].append(console_buffer_row)

    # result to memory memory for plot
    viseo = 0.
    plot_result = dict(
        name          = node.getName(),       #0
        position      = pos+pos_start,
        alfax         = alphaX,
        betax         = betaX,                #3
        emittx        = emittX,
        emittxn       = norm_emittX,
        alfay         = alphaY,
        betay         = betaY,                #7
        emitty        = emittY,
        emittyn       = norm_emittY,
        alfaz         = alphaZ,
        betaz         = betaZ,                #11
        emittz        = emittZ,
        phi_de_emittz = phi_de_emittZ,
        xrms          = x_rms,                #14
        yrms          = y_rms,
        zrms_deg      = z_rms_deg,
        ekin          = eKin,                 #17
        nparts        = nParts,
        viseo         = viseo)
    CONF['plot_data'].append(plot_result)

# def action_exit(paramsDict):
#         action_entrance(paramsDict)

def main():
    random.seed(100)

    # section list
    names = ["HE"]

    #---- the XML input file name with the linac structure
    xml_file_name = "../lattice.xml"

    #---- create the factory instance
    alceli_linac_factory = ALCELI_LinacLatticeFactory()
    alceli_linac_factory.setMaxDriftLength(0.01)
    accLattice = alceli_linac_factory.getLinacAccLattice(names,xml_file_name)
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
        DEBUG_MAIN(lineno(),rf_gap)
        rf_gap.setCppGapModel(cppGapModel)

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

    # prepare results file
    result_file = open(CONF['result_filename'],"w")
    CONF['result_file'] = result_file
    s = "%35s %4s "%('Node','position')
    s += "%6s %6s %6s %6s "%('alphaX','betaX','emittX','normEmittX')
    s += "%6s %6s %6s %6s "%('alphaY','betaY','emittY','normEmittY')
    s += "%6s %6s %6s %6s "%('alphaZ','betaZ','emittZ','emittZphiMeV')
    s += "%5s %5s %5s "%('sizeX','sizeY','sizeZ_deg')
    s += "%10s %8s "%('eKin','Nparts')
    result_file.write(s+"\n")

    # prepare for plot data
    CONF['plot_data'] = []

    # BUNCH generation
    print "Start Bunch Generation",
    bunch_gen  = ALCELI_Linac_BunchGenerator(twissX,twissY,twissZ)
    paramsDict = {'BunchGenerator':bunch_gen}
    #set the initial kinetic energy in GeV
    bunch_gen.setKinEnergy(Tkin)
    #set the beam peak current in mA
    bunch_gen.setBeamCurrent(1.0)
    bunch_in = bunch_gen.getBunch(nParticles = 100000, distributorClass = GaussDist3D)

    if CONF['dumpBunchIN']:
        bunch_in.dumpBunch(CONF['bunchIn_filename'])
    print " - finished"

    # DESIGN tracking
    print "Design tracking started",
    #set up design
    accLattice.trackDesignBunch(bunch_in)
    print " - finished"

    # BUNCH tracking preparation
    actionContainer         = AccActionsContainer("Test Design Bunch Tracking")
    paramsDict["old_pos"]   = -1.
    paramsDict["count"]     = 0
    paramsDict["pos_step"]  = 0.1
    paramsDict['pos_start'] = 0.

    # console buffer preparation
    console_buffer_header = ['N','node','pos','sizeX','sizeY','sizeZ[deg]','eKin','Nparts']
    # console_buffer_header = ['N','node','pos','sizeX']
    console_buffer_rows   = []
    CONF['console_buffer_rows'] = console_buffer_rows

    # registration of actions
    # actionContainer.addAction(action_entrance, AccActionsContainer.ENTRANCE)
    actionContainer.addAction(action_exit, AccActionsContainer.EXIT)

    # bunch tracking
    print "Bunch tracking started",
    time_start = time.clock()
    accLattice.trackBunch(bunch_in, paramsDict = paramsDict, actionContainer = actionContainer)
    time_exec = time.clock() - time_start
    print " - finished in {:4.2f} [sec]".format(time_exec)
    result_file.close()

    # dump results for console and plots
    print "RESULTS - RESULTS - RESULTS - RESULTS - RESULTS - RESULTS - RESULTS - RESULTS - RESULTS - RESULTS "
    print tblprnt(console_buffer_header,console_buffer_rows)
    if CONF['twissPlot']:
        #save structures data with json: here track results for plotting
        with open(CONF["twiss_filename"],"w") as file:
            json.dump(CONF['plot_data'], file)

if __name__ == '__main__':
    main()


    #------------------------------------------------------------------
    #---- BaseRF_Gap to  AxisFieldRF_Gap direct replacement
    #---- in the MEBT, CCL, SCLMed,SCLHigh  it could be done directly
    #---- because rf fields cover drifts only.
    #---- The DTL needs a special treatment.
    #------------------------------------------------------------------

    #---- axis fields files location
    # dir_location = ""
    # z_step = 4.e-4
    # Replace_BaseRF_Gap_to_AxisField_Nodes(accLattice,z_step,dir_location,names)

    # print "Linac lattice has been modified. New L[m] = ",accLattice.getLength()

    #-----------------------------------------------------
    # Set up Space Charge Acc Nodes
    #-----------------------------------------------------
    # from orbit.space_charge.sc3d import setSC3DAccNodes, setUniformEllipsesSCAccNodes
    # from spacecharge import SpaceChargeCalcUnifEllipse, SpaceChargeCalc3D
    # sc_path_length_min = 0.02
    #
    # print "Set up Space Charge nodes. "

    # set of uniformly charged ellipses Space Charge
    # nEllipses = 1
    # calcUnifEllips = SpaceChargeCalcUnifEllipse(nEllipses)
    # space_charge_nodes = setUniformEllipsesSCAccNodes(accLattice,sc_path_length_min,calcUnifEllips)


    # set FFT 3D Space Charge
    # sizeX = 64
    # sizeY = 64
    # sizeZ = 64
    # calc3d = SpaceChargeCalc3D(sizeX,sizeY,sizeZ)
    # space_charge_nodes =  setSC3DAccNodes(accLattice,sc_path_length_min,calc3d)

    # max_sc_length = 0.
    # min_sc_length = accLattice.getLength()
    # for sc_node in space_charge_nodes:
        # scL = sc_node.getLengthOfSC()
        # if(scL > max_sc_length): max_sc_length = scL
        # if(scL < min_sc_length): min_sc_length = scL
    # print "maximal SC length =",max_sc_length,"  min=",min_sc_length

    # print "===== Aperture Nodes START  ======="
    # aprtNodes = Add_quad_apertures_to_lattice(accLattice)
    # aprtNodes = Add_rfgap_apertures_to_lattice(accLattice,aprtNodes)
    # aprtNodes = AddMEBTChopperPlatesAperturesToSNS_Lattice(accLattice,aprtNodes)
    #
    # x_size = 0.042
    # y_size = 0.042
    # aprtNodes = AddScrapersAperturesToLattice(accLattice,"MEBT_Diag:H_SCRP",x_size,y_size,aprtNodes)
    #
    # x_size = 0.042
    # y_size = 0.042
    # aprtNodes = AddScrapersAperturesToLattice(accLattice,"MEBT_Diag:V_SCRP",x_size,y_size,aprtNodes)


    # for node in aprtNodes:
        # print "aprt=",node.getName()," pos =",node.getPosition()

    # print "===== Aperture Nodes Added ======="
