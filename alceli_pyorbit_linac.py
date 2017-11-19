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

import sys
import math
import random
import time

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

from alceli_linac_bunch_generator import ALCELI_Linac_BunchGenerator
from alceli_linac_lattice_factory import ALCELI_LinacLatticeFactory

def action_entrance(paramsDict):
    node = paramsDict["node"]
    bunch = paramsDict["bunch"]
    pos = paramsDict["path_length"]
    if(paramsDict["old_pos"] == pos): return
    if(paramsDict["old_pos"] + paramsDict["pos_step"] > pos): return
    paramsDict["old_pos"] = pos
    paramsDict["count"] += 1
    gamma = bunch.getSyncParticle().gamma()
    beta = bunch.getSyncParticle().beta()
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
    norm_emittX = emittX*gamma*beta
    norm_emittY = emittY*gamma*beta
    #---- phi_de_emittZ will be in [pi*deg*MeV]
    phi_de_emittZ = z_to_phase_coeff*emittZ
    eKin = bunch.getSyncParticle().kinEnergy()*1.0e+3

#     result into file
    s = " %35s %4.5f "%(node.getName(),pos+pos_start)
    s += "%6.4f %6.4f %6.4f %6.4f "%(alphaX,betaX,emittX,norm_emittX)
    s += "%6.4f %6.4f %6.4f %6.4f "%(alphaY,betaY,emittY,norm_emittY)
    s += "%6.4f %6.4f %6.4f %6.4f "%(alphaZ,betaZ,emittZ,phi_de_emittZ)
    s += "%5.3f %5.3f %5.3f "%(x_rms,y_rms,z_rms_deg)
    s += "%10.6f %8d"%(eKin,nParts)
    file_out.write(s +"\n")
    file_out.flush()

#     result to console buffer
    console_buffer_row=[]
    console_buffer_row.append(  '%5d'%paramsDict["count"])
    console_buffer_row.append(  '%25s'%node.getName())
    console_buffer_row.append( '%4.5f'%(pos+pos_start))
    console_buffer_row.append( '%5.3f'%x_rms)
    console_buffer_row.append( '%5.3f'%y_rms)
    console_buffer_row.append( '%5.3f'%z_rms_deg)
    console_buffer_row.append('%10.6f'%eKin)
    console_buffer_row.append(   '%8d'%nParts)
    console_buffer_rows.append(console_buffer_row)

def action_exit(paramsDict):
    action_entrance(paramsDict)

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

# =========================================================================================
random.seed(100)

# section list
names = ["HE"]

#---- the XML input file name with the linac structure
# xml_file_name = "./alceli.xml"
xml_file_name = "../SIMULINAC/25_09_2017_versuche_70_200MeV.xml"

#---- create the factory instance
alceli_linac_factory = ALCELI_LinacLatticeFactory()
alceli_linac_factory.setMaxDriftLength(0.01)

#---- make lattice from XML file
accLattice = alceli_linac_factory.getLinacAccLattice(names,xml_file_name)
print "Linac lattice is ready. L=",accLattice.getLength()

#----set up RF Gap Model -------------
#---- There are three available models at this moment
#---- BaseRfGap  uses only E0TL*cos(phi)*J0(kr) with E0TL = const
#---- MatrixRfGap uses a matrix approach like envelope codes
#---- RfGapTTF uses Transit Time Factors (TTF) like PARMILA
cppGapModel = BaseRfGap
#cppGapModel = MatrixRfGap
# cppGapModel = RfGapTTF
rf_gaps = accLattice.getRF_Gaps()
for rf_gap in rf_gaps:
    rf_gap.setCppGapModel(cppGapModel())

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

#-----TWISS Parameters at the entrance of MEBT ---------------
# transverse emittances are unnormalized and in pi*mm*mrad
# longitudinal emittance is in pi*eV*sec
e_kin_ini = 70.e-3 # in [GeV]
mass      = 0.939294    # in [GeV]
frequency = 816.e+6
v_light   = 2.99792458e+8  # in [m/sec]
gamma     = (mass + e_kin_ini)/mass
beta      = math.sqrt(gamma*gamma - 1.0)/gamma
print "relat. gamma=",gamma
print "relat.  beta=",beta

#------ emittances are normalized - transverse by gamma*beta and long. by gamma**3*beta
(alphaX,betaX,emittX) = (-1.9620, 0.1831, 0.21)
(alphaY,betaY,emittY) = ( 1.7681, 0.1620, 0.21)
(alphaZ,betaZ,emittZ) = ( 0.0196, 0.5844, 0.24153)
alphaZ = -alphaZ

#---make emittances un-normalized XAL units [m*rad]
emittX = 1.0e-6*emittX/(gamma*beta)
emittY = 1.0e-6*emittY/(gamma*beta)
emittZ = 1.0e-6*emittZ/(gamma**3*beta)
print " ========= XAL Twiss ==========="
print " aplha beta emitt[mm*mrad] X= %6.4f %6.4f %6.4f "%(alphaX,betaX,emittX*1.0e+6)
print " aplha beta emitt[mm*mrad] Y= %6.4f %6.4f %6.4f "%(alphaY,betaY,emittY*1.0e+6)
print " aplha beta emitt[mm*mrad] Z= %6.4f %6.4f %6.4f "%(alphaZ,betaZ,emittZ*1.0e+6)

#---- long. size in mm
sizeZ = math.sqrt(emittZ*betaZ)*1.0e+3

#---- transform to pyORBIT emittance[GeV*m]
emittZ = emittZ*gamma**3*beta**2*mass
betaZ = betaZ/(gamma**3*beta**2*mass)

print " ========= PyORBIT Twiss ==========="
print " aplha beta emitt[mm*mrad] X= %6.4f %6.4f %6.4f "%(alphaX,betaX,emittX*1.0e+6)
print " aplha beta emitt[mm*mrad] Y= %6.4f %6.4f %6.4f "%(alphaY,betaY,emittY*1.0e+6)
print " aplha beta emitt[mm*MeV] Z= %6.4f %6.4f %6.4f "%(alphaZ,betaZ,emittZ*1.0e+6)

twissX = TwissContainer(alphaX,betaX,emittX)
twissY = TwissContainer(alphaY,betaY,emittY)
twissZ = TwissContainer(alphaZ,betaZ,emittZ)

print "Start Bunch Generation",
bunch_gen = ALCELI_Linac_BunchGenerator(twissX,twissY,twissZ)

#set the initial kinetic energy in GeV
bunch_gen.setKinEnergy(e_kin_ini)

#set the beam peak current in mA
bunch_gen.setBeamCurrent(38.0)

bunch_in = bunch_gen.getBunch(nParticles = 10000, distributorClass = WaterBagDist3D)
#bunch_in = bunch_gen.getBunch(nParticles = 100000, distributorClass = GaussDist3D)
#bunch_in = bunch_gen.getBunch(nParticles = 10000, distributorClass = KVDist3D)
print " - finished"

print "Design tracking srated",
#set up design
accLattice.trackDesignBunch(bunch_in)
print " - finished"

#prepare bunch tracking through the lattice
paramsDict = {"old_pos":-1.,"count":0,"pos_step":0.1}
actionContainer = AccActionsContainer("Test Design Bunch Tracking")
pos_start = 0.
twiss_analysis = BunchTwissAnalysis()

file_out = open("pyorbit_twiss_sizes_ekin.dat","w")
s = "%35s %4s "%('Node','position')
s += "%6s %6s %6s %6s "%('alphaX','betaX','emittX','normEmittX')
s += "%6s %6s %6s %6s "%('alphaY','betaY','emittY','normEmittY')
s += "%6s %6s %6s %6s "%('alphaZ','betaZ','emittZ','emittZphiMeV')
s += "%5s %5s %5s "%('sizeX','sizeY','sizeZ_deg')
s += "%10s %8s "%('eKin','Nparts')
file_out.write(s+"\n")

# console buffer preparation
console_buffer_header = ['N','node','pos','sizeX','sizeY','sizeZ[deg]','eKin','Nparts']
console_buffer_rows   = []

# registration of actions
actionContainer.addAction(action_entrance, AccActionsContainer.ENTRANCE)
actionContainer.addAction(action_exit, AccActionsContainer.EXIT)

print "Bunch tracking started",
time_start = time.clock()
accLattice.trackBunch(bunch_in, paramsDict = paramsDict, actionContainer = actionContainer)
time_exec = time.clock() - time_start
file_out.close()
print " - finished in {:4.2f} [sec]".format(time_exec)

print "RESULTS - RESULTS - RESULTS - RESULTS - RESULTS - RESULTS - RESULTS - RESULTS - RESULTS - RESULTS "
print tblprnt(console_buffer_header,console_buffer_rows)
