"""
The ALCELI Linac Lattice Factory generates the Linac Accelerator Lattice from the information
inside of the XML input file. This structure of this file is specific for the ALCELI.
Users from other facilities can use the same XML files with the same structure, but if they
need something else they can create their own Factory for different structure.
Here we use XmlDataAdaptor to parse the XML file.
The ALCELI Linac Lattice Factory uses a predefined set of Linac Acc Elements.
"""

import os
import sys
import string
import math

# import the XmlDataAdaptor XML parser
from orbit.utils.xml import XmlDataAdaptor

from orbit.py_linac.lattice import LinacAccLattice
from orbit.py_linac.lattice import LinacAccNodes

from orbit.py_linac.lattice import BaseLinacNode, LinacNode, LinacMagnetNode, MarkerLinacNode, Drift, Quad, AbstractRF_Gap, Bend
from orbit.py_linac.lattice import DCorrectorH, DCorrectorV
from orbit.py_linac.lattice import RF_Cavity, Sequence
from orbit.py_linac.lattice import BaseRF_Gap

# import general accelerator elements
from orbit.lattice import AccNode

# import pyORBIT Python utilities classes for objects with names, types, and dictionary parameters
from orbit.utils import orbitFinalize

# DEBUG
from acDebugHelpers import caller_name, lineno, DEBUG_ON, DEBUG_OFF, DEXIT

DEBUG_FACTORY = DEBUG_OFF

class AcLinacLatticeFactory():
   """
   The ALCELI Linac Lattice Factory generates the Linac Accelerator Lattice
   from the XML file of the specific structure.
   """
   def __init__(self):
      #We need to compare positions, lengths etc. This is our delta in [m]
      self.zeroDistance = 0.00001
      #The maximal length of the drift. It will be devided if it is more than that.
      self.maxDriftLength = 1.

   def setMaxDriftLength(self, maxDriftLength = 1.0):
      """
      Sets the maximal drift length that is used for
      the purpose of the space charge calculations and diagnostics.
      """
      self.maxDriftLength = maxDriftLength
      # DEBUG_FACTORY(__file__,lineno(),'caller name: {}'.format(caller_name(skip=2)))

   def getMaxDriftLength(self):
      """
      Returns the maximal drift length that is used for the purpose
      of the space charge calculations and diagnostics.
      """
      return self.maxDriftLength

   def getLinacAccLattice(self,names,xml_file_name):
      """
      Returns the linac accelerator lattice for specified sequence names and for a specified XML file.
      """
      if(len(names) < 1):
         msg = "The AcLinacLatticeFactory method getLinacAccLattice(names,xml_file_name): you have to specify the names array!"
         msg = msg + os.linesep
         msg = msg + "Stop."
         msg = msg + os.linesep
         orbitFinalize(msg)
      #----- let's parse the XML file
      acc_da = XmlDataAdaptor.adaptorForFile(xml_file_name)
      # DEBUG_FACTORY(__file__,lineno(),acc_da)
      # DEBUG_FACTORY(__file__,lineno(),acc_da.__dict__)
      lattice = self.getLinacAccLatticeFromDA(names,acc_da)
      return (lattice,acc_da)

   def getLinacAccLatticeFromDA(self,names,acc_da):
      """
      Returns the linac accelerator lattice for specified sequence names.
      """
      if(len(names) < 1):
         msg = "The AcLinacLatticeFactory method getLinacAccLatticeFromDA(names,): you have to specify the names array!"
         msg = msg + os.linesep
         msg = msg + "Stop."
         msg = msg + os.linesep
         orbitFinalize(msg)
      # DEBUG_FACTORY(__file__,lineno(),acc_da.getName())

      #----- let's parse the XML DataAdaptor
      accSeq_da_arr = acc_da.childAdaptors()
      # DEBUG_FACTORY(__file__,lineno(),accSeq_da_arr[0].makeXmlText())

      #-----let's filter and check that the names in good order
      accSeq_da_arr = self.filterSequences_and_OptionalCheck(accSeq_da_arr,names)
      # DEBUG_FACTORY(__file__,lineno(),accSeq_da_arr)
      # DEBUG_FACTORY(__file__,lineno(),string.join(['{}'.format(i.getAttributes()) for i in accSeq_da_arr]))

      #----make linac latticeaccSeq
      linacAccLattice = LinacAccLattice(acc_da.getName())
      # DEBUG_FACTORY(__file__,lineno(),linacAccLattice.__dict__)

      #There are the folowing possible types of elements in the linac tree:
      #QUAD - quadrupole
      #RFGAP - RF Gap
      #DCH - horizontal dipole corrector
      #DCV - vertical dipole corrector
      #Marker - anything else with the length equals to 0
      #Before putting enerything into the linacAccLattice we will create sequences
      # with all nodes.
      #----------------------------------------------------------------------
      # The DRIFTS will be generated additionally and put into right places
      #----------------------------------------------------------------------
      def positionComp(node1_da,node2_da):
         if(node1_da.getParam("pos") > node2_da.getParam("pos")):
            return 1
         else:
            if(node1_da.getParam("pos") == node2_da.getParam("pos")):
               return 0
         return -1

      # loop sequences
      accSeqs = []
      accRF_Cavs = []
      seqPosition = 0.
      for seq_da in accSeq_da_arr:
         accSeq = Sequence(seq_da.getName())
         # DEBUG_FACTORY(__file__,lineno(),'seq: {}'.format(accSeq.getName()))
         accSeq.setLinacAccLattice(linacAccLattice)
         accSeq.setLength(seq_da.doubleValue("length"))
         accSeq.setPosition(seqPosition)
         seqPosition = seqPosition + accSeq.getLength()
         #---- BPM frequnecy for this sequence ----
         if(seq_da.hasAttribute("bpmFrequency")):
            bpmFrequency = seq_da.doubleValue("bpmFrequency")
            accSeq.addParam("bpmFrequency",bpmFrequency)
         #-----------------------------------------
         accSeqs.append(accSeq)
         # DEBUG_FACTORY(__file__,lineno(),'sequences: '+string.join(['{}'.format(i.getName()) for i in accSeqs]))

         #---- create RF Cavities
         if(len(seq_da.childAdaptors("Cavities")) == 1):
            cavs_da = seq_da.childAdaptors("Cavities")[0]
            cav_da_arr = cavs_da.childAdaptors("Cavity")
            # DEBUG_FACTORY(__file__,lineno(),'Cavities: '+string.join(['{}'.format(i.stringValue('name')) for i in cav_da_arr]))
            # loop cavities
            for cav_da in cav_da_arr:
               frequency = cav_da.doubleValue("frequency")
               cav_amp = cav_da.doubleValue("ampl")
               cav_name = cav_da.stringValue("name")
               cav_pos = cav_da.doubleValue("pos")
               cav = RF_Cavity(cav_name)
               cav.setAmp(cav_amp)
               cav.setFrequency(frequency)
               cav.setPosition(cav_pos)
               accSeq.addRF_Cavity(cav)
               # DEBUG_FACTORY(__file__,lineno(),cav.__dict__)
         #----------------------------
         #node_da_arr - array of accElements. These nodes are not AccNodes. They are XmlDataAdaptor class instances
         node_da_arr = seq_da.childAdaptors("accElement")
         #put nodes in order according to the position in the sequence
         for node_da in node_da_arr:
            node_da.setParam("pos",node_da.doubleValue("pos"))
         node_da_arr.sort(positionComp)
         #thinNodes - array of accNode nodes with zero length
         #They can be positioned inside the thick nodes, and this will be done at the end of this method
         thinNodes = []
         
         # node  loop
         for node_da in node_da_arr:
            # DEBUG_FACTORY(__file__, lineno(), node_da.__dict__)
            params_da    = node_da.childAdaptors("parameters")[0]
            # DEBUG_FACTORY(__file__, lineno(), params_da.__dict__)
            node_tagname = node_da.getName()
            node_name    = node_da.stringValue('name')
            node_type    = node_da.stringValue("type")
            node_length  = node_da.doubleValue("length")
            node_pos     = node_da.getParam("pos")
            # DEBUG_FACTORY(__file__,lineno(),'node_da: {} {} {} len= {} pos={}'.format( node_tagname,node_name,node_type,node_length,node_pos))
            #------------QUAD-----------------
            if(node_type == "QUAD"):
               accNode = Quad(node_da.stringValue("name"))
               accNode.setParam("dB/dr",params_da.doubleValue("field"))
               accNode.setParam("field",params_da.doubleValue("field"))
               accNode.setLength(node_length)
               if(params_da.hasAttribute("poles")):
                  #accNode.setParam("poles",[int(x) for x in eval(params_da.stringValue("poles"))])
                  accNode.setParam("poles",params_da.intArrayValue("poles"))
               if(params_da.hasAttribute("kls")):
                  #accNode.setParam("kls", [x for x in eval(params_da.stringValue("kls"))])
                  accNode.setParam("kls",params_da.doubleArrayValue("kls"))
               if(params_da.hasAttribute("skews")):
                  #accNode.setParam("skews",[int(x) for x in eval(params_da.stringValue("skews"))])
                  accNode.setParam("skews",params_da.intArrayValue("skews"))
               if(0.5*accNode.getLength() > self.maxDriftLength):
                  accNode.setnParts(2*int(0.5*accNode.getLength()/self.maxDriftLength  + 1.5 - 1.0e-12) )
               if(params_da.hasAttribute("aperture") and params_da.hasAttribute("aprt_type")):
                  accNode.setParam("aprt_type",params_da.intValue("aprt_type"))
                  accNode.setParam("aperture",params_da.doubleValue("aperture"))
               #---- possible parameters for PMQ description of the in Trace3D style
               if(params_da.hasAttribute("radIn") and params_da.hasAttribute("radOut")):
                  accNode.setParam("radIn",params_da.doubleValue("radIn"))
                  accNode.setParam("radOut",params_da.doubleValue("radOut"))
               accNode.setParam("pos",node_pos)
               accSeq.addNode(accNode)
               # DEBUG_FACTORY(__file__,lineno(),'maxDriftLength {}'.format(self.maxDriftLength))
               # DEBUG_FACTORY(__file__,lineno(),accNode.__dict__)
               
            #------------BEND-----------------
            elif(node_type == "BEND"):
               accNode = Bend(node_da.stringValue("name"))
               if(params_da.hasAttribute("poles")):
                  #accNode.setParam("poles",[int(x) for x in eval(params_da.stringValue("poles"))])
                  accNode.setParam("poles",params_da.intArrayValue("poles"))
               if(params_da.hasAttribute("kls")):
                  #accNode.setParam("kls", [x for x in eval(params_da.stringValue("kls"))])
                  accNode.setParam("kls",params_da.doubleArrayValue("kls"))
               if(params_da.hasAttribute("skews")):
                  #accNode.setParam("skews",[int(x) for x in eval(params_da.stringValue("skews"))])
                  accNode.setParam("skews",params_da.intArrayValue("skews"))
               accNode.setParam("ea1",params_da.doubleValue("ea1"))
               accNode.setParam("ea2",params_da.doubleValue("ea2"))
               accNode.setParam("theta",params_da.doubleValue("theta"))
               if(params_da.hasAttribute("aperture_x") and params_da.hasAttribute("aperture_y") and params_da.hasAttribute("aprt_type")):
                  accNode.setParam("aprt_type",params_da.intValue("aprt_type"))
                  accNode.setParam("aperture_x",params_da.doubleValue("aperture_x"))
                  accNode.setParam("aperture_y",params_da.doubleValue("aperture_y"))
               accNode.setLength(node_length)
               if(accNode.getLength() > self.maxDriftLength):
                  accNode.setnParts(2*int(accNode.getLength()/self.maxDriftLength  + 1.5 - 1.0e-12))
               accNode.setParam("pos",node_pos)
               accSeq.addNode(accNode)
            #------------RF_Gap-----------------
            elif(node_type == "RFGAP"):
               accNode = BaseRF_Gap(node_da.stringValue("name"))
               accNode.setLength(0.)
               accNode.setParam("E0TL",params_da.doubleValue("E0TL"))
               accNode.setParam("E0L",params_da.doubleValue("E0L"))
               accNode.setParam("mode",params_da.doubleValue("mode"))
               accNode.setParam("gap_phase",params_da.doubleValue("phase")*math.pi/180.)
               accNode.setParam("EzFile",params_da.stringValue("EzFile"))
               cav_name = params_da.stringValue("cavity")
               cav = accSeq.getRF_Cavity(cav_name)
               cav.addRF_GapNode(accNode)
               if(accNode.isFirstRFGap()):
                  cav.setPhase(accNode.getParam("gap_phase"))
               #---- TTFs parameters
               ttfs_da = node_da.childAdaptors("TTFs")[0]
               accNode.setParam("beta_min",ttfs_da.doubleValue("beta_min"))
               accNode.setParam("beta_max",ttfs_da.doubleValue("beta_max"))
               (polyT,polyS,polyTp,polySp) = accNode.getTTF_Polynimials()
               polyT_da = ttfs_da.childAdaptors("polyT")[0]
               polyS_da = ttfs_da.childAdaptors("polyS")[0]
               polyTp_da = ttfs_da.childAdaptors("polyTP")[0]
               polySp_da = ttfs_da.childAdaptors("polySP")[0]
               polyT.order(polyT_da.intValue("order"))
               polyS.order(polyS_da.intValue("order"))
               polyTp.order(polyTp_da.intValue("order"))
               polySp.order(polySp_da.intValue("order"))
               coef_arr = polyT_da.doubleArrayValue("pcoefs")
               for coef_ind in range(len(coef_arr)):
                  polyT.coefficient(coef_ind,coef_arr[coef_ind])
               coef_arr = polyS_da.doubleArrayValue("pcoefs")
               for coef_ind in range(len(coef_arr)):
                  polyS.coefficient(coef_ind,coef_arr[coef_ind])
               coef_arr = polyTp_da.doubleArrayValue("pcoefs")
               for coef_ind in range(len(coef_arr)):
                  polyTp.coefficient(coef_ind,coef_arr[coef_ind])
               coef_arr = polySp_da.doubleArrayValue("pcoefs")
               for coef_ind in range(len(coef_arr)):
                  polySp.coefficient(coef_ind,coef_arr[coef_ind])
               if(params_da.hasAttribute("aperture") and params_da.hasAttribute("aprt_type")):
                  accNode.setParam("aprt_type",params_da.intValue("aprt_type"))
                  accNode.setParam("aperture",params_da.doubleValue("aperture"))
               accNode.setParam("pos",node_pos)
               accSeq.addNode(accNode)
               # DEBUG_FACTORY(__file__,lineno(),'accNode: '+string.join(['\n\t{} : {}'.format(k,v) for k,v in accNode.__dict__.items()]))
            else:
               if(node_length != 0.):
                  msg = "The LinacLatticeFactory method getLinacAccLattice(names): there is a strange element!"
                  msg = msg + os.linesep
                  msg = msg + "name=" + node_da.stringValue("name")
                  msg = msg + os.linesep
                  msg = msg + "type="+node_type
                  msg = msg + os.linesep
                  msg = msg + "length(should be 0.)="+str(node_length)
                  orbitFinalize(msg)
               #------ thin nodes analysis
               accNode = None
               if(node_type == "DCV" or node_type == "DCH"):
                  if(node_type == "DCV"): accNode = DCorrectorV(node_da.stringValue("name"))
                  if(node_type == "DCH"): accNode = DCorrectorH(node_da.stringValue("name"))
                  accNode.setParam("effLength",params_da.doubleValue("effLength"))
                  if(params_da.hasAttribute("B")):
                     accNode.setParam("B",params_da.doubleValue("B"))
               else:
                  accNode = MarkerLinacNode(node_da.stringValue("name"))
               accNode.setParam("pos",node_pos)
               thinNodes.append(accNode)
         #----- assign the thin nodes that are inside the thick nodes
         unusedThinNodes = []
         for thinNode in thinNodes:
            thinNode_pos = thinNode.getParam("pos")
            isInside = False
            for accNode in accSeq.getNodes():
               length = accNode.getLength()
               if(length > 0.):
                  pos = accNode.getParam("pos")
                  if(thinNode_pos >= (pos-length/2) and thinNode_pos <= (pos+length/2)):
                     isInside = True
                     delta_pos = thinNode_pos - (pos-length/2)
                     s_path = 0.
                     part_ind_in = -1
                     for part_ind in range(accNode.getnParts()):
                        part_ind_in = part_ind
                        s_path += accNode.getLength(part_ind)
                        if(delta_pos <= s_path + self.zeroDistance):
                           break
                     accNode.addChildNode(thinNode, place = AccNode.BODY, part_index = part_ind_in , place_in_part = AccNode.AFTER)
                     thinNode.setParam("pos",(pos-length/2)+s_path)
            if(not isInside):
               unusedThinNodes.append(thinNode)
         thinNodes = unusedThinNodes
         newAccNodes = accSeq.getNodes()[:] + thinNodes
         newAccNodes.sort(positionComp)
         accSeq.setNodes(newAccNodes)

         #insert the drifts ======================start ===========================
         #-----now check the integrity quads and rf_gaps should not overlap
         #-----and create drifts
         copyAccNodes = accSeq.getNodes()[:]
         # DEBUG_FACTORY(__file__,lineno(),copyAccNodes)
         firstNode = copyAccNodes[0]
         lastNode = copyAccNodes[len(copyAccNodes)-1]
         driftNodes_before = []
         driftNodes_after = []
         #insert the drift before the first element if its half length is less than its position
         if(math.fabs(firstNode.getLength()/2.0 - firstNode.getParam("pos")) > self.zeroDistance):
            if(firstNode.getLength()/2.0 > firstNode.getParam("pos")):
               msg = "The LinacLatticeFactory method getLinacAccLattice(names): the first node is too long!"
               msg = msg + os.linesep
               msg = msg + "name=" + firstNode.getName()
               msg = msg + os.linesep
               msg = msg + "type=" + firstNode.getType()
               msg = msg + os.linesep
               msg = msg + "length=" + str(firstNode.getLength())
               msg = msg + os.linesep
               msg = msg + "pos=" + str(firstNode.getParam("pos"))
               orbitFinalize(msg)
            else:
               driftNodes = []
               driftLength = firstNode.getParam("pos") - firstNode.getLength()/2.0
               nDrifts = int(driftLength/self.maxDriftLength) + 1
               driftLength = driftLength/nDrifts
               for idrift in range(nDrifts):
                  drift = Drift(accSeq.getName()+":START:"+str(idrift+1)+":drift")
                  drift.setLength(driftLength)
                  drift.setParam("pos",0.+drift.getLength()*(idrift+0.5))
                  driftNodes.append(drift)
               driftNodes_before = driftNodes
         #insert the drift after the last element if its half length less + position is less then the sequence length
         if(math.fabs(lastNode.getLength()/2.0 + lastNode.getParam("pos") - accSeq.getLength()) > self.zeroDistance):
            if(lastNode.getLength()/2.0 + lastNode.getParam("pos") > accSeq.getLength()):
               msg = "The LinacLatticeFactory method getLinacAccLattice(names): the last node is too long!"
               msg = msg + os.linesep
               msg = msg + "name=" + lastNode.getName()
               msg = msg + os.linesep
               msg = msg + "type=" + lastNode.getType()
               msg = msg + os.linesep
               msg = msg + "length=" + str(lastNode.getLength())
               msg = msg + os.linesep
               msg = msg + "pos=" + str(lastNode.getParam("pos"))
               msg = msg + os.linesep
               msg = msg + "sequence name=" + accSeq.getName()
               msg = msg + os.linesep
               msg = msg + "sequence length=" + str(accSeq.getLength())
               orbitFinalize(msg)
            else:
               driftNodes = []
               driftLength = accSeq.getLength() - (lastNode.getParam("pos") + lastNode.getLength()/2.0)
               nDrifts = int(driftLength/self.maxDriftLength) + 1
               driftLength = driftLength/nDrifts
               for idrift in range(nDrifts):
                  drift = Drift(accSeq.getName()+":"+lastNode.getName()+":"+str(idrift+1)+":drift")
                  drift.setLength(driftLength)
                  drift.setParam("pos",lastNode.getParam("pos")+lastNode.getLength()/2.0 + drift.getLength()*(idrift+0.5))
                  driftNodes.append(drift)
               driftNodes_after = driftNodes
         #now move on and generate drifts between (i,i+1) nodes from copyAccNodes
         newAccNodes = driftNodes_before
         for node_ind in range(len(copyAccNodes)-1):
            accNode0 = copyAccNodes[node_ind]
            newAccNodes.append(accNode0)
            accNode1 = copyAccNodes[node_ind+1]
            dist = accNode1.getParam("pos") - accNode1.getLength()/2 - (accNode0.getParam("pos") + accNode0.getLength()/2)
            # DEBUG_FACTORY(__file__,lineno(),'distance from {},to {}, {:8.4f}[m]'.format(accNode0.getName(),accNode1.getName(),dist))
            if(abs(dist)<1.e-10): dist = 0.
            if(dist < 0.):
               msg = "The LinacLatticeFactory method getLinacAccLattice(names): two nodes are overlapping!"
               msg = msg + os.linesep
               msg = msg + "sequence name=" + accSeq.getName()
               msg = msg + os.linesep
               msg = msg + "node 0 name=" + accNode0.getName() + " pos="+ str(accNode0.getParam("pos")) + " L="+str(accNode0.getLength())
               msg = msg + os.linesep
               msg = msg + "node 1 name=" + accNode1.getName() + " pos="+ str(accNode1.getParam("pos")) + " L="+str(accNode1.getLength())
               msg = msg + os.linesep
               orbitFinalize(msg)
            elif(dist > self.zeroDistance):
               driftNodes = []
               nDrifts = int(dist/self.maxDriftLength) + 1
               driftLength = dist/nDrifts
               for idrift in range(nDrifts):
                  drift = Drift(accSeq.getName()+":"+accNode0.getName()+":"+str(idrift+1)+":drift")
                  drift.setLength(driftLength)
                  drift.setParam("pos",accNode0.getParam("pos")+accNode0.getLength()*0.5+drift.getLength()*(idrift+0.5))
                  driftNodes.append(drift)
               newAccNodes += driftNodes
            else:
               pass
         newAccNodes.append(lastNode)
         newAccNodes += driftNodes_after
         accSeq.setNodes(newAccNodes)
         #insert the drifts ======================stop ===========================
         #add all AccNodes to the linac lattice
         for accNode in accSeq.getNodes():
            linacAccLattice.addNode(accNode)
      #------- finalize the lattice construction
      linacAccLattice.initialize()
      nodes = linacAccLattice.getNodes()
      if DEBUG_FACTORY == DEBUG_ON:
         DEBUG_FACTORY(__file__,lineno(),'LinacAccLattice initalized')
         pass
         # for i in range(0,len(nodes)-1):
         #    l  = nodes[i].getLength()
         #    s0 = nodes[i].getParam("pos")
         #    si = s0-l/2.
         #    sf = s0+l/2.
         #    print('{} \t(si,s0,sf) ({},{},{})'.format(nodes[i],si,s0,sf))
      return linacAccLattice

   def filterSequences_and_OptionalCheck(self,accSeq_da_arr,names):
      """
      This method will filter the sequences according to names list
      and check the order of sequences in names.
      All sequences should be in the right order. For ALCELI linac is
      just a linac. It returns the filtered array with data adapters
      with the names in the names array.
      """
      seqencesLocal = accSeq_da_arr
      seqencesLocalNames = []
      for seq_da in seqencesLocal:
         seqencesLocalNames.append(seq_da.getName())
      ind_old = -1
      count = 0
      for name in names:
         ind = seqencesLocalNames.index(name)
         if(ind < 0 or (count > 0 and ind != (ind_old + 1))):
            msg = "The LinacLatticeFactory method getLinacAccLattice(names): sequence names array is wrong!"
            msg = msg + os.linesep
            msg = msg + "existing names=" + str(seqencesLocalNames)
            msg = msg + os.linesep
            msg = msg + "sequence names="+str(names)
            orbitFinalize(msg)
         ind_old = ind
         count += 1
      ind_start = seqencesLocalNames.index(names[0])
      accSeq_da_arr = accSeq_da_arr[ind_start:ind_start+len(names)]
      return accSeq_da_arr

   def makeDataAdaptorforLinacLattice(self,linacAccLattice):
      """
      This method generates and returns the accelerator data adaptor
      for the accelerator lattice.
      """
      acc_da = XmlDataAdaptor(linacAccLattice.getName())
      accSeqs = linacAccLattice.getSequences()
      for accSeq in accSeqs:
         self.makeDAforAccSeq(acc_da,accSeq)
      return acc_da

   def makeDAforAccSeq(self,acc_da,accSeq):
      """
      This method adds the accelerator sequence data to accelerator data adaptor.
      """
      seq_da =  acc_da.createChild(accSeq.getName())
      seq_da.setValue("length",accSeq.getLength())
      if(accSeq.hasParam("bpmFrequency")):
         seq_da.setValue("bpmFrequency",accSeq.getParam("bpmFrequency"))
      #----------------------------------
      nodes = accSeq.getNodes()
      for node in nodes:
         if(isinstance(node,Quad)):
            self.makeDA_quad(seq_da,node)
            continue
         if(isinstance(node,Bend)):
            self.makeDA_bend(seq_da,node)
            continue
         if(isinstance(node,DCorrectorH)):
            self.makeDA_dcorrH(seq_da,node)
            continue
         if(isinstance(node,DCorrectorV)):
            self.makeDA_dcorrV(seq_da,node)
            continue
         if(isinstance(node,BaseRF_Gap)):
            self.makeDA_rf_gap(seq_da,node)
            continue
         if(isinstance(node,MarkerLinacNode)):
            self.makeDA_marker(seq_da,node)
            continue
      #----------------------------------
      rf_cavs = accSeq.getRF_Cavities()
      cavs_da = seq_da.createChild("Cavities")
      for cav in rf_cavs:
         cav_da = cavs_da.createChild("Cavity")
         cav_da.setValue("ampl",cav.getAmp())
         cav_da.setValue("frequency",cav.getFrequency())
         cav_da.setValue("name",cav.getName())
         cav_da.setValue("pos",cav.getPosition())

   def makeDA_quad(self,seq_da,quad):
      """
      This method will fill out the DataAdaptor for a quad linac accelerator node.
      """
      quad_da =  seq_da.createChild("accElement")
      params_da = quad_da.createChild("parameters")
      quad_da.setValue("length",quad.getLength())
      quad_da.setValue("name",quad.getName())
      quad_da.setValue("pos",("%9.6f"%quad.getPosition()).strip())
      quad_da.setValue("type","QUAD")
      params_da.setValue("field",quad.getParam("dB/dr"))
      if(quad.hasParam("poles") and len(quad.getParam("poles")) > 0):
         params_da.setValue("poles",quad.getParam("poles"))
         params_da.setValue("kls",quad.getParam("kls"))
         params_da.setValue("skews",quad.getParam("skews"))
      if(quad.hasParam("aprt_type")):
         params_da.setValue("aprt_type",quad.getParam("aprt_type"))
         params_da.setValue("aperture",quad.getParam("aperture"))

   def makeDA_bend(self,seq_da,bend):
      """
      This method will fill out the DataAdaptor for a bend linac accelerator node.
      """
      bend_da =  seq_da.createChild("accElement")
      params_da = bend_da.createChild("parameters")
      bend_da.setValue("type","BEND")
      bend_da.setValue("name",bend.getName())
      bend_da.setValue("pos",("%9.6f"%bend.getPosition()).strip())
      bend_da.setValue("length",bend.getLength())
      params_da.setValue("theta",bend.getParam("theta"))
      if(bend.hasParam("ea1")):
         params_da.setValue("ea1",bend.getParam("ea1"))
      if(bend.hasParam("ea2")):
         params_da.setValue("ea2",bend.getParam("ea2"))
      if(bend.hasParam("poles") and len(bend.getParam("poles")) > 0):
         params_da.setValue("poles",bend.getParam("poles"))
         params_da.setValue("kls",bend.getParam("kls"))
         params_da.setValue("skews",bend.getParam("skews"))
      if(bend.hasParam("aprt_type")):
         params_da.setValue("aprt_type",bend.getParam("aprt_type"))
         params_da.setValue("aperture_x",bend.getParam("aperture_x"))
         params_da.setValue("aperture_y",bend.getParam("aperture_y"))

   def makeDA_dcorrH(self,seq_da,dcorr):
      """
      This method will fill out the DataAdaptor for a DCorrectorH linac accelerator node.
      """
      dcorr_da =  seq_da.createChild("accElement")
      params_da = dcorr_da.createChild("parameters")
      dcorr_da.setValue("length",0.)
      dcorr_da.setValue("name",dcorr.getName())
      dcorr_da.setValue("pos",("%9.6f"%dcorr.getPosition()).strip())
      dcorr_da.setValue("type","DCH")
      params_da.setValue("B",dcorr.getParam("B"))
      params_da.setValue("effLength",dcorr.getParam("effLength"))

   def makeDA_dcorrV(self,seq_da,dcorr):
      """
      This method will fill out the DataAdaptor for a DCorrectorV linac accelerator node.
      """
      dcorr_da =  seq_da.createChild("accElement")
      params_da = dcorr_da.createChild("parameters")
      dcorr_da.setValue("length",0.)
      dcorr_da.setValue("name",dcorr.getName())
      dcorr_da.setValue("pos",("%9.6f"%dcorr.getPosition()).strip())
      dcorr_da.setValue("type","DCV")
      params_da.setValue("B",dcorr.getParam("B"))
      params_da.setValue("effLength",dcorr.getParam("effLength"))

   def makeDA_marker(self,seq_da,marker):
      marker_da =  seq_da.createChild("accElement")
      params_da = marker_da.createChild("parameters")
      marker_da.setValue("length",0.)
      marker_da.setValue("name",marker.getName())
      marker_da.setValue("pos",("%9.6f"%marker.getPosition()).strip())
      marker_da.setValue("type","MARKER")

   def makeDA_rf_gap(self,seq_da,rf_gap):
      """
      This method will fill out the DataAdaptor for a RF Gap linac accelerator node.
      """
      rf_gap_da =  seq_da.createChild("accElement")
      rf_gap_da.setValue("length",0.)
      rf_gap_da.setValue("name",rf_gap.getName())
      rf_gap_da.setValue("pos",("%9.6f"%rf_gap.getPosition()).strip())
      rf_gap_da.setValue("type","RFGAP")
      params_da = rf_gap_da.createChild("parameters")
      params_da.setValue("E0L",rf_gap.getParam("E0L"))
      params_da.setValue("E0TL",rf_gap.getParam("E0TL"))
      params_da.setValue("EzFile",rf_gap.getParam("EzFile"))
      params_da.setValue("cavity",rf_gap.getRF_Cavity().getName())
      params_da.setValue("mode",int(rf_gap.getParam("mode")))
      params_da.setValue("phase","%7.3f"%(rf_gap.getParam("gap_phase")*180.0/math.pi))
      if(rf_gap.hasParam("aprt_type")):
         params_da.setValue("aprt_type",rf_gap.getParam("aprt_type"))
         params_da.setValue("aperture",rf_gap.getParam("aperture"))
      rf_gap_ttfs_da = rf_gap_da.createChild("TTFs")
      rf_gap_ttfs_da.setValue("beta_max",rf_gap.getParam("beta_max"))
      rf_gap_ttfs_da.setValue("beta_min",rf_gap.getParam("beta_min"))
      #---------------------------------------------
      def make_poly_da(node_name,ttfs_da,poly):
         poly_da = ttfs_da.createChild(node_name)
         poly_da.setValue("order",poly.order())
         coefs = []
         for ind in range(poly.order()+1):
            coefs.append(poly.coefficient(ind))
         poly_da.setValue("pcoefs",coefs)
      #---------------------------------------------
      (polyT,polyS,polyTp,polySp) = rf_gap.getTTF_Polynimials()
      make_poly_da("polyT",rf_gap_ttfs_da,polyT)
      make_poly_da("polyS",rf_gap_ttfs_da,polyS)
      make_poly_da("polyTP",rf_gap_ttfs_da,polyTp)
      make_poly_da("polySP",rf_gap_ttfs_da,polySp)


