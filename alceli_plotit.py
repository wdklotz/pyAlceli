# -*- coding: utf-8 -*-
import sys
import math
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import json

from alceli_conf import CONF

def display1(track_results):
   z=[]
   betaX=[]
   betaY=[]
   betaZ=[]
   Xrms=[]
   Yrms=[]
   Zrms_deg=[]
   emittX=[]
   emittXn=[]
   emittY=[]
   emittYn=[]
   viseo=[]     # the lattice elements plot
   zero=[]
   y1max=y2max=0.0

   for record in track_results:
      z.append(record['position'])
      betaX.append(record['betax'])
      betaY.append(record['betay'])
      betaZ.append(record['betaz']*1e-3)
      Xrms.append(record['xrms'])
      Yrms.append(record['yrms'])
      Zrms_deg.append(record['zrms_deg'])
      emittX.append(record['emittx'])
      emittXn.append(record['emittxn'])
      emittY.append(record['emitty'])
      emittYn.append(record['emittyn'])
      viseo.append(record['viseo'])
      zero.append(0)
      y1max = max(y1max,betaX[-1],betaY[-1])
      y2max = max(y2max,Xrms[-1],Yrms[-1])
   viseo1 = [0.2*y1max*x for x in viseo]     #scaling of lattice element plot
   viseo2 = [0.2*y2max*x for x in viseo]     #scaling of lattice element plot

   width=20; height=10
   plt.figure(CONF['title'],figsize=(width,height))
   splot = plt.subplot(221)
   splot.set_title('twiss: beta(X), beta(Y)')
   plt.plot(z,betaX,label=r'$\beta$-x [mm*mrad]',color='blue',linestyle='-')
   plt.plot(z,betaY,label=r'$\beta$-y [mm*mrad]',color='red',linestyle='-')
   plt.step(z,viseo1,color='black',linestyle='-')
   plt.plot(z,zero,color='black')
   plt.legend(loc='lower right',fontsize='large')

   splot = plt.subplot(222)
   splot.set_title('beam size(rms): X, Y')
   plt.plot(z,Xrms,label=r'$\sigma$-x [mm]',color='blue',linestyle='-')
   plt.plot(z,Yrms,label=r'$\sigma$-y [mm]',color='red',linestyle='-')
#    plt.step(z,viseo2,color='black',linestyle='-')
#    plt.plot(z,zero,color='black')
   plt.legend(loc='lower right',fontsize='large')

   splot = plt.subplot(223)
   splot.set_title('twiss: emittance(X), emittance(X)Norm')
   plt.plot(z,emittX,label=r'$\epsilon$-x  [mm*mrad]',color='blue',linestyle='-')
   plt.plot(z,emittXn,label=r'$\epsilon$-xN [mm*mrad]',color='red',linestyle='-')
   plt.legend(loc='lower right',fontsize='large')

   splot = plt.subplot(224)
   splot.set_title('long. size Z(rms) and beta(Z)')
   plt.plot(z,Zrms_deg,label=r'$\sigma$-z',color='blue',linestyle='-')
   plt.plot(z,betaZ,label=r'$\beta$-z',color='red',linestyle='-')
   plt.legend(loc='lower right',fontsize='large')

   plt.draw()

def make_scatter(axScatter,x,y,whazit):
#    max values
   xmax = np.max(np.fabs(x))
   ymax = np.max(np.fabs(y))

   # the scatter plot
   axScatter.scatter(x,y,s=1)

   # set tick label size
   axScatter.tick_params(labelsize='xx-small')

   # place a text box in upper left in axes coords
   props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)  # these are matplotlib.patch.Patch properties
   axScatter.text(0.05, 0.95, whazit, transform=axScatter.transAxes, fontsize=10, verticalalignment='top', bbox=props)

   #    create new axes on the right and on the top
   divider = make_axes_locatable(axScatter)
   axHistx = divider.append_axes('top',   1.2, pad=0.2, sharex=axScatter)
   axHisty = divider.append_axes('right', 1.2, pad=0.2, sharey=axScatter)
   axHistx.tick_params(labelsize='xx-small')
   axHisty.tick_params(labelsize='xx-small')

   # make some labels invisible
   plt.setp(axHistx.get_xticklabels() + axHisty.get_yticklabels(), visible=False)

   # now determine nice binning limits by hand
   binwidthx=xmax/100.
   binwidthy=ymax/100.
   limx = (int(xmax/binwidthx)+1)*binwidthx
   limy = (int(ymax/binwidthy)+1)*binwidthy
   binsx = np.arange(-limx, limx + binwidthx, binwidthx)
   binsy = np.arange(-limy, limy + binwidthy, binwidthy)

   # do the histograms
   axHistx.hist(x, bins=binsx)
   axHisty.hist(y, bins=binsy, orientation='horizontal')

#   axHistx.axis['bottom'].major_ticklabels.set_visible(False)
   for tl in axHistx.get_xticklabels():
      tl.set_visible(False)
   axHistx.set_yticks([0,50,100])

#   axHisty.axis['left'].major_ticklabels.set_visible(False)
   for tl in axHisty.get_yticklabels():
      tl.set_visible(False)
   axHisty.set_xticks([0,50,100])

def display2(bunch,whazit):
   x=[]; px=[]; y=[]; py=[];z=[]; pz=[]
   for item in bunch:
#       print item
      x.append (item['x']*1.e3)     #[mm]
      px.append(item['px']*1.e3)    #[mrad]?
      y.append (item['y']*1.e3)     #[mm]
      py.append(item['py']*1.e3)    #[mrad]?
      z.append (item['z']*1.e3)     #[mm]
      pz.append(item['pz']*1.e6)    #[??]

   width= 9.;   height = 9.
   fig = plt.figure(CONF['title']+", scatter plots@"+whazit,figsize=(width,height))
   ax1 = plt.subplot(221)
   make_scatter(ax1,x,px,'x,px')     #x,px
   ax2 = plt.subplot(222)
   make_scatter(ax2,y,py,'y,py')     #y,py
   ax3 = plt.subplot(223)
   make_scatter(ax3,x,y,'x,y')       #x,y
   ax4 = plt.subplot(224)
   make_scatter(ax4,z,pz,'z,pz')     #z,pz
   plt.draw()

## ============ MAIN ============
with open(CONF['plot_filename'],"r") as f:
   twiss_data = json.load(f)     # get the whole file in ram
   display1(twiss_data)

if CONF['dumpBunchIN']:
   bunch_in=[]
   with open(CONF['bunchIn_filename'],'r') as file:
      for line in file:
         if line[0] == '%':
            continue
         else:
            line = line[:-2].split(' ')
            bunch_in.append(
            dict(x=float(line[0]), px=float(line[1]), y=float(line[2]), py=float(line[3]), z=float(line[4]), pz=float(line[5])))
   display2(bunch_in,'ENTRANCE')

if CONF['dumpBunchOUT']:
   bunch_out=[]
   with open(CONF['bunchOut_filename'],'r') as file:
      for line in file:
         if line[0] == '%':
            continue
         else:
            line = line[:-2].split(' ')
            bunch_out.append(
            dict(x=float(line[0]), px=float(line[1]), y=float(line[2]), py=float(line[3]), z=float(line[4]), pz=float(line[5])))
   display2(bunch_out,'EXIT')

plt.show()

