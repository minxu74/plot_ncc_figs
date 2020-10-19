#!/usr/bin/env python

import numpy as np
import netCDF4 as nc4
import numpy.ma as ma

import matplotlib.pyplot as plt
import matplotlib.ticker as mtick


#dirname='/compyfs/dtn/yang954/fromCADES/PR-LUQ/181211_EVV_PR-LUQ_ICB20TRCNPRDCTCBC/run'
dirname='/compyfs/yang954/e3sm_scratch/20191208_testing_PR-LUQ_ICB20TRCNPRDCTCBC/run'

#casname='181211_EVV_PR-LUQ_ICB20TRCNPRDCTCBC'
casname='20191208_testing_PR-LUQ_ICB20TRCNPRDCTCBC'

fig, axs = plt.subplots(1, 3)


fig.set_size_inches(11, 8.5)
#plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
plt.ticklabel_format(style='sci', axis='x', scilimits=(0,0))

#axs[0].set_ylim(1.1e-8*86400,1.5e-8*86400)
#axs[1].set_ylim(1.1e-8*86400,1.5e-8*86400)
#axs[2].set_ylim(1.1e-8*86400,1.5e-8*86400)

#axs[0].set_ylim(0.3,0.6)
#axs[1].set_ylim(0.3,0.6)
#axs[2].set_ylim(0.3,0.6)

axs[2].tick_params(labelleft=False)
axs[1].tick_params(labelleft=False)

y3max=[]
y3min=[]
for iy in range(1850, 2003):
    print (iy)
    cy = f'{iy:04}'

    #for im in range(1, 13):
    #    cm = f'{im:02}'

    #    filename='181211_EVV_PR-LUQ_ICB20TRCNPRDCTCBC.clm2.h0.'+str(iy)+'-'+cm+'-01-00000.nc'
    #    with nc4.Dataset(dirname+'/'+filename, 'r') as f:
    #         print(f)
    filename=casname+'.clm2.h0.'+str(iy)+'-01-01-00000.nc'
    with nc4.Dataset(dirname+'/'+filename, 'r') as f:
         xx = f.variables['BIOCHEM_PMIN'][:,0] * 86400
         
         y1 = f.variables['TOTSOMP'][:,0] 
         y2 = f.variables['TOTSOMP_1m'][:,0]
         y3 = f.variables['SOLUTIONP'][:,0]

         ux = f.variables['BIOCHEM_PMIN'].units
         u1 = f.variables['TOTSOMP'].units
         u2 = f.variables['TOTSOMP_1m'].units
         u3 = f.variables['SOLUTIONP'].units

         #axs[0].xaxis.set_major_formatter(mtick.FormatStrFormatter('%.2e'))
         #axs[1].xaxis.set_major_formatter(mtick.FormatStrFormatter('%.2e'))
         #axs[2].xaxis.set_major_formatter(mtick.FormatStrFormatter('%.2e'))
#-
         axs[0].scatter(ma.average(y1), xx.sum(), s=1.6, label=cy)#, c, marker=verts)
         axs[1].scatter(ma.average(y2), xx.sum(), s=1.6, label=cy)#, c, marker=verts)
         axs[2].scatter(ma.average(y3), xx.sum(), s=1.6, label=cy)#, c, marker=verts)

         yy3=ma.average(y3)

         ux = "gP/m2/yr"
         axs[2].set_ylabel('BIOCHEM_PMIN ('+ux+')')
         axs[1].set_ylabel('BIOCHEM_PMIN ('+ux+')')
         axs[0].set_ylabel('BIOCHEM_PMIN ('+ux+')')

         #u3 = "gP/m2/day"
         #u2 = "gP/m2/day"
         #u1 = "gP/m2/day"
         axs[2].set_xlabel('SOLUTIONP ('+u3+')')
         axs[1].set_xlabel('TOTSOMP_1m ('+u2+')')
         axs[0].set_xlabel('TOTSOMP ('+u1+')')

         y3max.append(yy3.max())
         y3min.append(yy3.min())

print(min(y3min), max(y3max))
axs[2].set_xlim(min(y3min), max(y3max))
#axs[0].legend(scatterpoints=2, markerscale=3)


fig.savefig('test1.png', dpi=300)
