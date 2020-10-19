#!/usr/bin/env python


import matplotlib.pyplot as plt
import netCDF4 as nc4
import numpy as np
import matplotlib.cm as cmp
import matplotlib.colors as mcolors
import cartopy.crs as ccrs
import numpy.ma as ma



def hex_to_rgb(value):
    '''
    Converts hex to rgb colours
    value: string of 6 characters representing a hex colour.
    Returns: list length 3 of RGB values'''
    value = value.strip("#") # removes hash symbol if present
    lv = len(value)
    return tuple(int(value[i:i + lv // 3], 16) for i in range(0, lv, lv // 3))


def rgb_to_dec(value):
    '''
    Converts rgb to decimal colours (i.e. divides each value by 256)
    value: list (length 3) of RGB values
    Returns: list (length 3) of decimal values'''
    return [v/256 for v in value]


def get_continuous_cmap(hex_list, float_list=None):
    ''' creates and returns a color map that can be used in heat map figures.
        If float_list is not provided, colour map graduates linearly between each color in hex_list.
        If float_list is provided, each color in hex_list is mapped to the respective location in float_list. 
        
        Parameters
        ----------
        hex_list: list of hex code strings
        float_list: list of floats between 0 and 1, same length as hex_list. Must start with 0 and end with 1.
        
        Returns
        ----------
        colour map'''
    rgb_list = [rgb_to_dec(hex_to_rgb(i)) for i in hex_list]
    if float_list:
        pass
    else:
        float_list = list(np.linspace(0,1,len(rgb_list)))
        
    cdict = dict()
    for num, col in enumerate(['red', 'green', 'blue']):
        col_list = [[float_list[i], rgb_list[i][num], rgb_list[i][num]] for i in range(len(float_list))]
        cdict[col] = col_list
    cmp = mcolors.LinearSegmentedColormap('my_cmp', segmentdata=cdict, N=256)
    return cmp


CasNam1_CNP = "/compyfs/yang954/e3sm_scratch/20190912_hcru_hcru_ICB20TRCNPRDCTCBC/run/20190912_hcru_hcru_ICB20TRCNPRDCTCBC.clm2.h0."
CasNam2_CNP = "/compyfs/yang954/e3sm_scratch/EXP2CO2_20190912_hcru_hcru_ICB20TRCNPRDCTCBC/run/EXP2CO2_20190912_hcru_hcru_ICB20TRCNPRDCTCBC.clm2.h0."

DaysPerMon = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]


Biomass_ABG1 = ma.zeros((360, 720))
Biomass_ABG2 = ma.zeros((360, 720))

LAI1 = ma.zeros((360, 720))
LAI2 = ma.zeros((360, 720))

for iy in range(2000, 2010, 1):
    cy = '{:04d}'.format(iy)
    for im in range(1, 13, 1):
        cm = '{:02d}'.format(im)
        print (CasNam1_CNP + cy + cm + '.nc')
        with nc4.Dataset(CasNam1_CNP + '{}-{}.nc'.format(cy, cm), 'r') as ncf1, nc4.Dataset(CasNam2_CNP + '{}-{}.nc'.format(cy, cm), 'r' ) as ncf2:
             # g C/m2

             B_BLW1 = ncf1.variables['TOTVEGC'][0,:,:] - ncf1.variables['TOTVEGC_ABG'][0,:,:] - ncf1.variables['CPOOL'][0,:,:]
             B_BLW2 = ncf2.variables['TOTVEGC'][0,:,:] - ncf2.variables['TOTVEGC_ABG'][0,:,:] - ncf2.variables['CPOOL'][0,:,:]

             Ratio1 = ncf1.variables['TOTVEGC_ABG'][0,:,:] / (B_BLW1 + ncf1.variables['TOTVEGC_ABG'][0,:,:])
             Ratio2 = ncf2.variables['TOTVEGC_ABG'][0,:,:] / (B_BLW2 + ncf2.variables['TOTVEGC_ABG'][0,:,:])

             

             Biomass_ABG1 += (ncf1.variables['TOTVEGC_ABG'][0,:,:]+ncf1.variables['CPOOL'][0,:,:] * Ratio1)*DaysPerMon[im-1] / 365.
             Biomass_ABG2 += (ncf2.variables['TOTVEGC_ABG'][0,:,:]+ncf2.variables['CPOOL'][0,:,:] * Ratio2)*DaysPerMon[im-1] / 365.

             LAI1 += ncf1.variables['TLAI'][0,:,:]*DaysPerMon[im-1] / 365.
             LAI2 += ncf2.variables['TLAI'][0,:,:]*DaysPerMon[im-1] / 365.

             if iy == 2000 and im == 1:
                lat = ncf1.variables['lat'][:]
                lon = ncf1.variables['lon'][:]
                area = ncf1.variables['area'][:]    # km
                lfrc = ncf1.variables['landfrac'][:]  # fraction


# 10-year average 
Biomass_ABG1 = Biomass_ABG1 / 10.   # g C/m2
Biomass_ABG2 = Biomass_ABG2 / 10.

print (Biomass_ABG1.max())
print (Biomass_ABG1.min())

LAI1 = LAI1 / 10.
LAI2 = LAI2 / 10.

GBLSUM_Biomass_ABG1 = ma.sum(Biomass_ABG1 * area * 1e6 * lfrc * 1.e-15)   # Pg C
GBLSUM_Biomass_ABG2 = ma.sum(Biomass_ABG2 * area * 1e6 * lfrc * 1.e-15)

Beta_Biomass_ABG1 = (GBLSUM_Biomass_ABG2 - GBLSUM_Biomass_ABG1) / GBLSUM_Biomass_ABG1 * 100.


print ('CNTL total biomass (PgC)', GBLSUM_Biomass_ABG1)
print ('CO2E total biomass (PgC)', GBLSUM_Biomass_ABG2)
print ('Beta=', Beta_Biomass_ABG1/2.)

lons, lats = np.meshgrid(lon, lat)

# masked when < 2000 g C/m2
AdifBiomass_ABG = ma.masked_where(Biomass_ABG1 < 2000, Biomass_ABG2 - Biomass_ABG1)
RdifBiomass_ABG = ma.masked_where(Biomass_ABG1 < 2000, AdifBiomass_ABG / Biomass_ABG1 * 100.)

AdifLAI         = ma.masked_where(Biomass_ABG1 < 2000, LAI2 - LAI1)
RdifLAI         = ma.masked_where(Biomass_ABG1 < 2000, AdifLAI / LAI1 * 100.)



Lat_AdifLAI     = ma.sum(AdifLAI * area [:,:] * lfrc [:,:], axis=1) / ma.sum(area [:,:] * lfrc [:,:], axis=1)
Lat_RdifLAI     = ma.sum(RdifLAI * area [:,:] * lfrc [:,:], axis=1) / ma.sum(area [:,:] * lfrc [:,:], axis=1)
#Lat_RdifLAI     = ma.sum(AdifLAI * area [:,:] * lfrc [:,:], axis=1) * 100. / ma.sum(LAI1 * area [:,:] * lfrc [:,:], axis=1)
#AdifBiomass_ABG = Biomass_ABG2 - Biomass_ABG1

Lat_AdifBiomass_ABG     = ma.sum(AdifBiomass_ABG * area [:,:] * 1.e6 * lfrc [:,:] * 1.e-12, axis=1)   # Tg C
Lat_RdifBiomass_ABG     = ma.sum(RdifBiomass_ABG * area [:,:] * lfrc [:,:], axis=1) / ma.sum(area [:,:] * lfrc [:,:], axis=1)

#Lat_RdifBiomass_ABG     = ma.sum(AdifBiomass_ABG * 100. * area [:,:] * 1.e6 * lfrc [:,:] * 1.e-12, axis=1) * 100./  \
#                          ma.sum(Biomass_ABG1    * 100. * area [:,:] * 1.e6 * lfrc [:,:] * 1.e-12, axis=1)


fig = plt.figure(figsize=[30,10])
fig, ax = plt.subplots(1, 2, constrained_layout = True)
fig.set_figheight(9)
fig.set_figwidth(15)

x  = lat
y1 = Lat_RdifBiomass_ABG
ax[0].plot(x, y1, color='k', linewidth=0.3, label="LAI")
ax[0].fill_between(x, y1, 0,
                   color='C0',       # The outline color
                   alpha=0.3)         # Transparency of the fill
y2 = Lat_RdifLAI
ax02 = ax[0].twinx()  # instantiate a second axes that shares the same x-axis
ax02.plot(x, y2,  color='k', linewidth=0.3, label='4xCO2')
ax02.fill_between(x, y2, 0,
                  color='lightgreen',       # The outline color
                  alpha=0.3)         # Transparency of the fill
ax02.set_ylim(0,16)
ax[0].set_ylim(0, 16)

ax[0].set_xlim(-60, 90)
ax02.set_xlim(-60, 90)
ax[0].tick_params(labelsize=18)
ax02.tick_params(labelsize=18)

ax[0].set_xlabel('Latitude', fontsize=20)
ax[0].set_ylabel('Relative change in biomass (%)', fontsize=20)
ax02.set_ylabel('Relative change in LAI (%)', fontsize=20, rotation=270, labelpad=20)

x  = lat
y1 = Lat_AdifBiomass_ABG
ax[1].plot(x, y1, color='k', linewidth=0.3, label="Biomass")
p1=ax[1].fill_between(x, y1, 0,
                   color='C0',       # The outline color
                   alpha=0.3)         # Transparency of the fill

ax12 = ax[1].twinx()  # instantiate a second axes that shares the same x-axis

y2 = Lat_AdifLAI
ax12.plot(x, y2*100, color='k', linewidth=0.3, label='LAI')
ax[1].tick_params(labelsize=18)
p2=ax12.fill_between(x, y2*100, 0,
                  color='lightgreen',       # The outline color
                  alpha=0.3)         # Transparency of the fill
ax12.tick_params(labelsize=18)

ax[1].set_ylim(0, 900)
ax12.set_ylim(0, 90)
ax[1].set_xlim(-60, 90)
ax12.set_xlim(-60, 90)



ax[1].legend(handles=[p1, p2], labels=['Biomass', 'LAI'], prop={"size":16})

ax[1].set_xlabel('Latitude', fontsize=20)
ax[1].set_ylabel('Absolute change in biomass (Tg C)', fontsize=20)
ax12.set_ylabel('Absolute change in LAI (dm2/m2)', fontsize=20, rotation=270, labelpad=20)



#ax1.set_yticks(np.arange(-80, 100, 20))
#ax1.tick_params( labelsize=18)
#ax1.set_xticks(np.arange(-2, 3, 1))


#np.save('adif_biomass_abg.npy', AdifBiomass_ABG)
#np.save('rdif_biomass_abg.npy', RdifBiomass_ABG)
plt.savefig('fig3.png')



fig = plt.figure(figsize=[30,10])
fig, ax = plt.subplots(1, 2, constrained_layout = True, subplot_kw=dict(projection=ccrs.Robinson()))
fig.set_figheight(10)
fig.set_figwidth(20)


hex_list=["#B477A9", "#B08AA0", "#AC9C97", "#A7AC8C", "#9BB67B", "#98BE6C", "#AEC16B", 
          "#C5C368", "#D4BE5F", "#E3B956", "#E8AD4E", "#E29A48", "#DC8742", "#D7723D",
          "#D15F37", "#CD5337", "#C94936", "#C63E36", "#C23437", "#BF2937"]

mycmap = get_continuous_cmap(hex_list)



p=ax[0].contourf(lons, lats, RdifBiomass_ABG,  levels =  np.arange(0, 31, 1), transform=ccrs.PlateCarree(), cmap=cmp.viridis_r, extend='both')
#p=ax[0].contourf(lons, lats, RdifBiomass_ABG,  levels = np.arange(0, 31, 1), transform=ccrs.PlateCarree(), cmap=mycmap, extend='both')
ax[0].coastlines()
ax[0].gridlines()
ax[0].set_global()
ax[0].set_title("(a) Relative increase in biomass (%)", fontsize=28)
ax[0].stock_img()
cbar=fig.colorbar(p, ax=ax[0], location='bottom',shrink=0.6, extend='both', extendrect=True, ticks=[0,10,20,30])
cbar.ax.tick_params(labelsize=28)
 

#1.e2 from g C/m2 to Mg C/m2
p=ax[1].contourf(lons, lats, AdifBiomass_ABG/1.e2,  levels = np.arange(0, 31, 1), transform=ccrs.PlateCarree(), cmap=cmp.viridis_r, extend='both')
#p=ax[1].contourf(lons, lats, AdifBiomass_ABG,  transform=ccrs.PlateCarree(), cmap=mycmap, extend='both')
ax[1].coastlines()
ax[1].gridlines()
ax[1].set_global()
ax[1].stock_img()
ax[1].set_title("(b) Absolute increase in biomass (MgC ha-1)", fontsize=28)
 
cbar=fig.colorbar(p, ax=ax[1], location='bottom',shrink=0.6, extend='both', extendrect=True, ticks=[0,10,20,30])
cbar.ax.tick_params(labelsize=28)


plt.savefig('fig2.png')
