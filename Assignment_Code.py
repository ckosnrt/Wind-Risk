
"""
Template code for looping over dates, reading multiple ERA5 data files in netCDF format,
sub-setting the data from a particular location and creating a time series.

Author: 2020, John Methven
"""

import matplotlib.pyplot as plt
from netCDF4 import Dataset
import cartopy.crs as ccrs
import numpy as np
import xarray as xr
import seaborn as sns
import datetime
import os
from datetime import timedelta

def read_data(filename, fyear, model):
    """
    Read in the ground variable data from the netCDF file.
    Input: name of file to read.
    Output: 
    :longitude  - degrees
    :latitude   - degrees
    :outdates   - calendar dates for data points
    :u          - zonal component of the wind (m/s)
    :v          - meridional component of the wind (m/s)
    :mslp       - mean sea level pressure (hPa)

    """
    
    # Read netCDF file into an xarray DataSet - a form of catalogue
    ds = xr.open_dataset(filename)

    ftime = ds['time']

    if model.lower() == 'era5':
        alon = ds['lon'][:]
        alat = ds['lat'][:]
        u = ds['u10'][:,:,:]
        v = ds['v10'][:,:,:]
        if any((fyear < 2000, fyear >= 2018)):
            mslp = 0.01*ds['msl']
        else:
            mslp = 0.01*ds['sp']
    else:
        alon = ds['longitude'][:]
        alat = ds['latitude'][:]
        u = ds['u10'][:,:,:].mean(dim='number')
        v = ds['v10'][:,:,:].mean(dim='number')
        mslp = 0.01*ds['msl'][:,:,:].mean(dim='number')

    ds.close()
    return alon, alat, ftime, u, v, mslp, ds


def extract_series(fpath, fstem, lonpick, latpick, dstart, dend):
    '''
    High level function controlling extraction of runoff time series 
    for chosen location.
    Input: fpath, fstem determine the name of file to read
    :lonpick    - longitude of chosen location
    :latpick    - latitude of chosen location
    :dstart     - start date in datetime.date format
    :dend       - end date in datetime.date format
    Output: 
    :dayarr     - time in days since start
    :timarr     - time series in datetime format
    :windarr    - wind speed (m/s) time series at chosen location
    '''   
    #
    # Set end date and start date of required time series
    #
    dendp = dend+timedelta(days=1)
    tinterval = dendp-dstart
    ndays = tinterval.days
    #
    # Plot the data for the first date in the interval
    #
    fdate = dstart.strftime("%Y_%m")
    fyear = dstart.year
    
    #Read the data
    filename = str(fpath+fstem+fdate+'_DET.nc')
    # Note that the str() function is included to ensure that these
    # variables are interpreted as character strings.
    
    alon, alat, ftime, u, v, mslp, ds = read_data(filename, fyear, model='ERA5')
    
    # Find the indices of the grid box centred closest to the chosen location
    intlon, intlat = subset_field2(alon, alat, lonpick, latpick)
    
    # Setup arrays to save time series data
    dayarr = np.arange(ndays)
    timarr = np.arange(np.datetime64(str(dstart)), np.datetime64(str(dendp)))
    windarr = np.zeros(ndays)
    #
    # Loop over months, reading files and saving daily data
    #
    icarryon = 1
    dcur = dstart
    n=0
    while icarryon == 1:
        fdate = dcur.strftime("%Y_%m")
        fyear = dcur.year
        #Read the data
        filename = str(fpath+fstem+fdate+'_DET.nc')
        # Note that the str() function is included to ensure that these
        # variables are interpreted as character strings.
        alon, alat, outdates, u, v, mslp, ds = read_data(filename, fyear, model='ERA5')
        npts = len(outdates)
        #
        # Check whether the last day requested is within this month
        #
        if n+npts >= ndays:
            npts = ndays-n
            icarryon = 0
        #
        # Save the daily data required from this file
        #
        for i in range(npts):
            windarr[n+i] = np.sqrt(u[i, intlat, intlon]**2 + v[i, intlat, intlon]**2)
        # Increment the date variable by number of days in this month
        n = n+npts
        dcur=dcur+timedelta(days=npts)
    
    return dayarr, timarr, windarr


def read_S2S(fstem, fname):
    '''
    Read in the S2S variable data from the netCDF file.
    Input: name of file to read.
    Output: 
    :longitude  - degrees
    :latitude   - degrees
    :newcal     - calendar dates for data points
    :number     - ensemble member
    :u          - zonal component of the wind (m/s)
    :v          - meridional component of the wind (m/s)
    :mslp       - mean sea level pressure (hPa)
    '''
    filename = str(fstem+fname)
    data = Dataset(filename, 'r')
    
    rtime = data.variables['time'][:]
    alon = data.variables['longitude'][:]
    alat = data.variables['latitude'][:]
    member = data.variables ['number'][:]
    u = data.variables['u10'][:,:,:,:]
    v = data.variables['v10'][:,:,:,:]
        
    data.close()
   
    ftime = float(rtime[0])
    dtime = timedelta(hours=ftime)
    
    startcal = datetime.datetime(1900, 1, 1)
    newcal = startcal+dtime
    print(newcal)
    leadtime = np.arange(float(len(rtime)))

    return  newcal, leadtime, member, alon, alat, u, v


def extract_S2S (fdir, lonpick, latpick, dir_type):
    '''
    extract for model data.
    Input: fdir determine the name of file to read
    :lonpick    - longitude of chosen location
    :latpick    - latitude of chosen location
    dir_type    - type of directory
    '''
    # Make a Python list containing the names of files in the chosen directory
    #    
    file_list = os.listdir(fdir)
    n_forc = len(file_list)
    
    # Read the first file in the list and find the data dimensions
    filename = file_list[0]
    newcal, leadtime, member, alon, alat, u, v = read_S2S(fdir, filename)
    nlead = len(leadtime)
    nmem = len(member)
    # nlon = len(alon)
    # nlat = len(alat)
    

    intlon, intlat = subset_field2(alon, alat, lonpick, latpick)
    windarr_fc = np.zeros((nlead, nmem, n_forc))
   
    # Loop over all the forecast files in the directory
    icount = -1
    for filename in file_list:
        icount = icount+1
        newcal, leadtime, member, alon, alat, u, v = read_S2S(fdir, filename)
        for i in range(0, nlead):
            for j in range(0, nmem):
                windarr_fc[i,j,icount] = np.sqrt ((u [i,j, intlat, intlon])**2 + (v[i,j, intlat, intlon]**2))
    if dir_type.lower() == 'ecmwf':
        # loop transfer the forecast timelead into it its verifiction date 
        windarr_fc_tlead =np.zeros(((nlead, nmem, (n_forc-47))))
        for i in range(nlead):
            windarr_fc_tlead[i] = windarr_fc[i,:,47-i:n_forc-i]
    else:
        windarr_fc_tlead =np.zeros(((nlead, nmem, (n_forc-45))))
        for i in range(nlead):
            windarr_fc_tlead[i] = windarr_fc[i,:,45-i:n_forc-i]
    return  leadtime, windarr_fc, windarr_fc_tlead


def subset_field(alon, alat, lonpick, latpick):
    """
    Find the indices of the grid point centred closest to chosen location.
    Input: 
    :alon       - longitude points
    :alat       - latitude points
    :lonpick    - longitude of chosen location
    :latpick    - latitude of chosen location
    Output:
    :intlon     - index of longitude for chosen point
    :intlat     = index of latitude for chosen point
    """
    
    dlon = alon[1]-alon[0]
    dlat = alat[1]-alat[0]
    lonwest = alon[0]-0.5*dlon
    latnorth = alat[0]-0.5*dlat
    intlon = int(round((lonpick-lonwest)/dlon))
    intlat = int(round((latpick-latnorth)/dlat))
    return intlon, intlat


def subset_field2(alon, alat, lonpick, latpick):
    """
    Find the indices of the grid point centred closest to chosen location.
    Input: 
    :alon       - longitude points
    :alat       - latitude points
    :lonpick    - longitude of chosen location
    :latpick    - latitude of chosen location
    Output:
    :intlon     - index of longitude for chosen point
    :intlat     = index of latitude for chosen point
    """
    
    dlon = (alon[1] - alon[0]).item()
    dlat = (alat[1] - alat[0]).item()
    lonwest = (alon[0] - 0.5 * dlon).item()
    latnorth = (alat[0] - 0.5 * dlat).item()

    intlon = int(round((lonpick - lonwest) / dlon))
    intlat = int(round((latpick - latnorth) / dlat))
    return intlon, intlat


def calc_distribution(fpath, fstem, lonpick, latpick, dstart, dend, model):
    '''
    Calculates and optionally plots the wind speed distribution 
    for a specified location and time period.

    Parameters
    ----------
    fpath, fstem : str, determine the name of file to read.
    lonpick : the longitude of the location.
    latpick : the latitude of the location.
    dstart : start date in datetime.date format
    dend : end date in datetime.date format
    model : the name of the model.
   
    Returns
    -------
    ws_data : wind speed data collected over the specified period and location.

    '''
    ws_data = []
    current_date = dstart
    
    while current_date <= dend:
        if model.lower() == 'era5':
            fdate = current_date.strftime("%Y_%m")
            filename = f"{fpath}{fstem}{fdate}_DET.nc"
            next_date_func = lambda date: (date.replace(day=28) + timedelta(days=4)).replace(day=1)
        else:
            fdate = current_date.strftime("%Y%m%d")
            filename = f"{fpath}{fstem}{fdate}.nc"
            next_date_func = lambda date: date + timedelta(days=1)
            
        if os.path.exists(filename):
            # Read data once for the current date
            alon, alat, outdates, u, v, mslp, ds = read_data(filename, current_date.year, model)
            
            # Calculate the subset field for the location
            intlon, intlat = subset_field(alon.values, alat.values, lonpick, latpick)
            
            # Calculate wind speed
            wind_speed = np.sqrt(u[:, intlat, intlon]**2 + v[:, intlat, intlon]**2)
            ws_data.extend(wind_speed.values)
            
        current_date = next_date_func(current_date)
    return ws_data
    

def set_histograms(ws_data, locations, model):
    '''
    Plot the wind speed histograms for specified locations and time period.

    Parameters
    ----------
    ws_data : dict of wind speed data for each location
    locations : dict with the name of the location and its coordinates
    model : the name of the model

    Returns
    -------
    None.
    '''
    n_locations = len(ws_data)
    # Create a figure with subplots for histograms
    fig, axs = plt.subplots(1, n_locations, figsize=(20, 6), sharey=True)
    fig.suptitle(f'Wind Speed Analysis Distribution using {model}', fontsize=20, y=1.05)
     
    for i, (location, data) in enumerate(ws_data.items()):
        percentiles = np.percentile(data, [50, 75, 98])
        sns.histplot(data, bins=20, kde=True, color='green', stat='density', alpha=0.5, label=f'{location} Wind Speed Distribution', ax=axs[i])
        axs[i].axvline(percentiles[0], color='red', linestyle='--', alpha=0.7, label=f'50th Percentile: {percentiles[0]:.2f} m/s')
        axs[i].axvline(percentiles[1], color='blue', linestyle='--', alpha=0.7, label=f'75th Percentile: {percentiles[1]:.2f} m/s')
        axs[i].axvline(percentiles[2], color='purple', linestyle='--', alpha=0.7, label=f'98th Percentile: {percentiles[2]:.2f} m/s')
        axs[i].set_xlabel('Wind Speed (m/s)')
        axs[i].set_ylabel('Density')
        axs[i].set_title(f'{location}')
        axs[i].legend(loc='upper right')

    plt.show()
  
    
def set_boxplot(ws_data, locations, model):
    '''
    Plot a combined boxplot for the wind speed data for specified locations and time period.

    Parameters
    ----------
    ws_data : dict of wind speed data for each location
    locations : dict with the name of the location and its coordinates
    model : the name of the model

    Returns
    -------
    None.
    '''
    plt.figure(figsize=(5, 5))
    boxplot_data = [data for location, data in ws_data.items()]
    boxplot_labels = [location for location in ws_data.keys()]
    boxplot_elements = plt.boxplot(boxplot_data, vert=True, labels=boxplot_labels)
    legend_labels = []
    
    for i, (location, data) in enumerate(ws_data.items()):
        percentile_98 = np.percentile(data, 98)
        # Draw a line at the 98th percentile for each box
        plt.hlines(percentile_98, i + 1 - 0.25, i + 1 + 0.25, color='red', linestyles='dashed')
        legend_labels.append(f"{location}: {percentile_98:.2f} m/s")
        
    plt.legend(legend_labels, loc='upper right')

    plt.title(f'Comparative Wind Speed Boxplot Analysis using {model}')
    plt.ylabel('Wind Speed (m/s)')
    plt.grid(True)
    plt.tight_layout()
    plt.show()



def plot_distribution(fpath, fstem, dstart, dend, locations, model):
    '''
    Plots the wind speed distribution for specified locations and time period.

    Parameters
    ----------
    fpath, fstem : str, determine the name of file to read.
    dstart, dend : start and end date in datetime.date format.
    locations : dict with the name of the location and its coordinates.
    model : the name of the model.

    Returns
    -------
    None.
    
    '''
    ws_data = {}
    for location, coords in locations.items():
        lonpick, latpick = coords['lon'], coords['lat']
        data = calc_distribution(fpath, fstem, lonpick, latpick, dstart, dend, model)
        ws_data[location] = data
        
    set_histograms(ws_data, locations, model)
    set_boxplot(ws_data, locations, model)
    

def plot_leadtime(leadtime, percentile_data, yobs_data, locations, colors):
    '''
    Plot the 98th percentile of S2S data from different leadtime and compare to 
    the obs/rean
    Inputs:
        leadtime        - number of leadtime,
        percentile_data - 98th percentile from different leadtime 
        yobs_data       - 98th percentile value from Rean/obs
        locations       - locations
        colors          -color for each location
    '''
    plt.figure(figsize=(10, 6))
    for name in locations:
        plt.plot(leadtime, percentile_data[name], label=f"{name} forecast", color=colors[name])
        plt.axhline(y=yobs_data[name], linestyle='--', label=f"{name} reanalysis", color=colors[name])
    
    plt.xlabel("Lead time")
    plt.ylabel("windspeed (m/s)")
    plt.title('98th Percentile of Wind Speed')
    plt.legend()
    plt.show()
    return


def leadtime_process(fpath, fstem, dstart, dend, locations, fdir, dir_type):
    """
    Process and plot the 98th percentile wind speed forecast and reanalysis data 
    for multiple locations over specified lead times.

    Parameters
    ----------
    .
    fpath, fstem, fdir : str, determine the name of file to read.
    dstart : start date in datetime.date format.
    dend : end date in datetime.date format.
    locations : the name of the location
    
    Returns
    -------
    None.

    """
    percentile_data = {}
    yobs_data = {}
    for name, coords in locations.items():
        dayarr, timarr, windarr = extract_series(fpath, fstem, coords['lon'], coords['lat'], dstart_s2s, dend_s2s)
        p98_reanalysis = np.percentile(windarr, 98)
        leadtime, windarr_fc, windarr_fc_tlead = extract_S2S(fdir, coords['lon'], coords['lat'], dir_type)
        windarr_fc_tlead_avg = np.mean(windarr_fc_tlead, axis=1)
        p98_fc = np.percentile(windarr_fc_tlead_avg, 98, axis=1)
        percentile_data[name] = p98_fc
        yobs_data[name] = p98_reanalysis  
    plot_leadtime(leadtime, percentile_data, yobs_data, locations.keys(), colors)
    return


def calc_ws_percentile(fpath, fstem, dstart, dend, model, percentile):
    '''
    Calculates plots the specified percentile of wind speed 
    over a defined region and time period.

    Parameters
    ----------
    fpath, fstem : str, determine the name of file to read.
    dstart : start date in datetime.date format.
    dend : end date in datetime.date format.
    model : the name of the model.
    percentile : the percentile to calculate.

    Returns
    -------
    v_percentile : the spatial distribution of the specified percentile 
                    of wind speed.
    all_u: distribution of u speed
    all_v: distribution of v speed
    all_mslp: distribution of mslp
    all_ws: distribution of wind magnitude
    alon, alat: longitude, latitude

    '''
    ws_data = []
    u_data = []
    v_data = []
    mslp_data =[]
    current_date = dstart
    
    while current_date <= dend:
        if model.lower() == 'era5':
            fdate = current_date.strftime("%Y_%m")
            filename = f"{fpath}{fstem}{fdate}_DET.nc"
            next_date_func = lambda date: (date.replace(day=28) + timedelta(days=4)).replace(day=1)
        else:
            fdate = current_date.strftime("%Y%m%d")
            filename = f"{fpath}{fstem}{fdate}.nc"
            next_date_func = lambda date: date + timedelta(days=1)
        
        if os.path.exists(filename):
            # Read data once for the current date
            alon, alat, outdates, u, v, mslp, ds = read_data(filename, current_date.year, model)
            wind_speed = np.sqrt(u**2 + v**2)
            ws_data.append(wind_speed)
            u_data.append(u)
            v_data.append(v)
            mslp_data.append(mslp)
            
        current_date = next_date_func(current_date)
        if not ws_data:
            raise ValueError('No wind speed data')
    
    all_ws = xr.concat(ws_data, dim='time')
    all_u = xr.concat(u_data, dim='time').mean(dim='time')
    all_v = xr.concat(v_data, dim='time').mean(dim='time')
    all_mslp = xr.concat(mslp_data, dim='time').mean(dim='time')  
    v_percentile = all_ws.quantile(percentile/100, dim='time')
    return v_percentile, all_u, all_v, all_mslp, all_ws, alon, alat


def plot_ws_percentile(fpath, fstem, dstart, dend, regions, locations, model, percentile):
    '''
    plot the specified percentile of wind speed 
    over a defined region and time period.

    Parameters
    ----------
    fpath, fstem : str, determine the name of file to read.
    dstart : start date in datetime.date format.
    dend : end date in datetime.date format.
    regions : the bounding box of the region.
    locations : the name of the location.
    model : the name of the model.
    percentile : the percentile to calculate.
    plot : whether to plot the wind speed at the specified percentile.

    Returns
    -------
    None.

    '''
    v_percentile, all_u, all_v, all_mslp, all_ws, alon, alat = calc_ws_percentile(fpath, fstem, dstart, dend, model, percentile)
    plt.figure(figsize=(10, 4), layout='constrained')   
    ax = plt.axes(projection=ccrs.PlateCarree())
    ax.coastlines(resolution='50m')
    ax.gridlines(draw_labels=True, color='gray', alpha=0.5)
    # ax.set_extent([regions['lon_min'], regions['lon_max'], regions['lat_min'], 
    #                 regions['lat_max']], crs=ccrs.PlateCarree())
    
    # Plotting MSLP contours and streamline
    if model.lower() == 'era5':
        contour = ax.contour(all_mslp.lon, all_mslp.lat, all_mslp, levels=7, 
                             colors='black', transform=ccrs.PlateCarree())
        v_contour = ax.contourf(v_percentile.lon, v_percentile.lat, v_percentile, 
                                levels=10, cmap='Greens', transform=ccrs.PlateCarree(), extend='both')
    else:
        contour = ax.contour(all_mslp.longitude, all_mslp.latitude, all_mslp, 
                            levels=7, colors='black', transform=ccrs.PlateCarree())   
        v_contour = ax.contourf(v_percentile.longitude, v_percentile.latitude, v_percentile, 
                            levels=10, cmap='Greens', transform=ccrs.PlateCarree(), extend='both')
    
    plt.colorbar(v_contour, ax=ax, label='Wind Speed (m/s)')   
    plt.clabel(contour, inline=True, fontsize=8, fmt='%1.0f')
    ax.streamplot(alon, alat, all_u, all_v, transform=ccrs.PlateCarree(), color='blue', density=1.2, linewidth=1)
    ax.set_title(f'{percentile}th Percentile of Wind Speed ({model})')
    
    for location, coords in locations.items():
        ax.plot(coords['lon'], coords['lat'], 'or', label=location)
        ax.text(coords['lon'], coords['lat'], location, transform=ccrs.PlateCarree(), fontsize=15)

    plt.show()
    return


def calc_alp(v_percentile, all_ws):
    '''
    Calculates the Accumulated Loss Potential (ALP) based on wind speed 
    exceeding a specified percentile threshold over a defined region and 
    time period.

    Parameters
    ----------
    v_percentile : the spatial distribution of the specified percentile 
                    of wind speed
    all_ws: distribution of wind magnitude
    
    Returns
    -------
    ssi_daily: the spatial distribution of the ALP
    alp : the spatial distribution of the specified percentile 
          of potential loss.

    '''
   
    ssi_daily = xr.where(all_ws > v_percentile, ((all_ws / v_percentile) - 1) ** 3, 0)
    alp = ssi_daily.sum(dim="time")
    return ssi_daily, alp


def plot_alp(v_percentile, all_ws, range_cb, model):
    """
    Plots the Accumulated Loss Potential (ALP) based on wind speed 
    exceeding a specified percentile threshold over a defined region and 
    time period.

    Parameters
    ----------
    v_percentile : the spatial distribution of the specified percentile 
                    of wind speed
    all_ws: distribution of wind magnitude
    range_cb: range for colorbar
    model : the name of the model.
    
    Returns
    -------
    None.

    """
    ssi_daily, alp = calc_alp(v_percentile, all_ws)
                     
    plt.figure(figsize=(10, 4), layout='constrained')   
    ax = plt.axes(projection=ccrs.PlateCarree())
    if model.lower() == 'era5':
        contour = ax.contourf(alp.lon, alp.lat, alp, 
                              transform=ccrs.PlateCarree(), cmap='Reds', extend='both')
    else:
        contour = ax.contourf(alp.longitude, alp.latitude, alp, range_cb,
                              transform=ccrs.PlateCarree(), cmap='Reds', extend='both')
    
    plt.colorbar(contour, ax=ax, label='Storm Severity Index')
    ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True, color='gray', alpha=0.5) 
    # ax.set_extent([regions['lon_min'], regions['lon_max'], regions['lat_min'], 
    #                regions['lat_max']], crs=ccrs.PlateCarree())
    ax.coastlines()
    ax.set_title(f'Accumulated Loss Potential ({model})')
    
    for location, coords in locations.items():
        ax.plot(coords['lon'], coords['lat'], 'or', label=location)
        ax.text(coords['lon'], coords['lat'], location, 
                transform=ccrs.PlateCarree(), fontsize=15)

    plt.show()
    return


def fc_storm(v_percentile, regions, locations, time_storm, model):
    """
    Forecasts storm intensity based on wind speed exceeding 
    a specified percentile threshold  using a specific atmospheric model.

    Parameters
    ----------
    v_percentile : the spatial distribution of the specified percentile 
                    of wind speed
    regions : the bounding box of the region.
    locations : the name of the location.
    time_storm : time of storm
    model : the name of the model.

    Returns
    -------
    None.

    """
    if model.lower() == 'era5':
        fpath = f'../Assignment2_Climate Services/europe_1halfx1half_ERA5_winds/eur_remap_bilinear_1halfx1half_ERA5_3hr_'
        list_name = ['1999_12_DET.nc']
    else:
        fpath = f'../Assignment2_Climate Services/europe_{model}_inst_energy/eur_inst_energy_hc_'
        list_name = ['19991201.nc', '19991208.nc', '19991212.nc',
                     '19991215.nc', '19991219.nc', '19991222.nc']

    for data_name in list_name:
        ds = xr.open_dataset(fpath+data_name).sel(time=time_storm) # , method='nearest')
        ftime = ds['time']
        if model.lower() == 'era5':
            alon = ds['lon']
            alat = ds['lat']
            u = ds['u10']
            v = ds['v10']
        else:
            alon = ds['longitude']
            alat = ds['latitude']
            u = ds['u10'].mean(dim='number')
            v = ds['v10'].mean(dim='number')
        
        
        ws = np.sqrt(u**2 + v**2)
        ssi = xr.where(ws > v_percentile, ((ws / v_percentile) - 1) ** 3, 0)
        plt.figure(figsize=(10, 4), layout='constrained')   
        ax = plt.axes(projection=ccrs.PlateCarree())
        
        range_cb = np.arange(0, 1, 0.1)
        contour = ax.contourf(alon, alat, ssi, range_cb,
                              transform=ccrs.PlateCarree(), cmap='Reds', extend='both')
        
        ax.streamplot(alon, alat, u, v, transform=ccrs.PlateCarree(), color='blue', density=1.2, linewidth=1)
        plt.colorbar(contour, ax=ax, label='Storm Severity Index')
        ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True, color='gray', alpha=0.5) 
        # ax.set_extent([regions['lon_min'], regions['lon_max'], regions['lat_min'], 
        #                 regions['lat_max']], crs=ccrs.PlateCarree())
        ax.coastlines()
        ax.set_title(f'Accumulated Loss Potential ({model}) {data_name[:8]} for {time_storm}')

        for location, coords in locations.items():
            ax.plot(coords['lon'], coords['lat'], 'or', label=location)
            ax.text(coords['lon'], coords['lat'], location, 
                    transform=ccrs.PlateCarree(), fontsize=15)
    
        plt.show()
        
    return 


def extract_ssi_for_location(ssi_data, lon, lat):
    """
    Extracts the Storm Severity Index (SSI) data for a specified location 
    based on its nearest grid point.

    Parameters
    ----------
    ssi_data : A data array containing SSI values
    lon : The longitude of the location.
    lat : The latitude of the location.

    Returns
    -------
    ssi_location : The SSI values for the nearest grid point

    """
    
    nearest_lon = ssi_data.lon.sel(lon=lon, method="nearest").values
    nearest_lat = ssi_data.lat.sel(lat=lat, method="nearest").values
    ssi_location = ssi_data.sel(lon=nearest_lon, lat=nearest_lat)
    return ssi_location

def calculate_probability_above_ssi(ssi_daily, ssi_storm_value):
    """
    Calculates the probability that the SSI values on a daily basis 
    exceed the SSI value for a specific storm event.

    Parameters
    ----------
    ssi_daily : Daily SSI values over a time period
    ssi_storm_value : The SSI value for the storm event
    
    Returns
    -------
    probability : The probability (in percentage) that daily SSI values 
    exceed the SSI value of the storm event

    """
    count_above = (ssi_daily > ssi_storm_value).sum()
    total_days = ssi_daily.size
    probability = (count_above / total_days) * 100
    return probability

def plot_ssi_time_series(ssi_daily_era, ssi_storm, locations, time_storm):
    """"
    Plots the time series of daily SSI values for specified locations and 
    highlights the SSI value during a specific storm event.

    Parameters
    ----------
    ssi_daily_era : Daily SSI values over a time period for model data
    ssi_storm : The SSI value for the storm event
    locations : the name of the location.
    time_storm : time of storm

    Returns
    -------
    None.

    """
    for location, coords in locations.items():
        ssi_daily_location = extract_ssi_for_location(ssi_daily_era, coords['lon'], coords['lat'])
        ssi_storm_location = extract_ssi_for_location(ssi_storm, coords['lon'], coords['lat'])

        probability = calculate_probability_above_ssi(ssi_daily_location, ssi_storm_location)
        print(f'Probability SSI Storm > SSI Daily on {time_storm} in {location}: {probability:.10f}%')
        plt.figure(figsize=(10, 4))
        ssi_daily_location.plot(label=f'SSI Daily for {location}', color='skyblue')
        plt.axhline(y=ssi_storm_location, color='r', linestyle='--', label=f'SSI on {time_storm}')
        plt.title(f'SSI Time Series for {location}')
        plt.xlabel('Time')
        plt.ylabel('SSI')
        plt.legend()
        plt.show()  


#%%
'''Main'''
# Select the start and end date required for the time series
dstart_era = datetime.date(1979, 1, 1)
dend_era = datetime.date(2018, 12, 31)

dstart_s2s = datetime.date(1999, 12, 1)
dend_s2s = datetime.date(2010, 12, 28)

dstart_storm = datetime.date(1999, 12, 22)
dend_storm = datetime.date(1999, 12, 29)

# Set the path and filename stem for data files.   
fpath = '../Assignment2_Climate Services/'
fstem_era = 'europe_1halfx1half_ERA5_winds/eur_remap_bilinear_1halfx1half_ERA5_3hr_'
fstem_ecmwf = 'europe_ECMWF_inst_energy/eur_inst_energy_hc_'   
fstem_ncep = 'europe_NCEP_inst_energy/eur_inst_energy_hc_' 
fdir_ecmwf = '../Assignment2_Climate Services/europe_ECMWF_inst_energy/'
fdir_ncep = '../Assignment2_Climate Services/europe_NCEP_inst_energy/'

regions = {'lat_min':  47, 
            'lat_max': 60, 
            'lon_min': -8, 
            'lon_max': 15}

locations = {
    'London': {'lon':-0.1278, 'lat':51.5074},
    'Paris': {'lon':2.3522, 'lat':48.8566},
    'Hamburg': {'lon':9.9937, 'lat':53.5511}}

colors = {
    'London': 'blue',
    'Paris': 'green',
    'Hamburg': 'red'}

v_percentile_era, all_u_era, all_v_era, all_mslp_era, all_ws_era, alon, alat = calc_ws_percentile(fpath, fstem_era, dstart_era, dend_era, model='ERA5', percentile=98)
v_percentile_nce, all_u_nce, all_v_nce, all_mslp_nce, all_ws_nce, alon, alat = calc_ws_percentile(fpath, fstem_ncep, dstart_s2s, dend_s2s, model='NCEP', percentile=98)
v_percentile_ecm, all_u_ecm, all_v_ecm, all_mslp_ecm, all_ws_ecm, alon, alat = calc_ws_percentile(fpath, fstem_ecmwf, dstart_s2s, dend_s2s, model='NCEP', percentile=98)

v_percentile_sera, all_u_era, all_v_era, all_mslp_era, all_ws_sera, alon, alat = calc_ws_percentile(fpath, fstem_era, dstart_storm, dend_storm, model='ERA5', percentile=98)
v_percentile_snce, all_u_nce, all_v_nce, all_mslp_nce, all_ws_snce, alon, alat = calc_ws_percentile(fpath, fstem_ncep, dstart_storm, dend_storm, model='NCEP', percentile=98)
v_percentile_secm, all_u_ecm, all_v_ecm, all_mslp_ecm, all_ws_secm, alon, alat = calc_ws_percentile(fpath, fstem_ecmwf, dstart_storm, dend_storm, model='NCEP', percentile=98)

time_storm = '1999-12-25'
ssi_daily_era, alp_era = calc_alp(v_percentile_era, all_ws_era)
ssi_daily_sera, alp_sera = calc_alp(v_percentile_sera, all_ws_sera)
ssi_storm = ssi_daily_sera.sel(time=time_storm)


'''Question 1'''
plot_distribution(fpath, fstem_era, dstart_era, dend_era, locations, model='ERA5')

'''Question 2'''
plot_ws_percentile(fpath, fstem_era, dstart_era, dend_era, regions, locations, model='ERA5', percentile=98)

'''Question 3'''
range_cb = np.arange(0, 7, 1)
plot_alp(v_percentile_era, all_ws_era, range_cb, model='ERA5')

'''Question 4'''
plot_distribution(fpath, fstem_ecmwf, dstart_s2s, dend_s2s, locations, model='ECMWF')
plot_distribution(fpath, fstem_ncep, dstart_s2s, dend_s2s, locations, model='NCEP')
leadtime_process(fpath, fstem_era,dstart_s2s, dend_s2s, locations, fdir_ecmwf, dir_type='ecmwf')
leadtime_process(fpath, fstem_era,dstart_s2s, dend_s2s, locations, fdir_ncep, dir_type='ncep')

'''Question 5'''
range_cb = np.arange(0, 10, 1)
plot_alp(v_percentile_sera, all_ws_sera, range_cb, model='ERA5')
plot_alp(v_percentile_secm, all_ws_secm, range_cb, model='ECMWF')
plot_alp(v_percentile_snce, all_ws_snce, range_cb, model='NCEP')

'''Question 6(i)'''
fc_storm(v_percentile_era, regions, locations, time_storm, model='ERA5')
fc_storm(v_percentile_nce, regions, locations, time_storm, model='NCEP')
fc_storm(v_percentile_ecm, regions, locations, time_storm, model='ECMWF')

'''Question 6 (ii)'''
plot_ssi_time_series(ssi_daily_era, ssi_storm, locations, time_storm) 