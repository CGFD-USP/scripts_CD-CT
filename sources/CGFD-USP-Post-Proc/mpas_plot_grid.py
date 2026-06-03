#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct  7 17:17:43 2022

@author: daniloceano
"""

import math
import os

import xarray as xr
import numpy as np

import matplotlib as mpl
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from geopy.distance import distance

from matplotlib.collections import PolyCollection
from multiprocessing import Pool

import argparse

def add_mpas_mesh_variables(ds, full=True, **kwargs):
    for v in ds.data_vars:
        if v not in derived_variables:
            continue

        newvs = derived_variables[v]

        for newv in newvs:
            if newv in ds:
                #print(newv + ' already here')
                continue

            if 'lat' in v or 'lon' in v:
                ds[newv] = xr.apply_ufunc(np.rad2deg, ds[v])
                ds[newv] = ds[newv].where(ds[newv] <= 180.0, ds[newv] - 360.0)
                ds[newv].attrs['units'] = 'degrees'

            elif newv == 'area':
                radius_circle = ds.attrs.get('sphere_radius', 1.0)
                if radius_circle == 1:
                    # need to correct to earth radius
                    correction_rad_earth = 6371220.0
                else:
                    correction_rad_earth = 1

                ds[newv] = (ds[v] / 10 ** 6) * correction_rad_earth**2
                ds[newv].attrs['units'] = 'km^2 (assuming areaCell in m^2)'
                ds[newv].attrs['long_name'] = 'Area of the cell in km^2'

            elif newv == 'resolution':
                radius_circle = ds.attrs.get('sphere_radius', 1.0)
                if radius_circle == 1.0:
                    #print('need to correct to earth radius!!')
                    correction_rad_earth = 6371220.0
                else:
                    correction_rad_earth = 1

                # km^2 (assuming areaCell in m^2)
                area = (ds[v] / 10 ** 6) * correction_rad_earth**2

                ds[newv] = 2 * (xr.apply_ufunc(np.sqrt, area / math.pi))
                ds[newv].attrs['units'] = 'km'
                ds[newv].attrs['long_name'] = 'Resolution of the cell (approx)'

    return ds


def open_mpas_file(file, **kwargs):
    ds = xr.open_dataset(file)
    ds = add_mpas_mesh_variables(ds, **kwargs)
    return ds

def start_cartopy_map_axis(zorder=1):
    ax = plt.axes(projection=ccrs.PlateCarree())  # projection type
    add_cartopy_details(ax, zorder=zorder)
    return ax

def add_cartopy_details(ax, zorder=1, tick_labelsize=16):
    ax.add_feature(cfeature.BORDERS, linestyle=':', zorder=zorder)
    ax.coastlines(resolution='10m', zorder=zorder+1)

    gl = ax.gridlines(draw_labels=True, alpha=0.5, linestyle='--',
                      zorder=zorder+2)
    gl.top_labels = False
    gl.right_labels = False

    # Set the font size for gridline labels
    gl.xlabel_style = {'size': tick_labelsize}
    gl.ylabel_style = {'size': tick_labelsize}
    
def set_plot_kwargs(da=None, list_darrays=None, **kwargs):
    plot_kwargs = {k: v for k, v in kwargs.items()
                   if k in ['cmap', 'vmin', 'vmax']
                   and v is not None}

    if 'cmap' not in plot_kwargs:
        plot_kwargs['cmap'] = 'summer' #'viridis' #'Spectral'

    vmin = plot_kwargs.get('vmin', None)
    if vmin is None:
        if da is not None:
            vmin = np.min(da)
        elif list_darrays is not None:
            vmin = np.min([v.min() for v in list_darrays if v is not None])
    if vmin is not None:
        plot_kwargs['vmin'] = vmin

    vmax = plot_kwargs.get('vmax', None)
    if vmax is None:
        if da is not None:
            vmax = np.max(da)
        elif list_darrays is not None:
            vmax = np.max([v.max() for v in list_darrays if v is not None])

    if vmax is not None:
        plot_kwargs['vmax'] = vmax

    return plot_kwargs

def colorvalue(val, cmap='summer', vmin=None, vmax=None):
    """
    Given a value and the range max, min, it returns the associated
    color of the desired cmap.
    :param val: float
    :param cmap: str
    :param vmin: float (default None)
    :param vmax: float (default None)
    :return: cm(norm_val): color
    """
    # Get a colormap instance, defaulting to rc values if name is None.
    cm = mpl.colormaps[cmap] #cm.get_cmap(cmap, None)
    if vmin is None:
        vmin = xr.DataArray.min().values  # min value of the array
    if vmax is None:
        vmax = xr.DataArray.max().values  # max value of the array
    if vmin == vmax:
        # A class which, when called, linearly normalizes data into the
        # [0.0, 1.0] interval.
        norm_val = mpl.colors.Normalize(vmin=vmin - 1, vmax=vmax + 1,
                                        clip=True)(val)
    else:
        norm_val = mpl.colors.Normalize(vmin=vmin, vmax=vmax,
                                        clip=True)(val)
    return cm(norm_val)

def plot_cells_chunk(chunk, lats, lons, vertices, num_sides, values, plot_kwargs):
    polygons = []
    colors = []
    for i in chunk:
        cell_vertices = vertices[i, :num_sides[i]]
        cell_lats = lats[cell_vertices]
        cell_lons = lons[cell_vertices]

        if np.all(cell_lons >= -179) and np.all(cell_lons <= 179):
            polygons.append(np.column_stack((cell_lons, cell_lats)))
            colors.append(colorvalue(values[i], **plot_kwargs))
    return polygons, colors

def plot_cells_mpas(ds, vname, ax, num_chunks=128, **plot_kwargs):
    # Preload data into memory
    values = ds[vname].values
    vertices = ds['verticesOnCell'].values - 1
    num_sides = ds['nEdgesOnCell'].values
    lats = ds['latitudeVertex'].values
    lons = ds['longitudeVertex'].values

    indices = np.arange(len(ds['nCells']))
    chunks = np.array_split(indices, num_chunks * 4)  # Smaller chunks for better load balancing

    with Pool(processes=num_chunks) as pool:
        results = pool.starmap(
            plot_cells_chunk,
            [(chunk, lats, lons, vertices, num_sides, values, plot_kwargs) for chunk in chunks]
        )

    polygons = []
    colors = []
    for chunk_polygons, chunk_colors in results:
        polygons.extend(chunk_polygons)
        colors.extend(chunk_colors)

    poly_collection = PolyCollection(polygons, facecolors=colors, edgecolor="black", linewidth=0.01)
    ax.add_collection(poly_collection)

def plot_cells_mpas_old(ds, vname, ax, **plot_kwargs):
    
    # ax = start_cartopy_map_axis()
             
    
    for i, cell in enumerate(ds['nCells'].values):
        value = ds[vname].sel(nCells=cell)

        vals = ds['verticesOnCell'].sel(nCells=cell).values
        num_sides = int(ds['nEdgesOnCell'].sel(nCells=cell))
        vals = vals[:num_sides] - 1
        lats = ds['latitudeVertex'].sel(nVertices=vals)
        lons = ds['longitudeVertex'].sel(nVertices=vals)
    
        
        color = colorvalue(value, **plot_kwargs)
        
        ## For some reason, when plotting from -180 to 180 
        # gives a strange mesh plot (actually, any value higher than
        # 179 gives the strange result)
        
        if all(j for j in lons >= -179) and all(j for j in lons <= 179):
            
            ax.fill(lons, lats, edgecolor=None, linewidth=0.0,
                    facecolor=color)
        
        
def plot_dual_mpas(ds, vname, ax, **plot_kwargs):
    for vertex in ds['nVertices'].values:
        value = ds[vname].sel(nVertices=vertex)

        vals = ds['cellsOnVertex'].sel(nVertices=vertex).values
        if 0 in vals:
            # Border triangle
            continue
        vals = vals - 1
        lats = ds['latitude'].sel(nCells=vals)
        lons = ds['longitude'].sel(nCells=vals)

        color = colorvalue(value, **plot_kwargs)

        ax.fill(lons, lats, edgecolor=None, linewidth=0.0,
                facecolor=color)
        
def add_colorbar(axs, fig=None, label=None, **plot_kwargs):
    if fig is None:
        fig = plt.gcf()

    try:
        x = axs[0, 0]
    except:
        try:
            x = axs[0]
            n = len(axs)
        except:
            axs = np.array([axs]).reshape([1, 1])
        else:
            axs = axs.reshape([n, 1])

    cbar = fig.colorbar(
        mpl.cm.ScalarMappable(
            norm=mpl.colors.Normalize(vmin=plot_kwargs['vmin'],
                                      vmax=plot_kwargs['vmax'], clip=True),
            cmap=plot_kwargs['cmap']),
        ax=axs[:, :], shrink=0.6)
    cbar.ax.locator_params(nbins=10)

    if plot_kwargs['vmin'] is not None and plot_kwargs['vmax'] is not None:
        # Explicitly set ticks to include the start and end of the color range
        ticks = np.linspace(plot_kwargs['vmin'], plot_kwargs['vmax'], num=7)
        cbar.set_ticks(ticks)

    if label is not None:
        cbar.set_label(label, fontsize=16)
        # Increase the font size of the colorbar ticks
        cbar.ax.tick_params(labelsize=16)

    return

def close_plot(fig=None, size_fig=None, pdf=None, outfile=None,
               force_show=False):

    if size_fig is None:
        size_fig = [10, 8]

    if fig is None:
        fig = plt.gcf()
    fig.set_size_inches(size_fig)

    if outfile is not None:
        plt.savefig(outfile, dpi=3200)

    if pdf is not None:
        pdf.savefig(fig, dpi=800)

    if (outfile is None and pdf is None) or force_show:
        plt.show()

    plt.close()


def plot_mpas_darray(ds, vname, ax=None, outfile=None, **kwargs):
    
    ## plot_mpas_darray
    if vname not in ds.data_vars:
            print('Unplottable Data Array ' + vname)
            print(ds)
            
    da = ds[vname]
    for coord in ['time', 'lev']:
        if coord in da.dims:
            print('Selecting first slice for ' + coord + '.')
            da = da.isel({coord: 0})      
            
    # final = False
    # if ax is None:
    #     final = True
    #     ax = start_cartopy_map_axis()
    
    # Apply zoom here if provided
    zoom_extent = kwargs.get('zoom_extent', [-180.0, 180, -90.0, 90.0])
    ax.set_extent(zoom_extent, crs=ccrs.PlateCarree())
    
    plot_kwargs = set_plot_kwargs(da=da, **kwargs)
    
    if 'nCells' in ds[vname].dims:
        plot_cells_mpas(ds, vname, ax, **plot_kwargs)
    # elif 'nVertices' in ds[vname].dims:
    #     plot_dual_mpas(ds, vname, ax, **plot_kwargs)
    else:
        print('WARNING  Impossible to plot!')
    
        
    units = da.attrs.get('units', '')
    name = kwargs.get('name', '')
    ncells = str(len(da.values.flatten()))
    title = kwargs.get('title', '')
    title = title.replace('<VAR>', vname).replace('<UNITS>', units)
    title = title.replace('<NAME>', name).replace('<NCELLS>', ncells)
    title = vname
    #ax.set_title(title)
    
    # if final:
    #     title_legend = kwargs.get('title_legend', '<VAR>: <UNITS>')
    #     title_legend = title_legend.replace('<VAR>', vname)
    #     title_legend = title_legend.replace('<UNITS>', units)
    #     add_colorbar(ax, label=title_legend, **plot_kwargs)
    
    #     close_plot(outfile=outfile)

def view_mpas_mesh(mpas_grid_file, outfile=None,
                            do_plot_resolution_rings=True,
                            vname='resolution',
                            border_radius=None,
                            **kwargs):
    
    ds = open_mpas_file(mpas_grid_file)
    
    units = ds[vname].attrs.get('units', '')
    ncells = str(len(ds[vname].values.flatten()))
    name = os.path.basename(mpas_grid_file)
    
    ax = start_cartopy_map_axis(zorder=2)
    plot_kwargs = set_plot_kwargs(da=ds[vname])
    
    # Add vmin and vmax from kwargs if they exist
    if 'vmin' in kwargs and kwargs['vmin'] is not None:
        plot_kwargs['vmin'] = kwargs['vmin']
    else:
        plot_kwargs['vmin'] = ds[vname].min().values
    if 'vmax' in kwargs and kwargs['vmax'] is not None:
        plot_kwargs['vmax'] = kwargs['vmax']
    else:
        plot_kwargs['vmax'] = ds[vname].max().values

     # --------
    tit = vname + ': ' + name + ' (' + str(ncells) + ')'
    array_plot_kwgs = {**plot_kwargs}
    if 'border_radius' in kwargs:
        if kwargs['border_radius'] is not None:
            array_plot_kwgs['border_radius'] = kwargs['border_radius']

    plot_mpas_darray(ds, vname, ax=ax, title=tit,
                     border_radius=border_radius, **array_plot_kwgs)
    
    add_colorbar(ax, label=vname + ' (' + units + ')', **plot_kwargs)
    
    # Apply zoom here if provided
    print ("Applying zoom extent in view_mpas_mesh:")
    zoom_extent = kwargs.get('zoom_extent', [-180.0, 180, -90.0, 90.0])
    print (zoom_extent)
    ax.set_extent(zoom_extent, crs=ccrs.PlateCarree())
    
    close_plot(outfile=outfile)
        
    
if __name__ == "__main__":
    
    derived_variables = {
        'latCell': ['latitude'],
        'lonCell': ['longitude'],
        'latVertex': ['latitudeVertex'],
        'lonVertex': ['longitudeVertex'],
        'areaCell': ['area', 'resolution'],
    }
    
    parser = argparse.ArgumentParser(
    description=__doc__,
    formatter_class=argparse.RawTextHelpFormatter)

    parser.add_argument(
        "-g", "--grid", type=str, required=True,
        help="Name of an MPAS grid.nc",
    )
    
    parser.add_argument(
        "-o", "--outfile", type=str, default=None,
        help="File to save the MPAS plot",
    )

    parser.add_argument(
        "-nc", type=int, default=1, 
        help="Number of chunks for parallel processing"
    )

    parser.add_argument(
        "-vmin", type=float, default=None,
        help="Minimum value for colorbar"
    )

    parser.add_argument(
        "-vmax", type=float, default=None,
        help="Maximum value for colorbar"
    )

    # Zoom arguments
    parser.add_argument("-lat_min", type=float, default=None, help="Minimum latitude for zoom")
    parser.add_argument("-lat_max", type=float, default=None, help="Maximum latitude for zoom")
    parser.add_argument("-lon_min", type=float, default=None, help="Minimum longitude for zoom")
    parser.add_argument("-lon_max", type=float, default=None, help="Maximum longitude for zoom")

    args = parser.parse_args()

    # Determine zoom extent
    zoom_extent = [
        args.lon_min if args.lon_min is not None else -180,
        args.lon_max if args.lon_max is not None else 180,
        args.lat_min if args.lat_min is not None else -90,
        args.lat_max if args.lat_max is not None else 90,
    ]

    if not os.path.exists(args.grid):
        raise IOError('File does not exist: ' + args.grid)
    
    view_mpas_mesh(args.grid, outfile=args.outfile, num_chunks=args.nc, vmin=args.vmin, vmax=args.vmax,
                   zoom_extent=zoom_extent)
