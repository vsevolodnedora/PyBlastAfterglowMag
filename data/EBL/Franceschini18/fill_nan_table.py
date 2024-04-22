import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.colors as pltcol
from scipy import interpolate
import h5py
from scipy.interpolate import griddata
from scipy.spatial import cKDTree

def find_nearest_index(array, value):
    ''' Finds index of the value in the array that is the closest to the provided one '''
    idx = (np.abs(array - value)).argmin()
    return idx


def fill_nan_by_extrapolation(z, x_grid, y_grid):
    # Create meshgrid for interpolation
    X, Y = np.meshgrid(x_grid, y_grid)

    # Flatten the grids and data for processing with griddata
    points = np.column_stack((X.ravel(), Y.ravel()))
    values = z.ravel()

    # Mask for valid (non-NaN) and invalid (NaN) points
    valid_mask = ~np.isnan(values)
    missing_mask = np.isnan(values)

    # Only keep valid points and values
    valid_points = points[valid_mask]
    valid_values = values[valid_mask]

    # Initial interpolation (linear or cubic) - change method if needed
    filled_values = griddata(valid_points, valid_values, points, method='linear', fill_value=np.nan)

    # Check if there are still nans after griddata
    if np.isnan(filled_values).any():
        # Find nearest non-NaN points using a KDTree and fill NaNs
        kdtree = cKDTree(valid_points)
        distances, indices = kdtree.query(points[missing_mask])
        filled_values[missing_mask] = valid_values[indices]

    # Reshape back to the original shape of the input data array
    return filled_values.reshape(z.shape)

if __name__ == '__main__':

    table = np.loadtxt("./table.txt")

    zs = table[0, 1:]
    en = table[1:, 0]
    freqs = (en * 1e12) * (2.417989242e14) # TeV -> Hz
    new_table = table[1:,1:]

    new_table=fill_nan_by_extrapolation(new_table,zs,en)

    freq_ = 4.17852721e+25
    z_ = 0.5
    print(freqs)
    print(zs)
    val = new_table[find_nearest_index(freqs,freq_)][find_nearest_index(zs,z_)]
    print(val)

    freq_ = 4.46885518e+27
    z_ = 2.
    print(freqs)
    print(zs)
    val = new_table[find_nearest_index(freqs,freq_)][find_nearest_index(zs,z_)]
    print(val)

    # new_table[~np.isfinite(new_table)] = -1.
    new_table=np.exp(-new_table)
    fig,ax = plt.subplots(ncols=1,nrows=1,sharex='all',sharey='all',layout='constrained')
    norm = pltcol.LogNorm(vmin=1e-3,#np.min(new_table[np.isfinite(new_table) & (new_table > 0)]),
                          vmax=np.max(new_table[np.isfinite(new_table) & (new_table > 0)]))

    cntr1 = ax.pcolormesh(freqs, zs, new_table.T, cmap="jet", norm=norm)
    ax.set(xscale=('log'),yscale=('log'),ylim=(min(zs),max(zs)))
    fig.colorbar(cntr1, ax=ax)
    plt.show()

    dfile = h5py.File("./table.h5",'w')
    dfile.create_dataset(name="tau",data=new_table)
    dfile.create_dataset(name="freq",data=freqs)
    dfile.create_dataset(name="z",data=zs)
    dfile.close()

    #
    # fig,axes = plt.subplots(ncols=2,nrows=1,sharex='all',sharey='all',layout='constrained')
    # # ax.contour(en, zs, new_table.T, levels=14, linewidths=0.5, colors='k')
    # # cntr1 = ax.contourf(en, zs, new_table.T, levels=6, cmap="RdBu_r")
    # norm = pltcol.LogNorm(vmin=np.min(new_table[np.isfinite(new_table) & (new_table > 0)]),
    #                       vmax=np.max(new_table[np.isfinite(new_table) & (new_table > 0)]))
    # _ = axes[0].pcolormesh(freqs, zs, new_table.T, cmap="jet", norm=norm)
    #
    # new_table_=fill_nan_by_extrapolation(new_table.T,en,zs).T
    # new_table_=np.exp(-new_table_)
    # cntr1 = axes[1].pcolormesh(freqs, zs, new_table_.T, cmap="jet", norm=norm)
    #
    # tbl = new_table / new_table_
    #
    # fig.colorbar(cntr1, ax=axes[-1])
    # # ax.set(xlim=(0, 2), ylim=(-2, 2))
    # for ax in axes: ax.set(xscale=('log'),yscale=('log'),ylim=(min(zs),max(zs)))
    #
    #
    #
    # plt.show()