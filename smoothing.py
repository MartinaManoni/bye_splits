import numpy as np
import h5py
from airflow.airflow_dag import fill_path
from copy import copy

def valid1(energies, infile, outfile, nbinsRz, nbinsPhi):
    """
    compares all values of 2d histogram between local and CMSSW versions
    """
    with open(infile, 'w') as flocal, open(outfile, 'r') as fremote :
        lines = fremote.readlines()
         
        for line in lines:
            l = line.split('\t')
            if l[0]=='\n' or '#' in l[0]:
                continue
            bin1 = int(l[0])
            bin2 = int(l[1])
            val_remote = float(l[2].replace('\n', ''))
            val_local = energies[bin1,bin2]
            if abs(val_remote-val_local)>0.001:
                print('Diff found! Bin1={}\t Bin2={}\tRemote={}\tLocal={}'.format(bin1, bin2, val_remote, val_local))
        for bin1 in range(nbinsRz):
            for bin2 in range(nbinsPhi):
                flocal.write('{}\t{}\t{}\n'.format(bin1, bin2, np.around(energies[bin1,bin2], 6)))
                

def smoothAlongRz(arr, nbinsRz, nbinsPhi):
    """
    Smoothes the energy distribution of cluster energies deposited on trigger cells.
    Works along the Rz direction ("horizontal")
    """
    weights = 2 * np.ones_like(arr)
    weights[[0,nbinsRz-1],:] = 1.5

    # add top and bottom phi rows for boundary conditions
    phiPad = np.zeros((1,nbinsPhi))
    arr = np.concatenate( (phiPad,arr,phiPad) )

    arr_new = arr + ( np.roll(arr, shift=1, axis=0) + np.roll(arr, shift=-1, axis=0) ) * 0.5

    # remove top and bottom phi rows
    arr_new = arr_new[1:arr_new.shape[0]-1]
    assert(arr_new.shape[0] == nbinsRz)
    return arr_new / weights

def smoothAlongPhi(arr, binSums, nbinsRz, nbinsPhi, seedsNormByArea,
                   minROverZ, maxROverZ, areaPerTriggerCell):
    """
    Smoothes the energy distribution of cluster energies deposited on trigger cells.
    Works along the Phi direction ("horizontal")
    """
    arr_new = np.zeros_like(arr)

    nBinsSide = (np.array(binSums, dtype=np.int32) - 1) / 2; # one element per Rz bin
    assert(nBinsSide.shape[0] == nbinsRz)
    area = (1 + 2.0 * (1 - 0.5**nBinsSide)) # one element per Rz bin

    if seedsNormByArea:
        R1 = minROverZ + bin1 * (maxROverZ - minROverZ) / nbinsRz
        R2 = R1 + ((maxROverZ - minROverZ) / nbinsRz)
        area = area * ((np.pi * (R2**2 - R1**2)) / nbinsPhi);
    else:
        #compute quantities for non-normalised-by-area histoMax
        #The 0.1 factor in bin1_10pct is an attempt to keep the same rough scale for seeds.
        #The exact value is arbitrary.
        bin1_10pct = int(0.1) * nbinsRz
        R1_10pct = minROverZ + bin1_10pct * (maxROverZ - minROverZ) / nbinsRz
        R2_10pct = R1_10pct + ((maxROverZ - minROverZ) / nbinsRz)
        area_10pct_ = ((np.pi * (R2_10pct**2 - R1_10pct**2)) / nbinsPhi)
        area = area * area_10pct_;

    # loop per chunk of (equal Rz) rows with a common shift to speedup
    # unfortunately np.roll's 'shift' argument must be the same for different rows
    for idx in np.unique(nBinsSide):
        roll_indices = np.where(nBinsSide == idx)[0]
        arr_copy = arr[roll_indices,:]
        arr_smooth = arr[roll_indices,:]
        for nside in range(1, int(idx)+1):
            arr_smooth += ( (np.roll(arr_copy, shift=nside,  axis=1) +
                             np.roll(arr_copy, shift=-nside, axis=1))
                            / (2**nside) )
        arr_new[roll_indices,:] = arr_smooth / np.expand_dims(area[roll_indices], axis=-1)

    return arr_new * areaPerTriggerCell

def printHistogram(arr):
    for i in range(arr.shape[0]):
        for j in range(arr.shape[1]):
            if arr[i,j] == 0:
                print('-', end='|')
            else:
                print('X', end='|')
        print()

def createHistogram(bins, nbinsRz, nbinsPhi):
    """
    Creates a 2D histogram with fixed (R/z vs Phi) size.
    The input event must be a 2D array where the inner axis encodes, in order:
    - 1: R/z bin index
    - 2: Phi bin index
    - 3: Value of interest ("counts" of the histogram, "z axis")
    """
    arr = np.zeros((nbinsRz, nbinsPhi))

    for bin in bins[:]:
        assert(bin[0] >= 0)
        assert(bin[1] >= 0)
        rzbin = int(bin[0])
        phibin = int(bin[1])
        arr[rzbin,phibin] = bin[2]

    return arr

# Event by event smoothing
def smoothing(param, selection, **kwargs):
    insmoothing = fill_path(kwargs['SmoothingIn'], param=param, selection=selection)
    outsmoothing = fill_path(kwargs['SmoothingOut'], param=param, selection=selection) 

    with h5py.File(insmoothing,  mode='r') as storeIn, h5py.File(outsmoothing, mode='w') as storeOut :

        for falgo in kwargs['FesAlgos']:
            keys_old = [x for x in storeIn.keys()
                        if falgo in x and '_group_old' in x]
            keys_new = [x for x in storeIn.keys()
                        if falgo in x and '_group_new' in x]
            
            for kold,knew in zip(keys_old,keys_new):
                opts = (kwargs['NbinsRz'], kwargs['NbinsPhi'])
                energies_old = createHistogram(storeIn[kold][:,[0,1,2]],
                                               *opts)
                energies_new = createHistogram(storeIn[knew][:,[0,1,2]],
                                               *opts)
                wght_x_old = createHistogram(storeIn[kold][:,[0,1,3]],
                                             *opts)
                wght_y_old = createHistogram(storeIn[kold][:,[0,1,4]],
                                             *opts)
                wght_x_new = createHistogram(storeIn[knew][:,[0,1,3]],
                                             *opts)
                wght_y_new = createHistogram(storeIn[knew][:,[0,1,4]],
                                             *opts)

                # if '44317' in key:
                #     print(key)
                #     print(energies)
                #     breakpoint()

                # if '187544' in key:
                #     valid1(energies,
                #            infile='outLocalBeforeSmoothing.txt',
                #            outfile='outCMSSWBeforeSmoothing.txt')
         
                #printHistogram(ev)

                phi_opt = (kwargs['BinSums'],
                           kwargs['NbinsRz'],
                           kwargs['NbinsPhi'],
                           kwargs['SeedsNormByArea'],
                           kwargs['MinROverZ'],
                           kwargs['MaxROverZ'],
                           kwargs['AreaPerTriggerCell'])

                energies_old = smoothAlongPhi(
                    energies_old,
                    *phi_opt
                    )
                energies_new = smoothAlongPhi(
                    energies_new,
                    *phi_opt,
                    )
            
                # if '187544' in key:
                #     valid1(energies, '187544',
                #            infile='outLocalHalfSmoothing.txt',
                #            outfile='outCMSSWHalfSmoothing.txt')
         
                #printHistogram(ev)

                rz_opt = (kwargs['NbinsRz'],
                          kwargs['NbinsPhi'],)
                energies_old = smoothAlongRz(
                    energies_old,
                    *rz_opt,
                )
                energies_new = smoothAlongRz(
                    energies_new,
                    *rz_opt,
                    )
         
                #printHistogram(ev)
                # 'wght_x_new', 'wght_y_new'
                cols_old = [ 'energies_old',
                             'wght_x_old', 'wght_y_old' ] 
                cols_new = [ 'energies_new',
                             'wght_x_new', 'wght_y_new' ] 

                storeOut[kold] = (energies_old,
                                 wght_x_old, wght_y_old )
                storeOut[knew] = (energies_new,
                                 wght_x_new, wght_y_new )
                
                storeOut[kold].attrs['columns'] = cols_old
                storeOut[knew].attrs['columns'] = cols_new                
                doc_m = ( 'Energies (post-smoothing) ' +
                          'and projected bin positions' )
                doc_message = doc_m
                storeOut[kold].attrs['doc'] = doc_message
                storeOut[knew].attrs['doc'] = doc_message

if __name__ == "__main__":
    from airflow.airflow_dag import smoothing_kwargs
    smoothing( **smoothing_kwargs )
