# Generate new data at: C:\Users\Dvir\Desktop\Doppler Frequencies Mapping\Doppler\Simulations\Auto-encoding denoising simulation\Reconstruction_network\Gen_data.m

import torch
import numpy as np
import mat73
import os

import DopplerDataSet


import flags as _

def tr_gen():
    if _.FLAGS['noise']:
        name = 'Piezo_samps_' + str(_.FLAGS['NUM_SIM_SAMPS']) + '_noisy.mat'
        
    else:
        name = 'Piezo_samps_' + str(_.FLAGS['NUM_SIM_SAMPS']) + '.mat'
        
    train_data_struct = mat73.loadmat(os.path.join(_.FLAGS['rootdir'], name))['data_struct']
    tr_piezo_targets = torch.tensor(np.array(train_data_struct['piezo_freqs']),
                                    dtype=torch.float32, requires_grad=False)
                                    
    # tr_piezo_targets = torch.nn.functional.normalize(tr_piezo_targets, p=2, dim=0).unsqueeze(1)
    tr_piezo_targets = tr_piezo_targets.unsqueeze(1)
    tr_measures = torch.tensor(np.array(train_data_struct['meaures']), dtype=torch.float32,
                               requires_grad=True)
    
    tr_measure_set = DopplerDataSet.DopplerDataSet(tr_measures, tr_piezo_targets)
    tr_params = {'batch_size': _.FLAGS['batch_size'], 'shuffle': _.FLAGS['shuffle']}
    tr_measure_gen = torch.utils.data.DataLoader(tr_measure_set, **tr_params)
    return tr_measure_gen


def valid_gen():
    if _.FLAGS['noise']:
        name_valid = 'Piezo_samps_' + str(_.FLAGS['NUM_SIM_VALID_SAMPS']) + '_noisy.mat'
        
    else:
        name_valid = 'Piezo_samps_' + str(_.FLAGS['NUM_SIM_VALID_SAMPS']) + '.mat'
        
    valid_data_struct = mat73.loadmat(os.path.join(_.FLAGS['valid_rootdir'], name_valid))['data_struct']
    valid_piezo_targets = torch.tensor(np.array(valid_data_struct['piezo_freqs']),
                                    dtype=torch.float32, requires_grad=False)
                                    
    # valid_piezo_targets = torch.nn.functional.normalize(valid_piezo_targets, p=2, dim=0).unsqueeze(1)
    valid_piezo_targets = valid_piezo_targets.unsqueeze(1)
    valid_measures = torch.tensor(np.array(valid_data_struct['meaures']), dtype=torch.float32,
                               requires_grad=True)
                               
    valid_measure_set = DopplerDataSet.DopplerDataSet(valid_measures, valid_piezo_targets)
    valid_params = {'batch_size': _.FLAGS['valid_batch_size'], 'shuffle': False}
    valid_measure_gen = torch.utils.data.DataLoader(valid_measure_set, **valid_params)
    return valid_measure_gen
