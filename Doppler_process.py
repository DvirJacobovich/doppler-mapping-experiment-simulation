import os
import math
import time
import copy
import pathlib
import numpy as np
import torch
import torch.nn as nn
import torch.optim as optim

import Doppler_model
import Doppler_plots
import Doppler_DataLoader as data

import flags as _

from matplotlib import pyplot as plt
import scipy.io as sio

# CUDA for PyTorch
use_cuda = torch.cuda.is_available()
device = torch.device("cuda:0" if use_cuda else "cpu")
torch.backends.cudnn.benchmark = True


tr_measure_gen = data.tr_gen()
valid_measure_gen = data.valid_gen()

print("Finished loading data into torch tensors")

# Learning rate for optimizers
lr = 1e-4  # 0.0002
# Beta1 hyper-param for Adam optimizers
beta1 = 0.9

training_losses, validation_losses = [], []
breaking_points = [20, 30]

def train(model: Doppler_model, criterion=nn.MSELoss()):
    min_valid_model = copy.deepcopy(model)
    # global net_rec_batch, target_batch
    optimizer = optim.Adam(model.parameters(), lr=_.FLAGS['lr'], betas=(_.FLAGS['beta1'], _.FLAGS['beta2']))
    # Loop over epochs
    for epoch in range(1, _.FLAGS['epochs'] + 1):
        # Training
        for batch_num, (measure_batch, target_batch) in enumerate(tr_measure_gen):
            optimizer.zero_grad()

            # Transfer to GPU
            measure_batch, target_batch = measure_batch.to(_.FLAGS['device']), target_batch.to(_.FLAGS['device'])

            net_rec_batch = model(measure_batch)  # image that is being outputted by the whole net
            err = criterion(net_rec_batch, target_batch)

            err.backward(retain_graph=True)  # perform back-propagation
            optimizer.step()
 
            if epoch == _.FLAGS['epochs'] - 1:
                rand_samp = torch.randint(net_rec_batch.size()[0], (1,))
                Doppler_plots.show_diff_piezo(net_rec_batch[rand_samp, 0, :, :], target_batch[rand_samp, 0, :, :])
            
            loss_update(batch_num, epoch, err, net_rec_batch, target_batch)
            
            if batch_num % _.FLAGS['VERBOSE'] == 0 and batch_num != 0:
                training_losses.append(err.item())
                print('Training MSE Loss: ', training_losses[-1])
                min_valid_model = validation_test(model, min_valid_model, epoch, validation=True)
                model.train()
            
            if batch_num == int(_.FLAGS['NUM_SIM_SAMPS'] / _.FLAGS['batch_size']):
                rand_samp = torch.randint(net_rec_batch.size()[0], (1,))
                Doppler_plots.show_diff_piezo(net_rec_batch[rand_samp, 0, :, :], target_batch[rand_samp, 0, :, :], 'Training')
                
            
        if epoch in _.FLAGS['breaking_points']:
            Doppler_plots.plot_losses(training_losses, validation_losses, _.FLAGS['epochs'])
            min_valid_name = str(epoch) + '_epochs' + str(format(min(validation_losses), ".4f")).replace('.', '_') 
            torch.save(min_valid_model.state_dict(), min_valid_name)
            
            
    min_valid_name = str(_.FLAGS['epochs']) + str(min(validation_losses)).replace('.', '_')
    torch.save(min_valid_model.state_dict(), min_valid_name)
    Doppler_plots.plot_losses(training_losses, validation_losses, _.FLAGS['epochs'])


def validation_test(model, min_valid_model, epoch, validation):
    """
    validation test
    """
    if validation:
        with torch.set_grad_enabled(False):
            model.eval()
            valid_criterion = nn.MSELoss()
            for i, (sig_meas_valid, target_valid) in enumerate(valid_measure_gen):
                sig_meas_valid = sig_meas_valid.to(device)
                target_valid = target_valid.to(device)
                
                net_rec = model(sig_meas_valid) 
                curr_valid_err = valid_criterion(net_rec, target_valid)

                validation_losses.append(curr_valid_err.item())
                print('Validation loss: ', validation_losses[-1], 'current min valid value: ', min(validation_losses))
    
                if curr_valid_err.item() == min(validation_losses):
                    min_valid_model = copy.deepcopy(model)
                    min_valid = curr_valid_err.item()
                    print('Entered min valid scope, new min valid is now: ', min_valid, '\n\n')
                
                else:
                    print('\n')
                    
                if epoch in _.FLAGS['breaking_points']:
                    rand_samp = torch.randint(net_rec.size()[0], (1,))
                    Doppler_plots.show_diff_piezo(net_rec[rand_samp, 0, :, :], target_valid[rand_samp, 0, :, :], 'Validation')
                
    return min_valid_model


def loss_update(i, epoch, primal_loss, net_rec, target_batch):
    print('Epoch num %d out of %d, Batch num %d out of %d:' %
          (epoch, _.FLAGS['epochs'], i + 1, math.ceil(int(_.FLAGS['NUM_SIM_SAMPS'] / _.FLAGS['batch_size']) + 1)))

    print('MSE Training Loss: ', primal_loss.item(), '\n')


if __name__ == '__main__':
    model = Doppler_model.Doppler_model().to(device)
    model.apply(Doppler_model.weights_init)

    torch.autograd.set_detect_anomaly(True)

    start_time = time.time()
    train(model)
    print(time.time() - start_time)
    torch.save(model.state_dict(), 'full_trained_model')
