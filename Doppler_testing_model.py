import torch
import Doppler_model
import Doppler_plots
import numpy as np
from torch.nn import functional as F
import DopplerDataSet
import mat73
import os


# use_cuda = torch.cuda.is_available()
# device = torch.device("cuda:0" if use_cuda else "cpu")

device = "cuda:0"

# folder = 'pre-trained/with adaptive average pooling as output layer and all dropouts are 0.5'
# pretrained_model = 'min_valid'
# full_path = os.path.join(folder, pretrained_model)

full_path = 'net'  
# full_path = 'pre-trained/with_adaptive_average_pooling_as_output_layer_and_all_dropouts_are_0_5/min_valid'


model = Doppler_model.Doppler_model().to(device)
model.load_state_dict(torch.load(full_path))
# print('model device: ', model.device)


name = 'Piezo_samps_' + str(1000) + '.mat'
rootdir = './data'
train_data_struct = mat73.loadmat(os.path.join(rootdir, name))['data_struct']


tr_measures = torch.tensor(np.array(train_data_struct['meaures']), dtype=torch.float32,
                           requires_grad=True)


tr_piezo_targets = torch.tensor(np.array(train_data_struct['piezo_freqs']),
                                dtype=torch.float32, requires_grad=False)

tr_piezo_targets = torch.nn.functional.normalize(tr_piezo_targets, p=2, dim=0).unsqueeze(1)



tr_measure_set = DopplerDataSet.DopplerDataSet(tr_measures, tr_piezo_targets)


tr_params = {'batch_size': 1, 'shuffle': False}

tr_measure_gen = torch.utils.data.DataLoader(tr_measure_set, **tr_params)





for batch_num, (sig_meas_batch, target_batch) in enumerate(tr_measure_gen):
    if batch_num < 100:
        sig_meas_batch = sig_meas_batch.to(device)
        target_batch = target_batch.to(device)
        
        target_batch = torch.transpose(target_batch, dim0=2, dim1=1)
        
        model_output = model(sig_meas_batch)
        tst_err = F.mse_loss(target_batch, model_output)
        if batch_num == 20:
            Doppler_plots.show_diff_piezo(model_output, target_batch)
            
        print(tst_err.item())
    
    
        

 
 
#  def code_to_debug():
#     import pdb; pdb.set_trace()

#     device = "cuda:0"
    
#     # folder = 'pre-trained/with adaptive average pooling as output layer and all dropouts are 0.5'
#     # pretrained_model = 'min_valid'
#     # full_path = os.path.join(folder, pretrained_model)
    
#     full_path = 'pre-trained/with_adaptive_average_pooling_as_output_layer_and_all_dropouts_are_0_5/min_valid'
    
    
#     model = nv_model.NV_simulateed_data_model()
#     model.load_state_dict(torch.load(full_path))
#     _measure, _sampling_mtx, _magnetic, smp_freqs = testing_loader.testing_loader()
    
#     tst_generator = zip(_measure, _sampling_mtx, _magnetic)
    
    
#     # dataLoader = [[[img1_part1],[img1_part2],[img1_part3], label1], [[img2_part1],[img2_part2],[img2_part3], label2]]
    
    
#     for batch_num, ((sig_meas_batch, target_batch), (smp_mtx_batch, _), (physics_batch, _)) in enumerate(
#                     tst_generator):
#             sig_meas_batch, smp_mtx_batch, physics_batch = sig_meas_batch.to(device), smp_mtx_batch.to(
#                 device), physics_batch.to(device)
#             target_batch = target_batch.to(device)
#             target_batch = torch.transpose(target_batch, dim0=2, dim1=1)
    
#             model_output = model(sig_meas_batch, smp_mtx_batch, physics_batch)
#             tst_err = F.mse_loss(target_batch, model_output)
#             NV_plots.model_rec_plots(model_output, target_batch, smp_freqs, 'Testing')
#             print(tst_err.item())
        
        
        
# code_to_debug()
    