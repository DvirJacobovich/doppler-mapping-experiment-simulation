import torch
import torch.nn as nn

from Doppler_Transformer_folder import Doppler_transformer as dt

import flags as _


class Doppler_model(nn.Module):
    def __init__(self):
        super(Doppler_model, self).__init__()
        self._transformer = dt.Doppler_transformer()
        self.decoding_proj = nn.Sequential(
            # input vector size: batct_size x 1024
            nn.AlphaDropout(0.5),
            nn.Linear(_.FLAGS['measures'], 4096), # 4096 = 64 x 8 x 8
            nn.PReLU(num_parameters=1, init=0.25, device=None, dtype=None),
            nn.BatchNorm1d(1),
        
        )
        
        self.decoder2 = nn.Sequential( # input is 128 x 128
            nn.Conv2d(1, 64, kernel_size=(3, 3), stride=(1, 1), padding=1, padding_mode='reflect'), # 64 x 64 x 64
            nn.PReLU(num_parameters=64, init=0.25, device=None, dtype=None), 
            nn.BatchNorm2d(64),
            
            nn.AlphaDropout(_.FLAGS['in_between_prob']),
            
            nn.Conv2d(64, 128, kernel_size=(3, 3), stride=(1, 1), padding=1, padding_mode='reflect'), # 128 x 64 x 64
            nn.PReLU(num_parameters=128, init=0.25, device=None, dtype=None),
            nn.BatchNorm2d(128),
            
            nn.AlphaDropout(_.FLAGS['in_between_prob']),
            
            nn.Conv2d(128, 128, kernel_size=(3, 3), stride=(1, 1), padding=1, padding_mode='reflect'), # 128 x 64 x 64
            nn.PReLU(num_parameters=128, init=0.25, device=None, dtype=None),
            nn.BatchNorm2d(128),
            
            nn.AlphaDropout(_.FLAGS['in_between_prob']),
            
            nn.Conv2d(128, 128, kernel_size=(3, 3), stride=(1, 1), padding=1, padding_mode='reflect'), # 128 x 64 x 64
            nn.PReLU(num_parameters=128, init=0.25, device=None, dtype=None),
            nn.BatchNorm2d(128),
            
            nn.AlphaDropout(_.FLAGS['in_between_prob']),
            
            nn.Conv2d(128, 128, kernel_size=(3, 3), stride=(1, 1), padding=1, padding_mode='reflect'), # 128 x 64 x 64
            nn.PReLU(num_parameters=128, init=0.25, device=None, dtype=None),
            nn.BatchNorm2d(128),
            
            nn.AlphaDropout(_.FLAGS['in_between_prob']),
            
            nn.Conv2d(128, 64, kernel_size=(3, 3), stride=(1, 1), padding=1, padding_mode='reflect'), # 64 x 64 x 64
            nn.PReLU(num_parameters=64, init=0.25, device=None, dtype=None),
            nn.BatchNorm2d(64),
            
            nn.AlphaDropout(_.FLAGS['in_between_prob']),
            
            nn.Conv2d(64, 32, kernel_size=(3, 3), stride=(1, 1), padding=1, padding_mode='reflect'), # 1 x 64 x 64
            nn.PReLU(num_parameters=32, init=0.25, device=None, dtype=None),
            nn.BatchNorm2d(32),
            
            nn.AlphaDropout(_.FLAGS['in_between_prob']),
            
            nn.Conv2d(32, 16, kernel_size=(3, 3), stride=(1, 1), padding=1, padding_mode='reflect'), # 1 x 64 x 64
            nn.PReLU(num_parameters=16, init=0.25, device=None, dtype=None),
            nn.BatchNorm2d(16),
            
            nn.AlphaDropout(_.FLAGS['in_between_prob']),
            
            nn.Conv2d(16, 8, kernel_size=(3, 3), stride=(1, 1), padding=1, padding_mode='reflect'), # 1 x 64 x 64
            nn.PReLU(num_parameters=8, init=0.25, device=None, dtype=None),
            nn.BatchNorm2d(8),
            
            nn.Conv2d(8, 1, kernel_size=(3, 3), stride=(1, 1), padding=1, padding_mode='reflect'), # 1 x 64 x 64
           
        )
        
        self._bilinear = nn.Bilinear(64, 1024, 4096)
        
       
    def forward(self, measure_batch):
        x1 = self.decoding_proj(measure_batch[:, None, :]) # batch x 1 x 1024 
        x2 = self._transformer(x1)[:, None, :, :] # batch x 1x 64 x 64
        # x3 = self._bilinear(x2, x1)
        return self.decoder2(x2) 


def weights_init(m):
    classname = m.__class__.__name__
    if classname.find('Conv') != -1:
        nn.init.normal_(m.weight.data, 0.0, 0.02)

    elif classname.find('Linear') != -1:
        # apply a uniform distribution to the weights and a bias=0
        # was: m.weight.data.uniform_(0.0, 1.0)
        m.weight.data.uniform_(-0.5, 0.5)
        m.bias.data.fill_(0)

    elif classname.find('BatchNorm') != -1:
        nn.init.normal_(m.weight.data, 1.0, 0.02)
        nn.init.constant_(m.bias.data, 0)



#  self.decoder = nn.Sequential(
#             # input size: batch_size x 1 x 64 x 64
          
#             nn.Conv2d(1, _.FLAGS['channels'], kernel_size=(3, 3), stride=(1, 1), padding=1, padding_mode=_.FLAGS['PADDING_MODE']),
#             nn.PReLU(num_parameters=_.FLAGS['channels'], init=0.25, device=None, dtype=None),
#             nn.BatchNorm2d(_.FLAGS['channels']),
  
            
#             nn.Dropout(_.FLAGS['DROPOUT_PROB']),
            
#             nn.Conv2d(_.FLAGS['channels'], _.FLAGS['channels'], kernel_size=(3, 3), stride=(1, 1), padding=1, padding_mode=_.FLAGS['PADDING_MODE']),
#             nn.PReLU(num_parameters=_.FLAGS['channels'], init=0.25, device=None, dtype=None),
#             nn.BatchNorm2d(_.FLAGS['channels']),
            
            
#             nn.Dropout(_.FLAGS['DROPOUT_PROB']),
            
#             nn.Conv2d(_.FLAGS['channels'], _.FLAGS['channels'], kernel_size=(3, 3), stride=(1, 1), padding=1, padding_mode=_.FLAGS['PADDING_MODE']),
#             nn.PReLU(num_parameters=_.FLAGS['channels'], init=0.25, device=None, dtype=None),
#             nn.BatchNorm2d(_.FLAGS['channels']),
            
            
#             nn.Dropout(_.FLAGS['DROPOUT_PROB']),
            
#             nn.Conv2d(_.FLAGS['channels'], _.FLAGS['channels'], kernel_size=(3, 3), stride=(1, 1), padding=1, padding_mode=_.FLAGS['PADDING_MODE']),
#             nn.PReLU(num_parameters=_.FLAGS['channels'], init=0.25, device=None, dtype=None),
#             nn.BatchNorm2d(_.FLAGS['channels']),
            
            
#             nn.Dropout(_.FLAGS['DROPOUT_PROB']),
            
#             nn.Conv2d(_.FLAGS['channels'], _.FLAGS['channels'], kernel_size=(3, 3), stride=(1, 1), padding=1, padding_mode=_.FLAGS['PADDING_MODE']),
#             nn.PReLU(num_parameters=_.FLAGS['channels'], init=0.25, device=None, dtype=None),
#             nn.BatchNorm2d(_.FLAGS['channels']),
            
            
#             nn.Dropout(_.FLAGS['DROPOUT_PROB']),
            
#             nn.Conv2d(_.FLAGS['channels'], _.FLAGS['channels'], kernel_size=(3, 3), stride=(1, 1), padding=1, padding_mode=_.FLAGS['PADDING_MODE']),
#             nn.PReLU(num_parameters=_.FLAGS['channels'], init=0.25, device=None, dtype=None),
#             nn.BatchNorm2d(_.FLAGS['channels']),
            
            
#             nn.Dropout(_.FLAGS['DROPOUT_PROB']),
            
#             nn.Conv2d(_.FLAGS['channels'], _.FLAGS['channels'], kernel_size=(3, 3), stride=(1, 1), padding=1, padding_mode=_.FLAGS['PADDING_MODE']),
#             nn.PReLU(num_parameters=_.FLAGS['channels'], init=0.25, device=None, dtype=None),
#             nn.BatchNorm2d(_.FLAGS['channels']),
            
            
#             nn.Dropout(_.FLAGS['DROPOUT_PROB']),
            
#             nn.Conv2d(_.FLAGS['channels'], _.FLAGS['channels'], kernel_size=(3, 3), stride=(1, 1), padding=1, padding_mode=_.FLAGS['PADDING_MODE']),
#             nn.PReLU(num_parameters=_.FLAGS['channels'], init=0.25, device=None, dtype=None),
#             nn.BatchNorm2d(_.FLAGS['channels']),
            
            
#             nn.Dropout(_.FLAGS['DROPOUT_PROB']),
         
#             nn.Conv2d(_.FLAGS['channels'], _.FLAGS['channels'], kernel_size=(3, 3), stride=(1, 1), padding=1, padding_mode=_.FLAGS['PADDING_MODE']),
#             nn.PReLU(num_parameters=_.FLAGS['channels'], init=0.25, device=None, dtype=None),
#             nn.BatchNorm2d(_.FLAGS['channels']),
            
            
#             nn.Dropout(_.FLAGS['DROPOUT_PROB']),
            
#             nn.Conv2d(_.FLAGS['channels'], _.FLAGS['channels'], kernel_size=(3, 3), stride=(1, 1), padding=1, padding_mode=_.FLAGS['PADDING_MODE']),
#             nn.PReLU(num_parameters=_.FLAGS['channels'], init=0.25, device=None, dtype=None),
#             nn.BatchNorm2d(_.FLAGS['channels']),
            
            
#             nn.Dropout(_.FLAGS['DROPOUT_PROB']),
            
#             nn.Conv2d(_.FLAGS['channels'], _.FLAGS['channels'], kernel_size=(3, 3), stride=(1, 1), padding=1, padding_mode=_.FLAGS['PADDING_MODE']),
#             nn.PReLU(num_parameters=_.FLAGS['channels'], init=0.25, device=None, dtype=None),
#             nn.BatchNorm2d(_.FLAGS['channels']),
            
            
#             nn.Dropout(_.FLAGS['DROPOUT_PROB']),
            
#             nn.Conv2d(_.FLAGS['channels'], _.FLAGS['channels'], kernel_size=(3, 3), stride=(1, 1), padding=1, padding_mode=_.FLAGS['PADDING_MODE']),
#             nn.PReLU(num_parameters=_.FLAGS['channels'], init=0.25, device=None, dtype=None),
#             nn.BatchNorm2d(_.FLAGS['channels']),
            
            
#             nn.Dropout(_.FLAGS['DROPOUT_PROB']),
            
#             nn.Conv2d(_.FLAGS['channels'], _.FLAGS['channels'], kernel_size=(3, 3), stride=(1, 1), padding=1, padding_mode=_.FLAGS['PADDING_MODE']),
#             nn.PReLU(num_parameters=_.FLAGS['channels'], init=0.25, device=None, dtype=None),
#             nn.BatchNorm2d(_.FLAGS['channels']),
            
            
#             nn.Dropout(_.FLAGS['DROPOUT_PROB']),

#             nn.Conv2d(_.FLAGS['channels'], 1, kernel_size=(3, 3), stride=(1, 1), padding=1, padding_mode=_.FLAGS['PADDING_MODE']),
#             nn.PReLU(num_parameters=1, init=0.25, device=None, dtype=None),
#             nn.Linear(64, 64),

#         )
        