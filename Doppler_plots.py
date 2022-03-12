import torch
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import numpy as np

import flags as _


import torch
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable


def show_diff_piezo(net_rec, original, type):
    a = net_rec
    b = original
    a = a.squeeze()
    b = b.squeeze()
    a = a.cpu().detach().numpy()
    b = b.cpu().detach().numpy()
    fig = plt.figure(figsize=(6, 6), dpi=90)
    ax1 = fig.add_subplot(1, 2, 1)
    plt.title('Piezos Rec (Normalized) 25%')
    b1 = ax1.imshow(a)

    divider1 = make_axes_locatable(ax1)
    cax1 = divider1.append_axes('right', size='5%', pad=0.05)
    fig.colorbar(b1, cax=cax1, orientation='vertical')

    ax2 = fig.add_subplot(1, 2, 2)
    plt.title('Original Piezo (Normalized)')
    b2 = ax2.imshow(b)

    divider2 = make_axes_locatable(ax2)
    cax2 = divider2.append_axes('right', size='5%', pad=0.05)
    fig.colorbar(b2, cax=cax2, orientation='vertical')

    plt.tight_layout()
    suptitle = type + 'Comparison'
    plt.suptitle(suptitle)
    
    plt.savefig(suptitle)
    # plt.show()
    



def plot_losses(train_err, valid_err, max_epochs):
    max_epochs = 10
    fig = plt.figure(figsize=(7, 7))
    ax1 = fig.add_subplot(111)
    ax2 = ax1.twiny()

    ax1.plot(range(len(train_err)), train_err, color='red', label='Training loss')
    ax1.legend(loc=0)
    # plt.legend()
    ax1.plot(range(len(valid_err)), valid_err, color='green', label='Validation loss')
    ax1.set_xlim(left=1)
    ax1.set_ylim(0, 0.1)

    ax1.legend(loc=1)
    ax1.grid(axis='y', linestyle='-')

    plt.xticks(np.arange(len(train_err)), [str(i) for i in range(len(train_err))])
    ax1.set_xlabel('Total batches Num')

    min_valid = min(valid_err)
    min_label = 'Min validation err: ' + format(min_valid, ".6f")
    ax2.plot(range(max_epochs), np.ones(max_epochs) * min_valid, linestyle='--', label='Min validation err: ' +
                                                                        format(min_valid, ".6f"))
    ax2.set_xlabel('Number of epochs')
    ax2.set_ylabel('MSE loss value')
    # ax2.cla()
    ax2.legend(loc=2)

    ax2.grid(axis='x', linestyle='-')
    plt.title('Training and Validation MSE losses')
    name = 'Losses_plots' + str(max_epochs) + '.png'
    plt.savefig(name)

    # plt.show()







