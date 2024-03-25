#!/usr/bin/env python
# coding: utf-8

# # Widgets

# In[ ]:


btn_0 = widgets.Button(description='Save')

def save(b):
    fig_glob.savefig('spectroscopic map at {} (meV)'.format(energy_glob), format = 'jpg',bbox_inches='tight')
btn_0.on_click(save)

btn_1 = widgets.Button(description='Save')
def save(b):
    fig_glob.savefig('Filtered map at {} (meV)'.format(energy_glob), format = 'jpg',bbox_inches='tight')
btn_1.on_click(save)

btn_2 = widgets.Button(description='Save')

def save(b):
    fig_glob.savefig('dI_dV map at {} (meV) averaged over r = {} at x = {} y = {}'.format(energy_glob,r_glob, x_glob, y_glob),
                     format = 'jpg',bbox_inches='tight')
btn_2.on_click(save)

btn_3 = widgets.Button(description='Save')
def save(b):
    fig_glob.savefig('dI_dV map averaged at {} (meV)'.format(energy_glob), format = 'jpg', bbox_inches='tight')
btn_3.on_click(save)

btn_4 = widgets.Button(description='Save')
def save(b):
    fig_glob.savefig('dI_dV map averaged at {} (meV)'.format(energy_glob), format = 'jpg', bbox_inches='tight')
btn_4.on_click(save)

btn_6 = widgets.Button(description='Save')
def save(b):
    fig_glob.savefig('dI_dV average', format = 'jpg', bbox_inches='tight')
btn_6.on_click(save)

x_widget = FloatSlider(min=0, max=len_x, step = len_x/100, continuous_update=False)

y_widget = FloatSlider(min=0, max=abs(len_y), step = len_y/100, continuous_update=False)

width_widget = FloatSlider(min = 1, max=abs(len_y)/3, step = len_y/100, continuous_update=False)

upper_cut_widget = FloatSlider(min = (bias_end+bias_start)/2, max = max(bias_end,bias_start),
                               step = abs((bias_end-bias_start))/num_points ,value =  bias_end, continuous_update=False)

lower_cut_widget = FloatSlider(min = min(bias_start,bias_end), max = (bias_end+bias_start)/2,
                               step = abs((bias_end-bias_start))/num_points ,value = bias_start, continuous_update=False)

energy_widget = FloatSlider(min= lower_cut_widget.value, max=upper_cut_widget.value, step =(bias_end-bias_start)/num_points,
                            continuous_update = False)


def update_energy_range(*args):
    energy_widget.max = upper_cut_widget.value
    energy_widget.min = lower_cut_widget.value
    

upper_cut_widget.observe(update_energy_range, 'value')
lower_cut_widget.observe(update_energy_range, 'value')

