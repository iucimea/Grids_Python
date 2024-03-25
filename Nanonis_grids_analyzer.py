### Importing the Libraries

import nanonis_library
import main_functions
import median_filtering

#%% Extracting spectroscopic Data
 
#

file = "grid-spectroscopy004"
file_name = '{}.3ds'.format(file)
extract_grid = Grid(file_name)
datas = extract_grid._load_data()   #opens both the topography and the spectroscopic images
topo = datas['params']              #extracts parameters taken at the grid points 
grid_spec = datas['Current (A)']    #extracts the spectroscopic data
divider = 1
num_points=len(grid_spec[0][0])
x_points = len(grid_spec[0,:,0])
y_points = len(grid_spec[:,0,0])
grid_spec=flip(grid_spec,0)
bias_start = (topo [0][0][0] *1000)/divider  #starting point in meV
bias_end = (topo [0][0][1] *1000)/divider   #ending point in meV
len_x = (-topo[0][0][2]+topo[0][-1][2]) * 10**9 # map size in X
len_y = (-topo[0][0][3]+topo[-1][0][3])* 10**9 #map size in Y

energy_axis = linspace(bias_start, bias_end, num_points)               #finding all the energies
print("number of points in spectroscopy:", num_points,"\nBias_Start:",bias_start,"\nBias_End:", bias_end, "\nnumber of points in x:",x_points,"\nnumber of points in y:", y_points, "\nLenght in x:",len_x,'nm',"\nLenght in y:",len_y,'nm')


#%% Import widgets

import widgets


#%% Read params


params = {'len_x':len_x, 'len_y': len_y, 'x_points':x_points , 'y_points': y_points, 'bias_start':bias_start, 'bias_end':bias_end, 'num_points': num_points}


#%% Median filtering


def median_filtering_map(data, energy, n, m_tresh, params):
    data_median = median_filtering(data[index_of_energy(energy, params)], n, m_tresh)
    global m_glob
    m_glob = m_tresh
    figure(figsize = (10,10))

    imshow(data_median, extent=[0, len_x, 0, len_y])
    xlabel('X', fontsize=16)
    ylabel('Y', fontsize=16)


# In[19]:


grid_map =(asarray(make_map(grid_spec, params)))


# In[20]:


interact(median_filtering_map, data = fixed(grid_map), 
         energy = FloatSlider(min = min(bias_start, bias_end), max = max(bias_start ,bias_end),
         step =0.1, continuous_update = False),
         n = fixed(3),
         m_tresh = FloatSlider(min=0.01, max=2,step =0.01, continuous_update=False), 
         params = fixed(params))


# In[ ]:


all_maps_median_filtered = median_filtering_all_maps(grid_map, 3, 1)

np.save('All_maps_median_filtered_{}.npy'.format(file),all_maps_median_filtered)
all_maps_median_filtered = np.load('All_maps_median_filtered_{}.npy'.format(file))

imshow(all_maps_median_filtered[250,:,:])


# # Fourier transform

# In[ ]:


def fourier(data):
    f=np.fft.fft2(data)
    fshift = np.fft.fftshift(f)         
    return fshift

    
def image_fourier (data, sizemap, m_tresh):
    fourier_ =  median_filtering(log(absolute(fourier(data))), 3, m_tresh)
    fig = figure()
    imshow(fourier_)


def map_in_energy_fourier(data, n, energy, size, _range, x_points, y_points, m_tresh):
    fig=figure(figsize = (size, size))
    data_1 = data[index_of_energy(energy, params)]
    e_plus = index_of_energy(energy + _range, params)
    e_minus = index_of_energy(energy - _range, params)
    average = average_data(data, e_minus, e_plus)
    median = median_filtering (average, n, m_tresh)
    
    derived_fourier = log(absolute(fourier(data_1)))
    derived_averaged_fourier = log(absolute(fourier(average)))
    derived_averaged_fourier_median = log(absolute(fourier(median)))
    
    subplot(1,2,1)
    imshow(data_1)
    xlabel('X', fontsize=16)
    ylabel('Y', fontsize=16)
    title('derived map', fontsize=16) 
    
    subplot(1,2,2)
    imshow(derived_fourier)
    xlabel('X', fontsize=16)
    ylabel('Y', fontsize=16)
    title('derived map fourier', fontsize=16) 
    
    
    global energy_glob
    nergy_glob = energy
    global fig_glob 
    fg_glob = fig       


interact(map_in_energy_fourier, data = fixed(grid_map), n = fixed(3),
         energy = FloatSlider(min = bias_end, max = bias_start, step =0.1, continuous_update = False), 
         size= IntSlider (min=20, max=30, step = 1, continuous_update=False), 
         x_points= IntSlider(min=0, max=300, step = 1, continuous_update=False), 
         y_points=IntSlider(min=0, max=300, step = 1, continuous_update=False), 
         _range = IntSlider(min=1, max=50, step = 1, continuous_update=False) , 
         m_tresh = FloatSlider(min=0.1, max=1, step = 0.1, continuous_update=False),
         params = fixed(params))


# In[ ]:


def map_in_energy_derive_fourier (data, energy, figure_size, m_tresh, params):  
    e = index_of_energy (energy, params)
    fig = figure(figsize = (figure_size, figure_size))
    median = median_filtering(data[e], 3, m_tresh)
    mean = average(data, (1,2))
    #mean = mean*(10**12)
    
    median_averaged = substract_average(median)
    standard = np.std(median_averaged)
    
    fourier_ = (median_filtering(log(absolute(fourier(median_averaged))), 3, m_tresh))
    fourier_1 = fourier_- np.min(fourier_)
    standard_fourier = np.std(fourier_1)
    
    fig = figure()
    fig.set_figheight(figure_size)
    fig.set_figwidth(figure_size*3)
    fig.subplots_adjust(hspace=0.4, wspace=0.4)
    
    imshow(fourier_1 ,extent=[-2/params['len_x'],2/params['len_x'],-2/params['len_y'],2/params['len_y']])
    xlabel('X (1/nm)', fontsize=50,y=1)
    ylabel('Y (1/nm)', fontsize=50)
    #title(''.join('fourier transform of derived energy map at {} meV'.format(round(energy))), fontsize=50, y= 1.07)
    xticks(fontsize=50,rotation = 45)
    yticks(fontsize=50)

    global fig_glob 
    fig_glob = fig
    global energy_glob
    energy_glob = energy 
    global fourier_glob
    fourier_glob = fourier_ 


interact(map_in_energy_derive_fourier, data = fixed(grid_map), 
         energy = FloatSlider(min= bias_end, max=bias_start, step =0.1, continuous_update = False),
         figure_size= IntSlider (min=20, max=30, step = 1, continuous_update=False), 
         x= IntSlider(min=0, max=300, step = 1, continuous_update=False), 
         y=IntSlider(min=0, max=300, step = 1, continuous_update=False), 
         r = IntSlider(min=1, max=50, step = 1, continuous_update=False) , 
         m_tresh = FloatSlider(min=0.1, max=1, step = 0.1, continuous_update=False),
         params = fixed(params))


# # Average spectrum

# In[ ]:


grid_map =(asarray(make_map(grid_spec, params)))


# In[ ]:


average_spectrum = np.mean(grid_map, axis = (1,2))
energy_axis = linspace(params['bias_start'],params['bias_end'],params['num_points'])

plot(energy_axis, average_spectrum)


# # Derivative average spectrum 

# In[ ]:


plot(energy_axis[10:-10], derivative_1d(average_spectrum, params)[10:-10])


# # Gaussian filtering

# In[ ]:


spec_filtered = make_spec(all_maps_filtered, params)


# In[ ]:


interact(Gaussian_filtering_plot, data = fixed(grid_spec), 
         x= FloatSlider(min=0, max = len_x, step = len_x/100, continuous_update=False), 
         y=FloatSlider(min=0, max=len_y, step = len_y/100, continuous_update=False), 
         sigma=IntSlider(min=1, max=15, step = 1, continuous_update=False), num_points = fixed(256), params = fixed(params))


# In[ ]:


Filtered_data = all_gaussian(grid_spec, sigma_global)
np.save('Gaussian_filtered_Data_{}.npy'.format(file), Filtered_data)
Filtered_data = np.load('Gaussian_filtered_Data_{}.npy'.format(file))
Filtered_map = make_map(Filtered_data, params)


# # Shows the Gaussian filtered map and median filters it again

# In[ ]:


interact(median_filtering_map, data = fixed(Filtered_map), size= IntSlider (min=15, max=20, step = 1, continuous_update=False),
         energy = FloatSlider(min = min(bias_start, bias_end), max = max(bias_start,bias_end),
         step =0.1, continuous_update = False), n = fixed(3),
         m_tresh = FloatSlider(min=0.01, max=2,step =0.01, continuous_update=False), params = fixed(params))


# # First derivative

# In[ ]:


derived_data = all_derive(Filtered_data, params)


# In[ ]:


np.save('All_drive_{}.npy'.format(file), derived_data)


# In[ ]:


derived_data = np.load('All_Drive_{}.npy'.format(file))


# In[ ]:


derived_map = make_map(derived_data, params)
plot(energy_axis[40:-40], np.mean(derived_map, (1,2))[40:-40])


# In[ ]:


interact(map_in_energy_derived_with_average_spectra, data = fixed(derived_map), 
         energy = FloatSlider(min= min(bias_start,bias_end), 
         max=max(bias_start,bias_end)), figure_size = IntSlider (min=20, max=30, step = 1, continuous_update=False), 
         x = x_widget, y = y_widget,
         r = IntSlider(min=1, max=50, step = 1, continuous_update=False),
         m_tresh = FloatSlider(min=0.1, max=1, step = 0.1, continuous_update=False),
         lower_cut = lower_cut_widget, upper_cut = upper_cut_widget, 
         bar_lower_percentage= FloatSlider(min=0, max=20, step = 0.1, continuous_update=False), 
         bar_higher_percentage= FloatSlider(min=0, max=20, step = 0.1, continuous_update=False), params = fixed(params))


# # Second derivative

# In[ ]:


second_derivativeData = all_derive(derived_data, params)
second_derivativeMap = make_map(second_derivativeData, params)
plot(energy_axis[40:-40], np.mean(second_derivativeMap, (1,2))[40:-40])


# In[ ]:


interact(map_in_energy_derived_with_average_spectra, data = fixed(second_derivativeMap), 
         energy = FloatSlider(min = min(bias_start,bias_end), 
         max=max(bias_start,bias_end)), figure_size = IntSlider (min=20, max=30, step = 1, continuous_update=False), 
         x = x_widget, y = y_widget,
         r = IntSlider(min=1, max=50, step = 1, continuous_update=False),
         m_tresh = FloatSlider(min=0.1, max=1, step = 0.1, continuous_update=False),
         lower_cut = lower_cut_widget, upper_cut = upper_cut_widget, 
         bar_lower_percentage= FloatSlider(min=0, max=20, step = 0.1, continuous_update=False), 
         bar_higher_percentage= FloatSlider(min=0, max=20, step = 0.1, continuous_update=False), params = fixed(params))


# In[ ]:


interact(map_in_energy_derived, data = fixed(derived_map), energy = FloatSlider(min = min(bias_start, bias_end), 
         max = max(bias_start, bias_end), step =0.05, continuous_update = False), 
         figure_size= IntSlider (min=20, max=30, step = 1, continuous_update=False), 
         x= IntSlider(min=0, max=300, step = 1, continuous_update=False), 
         y=IntSlider(min=0, max=300, step = 1, continuous_update=False), 
         r = IntSlider(min=1, max=50, step = 1, continuous_update=False) , 
         m_tresh = FloatSlider(min=0.1, max=1, step = 0.1, continuous_update=False), 
         bar_lower_percentage= FloatSlider(min=0, max=20, step = 0.1, continuous_update=False), 
         bar_higher_percentage= FloatSlider(min=0, max=20, step = 0.1, continuous_update=False), 
         params = fixed(params))


# # Spectroscopic map after filtering and derivating

# In[ ]:


interact (map_in_energy_spectra_mean, data = fixed(derived_map), energy = energy_widget, 
          figure_size= IntSlider (min=20, max=30, step = 1, continuous_update=False), 
          x= x_widget, y=y_widget, r = FloatSlider(min=len_x/x_points+0.1, max=min(len_x,len_y), 
          step = len_x/100, continuous_update=False) , 
          m_tresh = FloatSlider(min=0.1, max=1, step = 0.1, continuous_update=False), 
          sigma = IntSlider(min=1, max=30, step = 1, continuous_update=False), 
          lower_cut=lower_cut_widget, upper_cut=upper_cut_widget, 
          bar_lower_percentage= FloatSlider(min=0, max=20, step = 0.1, continuous_update=False), 
          bar_higher_percentage= FloatSlider(min=0, max=20, step = 0.1, continuous_update=False),params = fixed(params))


# In[ ]:


def map_in_energy_spectra_line_profile (data, energy, figure_size, x, y, r, m_tresh,sigma,lower_cut,upper_cut,bar_lower_percentage,bar_higher_percentage,width, params):
    energy_axis = linspace(params['bias_start'],params['bias_end'],params['num_points'])                             
    x_1 = x_i(x,params)
    y_1 = y_i(y,params)
    r_1 = x_i(r,params)
    y_1_min = y_i(y - width/2 ,params)
    y_1_max = y_i(y + width/2 ,params)
    mean = mean_spectrum_circ(data,x_1,y_1,r_1, params['x_points'],params['y_points'],)
    filtered_mean = gaussian_filter1d(mean, sigma)
    e = index_of_energy(energy,params)
    theta = np.linspace(0, 2*np.pi, 100)
    radius = r
    a = radius*np.cos(theta)
    b = radius*np.sin(theta)
    median = median_filtering (data[e], 3, m_tresh)
    average_ = average(median[y_1_max:y_1_min],axis=0)
    fig = figure()
    fig.subplots_adjust(hspace=0.4, wspace=0.4) 
    fig.set_figheight(figure_size)
    fig.set_figwidth(figure_size*3)
    upper_index = index_of_energy(upper_cut,params)
    lower_index = index_of_energy(lower_cut,params)
  
    
    subplot(2,2,1)
    im = imshow(median, extent=[0,params['len_x'],0,params['len_y']])
    plot(a+x,b+y,color="red", linewidth=5)
    axhline(y = y, color = 'w', linestyle = '-', linewidth=5) 
    xlabel('X (nm)', fontsize=50, y=1.01)
    ylabel('Y (nm)', fontsize=50, x = 1.01)
    title(''.join('derived energy map at {} meV'.format(energy)), fontsize=50, y= 1.1)
    xticks(fontsize=50)
    yticks(fontsize=50)
    cbar =colorbar(im, boundaries=linspace(np.min(median)*bar_lower_percentage/100,np.max(median)*(100-bar_higher_percentage)/100,1000),shrink=1,ticks=linspace(np.min(median)*bar_lower_percentage/100,np.max(median)*(100-bar_higher_percentage)/100,10))
    cbar.minorticks_on()
    cbar.ax.tick_params(labelsize=50)
    
    subplot(2,2,2)
    plot(energy_axis[min(lower_index,upper_index):max(lower_index,upper_index)],filtered_mean[min(lower_index,upper_index):max(lower_index,upper_index)], linewidth=6)
    plot(energy_axis[e], filtered_mean[e], 'ro', markersize=25)
    xlabel('E (meV)', fontsize=50)
    ylabel('dI/dV (pS)', fontsize=50)
    title(''.join('Average Spectra of the map at the circle'.format(energy)), fontsize=50, y= 1.07)
    xticks(fontsize=50)
    yticks(fontsize=50)
    
    subplot(2,2,3)
    plot(arange(params['x_points']), average_)

    xlabel('X (nm)', fontsize=50)
    ylabel('dI/dV (pS)', fontsize=50)
    title(''.join('Average Spectra of the map along the line'.format(energy)), fontsize=50, y= 1.07)
    xticks(fontsize=50)
    yticks(fontsize=50)    
    plt.show()


# In[ ]:


def map_in_energy_spectra_mean_3circs (data, energy, figure_size, x_1, y_1, r_1,x_2,y_2,r_2,x_3,y_3,r_3, m_tresh,sigma,lower_cut,upper_cut,bar_lower_percentage,bar_higher_percentage,n_1,n_2,n_3,params):
    
    x_1_i = x_i(x_1,params)
    y_1_i = y_i(y_1,params)
    r_1_i = x_i(r_1,params)
    x_2_i = x_i(x_2,params)
    y_2_i = y_i(y_2,params)
    r_2_i = x_i(r_2,params)
    x_3_i = x_i(x_3,params)
    y_3_i = y_i(y_3,params)
    r_3_i = x_i(r_3,params)
    
    e = index_of_energy(energy,params)
    theta = np.linspace(0, 2*np.pi, 100)
    
    
    mean_1 = mean_spectrum_circ(data,x_1_i,y_1_i,r_1_i,params['x_points'], params['y_points'])
    filtered_mean_1 = gaussian_filter1d(mean_1,sigma)
    fig = figure(figsize = (figure_size,figure_size))
    radius_1 = r_1
    a_1 = radius_1*np.cos(theta)
    b_1= radius_1*np.sin(theta)
    
    mean_2 = mean_spectrum_circ(data,x_2_i,y_2_i,r_2_i,params['x_points'], params['y_points'])
    filtered_mean_2 = gaussian_filter1d(mean_2,sigma)
    radius_2 = r_2
    a_2 = radius_2*np.cos(theta)
    b_2 = radius_2*np.sin(theta)    
 
    mean_3 = mean_spectrum_circ(data,x_3_i,y_3_i,r_3_i,params['x_points'], params['y_points'])
    filtered_mean_3= gaussian_filter1d(mean_3,sigma)
    radius_3 = r_3
    a_3 = radius_3*np.cos(theta)
    b_3 = radius_3*np.sin(theta)

    median = median_filtering (data[e], 3, m_tresh)
    fig = figure(figsize = (60,20))
    fig.subplots_adjust(hspace=0.1, wspace=0.1) 
    upper_index = index_of_energy(upper_cut,params)
    lower_index = index_of_energy(lower_cut,params)
  
    
    subplot(1,2,1)
    im = imshow(median, extent=(0,params['len_x'],0,params['len_y']))
    plot(a_1+x_1,b_1+y_1,color="aqua", linewidth=8)
    plot(a_2+x_2,b_2+y_2,color="red", linewidth=8)
    plot(a_3+x_3,b_3+y_3,color="lightgreen", linewidth=8)
    xlabel('X (nm)', fontsize=50, y = 1.01)
    ylabel('Y (nm)', fontsize=50, x = 1.01)
    title(''.join('dI/dV map at {} meV'.format("{:.2f}".format(energy))), fontsize=50, y= 1.1)
    xticks(fontsize=50)
    yticks(fontsize=50)
    cbar =colorbar(im, boundaries=linspace(np.min(median)*bar_lower_percentage/100,np.max(median)*(100-bar_higher_percentage)/100,1000),shrink=0.3,ticks=linspace(np.min(median)*bar_lower_percentage/100,np.max(median)*(100-bar_higher_percentage)/100,5).astype(int))
    cbar.minorticks_on()
    cbar.ax.tick_params(labelsize=50)
  
        
    subplot(1,2,2)
    plot(energy_axis[min(lower_index,upper_index):max(lower_index,upper_index)],n_1*asarray(filtered_mean_1[min(lower_index,upper_index):max(lower_index,upper_index)]),color="aqua", linewidth=8)
    plot(energy_axis[min(lower_index,upper_index):max(lower_index,upper_index)],n_2*asarray(filtered_mean_2[min(lower_index,upper_index):max(lower_index,upper_index)]),color="red", linewidth=8)    
    plot(energy_axis[min(lower_index,upper_index):max(lower_index,upper_index)],n_3*asarray(filtered_mean_3[min(lower_index,upper_index):max(lower_index,upper_index)]),color="lightgreen", linewidth=8)    
    margins(x=0)
    xlabel('E (meV)', fontsize=50)
    ylabel('dI/dV (pS)', fontsize=50)
    title(''.join('Average Spectra of the map at the circle'.format(energy)), fontsize=50, y= 1.07)
    xticks(fontsize=50)
    yticks(fontsize=50)
    ylim(0)
    subplots_adjust(left=0.1,bottom=0.1, right=0.9, top=0.9, wspace=0.4, hspace=0.4)
    global fig_glob 
    fig_glob = fig
    global energy_glob
    energy_glob = energy


# In[ ]:


interact(map_in_energy_spectra_mean_3circs, data = fixed(derived_map), energy = energy_widget, 
         figure_size= IntSlider (min=20, max=30, step = 1, continuous_update=False), 
         x_1= FloatSlider(min=0, max=len_x, step = len_x/100, continuous_update=False), 
         y_1=FloatSlider(min=0, max=len_y, step = len_y/100, continuous_update=False), 
         r_1 = FloatSlider(min=len_x/x_points+0.5, max=min(len_x,len_y), 
         step = len_x/100, continuous_update=False),
         x_2= FloatSlider(min=0, max=len_x, step = len_x/100, continuous_update=False),
         y_2=FloatSlider(min=0, max=len_y, step = len_y/100, continuous_update=False),
         r_2 = FloatSlider(min=len_x/x_points+0.5, max=min(len_x,len_y),
         step = len_x/100, continuous_update=False),
         x_3= FloatSlider(min=0, max=len_x, step = len_x/100, continuous_update=False), 
         y_3=FloatSlider(min=0, max=len_y, step = len_y/100, continuous_update=False), 
         r_3 = FloatSlider(min=len_x/x_points+0.5, max=min(len_x,len_y), 
         step = len_x/100, continuous_update=False), 
         m_tresh = FloatSlider(min=0.1, max=10, step = 0.1, continuous_update=False), 
         sigma=IntSlider(min=1, max=30, step = 1, continuous_update=False),
         lower_cut=lower_cut_widget,upper_cut=upper_cut_widget, 
         bar_lower_percentage= FloatSlider(min=0, max=20, step = 0.1, continuous_update=False), 
         bar_higher_percentage= FloatSlider(min=0, max=20, step = 0.1, continuous_update=False), 
         n_1= FloatSlider(min=0.8, max=1, step = 0.01, continuous_update=False), 
         n_2= FloatSlider(min=0.8, max=1, step = 0.01, continuous_update=False), 
         n_3= FloatSlider(min=0.8, max=1, step = 0.01, continuous_update=False),params = fixed(params))

#btn_5


# In[ ]:


interact(map_in_energy_derived, data = fixed(derived_map), energy = FloatSlider(min= bias_start, max=bias_end, step =0.05, continuous_update = False), figure_size= IntSlider (min=20, max=30, step = 1, continuous_update=False), 
         x= IntSlider(min=0, max=300, step = 1, continuous_update=False), y=IntSlider(min=0, max=300, step = 1, continuous_update=False), r = IntSlider(min=1, max=50, step = 1, continuous_update=False) , m_tresh = FloatSlider(min=0.1, max=1, step = 0.1, continuous_update=False),bar_lower_percentage= FloatSlider(min=0, max=20, step = 0.1, continuous_update=False),bar_higher_percentage= FloatSlider(min=0, max=20, step = 0.1, continuous_update=False),params = fixed(params))

btn_4


# In[ ]:


interact(map_in_energy_derived_with_average_spectra, data = fixed(derived_map), 
         energy = energy_widget, figure_size = IntSlider (min=20, max=30, step = 1, continuous_update=False), 
         x = x_widget, y = y_widget,
         r = IntSlider(min=1, max=50, step = 1, continuous_update=False),
         m_tresh = FloatSlider(min=0.1, max=1, step = 0.1, continuous_update=False),
         lower_cut = lower_cut_widget, upper_cut = upper_cut_widget, 
         bar_lower_percentage= FloatSlider(min=0, max=20, step = 0.1, continuous_update=False), 
         bar_higher_percentage= FloatSlider(min=0, max=20, step = 0.1, continuous_update=False), params = fixed(params))


# In[ ]:


interact(map_in_energy_derive_fourier, data = fixed(derived_map), energy = FloatSlider(min= bias_start, max=bias_end, step =0.1, continuous_update = False), figure_size= IntSlider (min=20, max=30, step = 1, continuous_update=False),
         x= IntSlider(min=0, max=300, step = 1, continuous_update=False), y=IntSlider(min=0, max=300, step = 1, continuous_update=False), r = IntSlider(min=1, max=50, step = 1, continuous_update=False) , m_tresh = FloatSlider(min=0.1, max=1, step = 0.1, continuous_update=False),params = fixed(params))


# In[ ]:


def normalizing_map (data):
    norm = []
    for i in range (len(data)):
        norm.append(data[i]/sum(data[i]))
    return norm


# In[ ]:


normalized_map = asarray(normalizing_map (derived_map))


# # Treating topographic data

# In[15]:


topography = topo[:,:,4]
topography = flip(topography,0)
slope_corrected_topography = substract_average(topography)
figure(figsize = (10,10))

params = {'len_x':len_x, 'len_y': len_y, 'x_points':x_points , 'y_points': y_points, 'bias_start':bias_start,
          'bias_end':bias_end, 'num_points': num_points, 'slope_corrected_topography' : slope_corrected_topography}

imshow(slope_corrected_topography, extent=[0, len_x, 0, len_y])


# # Polynomyal slope

# In[21]:


interact(substract_polynomyal_slope,
         data = fixed(slope_corrected_topography), 
         num_points_x = fixed(x_points), 
         num_points_y = fixed(y_points),
         order = IntSlider(min=1,max=10,step=1, continuous_update=False),
         size = IntSlider (min=5, max=15, step = 1, continuous_update=False),
         len_x = fixed(len_x), len_y = fixed(len_y),params = fixed(params))

