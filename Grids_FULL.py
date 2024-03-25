#!/usr/bin/env python
# coding: utf-8

# In[ ]:


# Importing the Libraries

get_ipython().run_line_magic('run', 'nanonis_library')
get_ipython().run_line_magic('run', 'main_functions.ipynb')
get_ipython().run_line_magic('run', 'median_filtering.ipynb')

# Extracting spectroscopic Data

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


# Import widgets

get_ipython().run_line_magic('run', 'widgets.ipynb')

# Treating topographic data

topography = topo[:,:,4]
topography = flip(topography,0)
slope_corrected_topography = substract_average(topography)
figure(figsize = (10,10))

params = {'len_x':len_x, 'len_y': len_y, 'x_points':x_points , 'y_points': y_points, 'bias_start':bias_start,
          'bias_end':bias_end, 'num_points': num_points, 'slope_corrected_topography' : slope_corrected_topography}

imshow(slope_corrected_topography, extent=[0, len_x, 0, len_y])

# Median filtering

def median_filtering_map(data, energy, n, m_tresh, params):
    data_median = median_filtering(data[index_of_energy(energy, params)], n, m_tresh)
    global m_glob
    m_glob = m_tresh
    figure(figsize = (10,10))

    imshow(data_median, extent=[0, len_x, 0, len_y])
    xlabel('X', fontsize=16)
    ylabel('Y', fontsize=16)

interact(median_filtering_map, data = fixed(grid_map), 
         energy = FloatSlider(min = min(bias_start, bias_end), max = max(bias_start ,bias_end),
         step =0.1, continuous_update = False),
         n = fixed(3),
         m_tresh = FloatSlider(min=0.01, max=2,step =0.01, continuous_update=False), 
         params = fixed(params))


all_maps_median_filtered = median_filtering_all_maps(grid_map, 3, 1)

np.save('All_maps_median_filtered_{}.npy'.format(file),all_maps_median_filtered)
all_maps_median_filtered = np.load('All_maps_median_filtered_{}.npy'.format(file))

imshow(all_maps_median_filtered[250,:,:])

# Fourier transform

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


# Average spectrum

grid_map =(asarray(make_map(grid_spec, params)))

average_spectrum = np.mean(grid_map, axis = (1,2))
energy_axis = linspace(params['bias_start'],params['bias_end'],params['num_points'])

plot(energy_axis, average_spectrum)

# Derivative average spectrum 

plot(energy_axis[10:-10], derivative_1d(average_spectrum, params)[10:-10])

# Gaussian filtering

spec_filtered = make_spec(all_maps_filtered, params)

interact(Gaussian_filtering_plot, data = fixed(grid_spec), 
         x= FloatSlider(min=0, max = len_x, step = len_x/100, continuous_update=False), 
         y=FloatSlider(min=0, max=len_y, step = len_y/100, continuous_update=False), 
         sigma=IntSlider(min=1, max=15, step = 1, continuous_update=False), num_points = fixed(256), params = fixed(params))

Filtered_data = all_gaussian(grid_spec, sigma_global)
np.save('Gaussian_filtered_Data_{}.npy'.format(file), Filtered_data)
Filtered_data = np.load('Gaussian_filtered_Data_{}.npy'.format(file))
Filtered_map = make_map(Filtered_data, params)

# Shows the Gaussian filtered map and median filters it again

interact(median_filtering_map, data = fixed(Filtered_map), size= IntSlider (min=15, max=20, step = 1, continuous_update=False),
         energy = FloatSlider(min = min(bias_start, bias_end), max = max(bias_start,bias_end),
         step =0.1, continuous_update = False), n = fixed(3),
         m_tresh = FloatSlider(min=0.01, max=2,step =0.01, continuous_update=False), params = fixed(params))

# First derivative

derived_data = all_derive(Filtered_data, params)

np.save('All_drive_{}.npy'.format(file), derived_data)

derived_data = np.load('All_Drive_{}.npy'.format(file))

derived_map = make_map(derived_data, params)
plot(energy_axis[40:-40], np.mean(derived_map, (1,2))[40:-40])

interact(map_in_energy_derived_with_average_spectra, data = fixed(derived_map), 
         energy = FloatSlider(min= min(bias_start,bias_end), 
         max=max(bias_start,bias_end)), figure_size = IntSlider (min=20, max=30, step = 1, continuous_update=False), 
         x = x_widget, y = y_widget,
         r = IntSlider(min=1, max=50, step = 1, continuous_update=False),
         m_tresh = FloatSlider(min=0.1, max=1, step = 0.1, continuous_update=False),
         lower_cut = lower_cut_widget, upper_cut = upper_cut_widget, 
         bar_lower_percentage= FloatSlider(min=0, max=20, step = 0.1, continuous_update=False), 
         bar_higher_percentage= FloatSlider(min=0, max=20, step = 0.1, continuous_update=False), params = fixed(params))

# Second derivative

second_derivativeData = all_derive(derived_data, params)
second_derivativeMap = make_map(second_derivativeData, params)
plot(energy_axis[40:-40], np.mean(second_derivativeMap, (1,2))[40:-40])

interact(map_in_energy_derived_with_average_spectra, data = fixed(second_derivativeMap), 
         energy = FloatSlider(min = min(bias_start,bias_end), 
         max=max(bias_start,bias_end)), figure_size = IntSlider (min=20, max=30, step = 1, continuous_update=False), 
         x = x_widget, y = y_widget,
         r = IntSlider(min=1, max=50, step = 1, continuous_update=False),
         m_tresh = FloatSlider(min=0.1, max=1, step = 0.1, continuous_update=False),
         lower_cut = lower_cut_widget, upper_cut = upper_cut_widget, 
         bar_lower_percentage= FloatSlider(min=0, max=20, step = 0.1, continuous_update=False), 
         bar_higher_percentage= FloatSlider(min=0, max=20, step = 0.1, continuous_update=False), params = fixed(params))

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


# Spectroscopic map after filtering and derivating

interact (map_in_energy_spectra_mean, data = fixed(derived_map), energy = energy_widget, 
          figure_size= IntSlider (min=20, max=30, step = 1, continuous_update=False), 
          x= x_widget, y=y_widget, r = FloatSlider(min=len_x/x_points+0.1, max=min(len_x,len_y), 
          step = len_x/100, continuous_update=False) , 
          m_tresh = FloatSlider(min=0.1, max=1, step = 0.1, continuous_update=False), 
          sigma = IntSlider(min=1, max=30, step = 1, continuous_update=False), 
          lower_cut=lower_cut_widget, upper_cut=upper_cut_widget, 
          bar_lower_percentage= FloatSlider(min=0, max=20, step = 0.1, continuous_update=False), 
          bar_higher_percentage= FloatSlider(min=0, max=20, step = 0.1, continuous_update=False),params = fixed(params))


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

interact(map_in_energy_derived, data = fixed(derived_map), energy = FloatSlider(min= bias_start, max=bias_end, step =0.05, continuous_update = False), figure_size= IntSlider (min=20, max=30, step = 1, continuous_update=False), 
         x= IntSlider(min=0, max=300, step = 1, continuous_update=False), y=IntSlider(min=0, max=300, step = 1, continuous_update=False), r = IntSlider(min=1, max=50, step = 1, continuous_update=False) , m_tresh = FloatSlider(min=0.1, max=1, step = 0.1, continuous_update=False),bar_lower_percentage= FloatSlider(min=0, max=20, step = 0.1, continuous_update=False),bar_higher_percentage= FloatSlider(min=0, max=20, step = 0.1, continuous_update=False),params = fixed(params))

btn_4

interact(map_in_energy_derived_with_average_spectra, data = fixed(derived_map), 
         energy = energy_widget, figure_size = IntSlider (min=20, max=30, step = 1, continuous_update=False), 
         x = x_widget, y = y_widget,
         r = IntSlider(min=1, max=50, step = 1, continuous_update=False),
         m_tresh = FloatSlider(min=0.1, max=1, step = 0.1, continuous_update=False),
         lower_cut = lower_cut_widget, upper_cut = upper_cut_widget, 
         bar_lower_percentage= FloatSlider(min=0, max=20, step = 0.1, continuous_update=False), 
         bar_higher_percentage= FloatSlider(min=0, max=20, step = 0.1, continuous_update=False), params = fixed(params))

interact(map_in_energy_derive_fourier, data = fixed(derived_map), energy = FloatSlider(min= bias_start, max=bias_end, step =0.1, continuous_update = False), figure_size= IntSlider (min=20, max=30, step = 1, continuous_update=False),
         x= IntSlider(min=0, max=300, step = 1, continuous_update=False), y=IntSlider(min=0, max=300, step = 1, continuous_update=False), r = IntSlider(min=1, max=50, step = 1, continuous_update=False) , m_tresh = FloatSlider(min=0.1, max=1, step = 0.1, continuous_update=False),params = fixed(params))


def normalizing_map (data):
    norm = []
    for i in range (len(data)):
        norm.append(data[i]/sum(data[i]))
    return norm


normalized_map = asarray(normalizing_map (derived_map))

# Polynomyal slope

interact(substract_polynomyal_slope,
         data = fixed(slope_corrected_topography), 
         num_points_x = fixed(x_points), 
         num_points_y = fixed(y_points),
         order = IntSlider(min=1,max=10,step=1, continuous_update=False),
         size = IntSlider (min=5, max=15, step = 1, continuous_update=False),
         len_x = fixed(len_x), len_y = fixed(len_y),params = fixed(params))


# Cross Correlation

start = -500
end = -510
start_E = -1
points = absolute((start-end))
range_ = 20


def cross_multiplication (data, E_1, E_2, params):
    return mean((data[index_of_energy(E_1, params)])*(data[index_of_energy(E_2, params)]))
    
def cross_definition_simple (data,E_1, E_2, points,params):
    array_E_1_E_2 = []
    cross_array = []
    E_array = linspace(E_1, E_2, points)
    for i in E_array:
        array_E_1_E_2.append(np.absolute((E_1-i)))    
        cross_array.append(cross_multiplication (data, E_1, i,params))
    return array_E_1_E_2, cross_array
    print(E_array)
    
def cross_definition (data, E_1, E_2, points,range_, start_E,params):
    array_E_1_E_2 = []
    cross_array = []
    E_array = linspace(E_1, E_2, points)
    for i in E_array:
        if i - min(E_array) > start_E :
            array_E_1_E_2.append(np.absolute((E_1-i)))
            for j in range (range_):
                cross_array.append(cross_multiplication (data, E_1+j, i+j,params))
    return array_E_1_E_2, cross_array

def cross_correlation_map (data, E_1_start, E_1_end, E_2_start, E_2_end, points_1, points_2,params):
    cross_array = []
    E_1_array = linspace(E_1_start, E_1_end, points_1)
    E_2_array = linspace(E_2_start, E_2_end, points_2)
    for i in E_1_array:
        for j in E_2_array:
            cross_array.append(cross_multiplication (data, i, j, params))
    cross_map = reshape(cross_array,(len( E_1_array),len( E_2_array)))
    return E_1_array, E_2_array, cross_map



a = linspace(-500,-400, 100)

cross_correlation = cross_correlation_map (normalized_map, -400, -500, -400, -500, 50, 50, params) 

cross = cross_definition_simple (normalized_map, start, end, points, params)

E_start = -400
E_end = -600
steps = -10 
points = 10
D_2_array = []
E_array = []

for i in range(absolute(E_start), absolute(E_end-steps)):
    cross =cross_definition_simple (normalized_map, -i, -i+steps, points,params)
    D_2_array.append(2+((linregress(log(cross[0][1:-1]),log(cross[1][1:-1])))[0])*2)
    E_array.append(-i)

plot(E_array,D_2_array)

a = cross_definition(normalized_map, start, end, points,range_, start_E,params)
b = asarray(a[1]).reshape(int(absolute(start-end)), range_)
#b = asarray(a[1]).reshape(int(absolute(start-end)-start_E), range_)
c = sum(b, axis=1)

plot((a[0][1:-1]), ((c)[1:-1]))
yscale("log")
xscale("log")

from scipy.stats import linregress
(linregress(log(cross[0][1:-1]),log(cross[1][1:-1])))[0]

def sum_eta (data, energy, L_array, q, params) :
    sum_array = []
    for i in L_array:
        sum_array.append(sum_box_mu(data[index_of_energy(energy,params)], i, q, params))
    return sum_array

def sum_box_mu(data, L, q,params):
    mu_matrix = []
    temp = []
    for i in range (0, len(data[0]), L):
        for j in range(0, len(data[1]), L):
            if i < len(data[0]):
                if j < len(data[1]):             
                    temp.append((absolute(sum(data[i:i+L, j:j+L])))**q)
                    
    mu = sum(temp)
    return mu

L_array = [4,8,16,32,64]

log(asarray(eta_array[1]))

eta_array = sum_eta(asarray(derived_map), -800, L_array, 2, params)

plot(log(asarray(L_array)), log(asarray(eta_array)))

from scipy.stats import linregress
linregress(log(L_array), log(eta_array))[0]


D_2 = []
for i in energy_axis:
    eta_array = sum_eta(asarray(normalized_map), i, L_array, 2, params)
    D_2.append(linregress(log(L_array), log(eta_array))[0])  

np.save(''.join('D_2_{}.npy'.format(file)),D_2)

D_2_Array = np.load(''.join('D_2_{}.npy'.format(file)))

E_starting = -600
E_ending = -400

plot(energy_axis[index_of_energy(E_starting,params):index_of_energy(E_ending,params)], D_2[index_of_energy(E_starting,params):index_of_energy(E_ending,params)])

fig = figure()
fig.set_figheight(20)
fig.set_figwidth(40)
plot(energy_axis[index_of_energy(E_starting, params):index_of_energy(E_ending, params)], D_2[index_of_energy(E_starting, params):index_of_energy(E_ending, params)], linewidth=6, label = 'Box counting')
plot(E_array,D_2_array, linewidth=6, label = 'Cross correlation')
xlabel('E (meV)', fontsize=70)
ylabel('$D_2$', fontsize=70)
xticks(fontsize=40)
yticks(fontsize=40)
fig.savefig(''.join('D_2 cross and box'), format = 'jpg',bbox_inches='tight')
plt.legend(fontsize = 50)

D_2[index_of_energy(-500,params)]

linregress(energy_axis[index_of_energy(E_starting,params):index_of_energy(E_ending,params)],D_2 [index_of_energy(E_starting,params):index_of_energy(E_ending,params)])[0]

sum_box(normalized_map[index_of_energy(-500,params)], 16, 2)

# Plotting


def f_alpha_not_curve (data, E, L, params):
    index_E = []
    for i in E:
        index_E.append(index_of_energy(i, params))
    matrix_alpha_not = []
    matrix_f_alpha_not = []
    for i in index_E:
        f_alpha_q = asarray(f_q_alpha_q(substract_average(data[i])-np.min(substract_average(data[i])),0, L))
        matrix_alpha_not.append(f_alpha_q[1])
        matrix_f_alpha_not.append(f_alpha_q[0]) 
    return (matrix_alpha_not, matrix_f_alpha_not,index_E )

a =  f_alpha_not_curve (derived_map, energy_axis, 16, params)

np.save(''.join('f_alpha_not_curve_{}.npy'.format(file)),a)

a_array = np.load(''.join('f_alpha_not_curve_{}.npy'.format(file)))

mean = average(derived_map, (1,2))

fig, ax = plt.subplots(figsize = (15,5))
energy_axis_modif = energy_axis[index_of_energy(-950, params):index_of_energy(950, params)]
mean_modif = mean[index_of_energy(-950, params):index_of_energy(950, params)]
ax.plot(energy_axis_modif, mean_modif, color = 'black', linewidth=5)
ax.set_xlabel('E (mV)', fontsize=30)
ax.set_ylabel('dI/dV (pS)', fontsize=30,color = 'black')
xticks(arange(-950, 950, step=200), rotation= 45, fontsize=20)
yticks(fontsize=20)
fig.savefig(''.join('Average MAP'), format = 'jpg',bbox_inches='tight')


fig, ax = plt.subplots(figsize = (15,5))
energy_axis_modif = energy_axis[index_of_energy(-800, params):index_of_energy(900, params)]
mean_modif = mean[index_of_energy(-800, params):index_of_energy(900, params)]
ax.plot(energy_axis_modif, mean_modif, color = 'blue', linewidth=5)
ax.set_xlabel('E (meV)', fontsize=30)
ax.set_ylabel('dI/dV (pS)', fontsize=30,color = 'blue')
xticks(fontsize=15, rotation= 45)
yticks(fontsize=15,color = 'blue')
ax2=ax.twinx()
ax2.plot(E_1, gaussian_filter1d((float64(a[0])-2),3), color = 'red', linewidth=5)

ax2.set_ylabel('$\u03B1_{0}-2$',fontsize=30,x =2, color = 'red')
xticks(arange(-800, 900, step=400),fontsize=15)
yticks(fontsize=15, color = 'red')


fig.savefig(''.join('comparision'), format = 'jpg',bbox_inches='tight')



fig, ax = plt.subplots(figsize = (15,5))
energy_axis_modif = energy_axis[index_of_energy(-800, params):index_of_energy(900, params)]
mean_modif = mean[index_of_energy(-800, params):index_of_energy(900, params)]
ax.plot(energy_axis, D_2, color = 'blue', linewidth=5)
ax.set_xlabel('E (meV)', fontsize=30)
ax.set_ylabel('$D_2$', fontsize=30,color = 'blue')
xticks(fontsize=15, rotation= 45)
yticks(fontsize=15,color = 'blue')

ax2=ax.twinx()
ax2.plot(E_1,gaussian_filter1d((float64(a[0])-2),3), color = 'red', linewidth=5)
ax2.set_ylabel('$\u03B1_{0}-2$',fontsize=30,x =2, color = 'red')
xticks(arange(-800, 900, step=400),fontsize=15)
yticks(fontsize=15, color = 'red')


fig, ax = plt.subplots(figsize = (15,5))
energy_axis_modif = energy_axis[index_of_energy(-800, params):index_of_energy(-400, params)]
mean_modif = mean[index_of_energy(-800, params):index_of_energy(-400, params)]
D_2_modif = D_2[index_of_energy(-800, params):index_of_energy(-400, params)]
delta_2_modif = 2 - asarray(D_2_modif)
a_array_modif = a_array[0][index_of_energy(-800, params):index_of_energy(-400, params)]
ax.plot(energy_axis_modif, delta_2_modif, color = 'blue', linewidth=5)
plt.ylim(0, 0.1)
ax.set_xlabel('E (meV)', fontsize=30)
ax.set_ylabel('$\u0394_2$', fontsize=30,color = 'blue')
xticks(fontsize=15, rotation= 45)
yticks(fontsize=15,color = 'blue')

ax2=ax.twinx()
ax.plot(energy_axis_modif,gaussian_filter1d((2*float64(a_array_modif)-4),3), color = 'red', linewidth=5)
ax2.set_ylabel('$2(\u03B1_{0}-2)$',fontsize=30,x =2, color = 'red')
plt.ylim(0, 0.1)
yticks(fontsize=15, color = 'red')

fig.savefig(''.join('D_2 and Delta_2 curves of {}').format(file_name), format = 'jpg',bbox_inches='tight')



a_array_modif-2

from scipy.stats import linregress

linregress(L_points_log, f_q_points)[0]

interact(min(derived_map[i]), i = intslider(min=0, max = 255, step = 1))

map_1d = derived_map[index_of_energy(-500)].flatten()
count, bins, ignored = plt.hist(map_1d, bins=300) 
plot(bins[0:-1],count)

E_log = [-700,-600,-550,-500,-480,-450,-420]
count_1 = [[],[],[],[],[],[],[],[]]
bins_1 =  [[],[],[],[],[],[],[],[]]
bins_list = linspace(-0,10, 1000)
for i in range(0,len(E_log)):
    map_1d = ((derived_map[index_of_energy(E_log[i])])/average(derived_map[index_of_energy(E_log[i])])).flatten()
    count_1[i], bins_1[i], ignored = plt.hist(map_1d,bins= bins_list)

fig, ax = plt.subplots(figsize = (15,5))
for i in range(0,len(E_log)-1):
    ax.plot(bins_1[i][0:-1], gaussian_filter1d(count_1[i],10),  linewidth = 5, label=''.join('{} (meV)'.format(E_log[i])),)
    ax.legend(bbox_to_anchor =(1, 0.1), ncol = 1,fontsize=15)
    xlabel('dI/dV/<dI/dV>', fontsize=30)
    ylabel('Normalized count', fontsize=25)
    xticks(fontsize=15)
    yticks(fontsize=15)
    xlim(-0.0,3)
    
fig.savefig(''.join('Lognormal_Comparison'), format = 'jpg',bbox_inches='tight')

map_1d = derived_map[index_of_energy(-500)].flatten()
count, bins, ignored = plt.hist(map_1d, bins=300) 

s = lognormal(1.2, 0.02, 65536)*(max(map_1d)-min(map_1d))
count_1, bins_1, ignored = plt.hist(s, bins=300) 

from scipy.stats import *
lognorm.fit(bins[0:-1], loc=0)

count_1 = count_1/max(count_1)
count = count/max(count)
plot(bins_1[0:-1]-average(bins_1[0:-1]), count_1)
plot(bins[0:-1]-average(bins[0:-1])+1.15*(10**-13), count)

# Multifractal analysis of energy maps

L_points  = asarray([4], dtype=int)

E = [-410,-420,-430]

a_x = f_alpha_different_energies(derived_map,E,-1, 1, 99, L_points, 0.9)

f_alpha_different_energies(derived_map,[-900],-10, 10, 99, L_points, 0.9)

fig.savefig(''.join('multifractale analysis between {} and {} for {} data (meV), -400 for the transition'.format(E[0],E[-1],'Unnamed011')), format = 'jpg')

# Multi_fractal_analysis_different q

E = [-600, -550, -500, -480,-450, -420]

a_x_prin = [[],[],[],[],[],[]]

a_x_prin[0]=(f_alpha_different_energies(derived_map,[-600],-5, 5, 99, L_points, 0.9))

a_x_prin[1]=(f_alpha_different_energies(derived_map,[-550],-5, 5, 99, L_points, 0.9))

a_x_prin[2]=(f_alpha_different_energies(derived_map,[-500],-5, 5, 99, L_points, 0.9))

a_x_prin[3]=(f_alpha_different_energies(derived_map,[-470],-5, 5, 99, L_points, 0.9))

a_x_prin[4]=(f_alpha_different_energies(derived_map,[-450],-5, 5, 99, L_points, 0.9))

a_x_prin[5]=(f_alpha_different_energies(derived_map,[-420],-5, 5, 99, L_points, 0.9))

a_x_array = asarray(a_x_prin)
fig, ax = plt.subplots(figsize = (15,5))

for i in range(0,len(E)):
    ax.plot(a_x_array[i][0][1], (a_x_array[i][0][0]),  linewidth = 5, label=''.join('{} (meV)'.format(E[i])),)
    ax.legend(bbox_to_anchor =(1, 0.7), ncol = 1,fontsize=20)
    xlabel('\u03B1', fontsize=50)
    ylabel('$(f(\u03B1)-2){}*10^{5}$', fontsize=35)
    #xticks(arange(2, 2.09, step=0.01),fontsize=20)
    #yticks(arange(-5, 1, step=1),fontsize=20)
fig.savefig(''.join('Multifracral analysis'), format = 'jpg',bbox_inches='tight')

a_x_array = asarray(a_x_prin)
fig, ax = plt.subplots(figsize = (15,5))
ax.scatter(a_x_array[0][0][1], (a_x_array[0][0][0]),  linewidth = 5, label=''.join('{} (meV)'.format(E[i])),)
ax.scatter(a_x_array[4][0][1], (a_x_array[4][0][0]),  linewidth = 5, label=''.join('{} (meV)'.format(E[i])),)
ax.legend(bbox_to_anchor =(1, 0.7), ncol = 1,fontsize=20)
xlabel('\u03B1', fontsize=50)
ylabel('$(f(\u03B1)-2){}*10^{5}$', fontsize=35)


# Any map

def f_alpha(data, q_min, q_max, num_q, L_points):
    data_f = float64(data)
    q_points = asarray(linspace(-q_max, q_max, num=num_q))
    L_points_log = log(L_points)
    f = []
    alpha=[]
    for k in range (len (q_points)):
        points = f_q_alpha_q_points(data_f, q_points[k], L_points)
        f_q_points=[]
        alpha_points=[]
        for i in range(len(points)):
            f_q_points.append(points[i][0])
            alpha_points.append(points[i][1])
        f.append(linregress(L_points_log, f_q_points)[0])
        alpha.append(linregress(L_points_log, alpha_points)[0])
    return f , alpha

hop = float64(hopping_matrix(80,1,0.1,0.5,0,5,1))

a = f_alpha(hop,-0.1,0.1, 30, asarray([3,6,10,5]))

plot(a[1],a[0])



fig.savefig(''.join('multifractale analysis between {} and {} for {} data (meV)'.format(E[0],E[-1],'Unnamed011')), format = 'jpg')

# Hopping Matrix

def hopping_matrix(n,t,v,thresh,l,k, e):
    potential = v
    a  = zeros((n*n,n*n))
    b = potential*(np.random.rand(n*n))+e
    b[b < thresh]=0
    for i in range(n*n):
        a[i][i-1]=t
        a[i-n][i]=t
        a[i][i-n]=t
        if i<n*n-1:
            a[i][i+1]=t       
        if i < n*n-n:
            a[i+n][i]=t 
       
        if ((i/n).is_integer()==True):
            a[i][i-1]=0
            a[i-1][i]=0
            if i < n*n-n+1:
                a [i+n-1][i]=t
                a [i][i+n-1]=t
        a[i][i] = b[i] 
    #return a
    ander = eigenvector_matrix (a, n, l,k)
    fig = figure(figsize = (10,10))
    subplot(1,2,1)
    imshow(ander)
    
def hopping_matrix(n,t,v,thresh,e,k,l):
    potential = v
    a  = zeros((n*n,n*n))
    b = potential*(np.random.rand(n*n))+e
    b[b < thresh]=0
    for i in range(n*n):
        a[i][i-1]=t
        a[i-n][i]=t
        a[i][i-n]=t
        if i<n*n-1:
            a[i][i+1]=t       
        if i < n*n-n:
            a[i+n][i]=t 
       
        if ((i/n).is_integer()==True):
            a[i][i-1]=0
            a[i-1][i]=0
            if i < n*n-n+1:
                a [i+n-1][i]=t
                a [i][i+n-1]=t
        a[i][i] = b[i] 
    ander = eigenvector_matrix (a, n, l,k)
    fig = figure(figsize = (10,10))
    subplot(1,2,1)
    imshow(ander)
    return(ander)
    
def hopping_matrix_sum(n,t,v,thresh,e,k,l,n_sum):
    potential = v
    a  = zeros((n*n,n*n))
    b = potential*(np.random.rand(n*n))+e
    b[b < thresh]=0
    for i in range(n*n):
        a[i][i-1]=t
        a[i-n][i]=t
        a[i][i-n]=t
        if i<n*n-1:
            a[i][i+1]=t       
        if i < n*n-n:
            a[i+n][i]=t 
       
        if ((i/n).is_integer()==True):
            a[i][i-1]=0
            a[i-1][i]=0
            if i < n*n-n+1:
                a [i+n-1][i]=t
                a [i][i+n-1]=t
        a[i][i] = b[i] 
    #return a
    ander = eigenvector_matrix (a, n, l,k)
    for i in range (n_sum):
        ander_sum =+ eigenvector_matrix (a, n, i, k)
    fig = figure(figsize = (10,10))    
    subplot(1,2,1)
    imshow(ander)
    subplot(1,2,2)
    imshow(ander)
    
def eigenvector_matrix (ander, n, l, k):
    eigenvectors = eigs(ander, k)
    eigenvector_number = l
    eigenvectors_matrix = eigenvectors[1][:,eigenvector_number].reshape(n,n)
    return (abs(eigenvectors_matrix)**2)

"hopping_matrix(60,1,0.1,0.5,0,5,1)"

interact (hopping_matrix, n = IntSlider(min=4, max = 100, step = 1, continuous_update=False),  t = fixed(1) , v = FloatSlider(min=0, max = 1, step = 0.1, continuous_update=False), k = IntSlider(min=4, max = 20, step = 1, continuous_update=False), thresh = FloatSlider(min=0.2 ,max = 0.8, step = 0.05,ontinuous_update=False), l = IntSlider(min=0, max = 10, step = 1, continuous_update=False),e =IntSlider(min=0, max = 10, step = 1, continuous_update=False ))

interact (hopping_matrix_sum, n = IntSlider(min=4, max = 100, step = 1, continuous_update=False),  t = fixed(1) , v = FloatSlider(min=0, max = 1, step = 0.1, continuous_update=False), k = IntSlider(min=4, max = 20, step = 1, continuous_update=False), thresh = FloatSlider(min=0.2 ,max = 0.8, step = 0.05,ontinuous_update=False),l = IntSlider(min=0, max = 10, step = 1, continuous_update=False),e =IntSlider(min=0, max = 10, step = 1, continuous_update=False ), n_sum = IntSlider(min=0, max = 10, step = 1, continuous_update=False))

#grid_map

grid_spec
grid_matrix = []
for i in range(len(grid_spec)):
    for j in range(len(grid_spec[0])):
        grid_matrix.append(grid_spec[i,j])
a = grid_spec [0,0]
b = hankel(a)
U, s, Vh = linalg.svd(b)
plot(energy_axis, (U.transpose())[20])


def sing_value_decompose(data, x , y, n):
    #a = grid_spec [0,0]
    #b = hankel(a)
    U, s, Vh = linalg.svd(data)
    plot(energy_axis, (U.transpose())[n])

interact(sing_value_decompose, data = fixed(grid_matrix), x= IntSlider(min= 0, max = x_points, step=1,continuous_update = False ), y = IntSlider(min= 0, max = y_points, step=1,continuous_update = False), n= IntSlider(min= 1, max = 250, step=1,continuous_update = False ) )

mean(average(derived_map,(1,2)))

# Percolation

def thresholding(data, threshold, params, energy):
    data_e = data[index_of_energy(energy,params)] 
    thresholding = np.zeros((params['x_points'], params['y_points']))
    for i in range(params['x_points']):
        for j in range(params['y_points']):
            if  data_e [i][j] > threshold:
                 thresholding[i][j] = 1
            if  data_e [i][j] < threshold:
                 thresholding[i][j] = 0
    imshow(thresholding, interpolation = 'none')
    
def thresholding_energy(data, threshold, params):
    data_e = data
    thresholding = zeros((len(data), len(data[0])))
    for i in range(len(data[0])):
        for j in range(len(data[1])):
            if  data_e [i][j] > threshold:
                 thresholding[i][j] = 1
            if  data_e [i][j] < threshold:
                 thresholding[i][j] = 0
    return thresholding

def thresholding_mat(data, threshold, params, energy):
    data_e = data[index_of_energy(energy,params)] 
    thresholding = np.zeros((params['x_points'], params['y_points']))
    for i in range(params['x_points']):
        for j in range(params['y_points']):
            if  data_e [i][j] > threshold:
                 thresholding[i][j] = 1
            if  data_e [i][j] < threshold:
                 thresholding[i][j] = 0
    return thresholding      

a = thresholding (derived_map,0.0002, params,-500)

interact(thresholding, data = fixed(derived_map), threshold = FloatSlider(min = 0, max =0.01, step=0.001/100, continuous_update=False,readout_format = '.5s'), params = fixed(params), energy =  energy_widget)


def square_threshold_counting (data, points_cut, L, threshold, params, energy):
    data_e = data[index_of_energy(energy,params)] 
    data_e_cut = data_e[0: points_cut-1,0: points_cut-1]
    data_thresh = thresholding_energy(data_e_cut, threshold, params)
    matrix = np.zeros((points_cut, points_cut))
    for i in range (0, points_cut, L):
        for j in range(0, points_cut, L):
            if sum(data_thresh[i:i+L,j:j+L]) > 0:
                matrix [i:i+L,j:j+L] = 1
            if sum(data_thresh[i:i+L,j:j+L]) == 0:
                matrix [i:i+L,j:j+L] = 0
    imshow(matrix, interpolation = 'none')


L_square = [2,3,6,8,9,12,16]

interact (square_threshold_counting, data = fixed(derived_map), points_cut = fixed(216), L = IntSlider(min = 2, max= 24, steps = 1, continuous_update = False) , threshold = FloatSlider(min = 0, max =0.01, step=0.01/100, continuous_update=False,readout_format = '.5s'), params = fixed(params), energy =  energy_widget)


def square_threshold_counting_matrix (data, points_cut, L, threshold, params, energy):
    data_e = data[index_of_energy(energy,params)] 
    data_e_cut = data_e[0: points_cut-1,0: points_cut-1]
    data_thresh = thresholding_energy(data_e_cut, threshold, params)
    matrix = np.zeros((points_cut, points_cut))
    for i in range (0, points_cut, L):
        for j in range(0, points_cut, L):
            if sum(data_thresh[i:i+L,j:j+L]) > 0:
                matrix [i:i+L,j:j+L] = 1
            if sum(data_thresh[i:i+L,j:j+L]) == 0:
                matrix [i:i+L,j:j+L] = 0
    return matrix
def square_threshold_different_size (data, points_cut, L_array, threshold, params, energy):
    square_array = []
    L_array_reverse = []
    for i in L_array:
        square_array.append(np.sum(square_threshold_counting_matrix (data, points_cut, i, threshold, params, energy)/(i*i)))
        L_array_reverse .append(1/i)
    return square_array, L_array_reverse

a = square_threshold_different_size (derived_map, 216, L_square, 0.0018, params, -502)

plot(log(a[1]),log(a[0]))

linregress(log(a[1]),log(a[0]))[0]

a = square_threshold_different_size (derived_map, 216, L_square, 0.0030, params, -550)

plot(log(a[1]),log(a[0]))

linregress(log(a[1]),log(a[0]))

plot([-850, -750, -600, -555, -552,-502,-485,-465,-435,-400], [6.13,5.73,4.33,3.2,3.1,1.8,1.4,0.97,0.55, 0.00035])

linregress([-600, -555, -552,-502,-485,-465,-435], [4.33,3.2,3.1,1.8,1.4,0.96,0.55])

def find_D (data, points_cut, L_array, params, energy, init_value, fractal):
        init_matrix = square_threshold_different_size (data, points_cut, L_array, init_value, params, energy)
        D_2_init = linregress(log(init_matrix[1]),log(init_matrix[0]))
        
        while absolute (D_2_init[0] - fractal) > 0.001:
            if D_2_init[0] - fractal > 0:
                init_value = init_value + absolute (D_2_init[0] - fractal)/2000
            if D_2_init[0] - fractal < 0:
                init_value = init_value - absolute (D_2_init[0] - fractal)/2000
            init_matrix = square_threshold_different_size (data, points_cut, L_array, init_value, params, energy)
            D_2_init = linregress(log(init_matrix[1]),log(init_matrix[0]))
            print(D_2_init[0])
        return init_value

def find_D_different_energies (data, points_cut, L_array, params, E_array, init_value, fractal):
    D_2_array = []
    for i in E_array:
        D_2_array.append(find_D (data, points_cut, L_array, params, i, init_value,fractal)) 
    return D_2_array, E_array
    
    

E_array = linspace(-450, -600, 50)

array_of_D_energies = find_D_different_energies (derived_map , 216, L_square, params, E_array, 0.001)

plot( array_of_D_energies[1],array_of_D_energies[0],)

linregress(array_of_D_energies[0],array_of_D_energies[1])

j = []
for i in range(len(array_of_D_energies[1])):
    j.append(thresholding_mat(derived_map, array_of_D_energies[0][i], params, array_of_D_energies[1][i]))



def show(data, energy):
    i = int((energy-(-450))*((50-1)/(-600-(-450))))
    imshow(data[i], interpolation = 'none')

interact(show, data = fixed(j), energy = FloatSlider(min= -600, max = -450, step=1, continuous_update = False))


imshow(asarray(j)[:,4,:], interpolation = 'none')

def show_energy(data, i):
    imshow(asarray(data)[:,i,:], interpolation = 'none')

interact(show_energy, data = fixed(j), i = IntSlider(min = 0, max = 255, step = 1))

E_array_2 = linspace(-400, -800, 250)

array_of_D_energies_2 = find_D_different_energies (derived_map , 216, L_square, params, E_array_2, 0.001)



def show_2(data, energy):
    i = int((energy-(-400))*((250-1)/(-800-(-400))))
    imshow(data[i], interpolation = 'none')

j_2 = []
for i in range(len(array_of_D_energies_2[1])):
    j_2.append(thresholding_mat(derived_map, array_of_D_energies_2[0][i], params, array_of_D_energies_2[1][i]))

interact(show_2, data = fixed(j_2), energy = FloatSlider(min= -800, max = -400, step=1, continuous_update = False))

def show_energy_2(data, i):
    imshow(asarray(data)[:,i,:], interpolation = 'none', extent=[0,150,-800,-400])

interact(show_energy_2, data = fixed(j_2), i = IntSlider(min = 0, max = 255, step = 1))

plot( array_of_D_energies_2[1],array_of_D_energies_2[0])

summer = np.sum(j,axis = 0)

imshow(summer)

import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.ticker import LinearLocator
import numpy as np

fig, ax = plt.subplots(subplot_kw={"projection": "3d"})

# Make data.
X = np.arange(0, 256, 1)
Y = np.arange(0, 256, 1)
X, Y = np.meshgrid(X, Y)
Z = summer

# Plot the surface.
surf = ax.plot_surface(X, Y, Z, cmap=cm.coolwarm,
                       linewidth=0, antialiased=False)

# Customize the z axis.
ax.set_zlim(0, 50)
ax.zaxis.set_major_locator(LinearLocator(10))
# A StrMethodFormatter is used automatically
ax.zaxis.set_major_formatter('{x:.02f}')

# Add a color bar which maps values to colors.
fig.add_collection3d(surf, shrink=0.5, aspect=5)

plt.show()

np.max(summer)

# Monte Carlo

random_ising = np.random.randint(2, size=(256,256))

def ising_iteration (data, params, betha, iteration):
    zero_one = data
    random_prob = np.random.random((256, 256))
    for k in range(iteration):
        X = zeros((256,256))
        for i in range ((params['x_points'])):
            for j in range ((params['y_points'])):
                w_up = exp(betha*(zero_one[i-2][j-1]+zero_one[i][j-1]+zero_one[i-1][j-2]+zero_one[i-1][j]))
                w_down = exp(-betha*(zero_one[i-2][j-1]+zero_one[i][j-1]+zero_one[i-1][j-2]+zero_one[i-1][j]))
                P_up = w_up / (w_up + w_down)
                if P_up < random_prob[i-1][j-1]:
                    X [i-1][j-1] = 1
                else:
                    X [i-1][j-1] = -1
        zero_one = X
    return X


a = ising_iteration (random, params, 0.001, 5)

imshow (a, interpolation = 'none') 

b = ising_iteration (random, params, 100, 10)

imshow (a, interpolation = 'none') 

imshow(random,interpolation = 'none')

average(a)

average(a)



array_of_D_energies_1_74 = find_D_different_energies (derived_map , 216, L_square, params, E_array_2, 0.001, 1.74)

j_1_74 = []
for i in range(len(array_of_D_energies_1_74[1])):
    j_1_74.append(thresholding_mat(derived_map, array_of_D_energies_2[0][i], params, array_of_D_energies_2[1][i]))

array_of_D_energies_1_73 = find_D_different_energies (derived_map , 216, L_square, params, E_array_2, 0.001, 1.73)

j_1_73 = []
for i in range(len(array_of_D_energies_1_73[1])):
    j_1_73.append(thresholding_mat(derived_map, array_of_D_energies_2[0][i], params, array_of_D_energies_2[1][i]))

array_of_D_energies_1_77 = find_D_different_energies (derived_map , 216, L_square, params, E_array_2, 0.001, 1.77)

j_1_77 = []
for i in range(len(array_of_D_energies_1_77[1])):
    j_1_77.append(thresholding_mat(derived_map, array_of_D_energies_2[0][i], params, array_of_D_energies_2[1][i]))

array_of_D_energies_1_76 = find_D_different_energies (derived_map , 216, L_square, params, E_array_2, 0.001, 1.76)

j_1_76 = []
for i in range(len(array_of_D_energies_1_76[1])):
    j_1_76.append(thresholding_mat(derived_map, array_of_D_energies_2[0][i], params, array_of_D_energies_2[1][i]))

# Multifractal functions

def sum_box(data, L, q):
    temp = []
    for i in range (0, len(data[0]), L):
        for j in range(0, len(data[1]), L):
            if i < len(data[0]):
                if j < len(data[1]):             
                    temp.append((absolute(sum(data[i:i+L, j:j+L])))**q)
    mu = sum(temp)
    return mu, temp

def l_q_differentiate (data,q,L):
    khi_q = []
    khi_q_sum =[]
    eta_q = []
    #for i in range(1,10,1):
    for j in range(-q, q, 1):
        khi_q.append(sum_box(data, L, j))
    eta_q.append(log(khi_q)/log(L))
    return eta_q

def different_L_array_least_square (data, q, L_max):
    L_array=[]
    etq_q_array=[]
    for i in range(2, L_max):
        L_array.append(i)
        etq_q_array.append(l_q_differentiate (data, q, i))
    return L_array,  etq_q_array

def log_khi_log_L(data,q,L):
    L_array=[]
    khi_q = []
    L_array_log=[]
    khi_q_log = []
    tau_q_array=[]
    q_array = []
    for j in range (-q, q):
        for i in range (2, L):
            L_array.append(i) 
            khi_q.append(sum_box(data, i, j))
            L_array_log.append(log(i))
            khi_q_log.append(log(sum_box(data, i, j)))
        tau_q = scipy.stats.linregress(L_array_log, khi_q_log)[0]
        tau_q_array.append(tau_q)
        q_array.append(j)
    khi_correct = asarray(tau_q_array) * -1 
    return q_array, khi_correct


def least_square (x_data,y_data):
    def func(params, xdata, ydata):
        return (ydata - numpy.dot(xdata, params))

def f_q_alpha_q (data_f, q, L):
    data = float64(data_f)
    mu_matrix = sum_box(data, L, q)
    mu_matrix_q_1 = sum_box(data, L, 1)
    f_q = sum (((mu_matrix[1]/mu_matrix[0])*log(mu_matrix[1]/mu_matrix[0]))/log(L/256))
    alpha_q = sum (((mu_matrix[1]/mu_matrix[0])*log(mu_matrix_q_1[1]/mu_matrix_q_1[0]))/log(L/256))
    return f_q, alpha_q

def f_q_alpha_q_points (data, q, L_points):
    f_q_alpha_q_points =[]
    for i in L_points:
        f_q_alpha_q_points.append(f_q_alpha_q(data, q, i))
    return f_q_alpha_q_points

def f_alpha(data, q_min, q_max, num_q, L):
    
    q_points = asarray(linspace(-q_max, q_max, num=num_q))
    f = []
    alpha=[]
    for k in range (len (q_points)):
        points = f_q_alpha_q_points(data, q_points[k], L)
        f_q_points=[]
        alpha_points=[]
        for i in range(len(points)):
            f_q_points.append(points[i][0])
            alpha_points.append(points[i][1])
        f.append( f_q_points[0])
        alpha.append(alpha_points[0])
    return f , alpha

def f_alpha_different_energies(data, E, q_min, q_max, num_q, L_points,m_tresh):
    data = absolute(float64(data))
    index_E = asarray(list(map(index_of_energy,E)))
    matrix = []
    for i in index_E:
        a = f_alpha(data[i],-q_min, q_max,num_q, L_points)
        matrix.append(a)
    return (matrix)

