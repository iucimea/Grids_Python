#!/usr/bin/env python
# coding: utf-8

# # Median filtering

# In[ ]:


def median_filtering(data,n,m_tresh): #circular_bc
    if(n==1):return data
    temp = []
    indexer = n // 2  #indexer is used to shift the kernel from 0; filter_size to -indexer; +indexer
    data_final = []
    data_final = zeros((len(data),len(data[0])))
    for i in range(len(data)):
        for j in range(len(data[0])):
            for z in range (n):
                for w in range (n):
                    x=i+z-indexer
                    y=j+w-indexer
                    temp.append(data[x%len(data),y%len(data[0])])#circular bc
            temp.sort()
            data_final[i][j] = temp[len(temp) // 2]
            temp = []
    G = data_final
    I = data
    m_tresh*=std(data)
    return G-m_tresh*tanh((G-I)/m_tresh)


def median_filtering_2 (data,n,m_tresh):
    if(n==1):return data
    temp = []
    indexer = n // 2  #indexer is used to shift the kernel from 0;filter_size to -indexer;+indexer
    data_final = []
    data_final = zeros((len(data),len(data[0])))
    for i in range(len(data)):
        for j in range(len(data[0])):
            for z in range (n):
                for w in range (n):
                    x=i+z-indexer
                    y=j+w-indexer
                    temp.append(data[x%len(data),y%len(data[0])])#circular bc
            temp.sort()
            data_final[i][j] = temp[len(temp) // 2]
            temp = []
    G = data_final
    I = data
    m_tresh*=std(data)
    return G  
   

def median_filtering_no (data, n, m_tresh, x_points, y_points):
    data_1 = absolute(data)
    image = zeros(x_points*y_points).reshape(x_points,y_points)
    for i in range(0, x_points, n):
        for j in range (0, y_points, n):
            stds = np.std(data_1[i:i+n][j:j+n])
            med = median(data_1[i:i+n][j:j+n])
            filtered = med+(m_tresh*stds)*tanh((abs(data_1[i:i+n][j:j+n]-med))/(m_tresh*stds))
            image[i:i+n][j:j+n]=filtered
    return image


    
def median_filtering_all_maps(data, n, m_tresh):
    filtered_maps = []
    for i in range(len(data)):
        filtered_maps.append(median_filtering(data[i],n, m_tresh))
    return asarray(filtered_maps)


def median_filtering_all_maps_2(data, n, m_tresh):
    filtered_maps = []
    for i in range(len(data)):
        filtered_maps.append(median_filtering_2(data[i],n, m_tresh))
    return asarray(filtered_maps)


def map_in_energy_spectra_mean_range (data, energy, figure_size, x, y, r, _range, m_tresh):  
    mean = mean_spectrum_circ(data,x,y,r,len(data[0]), len(data))
    e = index_of_energy(energy,params)
    e_plus = index_of_energy(energy+_range,params)
    e_minus = index_of_energy(energy-_range,params)
    average = average_data(data, e_minus, e_plus)
    fig=figure(figsize = (figure_size,figure_size))
    theta = np.linspace(0, 2*np.pi, 100)
    radius = r
    a = radius*np.cos(theta)
    b = radius*np.sin(theta)
    x1 = [energy_axis[e], energy_axis[e]]
    y1 = [min(mean), max(mean)]
    x2 = [min(energy_axis), max(energy_axis)]
    y2 = [mean[e], mean[e]]
    fig = figure(figsize = (figure_size,figure_size))
    subplot(1,3,1)
    imshow(data[e])
    plot(a+x,b+y,color="black", linewidth=1)
    
    subplot(1,3,2)
    imshow(average)
    
    subplot(1,3,3)
    imshow (median_filtering(average,3, m_tresh))
    
    subplot(1,2,1)
    plot(energy_axis[25:-25],mean[25:-25])
    plot(x1, y1, x2, y2, color="black", linewidth=1)
    
    global fig_glob 
    fig_glob = fig
    global energy_glob
    energy_glob = energy


def map_in_energy_spectra_mean (data, energy, figure_size, x, y, r, m_tresh,sigma,lower_cut,upper_cut,bar_lower_percentage,bar_higher_percentage,params):
    energy_axis = linspace(params['bias_start'],params['bias_end'],params['num_points'])                             
    x_1 = x_i(x,params)
    y_1 = y_i(y,params)
    r_1 = x_i(r,params)
    mean = mean_spectrum_circ(data,x_1,y_1,r_1, params['x_points'],params['y_points'],)
    filtered_mean = gaussian_filter1d(mean, sigma)
    e = index_of_energy(energy,params)
    fig = figure(figsize = (figure_size,figure_size))
    theta = np.linspace(0, 2*np.pi, 100)
    radius = r
    a = radius*np.cos(theta)
    b = radius*np.sin(theta)
    median = median_filtering (data[e], 3, m_tresh)
    fig = figure()
    fig.subplots_adjust(hspace=0.4, wspace=0.4) 
    fig.set_figheight(figure_size)
    fig.set_figwidth(figure_size*3)
    upper_index = index_of_energy(upper_cut,params)
    lower_index = index_of_energy(lower_cut,params)
    
    subplot(1,2,1)
    im = imshow(median, extent=[0,params['len_x'],0,params['len_y']])
    plot(a+x,b+y,color="red", linewidth=5)
    xlabel('X (nm)', fontsize=50, y=1.01)
    ylabel('Y (nm)', fontsize=50, x = 1.01)
    title(''.join('derived energy map at {} meV'.format(energy)), fontsize=50, y= 1.1)
    xticks(fontsize=50)
    yticks(fontsize=50)
    cbar =colorbar(im, boundaries=linspace(np.min(median)*bar_lower_percentage/100,np.max(median)*(100-bar_higher_percentage)/100,1000),shrink=1,ticks=linspace(np.min(median)*bar_lower_percentage/100,np.max(median)*(100-bar_higher_percentage)/100,10))
    cbar.minorticks_on()
    cbar.ax.tick_params(labelsize=50)
    
    subplot(1,2,2)
    plot(energy_axis[min(lower_index,upper_index):max(lower_index,upper_index)],filtered_mean[min(lower_index,upper_index):max(lower_index,upper_index)], linewidth=6)
    plot(energy_axis[e], filtered_mean[e], 'ro', markersize=25)
    xlabel('E (meV)', fontsize=50)
    ylabel('dI/dV (pS)', fontsize=50)
    title(''.join('Average Spectra of the map at the circle'.format(energy)), fontsize=50, y= 1.07)
    xticks(fontsize=50, rotation = 45)
    yticks(fontsize=50)
    global fig_glob 
    fig_glob = fig
    global energy_glob
    energy_glob = energy
    global circle_spectra
    circle_spectra = filtered_mean
    global x_glob
    x_glob = x
    global y_glob
    y_glob = y
    global r_glob
    r_glob = r  

