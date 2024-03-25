# # Main functions


def index_of_energy (energy, params):
    spec_N = int((abs(energy-params['bias_start']))*((params['num_points']-1)/abs((params['bias_end']-params['bias_start']))))
    return spec_N


def energy_of_index (index, params):
    energy = (index/((params['num_points']-1)/(params['bias_end']-params['bias_start'])))+params['bias_start']
    return energy


def x_i (x, params):
    x_index = int((x)*((params['x_points']-1)/(params['len_x'])))
    return x_index


def y_i (y, params):
    y_index = params['y_points']-int((y)*((params['y_points']-1)/(params['len_y'])))-1
    return y_index


def show_spectra(data, x, y,params):
    energy_axis = linspace(params['bias_start'],params['bias_end'],params['num_points'])                                
    plot(energy_axis, data[x][y])
    
    
def Gaussian_filtering (data, sigma):
    filtered_data = gaussian_filter1d(data, sigma)
    return filtered_data
    
    
def Gaussian_filtering_plot (data, x, y, sigma, params):
    energy_axis = linspace(params['bias_start'],params['bias_end'],params['num_points'])    
    x_1 = x_i(x, params)
    y_1 = y_i(y, params)
    data_1 = (data[y_1][x_1])
    filtered_data = gaussian_filter1d(data_1 ,sigma)
    figure(figsize = (15,5))
    
    subplot(1,2,1)
    plot(energy_axis, data_1)
    plot(energy_axis, filtered_data)
    figure(figsize = (15,15))
    
    subplot(1,2,2)
    imshow(params['slope_corrected_topography'], extent=[0,params['len_x'],0,params['len_y']]) 
    scatter(x,y,color='r')
    global sigma_global
    sigma_global = sigma

    
def all_gaussian (data, sigma): 
    filtered = apply_along_axis(Gaussian_filtering, 2 ,data, sigma)
    return filtered


def derivative_1d (data,params):
    energy_axis = linspace(params['bias_start'],params['bias_end'],params['num_points'])                                  
    matrix_1 = insert((delete(data, -1)), 0, 2*data[0]-data[1])
    matrix_2 = insert((delete(data,0)),-1, 2*data[-1]-data[-2])
    diff = matrix_2-matrix_1
    derive = diff/(energy_axis[2]-energy_axis[1])
    return derive


def all_derive (data,params):
    derive = asarray(apply_along_axis(derivative_1d, 2, data,params))
    return derive


def make_map (data, params):
    grid_spec_map = list(map(lambda x: (data[:,:,x]), arange(params['num_points'])))
    return (grid_spec_map)


def make_spec (data, params):
    q = []
    for i in range(len(data[0][0])):
        for j in range(len(data[0])):
            q.append(data[:,j,i])
    return asarray(q).reshape(len(data[0][0]),len(data[0]),params['num_points'])
     
    
def map_in_energy_new (data, energy, size, _range, n, m_tresh):
    fig=figure(figsize = (size,size))
    e_plus = index_of_energy(energy+_range,params)
    e_minus = index_of_energy(energy-_range,params)
    average = average_data(data, e_minus, e_plus)
    median_filtered = median_filtering (data[index_of_energy(energy,params)], n, m_tresh)
    median_filtered_average = median_filtering(average, n, m_tresh)
    
    subplot(2,2,1)
    imshow(data[index_of_energy(energy,params)])
    xlabel('X', fontsize=16)
    ylabel('Y', fontsize=16)
    title('energy map after gaussian filtering', fontsize=16) 
    
    subplot(2,2,2)
    imshow(median_filtered)
    xlabel('X', fontsize=16)
    ylabel('Y', fontsize=16)
    title('energy map after gaussian filtering with median filtering', fontsize=16) 
    
    subplot(2,2,3)
    imshow(average)
    xlabel('X', fontsize=16)
    ylabel('Y', fontsize=16)
    title('energy map for a given energy window', fontsize=16) 
    
    subplot(2,2,4)
    imshow(median_filtered_average)
    xlabel('X', fontsize=16)
    ylabel('Y', fontsize=16)
    title('energy map for a given energy window with median filtering', fontsize=16) 
    
    global energy_glob
    energy_glob = energy
    global fig_glob 
    fg_glob = fig
    
    
def substract_average (data):
    average = apply_along_axis(np.average, 1, data)
    average = transpose(average[newaxis])
    data_averaged = data-average
    return data_averaged


def substract_average_map (data):
    substracted = []
    for i in range(len(data)):
        substracted .append(substract_average(data[i]))
    return substracted


def points_coordinate (x,y,params):
    x1 = [x, x]
    y1 = [0, params('num_points')]
    x2 = [0, params('num_points')]
    y2 = [y, y]
    plot(x1, y1, x2, y2, color="black", linewidth=1)
    
    
def circle (data, x , y, r):
    figure(figsize = (10,10))
    theta = np.linspace(0, 2*np.pi, 100)
    radius = r
    a =radius*np.cos(theta)
    b = radius*np.sin(theta)
    figure(figsize = (15,15))
    subplot(1,2,1)
    plot(a+x,b+y,color="black", linewidth=1)
    imshow(topo)
    figure(figsize = (20,5))
    subplot(1,2,2)  
    plot(energy_axis, mean_spectrum_circ(data,x,y,r,len(data[0]), len(data)))

    
def mean_spectrum_circ (data,centerx,centery,radius,maxspecx,maxspecy):
    specpts = len(data) # number of points in each spectrum
    xpoints, ypoints = maxspecx, maxspecy # number of spectra in x and y directions
    output = []
    for i in range(specpts):
        output.append(0)
        number_of_spectra = 0
        for j in range(max(centerx-radius, 0), min(centerx+radius, xpoints-1)):
            for k in range(max(centery-radius, 0), min(centery+radius, ypoints-1)):
                if (j-centerx)**2+(k-centery)**2 < radius**2:
                    output[i] += data[i][k][j]
                    number_of_spectra += 1
        output[i] = output[i]/float(number_of_spectra)
    return output    


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
    filtered_mean = filtered_mean
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
    plot(energy_axis[min(lower_index,upper_index):max(lower_index,upper_index)],mean[min(lower_index,upper_index):max(lower_index,upper_index)], linewidth=6)
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
    
    
def map_in_energy_derived_with_average_spectra(data, energy, figure_size, m_tresh,lower_cut,upper_cut,bar_lower_percentage,bar_higher_percentage,params): 
    energy_axis = linspace(params['bias_start'],params['bias_end'],params['num_points'])                                    
    e = index_of_energy (energy,params)
    fig = figure(figsize = (figure_size,figure_size))
    median = median_filtering(data[e], 3, m_tresh)
    mean = average(data, (1,2))

    upper_index = index_of_energy(upper_cut,params)
    lower_index = index_of_energy(lower_cut,params)
    fig = figure()
    fig.set_figheight(figure_size)
    fig.set_figwidth(figure_size*3)
    fig.subplots_adjust(hspace=0.4, wspace=0.4) 
    
    subplot(1,2,1)
    im = imshow(median, extent=[0,params['len_x'],0,params['len_y']])
    xlabel('X (nm)', fontsize=50,y=1)
    ylabel('Y (nm)', fontsize=50)
    title(''.join('derived energy map at {} meV'.format(round(energy))), fontsize=50, y= 1.07)
    xticks(fontsize=50)
    yticks(fontsize=50)
    cbar =colorbar(im, boundaries=linspace(np.min(median)*bar_lower_percentage/100,np.max(median)*(100-bar_higher_percentage)/100,1000),shrink=1,ticks=linspace(np.min(median)*bar_lower_percentage/100,np.max(median)*(100-bar_higher_percentage)/100,10))
    cbar.minorticks_on()
    cbar.ax.tick_params(labelsize=50)
    
    subplot(1,2,2)
    plot(energy_axis[min(lower_index,upper_index):max(lower_index,upper_index)],mean[min(lower_index,upper_index):max(lower_index,upper_index)], linewidth=6)
    plt.plot(energy_axis[e], mean[e], 'ro', markersize=25)
    xlabel('E (meV)', fontsize=50)
    ylabel('dI/dV (pS)', fontsize=50)
    title(''.join('Average Spectra of the map'.format(energy)), fontsize=50, y= 1.07)
    xticks(fontsize=50, rotation=45)
    yticks(fontsize=50)

    global fig_glob 
    fig_glob = fig
    global energy_glob
    energy_glob = energy 
    
    
def average_spectra (data, figure_size, lower_cut,upper_cut,sigma,params):
    energy_axis = linspace(params['bias_start'],params['bias_end'],params['num_points'])   
    mean = average(data, (1,2))
    filtered_mean = gaussian_filter1d(mean, sigma)
    upper_index = index_of_energy(upper_cut,params)
    lower_index = index_of_energy(lower_cut,params)
    fig = figure()
    fig.set_figheight(figure_size)
    fig.set_figwidth(figure_size*3)
    fig.subplots_adjust(hspace=0.4, wspace=0.4) 
    fig = figure()
    fig.set_figheight(figure_size)
    fig.set_figwidth(figure_size)
    plot(energy_axis[lower_index:upper_index],filtered_mean[lower_index:upper_index], linewidth=6)
    xlabel('E (meV)', fontsize=50)
    ylabel('dI/dV (pS)', fontsize=50)
    xticks(fontsize=50, rotation=45)
    yticks(fontsize=50)
    margins(x=0)
    global fig_glob 
    fig_glob = fig
    
def average_spectra_data (data, figure_size, lower_cut,upper_cut,sigma, params):
    energy_axis = linspace(params['bias_start'],params['bias_end'],params['num_points'])   
    mean = average(data, (1,2))
    filtered_mean = gaussian_filter1d(mean, sigma)
    upper_index = index_of_energy(upper_cut,params)
    lower_index = index_of_energy(lower_cut,params)
    return energy_axis[lower_index:upper_index], filtered_mean[lower_index:upper_index]
    
def average_spectra_2maps (data, data_2, figure_size, lower_cut, upper_cut, sigma, params):
    energy_axis = linspace(params['bias_start'],params['bias_end'],params['num_points'])   
    mean = average(data, (1,2))
    mean_2 = average(data_2, (1,2))
    filtered_mean = gaussian_filter1d(mean, sigma)
    filtered_mean_2 = gaussian_filter1d(mean_2, sigma)
    upper_index = index_of_energy(upper_cut,params)
    lower_index = index_of_energy(lower_cut,params)
    fig = figure()
    fig.set_figheight(figure_size)
    fig.set_figwidth(figure_size*3)
    fig.subplots_adjust(hspace=0.4, wspace=0.4) 
    fig = figure()
    fig.set_figheight(figure_size)
    fig.set_figwidth(figure_size) 
    plot(energy_axis[lower_index:upper_index],filtered_mean[lower_index:upper_index], linewidth=6, label = 'At 2 K')
    plot(energy_axis[lower_index:upper_index],filtered_mean_2[lower_index:upper_index], linewidth=6, label = 'At 300 mK')
    legend(loc='best', frameon=False,prop={'size': 40})
    xlabel('E (meV)', fontsize=50)
    ylabel('dI/dV (pS)', fontsize=50)
    xticks(fontsize=50)
    yticks(fontsize=50)
    margins(x=0)
    global fig_glob 
    fig_glob = fig
    
    
def map_in_energy_derived (data, energy, figure_size, m_tresh,bar_lower_percentage,bar_higher_percentage,params):  
    e = index_of_energy (energy,params)
    fig = figure(figsize = (figure_size,figure_size))
    median = median_filtering (data[e], 3, m_tresh)
    mean = average(data, (1,2))
    mean = mean*(10**12)
    fig = figure()
    fig.set_figheight(figure_size)
    fig.set_figwidth(figure_size*3)
    fig.subplots_adjust(hspace=0.4, wspace=0.4) 
    im = imshow(median, extent=[0,params['len_x'],0,params['len_y']])
    xlabel('X (nm)', fontsize=50,y=1)
    ylabel('Y (nm)', fontsize=50)
    title(''.join('derived energy map at {} meV'.format((energy))), fontsize=50, y= 1.07)
    xticks(fontsize=50)
    yticks(fontsize=50)
    cbar =colorbar(im, boundaries=linspace(np.min(median)*bar_lower_percentage/100,np.max(median)*(100-bar_higher_percentage)/100,1000),shrink=1,ticks=linspace(np.min(median)*bar_lower_percentage/100,np.max(median)*(100-bar_higher_percentage)/100,10))
    cbar.minorticks_on()
    cbar.ax.tick_params(labelsize=50)
    global fig_glob 
    fig_glob = fig
    global energy_glob
    energy_glob = energy

    
def mesh (data):
    a = []
    b = []
    for i in range(len(data[0])):
        for j in range( len(data)):
            a.append([j,i])
            b.append(data[j][i])
    return a , b


def substract_polynomyal_slope (data, order, size, len_x, len_y, params):
    X = (mesh(data))[0]
    vector = (mesh(data))[1]
    poly = PolynomialFeatures(degree=order)
    X_ = poly.fit_transform(X)
    clf = linear_model.LinearRegression()
    clf.fit(X_, vector)
    predict_ = poly.fit_transform (X)
    pre = clf.predict(predict_)
    vector = asarray(vector)
    fitted_data = (vector - pre).reshape((len(data[0]), len(data)))
    fig = figure(figsize = (size,size))
    imshow (transpose(fitted_data), extent=[0,params['len_x'],0,params['len_y']])


def average_data(data, E_minus, E_plus):
    data_cut = asarray(data) [E_minus:E_plus,:,:]
    data_averaged = average(data_cut, axis=0)
    return data_averaged


def average_(data, percentage, num_points):
    percent = percentage/100
    average_ = mean(asarray(data[0:int(percent*num_points)]))
    new_data = (data/average_)*-1
    return new_data
    
    
def all_Normalize (data, percentage, num_points, params): 
    Normalized = apply_along_axis(average_, 2 ,data, percentage, params['num_points'])
    return Normalized    


def plotting(data, x, y,params):
    x_1 = x_i(x)
    y_1 = y_i(y)
    data_1 = (data[x_1][y_1])
    figure(figsize = (15,5))
    subplot(1,2,1)
    plt.plot(energy_axis, data_1)
    figure(figsize = (15,15))
    subplot(1,2,2)
    imshow(slope_corrected_topography, extent=[0,params['len_x'], 0, params['len_y']]) 
    scatter(x,y,color='r')
    
    
def normalization_check (data, percentage, lower_cut, upper_cut, sigma, figure_size, sigma_global):
    normalized_spec  = all_Normalize(data, percentage)
    Filtered_data_normalized = all_gaussian(normalized_spec , sigma_global)
    Filtered_map_normalized = make_map (Filtered_data_normalized, num_points)
    derived_data_normalized = all_derive(Filtered_data_normalized)
    derived_map_normalized = make_map (derived_data_normalized, num_points) 
    average_spectra (derived_map_normalized, figure_size, lower_cut, upper_cut, sigma)    
    
    
def mean_spectrum_circ(data,centerx,centery,radius,maxspecx,maxspecy):
    specpts = len(data) # number of points in each spectrum
    xpoints, ypoints = maxspecx, maxspecy # number of spectra in x and y directions
    output = []
    for i in range(specpts):
        output.append(0)
        number_of_spectra = 0
        for j in range(max(centerx-radius, 0), min(centerx+radius, xpoints-1)):
            for k in range(max(centery-radius, 0), min(centery+radius, ypoints-1)):
                if (j-centerx)**2+(k-centery)**2 < radius**2:
                    output[i] += data[i][k][j]
                    number_of_spectra += 1
        output[i] = output[i]/float(number_of_spectra)
    return output
#Selects and sum the spectra in a circular zone of radius "size" around the points
#centerx and centery on set of data "data" with a threshold "threshold" for the noisy spectra around an average conductance "avg_conductance"    


def compress (data, compress):
    array_data = array(data)
    sigma = std(array_data)
    percentage = compress
    compressed_data = log((array_data /sigma*100)/((exp(compress)-0.9))+1)
    return compressed_data


# In[ ]:




