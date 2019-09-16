"""
Dysarthric speech processing toolbox
Utilities

@functions:
    - read_wav_file(path)
    - make_directory(new_directory)
    - save_data(data, file_name)
    - load_data(file_name)
    - find_peaks(data, order, interpolate)
    - linear_interpolation(y1, y2, weight)
    - parabolic_interpolation(alpha, beta, gamma, x)
    
@author: Zoltan Galaz
@email:  z.galaz@phd.feec.vutbr.cz

Last revision: 04.05.2016
"""

# Imports
import os
import math
import pickle
import numpy as np
import scipy.io.wavfile as wavfile

# Functions
def read_wav_file(path): 
    """
    Load WAV file and return its properties
    
    Load WAV file specified by the path (string), and return a list of 1/more
    numpy array(s) containing the data samples for each channel, the sampling 
    frequency, the number of channels, and the number of samples per channel

    Args:
        path:         (string) string specifying the path to the WAV file
    
    Returns:
        data:         (list) list of 1/more numpy array(s)->data for channels
        fs:           (int) sampling frequency [Hertz]
        num_channels: (int) number of channels
        num_samples:  (int) number of samples per channel
        
    """
        
    # Parse the input file's extension
    extension = os.path.splitext(path)[1]
    
    # Load the WAV file and set the output parameters
    try:
        if extension.lower() == '.wav':
            [fs, x] = wavfile.read(path)
            num_samples = len(x)
            try: 
                num_channels = x.shape[1]
            except:
                num_channels = 1
            data = []    
            for channel in range(num_channels):
                if num_channels == 1:
                    data.append(x.astype(np.float32)/float(2**15))
                else:
                    data.append(x[0:,channel].astype(np.float32)/float(2**15))
        else:
            raise IOError("unknown file type")
            return (-1,-1,-1)
    except: 
        IOError("file not found")
        return (-1,-1,-1)
    
    # Return the output data (tuple)
    return (data, fs, num_channels, num_samples)

def make_directory(new_directory):
    """
    Make a new directory
    
    Source: GNU <http://code.activestate.com/recipes/82465-a-friendly-mkdir/>
    Make new directory using os.mrdir() function: new_directory (string). The
    function follows the rules standard for mkdir, i.e.:
        - if a directory already exists, silently complete
        - if a file with the specified name already exists raise an exception
        - parent directory(ies) does not exist, make them as well 
   
    Args:
        new_directory: (string) full path name of the directory to be created
    
    Returns:
        None
        
    """
    
    # Make new directory (if possible)
    if os.path.isdir(new_directory):
        pass
    elif os.path.isfile(new_directory):
        raise OSError("file with the same name exists")
    else:
        (head, tail) = os.path.split(new_directory)
        if head and not os.path.isdir(head):
            make_directory(head)
        if tail:
            os.mkdir(new_directory)

def save_data(data, file_name):
    """
    Save data structure as a Python pickle
   
    Args:
        data:      (writeable data structure) data to be saved
        file_name: (string) name of the file to be saved
    
    Returns:
        None
        
    """
    
    # Save the data structure
    fid = open(file_name, "w") 
    if fid:
        pickle.dump(data, fid)
        fid.close()
    else:
        raise Exception("unable to save data to file")

def load_data(file_name):
    """
    Load data from a Python pickle
   
    Args:       
        file_name: (string) name of the file to be loaded
    
    Returns:
        data:      (writeable data structure) data to be loaded
        
    """
    
    # Load the data structure
    fid = open(file_name, "w") 
    if fid:
        data = pickle.load(fid)
        fid.close()
        return data
    else:
        raise Exception("unable the data from file")
    
def find_peaks(data, sort=False, interpolate=True):
    """
    Find peaks in the data
    
    Find peaks (maximum values) in the provided data array. This function uses
    a range of iteration (start,end): point[index-1]..point[index+1] to search
    for the peaks. The optional variables enables a user to sort, interpolate
    the peaks.
   
    Args:
        data:        (numpy array) data array to find peaks in (1-dimensional)
        sort:        (boolean) sort the peaks according to values
        interpolate: (boolean) perform parabolic interpolation
    
    Returns:
        Tuple
        peaks_x:     (list) indices of the peaks
        peaks_y:     (list) values of the peaks
        
    """
    
    # Pefrorm initial check
    if type(data).__name__.strip() <> "ndarray":
        raise ValueError("data argument is not an instance of numpy.array")
    if len(data) < 1:
        raise ValueError("data array is empty")
    peaks_x = []
    peaks_y = []
    
    # Find peaks in the data
    for i in xrange(1, len(data)-1):
        if data[i] >= data[i-1] and data[i] >= data[i + 1]:
            x_pos_max = i
            value_max = data[i]
            
            # Interpolate (parabolic interpolation) if desired
            if interpolate:
                if x_pos_max > 0 and x_pos_max < len(data)-1:
                    alpha = data[x_pos_max-1]
                    beta  = data[x_pos_max]
                    gamma = data[x_pos_max+1]
                    denom = (alpha-beta*2+gamma)/2.0
                    if denom == 0.0: 
                        denom += 0.0001
                    x = (alpha-gamma)/denom
                    x_pos_max = x + x_pos_max
                    value_max = parabolic_interpolation(alpha, beta, gamma, x)
            peaks_x.append(x_pos_max)
            peaks_y.append(value_max)

    # Sort (ascending->according to peaks_y) if desired
    if sort:
        index = range(len(peaks_y))
        index.sort(key=peaks_y.__getitem__)
        peaks_x[:] = [peaks_x[i] for i in index]
        peaks_y[:] = [peaks_y[i] for i in index]
    
    # Return the peaks (positions, values)
    return (peaks_x, peaks_y)

def linear_interpolation(y1, y2, weight):
    """
    Perform linear interpolation
    
    Perform the linear interpolation between two equally space values (y1, y2) 
    and apply the weighting -> [0..1]: 0 = 100%y1, 1 = 100%y2.
   
    Args:
        y1:     (float) first data value
        y1:     (float) second data value
        weight: (float) weighting factor [0..1]
    
    Returns:
        lin:    (float) linearly interpolated data value
        
    """
    
    # Return linearly interpolated data value
    return y1*(1.0-weight)+y2*weight

def parabolic_interpolation(alpha, beta, gamma, x):
    """
    Perform parabolic interpolation
    
    Perform the parabolic interpolation between three equally space values and 
    apply relative position of read offset [-1..1]. The points: (alpha, beta, 
    gamma). 
   
    Args:
        alpha: (float) first data value
        beta:  (float) second data value
        gamma: (float) third data value
        x:     (float) relative position of read offset [-1..1]
    
    Returns:
        par:   (float) parabolically interpolated data value
        
    """
    
    # Perform initial check
    if x == 0:
        return beta
    else:
        offset = alpha
        if (beta < offset):
            offset = beta
        if (gamma < offset): 
            offset = gamma
    
    # Apply the offset
    offset = math.fabs(offset)+1
    alpha += offset;
    beta  += offset;
    gamma += offset;
    
    # Return parabolically interpolated data value
    a = (alpha-2.0*beta+gamma)/2.0
    if (a == 0):
        if (x > 1):
            return linear_interpolation(beta, gamma, x)-offset
        else:
            return linear_interpolation(alpha, beta, x+1)-offset
    else:
        c = (alpha-gamma)/(4.0*alpha)
        b = beta-a*(c**2)
        return (a*(x-c)*(x-c)+b)-offset



# Main
if __name__ == "__main__":
    from matplotlib.pyplot import figure, plot, show
    
    [x, Fs, num_channels, num_samples] = read_wav_file("data/alkoholik8k.wav")
    for channel in range(num_channels):
        data  = np.array(x[channel])
        peaks = find_peaks(data, sort=False, interpolate=False)
        
        figure()
        plot(data, 'b')
        plot(peaks[0], peaks[1], 'r')
        show()