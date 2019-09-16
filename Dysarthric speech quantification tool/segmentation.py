"""
Dysarthric speech processing toolbox
Segmentation

@functions:
    - segment(data, **kwargs)
    - _check_segmentation(data)

@author: Zoltan Galaz
@email:  z.galaz@phd.feec.vutbr.cz

Last revision: 04.05.2016
"""

# Imports
import numpy as np

# Constants
NUM_SEGMENTS = 100

# Functions
def segment(data, **kwargs):
    """
    Return input data into divided into short-time segments
    
    Input data can be can be 1D numpy array, or 2D numpy array The function  
    also handles additional arguments (**kwargs), which are used to specify 
    the length and overlap of the windowing function during segmentation.

    Args:
        data:           (numpy array) input data vector (one-dimensional)
        
        **kwargs
        window_length : (int) length of window (in samples)
        window_overlap: (int) window overlapp (in samples)

    Returns:
        M:              (numpy array) 2D numpy array of segmented signal. 
                         Rows: samples
                         Cols: segments
    
    Raises:
        ValueError:     handles the case: window overlap is >= window length
        
    """
    
    # Check if the segmentation is necessary
    if not _check_segmentation(data):
        return data
    
    # Get optional arguments
    window_length  = kwargs.get('window_length', \
                                 np.floor(data.size/NUM_SEGMENTS))
    window_overlap = kwargs.get('window_overlap', \
                                 np.floor(window_length/2.0))
                                
    window_length  = int(window_length)
    window_overlap = int(window_overlap)
                                
    # Check if the window length does not exceed length of input data
    if window_length > len(data):
        padding = True
    else:
        padding = False
        
    # Check if the window overlap does not exceed the window length
    if window_overlap >= window_length:
        raise ValueError("window overlap is >= window length")
    
    # Get the number of columns (segments)
    cols = (len(data)-window_overlap)/float((window_length-window_overlap))
    if padding:
        if (cols-np.fix(cols) > 0):
            data = data.extend(np.zeros((window_length,1), dtype=np.int))
            cols = np.ceil(cols)
    else:
        cols = np.fix(cols)
    
    # Prepare the segmentation parameters
    tmp_start = 0
    tmp_stop  = cols*(window_length-window_overlap)
    tmp_step  = int(window_length-window_overlap)
    step = [x for x in range(tmp_start,tmp_stop.astype(np.int64),tmp_step)]
    seq  = [x for x in range(tmp_start,window_length)]
    
    # Segment the sequence
    M = np.zeros((int(cols),window_length), dtype=float)
    for i in xrange(0,int(cols)):  
        M[i,:] = data[seq+np.array(step[i])].copy()
    M = M.transpose()
    
    # Weight by Hamming window
    for i in xrange(0,int(cols)):
        M[:,i] = M[:,i]*np.hamming(window_length)
        
    # Return the segmented data
    return M
    
def _check_segmentation(data):
    """
    Return if the segmentation is necessary
    
    The function tries to convert data into numpy array and checks if the data
    are one dimensional array. If it is, the segmentation is necessary. A user
    can therefore call segment(data) to segment the signal into segments of
    predefined length and overlap.

    Args:
        data:            (numpy array/list) input data (one/two-dimensional)

    Returns:
        do_segmentation: (boolean) segmentation is necessary
    
    Raises:
        ValueError:      handles the case: unsupported input file type'
        
    """
    
    # Convert input data to numpy array
    try:
        data = np.array(data)
    except ValueError:
        print "<Error> _check_segmentation(): unsupported input file type"

    # Check if the segmentation is necessary
    do_segmentation = False
    if isinstance(data, np.ndarray):
        if data.ndim == 1:
            do_segmentation = True
            
    # Return if the segmentation is necessary
    return do_segmentation
    
    
    
# Main
if __name__ == "__main__":
    from matplotlib.pyplot import figure, plot, subplot, tight_layout
    from matplotlib.pyplot import grid, title, xlabel, ylabel, show
    from utilities import read_wav_file
    
    file_name = "alkoholik8k.wav"
    [x, Fs, num_channels, num_samples] = read_wav_file("data/" + file_name)
    for channel in xrange(num_channels):
        data   = np.array(x[channel], dtype=np.float64)
        frames = segment(data, window_length=0.02*Fs, window_overlap=0.01*Fs)
        seg    = frames[:,1]
        
        figure()
        subplot(2,1,1)
        plot(data)
        grid()
        title("Original file: " + file_name)
        xlabel("samples")
        ylabel("amplitude")
        
        subplot(2,1,2)
        plot(seg)
        grid()
        title("1st segment: " + file_name)
        xlabel("samples")
        ylabel("amplitude")
        
        tight_layout()
        show()