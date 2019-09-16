"""
Dysarthric speech processing toolbox
Speech signal class

@attributes:
    - data             : (list/tuple of numpy arrays, or numpy array)
                          1D numpy array: raw data samples
                          2D numoy array: segmented [samples, segments]
                          list/tuple: elements store data for channels
    - label            : (string) speech signal's label
    - fs               : (float) sampling frequency [Hertz]
    - num_channels     : (int) number of channels of data
    - num_samples      : (int) number of samples per channel
    - t_start          : (float) starting time [seconds]
    - t_end            : (float) ending time [seconds]
    - duration         : (float) duration [seconds]

@methods:
    - __init__         : initialize a class (class constructor)
    - __repr__         : operator overloading (representation)
    - __str__          : operator overloading (string representation)
    - __eq__           : operator overloading (equal to)
    - __ne__           : operator overloading (not equal to)
    - normalize_data   : normalize the data min_value <= data <= max_value
    - standardize_data : standardize the data mean() = 0, std() = 1
    - get_maximum      : returns the maximum values in data (per channel)
    - get_minimum      : returns the minimum values in data (per channel)
    - unwrap           : returns unwrapped data channel (np.array) 
    
@author: Zoltan Galaz
@email:  z.galaz@phd.feec.vutbr.cz

Last revision: 09.05.2016
"""

# Imports
import numpy as np

# Classes
class Speech(object):
    """
    Speech signal class container
    
    @attributes:
        - data             : (list/tuple of numpy arrays, or numpy array)
                              1D numpy array: raw data samples
                              2D numoy array: segmented [samples, segments]
                              list/tuple: elements store data for channels
        - label            : (string) speech signal's label
        - fs               : (float) sampling frequency [Hertz]
        - num_channels     : (int) number of channels of data
        - num_samples      : (int) number of samples per channel
        - t_start          : (float) starting time [seconds]
        - t_end            : (float) ending time [seconds]
        - duration         : (float) duration [seconds]
    
    @methods:
        - __init__         : initialize a class (class constructor)
        - __repr__         : operator overloading (representation)
        - __str__          : operator overloading (string representation)
        - __eq__           : operator overloading (equal to)
        - __ne__           : operator overloading (not equal to)
        - normalize_data   : normalize the data min_value <= data <= max_value
        - standardize_data : standardize the data mean() = 0, std() = 1
        - get_maximum      : returns the maximum values in data (per channel)
        - get_minimum      : returns the minimum values in data (per channel)
        - unwrap           : returns unwrapped data channel (np.array)
        
    """
    
    def __init__(self, data, fs, **kwargs):
        """
        Initializes an object of Speech class
        
        Construct an object of speech signal class given input atributes. This
        function checks the properties of the input data, and initializes all
        unspecified atributes to default values.
        
        Args:
            data:    (numpy array/list or tuple of numpy arrays) data samples 
                      1D numpy array: raw data samples
                      2D numoy array: segmented [samples, segments]
                      list/tuple: elements store data for channels
            fs:      (float) sampling frequency [Hertz]
            
            **kwargs
            label:   (string) speech signal's label
            t_start: (float) time of the start [seconds]
            t_end:   (float) time of the end [seconds]
            
        Returns:
            self:    (Speech object) initialized Speech class object
            
        """
        
        # Conctruct the Speech class object
        data_type = type(data).__name__
        self.data = []
        if data_type in ['list', 'tuple']:
            for i in xrange(len(data)):
                if not type(data[i]).__name__ == 'ndarray':
                    raise ValueError("expected a list of numpy arrays")  
                temp_data = np.array(data[i].copy(), dtype=np.float64)
                temp_data = temp_data[~np.isnan(temp_data)]
                temp_data = temp_data[~np.isinf(temp_data)]
                self.data.append(temp_data)
            self.num_channels = len(self.data)
            self.num_samples  = self.data[0].size*self.num_channels
        elif data_type == 'ndarray':
            self.data = [data[(~np.isnan(data)) & (~np.isinf(data))]]
            self.num_channels = 1
            self.num_samples  = self.data.size
        else:
            raise ValueError("unexpected data type encountered")
            
        self.fs       = float(abs(fs))    
        self.duration = self.num_samples/self.fs
        self.t_start  = kwargs.get('t_start', 0.0)
        self.t_end    = kwargs.get('t_end', self.t_start+self.duration)
        self.label    = kwargs.get('label', None)
        
    def normalize_data(self, min_value=0.0, max_value=1.0):
        """
        Normalize data of Speech class object
        
        Normalize the data of a given speech signal. Optional function args
        determines the boundaries: minimum, maximum of normalized data. The
        function handles multichannel signals.
        
        Args:
            self:      (Speech object) Speech class object
            min_value: (float) minimum value allowed in the data
            max_value: (float) maximum value allowed in the data
        
        Returns:
            self:      (Speech object) Speech class object
            
        """
        
        # Return the data normalized
        data_type = type(self.data).__name__
        if data_type in ['list', 'tuple']:
            for i in range(len(self.data)):
                try: 
                    num_segments = self.data[i].shape[1]
                except:
                    num_segments = 1
                for seg in xrange(0,num_segments):
                    try:
                        segment = self.data[i][:,seg]
                    except:
                        segment = self.data[i][seg,:]
                    if (segment.min() < min_value) or \
                       (segment.max() > max_value):
                        v_max = segment.max()
                        v_min = segment.min()
                        segment -= float(v_min)
                        segment *= (max_value-min_value)/(v_max-v_min)
                        segment += min_value
        else: 
            ValueError("unknown data type: list or tuple expected")
    
    def standardize_data(self):
        """
        Standardize data of Speech class object
        
        Standardize the data of a given speech signal to have a zero mean and
        a standard deviation of one (mean = 0, std = 1). The function handles 
        multichannel signals.
        
        Args:
            self: (Speech object) Speech class object
        
        Returns:
            self: (Speech object) Speech class object
         
        """
        
        # Return the data standardized
        data_type = type(self.data).__name__
        if data_type in ['list', 'tuple']:
            for i in range(len(self.data)):
                try: 
                    num_segments = self.data[i].shape[1]
                except:
                    num_segments = 1
                for seg in xrange(0,num_segments):
                    try:
                        segment = self.data[i][:,seg]
                    except:
                        segment = self.data[i][seg,:]
                    mean_val = np.mean(self.data[i])
                    std_val  = np.std(self.data[i])
                    segment -= mean_val
                    segment /= std_val
        else: 
            ValueError("unknown data type: list or tuple expected")
                
    def get_maximum(self):
        """
        Get the maximum values of Speech class object
        
        Get the maximum value stored in the data. Since the data can represent
        list/tuple of numpy arrays (numpy arrays hold samples of the channels)
        the function returns either a scalar value or a list of maximum values
        that holds the maximum value per channel.
        
        Args:
            self:      (Speech object) Speech class object
            
        Returns:
            max_value: (list) list of maximum values per channel
            
        """
        
        # Return the maximum value/s
        max_value = []
        data_type = type(self.data).__name__
        if data_type in ['list', 'tuple']:
            for i in range(len(self.data)):
                max_vals = []
                try: 
                    num_segments = self.data[i].shape[1]
                except:
                    num_segments = 1
                for seg in xrange(0,num_segments):
                    try:
                        segment = self.data[i][:,seg]
                    except:
                        segment = self.data[i][seg,:]
                    max_vals.append(segment.max())        
                max_value.append(max_vals)
            return max_value
        else: 
            ValueError("unknown data type: list or tuple expected")
    
    def get_minimum(self):
        """
        Get the minimum values of Speech class object
        
        Get the minimum value stored in the data. Since the data can represent
        list/tuple of numpy arrays (numpy arrays hold samples of the channels)
        the function returns either a scalar value or a list of minimum values
        that holds the minimum value per channel.
        
        Args:
            self:      (Speech object) Speech class object
            
        Returns:
            min_value: (list) list of minimum values per channel
            
        """
        
        # Return the minimum value/s
        min_value = []
        data_type = type(self.data).__name__
        if data_type in ['list', 'tuple']:
            for i in range(len(self.data)):
                min_vals = []
                try: 
                    num_segments = self.data[i].shape[1]
                except:
                    num_segments = 1
                for seg in xrange(0,num_segments):
                    try:
                        segment = self.data[i][:,seg]
                    except:
                        segment = self.data[i][seg,:]
                    min_vals.append(segment.min())    
                min_value.append(min_vals)
            return min_value
        else: 
            ValueError("unknown data type: list or tuple expected")
    
    def unwrap(self, index):
        """
        Unwrap the data samples of indexed channel
        
        Returns the data samples of the indexed channel unwrapped, which means
        the function converts the data from a list to numpy array of the type
        numpy.float64. It also checks the validity of the index.
        
        Args:
            self:      (Speech object) Speech class object
            index:     (int) index of the channel to be unwraooed
            
        Returns:
            unwrapped: (numpy array) unwrapped data samples
            
        """
        
        # Unwrap indexed data samples
        data_type = type(self.data).__name__
        if data_type in ['list', 'tuple']:
            if (index < 0) or (index >= self.num_channels):
                raise ValueError("index out of range")
            else:
                return np.array(self.data[index].copy(), dtype=np.float64)
        else: 
            ValueError("unknown data type: list or tuple expected")
    
    def __repr__(self):
        return "  data = '{0}'\n".format(self.data) + \
               "  fs = '{0}'\n".format(self.fs)
    
    def __str__(self):
        return "  data = '{0}'\n".format(self.data) + \
               "  fs = '{0}'\n".format(self.fs) + \
               "  num_channels = '{0}'\n".format(self.num_channels) + \
               "  num_samples = '{0}'\n".format(self.num_samples) + \
               "  t_start = '{0}'\n".format(self.t_start) + \
               "  t_end = '{0}'\n".format(self.t_end) + \
               "  duration = '{0}'\n".format(self.duration)
    
    def __hash__(self):
        return hash(tuple(sorted(self.__dict__.items())))
    
    def __eq__(self, other):
        if isinstance(other, self.__class__):
            return self.__dict__ == other.__dict__
        return NotImplemented
    
    def __ne__(self, other):
        if isinstance(other, self.__class__):
            return not self == other
        return NotImplemented
    
if __name__ == "__main__":
    from matplotlib.pyplot import figure, plot, subplot, tight_layout, show
    from segmentation import segment
    from utilities import read_wav_file
    from features import spectral_centroid
    
    [x, Fs, num_channels, num_samples] = read_wav_file("data/alkoholik8k.wav")
    wav = Speech(x, Fs)
    
    for channel in xrange(num_channels):
        data = wav.unwrap(channel)
        seg  = segment(data, window_length=0.02*Fs, window_overlap=0.01*Fs)
        SC   = spectral_centroid(seg, Fs)
        
        figure()
        subplot(2,1,1)
        plot(data)
        subplot(2,1,2)
        plot(SC)
        tight_layout()
        show()