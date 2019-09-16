"""
Dysarthric speech processing toolbox
Features

@functions:
    - squared_energy_operator(frames)
    - low_squared_energy_operator(frames)
    - teager_energy_operator(frames)
    - zero_crossing_rate(frames)
    - high_zero_crossing_rate(frames)
    - vowel_articulation_index(F1a, F1i, F1u, F2a, F2i, F2u)
    - vowel_space_area(F1a, F1i, F1u, F2a, F2i, F2u)
    - jitter_variants(F0)
    - shimmer_variants(A0)
    - spectral_centroid(frames, Fs)
    - spectral_spread(frames, Fs)
    - spectral_flux(frames)
    - spectral_rolloff(frames, c=0.8)
    - spectral_distance(frames)
    - spectral_flatness(frames)
 
@author: Zoltan Galaz
@email:  z.galaz@phd.feec.vutbr.cz

Last revision: 09.05.2016
"""

# Imports
import numpy as np
import segmentation as seg

from scipy.fftpack import fft
from math import exp
from math import log
from math import log10

# Functions
def squared_energy_operator(frames):
    """
    Return the squared energy operator of an input signal
    
    The input signal can be an 1D numpy array (waveform), or 2D numpy array.
    If the input signal is one-dimensional, function performs a segmentation
    to automatically segment the signal into frames of predefined length and
    overlap. Otherwise, function computes squared energy operator (SEO) per  
    segment. The output values are stored in 1D numpy array (elements hold  
    the SEO of all segments).

    Args:
        frames: (numpy array) input data vector (one/two-dimensional)
        
    Returns:
        SEO:    (numpy array) squared energy operator vector
        
    """
    
    # Perform segmentation if necesary
    frames = seg.segment(frames)
    
    # Compute squared energy operator (per segment)
    SEO = []
    for i in xrange(0,frames.shape[1]):
        frame  = frames[:,i].copy()
        energy = np.sum(frame**2)/np.float64(len(frame))
        SEO.append(energy)
        
    # Return squared energy operator vector
    return np.array(SEO)

def low_squared_energy_operator(frames):
    """
    Return the low squared energy operator of an input signal
    
    The input signal can be an 1D numpy array (waveform), or 2D numpy array.
    If the input signal is one-dimensional, function performs a segmentation
    to automatically segment the signal into frames of predefined length and
    overlap. Otherwise, function computes low squared energy operator (LSEO)   
    per segment. The output values are stored in 1D numpy array (elements   
    hold the LSEO of all segments).

    Args:
        frames: (numpy array) input data vector (one/two-dimensional)
        
    Returns:
        LSEO:   (numpy array) low squared energy operator vector
        
    """
    
    # Perform segmentation if necesary
    frames = seg.segment(frames)
    
    # Compute and segment squared energy operator
    LEN = 100
    SEO = seg.segment(squared_energy_operator(frames), LEN, 0)
    
    # Compute low squared energy operator
    LSEO = np.sum(np.sign(0.5*np.mean(SEO, axis=0)-SEO)+1)/(2*LEN);
    
    # Return low squared energy operator vector
    return np.array(LSEO)

def teager_energy_operator(frame):
    """
    Return the teager energy operator of an input signal
    
    The input signal is an 1D numpy array (waveform). This function computes  
    teager energy operator (TEO) and stores the output values in an 1D numpy
    array. In this function, no segmentation is performed.

    Args:
        frame: (numpy array) input data vector (one-dimensional)
        
    Returns:
        TEO:   (numpy array) teager energy operator vector
        
    """
    
    # Compute teager energy operator vector
    TEO = np.zeros((len(frame)-1, 1), dtype=np.float64);
    for i in xrange(1,len(frame)-1):
        a_sample = float(frame[i].copy())
        p_sample = float(frame[i-1].copy())
        n_sample = float(frame[i+1].copy())
        TEO[i-1] = a_sample**2 - n_sample*p_sample
        
    # Return teager energy operator vector
    return np.array(TEO)

def zero_crossing_rate(frames):
    """
    Return the zero-crossing rate of an input signal
    
    The input signal can be an 1D numpy array (waveform), or 2D numpy array.
    If the input signal is one-dimensional, function performs a segmentation
    to automatically segment the signal into frames of predefined length and
    overlap. Otherwise, the function computes zero-crossing rate (ZCR) per 
    segment. The output values are stored in 1D numpy array (elements hold 
    the ZCR of all segments).

    Args:
        frames: (numpy array) input data vector (one/two-dimensional)
        
    Returns:
        ZCR:    (numpy array) zero-crossing rate vector
        
    """
    
    # Perform segmentation if necesary
    frames = seg.segment(frames)
    
    # Compute zero-crossing rate (per segment)
    ZCR = []
    for i in xrange(0,frames.shape[1]):
        frame = frames[:,i].copy()
        count = np.sum(np.abs(np.diff(np.sign(frame))))/2
        ZCR.append(np.float64(count)/np.float64(len(frame)-1.0))
        
    # Return zero-crossing rate vector
    return np.array(ZCR)

def high_zero_crossing_rate(frames):
    """
    Return the high zero-crossing rate of an input signal
    
    The input signal can be an 1D numpy array (waveform), or 2D numpy array.
    If the input signal is one-dimensional, function performs a segmentation
    to automatically segment the signal into frames of predefined length and
    overlap. Otherwise, the function computes high zero-crossing rate (HZCR)  
    per segment. The output values are stored in 1D numpy array (elements  
    hold the HZCR of all segments).

    Args:
        frames: (numpy array) input data vector (one/two-dimensional)
        
    Returns:
        HZCR:   (numpy array) high zero-crossing rate vector
        
    """
    
    # Perform segmentation if necesary
    frames = seg.segment(frames)
    
    # Compute and segment zero crossing rate
    LEN = 100
    ZCR = seg.segment(zero_crossing_rate(frames), LEN, 0)
    
    # Compute high zero-crossing rate
    HZCR = np.sum(np.sign(ZCR-1.5*np.mean(ZCR, axis=0))+1)/(2*LEN);
    
    # Return high zero-crossing rate vector
    return np.array(HZCR)

def vowel_articulation_index(F1a, F1i, F1u, F2a, F2i, F2u):
    """
    Return vowel articulation index

    Args:
        F1a: (float) the 1. formant frequency of the vowel /a [Hz]
        F1i: (float) the 1. formant frequency of the vowel /i [Hz]
        F1u: (float) the 1. formant frequency of the vowel /u [Hz]
        F2a: (float) the 1. formant frequency of the vowel /a [Hz]
        F2i: (float) the 1. formant frequency of the vowel /i [Hz]
        F2u: (float) the 1. formant frequency of the vowel /u [Hz]
        
    Returns:
        VAI: (float) vowel articulation index
        
    """
    
    # Return vowel articulation index
    return float((F2i+F1a)/(F1i+F1u+F2u+F2a))

def vowel_space_area(F1a, F1i, F1u, F2a, F2i, F2u):
    """
    Return vowel space area

    Args:
        F1a: (float) the 1. formant frequency of the vowel /a [Hz]
        F1i: (float) the 1. formant frequency of the vowel /i [Hz]
        F1u: (float) the 1. formant frequency of the vowel /u [Hz]
        F2a: (float) the 1. formant frequency of the vowel /a [Hz]
        F2i: (float) the 1. formant frequency of the vowel /i [Hz]
        F2u: (float) the 1. formant frequency of the vowel /u [Hz]
        
    Returns:
        VSA: (float) vowel space area
        
    """
    
    # Compute vowel space area
    EDiu = np.sqrt((F1i-F1u)**2+(F2i-F2u)**2)
    EDia = np.sqrt((F1i-F1a)**2+(F2i-F2a)**2)
    EDau = np.sqrt((F1a-F1u)**2+(F2a-F2u)**2)
    S    = (EDiu+EDia+EDau)/(2.0)
    VSA  = np.sqrt(S*(S-EDiu)*(S-EDia)*(S-EDau))

    # Return vowel space area
    return float(VSA)

def jitter_variants(F0):
    """
    Return jitter variants
    
    The input signal can be 1D numpy array of pitch periods estimations. The 
    function computes 4 variants of jitter: absolute, relative, rap, ppq5.

    Args:
        F0:            (numpy array) input pitch periods vector
        
    Returns:
        Tuple
        jitt_absolute: (float) jitter (absolute)
        jitt_relative: (float) jitter (relative)
        jitt_rap:      (float) jitter rap
        jitt_ppq:      (float) jitter ppq5
        
    """
    
    # Compute jitter (absolute)
    jitt_absolute = np.sum(np.abs(np.diff(F0)))/np.float64(len(F0)-1.0)
    
    # Compute jitter (relative)
    jitt_relative = jitt_absolute/np.mean(F0)*100
    
    # Compute jitter (rap)
    inn_rap  = 0.0
    for i in xrange(1, len(F0)-1):
        inn_rap += abs(F0[i]-sum(F0[i-1:i+1].copy())/3.0)
    jitt_rap = (np.sum(inn_rap)\
                /np.float64(len(F0)-1.0))\
                /np.mean(F0)*100
    
    # Compute jitter (ppq5)
    inn_ppq  = 0.0
    for i in xrange(2, len(F0)-2):
        inn_ppq += abs(F0[i]-sum(F0[i-2:i+2].copy())/5.0)
    jitt_ppq = (np.sum(inn_ppq)\
                /np.float64(len(F0)-1.0))\
                /np.mean(F0)*100
    
    # Return jitter variants
    return (jitt_absolute, jitt_relative, jitt_rap, jitt_ppq)

def shimmer_variants(A0):
    """
    Return shimmer variants
    
    The input signal can be 1D numpy array of amplitudes estimations. The 
    function computes 4 variants of shimmer: local, local dB, apq3, apq5.

    Args:
        A0:          (numpy array) input amplitudes vector
        
    Returns:
        Tuple
        shimm_local: (float) shimmer (local)
        shimm_db:    (float) shimmer (local, dB)
        shimm_apq3:  (float) shimmer apq3
        shimm_apq5:  (float) shimmer apq5
    
    """
    
    # Compute shimmer (local)
    shimm_local = np.sum(np.abs(np.diff(A0)))/np.float64(len(A0)-1.0)
    
    # Compute shimmer (local, dB)
    inter = []
    for i in xrange(1, len(A0)):
        inner = float(A0[i-1])/float(A0[i]+0.001)
        if inner == 0.0 or np.isnan(inner) or np.isinf(inner):
            continue
        else:
            inter.append(20*log10(inner))
    shimm_db = np.sum(inter)/np.float64(len(A0)-1.0)
    
    # Compute shimmer (apq3)
    inn_apq3 = 0.0
    for i in xrange(1, len(A0)-1):
        inn_apq3 += abs(A0[i]-sum(A0[i-1:i+1].copy())/3.0)
    shimm_apq3 = (np.sum(inn_apq3)\
                  /np.float64(len(A0)-1.0))\
                  /np.mean(A0)*100
    
    # Compute shimmer (apq5)
    inn_apq5 = 0.0
    for i in xrange(2, len(A0)-2):
        inn_apq5 += abs(A0[i]-sum(A0[i-2:i+2].copy())/5.0)
    shimm_apq5 = (np.sum(inn_apq5)\
                  /np.float64(len(A0)-1.0))\
                  /np.mean(A0)*100
    
    # Return shimmer variants
    return (shimm_local, shimm_db, shimm_apq3, shimm_apq5)
    
def spectral_centroid(frames, Fs):
    """
    Return the spectral centroid of an input signal
    
    The input signal can be an 1D numpy array (waveform), or 2D numpy array.
    If the input signal is one-dimensional, function performs a segmentation
    to automatically segment the signal into frames of predefined length and
    overlap. Otherwise, function computes spectral centroid (SC) per segment
    The output values are stored in 1D numpy array (elements hold the SC of
    all segments).
    
    Args:
        frames: (numpy array) input data vector (one/two-dimensional)
        Fs:     (int) sampling rate (frequency)
        
    Returns:
        SC:     (numpy array) spectral centroid vector
        
    """
    
    # Perform segmentation if necesary
    frames = seg.segment(frames)
    
    # Compute spectral centroid (per segment)
    SC = []
    ID = (np.arange(1,frames.shape[0]+1))*(Fs/(2.0*frames.shape[0]))
    for i in xrange(0,frames.shape[1]):
        frame  = frames[:,i].copy()
        spect  = abs(fft(frame)[0:frames.shape[0]])
        spect /= spect.max()
    
        NUM = np.sum(ID*spect)
        DEN = np.sum(spect)+0.010
        SC.append((NUM/DEN)/(Fs/2.0))

    # Return pectral centroid vector
    return np.array(SC)

def spectral_spread(frames, Fs):
    """
    Return the spectral spread of an input signal
    
    The input signal can be an 1D numpy array (waveform), or 2D numpy array.
    If the input signal is one-dimensional, function performs a segmentation
    to automatically segment the signal into frames of predefined length and
    overlap. Otherwise, function computes spectral spread (SS) per segment.  
    The output values are stored in 1D numpy array (elements hold the SS of 
    all segments).

    Args:
        frames: (numpy array) input data vector (one/two-dimensional)
        Fs:     (int) sampling rate (frequency)
        
    Returns:
        SS:     (numpy array) spectral spread vector
        
    """
    
    # Perform segmentation if necesary
    frames = seg.segment(frames)
    
    # Compute spectral spread (per segment)
    SS = []
    ID = (np.arange(1,frames.shape[0]+1))*(Fs/(2.0*frames.shape[0]))
    for i in xrange(0,frames.shape[1]):
        frame  = frames[:,i].copy()
        spect  = abs(fft(frame)[0:frames.shape[0]])
        spect /= spect.max()
    
        NUM = np.sum(ID*spect)
        DEN = np.sum(spect)+0.010
        SC  = ((NUM/DEN)/(Fs/2.0))
        SS.append(np.sqrt(np.sum(((ID-SC)**2)*spect)/DEN))

    # Return spectral spread vector
    return np.array(SS)

def spectral_flux(frames):
    """
    Return the spectral flux of an input signal
    
    The input signal can be an 1D numpy array (waveform), or 2D numpy array.
    If the input signal is one-dimensional, function performs a segmentation
    to automatically segment the signal into frames of predefined length and
    overlap. Otherwise, function computes spectral flux (SF) per segment.  
    The output values are stored in 1D numpy array (elements hold the SF of 
    all segments).

    Args:
        frames: (numpy array) input data vector (one/two-dimensional)
        
    Returns:
        SF:     (numpy array) spectral flux vector
        
    """
    
    # Perform segmentation if necesary
    frames = seg.segment(frames)
    
    # Compute spectral flux (per segments)
    SF = []
    for i in xrange(1,frames.shape[1]):
        p_frame  = frames[:,i-1].copy()
        p_spect  = abs(fft(p_frame)[0:frames.shape[0]])
        p_spect /= p_spect.max()
        
        a_frame  = frames[:,i].copy()
        a_spect  = abs(fft(a_frame)[0:frames.shape[0]])
        a_spect /= a_spect.max()
        
        a_overall = a_spect/np.sum(a_spect)
        p_overall = p_spect/np.sum(p_spect)
        SF.append(np.sum((a_overall - p_overall)**2))
        
    # Return spectral flux vector
    return np.array(SF)

def spectral_rolloff(frames, c=0.8):
    """
    Return the spectral roll-off of an input signal
    
    The input signal can be an 1D numpy array (waveform), or 2D numpy array.
    If the input signal is one-dimensional, function performs a segmentation
    to automatically segment the signal into frames of predefined length and
    overlap. Otherwise, function computes spectral roll-off (SR) per segment  
    The output values are stored in 1D numpy array (elements hold the SR of 
    all segments).

    Args:
        frames: (numpy array) input data vector (one/two-dimensional)
        c:      (float) roll-off factor (default=0.8)
        
    Returns:
        SR:     (numpy array) spectral roll-off vector
        
    """
    
    # Perform segmentation if necesary
    frames = seg.segment(frames)
    
    # Compute spectral roll-off (per segment)
    SR = []
    for i in xrange(0,frames.shape[1]):
        frame  = frames[:,i].copy()
        spect  = abs(fft(frame)[0:frames.shape[0]])
        spect /= spect.max()
        
        [pos,] = np.nonzero(np.cumsum(spect**2) > c*np.sum(spect**2))
        if len(pos) > 0:
            SR.append(np.float64(pos[0])/(float(len(spect))))
        else:
            SR.append(0.0)
            
    # Return spectral roll-off vector
    return np.array(SR)

def spectral_distance(frames):
    """
    Return the spectral distance of an input signal
    
    The input signal can be an 1D numpy array (waveform), or 2D numpy array.
    If the input signal is one-dimensional, function performs a segmentation
    to automatically segment the signal into frames of predefined length and
    overlap. Otherwise, function computes spectral distance (SD) per segment  
    The output values are stored in 1D numpy array (elements hold the SD of 
    all segments).

    Args:
        frames: (numpy array) input data vector (one/two-dimensional)
        
    Returns:
    Tuple composed of
        SD_mag: (numpy array) spectral distance based on module vector
        SD_phs: (numpy array) spectral distance based on phase vector
        
    """
    
    # Perform segmentation if necesary
    frames = seg.segment(frames)
    
    # Compute spectral distance (per segment)
    SD_mag = []
    SD_phs = []
    for i in xrange(1,frames.shape[1]):
        p_frame = frames[:,i-1].copy()
        ps_mag  = np.abs(fft(p_frame)[0:frames.shape[0]])
        ps_phs  = np.angle(fft(p_frame)[0:frames.shape[0]])
        ps_mag /= ps_mag.max()
        ps_phs /= ps_phs.max()
        
        a_frame = frames[:,i].copy()
        as_mag  = np.abs(fft(a_frame)[0:frames.shape[0]])
        as_phs  = np.angle(fft(a_frame)[0:frames.shape[0]])
        as_mag /= as_mag.max()
        as_phs /= as_phs.max()
        
        SD_mag.append(np.sum(abs(ps_mag-as_mag)))
        SD_phs.append(np.sum(abs(ps_phs-as_phs)))
        
    # Return pectral tuple of distance vectors
    return (np.array(SD_mag), np.array(SD_phs))

def spectral_flatness(frames):
    """
    Return the spectral flatness of an input signal
    
    The input signal can be an 1D numpy array (waveform), or 2D numpy array.
    If the input signal is one-dimensional, function performs a segmentation
    to automatically segment the signal into frames of predefined length and
    overlap. Otherwise, function computes spectral flatness (SF) per segment  
    The output values are stored in 1D numpy array (elements hold the SF of 
    all segments).

    Args:
        frames: (numpy array) input data vector (one/two-dimensional)
        
    Returns:
        SF:     (numpy array) spectral flatness vector
        
    """
    
    # Perform segmentation if necesary
    frames = seg.segment(frames)
    
    # Compute spectral flatness (per segment)
    SF = []
    for i in xrange(0,frames.shape[1]):
        frame  = frames[:,i].copy()
        spect  = abs(fft(frame)[0:frames.shape[0]])**2
        spect /= spect.max()
        
        gMean = np.float64(0)
        aMean = 0 
        for i in range(len(spect)):
            sample = np.float64(spect[i])
            gMean += np.float64(log(sample))
            aMean += sample
        
        gMean /= np.float64(len(spect))
        gMean  = exp(gMean)
        aMean /= float(len(spect))
        SF.append(gMean/aMean)
        
    # Return spectral flatness vector
    return np.array(SF)



# Main
if __name__ == "__main__":
    from matplotlib.pyplot import figure, plot, subplot, tight_layout
    from matplotlib.pyplot import grid, title, xlabel, ylabel, show
    from segmentation import segment
    from utilities import read_wav_file
    
    file_name = "alkoholik8k.wav"
    [x, Fs, num_channels, num_samples] = read_wav_file("data/" + file_name)
    for channel in xrange(num_channels):
        data   = np.array(x[channel], dtype=np.float64)
        frames = segment(data, window_length=0.02*Fs, window_overlap=0.01*Fs)
        
        # Speech features
        SEO = squared_energy_operator(frames)
        TEO = teager_energy_operator(data)
        ZCR = zero_crossing_rate(frames)
        SC  = spectral_centroid(frames, Fs)
        SS  = spectral_spread(frames, Fs)
        SF  = spectral_flux(frames)
        SR  = spectral_rolloff(frames)
        (SD_mag, SD_phs) = spectral_distance(frames)
        
        figure()
        subplot(2,2,1)
        plot(np.linspace(0, len(data)/Fs, num=len(data)), data)
        grid()
        title("Original file: " + file_name)
        xlabel("time [s]")
        ylabel("amplitude")
        
        subplot(2,2,2)
        plot(SF)
        grid()
        title("Spectral flux")
        xlabel("segments")
        ylabel("spectral flux")
        
        subplot(2,2,3)
        plot(SEO)
        grid()
        title("Squared energy operator")
        xlabel("segments")
        ylabel("squared energy oper.")
        
        subplot(2,2,4)
        plot(SC)
        grid()
        title("Spectral centroid")
        xlabel("segments")
        ylabel("spectral centroid")
        
        tight_layout()
        show()