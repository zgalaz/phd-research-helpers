"""
Dysarthric speech processing toolbox
Speaker class

@attributes:
    - name             : (string) speaker's name ('1st_name 2nd_name')
    - id               : (int) speaker's unique identifier
    - gender           : (string) speaker's gender ('male'/'female'}
    - age              : (int) speaker's age
    - speech           : (speech class object) speech signal
    - features         : (dict) speech features
    - condition        : (dict) medical examination details
                          .clinical_state: (string) ['HC'/'PD']
                          .disease_duration: (int)
                          .updrs_III: (float)
                          .updrs_IV: (float)
                          .rbdsq: (float)
                          .fog_q: (float)
                          .nmss: (float)
                          .mmse: (float)
                          .ace_r: (float)
                          .bdi: (float)
    
@methods:
    - __init__         : initialize a class (class constructor)
    - __repr__         : operator overloading (representation)
    - __str__          : operator overloading (string representation)
    - __eq__           : operator overloading (equal to)
    - __ne__           : operator overloading (not equal to)
    
@author: Zoltan Galaz
@email:  z.galaz@phd.feec.vutbr.cz

Last revision: 09.05.2016
"""

# Imports
from speech import Speech

# Classes
class Speaker(object):
    """
    Speaker class container
    
    @attributes:
        - name             : (string) speaker's name ('1st_name 2nd_name')
        - id               : (int) speaker's unique identifier
        - gender           : (string) speaker's gender ('male'/'female'}
        - age              : (int) speaker's age
        - speech           : (speech class object) speech signal
        - features         : (dict) speech features
        - condition        : (dict) medical examination details
                              .clinical_state: (string) ['HC'/'PD']
                              .disease_duration: (int)
                              .updrs_III: (float)
                              .updrs_IV: (float)
                              .rbdsq: (float)
                              .fog_q: (float)
                              .nmss: (float)
                              .mmse: (float)
                              .ace_r: (float)
                              .bdi: (float)
        
    @methods:
        - __init__         : initialize a class (class constructor)
        - __repr__         : operator overloading (representation)
        - __str__          : operator overloading (string representation)
        - __eq__           : operator overloading (equal to)
        - __ne__           : operator overloading (not equal to)
    
    """
    
    def __init__(self, **kwargs):
        """
        Initializes an object of Speaker class
        
        Construct an object of speaker class given input (**kwargs) atributes. 
        This function checks the properties of the input data, and initializes 
        all unspecified atributes to default values.
        
        Args:            
            **kwargs
            name:      (string) speaker's name ('1st_name 2nd_name')
            id:        (int) speaker's unique identifier
            gender:    (string) speaker's gender ('male'/'female'}
            age:       (int) speaker's age
            speech:    (speech class object) speech signal
            features:  (dict) speech features
            condition: (dict) medical examination details
                        .clinical_state: (string) ['HC'/'PD']
                        .disease_duration: (int) [months]
                        .updrs_III: (float)
                        .updrs_IV: (float)
                        .rbdsq: (float)
                        .fog_q: (float)
                        .nmss: (float)
                        .mmse: (float)
                        .ace_r: (float)
                        .bdi: (float)
            
        Returns:
            self:      (Speaker object) initialized Speaker class object
            
        """
        
        # Conctruct the Speaker class object
        self.name      = kwargs.get('name', None)
        self.id        = kwargs.get('id', None)
        self.gender    = kwargs.get('gender', None)
        self.age       = kwargs.get('age', None)
        self.speech    = kwargs.get('speech', None)
        self.features  = kwargs.get('features', None)
        self.condition = kwargs.get('condition', None)
        
    def __repr__(self):
        return "\nSpeaker:\n" + \
                "  name = '{0}'\n".format(self.name) + \
                "  id = '{0}'\n".format(self.id) + \
                "  age = '{0}'\n".format(self.age) + \
                "  gender = '{0}'\n".format(self.gender)
    
    def __str__(self):
        speaker_str = "\nSpeaker:\n" + \
            "  name = '{0}'\n".format(self.name) + \
            "  id = '{0}'\n".format(self.id) + \
            "  age = '{0}'\n".format(self.age) + \
            "  gender = '{0}'\n".format(self.gender)
        if self.speech is not None:
            speaker_str += self.speech.__str__()
        if self.condition is not None:
            for key, value in self.condition.items():
                speaker_str += "  {0} = '{1}'\n".format(key, value)
        if self.features is not None:
            for key, value in self.features.items():
                speaker_str += "  {0} = '{1}'\n".format(key, value)
        return speaker_str
    
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
    from utilities import read_wav_file
    
    [x, Fs, num_channels, num_samples] = read_wav_file("data/alkoholik8k.wav")
    wav = Speech(x, Fs)
    prs = Speaker(speech=wav)
    
    print prs