'''
Created on 25 Jan 2012

@author: garyt
'''

class Constants_table(object):
    '''
    classdocs
    '''

    FLAT_BOTTOM_LIMIT = 'flat_bottom_limit'
    FLAT_BOTTOM_CONSTANT  = 'flat_bottom_constant'
    END_HARMONIC = 'end_harmonic'
    SCALE_HARMONIC = 'scale_harmonic'
    WEIGHT = 'weight'
    TANH_AMPLITUDE = 'tanh_amplitude'
    
    def __init__(self, table):
        self._table = table
    
    def get_flat_bottom_constant(self):
        return self._table[self.FLAT_BOTTOM_CONSTANT]
    
    def get_target_atoms(self):
        return self._table[self.FLAT_BOTTOM_LIMIT].keys()
    
    def _check_target_atom(self, atom):
        
        atoms = self.get_target_atoms()
        if not atom in atoms:
            message = "atom %s not in target atoms %s" % (atom,', '.join(atoms))
            raise KeyError(message)
    
    def get_flat_bottom_limit(self, target_atom):
        self._check_target_atom(target_atom)
        return self._table[self.FLAT_BOTTOM_LIMIT][target_atom]
        
    def get_end_harmonic(self, target_atom):
        self._check_target_atom(target_atom) 
        return self._table[self.END_HARMONIC][target_atom]   

    def get_scale_harmonic(self, target_atom):
        self._check_target_atom(target_atom) 
        return self._table[self.SCALE_HARMONIC][target_atom]
    
    def get_weight(self,target_atom):
        
        self._check_target_atom(target_atom)
        
        return self._table[self.WEIGHT][target_atom]
    def get_tanh_amplitude(self,target_atom):
        
        self._check_target_atom(target_atom)
        return self._table[self.TANH_AMPLITUDE][target_atom]
    
    def get_tanh_elongation(self,target_atom):
        self._check_target_atom(target_atom)
        
        tanh_amplitude =self.get_tanh_amplitude(target_atom)
        end_harmonic = self.get_end_harmonic(target_atom)
        scale_harmonic  =  self.get_scale_harmonic(target_atom)
        
        return 2.0 * end_harmonic / ((scale_harmonic**2.0) * tanh_amplitude);
    
    def get_tanh_y_offset(self,target_atom):
        self._check_target_atom(target_atom)
        
        end_harmonic = self.get_end_harmonic(target_atom)
        scale_harmonic  =  self.get_scale_harmonic(target_atom)
        
        return (end_harmonic / scale_harmonic)**2.0
