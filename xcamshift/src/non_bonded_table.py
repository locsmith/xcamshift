'''
Created on 19 Jan 2012

@author: garyt
'''



class Non_bonded_table(object):



    TARGET_ATOMS = "target_atoms"
    SPHERE_1 = 'sphere_1'
    SPHERE_2 = 'sphere_2'
    EXPONENT = 'exponent'
    DATA = 'data'
    CHEM_TYPE_TRANSLATIONS = 'chem_type_translations'
    
    def __init__(self, table):
        self._table = table
    
    def get_exponent(self,sphere):
        self._check_sphere(sphere)
        return self._table[self.DATA][sphere][self.EXPONENT]
    
    def get_target_atoms(self):
        return self._table[self.TARGET_ATOMS]    

    def _check_target_atom(self, atom_name):
        atoms = self.get_target_atoms()
        if atom_name not in atoms:
            template = "atom_name %s is not in sidechain shift _table target atoms (%s)"
            message = template % (atom_name, ', '.join(atoms))
            raise KeyError(message)
    
    def get_spheres(self):
        spheres = self._table[self.DATA].keys()
        spheres.sort()
        return tuple(spheres)


    def _check_sphere(self,sphere):
        spheres = self.get_spheres()
        if sphere not in spheres:
            template = "sphere %s is not in non bonded table spheres (%s)"
            message = template % (sphere, ', '.join(spheres))
            raise KeyError(message)
    
    def _check_coefficient_key(self,sphere,coefficent_key):
        self._check_sphere(sphere)
        coefficent_keys = self.get_remote_atom_types(sphere)
        if coefficent_key not in coefficent_keys:
            template = "coefficient key %s is not in non bonded table coefficient keys (%s)"
            coefficent_keys_string = [`coefficent_key` for coefficent_key in coefficent_keys]
            message = template % (coefficent_key, ', '.join(coefficent_keys_string))
            raise KeyError(message)
    
    def get_remote_atom_types(self, sphere):
        self._check_sphere(sphere)
        coefficent_keys = self._table[self.DATA][sphere].keys()
        coefficent_keys.remove('exponent')
        coefficent_keys.sort()
        return tuple(coefficent_keys) 
    
    def get_chem_types(self):
        return self._table[self.CHEM_TYPE_TRANSLATIONS].keys()
    
    def _check_chem_type(self,chem_type):
        chem_types = self.get_chem_types()
        if chem_type not in chem_types:
            template = "chemical type %s is not in non bonded table chemical type translations (%s)"
            chem_types_strings = [`chem_type` for chem_type in chem_types]
            message = template % (chem_types_strings, ', '.join(chem_types_strings))
            raise KeyError(message)
        
    def get_chem_type_translation(self, chem_type):
        self._check_chem_type(chem_type)
        
        return self._table[self.CHEM_TYPE_TRANSLATIONS][chem_type]
    
    def get_non_bonded_coefficient(self,target_atom,sphere,remote_atom_type,hybridisation):
        self._check_target_atom(target_atom)
        self._check_sphere(sphere)
        
        coefficent_key = remote_atom_type,hybridisation
        self._check_coefficient_key(sphere,coefficent_key)
        
        return self._table[self.DATA][sphere][coefficent_key][target_atom]
    
