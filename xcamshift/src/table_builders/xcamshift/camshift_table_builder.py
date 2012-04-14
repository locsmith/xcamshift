'''
Created on 30 Mar 2012

@author: garyt




'''
import sys
import os
import glob
import re
from table_builders.xcamshift.Backbone_distance_extractor import BB_table_extractor
from table_builders.xcamshift.Extra_distance_extractor import XTRA_table_extractor
from table_builders.xcamshift.Sidechain_distance_extractor import SC_table_extractor
from table_builders.xcamshift.Dihdedral_distance_extractor import DIHEDRALS_table_extractor
from table_builders.xcamshift.Ring_table_extractor import RING_table_extractor
from table_builders.xcamshift.Nonbonded_table_extractor import Nonbonded_table_extractor
from common_constants import CAMSHIFT_SUB_POTENTIALS

def get_atom_type_for_filename(file_path):
    file_name = os.path.split(file_path)[1]
    file_pattern = re.compile("1([A-Z]{1,2})orgpairs.par.*")
    match = re.search(file_pattern, file_name)
    return match.group(1)
    
def check_args():
    
    args = []
    if len(sys.argv) > 1:
        args = sys.argv[1:]
    if len(args) != 1:
        template = "Camshift_table_builder requires a single directory as its argument but got '%s'"
        print >> sys.stderr, template % `args`
        print >> sys.stderr, "exiting..."
        sys.exit(-1)
        
def parameter_file_iterator(table_dir, template = '1[A-Z]*orgpairs.par*'):
    if not os.path.isabs(table_dir):
        table_dir = os.path.join(os.getcwd(), table_dir)
    query = os.path.join(table_dir, template)
    for filename in glob.glob(query):
        yield filename



class Camshift_table_reader(object):
    '''
    classdocs
    '''


    def __init__(self,table_dir):
        '''
        Constructor
        '''
        
        self.table_dir =table_dir
        self.files_by_type = {}
        self.data = {}
        self.file_types = set()

    def read(self):
        self._find_files(self.table_dir)
        
        for file_type in self.file_types:
            self._parse_files(file_type)

    def _split_files_by_extension(self,par_files):
        for par_filename in par_files:
            ext = os.path.splitext(par_filename)
            ext = ext[1][1:]
            if ext == 'par':
                ext = ''
            self.file_types.add(ext)
            self.files_by_type.setdefault(ext, []).append(par_filename)

    def _find_files(self,table_dir):
    
        par_files = [file_name for file_name in parameter_file_iterator(table_dir)]
        
        self._split_files_by_extension(par_files)
    
    
    def _parse_files(self,file_type):
        self.data.setdefault(file_type,{})
        for file_path in self.files_by_type[file_type]:
            atom_type = get_atom_type_for_filename(file_path)
            file_h = open(file_path)
            
            for line in file_h:
                
                line = line.strip()
                if line.startswith("#"):
                    continue
                
                items = line.split()
                
                table_type = items[0]
                
                values = [float(item) for item in items[1:]]
                self.data[file_type].setdefault(table_type, {})[atom_type] = values
        
    
    def get_file_types(self):
        return self.file_types


def _get_file_type_name(file_type):
    if file_type == "":
        file_type_name = 'base'
    else:
        file_type_name = file_type
    
    return file_type_name

def _read_data_into_reader(table_dir):
    reader = Camshift_table_reader(table_dir)
    reader.read()
    return reader

def _get_extractor_classes(sub_potential = None):
    if isinstance(sub_potential, str):
        sub_potential = [sub_potential]
    
    extractor_classes = (
                  BB_table_extractor, 
                  XTRA_table_extractor, 
                  SC_table_extractor, 
                  DIHEDRALS_table_extractor, 
                  RING_table_extractor,
                  Nonbonded_table_extractor
                  )
    result  = []
    
    for extractor in extractor_classes:
        if sub_potential == None:
            result.append(extractor)
        elif extractor.get_name() in sub_potential:
            result.append(extractor)

    return result

# TODO integrate with rest of module
def extract(data_dir, sub_potential = None, residue_types=None):
    reader = _read_data_into_reader(data_dir)
    
    extractor_classes = _get_extractor_classes(sub_potential)
    
    if residue_types == None:
        residue_types = reader.get_file_types()
    elif isinstance(residue_types, str):
        residue_types = [residue_types]
        
    result  = {}   
    for extractor_class in extractor_classes:
        extractor =  extractor_class(reader.data)
        for residue_type in residue_types:
            result[sub_potential,residue_type] = extractor.extract(residue_type)
        
    return result
    
    
    
    
if __name__ == '__main__':

    check_args()
    
    table_dir = sys.argv[1]
    
    reader = _read_data_into_reader(table_dir)
    
    file_types = reader.get_file_types()
    
    table_types = CAMSHIFT_SUB_POTENTIALS
    
    #TODO use extract function here
    extractor_classes = _get_extractor_classes()

    for extractor_class, table_type in zip(extractor_classes,table_types):
        extractor = extractor_class(reader.data)
        
        for file_type in file_types:
            file_type_name  = _get_file_type_name(file_type)
            
            print "\n--- %s %s ---\n" % (table_type, file_type_name)
    
#            extractor = BB_table_extractor(reader.data)
            print extractor.extract(file_type)
