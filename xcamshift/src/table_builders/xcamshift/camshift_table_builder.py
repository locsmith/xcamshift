'''
Created on 30 Mar 2012

@author: garyt




'''
import sys
import os
import re
import glob
from table_builders.xcamshift.Backbone_distance_extractor import BB_table_extractor
from table_builders.xcamshift.Extra_distance_extractor import XTRA_table_extractor
from table_builders.xcamshift.Sidechain_distance_extractor import SC_table_extractor
from table_builders.xcamshift.Dihdedral_distance_extractor import DIHEDRALS_table_extractor
from table_builders.xcamshift.Ring_table_extractor import RING_table_extractor
from table_builders.xcamshift.Nonbonded_table_extractor import Nonbonded_table_extractor
from common_constants import CAMSHIFT_SUB_POTENTIALS
import argparse
from yaml import load


FAILURE = -1
VERSION  = '1.0.0'

def get_atom_type_for_filename(file_path):
    file_name = os.path.split(file_path)[1]
    file_pattern = re.compile("1([A-Z]{1,2})orgpairs.par.*")
    match = re.search(file_pattern, file_name)
    return match.group(1)
    
        
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
        
        for residue_type in self.file_types:
            self._parse_files(residue_type)

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
    
    
    def _parse_files(self,residue_type):
        self.data.setdefault(residue_type,{})
        for file_path in self.files_by_type[residue_type]:
            atom_type = get_atom_type_for_filename(file_path)
            file_h = open(file_path)
            
            for line in file_h:
                
                line = line.strip()
                if line.startswith("#"):
                    continue
                
                items = line.split()
                
                table_type = items[0]
                
                values = [float(item) for item in items[1:]]
                self.data[residue_type].setdefault(table_type, {})[atom_type] = values
        

            
    def get_file_types(self):
        return self.file_types
    
    def get_number_files(self):
        result = 0
        for residue_type in self.get_file_types():
            result += len(self.files_by_type[residue_type])
            
        return result


def _get_file_type_name(residue_type):
    if residue_type == "":
        file_type_name = 'base'
    else:
        file_type_name = residue_type
    
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
    
    
    
def _build_output_name(sub_potential_name, residue_type, camshift_version, template='cams_%s_%s_%s.yaml', version_template = '%i_%i_%i'):
    version_string  = version_template % camshift_version
    
    sub_potential_name = sub_potential_name.lower().strip()
    
    if residue_type == '':
        residue_type='base'
    return template % (version_string, sub_potential_name,residue_type) 



def _build_args_list():
    message = 'Application to parse camshift data tables to yaml for xcamshift'
    
    parser = argparse.ArgumentParser(description=message,version='1.0')
    
    parser.add_argument('directory', action="store")
    parser.add_argument('--verbose', action="store_true", default=False)
    parser.add_argument('-o', '--output', default="-", help="where to put the files either a file system path or - for stdout") 
    parser.add_argument('-t', '--template', default='../xcamshift_templates', help='directory where template files are kept: %(default)s}')
    
    return  parser.parse_args()





def _make_output_directory_or_exit(output_dir):
    if not os.path.isdir(output_dir):
        if os.path.isfile(output_dir):
            print >> sys.stderr, 'output directory (%s) cannot be a file!' % output_dir
            print >> sys.stderr, 'exiting...'
            sys.exit(FAILURE)
        os.mkdir(output_dir)
    if not os.path.isdir(output_dir):
        print >> sys.stderr, 'cannot create output directory (%s)!' % output_dir
        print >> sys.stderr, 'exiting...'
        sys.exit(FAILURE)
        



def _write_file(output_data, output_path):
    file_handle = None
    try:
        file_handle = open(output_path, 'w')
        print >> file_handle, output_data
        file_handle.close()
    except:
        if file_handle != None:
            file_handle.close()
        else:
            print >> sys.stderr, "couldn't write to %s" % output_path
            print >> sys.stderr, 'exiting...'
            sys.exit(FAILURE)


def _get_human_readable_file_types(file_types):
    result = []
    for file_type in file_types:
        if file_type == '':
            result.append('base')
        else:
            result.append(file_type)
    return result       


    
def _read_version(table_dir):
    try:
        version_file = os.path.join(table_dir,'version.yaml')
        version_fp = open(version_file)
        version_info = load(version_fp)
        
    except Exception as detail:
        raise Exception("couldn't load version file %s" % version_file, detail)
    
    return version_info


def _get_template_filename(sub_potential, residue_type,version):
    if residue_type == '':
        residue_type = 'base'
    sub_potential_name = sub_potential.strip().lower()
    version_string = '%i_%i_%i' % version
    file_name = 'cams_%s_%s_%s_template.txt' % (version_string, sub_potential_name, residue_type)
    return file_name



def _build_dir_path(dir_path, input_path):
    if os.path.isabs(input_path):
        dir_path = [input_path]
    else:
        dir_path = [dir_path, input_path]
    return os.path.join(*dir_path)

def _build_filename_path(dir_path, relative_path, file_name):
    directory_path = _build_dir_path(dir_path, relative_path)
    
    template_path = [directory_path, file_name]
    
    file_path = os.path.join(*template_path)
    
    return file_path

def _read_template(dir_path,sub_potential,residue_type,version, path='../xcamshift_templates'):
    file_name = _get_template_filename(sub_potential, residue_type, version)
    
    file_path = _build_filename_path(dir_path, path, file_name)

    template = None
    if os.path.isfile(file_path):
        with  open(file_path) as template_file_h:
            template  = template_file_h.readlines()
            template  = ''.join(template)

    return template
             
        

if __name__ == '__main__':

    
    args = _build_args_list()
    
    table_dir = args.directory
    
    if args.verbose:
        print >> sys.stderr, 'camshift_table_builder, version %s' % VERSION
        print >> sys.stderr
        print >> sys.stderr, 'reading files from %s' % table_dir
    
    
    reader = _read_data_into_reader(table_dir)
    
    version_info = _read_version(table_dir)
    camshift_version = tuple(version_info['version'])
    camshift_version_string = '%i.%i.%i' % tuple(camshift_version)
    camshift_source = version_info['source']
    
    residue_types = reader.get_file_types()
    
    table_types = CAMSHIFT_SUB_POTENTIALS
    
    table_type_names = [table_type.strip().upper() for table_type in table_types]
    files_to_output = len(table_type_names) * len(residue_types)
    
    if args.verbose:
        print >> sys.stderr, '  read %i files' % reader.get_number_files()
        print >> sys.stderr, '  camshift version %s (source: %s)' % (camshift_version_string,camshift_source)
        print >> sys.stderr, '  residue types are % s' % ','.join(_get_human_readable_file_types(residue_types))
        print >> sys.stderr, '  sub potentials are %s' % ', '.join(table_type_names)
        print >> sys.stderr, '  output directory is %s' % _build_dir_path(args.directory, args.output)
        print >> sys.stderr, '  there are %i files to produce...' % (files_to_output)
        print >> sys.stderr, ''
    
    
    STDOUT = '-'
    if args.output != STDOUT:
        output_dir = _build_dir_path(args.directory, args.output)
        _make_output_directory_or_exit(output_dir)    
    
    #TODO use extract function here
    extractor_classes = _get_extractor_classes()

    count = 1
    for extractor_class, table_type in zip(extractor_classes,table_types):
        extractor = extractor_class(reader.data)
        
        for residue_type in residue_types:
            file_type_name  = _get_file_type_name(residue_type)
            
            sub_potential_name  = extractor.get_name()
            
            output_data = extractor.extract(residue_type)
            
            template = _read_template(table_dir, sub_potential_name, residue_type, camshift_version,path=args.template)
            template_filename = _get_template_filename(sub_potential_name, residue_type, camshift_version,)
            
            if template != None:
                output_data = template % output_data            

            if args.output == STDOUT:
                title = _build_output_name(sub_potential_name, residue_type, camshift_version, 
                                          template='camshift %s - %s - %s', version_template='%s.%s.%s')
                print '----- %s ------' % title
                print output_data
                
            else:
                output_filename = _build_output_name(sub_potential_name, residue_type, camshift_version)
                
                if args.verbose:
                    print >> sys.stderr, '%i of %i. %s' % (count, files_to_output,  output_filename), 
                    if template != None:
                        print >> sys.stderr, '(with template %s)'  % template_filename
                    else:
                        print >> sys.stderr
                    
                output_path = _build_filename_path(table_dir, args.output, output_filename)
                
                _write_file(output_data, output_path)
                count += 1
    if args.verbose:
        print >> sys.stderr, ''
        print >> sys.stderr, 'done.'
        
