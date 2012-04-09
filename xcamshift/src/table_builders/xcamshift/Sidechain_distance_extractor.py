'''
Created on 7 Apr 2012

@author: garyt
'''

from table_builders.common_constants import HA,CA,H,N,C,CB,O, h_keys
from ..yaml_patches import apply_ordered_dict_patch, apply_patch_float_format_with_nulls, \
     apply_tuple_patch
from collection_backport import OrderedDict
from table_builders.formatters import fixup_null_values,\
    fixup_decimal_spacing, fixup_convert_H_to_HN, fixup_replace_plus_with_space,\
    fixup_tuple_key_spacing, fixup_put_lonely_keys_on_new_line,\
    fixup_complex_key_question_mark,\
    global_fixup_colons_on_same_line_as_tuple_key,\
    global_fixup_data_dicts_on_same_line_as_key, fixup_data_key_spacing,\
    fixup_line_data_key_spacing
from table_builders.table_extractor import Table_extractor



RESIDUE_TYPES = ('ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLN', 'GLU', 'GLY', 'HIS', 'ILE', 'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'SER', 'THR', 'TRP', 'TYR', 'VAL')

RESIDUE_ATOM_TYPES = {
    'ALA'   :  ( 'CB', 'HB1', 'HB2', 'HB3'),
    'ARG'   :  ( 'CB', 'CG', 'CD', 'NE', 'CZ', 'NH1', 'NH2', 'HB1', 'HB2', 'HG1', 'HG2', 'HD1', 'HD2', 'HE', 'HH11', 'HH12', 'HH21', 'HH22'),
    'ASN'   :  ( 'CB', 'CG', 'OD1', 'ND2', 'HB1', 'HB2', 'HD21', 'HD22'),
    'ASP'   :  ( 'CB', 'CG', 'OD1', 'OD2', 'HB1', 'HB2'),
    'CYS'   :  ( 'CB', 'SG', 'HB1', 'HB2', 'HG1'),
    'GLN'   :  ( 'CB', 'CG', 'CD', 'OE1', 'NE2', 'HB1', 'HB2', 'HG1', 'HG2', 'HE21', 'HE22'),
    'GLU'   :  ( 'CB', 'CG', 'CD', 'OE1', 'OE2', 'HB1', 'HB2', 'HG1', 'HG2'),
    'GLY'   :  ( 'HA2',),
    'HIS'   :  ( 'CB', 'CG', 'ND1', 'CD2', 'CE1', 'NE2', 'HB1', 'HB2', 'HD1', 'HD2', 'HE1', 'HE2'),
    'ILE'   :  ( 'CB', 'CG1', 'CG2', 'CD', 'HB', 'HG11', 'HG12', 'HG21', 'HG22', 'HG23', 'HD1', 'HD2', 'HD3'),
    'LEU'   :  ( 'CB', 'CG', 'CD1', 'CD2', 'HB1', 'HB2', 'HG1', 'HD11', 'HD12', 'HD13', 'HD21', 'HD22', 'HD23'),
    'LYS'   :  ( 'CB', 'CG', 'CD', 'CE', 'NZ', 'HB1', 'HB2', 'HG1', 'HG2', 'HD1', 'HD2', 'HE1', 'HE2', 'HZ1', 'HZ2', 'HZ3'),
    'MET'   :  ( 'CB', 'CG', 'SD', 'CE', 'HB1', 'HB2', 'HG1', 'HG2', 'HE1', 'HE2', 'HE3'),
    'PHE'   :  ( 'CB', 'CG', 'CD1', 'CD2', 'CE1', 'CE2', 'CZ', 'HB1', 'HB2', 'HD1', 'HD2', 'HE1', 'HE2', 'HZ'),
    'PRO'   :  ( 'CB', 'CG', 'CD', 'HB1', 'HB2', 'HG1', 'HG2', 'HD1', 'HD2'),
    'SER'   :  ( 'CB', 'OG', 'HB1', 'HB2', 'HG1'),
    'THR'   :  ( 'CB', 'OG1', 'CG2', 'HB', 'HG1', 'HG21', 'HG22', 'HG23'),
    'TRP'   :  ( 'CB', 'CG', 'CD1', 'CD2', 'NE1', 'CE2', 'CE3', 'CZ2', 'CZ3', 'CH2', 'HB1', 'HB2', 'HD1', 'HE1', 'HE3', 'HZ2', 'HZ3', 'HH2'),
    'TYR'   :  ( 'CB', 'CG', 'CD1', 'CD2', 'CE1', 'CE2', 'CZ', 'OH', 'HB1', 'HB2', 'HD1', 'HD2', 'HE1', 'HE2', 'HH'),
    'VAL'   :  ( 'CB', 'CG1', 'CG2', 'HB', 'HG11', 'HG12', 'HG13', 'HG21', 'HG22', 'HG23'),
}

SC  = 'COSC'


class SC_table_extractor(Table_extractor):
    
    def __init__(self,data):
        self._data =  data
        
        apply_patch_float_format_with_nulls()
        apply_tuple_patch()
        apply_ordered_dict_patch()

    def serialize(self,data):
        out_data = OrderedDict()
        out_data['data'] =  OrderedDict()
        
        for residue_type in RESIDUE_TYPES:
            residue_data_label = SC + residue_type + "2"
            residue_data = data[residue_data_label]
            
            for i, atom_type in enumerate(RESIDUE_ATOM_TYPES[residue_type]):
                out_line = out_data['data'].setdefault(residue_type,OrderedDict()).setdefault(atom_type,OrderedDict())
                for h_key in h_keys:
                    value=0.0
                    if i < len(residue_data[h_key]):
                        value  =  residue_data[h_key][i]
                    out_line[h_key] = value


        return out_data
    

    def _get_data(self, file_type=''):
        return super(SC_table_extractor, self)._get_data(file_type)

    def format_lines(self,lines):
        result = []
        for line in lines.split('\n'):
            line = fixup_null_values(line)
            line = fixup_data_key_spacing(line)
            line = fixup_line_data_key_spacing(line)
            line = fixup_decimal_spacing(line)
            line = fixup_replace_plus_with_space(line)
            line = fixup_convert_H_to_HN(line)
            line = fixup_put_lonely_keys_on_new_line(line)
            result.append(line)
        return result
    
