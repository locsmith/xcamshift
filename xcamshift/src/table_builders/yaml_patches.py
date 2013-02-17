#-------------------------------------------------------------------------------
# Copyright (c) 2013 Gary Thompson & The University of Leeds.
# All rights reserved. This program and the accompanying materials
# are made available under the terms of the GNU Lesser Public License v3.0
# which accompanies this distribution, and is available at
# http://www.gnu.org/licenses/lgpl-3.0.html
# 
# Contributors:
#     gary thompson - initial API and implementation
#-------------------------------------------------------------------------------
'''
Created on 7 Apr 2012

@author: garyt
'''
import yaml
from collection_backport import OrderedDict
from yaml.representer import Representer
import collection_backport
from yaml.constructor import ConstructorError
import collections
from yaml.dumper import Dumper
from python_utils import tupleit

# TODO hide functions which are not parts of public interfaces and add tests 
def add_access_to_yaml_list_based_keys(loader=None):
    if loader ==  None:
        return yaml.add_constructor(u'tag:yaml.org,2002:map', _construct_yaml_map)
    else:
        return yaml.add_constructor(u'tag:yaml.org,2002:map', _construct_yaml_map,loader)
    
def _construct_yaml_map(loader, node):
    pairs = [(tupleit(key) if isinstance(key, list) else key, value)
             for (key, value) in loader.construct_pairs(node, deep=True)]
    return dict(pairs)   

def float_representer(dumper, data):
    if abs(data) < 1e-9:
        result =  dumper.represent_scalar(u'tag:yaml.org,2002:null', u'~')
    else:
        result =  dumper.represent_scalar(u'tag:yaml.org,2002:float', u'%+11.8f' % data)
    return result

def apply_patch_float_format_with_nulls():
    yaml.add_representer(float,float_representer)
    
    
def construct_ordered_mapping(self, node, deep=False):
    if not isinstance(node, yaml.MappingNode):
        raise ConstructorError(None, None,
                "expected a mapping node, but found %s" % node.id,
                node.start_mark)
    mapping = collection_backport.OrderedDict()
    for key_node, value_node in node.value:
        key = self.construct_object(key_node, deep=deep)
        if not isinstance(key, collections.Hashable):
            raise ConstructorError("while constructing a mapping", node.start_mark,
                    "found unhashable key", key_node.start_mark)
        value = self.construct_object(value_node, deep=deep)
        mapping[key] = value
    return mapping



def construct_yaml_map_with_ordered_dict(self, node):
    data = collection_backport.OrderedDict()
    yield data
    value = self.construct_mapping(node)
    data.update(value)



def represent_ordered_mapping(self, tag, mapping, flow_style=None):
    value = []
    node = yaml.MappingNode(tag, value, flow_style=flow_style)
    if self.alias_key is not None:
        self.represented_objects[self.alias_key] = node
    best_style = True
    if hasattr(mapping, 'items'):
        mapping = list(mapping.items())
    for item_key, item_value in mapping:
        node_key = self.represent_data(item_key)
        node_value = self.represent_data(item_value)
        if not (isinstance(node_key, yaml.ScalarNode) and not node_key.style):
            best_style = False
        if not (isinstance(node_value, yaml.ScalarNode) and not node_value.style):
            best_style = False
        value.append((node_key, node_value))
    if flow_style is None:
        if self.default_flow_style is not None:
            node.flow_style = self.default_flow_style
        else:
            node.flow_style = best_style
    return node


    
def apply_ordered_dict_patch():
    yaml.add_representer(OrderedDict,Representer.represent_dict)
    
    
    yaml.representer.BaseRepresenter.represent_mapping = represent_ordered_mapping
    
    yaml.representer.Representer.add_representer(collection_backport.OrderedDict,
            yaml.representer.SafeRepresenter.represent_dict)
    

    
    
def tuple_representer(dumper, data):
    return dumper.represent_list(data)

def apply_tuple_patch():
    yaml.add_representer(tuple, tuple_representer)
    
def ignore_aliases(self,_data):
    return True

def apply_no_aliases_patch():
    Dumper.ignore_aliases = ignore_aliases
