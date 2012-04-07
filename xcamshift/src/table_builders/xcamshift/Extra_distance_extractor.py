'''
Created on 7 Apr 2012

@author: garyt
'''

#xtra_v_keys = (   (H,0,HA,0),
#                  (H,0,C,0),
#                  (H,0,CB,0),
#                  
#                  (C,-1,HA,0),
#                  (C,-1,C,0),
#                  (C,-1,CB,0),
#                  
#                  (O,0,HA,0),
#                  (O,0,N,0),
#                  (O,0,CB,0),
#                  
#                  (N,1,HA,0),
#                  (N,1,HA,0),
#                  (N,1,HA,0),
#                  
#                  (O,-1,HA,-1),
#                  (O,-1,N,-1),
#                  (O,-1,CB,-1),
#                  
#                  (N,0,HA,-1),
#                  (N,0,N,-1),
#                  (N,0,CB,-1),
#                  
#                  (CG,0,HA,0),
#                  (CG,0,N,0),
#                  (CG,0,C,0),
#                  
#                  (CG,0,C,-1),
#                  (CG,0,N,+1),
#                  
#                  (CG,-1,CA,0),
#                  (CG,1,CA,0),
#                  
#                  (CA,-1,CA,1)
#)

#XTRA = 'XTRADISTS'
#class Extra_table_extractor(object):
#    
#    def __init__(self,data):
#        self._data =  data
#        
#        yaml.add_representer(float,BB_table_extractor.float_representer)
#        yaml.add_representer(OrderedDict,Representer.represent_dict)
#    
#    @staticmethod
#    def float_representer(dumper, data):
#        if abs(data) < 1e-9:
#            result =  dumper.represent_scalar(u'tag:yaml.org,2002:null', u'~')
#        else:
#            result =  dumper.represent_scalar(u'tag:yaml.org,2002:float', u'%+11.8f' % data)
#        return result
#    
#    
#    def extract(self, file_type  = ''):
#        
#        data =  self._data[file_type]
#        
#        serialized_data = self.serialize(data)
#        
##        print serialized_data
#        
#        lines = self.build_output_lines(serialized_data)
#        
##        print lines
##        
#        lines =  self.format_lines(lines)
##        
#        return "\n".join(lines)
#        
#    def serialize(self,data):
#        out_data = OrderedDict()
#        
##        for i,v_key in enumerate(bb_v_keys):
##            if v_key in missing_v_keys:
##                for h_key in h_keys:
##                    data[XTRA][h_key].insert(i, 0.0) 
#                
#                    
#        for i,v_key in enumerate(xtra_v_keys):
#    
##            print 'key',v_key[0:2],v_key[2:4]
#            out_line = out_data.setdefault(v_key[0:2],OrderedDict()).setdefault(v_key[2:4],OrderedDict())
#    
#    #        print v_key,
#    
#            for j,h_key in enumerate(h_keys):
##                print h_key,data[XTRA][h_key][i]
#                out_line[h_key] = data[XTRA][h_key][i]
##        print out_data[('H',0)][HA,0][HA]
##        sys.exit(-1)
#        print out_data
#        return out_data
#    
#    def build_output_lines(self,serialized_data):
#        return  yaml.dump(serialized_data,default_flow_style=None,width=1000,indent=6)
#    
#    @staticmethod
#    def fixup_key_spacing_callback(m):
#        values  =  m.group(1).split(',')
#        value_1 =  ("%s," % values[0]).ljust(3)
#        value_2 = ("%s" % values[1]).rjust(3)
#
#        return '[%s%s]' % (value_1,value_2)
#    
#    @staticmethod
#    def fixup_decimal_spacing_callback(m):
#        reduce_spaces_by  = len(m.group(3)) - 1
#        num_spaces =  len(m.group(1))
#        return " " * (num_spaces - reduce_spaces_by) +  m.group(2) +  m.group(3) + "."
#    
#    def format_lines(self,lines):
#        result = []
#        lines = re.sub("\]\n(\s*):","]:\n\g<1> ",lines)
#        lines = re.sub("\]:\n\s+{","] : {",lines)
#        for line in lines.split('\n'):
#            line = line.replace("? [","[")
#            line =  line.replace('~,'," ,%10s" % "")
#            line =  line.replace('~}'," ,%9s}" % "")
#            line = re.sub("\[(.*)\]",Extra_table_extractor.fixup_key_spacing_callback,line)
#            line = re.sub("(\s+)([+-])([0-9]+)\.",BB_table_extractor.fixup_decimal_spacing_callback,line)
##            line = re.sub("^([0-9])","\n\g<0>",line)
#            line = re.sub("\+"," ",line)
#            line = re.sub("H ","HN",line)
#            line = re.sub("^\[","\n[",line)
#            result.append(line)
#        return result
#    
