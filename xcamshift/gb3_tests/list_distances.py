from pymol import cmd
atoms = cmd.get_model("name HN")
for atom_1 in atoms.atom:
    for atom_2 in atoms.atom:
       id_1 =  atom_1.id
       id_2 =  atom_2.id
       name_1 = atom_1.name
       name_2 = atom_2.name
       resid_1 = atom_1.resi
       resid_2 = atom_2.resi
       distance = cmd.distance( "id %s" % id_1 , "id %s" % id_2)
       if distance != 0.0 and  distance <= 5.0:
           print "assign (resid %4i and name %3s) (resid %4i and name %3s) %4.3f  %4.3f %4.3f\n" % (int(resid_1),name_1,int(resid_2),name_2 ,distance,1.0,1.0)
