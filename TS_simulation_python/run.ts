# PREREQUISITES:
#
#     You must use moltemplate.sh to create 3 files:
#        system.data  system.in.init  system.in.settings
#     (Follow the instructions in README_setup.sh, 
#      or run the file as a script using ./README_setup.sh)

# ------------------------------- Initialization Section --------------------

 units real
 atom_style full

 read_data       "system.data_modified"

include         "system.in.init"

# ------------------------------- Atom Definition Section -------------------

#read_data       "system.data"

# ------------------------------- Settings Section --------------------------

#include         "system.in.settings"

# ------------------------------- Run Section -------------------------------


#displace_atoms all move 100 100 100 units box

# -- minimization protocol --

variable watMoleMass equal 18.0153 # /(g/mol)
variable nAvog equal 6.0221415e23 # Avogadro's number
variable watMoleculeMass equal (${watMoleMass}/${nAvog}) # /(g/molecule)
variable A3_in_cm3 equal 1e-24 # Angstrom^3 in cm^3
variable nAtoms equal atoms
variable nMolecules equal v_nAtoms/3

variable    T equal 298.15
variable    dt equal 1

#timestep ${dt}

#neighbor 2.0 nsq
neigh_modify delay 0 every 1 check yes
replicate 1 1 1 

velocity all create $T 1234 rot no mom yes dist gaussian

timestep        0.5

#thermo 20000
thermo          100
thermo_style    custom step temp pe ke lx ly lz density press pxx pyy pzz
dump 2 all xyz  1000 dump.xyz
dump 3 all xyz  2000 dump_2000.xyz
dump 4 all xyz  5000 dump_5000.xyz
dump 5 all xyz  20000 dump_20000.xyz
dump_modify 2 element C H N O S P Cl F
dump_modify 3 element C H N O S P Cl F
dump_modify 4 element C H N O S P Cl F
dump_modify 5 element C H N O S P Cl F

#replicate 2 2 2

#fix             qeq all qeq/reax 1 0.0 10 1e-6 reax/c

#fix 202 all deform 1 y final 0 25 z final 0 25 units box 
#fix 3333 all nvt temp 300.0 300.0 100
#run 10000
#unfix 3333
#unfix 202

#fix NPT all npt temp 2000 2000 1000 iso 0.0 0.0 1000 drag 1
#run 20000
#unfix NPT

#fix NPT all npt temp $T $T 1000 iso 0.0 0.0 1000 drag 1
#run 20000
#unfix NPT

fix             2020 all qeq/reax 1 0.0 10 1e-6 reax/c

#fix NVT all nvt temp $T $T 1000
#run 10000

#unfix NVT

thermo 2000

variable tmp equal "lx"
variable L0 equal ${tmp}
print "Initial Length, L0: ${L0}"

######################################
# DEFORMATION
reset_timestep	0
#
fix		1 all npt temp 300 300 1000 y 0 0 1000 z 0 0 1000 drag 1
variable srate equal 1.0e10/2
variable srate1 equal "v_srate / 1.0e15"
fix		2 all deform 2000 x erate ${srate1} units box remap x
#
# Output strain and stress info to file
# for units metal, pressure is in [bars] = 100 [kPa] = 1/10000 [GPa]
# p2, p3, p4 are in GPa
variable strain equal "(lx - v_L0)/v_L0"
variable p1 equal "v_strain"
variable p2 equal "-pxx/10000"
variable p3 equal "-pyy/10000"
variable p4 equal "-pzz/10000"
fix def1 all print 2000 "${p1} ${p2} ${p3} ${p4}" file Al_SC_100.def1.txt screen no

fix average all ave/time 1 401 2000 v_p1 v_p2 v_p3 v_p4 file ts.profile

#
# Use cfg for AtomEye
#dump 		1 all cfg 250 dump.tensile_*.cfg mass type xs ys zs c_csym c_peratom fx fy fz
#dump_modify 1 element Al
#
# Display thermo
#thermo 	1000
thermo_style	custom step v_strain temp f_average[1] f_average[2] f_average[3] f_average[4] ke pe press vol
#
run		2000000
#
#######################################
# SIMULATION DONE
print "All done"
