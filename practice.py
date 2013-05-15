# first, open PyMOL and run PyMOLPyRosettaServer.py
"""
copyright Ching Yao Yang
Polytechnic of NYU
for my lovely PTE and PJB
"""
from rosetta import *
init()
# init of PyRosetta

pose = pose_from_sequence("AAAAAALVLSFFFWWAAHHHHHPPP");
# create pose from pdb

"""
use pose = pose_from_sequence()
or  pose = pose_from_rcsb()

"""

# pose.isfullatom()
# check full atom model


print pose.residue(5)
print pose.phi(5)
print pose.psi(5)
print pose.omega(5)
# outputs the 5th amino acid details
# Phe104 in PTE(out.pdb and 1HZY is residue(70)

print pose.sequence();
get_secstruct(pose);
# print out secondary structure

scorefxn = create_score_function("score12");
print scorefxn;
# create score function and display its measures

scorefxn.show(pose);
# create a table for the socre of current pose

mutation = Pose();
mutation.assign(pose);
# create a test for the pose

kT = 1.0;
n_moves = 5;
movemap = MoveMap();
# movemap.set_bb_true_range(9, 15)
movemap.set_bb(True);
small_mover = SmallMover(movemap, kT, n_moves);
shear_mover = ShearMover(movemap, kT, n_moves);
# this is for small mover, set up kT=1 and n_moves=1

min_mover = MinMover()
min_mover.movemap()
min_mover.score_function(scorefxn)
# set up minimization
"""minimization is little bit weird here?
small/shear move -> then minimization from side chain - > monte Carlo?!

"""

seq_mover = SequenceMover()
seq_mover.add_mover(small_mover)
seq_mover.add_mover(shear_mover)
seq_mover.add_mover(min_mover)

"""small_mover.apply(mutation);
shear_mover.apply(mutation);
"""
#simulate small and shear mover
repeat_mover = RepeatMover(seq_mover, 100);
# repeat mover

observer = PyMOL_Observer(mutation, True);

mc = MonteCarlo (mutation,scorefxn,kT);
# Monte Carlo simulation

trial_mover = TrialMover(small_mover, mc)
trial_mover.apply(mutation)
# set up tiral mover, combining mc and small movers

pmm = PyMOL_Mover();
pose.pdb_info().name("PTE");
mutation.pdb_info().name("mutation");
pmm.apply(pose);
pmm.apply(mutation);
pmm.keep_history(True);
# PyMOL Mover
"""This way you may be able to generate two models in PyMOL
one is original one
the other one is mutation from small_mover
"""

task_pack = standard_packer_task(pose);


generate_resfile_from_pose(pose,"pose.pdb")
