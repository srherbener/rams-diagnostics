#!/usr/bin/env python3

import os
import sys
sys.path.append("{0:s}/etc/python/common".format(os.environ['HOME']))

import f90nml  # namelist parser

############### SUBROUTINES #######################
def FindListDiffs(List1, List2):
    # Routine to find common items in the two lists
    # and items unique to each list.
    # Use sets and instersection, etc. operators)

    Set1 = set(List1)
    Set2 = set(List2)
    Comm = Set1 & Set2
    Only1 = Set1 - (Comm)
    Only2 = Set2 - (Comm)

    # Sort each resulting list
    Only1Sorted = sorted(list(Only1))
    Only2Sorted = sorted(list(Only2))
    CommSorted = sorted(list(Comm))

    # return three lists
    return Only1Sorted, Only2Sorted, CommSorted

def ReportGroupVarDiffs(Group, Nml1, Nml2):
    # Routine to find common variables in the given
    # group for both namelists, plus the variable unique
    # to each namelist group.

    print("Variable check for group: {0:s}:".format(Group))
    Vars1 = list(Nml1[Group].keys())
    Vars2 = list(Nml2[Group].keys())

    # Find common variables
    Only1Vars, Only2Vars, CommVars = FindListDiffs(Vars1, Vars2)

    print("  Variables unique to File1:")
    for Var in Only1Vars:
        print("    {0:s}".format(Var))
    print("")

    print("  Variables unique to File2:")
    for Var in Only2Vars:
        print("    {0:s}".format(Var))
    print("")

    # while listing common vars, also report the differences
    print("  Variables common to both files:")
    for Var in CommVars:
        Val1 = Nml1[Group][Var]
        Val2 = Nml2[Group][Var]

        print("    {0:s} --> ".format(Var), end="")
        if (Val1 != Val2):
            print("DIFF, values: ", end="")
            print(Val1, Val2)
        else:
            VarStatus = "match"
            print("SAME, value: ", end="")
            print(Val1)

    print("")

    return



############### MAIN #######################
ScriptName = os.path.basename(sys.argv[0])

# Need two arguments which are the two files to be compared
if (len(sys.argv) != 3):
    print("{0:s}: ERROR: Must supply exactly two arguments".format(ScriptName))
    print("")
    print("    USAGE: {0:s} <File1> <File2>".format(ScriptName))
    print("    USAGE:     <File1> and <File2> are the paths to the two RAMSIN files being compared")
    sys.exit(1)

RamsFname1 = sys.argv[1]
RamsFname2 = sys.argv[2]

print("Comparing two RAMSIN files:")
print("  File1: {0:s}".format(RamsFname1))
print("  File2: {0:s}".format(RamsFname2))
print("")

# Apply the parser to the two files
Rnml1 = f90nml.read(RamsFname1)
Rnml2 = f90nml.read(RamsFname2)

# Data is stored in python dictionary form. The first level of keys are
# the namelist sections (model_grid, model_file_info, etc.), and the
# second level of keys are the varibles within each section.

# Get the group names
Rnml1Groups = list(Rnml1.keys())
Rnml2Groups = list(Rnml2.keys())

# Find the common groups and those unique to each file groups
Only1Groups, Only2Groups, CommGroups = FindListDiffs(Rnml1Groups, Rnml2Groups)

print("Namelist groups unique to File1:")
for Group in Only1Groups:
    print("  {0:s}".format(Group))
print("")

print("Namelist groups unique to File2:")
for Group in Only2Groups:
    print("  {0:s}".format(Group))
print("")

print("Namelist groups common to both files:")
for Group in CommGroups:
    print("  {0:s}".format(Group))
print("")

# For each of the common groups, repeat the FindListDiffs process on the lists
# of variables. Then for the variables in common compare their values and report
# the ones that differ.
for Group in CommGroups:
    ReportGroupVarDiffs(Group, Rnml1, Rnml2)

