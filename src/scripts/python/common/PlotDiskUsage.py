#!/usr/bin/env python3

import os
import sys
sys.path.append("{0:s}/etc/python/common".format(os.environ['HOME']))

import numpy as np
import h5py
import matplotlib.pyplot as plt

ScriptName = sys.argv[0]

if (len(sys.argv) != 3):
    print("{0:s}: ERROR: Must supply exactly two arguments".format(ScriptName))
    print("")
    print("    USAGE: {0:s} <Type> <Name>".format(ScriptName))
    print("    USAGE:     <Type> is one of: 'Server' or 'User'")
    print("    USAGE:     <Name> is the name of the type selected by <Type>")
    sys.exit(1)

Type = sys.argv[1]
Name = sys.argv[2]

# Check the arguments
#   Type must be one of "Server", "User"
if (Type == "Server"):
    GroupName = "/Servers/{0:s}".format(Name)
elif (Type == "User"):
    GroupName = "/Users/{0:s}".format(Name)
else:
    print("{0:s}: ERROR: <Type> argument must be one of: 'Server' or 'User'".format(ScriptName))
    sys.exit(2)

print("Creating pie chart for {0:s} named {1:s}".format(Type, Name))

InFname = "{0:s}/tmp/DiskUsage.h5".format(os.environ['HOME'])

InFile = h5py.File(InFname, "r")

TimeStamp = InFile['TimeStamp'][...]
Units = InFile['Units'][...]

# Create a dictionary for looking up plot colors
# If plotting a server, you need the user colors
# If plotting a user, you need the server colors
if (Type == "Server"):
    TopGroup = InFile['/Users']
elif (Type == "User"):
    TopGroup = InFile['/Servers']

ColorMap = { }
for (ObjName, Object) in TopGroup.items():
    if isinstance(Object, h5py.Group):
        ColorMap[ObjName] = Object.attrs['Pcolor']

# Walk through the group that is specified by the Name argument
# grabbing the element names and amounts for the plot.
Group = InFile[GroupName]
UsageList = { }
for (ObjName, Object) in Group.items():
    if isinstance(Object, h5py.Dataset):
        # Object is a dataset. Use "ObjName" instead of Object.name since
        # the former contains the basename part of the path and
        # the latter contains the full path from the top of the tree.
        UsageList[ObjName] = Object[...]

# Check the attributes for GroupName.
# A server will have "Size", "Used", "Free" attributes
# A user will not have any of these attributes.
# So, initialize the Grp* varialbes to the user case (size = 15 TB).
GrpSize = 15.0
GrpUsed = 0.0
GrpFree = 0.0
for (AttrName, AttrVal) in Group.attrs.items():
    if (AttrName == 'Size'):
        GrpSize = AttrVal
    elif (AttrName == 'Used'):
        GrpUsed = AttrVal
    elif (AttrName == 'Free'):
        GrpFree = AttrVal

InFile.close()

# Sort the usage list in descending order of the amounts
# In the sorted command:
#    the .items() method is creating a list of tuples of the form (key, value)
#       from the dictionary UsageList
#    the key= arguemnt is telling sort to use the second item in the tuple
#       (ie, the disk usage) to do the sorting.
SortedUsageList = sorted(UsageList.items(), key=lambda x : x[1], reverse=True)

# Split off the data in the SortedUsageList
Labels = [ ]
Values = [ ]
PctValues = [ ]
Colors = [ ]
for i in range(len(SortedUsageList)):
    Lab = SortedUsageList[i][0]
    Val = SortedUsageList[i][1]
    PctVal = (Val / GrpSize) * 100.0
    if (PctVal >= 0.1):
        Labels.append(Lab)
        Values.append(Val)
        PctValues.append(PctVal)
        Colors.append(ColorMap[Lab])

# Put the free space, if non-zero at end of lists
# Use white for the free space
if (GrpFree > 0):
    Labels.append('Free Space')
    Values.append(GrpFree)
    PctValues.append((GrpFree/GrpSize)*100.0)
    Colors.append('white')

# Walk through the above arrays and add the absolute values and percentages
# to the labels
for i in range(len(Labels)):
    Lab = Labels[i]
    Val = Values[i]
    PctVal = PctValues[i]
    Labels[i] = "{0:s} ({1:.1f} {2:s}, {3:.1f} %)".format(Lab, float(Val), str(Units), float(PctVal))

if (Type == "Server"):
    Ptitle = "Disk Usage: {0:s}\nServer: {1:s} (Size = {2:.1f} {3:s})".format(
              str(TimeStamp), str(Name), float(GrpSize), str(Units))
elif (Type == "User"):
    Ptitle = "Disk Usage: {0:s}\nUser: {1:s} (Allot. = {2:.1f} {3:s})".format(
              str(TimeStamp), str(Name), float(GrpSize), str(Units))

Fig = plt.figure()

# Need to get pie chart and legend to fit without overlapping.
#   Place the legend on the right side. Specify an anchor point that is in the center of
#   the right side of the plot frame (bbox_to_anchor=(1,0.5) ). The bbox_transform
#   setting is to make it so that bbox_to_anchor arguement is in figure
#   coordinates (x and y go from 0 to 1).
#
#   Place the pie chart on the left side. This is done with the subplots_adjust method.
Patches, Text = plt.pie(Values, colors=Colors)
plt.legend(Patches, Labels, bbox_to_anchor=(1,0.5), loc="center right", fontsize=12, 
           bbox_transform=plt.gcf().transFigure)
plt.subplots_adjust(left=0.1, bottom=0.1, right=0.5)

plt.axis('equal') # make the pie a circle
plt.title(Ptitle)

OutFname = "{0:s}/tmp/DiskUsage_{1:s}_{2:s}.png".format(os.environ['HOME'], Type, Name)
print("  Writing: {0:s}".format(OutFname))
Fig.savefig(OutFname)

