#!/usr/bin/env python3

import os
import sys
sys.path.append("{0:s}/etc/python/common".format(os.environ['HOME']))

import re

import requests as req
from bs4 import BeautifulSoup
import numpy as np
from matplotlib import pyplot as plt

ServerList = [
    "tasman",
    "snow-home",
    "blizzard",
    "frost-home",
    "dendrite-home",
    "ccn-home",
    "icehome",
    "dendriteuser1",
    "dendriteuser2",
    "dendriteuser3",
    ]

print("Looking at group disk usage:")

# Grab the group totals information, and load into a "soup" database
GroupTotalsUrl = "http://salix.atmos.colostate.edu/serverStatus/sueGroup.html"
NetworkDriveUrl = "http://salix.atmos.colostate.edu/serverStatus/networkDrives.html"

print("  Reading URL: {0:s}".format(GroupTotalsUrl))
Gtot = req.get(GroupTotalsUrl)

print("  Reading URL: {0:s}".format(NetworkDriveUrl))
Ndrive = req.get(NetworkDriveUrl)

GtotSoup = BeautifulSoup(Gtot.content, "lxml")
NdriveSoup = BeautifulSoup(Ndrive.content, "lxml")

# All of the individual user usage data is formatted as:
#
#   <b> UserName </b> Amount Units <br/>
#
# So, look for all the "b" tags plus their next siblings
UsageList = {}
for tag in GtotSoup.find_all("b"):
    # Continue the loop until the user name is "Top Three users"
    # (which is a summary listing the heaviest disk usage by user)
    UserName = tag.string.replace("\n", "")
    if (UserName == "Top Three users"):
        break

    # Grab the Amount and Units. Convert all amounts to Tb.
    UserUsage = tag.next_sibling.replace("\n", "")
    UserUsage = UserUsage.split()

    Amount = float(UserUsage[0])
    if (UserUsage[1] == "Gb"):
        Amount = Amount / 1024.0

    UsageList[UserName] = Amount

# Find user usage on each server
# Servers are in sections starting with:
#    <a name="server_name"> </a>
#
for server in ServerList:
    print("Checking usage on server: {0:s}".format(server))
    svr = NdriveSoup.find("a", {"name":server})

    tbl = svr.find_next_sibling("table")

    for tbl_item in tbl.find_all("td"):
        user_list = tbl_item.contents
        # skip over empty lists
        if (user_list):
            # split the usage info string into components
            user_list = user_list[0].replace("\n", "").split()
            # pull off the user name from the path in the first component
            user_list[0] = os.path.basename(user_list[0])
            print(user_list)

SortedUsageList = (sorted(UsageList.items(), key=lambda t: t[1]))

# In this case, the sorted function will create a list of tuples where each tuple is
# a key:value pair from the UsageList dictionary. Peel these off into a text list (keys)
# and a numpy array (values) in order to facilitate plotting.
Nitems = len(SortedUsageList)

Labels = np.zeros(Nitems, dtype='<U10')
DiskUsage = np.zeros(Nitems)

for i in range(Nitems):
    Labels[i] = SortedUsageList[i][0]
    DiskUsage[i] = SortedUsageList[i][1]

TotalUsage = DiskUsage.sum()

# For the plot, chop off the small usage values. This will help de-clutter
# the plot.
MinUsage = 1.0 # Tb
TopSelect = (DiskUsage > MinUsage)
TopUsage = DiskUsage[TopSelect]
TopLabels = Labels[TopSelect]
N = len(TopUsage)

# Create a horizontal bar plot showing the heaviest users first
Width = 0.8
Ind = np.arange(N)
Title = "Top Disk Usage by User\n(Total Usage = {0:.2f} Tb)".format(TotalUsage)
Color = 'skyblue'

Fig = plt.figure()

plt.barh(Ind, TopUsage, Width, color=Color, edgecolor=Color, align='center')

plt.title(Title)
plt.xlabel('Usage (Tb)')
plt.yticks(Ind, TopLabels)

OutFile = "{0:s}/tmp/DiskUsageByUser.png".format(os.environ['HOME'])
print("  Writing: {0:s}".format(OutFile))
Fig.savefig(OutFile)
