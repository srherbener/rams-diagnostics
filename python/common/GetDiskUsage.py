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

UserAliases = {
    "areyh" : "aryeh",
    "adele" : "aigel",
    "michalcl" : "clavner",
    }

print("Looking at group disk usage:")

# Grab the group totals information, and load into a "soup" database
NetworkDriveUrl = "http://salix.atmos.colostate.edu/serverStatus/networkDrives.html"

print("  Reading URL: {0:s}".format(NetworkDriveUrl))
Ndrive = req.get(NetworkDriveUrl)

print("")

NdriveSoup = BeautifulSoup(Ndrive.content, "lxml")

# Find user usage on each server
# Servers are in sections starting with:
#    <a name="server_name"> </a>
# The usage for each server is located in a <table></table>
# structure following the <a></a> with the server name
#
UsageByServer = {}
UserList = []
TotName = "TotalAmount"
for Server in ServerList:
    UsageByServer[Server] = {}

    print("  Checking usage on server: {0:s}".format(Server))
    svr = NdriveSoup.find("a", {"name":Server})

    tbl = svr.find_next_sibling("table")

    for tbl_item in tbl.find_all("td"):
        user_list = tbl_item.contents
        # skip over empty lists
        if (user_list):
            # split the usage info string into components
            user_list = user_list[0].replace("\n", "").split()

            # Pull off the user name from the path in the first component
            # and amount of disk usage in the second component.
            User = os.path.basename(user_list[0])
            if (User in UserAliases):
                User = UserAliases[User]

            Amount = float(user_list[1])
            if (user_list[2] == "Gb"):
                Amount = Amount / 1024.0  # convert to TB

            # Store the server, user, amount
            UserList.append(User)
            UsageByServer[Server][User] = Amount

UserList = set(UserList) # remove duplicates
print("")

# Reformat the usage by server data into usage by user.
# Sum up the total amount of each user and add that entry.
#
# A few users have different user names on different servers. For these
# users, add together the usage for each server, and place into just
# one amount under a single user name.
UsageByUser = {}
GroupTotal = 0.0
for User in UserList:
    UsageByUser[User] = {}
    UserTotal = 0.0

    for Server in ServerList:
        UsageByUser[User][Server] = 0.0

        if (User in UsageByServer[Server]):
            Amount = UsageByServer[Server][User]
            UsageByUser[User][Server] = Amount
            UserTotal = UserTotal + Amount

    UsageByUser[User][TotName] = UserTotal
    GroupTotal = GroupTotal + UserTotal

# Create a list of users for the plot.
#   Select users consuming >= 1.0 TB of file space
#   Sort the list of users into low to high usage
#
# In the sorted command:
#    the .items() method is creating a list of tuples of the form (key, value)
#       from the dictionary TopUsers
#    the key= arguemnt is telling sort to use the second item in the tuple
#       (ie, the total disk usage) to do the sorting.
MinUsage = 1.0
TopUsers = {}
for User in UserList:
    UserUsage = UsageByUser[User][TotName]
    if (UserUsage >= MinUsage):
        TopUsers[User] = UserUsage

SortedTopUsers = (sorted(TopUsers.items(), key=lambda t: t[1]))
N = len(SortedTopUsers)

TopLabels = np.zeros(N, dtype='<U10')
TopUsage = np.zeros(N)

for i in range(N):
    TopLabels[i] = SortedTopUsers[i][0]
    TopUsage[i]  = SortedTopUsers[i][1]

# Create a horizontal bar plot showing the heaviest users first
Width = 0.8
Ind = np.arange(N)
Title = "Top Disk Usage by User\n(Total Usage = {0:.2f} Tb)".format(GroupTotal)
Color = 'skyblue'

Fig = plt.figure()

plt.barh(Ind, TopUsage, Width, color=Color, edgecolor=Color, align='center')

plt.title(Title)
plt.xlabel('Usage (Tb)')
plt.yticks(Ind, TopLabels)

OutFile = "{0:s}/tmp/DiskUsageByUser.png".format(os.environ['HOME'])
print("  Writing: {0:s}".format(OutFile))
Fig.savefig(OutFile)
