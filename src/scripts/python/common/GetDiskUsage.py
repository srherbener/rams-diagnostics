#!/usr/bin/env python3

import os
import sys
sys.path.append("{0:s}/etc/python/common".format(os.environ['HOME']))

import requests as req
from bs4 import BeautifulSoup
import numpy as np
import h5py

import time

ServerList = {
    # name : [ alias from web page, size in TB ]
    "Tasman"           : [ "tasman", 30.0 ],
    "Avalanche"        : [ "snow-home", 34.0 ],
    "Blizzard"         : [ "blizzard", 33.0 ],
    "PermaFrost"       : [ "permafrost", 18.0 ],
    "DendriteHome"     : [ "dendrite-home", 1.8 ],
    "Cloudseed"        : [ "ccn-home", 13.0 ],
    "Icicle"           : [ "icicle", 28 ],
    "DendriteUser1"    : [ "dendriteuser1", 3.6 ],
    "DendriteUser2"    : [ "dendriteuser2", 3.6 ],
    "DendriteUser3"    : [ "dendriteuser3", 3.6 ],
    "vandenHeeverHome" : [ "vandenHeeverHome", 15.0 ],
    }

UserList = {
    # Name, [ list of aliases from web page ]
    "Minnie"  : [ "jpark" ],
    "Jennie"  : [ "jbukowski" ],
    "Leah"    : [ "ldgrant" ],
    "Peter"   : [ "pmarin" ],
    "Aryeh"   : [ "aryeh", "areyh" ],
    "Sean"    : [ "sfreeman" ],
    "Ben"     : [ "btoms" ],
    "Stacey"  : [ "skawecki" ],
    "SteveS"  : [ "smsaleeb" ],
    "SteveH"  : [ "sherbener" ],
    "Sue"     : [ "sue" ],

    "Adele"   : [ "adele", "aigel" ],
    "Matt"    : [ "mattigel" ],
    "Rachel"  : [ "storer" ],
    "Clay"    : [ "cjmcgee" ],
    "Amanda"  : [ "asheffie" ],

    "Other"   : [ ],
    }

# Create a map going from each of the aliases in the UserList to the
# primary name
Web2UserMap = { }
for User in UserList:
    for Alias in UserList[User]:
        Web2UserMap[Alias] = User

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
for Server in ServerList:
    ServerWebName = ServerList[Server][0]
    ServerSize    = ServerList[Server][1]

    UsageByServer[Server] = {}
    ServerUsed = 0.0

    print("  Checking usage on server: {0:s}".format(Server))
    svr = NdriveSoup.find("a", {"name":ServerWebName})

    tbl = svr.find_next_sibling("table")

    for tbl_item in tbl.find_all("td"):
        user_list = tbl_item.contents
        # skip over empty lists
        if (user_list):
            # split the usage info string into components
            user_list = user_list[0].replace("\n", "").split()

            # Pull off the user name from the path in the first component
            # and amount of disk usage in the second component.
            WebUser = os.path.basename(user_list[0])
            if (WebUser in Web2UserMap):
                User = Web2UserMap[WebUser]
            else:
                User = "Other"

            Amount = float(user_list[1])
            if (user_list[2] == "Gb"):
                Amount = Amount / 1024.0  # convert to TB

            # Store the server, user, amount
            if (User in UsageByServer[Server]):
                UsageByServer[Server][User] += Amount
            else:
                UsageByServer[Server][User] = Amount

            ServerUsed += Amount

    ServerFree = np.amax([ (ServerSize - ServerUsed), 0.0 ])
    UsageByServer[Server]['SIZE'] = ServerSize
    UsageByServer[Server]['USED'] = ServerUsed
    UsageByServer[Server]['FREE'] = ServerFree

print("")

# Reformat the usage by server data into usage by user.
# Sum up the total amount of each user and add that entry.
#
# A few users have different user names on different servers. For these
# users, add together the usage for each server, and place into just
# one amount under a single user name.
UsageByUser = {}
for User in UserList:
    UsageByUser[User] = {}

    for Server in ServerList:
        UsageByUser[User][Server] = 0.0

        if (User in UsageByServer[Server]):
            Amount = UsageByServer[Server][User]
            UsageByUser[User][Server] = Amount

# Create a timestamp
TimeStamp = time.ctime()

# Write out data into a single hdf5 file
OutFname = "{0:s}/tmp/DiskUsage.h5".format(os.environ['HOME'])
print("  Writing: {0:s}".format(OutFname))

OutFile = h5py.File(OutFname, mode='w')

OutFile['TimeStamp'] = TimeStamp
OutFile['Units'] = "TB"

for Server in UsageByServer:
    for User in UsageByServer[Server]:
        DsetName = "/Servers/{0:s}/{1:s}".format(Server, User)
        OutFile[DsetName] = UsageByServer[Server][User]

for User in UsageByUser:
    for Server in UsageByUser[User]:
        DsetName = "/Users/{0:s}/{1:s}".format(User, Server)
        OutFile[DsetName] = UsageByUser[User][Server]

OutFile.close()
