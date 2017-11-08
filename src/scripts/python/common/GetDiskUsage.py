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
    # name : [ alias from web page, size in TB, color for plotting ]
    "Tasman"     : [ "tasman",            30.0, "magenta"    ],
    "Avalanche"  : [ "snow-home",         34.0, "blue"       ],
    "Blizzard"   : [ "blizzard",          33.0, "green"      ],
    "PermaFrost" : [ "permafrost",        18.0, "gold"       ],
    "DendHome"   : [ "dendrite-home",      1.8, "red"        ],
    "CloudSeed"  : [ "ccn-home",          13.0, "violet"     ],
    "Icicle"     : [ "icicle",            28.0, "skyblue"    ],
    "DendUser1"  : [ "dendriteuser1",      3.6, "aquamarine" ],
    "DendUser2"  : [ "dendriteuser2",      3.6, "darkorange" ],
    "DendUser3"  : [ "dendriteuser3",      3.6, "brown"      ],
    "SharedHome" : [ "vandenHeeverHome",  15.0, "cyan"       ],
    "Squall"     : [ "squall",           174.0, "darkblue"   ],
    }

UserList = {
    # Name, [ [ list of aliases from web page ], alloted disk space, plot color ]
    "Minnie"  : [ [ "jpark" ],          15.0, "magenta"     ],
    "Jennie"  : [ [ "jbukowski" ],      15.0, "blue"        ],
    "Leah"    : [ [ "ldgrant" ],        15.0, "green"       ],
    "Peter"   : [ [ "pmarin" ],         15.0, "gold"        ],
    "Aryeh"   : [ [ "aryeh", "areyh" ], 15.0, "red"         ],
    "Sean"    : [ [ "sfreeman" ],       15.0, "violet"      ],
    "Ben"     : [ [ "btoms" ],          15.0, "skyblue"     ],
    "Stacey"  : [ [ "skawecki" ],       15.0, "aquamarine"  ],
    "SteveS"  : [ [ "smsaleeb" ],       15.0, "darkorange"  ],
    "SteveH"  : [ [ "sherbener" ],      15.0, "brown"       ],
    "Sue"     : [ [ "sue" ],            15.0, "cyan"        ],

    "Adele"   : [ [ "adele", "aigel" ], 15.0, "orange"      ],
    "Matt"    : [ [ "mattigel" ],       15.0, "forestgreen" ],
    "Rachel"  : [ [ "storer" ],         15.0, "yellow"      ],
    "Clay"    : [ [ "cjmcgee" ],        15.0, "purple"      ],
    "Amanda"  : [ [ "asheffie" ],       15.0, "darkblue"    ],

    "Other"   : [ [ ],                  15.0, "black"       ],
    }

# Create a map going from each of the aliases in the UserList to the
# primary name
Web2UserMap = { }
for User in UserList:
    AliasList = UserList[User][0]
    for Alias in AliasList:
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
    ServerPcolor  = ServerList[Server][2]

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
    UsageByServer[Server]['PCOLOR'] = ServerPcolor

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

    UserSize   = UserList[User][1]
    UserPcolor = UserList[User][2]

    UserUsed = 0.0

    for Server in ServerList:
        UsageByUser[User][Server] = 0.0

        if (User in UsageByServer[Server]):
            Amount = UsageByServer[Server][User]
            UsageByUser[User][Server] = Amount
            UserUsed += Amount

    UserFree = np.max([ (UserSize - UserUsed), 0.0 ])
    UsageByUser[User]['SIZE'] = UserSize
    UsageByUser[User]['USED'] = UserUsed
    UsageByUser[User]['FREE'] = UserFree
    UsageByUser[User]['PCOLOR'] = UserPcolor

# Create a timestamp
TimeStamp = time.ctime()

# Write out data into a single hdf5 file
OutFname = "{0:s}/tmp/DiskUsage.h5".format(os.environ['HOME'])
print("  Writing: {0:s}".format(OutFname))

OutFile = h5py.File(OutFname, mode='w')

OutFile['TimeStamp'] = TimeStamp
OutFile['Units'] = "TB"

for Server in UsageByServer:
    GroupName = "/Servers/{0:s}".format(Server)
    OutFile.create_group(GroupName)
    for User in UsageByServer[Server]:
        if (User == "SIZE"):
            OutFile[GroupName].attrs['Size'] = UsageByServer[Server][User]
        elif (User == "USED"):
            OutFile[GroupName].attrs['Used'] = UsageByServer[Server][User]
        elif (User == "FREE"):
            OutFile[GroupName].attrs['Free'] = UsageByServer[Server][User]
        elif (User == "PCOLOR"):
            OutFile[GroupName].attrs['Pcolor'] = UsageByServer[Server][User]
        else:
            DsetName = "{0:s}/{1:s}".format(GroupName, User)
            OutFile[DsetName] = UsageByServer[Server][User]

for User in UsageByUser:
    GroupName = "/Users/{0:s}".format(User)
    OutFile.create_group(GroupName)
    for Server in UsageByUser[User]:
        if (Server == "SIZE"):
            OutFile[GroupName].attrs['Size'] = UsageByUser[User][Server]
        elif (Server == "USED"):
            OutFile[GroupName].attrs['Used'] = UsageByUser[User][Server]
        elif (Server == "FREE"):
            OutFile[GroupName].attrs['Free'] = UsageByUser[User][Server]
        elif (Server == "PCOLOR"):
            OutFile[GroupName].attrs['Pcolor'] = UsageByUser[User][Server]
        else:
            DsetName = "{0:s}/{1:s}".format(GroupName, Server)
            OutFile[DsetName] = UsageByUser[User][Server]

OutFile.close()
