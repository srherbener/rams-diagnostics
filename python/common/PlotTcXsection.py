###################################################################
# Plot cross sections in a TC (azimuthally averaged quantities)

import matplotlib.pyplot as plt
import numpy as np
import h5py
import PlotUtils as plu
import ConfigTsd as conf

class StormXsection:
    '''Class for plotting storm cross sections, including factors'''

    def __init__(self, InFname, InVname, Ptitle, Tperiod, Cspecs, FactCspecs, OutFname):
        self.InFname = InFname
        self.InVname = InVname

        self.Title = Ptitle
        self.Tperiod = Tperiod

        self.Cmin = Cspecs[0]
        self.Cmax = Cspecs[1]
        self.Cnum = Cspecs[2]
        self.Cmap = 'nipy_spectral'

        self.FactCmin = FactCspecs[0]
        self.FactCmax = FactCspecs[1]
        self.FactCnum = FactCspecs[2]
        self.FactCmap = 'bwr'

        self.Cfilled = 1

        self.OutFname = OutFname


    def PlotXsection(self):

        # get labels for sims and factors
        LabelScheme = conf.SetLabelScheme()

        CaseList = [ 
            'TSD_NONSAL_NODUST',
            'TSD_SAL_NODUST',
            'TSD_NONSAL_DUST',
            'TSD_SAL_DUST'
            ]
        Ncases = len(CaseList)

        FacList = [
            'FAC_SAL',
            'FAC_DUST',
            'FAC_INT',
            ]
        Nfacs = len(FacList)

        # Collect input data in a 3D array: (x,z,case) 
        Titles = [ ]
        for i in range(0, Ncases):
            Case  = CaseList[i]
        
            InFname = self.InFname.replace("<CASE>", Case)
            InVname = self.InVname
            print("Reading {0:s} ({1:s})".format(InFname, InVname))
        
            InFile = h5py.File(InFname, mode='r')
        
            # grab the coordinate values, and initialize the wind array
            if (i == 0):
                X = InFile['/x_coords'][:] / 1000.0  # radius in km
                Z = InFile['/z_coords'][:] / 1000.0  # height in km
                Nx = len(X)
                Nz = len(Z)
                InVar = np.zeros( [ Nz, Nx, Ncases ], dtype=np.float_)
        
            InVar[ :, :, i ] = InFile[InVname][:,:]
            
            InFile.close()
        
            Titles.append(r'${0:s}\ ({1:s})$'.format(self.Title, LabelScheme[Case]))
        
        Pmarkers = [ 'a', 'c', 'e', 'g' ]
        
        # Build the factors from the input cross sections
        #   The third dimension in InVar is the case which
        #   is mapped as:
        #
        #      index     case
        #        0       NSND
        #        1        SND
        #        2        NSD
        #        3         SD
        #
        # Factors are in this order:
        #
        #      index     factor            formula
        #        0       FAC_SAL     SND - NSND
        #        1      FAC_DUST     NSD - NSND
        #        2       FAC_INT     SD - (SND + NSD) + NSND
        
        FactTitles = [ ]
        FactVar = np.zeros( [ Nz, Nx, Nfacs ], dtype=np.float_)
        for i in range(0, Nfacs):
            Factor = FacList[i]
        
            if (Factor == 'FAC_SAL'):
                # SND - NSND
                FactVar[:,:,i] = InVar[:,:,1] - InVar[:,:,0]
            elif (Factor == 'FAC_DUST'):
                # NSD - NSND
                FactVar[:,:,i] = InVar[:,:,2] - InVar[:,:,0]
            elif (Factor == 'FAC_INT'):
                # SD - (SND + NSD) + NSND
                FactVar[:,:,i] = InVar[:,:,3] - (InVar[:,:,1] + InVar[:,:,2]) + InVar[:,:,0]
            else:
                print("ERROR: PlotXsection: Attempted to calculate too many factors")
                print("ERROR: PlotXsection: Stopping!")
                sys.exit(1)
        
            FactTitles.append(r'${0:s}\ ({1:s})$'.format(self.Title, LabelScheme[Factor]))
        
        FacPmarkers = [ 'd', 'f', 'h' ]
        
        # 7 panel plot
        #   Sims on the left side, 4 panels
        #   Factors on right side, 3 panels
        FigPsap = plt.figure()
        
        # Place axes for all panels
        AxW = 0.40
        AxH = 0.16
        AxLeftLlx = 0.07
        AxRightLlx = 0.53
        AxLly = [ 0.10, 0.32, 0.54, 0.76 ]
        
        SimAx = [ ]
        SimAx.append(FigPsap.add_axes([ AxLeftLlx, AxLly[3], AxW, AxH ]))
        SimAx.append(FigPsap.add_axes([ AxLeftLlx, AxLly[2], AxW, AxH ]))
        SimAx.append(FigPsap.add_axes([ AxLeftLlx, AxLly[1], AxW, AxH ]))
        SimAx.append(FigPsap.add_axes([ AxLeftLlx, AxLly[0], AxW, AxH ]))
        
        FactAx = [ ]
        FactAx.append(FigPsap.add_axes([ AxRightLlx, AxLly[2], AxW, AxH ]))
        FactAx.append(FigPsap.add_axes([ AxRightLlx, AxLly[1], AxW, AxH ]))
        FactAx.append(FigPsap.add_axes([ AxRightLlx, AxLly[0], AxW, AxH ]))
        
        # for placing text in upper right (blank) space
        # create axes, turn off the axes, then place the text in that axes
        TxtAxPsap = FigPsap.add_axes([ AxRightLlx, AxLly[3], AxW, AxH ])
        TxtAxPsap.set_axis_off()
        
        Xaxis = plu.AxisConfig('x', [ 0, 500 ], "Radius (km)")
        Xaxis.fontsize = 12
        
        Yaxis = plu.AxisConfig('y', [ 0, 6 ], "Z (km)")
        Yaxis.ticks = [ 0, 2, 4, 6 ]
        Yaxis.fontsize = 12
        
        # Sim cross sections
        Clim = [ self.Cmin, self.Cmax ]
        Clevs = np.linspace(self.Cmin, self.Cmax, num=self.Cnum)
        Xshow = [ 0, 0, 0, 1 ]
        Yshow = [ 1, 1, 1, 1 ]
        for i in range(0, Ncases):
            Ptitle = plu.TitleConfig(Pmarkers[i], Titles[i])
            Ptitle.fontsize = 18
        
            Xaxis.show = Xshow[i]
            Yaxis.show = Yshow[i]
            plu.PlotContour(SimAx[i], X, Z, InVar[:,:,i], Ptitle, Xaxis, Yaxis, self.Cmap, Clim, Clevs, self.Cfilled)
        
        # factor cross sections
        FactClim = [ self.FactCmin, self.FactCmax ]
        FactClevs = np.linspace(self.FactCmin, self.FactCmax, num=self.FactCnum)
        Xshow = [ 0, 0, 1 ]
        Yshow = [ 0, 0, 0 ]
        for i in range(0, Nfacs):
            Ptitle = plu.TitleConfig(FacPmarkers[i], FactTitles[i])
            Ptitle.fontsize = 18
        
            Xaxis.show = Xshow[i]
            Yaxis.show = Yshow[i]
            plu.PlotContour(FactAx[i], X, Z, FactVar[:,:,i], Ptitle, Xaxis, Yaxis, self.FactCmap, FactClim, FactClevs, self.Cfilled)
        
        
        # transAxes says to use the Axes coordinates which always run from 0 to 1
        TxtAxPsap.text(1, 1, self.Tperiod, transform=TxtAxPsap.transAxes,
            fontsize=20, horizontalalignment='right', verticalalignment='bottom')
        
        print("Writing: {0:s}".format(self.OutFname))
        FigPsap.savefig(self.OutFname)
        
        print("")
