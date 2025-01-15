
# Copyright (c) 2020 Greg Pintilie - gregp@slac.stanford.edu

# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
# THE SOFTWARE.

import chimera
import os
import os.path
import Tkinter
import tkFont
from CGLtk import Hybrid
import VolumeData
import _multiscale
import MultiScale.surface
import _surface
import numpy
import _contour
import Matrix
import Surface
import VolumeViewer
import FitMap
from sys import stderr
from time import clock
import _contour
import chimera.match
import time


from axes import prAxes
import _multiscale
from CGLutil.AdaptiveTree import AdaptiveTree
import random
from VolumePath import Marker_Set, Marker, Link
from _contour import affine_transform_vertices as transform_vertices
from Matrix import xform_matrix, multiply_matrices, chimera_xform, identity_matrix, invert_matrix, shift_and_angle
import struct


from Rotamers import getRotamers
from chimera.resCode import protein1to3


try :
    from segment_dialog import current_segmentation, segmentation_map
    import molbuild
    reload (molbuild)
except :
    pass


#gSigma = 0.6
mapqVersion = "2.9.7"
#showDevTools = True

showDevTools = False

try :
    from Segger import showDevTools, timing, seggerVersion
    #from Segger import timing, seggerVersion
except :
    pass

import qscores
reload (qscores)

import mmcif
reload (mmcif)

OML = chimera.openModels.list

isModelZ = False

dlgName = "mapqdlg"
dlgTitle = "MapQ (v"+mapqVersion+")"
dlgHelp = 'https://github.com/gregdp/mapq'

if isModelZ :
    dlgName = "modelzdlg"
    dlgTitle = "ModelZ (v1.2)"
    dlgHelp = 'https://github.com/gregdp/modelz'


chargedIons = { "MG":2, "NA":1, "CL":-1, "CA":2, "ZN":2, "MN":2, "FE":3, "CO":2, "NI":2 }

atomColors = {'C' : chimera.MaterialColor (0.565,0.565,0.565),
            'Cbb' : chimera.MaterialColor (0.2,0.6,0.2),
            'S' : chimera.MaterialColor (1.000,1.000,0.188),
            'O' : chimera.MaterialColor (1.000,0.051,0.051),
            'N' : chimera.MaterialColor (0.188,0.314,0.973),
            'P' : chimera.MaterialColor (1.0, 0.502, 0.0),
            'H' : chimera.MaterialColor (0.9,.9,.9),
            ' ' : chimera.MaterialColor (0.2,1,.2),
            "MG" : chimera.MaterialColor (.4,.4,.6),
            "NA" : chimera.MaterialColor (.7,.4,.9),
            "CL" : chimera.MaterialColor (0,1,0),
            "CA" : chimera.MaterialColor (.4,.4,.6),
            "ZN" : chimera.MaterialColor (.4,.4,.6),
            "MN" : chimera.MaterialColor (.4,.4,.6),
            "FE" : chimera.MaterialColor (.4,.4,.6),
            "CO" : chimera.MaterialColor (.4,.4,.6),
            "NI" : chimera.MaterialColor (.4,.4,.6)
}


atomColors = {'C' : chimera.MaterialColor (0.565,0.565,0.565),
            'Cbb' : chimera.MaterialColor (0.2,0.6,0.2),
            'S' : chimera.MaterialColor (1.000,1.000,0.188),
            'O' : chimera.MaterialColor (1.000,0.051,0.051),
            'N' : chimera.MaterialColor (0.188,0.314,0.973),
            'P' : chimera.MaterialColor (1.0, 0.502, 0.0),
            'H' : chimera.MaterialColor (0.9,.9,.9),
            ' ' : chimera.MaterialColor (0.2,1,.2),
            "MG" : chimera.MaterialColor (0,1,0),
            "NA" : chimera.MaterialColor (.6,.3,.6),
            "CL" : chimera.MaterialColor (.2,.6,.2),
            "CA" : chimera.MaterialColor (.4,.4,.6),
            "ZN" : chimera.MaterialColor (.2,.8,.2),
            "MN" : chimera.MaterialColor (.4,.4,.6),
            "FE" : chimera.MaterialColor (.4,.4,.6),
            "CO" : chimera.MaterialColor (.4,.4,.6),
            "NI" : chimera.MaterialColor (.4,.4,.6)
}




def umsg ( txt ) :
    print txt
    status ( txt )

def status ( txt ) :
    txt = txt.rstrip('\n')
    msg.configure(text = txt)
    msg.update_idletasks()


class MapQ_Dialog ( chimera.baseDialog.ModelessDialog ) :

    name = dlgName
    if showDevTools :
        buttons = ( "Options", "Select", "Log", "Close" )
    else :
        buttons = ( "Options", "Log", "Close" )
    title = dlgTitle
    help = dlgHelp


    def fillInUI(self, parent):

        self.group_mouse_mode = None

        tw = parent.winfo_toplevel()
        self.toplevel_widget = tw
        tw.withdraw()

        self.parent = parent

        parent.columnconfigure(0, weight = 1)

        row = 0

        menubar = Tkinter.Menu(parent, type = 'menubar', tearoff = False)
        tw.config(menu = menubar)

        f = Tkinter.Frame(parent)
        f.grid(column=0, row=row, sticky='ew')

        #l = Tkinter.Label(f, text='  ')
        #l.grid(column=0, row=row, sticky='w')


        if 0 :
            menubar = Tkinter.Menu(parent, type = 'menubar', tearoff = False)
            tw.config(menu = menubar)

            file_menu_entries = (
                ('Open Model...', self.LoadModel),
                ('Save Model...', self.SaveModel)
                )
            fmenu = Hybrid.cascade_menu(menubar, 'File', file_menu_entries)

            from chimera.tkgui import aquaMenuBar
            aquaMenuBar(menubar, parent, row = 0, columnspan=3)



        # ---------------------------------------------------------------------------------

        self.InitVars()


        if 1 :
            #row += 1
            ff = Tkinter.Frame(f)
            ff.grid(column=0, row=row, sticky='nsew', pady=0, padx=0)

            Tkinter.Grid.columnconfigure(f, 0, weight=1)
            Tkinter.Grid.columnconfigure(ff, 0, weight=1)

            Tkinter.Grid.rowconfigure(f, row, weight=1)
            Tkinter.Grid.rowconfigure(ff, 0, weight=1)


            self.Canvas = Tkinter.Canvas(ff, height=80)
            self.Canvas.grid(column=0, row=0, sticky='nsew')

            self.modX = 10; self.modY = 10; self.modH = 30
            self.seqX = 10; self.seqY = 45; self.seqH = 30

            self.Canvas.bind("<ButtonPress-1>", lambda event : self.B1_Down ( event ) )
            self.Canvas.bind("<Control-ButtonPress-1>", lambda event : self.B1_Down_Ctrl ( event ) )
            self.Canvas.bind("<Shift-ButtonPress-1>", lambda event : self.B1_Down_Shift ( event ) )
            self.Canvas.bind("<Option-ButtonPress-1>", lambda event : self.B1_Down_Alt ( event ) )
            self.Canvas.bind("<Alt-ButtonPress-1>", lambda event : self.B1_Down_Alt ( event ) )
            self.Canvas.bind("<ButtonPress-2>", lambda event : self.B2_Down (event) )
            self.Canvas.bind("<ButtonPress-3>", lambda event : self.B3_Down (event) )
            self.Canvas.bind("<ButtonRelease-1>", lambda event : self.B1_Up ( event ) )
            self.Canvas.bind("<Control-ButtonRelease-1>", lambda event : self.B1_Up_Ctrl ( event ) )
            self.Canvas.bind("<Shift-ButtonRelease-1>", lambda event : self.B1_Up_Shift ( event ) )
            self.Canvas.bind("<Alt-ButtonRelease-1>", lambda event : self.B1_Up_Alt ( event ) )
            self.Canvas.bind("<Option-ButtonRelease-1>", lambda event : self.B1_Up_Alt ( event ) )

            self.Canvas.bind("<ButtonRelease-2>", lambda event : self.B2_Up (event) )
            self.Canvas.bind("<Option-ButtonRelease-2>", lambda event : self.B2_Up_Alt (event) )
            self.Canvas.bind("<Alt-ButtonRelease-2>", lambda event : self.B2_Up_Alt (event) )
            self.Canvas.bind("<Control-ButtonRelease-2>", lambda event : self.B2_Up_Ctrl (event) )
            self.Canvas.bind("<Command-ButtonRelease-2>", lambda event : self.B2_Up_Comm (event) )
            self.Canvas.bind("<Shift-ButtonRelease-2>", lambda event : self.B2_Up_Shift (event) )

            self.Canvas.bind("<ButtonRelease-3>", lambda event : self.B2_Up (event) )
            self.Canvas.bind("<Option-ButtonRelease-3>", lambda event : self.B2_Up_Alt (event) )
            self.Canvas.bind("<Alt-ButtonRelease-3>", lambda event : self.B2_Up_Alt (event) )
            self.Canvas.bind("<Control-ButtonRelease-3>", lambda event : self.B2_Up_Ctrl (event) )
            self.Canvas.bind("<Command-ButtonRelease-3>", lambda event : self.B2_Up_Comm (event) )
            self.Canvas.bind("<Shift-ButtonRelease-3>", lambda event : self.B2_Up_Shift (event) )

            self.Canvas.bind("<B1-Motion>", lambda event : self.B1_Drag ( event ) )
            self.Canvas.bind("<B2-Motion>", lambda event : self.B2_Drag ( event ) )
            self.Canvas.bind("<B3-Motion>", lambda event : self.B3_Drag ( event ) )
            self.Canvas.bind("<Motion>", lambda event : self.Mouse_Move ( event ) )
            self.Canvas.bind("<Configure>", lambda event : self.Canvas_Config (event) )
            self.Canvas.bind("<Leave>", lambda event : self.Canvas_Leave (event) )
            self.Canvas.bind("<MouseWheel>", lambda event : self.Canvas_Wheel (event) )


        row += 1
        ff = Tkinter.Frame(f)
        ff.grid(column=0, row=row, sticky='w', pady=0, padx=2)

        if 1 :
            ff = Tkinter.Frame(f)
            ff.grid(column=0, row=row, sticky='w', pady=1, padx=0)

            l = Tkinter.Label(ff, text=' Map:', anchor=Tkinter.W)
            l.grid(column=0, row=0, sticky='w')

            self.dmap = Tkinter.StringVar(parent)
            self.dmapMB  = Tkinter.Menubutton ( ff, textvariable=self.dmap, relief=Tkinter.RAISED, width=17 )
            self.dmapMB.grid (column=1, row=0, sticky='we', padx=1)
            self.dmapMB.menu  =  Tkinter.Menu ( self.dmapMB, tearoff=0, postcommand=self.MapMenu )
            self.dmapMB["menu"]  =  self.dmapMB.menu

            self.cur_dmap = None
            self.SetVisMap ()
            if self.cur_dmap == None :
                self.dmap.set ( "Open a map ..." )


            l = Tkinter.Label(ff, text='Model:', anchor=Tkinter.W)
            l.grid(column=2, row=0, sticky='w')

            self.struc = Tkinter.StringVar(parent)
            self.strucMB  = Tkinter.Menubutton ( ff, textvariable=self.struc, relief=Tkinter.RAISED, width=17 )
            self.strucMB.grid (column=3, row=0, sticky='we', padx=1)
            self.strucMB.menu  =  Tkinter.Menu ( self.strucMB, tearoff=0, postcommand=self.StrucMenu )
            self.strucMB["menu"]  =  self.strucMB.menu

            self.cur_mol = None
            self.cur_chains = []
            self.SetVisMol ()
            if self.cur_mol == None :
                self.struc.set ( "Open a model ..." )

            l = Tkinter.Label(ff, text=" Chain:" )
            l.grid(column=4, row=0, sticky='w')

            self.chain = Tkinter.StringVar(parent)
            self.chainMB  = Tkinter.Menubutton ( ff, textvariable=self.chain, relief=Tkinter.RAISED, width=5 )
            self.chainMB.grid (column=5, row=0, sticky='we', padx=1)
            self.chainMB.menu  =  Tkinter.Menu ( self.chainMB, tearoff=0, postcommand=self.ChainMenu )
            self.chainMB["menu"]  =  self.chainMB.menu

            b = Tkinter.Button(ff, text="...", command=self.LoadModel)
            b.grid (column=6, row=0, sticky='w', padx=1)

            if showDevTools :
                b = Tkinter.Button(ff, text="^", command=self.SaveModel)
                b.grid (column=7, row=0, sticky='w', padx=1)

            #b = Tkinter.Button(ff, text="S", command=self.SaveModel)
            #b.grid (column=7, row=0, sticky='w', padx=1)


            l = Tkinter.Label(ff, text=" Show:" )
            l.grid(column=12, row=0, sticky='w')


            b = Tkinter.Button(ff, text="Chain", command=self.AllChain)
            b.grid (column=13, row=0, sticky='w', padx=1)

            b = Tkinter.Button(ff, text="All", command=self.AllChains)
            b.grid (column=14, row=0, sticky='w', padx=1)

            b = Tkinter.Button(ff, text="Sel.", command=self.ShowOnlySel)
            b.grid (column=15, row=0, sticky='w', padx=1)

            b = Tkinter.Button(ff, text="At.", command=self.SetSelAtoms)
            b.grid (column=16, row=0, sticky='w', padx=1)

            b = Tkinter.Button(ff, text="Rib.", command=self.SetSelRibbon)
            b.grid (column=17, row=0, sticky='w', padx=1)

            b = Tkinter.Button(ff, text="SCs", command=self.ShowSCs)
            b.grid (column=18, row=0, sticky='w', padx=1)

            b = Tkinter.Button(ff, text="~SCs", command=self.HideSCs)
            b.grid (column=19, row=0, sticky='w', padx=1)

            b = Tkinter.Button(ff, text="W", command=self.Wire)
            b.grid (column=20, row=0, sticky='w', padx=1)

            b = Tkinter.Button(ff, text=" ", command=self.HideSel)
            b.grid (column=21, row=0, sticky='w', padx=0)



        if 1 :

            #l = Tkinter.Label(ff, text=' Go:', fg="#777")
            #l.grid(column=35, row=0, sticky='e')

            b = Tkinter.Button(ff, text="<", command=self.ZoomBegin)
            b.grid (column=38, row=0, sticky='w', padx=0)

            b = Tkinter.Button(ff, text=">", command=self.ZoomEnd)
            b.grid (column=39, row=0, sticky='w', padx=0)




        if 1 :
            row += 1
            op = Hybrid.Popup_Panel(f)
            ff = op.frame
            ff.grid(row = row, column = 0, sticky = 'news')
            ff.grid_remove()
            #ff.columnconfigure(0, weight=1)
            self.optionsPanel = op.panel_shown_variable
            self.optionsPanel.set ( showDevTools )

            b = Tkinter.Label(ff, text=" Resolution:")
            b.grid (column=0, row=0, sticky='w', padx=0, pady=1)

            self.mapRes = Tkinter.StringVar(f)
            self.mapRes.set ( "3" )
            e = Tkinter.Entry(ff, width=3, textvariable=self.mapRes)
            e.grid(column=1, row=0, sticky='w', padx=0, pady=1)

            b = Tkinter.Label(ff, text="A ")
            b.grid (column=2, row=0, sticky='w', padx=0, pady=1)

            b = Tkinter.Label(ff, text=" Sigma:")
            b.grid (column=3, row=0, sticky='w', padx=0, pady=1)

            self.sigma = Tkinter.StringVar(f)
            self.sigma.set ( "0.4" )
            e = Tkinter.Entry(ff, width=3, textvariable=self.sigma)
            e.grid(column=4, row=0, sticky='w', padx=0, pady=1)

            self.qmenu = Tkinter.StringVar(parent)
            self.qmenu.set ( "Q-scores ..." )
            self.qmenuMB  = Tkinter.Menubutton ( ff, textvariable=self.qmenu, relief=Tkinter.RAISED, width=11 )
            self.qmenuMB.grid (column=5, row=0, sticky='we', padx=1)
            self.qmenuMB.menu  =  Tkinter.Menu ( self.qmenuMB, tearoff=0, postcommand=self.QMenu )
            self.qmenuMB["menu"]  =  self.qmenuMB.menu

            self.cmenu = Tkinter.StringVar(parent)
            self.cmenu.set ( "Visualize ..." )
            self.cmenuMB  = Tkinter.Menubutton ( ff, textvariable=self.cmenu, relief=Tkinter.RAISED, width=11 )
            self.cmenuMB.grid (column=6, row=0, sticky='we', padx=1)
            self.cmenuMB.menu  =  Tkinter.Menu ( self.cmenuMB, tearoff=0, postcommand=self.CMenu )
            self.cmenuMB["menu"]  =  self.cmenuMB.menu

            #l = Tkinter.Label(ff, text='  |  Seq: ', fg="#000")
            #l.grid(column=21, row=0, sticky='ens')

            oft = Hybrid.Checkbutton(ff, 'Gaps', True)
            oft.button.grid(column = 22, row = 0, sticky = 'w')
            self.showGaps = oft.variable
            #self.showRibbon.set ( 1 )


            #l = Tkinter.Label(ff, text=' Select: ', fg="#000", font = 'TkCaptionFont')
            #l = Tkinter.Label(ff, text='  |  On Select: ', fg="#000")
            #l.grid(column=35, row=0, sticky='ens')

            oft = Hybrid.Checkbutton(ff, 'Extract', False)
            oft.button.grid(column = 37, row = 0, sticky = 'w')
            self.selExtract = oft.variable
            self.selExtract.set ( 1 )

            oft = Hybrid.Checkbutton(ff, 'Mesh', False)
            oft.button.grid(column = 38, row = 0, sticky = 'w')
            self.showMesh = oft.variable
            #self.showRibbon.set ( 1 )

            self.showLigands = Tkinter.IntVar()
            self.showLigands.set(True)
            oft = Tkinter.Checkbutton( ff, text="Ligands", variable=self.showLigands )
            oft.grid(column = 39, row = 0, sticky = 'w')

            #oft = Hybrid.Checkbutton(ff, 'Preserve', False, command=self.cb)
            #oft.button.grid(column = 39, row = 0, sticky = 'w')
            #self.preserveSel = oft.variable
            self.preserveSel = Tkinter.IntVar()
            oft = Tkinter.Checkbutton( ff, text="Keep", variable=self.preserveSel, command=self.preserveSelCb)
            oft.grid(column = 40, row = 0, sticky = 'w')
            #self.showRibbon.set ( 1 )

            self.preserveVol = Tkinter.IntVar()
            if 0 :
                oft = Tkinter.Checkbutton( ff, text="+Map", variable=self.preserveVol, command=self.preserveVolCb)
                oft.grid(column = 41, row = 0, sticky = 'w')
                #self.showRibbon.set ( 1 )

            self.showH = Tkinter.IntVar()
            oft = Tkinter.Checkbutton( ff, text="H", variable=self.showH)
            oft.grid(column = 42, row = 0, sticky = 'w')
            #self.showRibbon.set ( 1 )

            self.showW = Tkinter.IntVar()
            if 0 :
                oft = Tkinter.Checkbutton( ff, text="W", variable=self.showW)
                oft.grid(column = 43, row = 0, sticky = 'w')


            b = Tkinter.Button(ff, text="<", command=self.KeepBack)
            b.grid (column=45, row=0, sticky='w', padx=0)

            b = Tkinter.Button(ff, text="Zone", command=self.SelReLoad)
            b.grid (column=46, row=0, sticky='w', padx=0)

            if 0 and showDevTools :

                b = Tkinter.Button(ff, text="L", command=self.SelLoad)
                b.grid (column=47, row=0, sticky='w', padx=2)

            b = Tkinter.Button(ff, text="Near", command=self.ShowNear)
            b.grid (column=47, row=0, sticky='w', padx=2)

            b = Tkinter.Button(ff, text="W&I", command=self.ShowNearWI)
            b.grid (column=48, row=0, sticky='w', padx=2)

            if showDevTools :
                b = Tkinter.Button(ff, text="x", command=self.FindClashes)
                b.grid (column=49, row=0, sticky='w', padx=2)

            if 0 :
                b = Tkinter.Button(ff, text="Zone", command=self.Zone)
                b.grid (column=48, row=0, sticky='w', padx=1, pady=1)

                self.zoneRad = Tkinter.StringVar(ff)
                self.zoneRad.set ( "2" )
                e = Tkinter.Entry(ff, width=2, textvariable=self.zoneRad)
                e.grid(column=49, row=0, sticky='w', padx=1, pady=1)



            #l = Tkinter.Label(ff, text='   Zoom:', fg="#777")
            l = Tkinter.Label(ff, text=' Zoom:' )
            l.grid(column=60, row=0, sticky='e')

            b = Tkinter.Button(ff, text="-", command=self.ZoomMinus)
            b.grid (column=61, row=0, sticky='w', padx=0)

            b = Tkinter.Button(ff, text="+", command=self.ZoomPlus)
            b.grid (column=62, row=0, sticky='w', padx=0)


        #b = Tkinter.Button(ff, text="D", command=self.Domains)
        #b.grid (column=18, row=0, sticky='w', padx=2)

        #b = Tkinter.Button(ff, text="S", command=self.SS)
        #b.grid (column=19, row=0, sticky='w', padx=2)



        # ----------- select panel ----------------------------------

        if 1 and showDevTools :

            row += 1
            op = Hybrid.Popup_Panel(f)
            ff = op.frame
            ff.grid(row = row, column = 0, sticky = 'news')
            ff.grid_remove()
            #ff.columnconfigure(0, weight=1)
            self.selPanel = op.panel_shown_variable
            self.optionsPanel.set ( showDevTools )


            #ff = Tkinter.Frame(f)
            #ff.grid(column=0, row=row, sticky='w', pady=0, padx=2)

            #l = Tkinter.Label(ff, text=' Sel:', font = 'TkCaptionFont')
            #l.grid(column=1, row=0, sticky='w', pady=1)


            if 0 :
                #b = Tkinter.Button(ff, text="Asp", command=self.asp )
                #b.grid (column=1, row=0, sticky='w', padx=2)

                b = Tkinter.Button(ff, text="Extr", command=self.Extract )
                b.grid (column=2, row=0, sticky='w', padx=2)

                b = Tkinter.Button(ff, text="Al 1", command=self.AlignRes1 )
                b.grid (column=3, row=0, sticky='w', padx=2)

                b = Tkinter.Button(ff, text="Al 2", command=self.AlignRes2 )
                b.grid (column=4, row=0, sticky='w', padx=2)

                b = Tkinter.Button(ff, text="Avg", command=self.Avg )
                b.grid (column=5, row=0, sticky='w', padx=2)

                b = Tkinter.Button(ff, text="~Extr", command=self.CloseExtracted )
                b.grid (column=6, row=0, sticky='w', padx=2)


                #b = Tkinter.Button(ff, text="Sbb", command=self.BB_Sigma )
                #b.grid (column=8, row=0, sticky='w', padx=2)

                #b = Tkinter.Button(ff, text="Z", command=self.ZScoreSel )
                #b.grid (column=9, row=0, sticky='w', padx=2)

                #b = Tkinter.Button(ff, text="Zr", command=self.RotaZ1 )
                #b.grid (column=10, row=0, sticky='w', padx=2)

                #b = Tkinter.Button(ff, text="R1", command=self.R1 )
                #b.grid (column=11, row=0, sticky='w', padx=2)

                #b = Tkinter.Button(ff, text="ExA", command=self.ExCustA )
                #b.grid (column=12, row=0, sticky='w', padx=2)

                #b = Tkinter.Button(ff, text="ExB", command=self.ExCustB )
                #b.grid (column=13, row=0, sticky='w', padx=2)

                #b = Tkinter.Button(ff, text="ExC", command=self.ExCustC )
                #b.grid (column=14, row=0, sticky='w', padx=2)


            #b = Tkinter.Button(ff, text="S-sel", command=self.S_sel )
            #b.grid (column=20, row=0, sticky='w', padx=2)

            b = Tkinter.Button(ff, text="Q-sel", command=lambda:self.Q_sel(show=False, fitG=False) ) # lambda event : self.B1_Drag ( event )
            b.grid (column=21, row=0, sticky='w', padx=2)

            b = Tkinter.Button(ff, text="Q-sel-show", command=lambda:self.Q_sel(show=True, fitG=False) ) # lambda event : self.B1_Drag ( event )
            b.grid (column=22, row=0, sticky='w', padx=2)

            b = Tkinter.Button(ff, text="Q-sel-g", command=lambda:self.Q_sel(show=True, fitG=True) ) # lambda event : self.B1_Drag ( event )
            b.grid (column=23, row=0, sticky='w', padx=2)

            if 0 :
                b = Tkinter.Button(ff, text="Q-show", command=self.Q_show )
                b.grid (column=22, row=0, sticky='w', padx=2)

                b = Tkinter.Button(ff, text="SA-Q", command=self.SA_Q )
                b.grid (column=23, row=0, sticky='w', padx=2)


            #b = Tkinter.Button(ff, text="Ats", command=self.ShowAts)
            #b.grid (column=25, row=0, sticky='w', padx=10)

            if 1 :
                b = Tkinter.Button(ff, text="Alts", command=self.FindAlts)
                b.grid (column=28, row=0, sticky='w', padx=2)

                b = Tkinter.Button(ff, text="X-Alts", command=self.DelAlts)
                b.grid (column=29, row=0, sticky='w', padx=2)

            if 0 :
                b = Tkinter.Button(ff, text="APro", command=self.AProfs)
                b.grid (column=28, row=0, sticky='w', padx=2)

            #b = Tkinter.Button(ff, text="Ligs", command=self.Ligs)
            #b.grid (column=43, row=0, sticky='w', padx=2)

            #b = Tkinter.Button(ff, text="Scale", command=self.Scale)
            #b.grid (column=44, row=0, sticky='w', padx=2)



            b = Tkinter.Label(ff, text="   Str:")
            b.grid (column=30, row=0, sticky='w', padx=0, pady=1)

            self.selText = Tkinter.StringVar(f)
            self.selText.set ( "" )
            e = Tkinter.Entry(ff, width=20, textvariable=self.selText)
            e.grid(column=31, row=0, sticky='w', padx=2, pady=1)


            b = Tkinter.Button(ff, text="Sel", command=self.SelText)
            b.grid (column=32, row=0, sticky='w', padx=2)



            b = Tkinter.Label(ff, text="Rad:")
            b.grid (column=33, row=0, sticky='w', padx=0, pady=1)

            self.maskRad = Tkinter.StringVar(f)
            self.maskRad.set ( "2.5" )
            e = Tkinter.Entry(ff, width=3, textvariable=self.maskRad)
            e.grid(column=34, row=0, sticky='w', padx=2, pady=1)


            b = Tkinter.Button(ff, text="AddSel", command=self.AdSel)
            b.grid (column=35, row=0, sticky='w', padx=2)


            b = Tkinter.Button(ff, text="Ds", command=self.ShowDists)
            b.grid (column=41, row=0, sticky='w', padx=2)

            b = Tkinter.Button(ff, text="Inter", command=self.Inter)
            b.grid (column=42, row=0, sticky='w', padx=2)


            b = Tkinter.Button(ff, text="Occ", command=self.Occ)
            b.grid (column=43, row=0, sticky='w', padx=2)


            b = Tkinter.Button(ff, text="Rmsd", command=self.RMSD)
            b.grid (column=44, row=0, sticky='w', padx=2)

            b = Tkinter.Button(ff, text="RibD", command=self.RibD)
            b.grid (column=45, row=0, sticky='w', padx=2)

            b = Tkinter.Button(ff, text="Af", command=self.AfColor)
            b.grid (column=46, row=0, sticky='w', padx=2)

            b = Tkinter.Button(ff, text="T", command=self.TunnelVis)
            b.grid (column=47, row=0, sticky='w', padx=2)

            b = Tkinter.Button(ff, text="C", command=self.CompareChains)
            b.grid (column=48, row=0, sticky='w', padx=2)

            self.selPanel.set(True)





        dummyFrame = Tkinter.Frame(parent, relief='groove', borderwidth=1)
        Tkinter.Frame(dummyFrame).pack()
        dummyFrame.grid(row=row,column=0,columnspan=7, padx=1, sticky='we')
        row += 1


        global msg
        msg = Tkinter.Label(parent, width = 60, anchor = 'w', justify = 'left', fg="red", pady=1, padx=10)
        msg.grid(column=0, row=row, sticky='ew')
        self.msg = msg

        self.showingAtoms = False

        if len ( self.cur_chains ) > 0 :
            self.chain.set ( self.cur_chains[0] )
            #self.ShowCh ( self.cur_chains[0] )
            self.GetSeq ()
        else :
            self.chain.set ( '-' )

        #umsg ( 'Select one or more segmented regions then press "Place Points" to start' )

        callbacks = (self.mouse_down_cb, self.mouse_drag_cb, self.mouse_up_cb)
        #callbacks = (self.mouse_down_cb)
        from chimera import mousemodes
        mousemodes.addFunction('mark mapq', callbacks, self.mouse_mode_icon())

        if 0 :
            # bind, unbind in case it was left bound before...
            from chimera import mousemodes
            print " - unbinding mouse..."
            button, modifiers = ('3', ['Ctrl'])
            def_mode = mousemodes.getDefault(button, modifiers)
            mousemodes.setButtonFunction(button, modifiers, def_mode)
            self.bound_button = None


        chimera.openModels.addRemoveHandler(self.ModelClosed, None)


    def ModelClosed(self, trigger, n, mlist):

        # Clear menus that are showing closed models.
        if self.cur_dmap in mlist:
            self.cur_dmap = None
            self.dmap.set ( "Select a map ..." )

        if self.cur_mol in mlist:
            self.struc.set ( "Select a model ..." )
            self.cur_mol = None
            self.cur_chains = []
            self.chain.set ( "-" )
            self.RemoveSeq ()
            #self.UpdateSeqFont ()
            #self.UpdateSeq ()


    def bind_placement_button_cb(self) :

        if self.use_mouse.get() :
            print " - binding mouse..."
            button, modifiers = ('3', ['Ctrl'])
            from chimera import mousemodes
            mousemodes.setButtonFunction(button, modifiers, 'mark mapq')
            self.bound_button = (button, modifiers)
        elif self.bound_button:
            print " - unbinding mouse..."
            button, modifiers = self.bound_button
            from chimera import mousemodes
            def_mode = mousemodes.getDefault(button, modifiers)
            mousemodes.setButtonFunction(button, modifiers, def_mode)
            self.bound_button = None


    def mouse_mode_icon(self) :

        import os.path
        icon_path = os.path.join(os.path.dirname(__file__), 'marker.gif')
        from PIL import Image
        image = Image.open(icon_path)
        from chimera import chimage
        from chimera import tkgui
        icon = chimage.get(image, tkgui.app)
        return icon

    def mouse_down_cb(self, viewer, event) :

        print " mouse - "

        #print event.x, event.y
        if 0 :
            print dir(event)
            print event.char
            print event.keycode
            print event.keysym
            print event.keysym_num
            print event.num
            print event.state

        hits = []
        import VolumePath.tracer as tracer

        if 1 :
            from VolumeViewer import volume_list
            hits.extend(tracer.volume_maxima(event.x, event.y, volume_list()))
            print "vol"

        if 0 :
            from VolumeViewer import volume_list
            hits.extend(VolumePath.tracer.volume_plane_intercepts(event.x, event.y, volume_list()))

        if 0 :
            from Surface import surface_models
            hits.extend(tracer.surface_intercepts(event.x, event.y, surface_models()))
            print "surf"

        for C, vol in hits :
            print " --> ", vol.name, " --> %.1f, %.1f, %.1f" % (C[0], C[1], C[2])
            self.PlaceAt ( C, vol )





        #grabbed = (self.move_markers.get() and self.grab_marker(event.x, event.y))
        #if not grabbed:
        #    self.add_marker_at_screen_xy(event.x, event.y)



    def mouse_drag_cb(self, viewer, event):
        shift_mask = 1
        shift = (event.state & shift_mask)
        capslock_mask = 2
        capslock = (event.state & capslock_mask)
        #self.move_or_resize_marker(event.x, event.y, shift, capslock):


    def mouse_up_cb(self, viewer, event):
        #self.ungrab_marker()
        #self.pause_marker_placement = False
        #print "mouse up"
        pass



    def Options ( self ) :
        self.optionsPanel.set (not self.optionsPanel.get())


    def Select ( self ) :
        self.selPanel.set (not self.selPanel.get())


    def Log ( self ) :
        import Idle
        Idle.start_shell()


    def InitVars ( self ) :

        self.mag = 13
        self.seqt = []
        self.boldSeqT = None
        self.drag = ''

        #self.sheetBaseClr = numpy.array ( [50.0,205.0,50.0] )
        #self.sheetClr = numpy.array ( [204.0,255.0,204.0] )
        self.sheetBaseClr = numpy.array ( [55.0,55.0,150.0] )
        self.sheetClr = numpy.array ( [150.0,150.0,250.0] )
        self.sheetClrD = self.sheetClr - self.sheetBaseClr

        self.helixBaseClr = numpy.array ( [150.0,50.0,50.0] )
        self.helixClr = numpy.array ( [255.0,150.0,150.0] )
        self.helixClrD = self.helixClr - self.helixBaseClr

        c = self.helixBaseClr; self.helix1 = "#" + struct.pack('BBB',c[0],c[1],c[2]).encode('hex')
        c = self.helixClr;     self.helix2 = "#" + struct.pack('BBB',c[0],c[1],c[2]).encode('hex')

        self.switch = "#522"

        c = self.sheetBaseClr; self.strand1 = "#77F"
        c = self.sheetClr;     self.strand2 = "#77F"

        c = self.sheetBaseClr; self.sheet1 = "#" + struct.pack('BBB',c[0],c[1],c[2]).encode('hex')
        c = self.sheetClr;     self.sheet2 = "#" + struct.pack('BBB',c[0],c[1],c[2]).encode('hex')

        self.loop1 = "#999"

        self.selColor = "#7e7"


        self.font = tkFont.Font(family='Courier', size=(self.mag), weight='normal')
        #self.boldFont = tkFont.Font(family='Courier', size=(self.mag+4), weight='bold')
        self.tw = self.font.measure ( "a" )

        self.seq = ""

        #self.OrderMods ()


    def SetVisMap ( self ) :
        dmap = None
        mlist = OML(modelTypes = [VolumeViewer.volume.Volume])
        for m in mlist :
            if m.display and not "sel_masked" in m.name :
                dmap = m
                break

        if dmap == None :
            if len(mlist) > 0 :
                dmap = mlist[0]

        if dmap != None :
            self.dmap.set ( "[%d] %s" % (dmap.id,dmap.name) )
            self.cur_dmap = dmap


    def QMenu ( self ) :
        self.qmenuMB.menu.delete ( 0, 'end' )   # Clear menu
        options = []
        options.append ( "Calc. for selected atoms" )
        options.append ( "Calc. for all atoms in chain (single process)" )
        options.append ( "Calc. for all atoms in chain (2 processes)")
        options.append ( "Calc. for all atoms in chain (4 processes)")
        options.append ( "Calc. for all atoms in chain (8 processes)")
        options.append ( "Calc. for all atoms in chain (auto # processses)" )
        options.append ( "Load (and calculate stats for selected chain)" )

        if showDevTools :
            options.append ( "---" )
            options.append ( "Per-residue CC" )
            options.append ( "B-factors from Q-scores" )
            #options.append ( "Map to phenix-model-map CC" )
            #options.append ( "R-score for [selected atoms]" )
            options.append ( "sseQ - Show" )

        for op in options :
            self.qmenuMB.menu.add_radiobutton ( label=op, variable=self.qmenu,
                                command=lambda o=op: self.QMenuSelected(o) )

    def QMenuSelected ( self, op ) :

        self.qmenu.set ( "Q-scores ..." )
        print op
        if "R-score" in op :
            self.CalcSelR ()
        elif "B-factors from Q-scores" in op :
            self.CalcBFactors ()
        elif "Map to phenix-model-map CC" in op :
            self.PhenixCC ()
        elif "Calc." in op :
            if "selected atoms" in op.lower() :
                self.CalcSelQ ()
            elif "single" in op :
                self.CalcAllQ()
            elif "2" in op :
                self.CalcAllQp(2)
            elif "4" in op :
                self.CalcAllQp(4)
            elif "8" in op :
                self.CalcAllQp(8)
            elif "auto" in op :
                self.CalcAllQp()
        elif "Load" in op :
            self.GetQsFromFile ()
        elif "sseQ - Show" in op :
            self.CalcSseQSel ()
        elif op == "Fit Scores" :
            import fit
            reload ( fit )
            fit.fitScores ( float(self.mapRes.get()) )
        elif op == "Per-residue CC" :
            self.ResCC ()



    def CMenu ( self ) :
        self.cmenuMB.menu.delete ( 0, 'end' )   # Clear menu
        options = []
        options.append ( "Residue Q-scores" )
        options.append ( "Backbone Q-scores" )
        options.append ( "Sidechain Q-scores" )
        options.append ( "Atom Q-scores" )
        options.append ( "Random color all chains" )

        if showDevTools :
            options.append ( "------------" )
            options.append ( "Expand selection to all bonded atoms" )
            options.append ( "Plot" )
            options.append ( "Domains" )

        for op in options :
            self.cmenuMB.menu.add_radiobutton ( label=op, variable=self.cmenu,
                                command=lambda o=op: self.CMenuSelected(o) )


    def CMenuSelected ( self, op ) :

        self.cmenu.set ( "Visualize ..." )
        print op
        if "Residue" in op :
            self.DoColorRes()
        elif "Backbone" in op :
            self.DoColorBB()
        elif "Sidechain" in op :
            self.DoColorSC()
        elif "Atom" in op :
            self.DoColorAtoms ()
        elif "Random" in op :
            self.DoColorRandom ()

        elif "bonded" in op :
            self.ExpandSelBonded ()
        elif "Plot" in op :
            self.PlotQ ()
        elif op == "Domains" :
            self.ColorDomains ()


    def MapMenu ( self ) :
        #print "Map menu..."
        self.dmapMB.menu.delete ( 0, 'end' )   # Clear menu
        #self.cur_dmap = None
        #self.dmap.set("")
        mlist = OML(modelTypes = [VolumeViewer.volume.Volume])
        if len(mlist) == 0 :
            self.LoadModel()
            return
        for m in mlist :
            self.dmapMB.menu.add_radiobutton ( label="[%d] %s"%(m.id,m.name), variable=self.dmap,
                                command=lambda m=m: self.MapSelected(m) )


    def MapSelected ( self, dmap ) :

        self.cur_dmap = dmap
        print "Selected " + dmap.name
        self.dmap.set( "[%d] %s" % (dmap.id, dmap.name) )

        #self.GetSeq ()
        #self.ZoomBegin ()


    def GetChains ( self, mol ) :
        ct = {}
        for r in mol.residues:
            ct[r.id.chainId] = 1
        clist = ct.keys()
        clist.sort()
        return clist


    def SetVisMol ( self ) :
        mol = None
        mlist = OML(modelTypes = [chimera.Molecule])
        for m in mlist :
            if m.display :
                mol = m
                break

        if mol == None :
            if len(mlist) > 0 :
                mol = mlist[0]

        if mol != None :
            self.struc.set ( "[%d] %s" % (mol.id, mol.name) )
            self.cur_mol = mol
            self.cur_chains = self.GetChains ( mol )
            qscores.SetBBAts ( mol )


    def StrucSelected ( self, mol ) :

        self.cur_mol = mol
        print "Selected ", mol.name, " - ", mol.id
        if mol :

            #mlist = OML(modelTypes = [chimera.Molecule])
            #for m in mlist :
            #    m.display = False

            mol.display = True

            self.cur_chains = self.GetChains ( mol )

            if len(self.cur_chains) == 0 :
                self.chain.set ( "-" )
            elif self.chain.get() in self.cur_chains :
                print " - ch " + self.chain.get() + " already sel"
                #self.ShowCh ( self.chain.get() )
            else :
                self.chain.set ( self.cur_chains[0] )
                #self.ShowCh ( self.chain.get() )


        qscores.SetBBAts ( mol )
        self.parent.after(100, self.DoSeq)


    def DoSeq ( self ) :
        #print "after 100"

        self.GetSeq ()
        #self.ZoomBegin ()

        if self.cur_mol != None :
            self.ShowQScores ()


    def ChainSelected ( self, ch ) :
        print " - sel chain: ", ch, self.chain.get()
        #self.ShowCh ( ch )
        self.parent.after(100, self.DoSeq)



    def StrucMenu ( self ) :
        self.strucMB.menu.delete ( 0, 'end' )   # Clear menu
        mlist = OML(modelTypes = [chimera.Molecule])
        if len(mlist) == 0 :
            self.LoadModel()
            return
        for m in mlist :
            self.strucMB.menu.add_radiobutton ( label="[%d] %s"%(m.id,m.name), variable=self.struc,
                                           command=lambda m=m: self.StrucSelected(m) )

    def ChainMenu ( self ) :
        self.chainMB.menu.delete ( 0, 'end' )   # Clear menu
        #print " - chain menu"
        #print self.cur_chains
        for ch in self.cur_chains :
            self.chainMB.menu.add_radiobutton ( label=ch, variable=self.chain, command=lambda ch=ch: self.ChainSelected(ch) )

        self.chainMB.menu.add_radiobutton ( label="All", variable=self.chain, command=lambda ch="All": self.ChainSelected("All") )



    def loadFile ( self, path ) :

        print "Loading - %s" % path

        ext = os.path.splitext ( path )[1]
        mols = []

        if ext == ".cif" :
            start = time.time()
            mols = mmcif.LoadMol2 ( path, log=False )
            print "Loaded %d mols in %.1fs" % ( len(mols), time.time()-start)


        elif ext == ".mrc" or ext == ".map" or ext == ".ccp4" or ext == ".hdf" :
            om = chimera.openModels.open ( path )[0]
            chimera.runCommand ( "vol #%d style surface region all step 1" % om.id )
            for sp in om.surfacePieces :
                v, t = sp.geometry
                if len(v) == 8 and len(t) == 12 :
                    sp.display = False
                else :
                    #sp.displayStyle = sp.Mesh
                    sp.color = sp.color = (.7, .7, .7, .7)

            self.MapSelected (om)

        elif ext == ".pdb" or ext == ".ent" :
            mols = chimera.openModels.open ( path )

            for mol in mols :
                mmcif.ColorMol ( mol )

                #def nucleicOff():
                from NucleicAcids.cmd import sidechain
                sidechain("atoms", sel="#%d" % mol.id)
                from chimera.resCode import nucleic3to1
                for r in mol.residues :
                    if r.type in nucleic3to1 :
                        r.fillDisplay = False

            mol = mols[0]
            self.struc.set ( "[%d] %s" % (mol.id, mol.name) )
            self.cur_mol = mol
            self.cur_chains = self.GetChains ( mol )
            if len ( self.cur_chains ) > 0 :
                self.chain.set ( self.cur_chains[0] )
                self.GetSeq ()
                self.ShowQScores ()



    def load ( self, okay, dialog ):
        if okay:
            paths = dialog.getPaths ( )
            print "%d files" % len(paths)

            # load maps first...
            for path in paths :
                ext = os.path.splitext ( path )[1]
                if ext == ".mrc" or ext == ".map" or ext == ".ccp4" or ext == ".hdf":
                    self.loadFile ( path )

            # then models...
            for path in paths :
                ext = os.path.splitext ( path )[1]
                if ext == ".pdb" or ext == ".cif" or ext == ".ent" :
                    self.loadFile ( path )

    def LoadModel ( self ) :
        init = None
        mol = None
        for m in chimera.openModels.list() :
            if type(m) == chimera.Molecule and m.display == True and hasattr ( m, 'openedAs' ) :
                init = os.path.split ( m.openedAs[0] ) [0]
                break
            if type(m) == VolumeViewer.volume.Volume :
                init = os.path.split( m.data.path ) [0]

        if init == None :
            init = "/Users/greg/Box Sync/_data"

        print "init: %s" % init

        if 1 :
            from OpenSave import OpenModeless
            OpenModeless ( title = 'Open Model',
                           #filters = [('TXT', '*.txt', '.txt')],
                           filters = [],
                           initialfile = init, command = self.load )

        else :
            fpath = "/Users/greg/Box Sync/_data/problems/emd_30342/7cec.cif"



    def save ( self, okay, dialog ):
        print "save..."
        if okay:
            paths = dialog.getPaths ( )
            print "%d files" % len(paths)
            for path in paths :
                print path
                fname = os.path.split ( path )[1]
                fext = os.path.splitext ( fname )[1]
                print " ->", path

                if fext == ".cif" :
                    mmcif.WriteMol ( self.cur_mol, path, dmap = self.cur_dmap )
                    self.cur_mol.name = fname
                    self.cur_mol.openedAs = [ path, [] ]
                    self.struc.set ( "[%d] %s" % (self.cur_mol.id, fname) )

                elif fext == ".pdb" or fext == ".ent" :
                    if self.cur_dmap != None :
                        for at in self.cur_mol.atoms :
                            at.c0 = at.coord()
                            at.setCoord ( self.cur_dmap.openState.xform.inverse().apply ( at.xformCoord() ) )

                    print "."
                    chimera.PDBio().writePDBfile ( [self.cur_mol], path )

                    if self.cur_dmap != None :
                        for at in self.cur_mol.atoms :
                            at.setCoord ( at.c0 )
                            del at.c0

                    print "."
                    self.cur_mol.name = fname
                    self.cur_mol.openedAs = [ path, [] ]
                    self.struc.set ( "[%d] %s" % (self.cur_mol.id, fname) )

                    if 1 :
                        print " -< %s " % path
                        fin = open ( path, "r" )
                        fout = open ( os.path.splitext(path)[0] + "_AT0.pdb", "w" )
                        for line in fin :
                            if line[0:4] == "ATOM" or line[0:6] == "HETATM" :
                                fout.write ( line )
                        fin.close()
                        fout.close()



    def SaveModel ( self ) :

        if self.cur_mol == None :
            umsg ( "Select a model to save..." )
            return

        #initFile = os.path.splitext( self.cur_mol.name )[0] + ".cif"
        initFile = self.cur_mol.name

        from OpenSave import SaveModeless
        SaveModeless ( title = 'Save Model',
                       #filters = [('TXT', '*.txt', '.txt')],
                       filters = [],
                       initialfile = initFile, command = self.save )



    def DoColorBB ( self ) :

        self.UpdateModColor ( "bb" )
        # self.RandColorChains()

    def DoColorSC ( self ) :

        self.UpdateModColor ( "sc" )
        # self.RandColorChains()


    def DoColorRes ( self ) :

        self.UpdateModColor ( "res" )
        # self.RandColorChains()


    def DoColorAtoms ( self ) :

        self.UpdateModColor ( "ats" )
        # self.RandColorChains()

        ats = [at for at in self.cur_mol.atoms if not at.element.name == "H"]
        points = _multiscale.get_atom_coordinates ( ats, transformed = False )
        print " - search tree: %d/%d ats" % ( len(ats), len(self.cur_mol.atoms) )
        allAtTree = AdaptiveTree ( points.tolist(), ats, 2.0)

        doRess = self.GetCurRess()
        for r in doRess :
            for at in r.atoms :
                at.label = ""


        for at in self.cur_mol.atoms :
            at.label = ""

        if 1 :
            for at in chimera.selection.currentAtoms () :
                if at.display == True :
                    if 1 and hasattr (at, 'Q1') and hasattr (at, 'Q2') :
                        at.label = "(%.2f)" % ( (at.Q1+at.Q2)/2.0 )
                    elif hasattr (at, 'Q') :
                        at.label = "%.2f" % at.Q

                    at.labelColor = chimera.MaterialColor (0,0,0,1)

        else :
            doRess = chimera.selection.currentResidues()
            #if len(doRess) == 0 :
            #    doRess = self.GetCurRess()

            if len(doRess) > 0 :
                for r in doRess :
                    for at in r.atoms :
                        if at.display == True :
                            if 1 and hasattr (at, 'Q1') and hasattr (at, 'Q2') :
                                at.label = "(%.2f)" % ( (at.Q1+at.Q2)/2.0 )
                            elif hasattr (at, 'Q') :
                                at.label = "%.2f" % at.Q

                            at.labelColor = chimera.MaterialColor (0,0,0,1)
                            #at.labelOffset = chimera.Vector(0,0,0)

                    nats = self.AtsWithin ( r.atoms, 3.0, allAtTree )
                    for at in nats :
                        if at.display == True :
                            if 1 and hasattr (at, 'Q1') and hasattr (at, 'Q2') :
                                at.label = "(%.2f)" % ( (at.Q1+at.Q2)/2.0 )
                            elif hasattr (at, 'Q') :
                                at.label = "%.2f" % at.Q

                            at.labelColor = chimera.MaterialColor (0,0,0,1)
                            #at.labelOffset = chimera.Vector(0,0,0)



        # at.label, labelColor, labelCoord, labelOffset
        # at.label = "HI"
        # at.labelColor = chimera.MaterialColor (0,0,0,1)
        umsg ( "Labeled atoms" )



    def DoColorRandom ( self ) :

        self.RandColorChains ()


    def ExpandSelBonded ( self ) :

        print self.cur_mol.name, len(chimera.selection.currentAtoms())
        #print self.chain.get()

        bondedAtoms = {}
        for at in self.cur_mol.atoms :
            bondedAtoms[at] = []

        for b in self.cur_mol.bonds :
            at1, at2 = b.atoms
            bondedAtoms[at1].append ( at2 )
            bondedAtoms[at2].append ( at1 )

        visAts = {}
        Q = []
        for at in chimera.selection.currentAtoms() :
            visAts[at] = 1
            Q.append (at)

        while len(Q) > 0 :
            at = Q.pop(0)
            visAts[at] = 1
            if at in bondedAtoms :
                for at2 in bondedAtoms[at] :
                    if not at2 in visAts :
                        Q.append ( at2 )

        selAts = visAts.keys()
        umsg ( "Selected %d connected atoms" % len(selAts) )
        chimera.selection.clearCurrent()
        for at in selAts :
            chimera.selection.addCurrent ( at )


    def PlotQ ( self ) :

        print "https://stackoverflow.com/questions/64276513/draw-dotted-or-dashed-rectangle-from-pil"


    def UpdateSurfColor ( self ) :

        print " - surf of %s, by %s" % ( self.cur_dmap.name, self.cur_mol.name )

        numAt = 0
        for r in self.cur_mol.residues :
            for at in r.atoms :
                if at.element.name == "H" :
                    pass
                else :
                    numAt += 1

        allAtPos = numpy.zeros ( (numAt, 3) )
        allAts = [None] * numAt

        numAt = 0
        for r in self.cur_mol.residues :
            for at in r.atoms :
                if at.element.name == "H" :
                    pass
                else :
                    allAtPos[numAt] = at.coord().data()
                    allAts[numAt] = at
                    at.allPtI = numAt
                    numAt += 1


        print " - tree with %d ats" % numAt
        allAtTree = AdaptiveTree ( allAtPos.tolist(), allAts, 4.0)
        print " - done"





    def UpdateModColor ( self, colorMod ) :

        ress = []
        try :
            ress = self.seqRes
        except :
            pass

        if len ( ress ) == 0 :
            if self.cur_mol :
                ress = self.cur_mol.residues
            else :
                umsg ( "No molecule/chain selected?" )
                return

        if not hasattr (self, 'scores') :
            umsg ( "No scores... press Q, Qp, or Qf button first" )
            return

        foundScore = True
        for ri, r in enumerate ( ress ) :
            if r != None and not hasattr (r, 'Q') :
                foundScore = False

        if not foundScore :
            umsg ( "No scores... press Calc or Load button first" )
            return


        minScore, maxScore = 0,0
        if colorMod == "sc" :
            minScore, maxScore = self.minScore1, self.maxScore1
        else :
            minScore, maxScore = self.minScore2, self.maxScore2

        #cH = numpy.array( [0.19,0.53,0.87] ) # [0.0,1.0,0.0]
        #cL = numpy.array( [1.0,0.0,0.0] )

        cH = numpy.array ( [.33,.56,.88] )
        cL = numpy.array ( [.99,.99,.3] )

        if 0 : # qRedGreen :
            cH = numpy.array( [50.0/255.0,250.0/255.0,50.0/255.0] )
            cL = numpy.array( [250.0/255.0,50.0/255.0,50.0/255.0] )
        elif 1 : # qRedBlue :
            print " - red-blue"
            cH = numpy.array( [0.31, 0.46, 0.76] )
            #cH = numpy.array( [0.64, 0.92, 1.00] )
            cL = numpy.array( [0.99, 0.33, 0.33] )




        for ri, r in enumerate ( ress ) :
            sc = None
            if r == None :
                continue
            #sc = self.scores[ri] if colorSC else self.scores2[ri]
            if colorMod == "sc" :
                sc = r.qSC if hasattr (r, 'qSC') else 0
            elif colorMod == "bb" :
                sc = r.qBB if hasattr (r, 'qBB') else 0
            else :
                sc = r.Q if hasattr (r, 'Q') else 0

            if sc == None  :
                r.ribbonColor = chimera.MaterialColor ( .7, .7, .7, 1.0 )
                for at in r.atoms :
                    #at.color = r.ribbonColor
                    try :
                        at.color = atomColors[at.element.name.upper()]
                    except :
                        at.color = atomColors[' ']

            else :
                h = (sc - minScore) / (maxScore - minScore)
                if h > 1 : h = 1
                if h < 0 : h = 0
                c = h * cH + (1-h) * cL
                r.ribbonColor = chimera.MaterialColor ( c[0], c[1], c[2], 1.0 )
                for at in r.atoms :
                    #at.color = r.ribbonColor
                    try :
                        at.color = atomColors[at.element.name.upper()]
                    except :
                        at.color = atomColors[' ']





    def RandColorChains ( self ) :

        if self.cur_mol == None :
            umsg ("Select a molecule first")
            return

        m = self.cur_mol

        from random import random as rand

        ct = {}
        for r in m.residues: ct[r.id.chainId] = 1
        clist = ct.keys()
        clist.sort()
        chains_clrs = {}
        cnames = ""

        for ci, cid in enumerate ( clist ) :
            clr = ( rand()*.8+.1, rand()*.8+.1, rand()*.8+.1 )
            chains_clrs[cid] = chimera.MaterialColor ( clr[0], clr[1], clr[2], 1.0 )
            cnames = cnames + cid

        print "%s - color ribbon for %d chains -" % ( m.name, len(cnames) ), cnames

        # color atoms
        for r in m.residues :
            clr = chains_clrs[r.id.chainId]
            r.ribbonColor = clr
            for at in r.atoms :
                at.color = clr


    def AllChain ( self ) :

        if self.cur_mol == None :
            umsg ("Select a molecule first")
            return

        chainId = self.chain.get()
        if len(chainId) == 0 :
            umsg ("Select a chain first")
            return

        umsg ( "Showing mol %s chain %s" % (self.cur_mol.name, chainId) )

        qscores.SetBBAts ( self.cur_mol )
        #ct = {}
        #for r in self.cur_mol.residues: ct[r.id.chainId] = 1
        #clist = ct.keys()
        #clist.sort()

        for r in self.cur_mol.residues :
            if r.id.chainId == chainId :
                if r.isProt or r.isNA :
                    r.ribbonDisplay = True
                    r.ribbonDrawMode = 2
                else :
                    r.ribbonDisplay = False
                    for at in r.atoms :
                        at.drawMode = at.EndCap
                        at.display = True
            else :
                r.ribbonDisplay = False
                for at in r.atoms :
                    #at.drawMode = at.EndCap
                    at.display = False


    def ShowOnlySel ( self ) :
        if self.cur_mol == None :
            umsg ("Select a molecule first")
            return

        m = self.cur_mol

        rsel = chimera.selection.currentResidues ()

        if len(rsel) == 0 :
            umsg ("Show only selected residues - no residue found to be selected")
            return

        risel = {}
        for r in rsel :
            risel["%d.%d.%s" % (r.molecule.id, r.id.position, r.id.chainId)] = 1

        for r in m.residues :
            rid = "%d.%d.%s" % (m.id, r.id.position, r.id.chainId)
            if rid in risel :
                r.ribbonDisplay = not self.showingAtoms
                for at in r.atoms :
                    if at.element.name == "H" :
                        at.display = self.showH.get()
                    else :
                        at.display = True
            else :
                r.ribbonDisplay = False
                for at in r.atoms :
                    at.display = False



    def FindAlts ( self ) :
        if self.cur_mol == None :
            umsg ("Select a molecule first")
            return

        m = self.cur_mol

        atMap = {}
        for r in m.residues :

            hasAlt = False
            for at in r.atoms :
                if len(at.altLoc) > 0 :
                    hasAlt = True
                    break

            if hasAlt :
                r.ribbonDisplay = True
                for at in r.atoms :
                    if at.element.name == "H" :
                        at.display = False
                    else :
                        at.display = True
                        atMap[at] = 1
                        at.drawMode = at.EndCap
                        if at.element.name.upper() in atomColors :
                            at.color = atomColors[at.element.name.upper()]
                        else :
                            at.color = atomColors[" "]
            else :
                r.ribbonDisplay = True
                for at in r.atoms :
                    at.display = False

        for bond in m.bonds :
            #if bond.atoms[0] in atMap or bond.atoms[1] in atMap :
            bond.display = bond.Smart
            bond.drawMode = bond.Stick




    def DelAlts ( self ) :
        if self.cur_mol == None :
            umsg ("Select a molecule first")
            return

        m = self.cur_mol
        from random import random

        atMap = {}
        for r in m.residues :

            altScores = {}
            for at in r.atoms :
                if at.isSC :
                    alt = "_" if at.altLoc == '' else at.altLoc
                    if alt in altScores :
                        if hasattr ( at, 'Q' ) :
                            altScores[alt].append ( at.Q )
                        else :
                            altScores[alt].append ( random() )
                    else :
                        if hasattr ( at, 'Q' ) :
                            altScores[alt] = [at.Q]
                        else :
                            altScores[alt] = [random()]

            if len ( altScores.keys() ) > 1 :
                #print " - res %s %d.%s" % (r.type, r.id.position, r.id.chainId)
                keepAlt = ''
                maxScore = 0
                for alt, scores in altScores.iteritems() :
                    avg = numpy.mean(scores)
                    #print " %s: %.2f - %d" % (alt, avg, len(scores))
                    if avg > maxScore :
                        keepAlt = alt
                        maxScore = avg
                print " - %s %d.%s, keeping %s score %.2f" % (r.type, r.id.position, r.id.chainId, keepAlt, maxScore)

                for at in r.atoms :
                    if len(at.altLoc) > 0 :
                        if at.altLoc == keepAlt :
                            at.altLoc = ''
                        else :
                            m.deleteAtom ( at )




    def AllChains ( self ) :
        if self.cur_mol == None :
            umsg ("Select a molecule first")
            return

        m = self.cur_mol

        #ct = {}
        #for r in m.residues: ct[r.id.chainId] = 1
        #clist = ct.keys()
        #clist.sort()

        for r in m.residues :
            if ("CA" in r.atomsMap and "N" in r.atomsMap and "C" in r.atomsMap) or ("O3'" in r.atomsMap and "O5'" in r.atomsMap)  :
                r.ribbonDisplay = True
                r.ribbonDrawMode = 2
            else :
                r.ribbonDisplay = False
                for at in r.atoms :
                    #at.drawMode = at.Ball
                    at.display = True


    def GetCurRess ( self ) :
        if self.cur_mol == None :
            umsg ("Select a molecule first")
            return []

        chainId = self.chain.get()
        if len(chainId) == 0 :
            umsg ("Select a chain first")
            return []

        ress = []
        for r in self.cur_mol.residues :
            if r.id.chainId == chainId :
                    ress.append ( r )

        return ress


    def SetSelRibbon ( self ) :

        selRess = chimera.selection.currentResidues()
        if len(selRess) > 0 :
            self.SetDrawMode ( chimera.selection.currentResidues(), showRibbon = True )
        else :
            self.SetDrawMode ( self.GetCurRess(), showRibbon = True )

        self.showingAtoms = False


    def SetSelAtoms ( self ) :

        selRess = chimera.selection.currentResidues()
        if len(selRess) > 0 :
            self.SetDrawMode ( chimera.selection.currentResidues(), showRibbon = False )
            print " - on %d res" % len(selRess)
        else :
            self.SetDrawMode ( self.GetCurRess(), showRibbon = False )

        self.showingAtoms = True




    def SetDrawMode ( self, ress, showRibbon = None ) :

        #if showRibbon == None :
        #    showRibbon = segmod_dialog().showRibbon.get()

        #showRibbon = True

        #qscores.SetBBAts ( ress[0].molecule )

        for m in chimera.openModels.list() :
            if type(m) == chimera.Molecule :
                if not hasattr ( m, 'bbats' ) :
                    qscores.SetBBAts(m)
                    m.bbats = True


        #for at in ress[0].molecule.atoms :
        #    at.drawMode = at.EndCap
        #    at.display = False # not showRibbon

        #for res in ress[0].molecule.residues :
        #    res.ribbonDisplay = res.ribbonDisplay

        showH = self.showH.get()

        atMap = {}
        #atI = 0
        #c1 = (1.0,0.0,0.0,1)
        #c1 = (1.0,0.0,0.0,1)
        for res in ress :
            for at in res.atoms :

                if not hasattr (res, 'isProt') :
                    qscores.SetBBAts (res.molecule)

                if res.isProt or res.isNA :
                    at.drawMode = at.EndCap

                    if at.element.name == "H" :
                        at.display = showH
                    else :
                        at.display = True # not showRibbon

                    if at.element.name.upper() in atomColors :
                        at.color = atomColors[at.element.name.upper()]
                    else :
                        at.color = atomColors[" "]
                    atMap[at] = 1

                    res.ribbonDisplay, res.ribbonDrawMode = showRibbon, res.Ribbon_Round


            #f = float(atI) / float(len(ress)-1)
            #res.ribbonColor = chimera.MaterialColor( f*0.8+0.2, 0.02, (1-f)*0.8+0.2, 1.0 );
            #atI+=1

        for bond in ress[0].molecule.bonds :
            if bond.atoms[0] in atMap or bond.atoms[1] in atMap :
                bond.display = bond.Smart
                bond.drawMode = bond.Stick


    def Wire ( self ) :

        showH = self.showH.get()

        selRess = chimera.selection.currentResidues()
        if len(selRess) > 0 :

            atMap = {}
            for res in selRess :
                for at in res.atoms :
                    if res.isProt or res.isNA :
                        at.drawMode = at.EndCap
                        at.display = True # not showRibbon
                        if showH == False and at.element.name == "H" :
                            at.display = False
                        if at.element.name.upper() in atomColors :
                            at.color = atomColors[at.element.name.upper()]
                        else :
                            at.color = atomColors[" "]
                        atMap[at] = 1

                        res.ribbonDisplay, res.ribbonDrawMode = False, res.Ribbon_Round


            for bond in selRess[0].molecule.bonds :
                if bond.atoms[0] in atMap or bond.atoms[1] in atMap :
                    bond.display = bond.Smart
                    bond.drawMode = bond.Wire


    def ShowAts ( self ) :

        for mod in chimera.openModels.list() :
            if type(mod) == chimera.Molecule and mod.display == True :
                for res in mod.residues :
                    #if res.id.position in rs and res.id.chainId == cid :
                    if res.id.position in rs :
                        for at in res.atoms :
                            at.drawMode = at.EndCap
                            at.display = True
                            try :
                                at.color = atomColors[at.element.name.upper()]
                            except :
                                at.color = atomColors[" "]




    def HideSCs ( self ) :

        for mol in chimera.selection.currentMolecules() :
            if 1 or not hasattr ( mol, 'bbats' ) :
                qscores.SetBBAts(mol)
                mol.bbats = True

        ress = chimera.selection.currentResidues()
        if len(ress) == 0 :
            ress = self.GetCurRess()

        showH = self.showH.get()

        #print len(ress)

        for res in ress :
            #if res.id.position in rs and res.id.chainId == cid :
            for at in res.atoms :
                #at.drawMode = at.EndCap

                if not hasattr (at, 'isBB') :
                    qscores.SetBBAts (at.molecule)
                if not hasattr (at, 'isSugar') :
                    print "%s.%d.%s" % (at.name, at.residue.id.position, at.residue.type)

                wat = at
                if at.element.number == 1 :
                    if showH == False :
                        at.display = False
                        continue
                    wat = at.neighbors[0]

                if self.showingAtoms :
                    at.display = wat.isBB or wat.isSugar
                else :
                    at.display = wat.isBB

                #try :
                #    at.color = atomColors[at.element.name.upper()]
                #except :
                #    at.color = atomColors[" "]


    def ShowSCs ( self ) :

        for mol in chimera.selection.currentMolecules() :
            if 1 or not hasattr ( mol, 'bbats' ) :
                qscores.SetBBAts(mol)
                mol.bbats = True

        ress = chimera.selection.currentResidues()
        if len(ress) == 0 :
            ress = self.GetCurRess()

        showH = self.showH.get()

        for res in ress :
            for at in res.atoms :
                #at.drawMode = at.EndCap
                if at.element.name == "H" :
                    at.display = showH
                else :
                    at.display = True

                try :
                    at.color = atomColors[at.element.name.upper()]
                except :
                    at.color = atomColors[" "]


    def ShowNear ( self ) :

        selMol = {}
        if 1 :
            selMol = {self.cur_mol : 1}
            print " - using cur mol %s" % self.cur_mol.name
        elif 0 :
            for at in chimera.selection.currentAtoms() :
                selMol[at.molecule] = 1

        ats = []
        #for mol in chimera.selection.currentMolecules() :
        for mol in selMol.keys() :
            print mol.name
            ats.extend ( mol.atoms )
            for r in mol.residues :
                r.ribbonDisplay = False
                for at in r.atoms :
                    r.display = False

        print " - %d atoms " % len(ats)
        #ats = [at for at in self.cur_mol.atoms if not at.element.name == "H"]

        if 0 :
            points = _multiscale.get_atom_coordinates ( ats, transformed = False )
            print " - search tree: %d/%d ats" % ( len(ats), len(self.cur_mol.atoms) )
            start = time.time()
            allAtTree = AdaptiveTree ( points.tolist(), ats, 2.0)
            print " - tree: %.2f sec" % (time.time() - start)

        import gridm; reload ( gridm )
        start = time.time()
        agridm = gridm.Grid ()
        agridm.FromAtomsLocal ( ats, 5.0 )
        #print " - gridm: %.2f sec" % (time.time() - start)

        if 0 :
            start = time.time()
            agrid = grid.Grid ()
            agrid.FromAtomsLocal ( ats, 4.0 )
            print " - grid: %.2f sec" % (time.time() - start)


        #chimera.selection.clearCurrent ()

        start = time.time()

        for near_R in [5.0] :
            #print "within %.1fA" % R
            nearRes = {}
            nearAts = {}
            for at in chimera.selection.currentAtoms() :
                #selMol[at.molecule] = 1
                nats = agridm.AtsNearPtLocal ( at.coord() )
                for nat, v in nats :
                    if 1 or nat.residue.isProt or nat.residue.isNA :
                        if nat.name[0] == "H" :
                            continue
                        if nat.residue in nearRes :
                            if v.length < nearRes[nat.residue][0] :
                                nearRes[nat.residue] = [v.length,at.serialNumber]
                        else :
                            nearRes[nat.residue] = [v.length,at.serialNumber]
                        nearAts[nat] = 1


            #print "atoms:"
            #fp = open ( "/Users/greg/Desktop/_atoms_%.0fA.txt" % R, "w" )
            #for at in nearAts.keys() :
                #print "%s.%d.%s" % (at.name, at.residue.id.position, at.residue.id.chainId),
            #    fp.write ( "%s.%d.%s\n" % (at.name, at.residue.id.position, at.residue.id.chainId) )

            #print ""
            #fp.close()

            #print "residues: "
            #fp = open ( "/Users/greg/Desktop/_residues_%.0fA.txt" % R, "w" )
            for r in nearRes.keys() :
                #chimera.selection.mergeCurrent ( chimera.selection.EXTEND, chimera.selection.OSLSelection ("") )
                if r.molecule in selMol :
                    chimera.selection.addCurrent ( r )
                    #r.ribbonDisplay = False
                    # print "%s\t%d\t%s\t%d\t%f" % (r.id.chainId, r.id.position, r.type, nearRes[r][1], nearRes[r][0])
                    for at in r.atoms :
                        #at.drawMode = at.EndCap
                        at.display = True
                        if at.element.name.upper() in atomColors :
                            at.color = atomColors[at.element.name.upper()]


            #print ""
            #fp.close()

        #print " - grid: %.5f sec" % (time.time() - start)




    def ShowNearWI ( self ) :

        for mol in chimera.selection.currentMolecules() :
            if 1 or not hasattr ( mol, 'bbats' ) :
                qscores.SetBBAts(mol)
                mol.bbats = True

        ats = [at for at in self.cur_mol.atoms if not at.element.name == "H"]

        if 0 :
            points = _multiscale.get_atom_coordinates ( ats, transformed = False )
            print " - search tree: %d/%d ats" % ( len(ats), len(self.cur_mol.atoms) )
            start = time.time()
            allAtTree = AdaptiveTree ( points.tolist(), ats, 2.0)
            print " - tree: %.2f sec" % (time.time() - start)

        import gridm; reload ( gridm )
        start = time.time()
        agridm = gridm.Grid ()
        agridm.FromAtomsLocal ( ats, 4.0 )
        #print " - gridm: %.2f sec" % (time.time() - start)

        if 0 :
            start = time.time()
            agrid = grid.Grid ()
            agrid.FromAtomsLocal ( ats, 4.0 )
            print " - grid: %.2f sec" % (time.time() - start)

        #chimera.selection.clearCurrent ()
        start = time.time()

        for R in [4.0] :
            #print "within %.1fA" % R
            nearRes = {}
            nearAts = {}
            for at in chimera.selection.currentAtoms() :
                nats = agridm.AtsNearPtLocal ( at.coord() )
                for nat, v in nats :
                    if nat.residue.type.upper() in chargedIons or nat.residue.type.upper() == "HOH" :
                        nearRes[nat.residue] = 1
                        nearAts[nat] = 1

            for r in chimera.selection.currentResidues() :
                chimera.selection.addCurrent ( r )

            for r in nearRes.keys() :
                chimera.selection.addCurrent ( r )
                for at in r.atoms :
                    #at.drawMode = at.EndCap
                    at.display = True
                    if at.element.name.upper() in atomColors :
                        at.color = atomColors[at.element.name.upper()]




    def FindClashes ( self ) :

        if self.cur_mol == None :
            umsg ( "Select a molecule" )
            return

        for mol in chimera.openModels.list() :
            if type(mol) == chimera.Molecule :
                mol.display = False

        mol = self.cur_mol
        mol.display = True

        ats = [at for at in self.cur_mol.atoms if not at.element.name == "H"]

        import gridm; reload ( gridm )
        start = time.time()
        agridm = gridm.Grid ()
        agridm.FromAtomsLocal ( ats, 3.0 )
        #print " - gridm: %.2f sec" % (time.time() - start)

        for r in self.cur_mol.residues :
            if r.isProt or r.isNA :
                r.ribbonDisplay = False
            for at in r.atoms :
                at.display = False

        for ai, at in enumerate ( ats ) :
            nearRes = {}
            #nearAts = {}
            nats = agridm.AtsNearPtLocal ( at.coord() )
            bondedAts = at.bondsMap.keys()
            for nat, v in nats :
                if v.length < 1.8 and not nat in bondedAts and not at == nat :
                    #if nat.residue.type.upper() in chargedIons or nat.residue.type.upper() == "HOH" :
                    nearRes[nat.residue] = 1
                    #nearAts[nat] = 1
                    print " -x- %d.%s@%s,%d.%s@%s - %.2f" % (at.residue.id.position, at.residue.id.chainId, at.name, nat.residue.id.position, nat.residue.id.chainId, nat.name, v.length)


            if len(nearRes.keys()) > 0 :
                for r in [at.residue] + nearRes.keys() :
                    r.ribbonDisplay = True
                    for at in r.atoms :
                        #at.drawMode = at.EndCap
                        at.display = True
                        if at.element.name.upper() in atomColors :
                            at.color = atomColors[at.element.name.upper()]

            #if ai % 10000 == 0 :
            #    print "%d/%d" % (ai/10000, len(ats)/10000),

        print ""



    def Zone ( self ) :

        print "Zone:", self.zoneRad.get()

        try :
            rad = float ( self.zoneRad.get() )
        except :
            umsg ( "Enter a number for zone radius" )
            return

        atoms = chimera.selection.currentAtoms()
        if len(atoms) == 0 :
            umsg ( "Nothing selected" )
            return

        if self.cur_dmap == None :
            umsg ( "Select a Map" )
            return

        dmap = self.cur_dmap
        m = atoms[0].molecule

        from _multiscale import get_atom_coordinates
        points = get_atom_coordinates ( atoms, transformed = True )

        mods = {}
        for m in chimera.openModels.list() :
            mods[m.name] = m

        for i in range ( 10000 ) :
            nname = os.path.splitext(dmap.name)[0] + "_Z%.0f_%d" % (rad,i+1) + ".mrc"
            if not nname in mods :
                break

        cmap = self.PtsToMap ( points, dmap, rad, nname, showMesh=False, alpha=0.4 )
        #self.PtsToMap ( points, dmap, R, dmap.name + label, False, alpha=0.2 if self.showMesh.get() else 0.4 )

        umsg ( "Made zone map: " + nname )
        dmap.display = False

        chimera.runCommand ( "vol #%d style surface region all step 1" % cmap.id )



    def Inter ( self ) :

        for mol in chimera.selection.currentMolecules() :
            if 1 or not hasattr ( mol, 'bbats' ) :
                qscores.SetBBAts(mol)
                mol.bbats = True

        print ""
        print "Interactions for: %s, %d atoms" % ( self.cur_mol, len(self.cur_mol.atoms) )


        ats = [at for at in self.cur_mol.atoms if not at.element.name == "H"]
        points = _multiscale.get_atom_coordinates ( ats, transformed = False )
        print " - search tree: %d/%d ats" % ( len(ats), len(self.cur_mol.atoms) )
        allAtTree = AdaptiveTree ( points.tolist(), ats, 2.0)


        polar, hyd, wat, watm = {}, {}, {}, {}

        def setI_ (I, R1, R2) :
            if R1 in I :
                if R2 in I[R1] :
                    I[R1][R2] += 1
                else :
                    I[R1][R2] = 1
            else :
                I[R1] = {}
                I[R1][R2] = 1

        def setI (I, R1, R2) :
            setI_(I, R1, R2)
            setI_(I, R2, R1)


        def addI (at1, at2) :

            R1 = at1.residue.id.chainId
            R2 = at2.residue.id.chainId

            if (at1.element.name == "O" or at1.element.name == "N") and (at2.element.name == "O" or at2.element.name == "N") :
                setI ( polar, R1, R2 )
            else :
                setI ( hyd, R1, R2 )


        for at in self.cur_mol.atoms :

            nats = self.AtsWithin ( [at], 3.5, allAtTree )

            chains = {}
            if at.residue.type == "HOH" :
                for nat in nats :
                    if nat.residue.type == "HOH" :
                        continue
                    chains[nat.residue.id.chainId] = nat

            if len(chains.keys()) == 2 :
                c1, c2 = chains.keys()
                a1, a2 = chains[c1], chains[c2]

                if (a1.coord() - a2.coord()).length > 3.5 :
                    setI ( watm, chains.keys()[0], chains.keys()[1] )
                else :
                    setI ( wat, chains.keys()[0], chains.keys()[1] )

            if len(chains.keys()) > 2 :
                print "wat:", chains.keys()

            else :
                for nat in nats :
                    if nat.residue.type == "HOH" :
                        continue
                    if at.residue.id.chainId == nat.residue.id.chainId :
                        continue
                    addI ( at, nat )


        print "Polar: "
        print polar

        print "Hydrophobic: "
        print hyd

        print "Water: "
        print wat

        print "Water Mediated: "
        print watm




    def ShowDists ( self ) :

        m1 = [m for m in chimera.openModels.list() if m.display==True and type(m) == chimera.Molecule][0]
        print " - m1: %s" % m1.name

        avgB, avgN = 0.0, 0.0
        maxV = 0.0
        for r1 in m1.residues :
            avg, N = 0.0, 0.0
            for at in r1.atoms :
                avg += at.bfactor
                avgB += at.bfactor
                avgN += 1.0
                N += 1.0
            r1.avgB = avg / N
            if maxV < r1.avgB :
                maxV = r1.avgB

        print "avg B: %.3f" % (avgB/avgN)

        fp = open ( "/Users/greg/Desktop/bfactors.txt", "w" )
        for r1 in m1.residues :
            avg, N = 0.0, 0.0
            for at in r1.atoms :
                avg += at.bfactor
                N += 1.0
            fp.write (  "%d\t%.3f\t%.3f\n" % (r1.id.position, r1.avgB, (0.9-0.9*r1.avgB/maxV)) )
        fp.close()


    # per-residue distances
    def ShowDists_0 ( self ) :

        m1, m2 = [m for m in chimera.openModels.list() if m.display==True and type(m) == chimera.Molecule]
        print " - m1: %s" % m1.name
        print " - m2: %s" % m2.name

        qscores.SetBBAts ( m1 )
        qscores.SetBBAts ( m2 )

        rmap = {}
        for r in m1.residues :
            rmap[r.id.position] = r

        maxD = 0.0
        for r2 in m2.residues :
            r1 = rmap[r2.id.position]

            sumD, sumN = 0.0, 0.0
            for at2 in r2.atoms :
                at1 = r1.atomsMap[at2.name][0]
                v = at1.xformCoord() - at2.xformCoord()
                sumD += v.length
                sumN += 1.0

            r2.avgD = min ( sumD / sumN, 3.0 )
            if maxD < r2.avgD :
                maxD = r2.avgD

        fp = open ( "/Users/greg/Desktop/rdists.txt", "w" )
        for r2 in m2.residues :
            fp.write (  "%d\t%.3f\t%.3f\n" % (r2.id.position, r2.avgD, (1.0-r2.avgD/maxD)) )
            for at in r2.atoms :
                at.bfactor = r2.avgD
        fp.close()


    def ShowDists_0 ( self ) :

        m1, m2 = [m for m in chimera.openModels.list() if m.display==True and type(m) == chimera.Molecule]
        print " - m1: %s" % m1.name
        print " - m2: %s" % m2.name


        amap = {}
        for at in m2.atoms :
            atId = "%d.%s.%s.%s" % (at.residue.id.position,at.residue.id.chainId,at.name,at.altLoc)
            amap[atId] = at

        from chimera.resCode import protein3to1

        tt, tt2, nt = {}, {}, {}

        for at in m1.atoms :
            atId = "%d.%s.%s.%s" % (at.residue.id.position,at.residue.id.chainId,at.name,at.altLoc)
            if atId in amap :
                at2 = amap[atId]
                d = (at.coord()-at2.coord()).length
            else :
                print " - not found:", atId
                continue

            if at.display and not at.residue.type in protein3to1 :
                if not at.name in tt :
                    tt[at.name] = d; tt2[at.name] = d*d; nt[at.name] = 1.0
                else :
                    tt[at.name] += d; tt2[at.name] += d*d; nt[at.name] += 1.0

            if at.residue.type in protein3to1 :
                if at.isBB :
                    if not "BB" in tt :
                        tt["BB"] = d; tt2["BB"] = d*d; nt["BB"] = 1.0
                    else :
                        tt["BB"] += d; tt2["BB"] += d*d; nt["BB"] += 1.0
                else :
                    if not "SC" in tt :
                        tt["SC"] = d; tt2["SC"] = d*d; nt["SC"] = 1.0
                    else :
                        tt["SC"] += d; tt2["SC"] += d*d; nt["SC"] += 1.0


        for tp, D in tt.iteritems () :
            N, D2 = nt[tp], tt2[tp]
            rmsd = numpy.sqrt ( D2/N )
            avgd = D/N
            print "%s - %.0f atoms, avgd %.5f, rmsd: %.5f" % ( tp, N, avgd, rmsd )




    def ShowCh ( self, ch ) :

        if self.cur_mol == None :
            umsg ("Select a molecule first")
            return

        print " - showing chain:", ch

        qscores.SetBBAts ( self.cur_mol )

        m = self.cur_mol
        print " - cur mol:", m.name

        ct = {}
        for r in m.residues: ct[r.id.chainId] = 1
        clist = ct.keys()
        clist.sort()

        atsMap = {}
        for r in m.residues :
            show = True if r.id.chainId == ch else False

            if r.isProt or r.isNA :
                r.ribbonDisplay = show
                #r.ribbonDrawMode = 2
                for at in r.atoms :
                    at.display = False
            else :
                r.ribbonDisplay = False
                for at in r.atoms :
                    #at.drawMode = at.Ball
                    at.display = show
                    atsMap[at] = 1
                    if show :
                        at.drawMode = at.EndCap
                        if at.element.name.upper() in atomColors :
                            at.color = atomColors[at.element.name.upper()]
                        else :
                            at.color = atomColors[" "]
        for bond in m.bonds :
            #if bond.atoms[0] in atsMap or bond.atoms[1] in atsMap :
            bond.display = bond.Smart
            #else :
            #    bond.display = bond.Never




    def GetMod ( self, name ) :
        for m in chimera.openModels.list() :
            if name != None and len(name) > 0 :
                if m.name == name :
                    return m
            else :
                if m.display == True :
                    return m
        return None



    def GetSeq ( self ) :

        if self.cur_mol == None :
            umsg ( "No selected molecule" )
            return

        if len ( self.chain.get() ) == 0 :
            umsg ( "No selected chain" )
            return

        self.RemoveSeq ()
        print " - removed seq"

        try :
            print self.cur_mol.name
        except :
            print " - mol may have been closed"
            return

        self.GetSeqFromStruc ( self.cur_mol, self.chain.get() )

        if len(self.seq) > 0 :

            print "-- seq from open mol -- %d res" % len(self.seq)
            #print self.seq

            self.seqt = []
            self.seqSheetR = [None] * len(self.seq)
            self.seqHelixR = [None] * len(self.seq)
            self.seqScoreR = [None] * len(self.seq)
            self.seqScoreR2 = [None] * len(self.seq)
            self.scores2 = [None] * len(self.seq)
            self.scores = [None] * len(self.seq)

            self.UpdateSeqFont ()
            self.UpdateSeq ()

            return True

        return False



    def RemoveSeq  (self) :

        if self.seq == "" :
            return

        for si in range ( len(self.seq) ) :
            res = self.seq[si]
            pred = self.pred[si]
            conf = float ( self.conf[si] ) / 10.0

            if pred == 'E' :
                if self.seqSheetR[si] != None :
                    self.Canvas.delete ( self.seqSheetR[si] )

            elif pred == 'H' :
                if self.seqHelixR[si] != None :
                    self.Canvas.delete ( self.seqHelixR[si] )

            if self.seqScoreR[si] != None :
                self.Canvas.delete ( self.seqScoreR[si] )

            if self.seqScoreR2[si] != None :
                self.Canvas.delete ( self.seqScoreR2[si] )


        # box showing selected Residue
        if hasattr ( self, 'seqMouseR' ) :
            self.Canvas.delete ( self.seqMouseR )
            del self.seqMouseR

        if hasattr ( self, 'seqText' ) :
            self.Canvas.delete ( self.seqText )
            self.seqText = None
            del self.seqText

        self.seqSel = None
        self.seq = ""
        self.UpdateSeqSel ()



    def GetSeqFromStruc ( self, mol, chainId ) :

        print "Getting seq from %s, %s" % (mol.name, chainId)

        if self.showGaps.get() :
            print " - showing gaps"

        self.conf = ""
        self.pred = ""
        self.seq = ""
        self.seqRes = []
        self.seqRi = []

        if chainId == 'All' :
            return

        from chimera.resCode import protein3to1
        from chimera.resCode import nucleic3to1
        protein3to1['HSD'] = protein3to1['HIS']

        minri, maxri = None, None
        rids = {}
        for r in mol.residues :
            if r.id.chainId == chainId :
                if r.type in protein3to1 or r.type in nucleic3to1 :
                    rids[r.id.position] = r
                    if minri == None or r.id.position < minri : minri = r.id.position
                    if maxri == None or r.id.position > maxri : maxri = r.id.position


        ris = rids.keys()
        ris.sort()

        if maxri == None :
             return

        full_seq = ""
        for ri in range ( minri, maxri+1 ) :
            if ri in rids :
                r = rids[ri]
                if r.type in protein3to1 :
                    self.seq = self.seq + protein3to1[r.type]
                    self.conf = self.conf + "9"
                    predi = "C"
                    if r.isSheet : predi = "E"
                    if r.isHelix : predi = "H"
                    self.pred = self.pred + predi
                    self.seqRes.append ( r )
                    self.seqRi.append ( ri )
                    full_seq += protein3to1[r.type]
                elif r.type in nucleic3to1 :
                    self.seq = self.seq + nucleic3to1[r.type]
                    self.conf = self.conf + "9"
                    self.predi = "C"
                    self.pred = self.pred + self.predi
                    self.seqRes.append ( r )
                    self.seqRi.append ( ri )
                    full_seq += nucleic3to1[r.type]
            else :
                if self.showGaps.get() :
                    self.seq = self.seq + "."
                    self.conf = self.conf + "9"
                    self.pred = self.pred + "C"
                    self.seqRes.append ( None )
                    self.seqRi.append ( ri )
                    full_seq += "X"

        #print "seq:", full_seq





    def SSE ( self ) :

        print "sse"
        #self.GetFromMol ( mod, chainId )


    def CurRes ( self ) :

        #self.GetFromMol ( mod, chainId )

        if self.cur_mol == None :
            umsg ( "No selected molecule" )
            return []

        if self.cur_dmap == None :
            umsg ( "No selected map" )
            return []

        if len ( self.chain.get() ) == 0 :
            umsg ( "No selected chain" )
            return []

        from chimera.resCode import protein3to1
        protein3to1['HSD'] = protein3to1['HIS']

        rids = {}
        for r in self.cur_mol.residues :
            if r.id.chainId == self.chain.get() :
                if r.type in protein3to1 :
                    rids[r.id.position] = r

        print " - %d residues" % len(rids.values())
        #return [ rids[6] ]
        return rids.values ()



    def CalcZScores ( self ) :

        ress = []
        try :
            ress = self.seqRes
        except :
            pass

        if len ( ress ) == 0 :
            umsg ( "No molecule/chain selected?" )
            return

        self.scores2 = [None] * len(self.seqRes)
        scoreI = 0

        if self.cur_dmap == None :
            status ( "Select or open a map..." )
            return

        if self.cur_mol == None :
            status ( "Select or open a model..." )
            return

        status ( "Getting secondary structure elements..." )

        resolution = 3.0 * self.cur_dmap.data.step[0]
        #resolution = 3.0
        umsg ( "Calculating backbone Z-scores..." )

        zscores2 = []

        if 1 : # old
            sses = SSEs ( self.seqRes )
            print " - ", len(sses), "sse for ", len(ress), "res"

            atI = 1
            numH, numE, num_ = 0, 0, 0

            for el in sses :
                si, ei, ss, elRess = el

                if atI % 10 == 0 :
                    status ( "BB scores: %d/%d" % (atI,len(sses) ) )
                atI += 1

                #if 1 or (startRes < 129 and endRes > 129) :
                startResI, endResI, sseType, ress = el
                #print " : %d-%d, %s, %d res" % (startResI, endResI, sseType, len(ress))

                zscore, ccs = zBB ( self.cur_mol, ress, resolution, self.cur_dmap )
                #print ss, si, "-", ei, zscore
                if zscore != None :
                    zscores2.append ( zscore )

                for r in elRess :
                    r.bbZ = zscore
                    self.scores2[scoreI] = zscore
                    scoreI += 1

        else :

            bbs = BBsegs ( self.seqRes )
            W = 5
            atRes = 0

            for bb in bbs :
                print "%d res, %d-%d" % (len(bb),bb[0].id.position,bb[-1].id.position)

                for ri, r in enumerate ( bb ) :
                    firstRi = max ( 0, ri-(W-1)/2 )
                    lastRi = min ( len(bb)-1, ri+(W-1)/2 )
                    ress = bb[firstRi:lastRi+1]
                    zscore, ccs = zBB ( self.cur_mol, ress, resolution, self.cur_dmap )

                    #print "  %d : %d - %d, %.3f" % (ri, firstRi, lastRi, zscore)
                    if atRes % 50 == 0 :
                        status ( "Backbone - residue %d/%d" % (atRes,len(self.seqRes) ) )
                        #print "%d/%d" % (atRes,len(self.seqRes))
                        print "."

                    atRes += 1

                    if zscore != None :
                        zscores2.append ( zscore )

                    r.bbZ = zscore
                    r.CCS = ccs
                    self.scores2[scoreI] = zscore
                    scoreI += 1


        #print zscores2

        print " - %d res, min %.2f max %.2f, avg %.2f" % (len(ress), min(zscores2), max(zscores2), numpy.average(zscores2) )
        self.avgScore2 = numpy.average ( zscores2 )

        doRes = []

        doAllResInMol = False

        if doAllResInMol :
            for res in self.cur_mol.residues :
                if "CA" in res.atomsMap and "N" in res.atomsMap and "C" in res.atomsMap :
                    doRes.append ( res )

            print "++++ added all %d res from %s ++++" % (len(doRes), self.cur_mol.name)

        else :
            for r in self.seqRes :
                try :
                    blah
                    ra = r.scZ
                except :
                    doRes.append ( r )



        #doRes = self.seqRes
        #doRes = self.CurRes()
        print " - need score for %d res" % len(doRes)

        umsg ( "Calculating Side Chains / Bases Z-scores..." )

        sczScores = []
        if len(doRes) > 0 :
            sczScores = CalcRotaZ ( self.cur_dmap, self.cur_mol, doRes )
            #avgA, stdA = numpy.average ( A ), numpy.std ( A )
            #umsg ( "Avg side chain Z-score: %.3f" % ( avgA ) )

        if not doAllResInMol :
            doRes = self.seqRes

        self.scores = [None] * len(doRes)
        for ri, r in enumerate ( doRes ) :
            self.scores[ri] = r.scZ

        scores = [x for x in self.scores if x is not None]

        self.minScore = min ( scores )
        self.maxScore = max ( scores )
        self.avgScore = numpy.average ( scores )

        print " - %d res, min %.2f max %.2f, avg %.2f" % (len(doRes),self.minScore,self.maxScore, self.avgScore)

        self.minScore1, self.maxScore1 = 0,2
        self.minScore2, self.maxScore2 = 0,4

        bbRes = numpy.power ( numpy.e, (self.avgScore2 - 8.0334) / -4.128 ) # y = -4.128ln(x) + 8.0334
        scRes = numpy.power ( numpy.e, (self.avgScore - 4.8261) / -3.097 ) # y = -3.097ln(x) + 4.8261
        #scRes = (self.avgScore2 - 3.507) / -0.721
        #bbRes = (self.avgScore - 6.1234) / -0.9191

        umsg ( "Average BB Z-score: %.2f (%.1fA), Average Side Chain Z-score: %.2f (%.1fA)" % (self.avgScore2, bbRes, self.avgScore, scRes) )

        self.UpdateSeq ()



        sByType = {}
        rByType = {}
        for r in doRes :
            if r.scZ != None :
                if not r.type in sByType :
                    rByType[r.type] = []
                    sByType[r.type] = []
                rByType[r.type].append ( [r.scZ, r] )
                sByType[r.type].append ( [r.scZ] )

        avgs = []
        for rtype, ra in sByType.iteritems () :
            avgs.append ( [numpy.average (ra), rtype] )

        from chimera.resCode import protein3to1
        from chimera.resCode import nucleic3to1
        avgs.sort ( reverse=True, key=lambda x: x[0] )


        #mpath, mname = os.path.split ( dmap.data.path )
        dname, dext = os.path.splitext ( self.cur_dmap.data.path )
        #mfname = os.path.split ( self.cur_mol.openedAs[0] )[-1]
        #mname, mext = os.path.splitext ( mfname )

        avgm, numt = {}, {}
        for avgScore, rtype in avgs :

            rscores = rByType[rtype]
            rscores.sort ( reverse=True, key=lambda x: x[0] )
            hr = rscores[0]
            R = hr[1]
            highestScore = hr[0]
            numRes = len(rscores)

            if R.isProt :
                print "%s\t%s\t%d\t%f\t%d\t.%s\t%f" % (rtype, protein3to1[rtype], numRes, avgScore, R.id.position, R.id.chainId, highestScore)
            else :
                print "%s\t%s\t%d\t%f\t%d\t.%s\t%f" % (rtype, nucleic3to1[rtype], numRes, avgScore, R.id.position, R.id.chainId, highestScore)

            avgm[rtype] = avgScore
            numt[rtype] = numRes


        ofname = "%s__%s__scz_rtype.txt" % (dname, self.cur_mol.name)
        print " -> ", ofname
        fp = open ( ofname, "w" )

        for rt in ["PHE", "PRO", "ILE", "LEU", "VAL"] : # , "GLY", , "ALA"
            fp.write ( "%s\t%d\t%f\n" % (rt, numt[rt], avgm[rt]) )

        for rt in ["MET"] :
            fp.write ( "%s\t%d\t%f\n" % (rt, numt[rt], avgm[rt]) )

        for rt in ["HIS", "ARG", "LYS", "TRP", "CYS"] : #
            try :
                fp.write ( "%s\t%d\t%f\n" % (rt, numt[rt], avgm[rt]) )
            except :
                print " - no %s" % rt

        for rt in ["GLN", "ASN", "THR"] :
            fp.write ( "%s\t%d\t%f\n" % (rt, numt[rt], avgm[rt]) )

        for rt in ["TYR", "GLU", "ASP", "SER"] :
            fp.write ( "%s\t%d\t%f\n" % (rt, numt[rt], avgm[rt]) )

        fp.close()





    def CalcAllR (self) :

        ress = []
        try :
            ress = self.seqRes
        except :
            pass

        if len ( ress ) == 0 :
            umsg ( "No molecule/chain selected?" )
            return

        if self.cur_dmap == None :
            status ( "Select or open a map..." )
            return

        if self.cur_mol == None :
            status ( "Select or open a model..." )
            return

        cid = self.chain.get()

        CalcSCBBr ( self.cur_mol, cid, self.cur_dmap )


        self.scores, self.scores2 = [], []
        scBB, scSC = [], []

        for r in self.cur_mol.residues :
            if cid == None or r.id.chainId == cid :
                self.scores2.append ( r.SCBBr )
                self.scores.append ( r.SCBBr )
                r.scZ = r.SCBBr
                r.bbZ = r.SCBBr
                if r.SCBBr != None :
                    scBB.append ( r.SCBBr )
                if r.SCBBr != None :
                    scSC.append ( r.SCBBr )


        scMin, scMax, scAvg = min(scSC), max(scSC), numpy.average(scSC)
        bbMin, bbMax, bbAvg = min(scBB), max(scBB), numpy.average(scBB)


        print "Average R sc : %.2f - %.2f, avg %.2f" % (scMin, scMax, scAvg)
        print "Average R bb : %.2f - %.2f, avg %.2f" % (bbMin, bbMax, bbAvg)



        self.minSCscore, self.maxSCscore = 0.0,1
        self.minBBscore, self.maxBBscore = 0.0,1

        self.UpdateSeq ()



    def CalcAllSigma (self) :

        ress = []
        try :
            ress = self.seqRes
        except :
            pass

        if len ( ress ) == 0 :
            umsg ( "No molecule/chain selected?" )
            return

        if self.cur_dmap == None :
            status ( "Select or open a map..." )
            return

        if self.cur_mol == None :
            status ( "Select or open a model..." )
            return


        cid = self.chain.get()

        ats = [at for at in self.cur_mol.atoms if not at.element.name == "H"]
        points = _multiscale.get_atom_coordinates ( ats, transformed = False )
        print " - search tree: %d/%d ats" % ( len(ats), len(self.cur_mol.atoms) )
        allAtTree = AdaptiveTree ( points.tolist(), ats, 1.0)


        scBB, scSC = [], []

        for r in self.cur_mol.residues :
            if cid == None or r.id.chainId == cid :
                r.scores2 = r.bbZ
                r.scores1 = r.scZ
                if r.bbZ != None : scBB.append ( r.bbZ )
                if r.scZ != None : scSC.append ( r.scZ )


        #bbRes = numpy.power ( numpy.e, (self.avgScore2 - 8.0334) / -4.128 ) # y = -4.128ln(x) + 8.0334
        #scRes = numpy.power ( numpy.e, (self.avgScore - 4.8261) / -3.097 ) # y = -3.097ln(x) + 4.8261
        #scRes = (self.avgScore2 - 3.507) / -0.721
        #bbRes = (self.avgScore - 6.1234) / -0.9191

        scMin, scMax, scAvg = min(scSC), max(scSC), numpy.average(scSC)
        bbMin, bbMax, bbAvg = min(scBB), max(scBB), numpy.average(scBB)


        print "Average Sigma sc : %.2f - %.2f, avg %.2f | %.2f - %.2f, avg %.2f" % (scMin, scMax, scAvg, 1.0/scMin, 1.0/scMax, 1.0/scAvg)
        print "Average Sigma bb : %.2f - %.2f, avg %.2f | %.2f - %.2f, avg %.2f" % (bbMin, bbMax, bbAvg, 1.0/bbMin, 1.0/bbMax, 1.0/bbAvg)


        self.minScore1, self.maxScore1 = 0.0,0.5
        self.minScore2, self.maxScore2 = 0.0,0.2

        self.UpdateSeq ()




    def CalcAllQ (self) :

        ress = []
        try :
            ress = self.seqRes
        except :
            pass

        if len ( ress ) == 0 :
            umsg ( "No molecule/chain selected?" )
            #return

        if self.cur_dmap == None :
            status ( "Select or open a map..." )
            return

        if self.cur_mol == None :
            status ( "Select or open a model..." )
            return

        cid = self.chain.get()

        umsg ( "Calculating Q-scores - see bottom of main window for status or to cancel..." )

        res = float(self.mapRes.get())
        sigma = float ( self.sigma.get() )

        if sigma < 1.0 :
            Qavg = qscores.CalcQ (self.cur_mol, self.chain.get(), self.cur_dmap, res, sigma, log=True )
        else :
            Qavg = qscores.CalcQ (self.cur_mol, self.chain.get(), self.cur_dmap, res, sigma, log=True )

        qscores.SaveQStats ( self.cur_mol, self.chain.get(), self.cur_dmap.name, sigma, res )
        self.ShowQScores ()

        #umsg ( "Average Q-score for %s: %.2f" % (self.cur_mol.name, Qavg) )
        umsg ( "Done Q-scores for %s" % (self.cur_mol.name) )




    def CalcAllQp (self, numProc=None) :

        if self.cur_dmap == None :
            status ( "Select or open a map..." )
            return

        if self.cur_mol == None :
            status ( "Select or open a model..." )
            return

        cid = self.chain.get()

        if cid == "All" :
            cid = None

        #ats = [at for at in self.cur_mol.atoms if not at.element.name == "H"]
        #points = _multiscale.get_atom_coordinates ( ats, transformed = False )
        #print " - search tree: %d/%d ats" % ( len(ats), len(self.cur_mol.atoms) )
        #allAtTree = AdaptiveTree ( points.tolist(), ats, 1.0)

        if 0 :
            for r in self.cur_mol.residues :
                if hasattr ( r, 'Q' ) : del r.Q
                if hasattr ( r, 'qSC' ) : del r.qSC
                if hasattr ( r, 'qBB' ) : del r.qBB

        if len(self.cur_dmap.data.path) == 0 :
            umsg ( "No file for map - %s - must be saved first..." % (self.cur_dmap.name) )
            return

        if not os.path.isfile ( self.cur_dmap.data.path ) :
            umsg ( "Map file not found - %s - must be saved first..." % (self.cur_dmap.data.path) )
            return

        if not hasattr (self.cur_mol, 'openedAs') :
            umsg ( "No file for model %s - must be saved first..." % (self.cur_mol.name) )
            return

        if not os.path.isfile ( self.cur_mol.openedAs[0] ) :
            umsg ( "Model file not found - %s - must be saved first..." % (self.cur_mol.openedAs[0]) )
            return

        sigma = float(self.sigma.get())

        qscores.CalcQpn (self.cur_mol, cid, self.cur_dmap, sigma, numProc=numProc )
        qscores.SaveQStats ( self.cur_mol, self.chain.get(), self.cur_dmap.name, sigma, float(self.mapRes.get()) )

        self.ShowQScores ()



    def ShowQScores (self) :

        cid = self.chain.get()
        bbScores, scScores = [], []
        hasResCC = False
        for r in self.cur_mol.residues :
            #if cid == None or cid == "All" or r.id.chainId == cid :
            if r.id.chainId == cid :
                if hasattr ( r, "resCC" ) and r.resCC != None :
                    hasResCC = True
                if not hasattr ( r, 'Q2' ) or r.Q2 == None :
                    qscores.CalcResQ ( r )
                    if r.isProt or r.isNA :
                        r.score1 = r.qSC
                        r.score2 = r.qBB
                        if r.qSC != None : scScores.append ( r.qSC )
                        if r.qBB != None : bbScores.append ( r.qBB )
                    else :
                        r.score1 = r.Q
                        r.score2 = r.Q
                        if r.Q != None :
                            scScores.append ( r.Q )
                            bbScores.append ( r.Q )
                else :
                    if r.isProt :
                        r.score1 = None
                        r.score2 = r.Q2
                        if r.Q2 != None : scScores.append ( 0.0 )
                        if r.Q2 != None : bbScores.append ( r.Q2 )
                    else :
                        r.score1 = None
                        r.score2 = r.Q2
                        if r.Q2 != None :
                            scScores.append ( 0.0 )
                            bbScores.append ( r.Q2 )


        #bbRes = numpy.power ( numpy.e, (self.avgScore2 - 8.0334) / -4.128 ) # y = -4.128ln(x) + 8.0334
        #scRes = numpy.power ( numpy.e, (self.avgScore - 4.8261) / -3.097 ) # y = -3.097ln(x) + 4.8261
        #scRes = (self.avgScore2 - 3.507) / -0.721
        #bbRes = (self.avgScore - 6.1234) / -0.9191


        #try :
        if len(scScores) > 0 and len(bbScores) > 0 :
            scMin, scMax, scAvg = min(scScores), max(scScores), numpy.average(scScores)
            bbMin, bbMax, bbAvg = min(bbScores), max(bbScores), numpy.average(bbScores)
            print "Average Q sc : %.2f - %.2f, avg %.2f" % (scMin, scMax, scAvg )
            print "Average Q bb : %.2f - %.2f, avg %.2f" % (bbMin, bbMax, bbAvg )


        self.GetMaxScores()


        RES = float(self.mapRes.get())
        sigma = float(self.sigma.get())
        expQScore, lowQ, highQ, eqn = qscores.ExpectedQScore ( RES, sigma )
        #print " - res %.2f - expected Q-score: %.2f" % (RES, expQScore)

        if not hasResCC :
            self.minScore1, self.maxScore1 = 0.0,expQScore
            self.minScore2, self.maxScore2 = 0.0,expQScore
        else :
            self.minScore1, self.maxScore1 = 0.0,1.0
            self.minScore2, self.maxScore2 = 0.0,1.0

        #except :
        #    pass

        self.UpdateSeq ()



    def QuickQ (self) :

        if self.cur_dmap == None :
            status ( "Select or open a map..." )
            return

        if self.cur_mol == None :
            status ( "Select or open a model..." )
            return

        cid = self.chain.get()

        if cid == "All" :
            cid = None

        #ats = [at for at in self.cur_mol.atoms if not at.element.name == "H"]
        #points = _multiscale.get_atom_coordinates ( ats, transformed = False )
        #print " - search tree: %d/%d ats" % ( len(ats), len(self.cur_mol.atoms) )
        #allAtTree = AdaptiveTree ( points.tolist(), ats, 1.0)

        if 0 :
            for r in self.cur_mol.residues :
                if hasattr ( r, 'Q' ) : del r.Q
                if hasattr ( r, 'qSC' ) : del r.qSC
                if hasattr ( r, 'qBB' ) : del r.qBB


        sigma = float(self.sigma.get())
        CalcQp (self.cur_mol, cid, self.cur_dmap, sigma)



        scBB, scSC = [], []

        for r in self.cur_mol.residues :
            if cid == None or r.id.chainId == cid :
                if r.isProt or r.isNA :
                    r.score1 = r.qSC
                    r.score2 = r.qBB
                    if r.qBB != None : scBB.append ( r.qBB )
                    if r.qSC != None : scSC.append ( r.qSC )
                else :
                    r.score1 = r.Q
                    r.score2 = r.Q

        scMin, scMax, scAvg = min(scSC), max(scSC), numpy.average(scSC)
        bbMin, bbMax, bbAvg = min(scBB), max(scBB), numpy.average(scBB)


        print " - Average Q sc : %.2f - %.2f, avg %.2f" % (scMin, scMax, scAvg )
        print " - Average Q bb : %.2f - %.2f, avg %.2f" % (bbMin, bbMax, bbAvg )


        self.minScore1, self.maxScore1 = 0.0,1
        self.minScore2, self.maxScore2 = 0.0,1


        sigma = float(self.sigma.get())
        self.UpdateSeq ()
        qscores.SaveQStats ( self.cur_mol, self.chain.get(), self.cur_dmap.name, sigma, float(self.mapRes.get()) )





    def GetQsFromFile (self) :

        if self.cur_dmap == None :
            status ( "Select or open a map..." )
            return

        if self.cur_mol == None :
            status ( "Select or open a model..." )
            return

        chainId = self.chain.get()
        umsg ( "Loading Q-scores for chain %s..." % chainId )

        #molPath, molExt = os.path.splitext(self.cur_mol.openedAs[0])
        #mapName = os.path.splitext(self.cur_dmap.name)[0]

        nname, nname2 = qscores.QScoreFileName ( self.cur_mol, self.cur_dmap )
        print " -Q - ", nname
        print " -Q2- ", nname2
        if os.path.isfile ( nname ) :
            if hasattr ( self.cur_mol, 'cif' ) :
                qscores.QsFromCifFile ( self.cur_mol, nname )
            else :
                qscores.QsFromPdbFile ( self.cur_mol, nname )
        elif os.path.isfile ( nname2 ) :
            if hasattr ( self.cur_mol, 'cif' ) :
                qscores.QsFromCifFile ( self.cur_mol, nname2, Q2=True )
            else :
                qscores.QsFromPdbFile ( self.cur_mol, nname2, Q2=True )
        else :
            #print "-? -", nname
            umsg ( "Q scores not found for this map and file ---" )
            return


        if 0 :
            umsg ( "Saving files with Q-score B-factor" )
            self.SaveQsBfs ( self.cur_mol, 50.0 )
            self.SaveQsBfs ( self.cur_mol, 100.0 )
            self.SaveQsBfs ( self.cur_mol, 150.0 )
            self.SaveQsBfs ( self.cur_mol, 200.0 )
            self.SaveQsBfs ( self.cur_mol, 300.0 )

        umsg ( "Saving stats files for chain %s..." % chainId )
        res = float(self.mapRes.get())
        sigma = float(self.sigma.get())
        qscores.SaveQStats ( self.cur_mol, chainId, self.cur_dmap.name, sigma, res )

        #qscores.QStatsProt ( self.cur_mol, self.cur_dmap, chainId )
        #qscores.QStatsRNA ( self.cur_mol, self.cur_dmap, chainId )
        #qscores.QStats1 (self.cur_mol, chainId)
        avgq = qscores.QStats1 (self.cur_mol, res=res, sigma=sigma)
        print " - overall average Q: %.3f" % avgq

        umsg ( "Showing Q-scores for chain %s" % chainId )
        self.ShowQScores ()


    def SaveQsBfs ( self, mol, f ) :

        res = float ( self.mapRes.get() )
        sigma = float ( self.sigma.get() )
        expQ, lowQ, highQ, expF = qscores.ExpectedQScore (res, sigma)

        print ""
        print "B-factors -- expected Q at %.2fA, sigma %.2f: %.3f" % (res, sigma, expQ)
        print expF

        maxB = 0.0
        for at in mol.atoms :
            if not at.element.name == "H" :
                if hasattr (at, 'Q') :
                    #at.bfactor = max ( 1.0, f * max(0.0, (expQ - at.Q)) )
                    at.bfactor = f * max(0.0, (1.0 - at.Q))
                    maxB = max ( at.bfactor, maxB )
                else :
                    return None

        bondedAts = {}
        for b in mol.bonds :
            bondedAts[b.atoms[0]] = b.atoms[1]
            bondedAts[b.atoms[1]] = b.atoms[0]

        for at in mol.atoms :
            if at.element.name == "H" :
                bat = bondedAts[at]
                at.bfactor = bat.bfactor

        molPath = os.path.splitext(mol.openedAs[0])[0]

        nname = molPath + "__Bf%.0f__.pdb" % f
        print " - saving %s - max B: %.3f" % (nname, maxB)
        chimera.PDBio().writePDBfile ( [mol], nname )

        return nname


    def PhenixCC ( self ) :

        print "cc"
        if self.cur_dmap == None :
            status ( "Select or open a map..." )
            return

        dmap = self.cur_dmap

        if self.cur_mol == None :
            status ( "Select or open a model..." )
            return

        try :
            RES = float(self.mapRes.get())
        except :
            umsg ( "Please enter a numeric value for Resolution in Options" )


        umsg ( "Calculating Phenix-CC for %s in %s" % (self.cur_mol.name, self.cur_dmap.name) )

        print " - map res: %.2f" % RES

        self.phPath = "/Users/greg/_mol/phenix-installer-1.21.1-5286-mac-intel-osx-x86_64/build/bin/"
        self.phPath = "/Users/greg/_mol/phenix-installer-1.19.2-4158-mac-intel-osx-x86_64/build/bin/"
        self.phPath = "/Users/greg/_mol/phenix-1.20.1-4487/build/bin/"
        import subprocess

        print ""
        molPath = self.cur_mol.openedAs[0]

        mapf = molPath + "_fmodel.ccp4"
        if not os.path.exists(mapf) :
            self.MakePhMap ( molPath, RES )

        print " - loading map:", molPath + "_fmodel.ccp4"
        dm = chimera.openModels.open ( molPath + "_fmodel.ccp4" )[0]

        fpoints, fpoint_weights = fit_points_g ( dm.data, 0.01 )
        map_values = dmap.interpolated_values ( fpoints, dm.openState.xform )
        ov, cc, ccm = FitMap.overlap_and_correlation ( fpoint_weights, map_values )
        print " - olap: %.3f, CC: %.3f, CCm: %.3f" % (ov, cc, ccm)

        #chimera.openModels.close ( [dm] )

        return ccm, molPath




    def CalcBFactors ( self ) :

        if self.cur_dmap == None :
            status ( "Select or open a map..." )
            return

        if self.cur_mol == None :
            status ( "Select or open a model..." )
            return

        umsg ( "Calculating B-factors for %s in %s" % (self.cur_mol.name, self.cur_dmap.name) )

        try :
            RES = float(self.mapRes.get())
        except :
            umsg ( "Please enter a numeric value for Resolution in Options" )

        print " - map res: %.2f" % RES

        maxcc = -1.0
        maxf = 0.0
        self.phPath = "/Users/greg/_mol/phenix-installer-1.21.1-5286-mac-intel-osx-x86_64/build/bin/"
        self.phPath = "/Users/greg/_mol/phenix-installer-1.19.2-4158-mac-intel-osx-x86_64/build/bin/"
        self.phPath = "/Users/greg/_mol/phenix-1.20.1-4487/build/bin/"

        import subprocess

        molPaths = []
        #min, max = 0, 600
        #fs = range ( 0, 301, 20 )
        fs = [200]

        for f in fs :
            print " ---------------------------- %d -> 300 --------------------------------" % f
            ov, cc, ccm, ov2, cc2, ccm2, molPath = self.QB ( self.cur_mol, self.cur_dmap, f, RES)
            molPaths.append ( [f, cc, molPath] )

            fp = open ( os.path.splitext(self.cur_dmap.data.path)[0] + "_f_cc.txt", "a"  )
            fp.write ( "%f\t%f\t%f\t%f\t\t%f\t%f\t%f\t%f\n" % (f, ov, cc, ccm, f, ov2, cc2, ccm2) )
            fp.close()

            if ccm > maxcc :
                maxcc = ccm
                maxf = f

            #break

        print " - maxf: %.3f" % maxf

        if 0 :
            for f, cc, molPath in molPaths :
                print "%d\t%.4f" % (f,cc),
                if abs(f-maxf) > 0.1 :
                    print " - removing %s" % molPath
                    os.remove ( molPath)
                    #print " - removing %s" % (molPath + "_fmodel.ccp4")
                    os.remove ( molPath + "_fmodel.ccp4" )
                else :
                    print " - max"


    def CalcBFactors_BS ( self ) :

        if self.cur_dmap == None :
            status ( "Select or open a map..." )
            return

        if self.cur_mol == None :
            status ( "Select or open a model..." )
            return

        umsg ( "Calculating B-factors for %s in %s" % (self.cur_mol.name, self.cur_dmap.name) )

        try :
            RES = float(self.mapRes.get())
        except :
            umsg ( "Please enter a numeric value for Resolution in Options" )

        print " - map res: %.2f" % RES

        maxcc = -1.0
        maxf = 0.0
        self.phPath = "/Users/greg/_mol/phenix-installer-1.21.1-5286-mac-intel-osx-x86_64/build/bin/"
        self.phPath = "/Users/greg/_mol/phenix-installer-1.19.2-4158-mac-intel-osx-x86_64/build/bin/"
        self.phPath = "/Users/greg/_mol/phenix-1.20.1-4487/build/bin/"

        molPaths = []

        lowF = 50
        lowCC, molPath = self.QB ( self.cur_mol, self.cur_dmap, lowF, RES )
        molPaths.append ( [lowF, lowCC, molPath] )
        print "%.2f -> CC %.2f" % (lowF, lowCC)

        highF = 500
        highCC, molPath = self.QB ( self.cur_mol, self.cur_dmap, highF, RES )
        molPaths.append ( [highF, highCC, molPath] )
        print "%.2f -> CC %.2f" % (highF, highCC)

        while 1 :

            midF = (lowF + highF) / 2.0
            print "\n", midF, " -> ",
            midF = numpy.round ( midF/10.0 ) * 10.0
            print midF, " diff low:", abs(midF - lowF), "diff high:", abs(highF - midF)

            if abs(midF - lowF) < 10.01 or abs(highF - midF) < 10.01 :
                break

            midCC, molPath = self.QB ( self.cur_mol, self.cur_dmap, midF, RES )
            molPaths.append ( [midF, midCC, molPath] )

            #break


        print ""
        print ""
        print " - maxf: %.3f" % maxf

        for f, cc, molPath in molPaths :
            print f, "%.4f" % cc,
            if abs(f-maxf) > 0.1 :
                print " - removing %s" % molPath
                #os.remove ( molPath)
                #print " - removing %s" % (molPath + "_fmodel.ccp4")
                #os.remove ( molPath + "_fmodel.ccp4" )
            else :
                print " - max"


    def QB ( self, mol, dmap, f, RES ) :

            print ""
            umsg ( " __ f=%.0f __" % f )
            molPath = self.SaveQsBfs ( mol, f )
            if molPath == None :
                umsg ( "Make sure Q-scores have been calculated for all atoms" )
                return None

            if 1 :
                mapf = molPath + "_fmodel.ccp4"
                if not os.path.exists(mapf) :
                    self.MakePhMap ( molPath, RES )
                    print " - making %s" % mapf
                else :
                    print " - found _fmodel %s" % mapf

                #return


                print " - loading map:", molPath + "_fmodel.ccp4"
                dm = chimera.openModels.open ( molPath + "_fmodel.ccp4" )[0]

                # whole map to map
                fpoints, fpoint_weights = fit_points_g ( dm.data, 0.01 )
                map_values = dmap.interpolated_values ( fpoints, dm.openState.xform )
                ov, cc, ccm = FitMap.overlap_and_correlation ( fpoint_weights, map_values )
                print " - map / olap : %.3f, CC: %.3f, CCm: %.3f" % (ov, cc, ccm)

                # just at atom positions
                atoms = [at for at in mol.atoms if not at.element.name == "H"]
                xyz = _multiscale.get_atom_coordinates(atoms, transformed = False)

                phMapVals = dm.interpolated_values ( xyz, dm.openState.xform )
                mean, stdev = numpy.average (phMapVals), numpy.std (phMapVals)
                #print phMapVals
                #phMapVals = (phMapVals - mean)/stdev
                #print " - mean %f stdev %f" % (mean, stdev)
                #print phMapVals

                cmMapVals = dmap.interpolated_values ( xyz, dmap.openState.xform )
                mean, stdev = numpy.average (cmMapVals), numpy.std (cmMapVals)
                #cmMapVals = (cmMapVals - mean)/stdev
                #print " - mean %f stdev %f" % (mean, stdev)

                ov2, cc2, ccm2 = FitMap.overlap_and_correlation ( phMapVals, cmMapVals )
                print " - ats / olap : %.3f, CC: %.3f, CCm: %.3f" % (ov2, cc2, ccm2)

                chimera.openModels.close ( [dm] )

                return ov, cc, ccm, ov2, cc2, ccm2, molPath

            return 0, molPath


    def MakePhMap ( self, molPath, RES ) :

        print " -> fmodel"
        #args = [self.phPath+'phenix.fmodel', "high_resolution=%.1f"%RES, "scattering_table=electron", "generate_fake_p1_symmetry=True", molPath ]
        #args = [self.phPath+'phenix.fmodel', "high_resolution=%.1f"%RES, "low_resolution=5.0", "scattering_table=n_gaussian", "generate_fake_p1_symmetry=True", molPath ]
        args = [self.phPath+'phenix.fmodel', "high_resolution=%.1f"%RES, "scattering_table=electron", "generate_fake_p1_symmetry=True", molPath ]
        print " - : ", "] [".join ( args )
        fpath = os.path.split ( molPath )[0]
        #print " - in %s" % fpath

        logfp1 = os.path.splitext ( molPath )[0] + '_fmodel.log'; logf = open ( logfp1, "w" )
        logfp1e = os.path.splitext ( molPath )[0] + '_fmodel_err.log'; logfe = open ( logfp1e, "w" )
        from subprocess import Popen
        p = Popen(args, stdout=logf, stderr=logfe, cwd=fpath)
        p.wait()
        logf.close()
        logfe.close()

        print " -> mtz2map"
        args = [self.phPath+'phenix.mtz2map', "high_resolution=%.1f"%RES, "include_fmodel=true", "scattering_table=electron", molPath, molPath + ".mtz" ]
        #print " - : ", "] [".join ( args )
        logfp2 = os.path.splitext ( molPath )[0] + '_mtz2map.log'
        logf = open ( logfp2, "w" )
        p = Popen(args, stdout=logf, cwd=fpath)
        p.wait()
        logf.close()

        os.remove (logfp1)
        os.remove (logfp1e)
        os.remove (logfp2)
        os.remove (molPath + ".mtz")


    def SA_Q (self) :

        ress = []
        try :
            ress = self.seqRes
        except :
            pass

        if len ( ress ) == 0 :
            umsg ( "No molecule/chain selected?" )
            return

        if self.cur_dmap == None :
            status ( "Select or open a map..." )
            return

        if self.cur_mol == None :
            status ( "Select or open a model..." )
            return


        chainId = self.chain.get()



        #ats = [at for at in self.cur_mol.atoms if not at.element.name == "H"]
        #points = _multiscale.get_atom_coordinates ( ats, transformed = False )
        #print " - search tree: %d/%d ats" % ( len(ats), len(self.cur_mol.atoms) )
        #allAtTree = AdaptiveTree ( points.tolist(), ats, 1.0)

        umsg ( "Solvent Accessibility vs. Q... making surface for %d atoms..." % len(self.cur_mol.atoms) )
        print ".",


        # https://en.m.wikipedia.org/wiki/Van_der_Waals_radius
        vdwRadii = { 'H' : 1.2, # (1.09)[1]
                    'C' : 1.7,
                    'N' : 1.55,
                    'O' : 1.52,
                    'F' : 1.47,
                    'P' : 1.8,
                    'S' :  1.8   }

        #vdwRadii = { 'H' : 1.5, 'C' : 1.5, 'N' : 1.5, 'O' : 1.5, 'F' : 1.5, 'P' : 1.5, 'S' :  1.5 }



        if GetMod ( "Surface Pts" ) : chimera.openModels.close ( [GetMod ( "Surface Pts" )] )
        if GetMod ( "SA pts" ) : chimera.openModels.close ( [GetMod ( "SA pts" )] )
        if GetMod ( "SA- pts" ) : chimera.openModels.close ( [GetMod ( "SA pts" )] )
        if GetMod ( "ASP pts" ) : chimera.openModels.close ( [GetMod ( "SA pts" )] )


        surfPts = []
        for at in self.cur_mol.atoms :
            VWR = vdwRadii[at.element.name] if at.element.name in vdwRadii else 1.55
            apts = SpherePts ( at.coord(), VWR, 100 )
            #apts.extend ( SpherePts ( at.coord(), VWR/2.0, 50 ) )
            apts.append (  at.coord().data() )

            surfPts.extend ( apts )

            #AddSpherePts ( apts, atomColors2[at.element.name], 0.1, "Surface Pts" )
            #AddSpherePts ( apts, (.7,.7,.7,1), 0.1, "Surface Pts" )


        umsg ( "Solvent Accessibility vs. Q... making tree for %d atoms, %d points..." % (len(self.cur_mol.atoms), len(surfPts) ) )
        print ".",

        surfPtsTree = AdaptiveTree ( surfPts, surfPts, 2.0 )


        #AddSpherePts ( surfPts, atomColors2[at.element.name], 0.1, "Surface Pts" )



        molPath = os.path.splitext(self.cur_mol.openedAs[0])[0]
        mapName = os.path.splitext(self.cur_dmap.name)[0]

        nname = molPath + "__SA-Q__" + mapName + ".txt"
        fp = open ( nname, "w" )

        umsg ( "Solvent Accessibility vs. Q ... saving to file %s" % nname )
        print ".",

        doRess = []
        for r in self.cur_mol.residues :
            if r.id.chainId == chainId :
                doRess.append ( r )

        print " - calc for %d res..." % len (doRess)
        waterRad = 1.4
        waterRad2 = 1.4*1.4

        rt_sa = {}

        for ri, r in enumerate ( doRess ) :

            if 0 or r.type == "ASP" :

                showPts = r.type == "ASP"

                if 1 or not hasattr ( r, 'SAArea' ) :

                    numPtsOnSAS, tryPts = 0.0, 300.0
                    for at in r.scAtoms :
                        VWR = vdwRadii[at.element.name] if at.element.name in vdwRadii else 1.55
                        outPts = SpherePts ( at.coord(), VWR + waterRad, int(tryPts) )
                        #AddSpherePts ( outPts, (.9,.9,.2,1.0), 0.1, "ASP pts" )
                        for pt in outPts :
                            vPt = [pt[0], pt[1], pt[2]]; apt = numpy.array ( vPt )
                            opointsNear = surfPtsTree.searchTree ( vPt, waterRad )
                            onSurf = True
                            for npt in opointsNear :
                                v = apt - npt; r2 = numpy.sum ( v * v )
                                if r2 < waterRad2 :
                                    onSurf = False; break
                            if onSurf :
                                numPtsOnSAS += 1.0
                                if showPts :
                                    v = chimera.Point(pt[0], pt[1], pt[2]) - at.coord(); v.normalize()
                                    pt = at.coord() + v * vdwRadii[at.element.name]
                                    AddSpherePts ( [pt.data()], (.9,.2,.9,1.0), 0.1, "SA pts" )
                                    #AddSpherePts ( [vPt], (.2,.9,.9,0.8), 0.11, "SA- pts" )

                    r.SAArea = 4 * numpy.pi * numpy.power ( vdwRadii[at.element.name], 2.0  ) * numPtsOnSAS / tryPts

                if hasattr (r, 'qSC') and r.qSC != None :
                    fp.write ( "%s\t%d\t%f\t%f\n" % (r.type, r.id.position, r.qSC, r.SAArea) )
                elif hasattr (r, 'Q') and r.Q != None :
                    fp.write ( "%s\t%d\t%f\t%f\n" % (r.type, r.id.position, r.Q, r.SAArea) )

                if r.type in rt_sa :
                    rt_sa[r.type] += r.SAArea
                else :
                    rt_sa[r.type] = r.SAArea
            if ri % 10 == 0 :
                umsg ( "SA - res %d/%d" % (ri, len(doRess)) )

        fp.close()
        umsg ( "Solvent Accessibility vs. Q ... saved to file %s ... done" % nname )

        #nname = molPath + "__SAa-Q__" + mapName + ".txt"
        #fp = open ( nname, "w" )

        print "SA area by rtype:"
        totalSA = 0.0
        for rt, saa in rt_sa.iteritems () :
            print "%s\t%f" % (rt, saa)
            totalSA += saa

        print " - total SA:%f" % totalSA

        print "SA area / total area by rtype:"
        for rt, saa in rt_sa.iteritems () :
            print "%s\t%f" % (rt, saa/totalSA)





    def CalcAllRadZ (self) :

        ress = []
        try :
            ress = self.seqRes
        except :
            pass

        if len ( ress ) == 0 :
            umsg ( "No molecule/chain selected?" )
            return

        if self.cur_dmap == None :
            status ( "Select or open a map..." )
            return

        if self.cur_mol == None :
            status ( "Select or open a model..." )
            return

        cid = self.chain.get()

        ats = [at for at in self.cur_mol.atoms if not at.element.name == "H"]
        points = _multiscale.get_atom_coordinates ( ats, transformed = False )
        print " - search tree: %d/%d ats" % ( len(ats), len(self.cur_mol.atoms) )
        allAtTree = AdaptiveTree ( points.tolist(), ats, 1.0)
        #allAtTree = None


        CalcRadZ ( self.cur_mol, cid, self.cur_dmap, allAtTree, useOld=False, log=True )

        self.scores, self.scores2 = [], []
        scBB, scSC = [], []

        for r in self.cur_mol.residues :
            if cid == None or r.id.chainId == cid :
                self.scores2.append ( r.bbZ )
                self.scores.append ( r.scZ )
                if r.bbZ != None :
                    scBB.append ( r.bbZ )
                if r.scZ != None :
                    scSC.append ( r.scZ )


        #bbRes = numpy.power ( numpy.e, (self.avgScore2 - 8.0334) / -4.128 ) # y = -4.128ln(x) + 8.0334
        #scRes = numpy.power ( numpy.e, (self.avgScore - 4.8261) / -3.097 ) # y = -3.097ln(x) + 4.8261
        #scRes = (self.avgScore2 - 3.507) / -0.721
        #bbRes = (self.avgScore - 6.1234) / -0.9191

        scMin, scMax, scAvg = min(scSC), max(scSC), numpy.average(scSC)
        bbMin, bbMax, bbAvg = min(scBB), max(scBB), numpy.average(scBB)

        print "Average RadZ sc : %.2f - %.2f, avg %.2f" % (scMin, scMax, scAvg)
        print "Average RadZ bb : %.2f - %.2f, avg %.2f" % (bbMin, bbMax, bbAvg)

        umsg ( "Average Side Chain: %.2f, Backbone: %.2f" % (scAvg, bbAvg) )


        self.minSCscore, self.maxSCscore = 0.0,4
        self.minBBscore, self.maxBBscore = 0.0,4

        self.UpdateSeq ()




    def CalcAllRotaZ (self) :

        ress = []
        try :
            ress = self.seqRes
        except :
            pass

        if len ( ress ) == 0 :
            umsg ( "No molecule/chain selected?" )
            return

        if self.cur_dmap == None :
            status ( "Select or open a map..." )
            return

        if self.cur_mol == None :
            status ( "Select or open a model..." )
            return

        cid = self.chain.get()

        self.scores, self.scores2 = [], []
        scBB, scSC = [], []

        print "..."

        #CalcRotaZ ( self.cur_dmap, self.cur_mol, self.cur_mol.residues )


        for r in self.cur_mol.residues :



            for at in r.atoms :
                if not hasattr ( at, 'isBB' ) :
                    print " - noBB - atom %s, res %d.%s, chain %s" % (at.name, at.residue.id.position, at.residue.type, at.residue.id.chainId)

            if cid == None or r.id.chainId == cid :
                self.scores2.append ( r.bbZ )
                self.scores.append ( r.scZ )
                if r.bbS != None :
                    scBB.append ( r.bbZ )
                if r.scS != None :
                    scSC.append ( r.scZ )


        #bbRes = numpy.power ( numpy.e, (self.avgScore2 - 8.0334) / -4.128 ) # y = -4.128ln(x) + 8.0334
        #scRes = numpy.power ( numpy.e, (self.avgScore - 4.8261) / -3.097 ) # y = -3.097ln(x) + 4.8261
        #scRes = (self.avgScore2 - 3.507) / -0.721
        #bbRes = (self.avgScore - 6.1234) / -0.9191

        scMin, scMax, scAvg = min(scSC), max(scSC), numpy.average(scSC)
        bbMin, bbMax, bbAvg = min(scBB), max(scBB), numpy.average(scBB)


        umsg ( "Average Sigma sc : %.2f - %.2f, avg %.2f | %.2f - %.2f, avg %.2f" % (scMin, scMax, scAvg, 1.0/scMin, 1.0/scMax, 1.0/scAvg) )
        umsg ( "Average Sigma bb : %.2f - %.2f, avg %.2f | %.2f - %.2f, avg %.2f" % (bbMin, bbMax, bbAvg, 1.0/bbMin, 1.0/bbMax, 1.0/bbAvg) )


        self.minSCscore, self.maxSCscore = 0.0,2.0
        self.minBBscore, self.maxBBscore = 0.0,2.0

        self.UpdateSeq ()




    def RtypeOut ( self, avgScore, rtype, rByType, fout ) :
        pass




    def UpdateSeqFont ( self ) :
        # http://stackoverflow.com/questions/4296249/how-do-i-convert-a-hex-triplet-to-an-rgb-tuple-and-back

        if not hasattr ( self, 'seq' ) :
            print " - update seq font - no seq"
            return

        #print "seq len %d, text w %d" % ( len(self.seq), self.tw )

        # boxes for BBs
        x_at = self.seqX
        y_at = self.seqY + self.seqH/2

        y0 = self.seqY+5
        y1 = self.seqY+self.seqH-5

        for si in range ( len(self.seq) ) :
            res = self.seq[si]
            pred = self.pred[si]
            conf = float ( self.conf[si] ) / 10.0

            if pred == 'E' :
                x0 = self.seqX + si * self.tw
                x1 = x0 + self.tw
                #self.Canvas.coords ( self.seqMouseR, x0, y0, x1, y1 )
                #self.Canvas.itemconfigure ( self.seqMouseR, state=Tkinter.NORMAL )

                if self.seqSheetR[si] == None :
                    c = self.sheetBaseClr + self.sheetClrD * conf
                    clr = "#" + struct.pack('BBB',c[0],c[1],c[2]).encode('hex')
                    self.seqSheetR[si] = self.Canvas.create_rectangle(x0, y0, x1, y1, outline=clr, fill=clr)
                else :
                    self.Canvas.coords ( self.seqSheetR[si], x0, y0, x1, y1 )

            elif pred == 'H' :
                x0 = self.seqX + si * self.tw
                x1 = x0 + self.tw

                if self.seqHelixR[si] == None :
                    c = self.helixBaseClr + self.helixClrD * conf
                    clr = "#" + struct.pack('BBB',c[0],c[1],c[2]).encode('hex')
                    self.seqHelixR[si] = self.Canvas.create_rectangle(x0, y0, x1, y1, outline=clr, fill=clr)
                else :
                    self.Canvas.coords ( self.seqHelixR[si], x0, y0, x1, y1 )



        # box showing selected Residue
        if hasattr ( self, 'seqMouseR' ) :
            self.Canvas.coords ( self.seqMouseR, 0, 0, 0, 0 )
        else :
            self.seqMouseR = self.Canvas.create_rectangle(0, 0, 0, 0, outline="#aab", fill="#bbc", state=Tkinter.HIDDEN)



        x_at = self.seqX
        y_at = self.seqY + self.seqH/2

        if hasattr ( self, 'seqText' ) and self.seqText != None :
            #self.Canvas.coords ( self.seqText, x_at, y_at )
            #self.Canvas.itemconfigure ( self.seqText, font=self.font )
            #print " - has seq?"
            self.Canvas.delete ( self.seqText )
            #self.seqText = None
            #del self.seqText

        self.seqText = self.Canvas.create_text( x_at, y_at, text=self.seq, font=self.font, anchor='w')
        print " - created seq text - font"


        #self.UpdateSeqSel ()


    def GetMaxScores ( self ) :

        self.minScore1, self.maxScore1 = 0.0,1.0
        self.minScore2, self.maxScore2 = 0.0,1.0

        try :
            RES = float(self.mapRes.get())
        except :
            umsg ( "Please enter a numeric value for Resolution in Options" )
            return

        try :
            sigma = float(self.sigma.get())
        except :
            umsg ( "Please enter a numeric value for sigma in Options" )
            return

        #avgQrna = -0.1574 * RES + 1.0673 # rna
        #avgQprot = -0.1794 * RES + 1.1244 # protein
        #avgQIon =  -0.1103 * RES + 1.0795 # ion
        #avgQWater =  -0.0895 * RES + 1.0001 # water

        expQScore, lowQ, highQ, eqn = qscores.ExpectedQScore ( RES, sigma )

        print " - res %.2f - expected Q-score: %.2f" % (RES, expQScore)

        if 1 :
            self.minScore1, self.maxScore1 = 0.0,expQScore
            self.minScore2, self.maxScore2 = 0.0,expQScore
        else :
            self.minScore1, self.maxScore1 = 0.0,1.0
            self.minScore2, self.maxScore2 = 0.0,1.0

        #self.minScore1, self.maxScore1 = lowQ,expQScore
        #self.minScore2, self.maxScore2 = lowQ,expQScore


    def UpdateSeq ( self ) :

        if not hasattr ( self, 'seq' ) :
            print " - update seq - no seq"
            return

        if not hasattr ( self, 'maxScore1' ) :
            self.GetMaxScores ()

        x_at = self.seqX
        y_at = self.seqY + self.seqH/2

        if hasattr ( self, 'seqText' ) and self.seqText != None :
            self.Canvas.coords ( self.seqText, x_at, y_at )
        #else :
        #    self.seqText = self.Canvas.create_text( x_at, y_at, text=self.seq, font=self.font, anchor='w')
        #    print " - created seq text"

        if 1 :
            y0 = self.seqY+5
            y1 = self.seqY+self.seqH-5

            cH = numpy.array ( [.33*255.0,.56*255.0,.88*255.0] )
            cL = numpy.array ( [.99*255.0,.99*255.0,.3*255.0] )

            # qRedGreen
            if 1 :
                cH = numpy.array( [50,250,50] )
                cL = numpy.array( [250,50,50] )

            for si in range ( len(self.seq) ) :
                #if i >= len ( self.seqt ) :
                #    t = self.Canvas.create_text( x_at, y_at, text=self.seq[i], font=self.font)
                #    self.seqt.append ( t )
                #else :
                #    t = self.seqt [ i ]
                #    self.Canvas.coords ( t, x_at, y_at )
                # x_at += self.tw

                res = self.seqRes[si]

                pred = self.pred[si]
                if pred == 'E' :
                    if self.seqSheetR[si] != None :
                        x0 = self.seqX + si * self.tw
                        x1 = x0 + self.tw
                        self.Canvas.coords ( self.seqSheetR[si], x0, y0, x1, y1 )

                elif pred == 'H' :
                    if self.seqHelixR[si] != None :
                        x0 = self.seqX + si * self.tw
                        x1 = x0 + self.tw
                        self.Canvas.coords ( self.seqHelixR[si], x0, y0, x1, y1 )

                if res == None :
                    continue


                if not hasattr(res, 'score1') or res.score1 == None :
                    if self.seqScoreR[si] != None :
                        self.Canvas.delete ( self.seqScoreR[si] )
                    self.seqScoreR[si] = None
                else :
                    xx0 = self.seqX + si * self.tw + 2
                    xx1 = xx0 + self.tw - 2
                    h = (res.score1 - self.minScore1) / (self.maxScore1 - self.minScore1)
                    if h > 1 : h = 1
                    if h < 0 : h = 0
                    Y, H = self.modY, (self.modH/2 - 2)
                    yy0, yy1 = numpy.ceil(Y+H - H*h), numpy.floor(Y+H)
                    #c = self.helixBaseClr + self.helixClrD * conf
                    c = h * cH + (1-h) * cL
                    clr = "#" + struct.pack('BBB',c[0],c[1],c[2]).encode('hex')
                    if self.seqScoreR[si] != None :
                        self.Canvas.coords ( self.seqScoreR[si], xx0, yy0, xx1, yy1 )
                        self.Canvas.itemconfigure ( self.seqScoreR[si], outline=clr, fill=clr )
                    else :
                        self.seqScoreR[si] = self.Canvas.create_rectangle(xx0, yy0, xx1, yy1, outline=clr, fill=clr)

                if not hasattr(res, 'score2') or res.score2 == None :
                    if self.seqScoreR2[si] != None :
                        self.Canvas.delete ( self.seqScoreR2[si] )
                    self.seqScoreR2[si] = None
                else :
                    xx0 = self.seqX + si * self.tw + 2
                    xx1 = xx0 + self.tw - 2
                    h = (res.score2 - self.minScore2) / (self.maxScore2 - self.minScore2)
                    if h > 1 : h = 1
                    if h < 0 : h = 0
                    Y, H = self.modY, self.modH/2
                    #yy0, yy1 = Y+H, Y+H+H*h #upside down chart
                    yy0, yy1 = numpy.ceil(Y+H+H-H*h), numpy.floor(Y+H+H)
                    #c = self.helixBaseClr + self.helixClrD * conf
                    c = h * cH + (1-h) * cL
                    clr = "#" + struct.pack('BBB',c[0],c[1],c[2]).encode('hex')
                    if self.seqScoreR2[si] != None :
                        self.Canvas.coords ( self.seqScoreR2[si], xx0, yy0, xx1, yy1 )
                        self.Canvas.itemconfigure ( self.seqScoreR2[si], outline=clr, fill=clr )
                    else :
                        self.seqScoreR2[si] = self.Canvas.create_rectangle(xx0, yy0, xx1, yy1, outline=clr, fill=clr)


            self.UpdateSeqSel ()




    def SeqRec ( self, sel ) :
        y0 = self.seqY+5
        y1 = self.seqY+self.seqH-5

        x0 = self.seqX + sel[0] * self.tw
        x1 = self.seqX + (sel[1]+1) * self.tw

        return x0, y0, x1, y1


    def UpdateSeqSel ( self ) :

        if not hasattr ( self, 'seqSel' ) :
            return

        if self.seqSel == None :
            if hasattr(self, 'seqSelRect') :
                self.Canvas.delete ( self.seqSelRect )
                self.seqSelRect = None
            return

        x0, y0, x1, y1 = self.SeqRec ( self.seqSel )

        if hasattr(self, 'seqSelRect') and self.seqSelRect != None :
            self.Canvas.coords ( self.seqSelRect, x0, y0, x1, y1  )
        else :
            #c = self.helixBaseClr + self.helixClrD * conf
            #clr = "#" + struct.pack('BBB',c[0],c[1],c[2]).encode('hex')
            self.seqSelRect = self.Canvas.create_rectangle(x0, y0, x1, y1, outline=self.selColor,  width=3)





    def B1_Down (self, event):
        self.drag = ''

        #print "b1 _", event.x, event.y
        if self.isInSeq ( event.x, event.y ) :
            self.drag = 'seq'
        self.last_x = event.x
        self.last_y = event.y


    def B1_Down_Ctrl ( self, event ) :
        #print "b1 _ <ctrl>", event.x, event.y
        self.drag = ''

        if self.isInSeq ( event.x, event.y ) :
            self.drag = 'seqSel'

            if hasattr ( self, 'seqSel' ) and self.seqSel != None :
                self.prevSeqSel = self.seqSel
            else :
                self.prevSeqSel = None

            #print "sel seq..."
            seqI = ( event.x - self.seqX ) / self.tw

            resStartI = seqI
            try :
                resStartI = self.seqRes[seqI].id.position
            except :
                pass

            status ( "Start sequence sel at %d" % resStartI )

            self.seqSel = [seqI, seqI]
            self.UpdateSeqSel ()

        self.last_x = event.x
        self.last_y = event.y


    def B1_Down_Shift ( self, event ) :
        print "B1 down - shift"

        self.drag = ''

        if self.isInSeq ( event.x, event.y ) :
            if hasattr ( self, 'seqSel' ) and self.seqSel != None :
                seqI = ( event.x - self.seqX ) / self.tw
                if seqI >= self.seqSel[0] and seqI <= self.seqSel[1] :
                    self.drag = "con"
                    if not hasattr ( self, 'conLine' ) or self.conLine == None :
                        self.conLine = self.Canvas.create_line( event.x, event.y, event.x, event.y, fill="red", dash=(1, 1), width=2)
                    status ( "In selected sequence at %d" % seqI )


    def B1_Down_Alt ( self, event ) :
        print "B1 down - alt"

        self.drag = ''

        if self.isInMod ( event.x, event.y ) :
            self.dragMod = self.SelectedMod ( event.x, event.y )
            if self.dragMod != None :
                if self.dragMod.type == "Helix" :
                    self.drag = 'modRot'
                    self.dragStartX = event.x



    def B1_Up_Ctrl ( self, event ) :
        #print "b1 up - ctrl - ", event.x, event.y
        self.B1_Up ( event )


    def B1_Up_Shift ( self, event ) :
        #print "b1 up - shift - "
        self.B1_Up ( event )

    def B1_Up_Alt ( self, event ) :
        #print "b1 up - alt - "
        self.B1_Up ( event )


    def B1_Up (self, event) :
        #print "b1 up - ", event.x, event.y

        if self.drag == 'seqSel' and hasattr ( self, 'seqSel' ) :
            status ( "Selected: %d-%d" % (self.seqSel[0], self.seqSel[1]) )

            if hasattr ( self, 'prevSeqSel' ) and self.prevSeqSel != None :
                if self.seqSel[0] == self.seqSel[1] :
                    self.seqSel = None
                    self.prevSeqSel = None
                    self.UpdateSeqSel ()
                    status ( "Cleared sequence selection" )
                    chimera.selection.clearCurrent ()

            if self.seqSel != None :
                m, cid = self.cur_mol, self.chain.get()
                if m != None :

                    #startI = self.seqRes [ max(self.seqSel[0],0) ].id.position
                    #endI = self.seqRes [ min(self.seqSel[1],len(self.seqRes)-1) ].id.position

                    startI = self.seqRi [ max(self.seqSel[0],0) ]
                    endI = self.seqRi [ min(self.seqSel[1],len(self.seqRes)-1) ]

                    selStr = "#%d:%d-%d.%s" % (m.id,startI,endI,cid)

                    self.lastSelStr = selStr # "%d-%d.%s" % (startI,endI,cid)
                    sel = chimera.selection.OSLSelection ( selStr )
                    chimera.selection.setCurrent ( sel )

                    if hasattr ( self, 'prevSel' ) and self.preserveSel.get () :
                        for s in self.prevSel :
                            print " -s- adding to sel:", s
                            chimera.selection.mergeCurrent ( chimera.selection.EXTEND, chimera.selection.OSLSelection (s) )
                    else :
                        self.prevSel = []

                    if self.preserveSel.get () :
                        #self.prevSel.append ( "%d-%d.%s" % (startI,endI,cid) )
                        self.prevSel.append ( selStr )
                        print " - added to selection list: ", selStr

                    umsg ( "Selected: " + selStr )
                    #chimera.selection.addCurrent ( sel )

                    if self.selExtract.get () :
                        self.ShowSel ()

                else :
                    status ( "no model visible" )

            #else :
            #    print "cleared past sel"
            #    self.prevSel = []


        elif self.drag == 'modSel' :
            status ( 'Selected %d mods' % len(self.selMods) )

        elif self.drag == 'con' :
            selMod = None
            if hasattr ( self, 'selModPiece' ) and self.selModPiece != None :
                selMod = self.selModPiece
                self.selModPiece = None
            else :
                return

            if hasattr ( self, 'conLine' ) and self.conLine != None :
                self.Canvas.delete ( self.conLine )
                self.conLine = None

            status ( "connected to %s" % selMod.type )

            selMod.seq = self.seqSel
            selMod.numRes = (self.seqSel[1] - self.seqSel[0] + 1)
            selMod.MakeMod ()

            self.UpdateMod ()

        self.drag = ''
        #print "mod: ", self.modX, " seq:", self.seqX


    def KeepBack ( self ) :

        if not hasattr ( self, 'prevSel' ) or self.prevSel == None :
            umsg ( "Nothing selected previously... select something by Ctrl+Click+Drag on the sequence; this undoes the last selection, use Keep" )
            return

        if hasattr ( self, 'prevSel' ) and len(self.prevSel) > 0 :
            self.prevSel.pop()

            chimera.selection.clearCurrent()

            for s in self.prevSel :
                print " -s- adding to sel:", s
                chimera.selection.mergeCurrent ( chimera.selection.EXTEND, chimera.selection.OSLSelection (s) )

            if self.selExtract.get () :
                self.ShowSel ()



    def SelReLoad ( self ) :

        atoms = chimera.selection.currentAtoms()

        if len(atoms) > 0 :
            if self.cur_dmap == None :
                umsg ( "Select a Map" )
                return

            rad = 2.0
            dmap = self.cur_dmap
            m = atoms[0].molecule

            from _multiscale import get_atom_coordinates
            points = get_atom_coordinates ( atoms, transformed = True )

            if not self.preserveSel.get() :
                delMods = []
                for m in chimera.openModels.list() :
                    if "_mapq_zz_" in m.name :
                        delMods.append ( m )
                chimera.openModels.close ( delMods )

            nname = os.path.splitext(dmap.name)[0] + "_mapq_zz_%.0f_" % rad + ".mrc"
            cmap = self.PtsToMap ( points, dmap, rad, nname, showMesh=False, alpha=0.4 )
            #      self.PtsToMap ( points, dmap, R, dmap.name + label + "_mesh", True )

            #umsg ( "Made zone map: " + nname )
            dmap.display = False

            M = dmap.data.full_matrix()
            sdev = numpy.std(M)
            avg = numpy.average(M)

            #cmap.surface_levels = [avg + 3.0 * sdev]
            cmap.surface_levels = self.cur_dmap.surface_levels
            #chimera.runCommand ( "vol #%d style surface region all step 1 color grey transparency 0.8" % cmap.id )
            chimera.runCommand ( "vol #%d style surface region all step 1 color grey" % cmap.id )
            for sp in cmap.surfacePieces :
                v, t = sp.geometry
                if not (len(v) == 8 and len(t) == 12) :
                    sp.color = sp.color = (.7, .7, .7, .3)


        elif 0 :
            if not hasattr ( self, 'prevSel' ) or self.prevSel == None :
                umsg ( "Nothing selected previously... select something by Ctrl+Click+Drag on the sequence; this refreshes the selection" )
                return

            for s in self.prevSel :
                print " -s- adding to sel:", s
                chimera.selection.mergeCurrent ( chimera.selection.EXTEND, chimera.selection.OSLSelection (s) )

            if self.selExtract.get () :
                self.ShowSel ()

    def SelLoad ( self ) :

        self.prevSel = []

        if 1 :
            self.prevSel.append ( "#%d:735-735.A" % self.cur_mol.id )
            self.prevSel.append ( "#%d:796-796.A" % self.cur_mol.id )
            self.prevSel.append ( "#%d:799-799.A" % self.cur_mol.id )
            self.prevSel.append ( "#%d:137-137.L" % self.cur_mol.id )
            self.prevSel.append ( "#%d:108-108.C" % self.cur_mol.id )
            self.prevSel.append ( "#%d:789-789.A" % self.cur_mol.id )

        elif 0 :
            self.prevSel.append ( "#%d:41-41.A" % self.cur_mol.id )
            self.prevSel.append ( "#%d:7-7.A" % self.cur_mol.id )
            #self.prevSel.append ( "#%d:63-63.A" % self.cur_mol.id )
            self.prevSel.append ( "#%d:2-2.A" % self.cur_mol.id )
            self.prevSel.append ( "#%d:7-7.A" % self.cur_mol.id )



        for s in self.prevSel :
            print " -s- adding to sel:", s
            chimera.selection.mergeCurrent ( chimera.selection.EXTEND, chimera.selection.OSLSelection (s) )

        if self.selExtract.get () :
            self.ShowSel ()





    def SelText ( self ) :

        self.prevSel = []

        print self.selText.get()

        ls = self.selText.get().split(";")

        for l in ls :
            #self.prevSel.append ( "#%d:%s" % (self.cur_mol.id,l) )
            self.prevSel.append ( "%s" % l )

        for s in self.prevSel :
            print " -s- adding to sel:", s
            chimera.selection.mergeCurrent ( chimera.selection.EXTEND, chimera.selection.OSLSelection (s) )

        if self.selExtract.get () :
            self.ShowSel ()


        if 0 :
            fp = os.path.split ( self.cur_dmap.data.path )[0] + "/_sel.txt"
            found = False
            ls = []
            try :
                fo = open ( fp, "r" )
                for l in fo :
                    if self.selText.get() in l :
                        found = True
                fo.close()
            except :
                pass

            if found :
                print " - found sel text"
            else :
                fo = open ( fp, "a" )
                fo.write ( "%s\n" % self.selText.get() )
                fo.close()





    def preserveSelCb (self) :
        print "Preserve set to ", self.preserveSel.get()
        if self.preserveSel.get() :
            print " - setting current selection to preserve..."
            if hasattr ( self, 'lastSelStr' ) :
                self.prevSel = [ self.lastSelStr ]
        else :
            print " - clearing current"
            self.prevSel = []


    def preserveVolCb (self) :
        print "Preserve vol set to ", self.preserveVol.get()


    #def keepExMapCb (self) :
    #    print "Kep ex map set to ", self.keepExMap.get()


    def ClearSel ( self ) :
        self.prevSel = []
        self.seqSel = None
        self.prevSeqSel = None
        self.UpdateSeqSel ()
        status ( "Cleared sequence selection" )
        chimera.selection.clearCurrent ()




    def ExCustA ( self ) :

        if self.cur_dmap == None :
            umsg ("Select a map first")
            return

        if self.cur_mol == None :
            umsg ("Select a molecule first")
            return

        #selStr = "#%d:80-87.I,171-184.I,227-239.I,80-87.A,171-184.A,227-239.A,80-87.B,171-184.B,227-239.B,80-87.J,171-184.J,227-239.J,80-87.H,171-184.H,227-239.H" % self.cur_mol.id
        selStr = "#%d:80-87.I,171-184.I,227-239.I,80-87.A,171-184.A,227-239.A,80-87.J,171-184.J,227-239.J" % self.cur_mol.id

        umsg ( "Selected: " + selStr )
        sel = chimera.selection.OSLSelection ( selStr )
        chimera.selection.setCurrent ( sel )
        self.ShowSel()


    def ExCustB ( self ) :

        if self.cur_dmap == None :
            umsg ("Select a map first")
            return

        if self.cur_mol == None :
            umsg ("Select a molecule first")
            return

        selStr = "#%d:428-435.F,365-376.F,396-402.F,428-435.I,365-376.I,396-402.I" % self.cur_mol.id
        #selStr = "#%d:428-435.A,365-376.A,396-402.A,428-435.H,365-376.H,396-402.H" % self.cur_mol.id


        umsg ( "Selected: " + selStr )
        sel = chimera.selection.OSLSelection ( selStr )
        chimera.selection.setCurrent ( sel )
        self.ShowSel()

    def ExCustC ( self ) :

        if self.cur_dmap == None :
            umsg ("Select a map first")
            return

        if self.cur_mol == None :
            umsg ("Select a molecule first")
            return

        #selStr = "#%d:10:548-558.I,520-530.I,548-558.F,520-530.F" % self.cur_mol.id
        selStr = "#%d:428-435.F,365-376.F,396-402.F,428-435.I,365-376.I,396-402.I,548-558.I,520-530.I,548-558.F,520-530.F" % self.cur_mol.id


        umsg ( "Selected: " + selStr )
        sel = chimera.selection.OSLSelection ( selStr )
        chimera.selection.setCurrent ( sel )
        self.ShowSel()


    def AdSel ( self ) :

        atoms = chimera.selection.currentAtoms ()

        R = float ( self.maskRad.get() )

        if len(atoms) > 0 :

            from _multiscale import get_atom_coordinates
            points = get_atom_coordinates ( atoms, transformed = True )
            COM, U, S, V = prAxes ( points )

            #atomRad = 2.0 # float ( self.maskWithSelDist.get() )
            print " - %d selected atoms, mask at %.2f" % ( len(atoms), R )
            dmap = self.cur_dmap

            label = " %d_ats_%s.%d.%s mask" % (len(atoms), atoms[0].name, atoms[0].residue.id.position, atoms[0].residue.id.chainId )

            if len ( atoms ) > 0 and dmap != None :
                #points = get_atom_coordinates ( atoms, transformed = False )
                if self.showMesh.get () :
                    self.PtsToMap ( points, dmap, R, dmap.name + label, False, alpha=0.2 if self.showMesh.get() else 0.4 )
                    self.PtsToMap ( points, dmap, R, dmap.name + label + "_mesh", True )
                else :
                    self.PtsToMap ( points, dmap, R, dmap.name + label, False, alpha=0.4 if self.showMesh.get() else 0.4 )


    def HideSel ( self ) :

        for r in chimera.selection.currentResidues() :
            r.ribbonDisplay = False
            for at in r.atoms :
                at.display = False


    def ShowSel ( self ) :

        #showRibbon = self.showRibbon.get()
        showRibbon = not self.showingAtoms # self.showRibbon.get()
        showLigands = self.showLigands.get()
        showSC = True # self.showAtoms.get()
        showH = self.showH.get()
        showW = self.showW.get()

        atoms = []
        scores = []
        selResM = {}

        if len ( chimera.selection.currentResidues () ) == 0 :
            umsg ( "Nothing selected..." )
            return

        for r in chimera.selection.currentResidues () :
            rid = "%d.%s" % (r.id.position, r.id.chainId)
            selResM [rid] = 1

        if self.cur_mol == None :
            return

        if 1 or not hasattr ( self.cur_mol, 'bbats' ) :
            qscores.SetBBAts(self.cur_mol)
            self.cur_mol.bbats = True


        atMap = {}
        for r in self.cur_mol.residues :
            rid = "%d.%s" % (r.id.position, r.id.chainId)
            if rid in selResM :

                if hasattr (r, 'scZ') and r.scZ != None :
                    scores.append(r.scZ)

                r.ribbonDisplay = showRibbon

                for at in r.atoms :
                    atMap[at] = 1
                    if at.element.name == "H" :
                        at.display = showH
                    elif at.isSC :
                        if showSC :
                            at.display = True
                            atoms.append ( at )
                        else :
                            at.display = False
                    else :
                        at.display = True
                        atoms.append ( at )
                    if at.element.name in atomColors :
                        if at.element.name == "H" :
                            continue
                        if at.isBB :
                            at.color = atomColors[at.element.name.upper()]
                            #if at.element.name == "C" :
                            #    at.color = atomColors['Cbb']
                        else :
                            at.color = atomColors[at.element.name.upper()]

            else :
                r.ribbonDisplay = False
                for at in r.atoms :
                    at.display = False


        atTree = None
        if showLigands :
            points = _multiscale.get_atom_coordinates ( atoms, transformed = False )
            print " - search tree: %d/%d ats" % ( len(atoms), len(r.molecule.atoms) )
            atTree = AdaptiveTree ( points.tolist(), atoms, 2.0)
            from chimera.resCode import protein3to1

            ligAts = []
            for r in self.cur_mol.residues :
                rid = "%d.%s" % (r.id.position, r.id.chainId)

                #if r.type == "MG" or r.type == "HOH" :
                if not r.isProt and not r.isNA :
                    if len ( self.AtsWithin (r.atoms, 4.0, atTree) ) > 0 :
                        for at in r.atoms :
                            at.display = True
                            atMap[at] = 1

                            at.radius = 1.46
                            at.drawMode = at.EndCap # if at.element.name.lower() == "o" else at.Ball

                            if at.element.name.upper() in atomColors :
                                at.color = atomColors[at.element.name.upper()]
                            atoms.append ( at )
                            ligAts.append ( at )
                    else :
                        for at in r.atoms :
                            at.display = False

            #chimera.selection.clearCurrent ()
            print " - added %d ligand atoms to sel" % len(ligAts)
            chimera.selection.addCurrent ( ligAts )


        for bond in self.cur_mol.bonds :
            a1, a2 = bond.atoms
            if a1 in atMap and a2 in atMap :
                if showW :
                    bond.drawMode = bond.Wire
                else :
                    bond.drawMode = bond.Stick
                #bond.display = bond.Smart


        if len(atoms) > 0 :

            from _multiscale import get_atom_coordinates
            points = get_atom_coordinates ( atoms, transformed = True )
            COM, U, S, V = prAxes ( points )

            moveCam = 1
            if moveCam :
                p0 = numpy.array ( chimera.viewer.camera.center )
                p1 = numpy.array ( [ COM[0], COM[1], COM[2] ] )
                for i in range (10) :
                    f = float(i) / 9.0
                    f1, f2 = 2.0*f*f*f-3.0*f*f+1.0, 3*f*f-2*f*f*f
                    P = p0 * f1 + p1 * f2
                    chimera.viewer.camera.center = (P[0],P[1],P[2])
                    print ".",
                print ""


            atomRad = 2.5 # float ( self.maskWithSelDist.get() )
            print " - %d selected atoms, mask at %.2f" % ( len(atoms), atomRad )
            dmap = self.cur_dmap

            mlist = OML(modelTypes = [VolumeViewer.volume.Volume])

            at = 1
            for m in mlist :
                if "sel_masked" in m.name :
                    mname = m.name.split()[0]
                    if not hasattr (self, 'cLevels') :
                        self.cLevels = {}
                    if not "_mesh" in m.name :
                        self.cLevels[mname] = m.surface_levels[0]
                    if not self.preserveVol.get() :
                        chimera.openModels.close ( [m] )
                    else :
                        m.name = m.name.replace ( "sel_masked", "prev_masked" )

            if len ( atoms ) > 0 and dmap != None :
                #points = get_atom_coordinates ( atoms, transformed = False )
                self.PtsToMap ( points, dmap, atomRad, dmap.name + " sel_masked", False )
                if self.showMesh.get () :
                    self.PtsToMap ( points, dmap, atomRad, dmap.name + " sel_masked_mesh", True )

            if len(scores) > 0 :
                umsg ( "%d residue scores, avg score %.1f" % ( len(scores), numpy.average(scores) ) )

        else :
            umsg ( "no atoms selected, try reselecting the model and chain..." )


    def AtsWithin (self, ats, R, atTree) :

        nearAts = []
        R2 = R * R
        for at in ats :
            pt = at.coord()
            vPt = numpy.array ( pt.data() )
            opointsNear = atTree.searchTree ( [pt[0], pt[1], pt[2]], R )
            if len(opointsNear) > 0 :
                for p in opointsNear :
                    try :
                        v = vPt - p.coord().data()
                    except :
                        continue
                    sqSum = numpy.sum ( v * v )
                    if sqSum < R2 :
                        nearAts.append (p)

        return nearAts


    def AtsWithinPt (self, pt, R, atTree) :

        nearAts = []
        R2 = R * R

        vPt = numpy.array ( [pt[0], pt[1], pt[2]] )
        opointsNear = atTree.searchTree ( [pt[0], pt[1], pt[2]], R )
        if len(opointsNear) > 0 :
            for p in opointsNear :
                try :
                    v = vPt - p.coord().data()
                except :
                    continue
                sqSum = numpy.sum ( v * v )
                if sqSum < R2 :
                    nearAts.append ( [numpy.sqrt(sqSum), p] )

        return nearAts




    def PtsToMap0 ( self, points, dmap, atomRad, nname, neg = 1.0 ) :
        import _contour
        _contour.affine_transform_vertices ( points, Matrix.xform_matrix( dmap.openState.xform.inverse() ) )
        mdata = VolumeData.zone_masked_grid_data ( dmap.data, points, atomRad )

        gdata = VolumeData.Array_Grid_Data ( mdata.full_matrix()*neg, dmap.data.origin, dmap.data.step, dmap.data.cell_angles, name = nname )
        nv = VolumeViewer.volume.volume_from_grid_data ( gdata )
        nv.name = nname
        dmap.display = False
        nv.region = ( nv.region[0], nv.region[1], [1,1,1] )
        nv.surface_levels[0] = dmap.surface_levels[0]
        ro = VolumeViewer.volume.Rendering_Options()
        #ro.smoothing_factor = .5
        #ro.smoothing_iterations = 20
        #ro.surface_smoothing = True
        nv.update_surface ( False, ro )
        for sp in nv.surfacePieces :
            v, t = sp.geometry
            if len(v) == 8 and len(t) == 12 :
                sp.display = False
            else :
                sp.color = (0.7, 0.7, 0.7, 0.4)


    def PtsToMap ( self, points, dmap, atomRad, nname, showMesh = False, alpha=0.4 ) :

        #_contour.affine_transform_vertices ( points, Matrix.xform_matrix( dmap.openState.xform.inverse() ) )
        #mdata = VolumeData.zone_masked_grid_data ( dmap.data, points, atomRad )

        import _contour
        points1 = numpy.copy ( points )
        _contour.affine_transform_vertices ( points1, Matrix.xform_matrix( dmap.openState.xform.inverse() ) )
        points0 = numpy.copy ( points1 )
        _contour.affine_transform_vertices ( points1, dmap.data.xyz_to_ijk_transform )

        bound = 5
        li,lj,lk = numpy.min ( points1, axis=0 ) - (bound, bound, bound)
        hi,hj,hk = numpy.max ( points1, axis=0 ) + (bound, bound, bound)

        n1 = hi - li + 1
        n2 = hj - lj + 1
        n3 = hk - lk + 1

        #print " - bounds - %d %d %d --> %d %d %d --> %d %d %d" % ( li,lj,lk, hi,hj,hk, n1,n2,n3 )

        #nmat = numpy.zeros ( (n1,n2,n3), numpy.float32 )
        #dmat = dmap.full_matrix()

        nstep = (dmap.data.step[0], dmap.data.step[1], dmap.data.step[2] )
        #nstep = (fmap.data.step[0]/2.0, fmap.data.step[1]/2.0, fmap.data.step[2]/2.0 )

        nn1 = int ( round (dmap.data.step[0] * float(n1) / nstep[0]) )
        nn2 = int ( round (dmap.data.step[1] * float(n2) / nstep[1]) )
        nn3 = int ( round (dmap.data.step[2] * float(n3) / nstep[2]) )

        O = dmap.data.origin
        #print " - %s origin:" % dmap.name, O
        nO = ( O[0] + float(li) * dmap.data.step[0],
               O[1] + float(lj) * dmap.data.step[1],
               O[2] + float(lk) * dmap.data.step[2] )

        #print " - new map origin:", nO

        nmat = numpy.zeros ( (nn1,nn2,nn3), numpy.float32 )
        ndata = VolumeData.Array_Grid_Data ( nmat, nO, nstep, dmap.data.cell_angles )

        #print " - fmap grid dim: ", numpy.shape ( fmap.full_matrix() )
        #print " - new map grid dim: ", numpy.shape ( nmat )

        npoints = VolumeData.grid_indices ( (nn1, nn2, nn3), numpy.single)  # i,j,k indices
        _contour.affine_transform_vertices ( npoints, ndata.ijk_to_xyz_transform )

        dvals = dmap.interpolated_values ( npoints, dmap.openState.xform )
        #dvals = dmap.interpolated_values ( npoints, chimera.Xform.identity() )
        #dvals = dmap.interpolated_values ( npoints, dmap.openState.xform.inverse() )
        #dvals = numpy.where ( dvals > threshold, dvals, numpy.zeros_like(dvals) )
        #nze = numpy.nonzero ( dvals )

        nmat = dvals.reshape( (nn3,nn2,nn1) )
        #f_mat = fmap.data.full_matrix()
        #f_mask = numpy.where ( f_mat > fmap.surface_levels[0], numpy.ones_like(f_mat), numpy.zeros_like(f_mat) )
        #df_mat = df_mat * f_mask

        ndata = VolumeData.Array_Grid_Data ( nmat, nO, nstep, dmap.data.cell_angles )
        #try : nv = VolumeViewer.volume.add_data_set ( ndata, None )
        #except : nv = VolumeViewer.volume.volume_from_grid_data ( ndata )

        #nv.openState.xform = dmap.openState.xform

        mdata = VolumeData.zone_masked_grid_data ( ndata, points0, atomRad )
        gdata = VolumeData.Array_Grid_Data ( mdata.full_matrix(), nO, nstep, dmap.data.cell_angles, name = nname )
        nv = VolumeViewer.volume.volume_from_grid_data ( gdata )
        nv.openState.xform = dmap.openState.xform

        nv.name = nname
        dmap.display = False
        nv.region = ( nv.region[0], nv.region[1], [1,1,1] )

        if hasattr (self, 'cLevels') and dmap.name in self.cLevels :
            print "%s -- %.2f" % (dmap.name, self.cLevels[dmap.name])
            nv.surface_levels[0] = self.cLevels[dmap.name]
        else :
            nv.surface_levels[0] = dmap.surface_levels[0]

        M = dmap.data.full_matrix()
        sdev, avg, thr = numpy.std(M), numpy.average(M), nv.surface_levels[0]

        M = dmap.data.full_matrix()
        lsdev, lavg = numpy.std(nmat), numpy.average(nmat)

        #print "Avg: %.3f, sdev: %.3f, thr: %.4f [%.4f sdev above mean]" % (avg, sdev, thr, (thr-avg)/sdev)
        sigmaGlobal = (thr-avg)/sdev
        sigmaLocal = (thr-lavg)/lsdev
        umsg ( "Contour level: %.4f, %.2f/%.2f sigma above average global/local" % (thr, sigmaGlobal, sigmaLocal)  )
        #print sigmaGlobal, sigmaLocal


        ro = VolumeViewer.volume.Rendering_Options()
        ro.smoothing_factor = .3
        ro.smoothing_iterations = 2
        ro.surface_smoothing = False
        ro.square_mesh = True
        ro.line_thickness = 2
        nv.update_surface ( False, ro )
        setro (ro)
        for sp in nv.surfacePieces :
            v, t = sp.geometry
            if len(v) == 8 and len(t) == 12 :
                sp.display = False
            else :
                if showMesh :
                    sp.color = (.5, .5, .5, 1.0)
                    sp.displayStyle = sp.Mesh
                else :
                    sp.color = (0.7, 0.7, 0.7, alpha)

        return nv


    def B1_Drag (self, event):
        #print "b1m ", event.x, event.y

        if self.drag == 'seq' :
            d = event.x - self.last_x
            self.seqX += d
            #GetSegMod().seqX = self.seqX
            self.UpdateSeq ()
        elif self.drag == 'mod' :
            d = event.x - self.last_x
            self.modX += d
            #GetSegMod().modX = self.modX
            self.UpdateMod ()
        elif self.drag == 'seqSel' :
            if hasattr ( self, 'seqSel' ):
                seqI = ( event.x - self.seqX ) / self.tw
                if seqI > self.seqSel[0] :
                    self.seqSel[1] = seqI
                elif seqI < self.seqSel[1] :
                    self.seqSel[0] = seqI

                resStartI = self.seqSel[0]+1
                resEndI = self.seqSel[1]+1
                try :
                    resStartI = self.seqRes [ self.seqSel[0] ].id.position
                    resEndI = self.seqRes [ self.seqSel[1] ].id.position
                except :
                    pass
                status ( "Sequence selected %d - %d" % (resStartI,resEndI) )

                self.UpdateSeqSel ()
        elif self.drag == 'con' :
            x1, y1, x2, y2 = self.Canvas.coords ( self.conLine )
            self.Canvas.coords ( self.conLine, x1, y1, event.x, event.y )
            self.SelectedModClr ( event.x, event.y )
        elif self.drag == "modRot" :
            dx = event.x - self.dragStartX
            self.dragStartX = event.x
            self.dragMod.Rotate ( dx )



        self.last_x = event.x
        self.last_y = event.y


    def B2_Down (self, event):
        print "b2 - down"




    def B2_Up (self, event):
        print "b2  - up", event.x, event.y

        if hasattr ( self, 'selModPiece' ) and self.selModPiece != None :

            if self.selModPiece.type == "Loop" :
                self.selModPiece.MakeMod ()

            else :
                self.selModPiece.switch = not self.selModPiece.switch
                self.selModPiece.MakeMod ()
                self.UpdateMod ()


    def B2_Up_Ctrl (self, event):
        print "b2  - up - control", event.x, event.y
        if hasattr ( self, 'selModPiece' ) and self.selModPiece != None :
            if self.selModPiece.type == "Loop" :
                MakeLoopMod1 ( self.selModPiece )
                #MakeLoopMod ( self.selModPiece )



    def B2_Up_Alt (self, event):
        print "b2  - up - alt", event.x, event.y
        if hasattr ( self, 'selModPiece' ) and self.selModPiece != None :
            if self.selModPiece.type == "Loop" :
                LoopPathOpt ( self.selModPiece, self.refUseMap.get() )


    def B2_Up_Shift (self, event):
        print "b2  - up - alt", event.x, event.y
        if hasattr ( self, 'selModPiece' ) and self.selModPiece != None :
            if self.selModPiece.type == "Loop" :
                LoopPathOpt ( self.selModPiece, self.refUseMap.get() )



    def B2_Up_Comm (self, event):
        print "b2  - up - command", event.x, event.y




    def B2_Drag (self, event):
        #print "b2m ", event.x, event.y
        pass



    def B3_Down (self, event):

        print "b3 _", event.x, event.y




    def B3_Up (self, event):
        print "b3 ^", event.x, event.y
        self.B2_Up ( event )


    def B3_Drag (self, event):
        #print "b3m ", event.x, event.y
        pass


    def isInSeq ( self, x, y ) :
        if y >= self.seqY and y <= self.seqY + self.seqH :
            return True
        else :
            return False

    def isInMod ( self, x, y ) :
        if y >= self.modY and y <= self.modY + self.modH :
            return True
        else :
            return False


    def Mouse_Move (self, event):
        #print "mod m ", event.x, event.y
        #self.Canvas.coords ( self.seqMouseLine, event.x,self.seqY,event.x,self.seqY+self.seqH )

        if self.isInSeq ( event.x, event.y ) and hasattr ( self, 'seq') and len(self.seq) > 0 :

            if hasattr ( self, 'seqRec' ) and hasattr ( self, 'tw' ) and hasattr ( self, 'seqMouseR' ) :
                self.Canvas.itemconfigure ( self.seqRec, state=Tkinter.NORMAL )

                si = ( event.x - self.seqX ) / self.tw
                if si < 0 :
                    si = 0
                if si < len ( self.seq ) :
                    res = self.seqRes [ si ]
                    ri = self.seqRi [si]
                    resEnd = self.seqRes [ len(self.seqRes) - 1 ]
                    resStart = self.seqRes [ 0 ]

                    if res != None :
                        try :
                            status ( "Sequence: %s/%s %d/%d" % ( self.seq[si], res.type, res.id.position, resEnd.id.position ) )
                        except :
                            return
                    else :
                        try :
                            status ( "Sequence: ?/? %d/%d" % ( ri, resEnd.id.position ) )
                        except :
                            return

                    y0 = self.seqY+5
                    y1 = self.seqY+self.seqH-5
                    if event.y >= y0 and event.y <= y1 and hasattr ( self, 'seqMouseR' ) :
                        x0 = self.seqX + si * self.tw
                        x1 = x0 + self.tw
                        self.Canvas.coords ( self.seqMouseR, x0, y0, x1, y1 )
                        self.Canvas.itemconfigure ( self.seqMouseR, state=Tkinter.NORMAL )
                    else :
                        self.Canvas.itemconfigure ( self.seqMouseR, state=Tkinter.HIDDEN )

            else :
                self.Canvas.itemconfigure ( self.seqRec, state=Tkinter.HIDDEN )

                if hasattr ( self, 'seqMouseR' ) :
                    self.Canvas.itemconfigure ( self.seqMouseR, state=Tkinter.HIDDEN )


        self.last_x = event.x
        self.last_y = event.y


    def Canvas_Leave ( self, event ) :
        #self.Canvas.coords ( self.seqMouseLine, 0,0,0,0 )
        pass


    def Canvas_Config (self, event) :
        #print "mod cfg ", event.width, event.height
        self.W = event.width
        self.H = event.height

        #self.Canvas.delete("all")
        if 1 :
            if hasattr(self, 'backRec') :
                self.Canvas.coords (self.backRec, 0, 0, self.W, self.H)
            else :
                self.backRec = self.Canvas.create_rectangle(0, 0, self.W, self.H, outline="#eee", fill="#eee")
                #self.seqMouseLine = self.Canvas.create_line(0, 0, 0, 0, fill="#66a")

            if hasattr ( self, 'seqRec' ) :
                self.Canvas.coords ( self.seqRec, 0, self.seqY, self.W, self.seqY+self.seqH )
            else :
                self.seqRec = self.Canvas.create_rectangle(0, self.seqY, self.W, self.seqY+self.seqH, outline="#ddd", fill="#ddd" )

            self.Canvas.tag_lower(self.seqRec)
            self.Canvas.tag_lower(self.backRec)


    def Canvas_Wheel ( self, event ) :

        if self.isInSeq (self.last_x, self.last_y) :

            self.seqX += event.delta * 10

            if 0 :
                self.mag = self.mag + event.delta
                if self.mag > 15 : self.mag = 15
                if self.mag < 2 : self.mag = 2
                status ( "Mag: %d" % self.mag )

                self.font = tkFont.Font(family='Courier', size=(self.mag), weight='normal')
                #self.boldFont = tkFont.Font(family='Courier', size=(self.mag+4), weight='bold')
                self.tw = self.font.measure ( "a" )

                #GetSegMod().seqX = self.seqX
                #self.UpdateSeqFont ()

            #self.UpdateSeqFont ()
            self.UpdateSeq ()

            # ['__doc__', '__module__', 'char', 'delta', 'height', 'keycode', 'keysym', 'keysym_num', 'num', 'send_event', 'serial', 'state', 'time', 'type', 'widget', 'width', 'x', 'x_root', 'y', 'y_root']
            #print dir(event)
            #print event.delta



    def ZoomMinus ( self ) :
        self.mag = self.mag - 1
        if self.mag > 15 : self.mag = 15
        if self.mag < 2 : self.mag = 2
        #print "w ", event.delta, " mag: ", self.mag

        self.font = tkFont.Font(family='Courier', size=(self.mag), weight='normal')
        #self.boldFont = tkFont.Font(family='Courier', size=(self.mag+4), weight='bold')
        self.tw = self.font.measure ( "a" )

        self.UpdateSeqFont ()
        self.UpdateSeq ()
        status ( "Magnification: %d" % self.mag )



    def ZoomPlus ( self ) :
        self.mag = self.mag + 1
        if self.mag > 15 : self.mag = 15
        if self.mag < 2 : self.mag = 2
        #print "w ", event.delta, " mag: ", self.mag

        self.font = tkFont.Font(family='Courier', size=(self.mag), weight='normal')
        #self.boldFont = tkFont.Font(family='Courier', size=(self.mag+4), weight='bold')
        self.tw = self.font.measure ( "a" )

        self.UpdateSeqFont ()
        self.UpdateSeq ()
        status ( "Magnification: %d" % self.mag )


    def ZoomBegin ( self ) :
        self.seqX = 10
        self.UpdateSeqFont ()
        self.UpdateSeq ()

    def ZoomEnd ( self ) :
        self.seqX = - ( len(self.seq) - 50 ) * self.tw
        self.UpdateSeqFont ()
        self.UpdateSeq ()




    def isSelected ( self, fmap ) :
        for sp in fmap.surfacePieces :
            if sp in Surface.selected_surface_pieces() :
                return True
        return False




    def S_sel (self) :

        # show sigma for a side chain

        selAts = chimera.selection.currentAtoms()
        if len ( selAts ) == 0 :
            return

        dmap = self.cur_dmap

        selAtom = selAts[0]
        r = selAtom.residue
        print "Res: %s - %d.%s - %s - Atom: %s" % (r.type, r.id.position, r.id.chainId, r.molecule.name, selAtom.name)

        if 1 or not hasattr ( r.molecule, 'bbats' ) :
            qscores.SetBBAts(r.molecule)
            r.molecule.bbats = True

        removeMods = []
        for m in chimera.openModels.list() :
            if "RAD points" in m.name :
                removeMods.append ( m )
        chimera.openModels.remove ( removeMods )


        ats = [at for at in self.cur_mol.atoms if not at.element.name == "H"]
        points = _multiscale.get_atom_coordinates ( ats, transformed = False )
        print " - search tree: %d/%d ats" % ( len(ats), len(r.molecule.atoms) )
        allAtTree = AdaptiveTree ( points.tolist(), ats, 1.0)

        #allAtTree = None
        #print "-"

        import time
        start = time.time()

        #sigma = RadAts ( r.scAtoms, dmap, allAtTree=allAtTree, show=1, log=1, numPts=30, toRAD=2, dRAD=0.5 )

        if 0 :
            print "_sigma____________________________"
            sigma = RadAts ( [selAtom], dmap, allAtTree=allAtTree, show=1, log=1, numPts=30, toRAD=2, dRAD=0.1 )
            res = sigma * numpy.pi * numpy.sqrt(2.0)
            end = time.time()
            print "%s - sigma: %.3f, res: %.3f, time: %f" % ( selAtom.name, sigma, res, (end - start) )


        minD, maxD = qscores.MinMaxD ( dmap )
        print " - mind: %.3f, maxd: %.3f" % (minD, maxD)
        #sigma = 0.6

        sigma = float(self.sigma.get())

        start = time.time()
        qq = qscores.Qscore  ( [selAtom], dmap, sigma, allAtTree=allAtTree, show=1, log=1, numPts=20, toRAD=3.0, dRAD=0.5, minD=minD, maxD=maxD, fitg=0 )
        end = time.time()
        print " - time: %f" % ( (end - start) )

        start = time.time()
        qq = qscores.Qscore  ( [selAtom], dmap, sigma, allAtTree=allAtTree, show=0, log=1, numPts=50, toRAD=2.0, dRAD=0.1, minD=minD, maxD=maxD, fitg=1 )
        end = time.time()
        print " - time: %f" % ( (end - start) )

        #start = time.time()
        #qq = qscores.Qscore ( [selAtom], dmap, sigma, allAtTree=allAtTree, show=0, log=1, numPts=20, toRAD=3.0, dRAD=0.1, minD=minD, maxD=maxD, fitg=1 )
        #end = time.time()
        #print " - time: %f" % ( (end - start) )

        #CC, CCm, yds, err = rr

        print r


    def ResCC ( self ) :

        if self.cur_mol == None :
            umsg ( "Select a model" )
            return

        if self.cur_dmap == None :
            umsg ( "Select a map" )
            return

        res = float ( self.mapRes.get() )
        chainId = self.chain.get()

        from qscores import ResCC
        ResCC ( self.cur_dmap, self.cur_mol, chainId, res )

        self.ShowQScores()

    def Q_sel (self, show=False, fitG=False) :

        # show sigma for a side chain

        selAts = chimera.selection.currentAtoms()
        if len ( selAts ) == 0 :
            return

        dmap = self.cur_dmap

        selAtom = selAts[0]
        selRes = selAtom.residue
        selMol = selAtom.molecule

        print ""
        print "Res: %s - %d.%s - %s - Atom: %s" % (selRes.type, selRes.id.position, selRes.id.chainId, selMol.name, selAtom.name)
        print " - in map: %s" % dmap.name

        qscores.SetBBAts(selMol)

        removeMods = []
        for m in chimera.openModels.list() :
            if "RAD points" in m.name :
                removeMods.append ( m )
            if "qp points" in m.name :
                removeMods.append ( m )
        chimera.openModels.remove ( removeMods )

        sigma = float(self.sigma.get())

        ats = None
        points = []

        if sigma < 1.0 :
            ats = [at for at in selMol.atoms if not at.element.name == "H"]
            print " - %d all non-H atoms (for sigma < 1)" % len(ats)
            points = _multiscale.get_atom_coordinates ( ats, transformed = False )

        else :
            ats = [at for at in selMol.atoms if (at.residue.isProt and at.name == "CA")]
            print " - %d protein CA atoms for sigma > 1" % len(ats)
            if len(ats) > 0 :
                points = _multiscale.get_atom_coordinates ( ats, transformed = False )
                pmod = qscores.AddSpherePts ( points, (.9,.9,.3,.4), 1.0, "%s - prot qp points" % (selRes.molecule.name) )
                pmod.openState.xform = selAtom.molecule.openState.xform

            pts1 = [res.qp_P for res in selMol.residues if (res.isNA and res.qp_P != None) ]
            pts2 = [res.qp_Sugar for res in selMol.residues if (res.isNA and res.qp_Sugar != None) ]
            pts3 = [res.qp_Base for res in selMol.residues if (res.isNA and res.qp_Base != None) ]
            print " - %d/%d/%d P/Sugar/Base nucleic atoms for sigma > 1" % ( len(pts1), len(pts2), len(pts3) )

            if len(pts1) > 0 :
                points.extend ( pts1 )
                pmod = qscores.AddSpherePts ( pts1, (.9,.3,.3,.4), 1.0, "%s - na qp points" % (selRes.molecule.name) )
                pmod.openState.xform = selAtom.molecule.openState.xform
            if len(pts2) > 0 :
                points.extend ( pts2 )
                pmod = qscores.AddSpherePts ( pts2, (.3,.9,.3,.4), 1.0, "%s - qp points" % (selRes.molecule.name) )
                pmod.openState.xform = selAtom.molecule.openState.xform
            if len(pts3) > 0 :
                points.extend ( pts3 )
                pmod = qscores.AddSpherePts ( pts3, (.3,.3,.9,.4), 1.0, "%s - qp points" % (selRes.molecule.name) )
                pmod.openState.xform = selAtom.molecule.openState.xform

            print ( "%d total points" % len(points) )

        atPt = None
        if sigma < 1.0 :
            atPt = selAtom.coord().data()
        else :
            if selRes.isProt :
                if "CA" in selRes.atomsMap :
                    atPt = selRes.atomsMap["CA"][0].coord().data()
            elif selRes.isNA :
                if selAtom.isBB :
                    atPt = selRes.qp_P
                elif selAtom.isSugar :
                    atPt = selRes.qp_Sugar
                elif selAtom.isBase :
                    atPt = selRes.qp_Base

        if atPt == None :
            print " - no selected point??"
            return

        gridD = 3.0 if sigma < 1.0 else 6.0

        ptGrid, atGrid = None, None
        import gridm
        reload(gridm)
        ptGrid = gridm.Grid()
        ptGrid.FromPoints ( points, gridD )
        print " - %d pts grid - %.2f" % ( len(points), gridD )

        #atGrid = gridm.Grid()
        #atGrid.FromAtomsLocal ( ats, 3.0 )
        #print " - %d ats grid in %.3f sec" % (len(ats), 0.0)

        minD, maxD = qscores.MinMaxD ( dmap )
        print " - mind: %.3f, maxd: %.3f" % (minD, maxD)

        import time
        start = time.time()

        print ""
        print "_Q_score____________________________"

        xfI = selAtom.molecule.openState.xform
        start = time.time()

        toRad = 3.0 if sigma < 1.0 else 6.0
        dRad = 0.1 if sigma < 1.0 else 0.2
        dRadShow = 0.5 if sigma < 1.0 else 1.0

        if not show :
            qs = qscores.QscorePt3 ( atPt, xfI, dmap, sigma, ptGrid=ptGrid, log=0, numPts=8, toRAD=toRad, dRAD=dRad, minD=minD, maxD=maxD, fitg=0 )
            end = time.time()
            print " - sigma (pt): %.3f, Q-score: %.3f, time: %f" % ( sigma, qs, (end - start) )
        elif fitG :
            qs = qscores.QscorePt3 ( atPt, xfI, dmap, sigma, ptGrid=ptGrid, log=1, numPts=10, toRAD=toRad, dRAD=dRad, minD=minD, maxD=maxD, fitg=1, show=0 )
            end = time.time()
            print " - sigma (pt): %.3f, Q-score: %.3f, time: %f" % ( sigma, qs[0], (end - start) )
        else :
            qs = qscores.QscorePt3 ( atPt, xfI, dmap, sigma, ptGrid=ptGrid, log=1, numPts=100, toRAD=toRad, dRAD=dRadShow, minD=minD, maxD=maxD, fitg=0, show=1 )
            end = time.time()
            print " - sigma (pt): %.3f, Q-score: %.3f, time: %f" % ( sigma, qs, (end - start) )





    def CalcSseQAll (self) :

        if self.cur_mol == None :
            umsg ( "Select a model" )
            return

        if self.cur_dmap == None :
            umsg ( "Select a map" )
            return

        #qscores.sseQscores ( self.cur_mol, self.cur_dmap, 3.0 )

        qscores.sseQscores2 ( self.cur_mol, self.cur_dmap, 3.0 )


    def CalcSseQAllMaps (self) :

        if self.cur_mol == None :
            umsg ( "Select a model" )
            return

        if self.cur_dmap == None :
            umsg ( "Select a map" )
            return

        for m in OML(modelTypes = [VolumeViewer.volume.Volume]) :
            if m.display == True :
                ts = m.name.split("_")[-1]
                res = ts.split (".")[0][0:-1]
                print ""
                print m.name, res
                print ""
                ssQ, ssQMax = qscores.sseQscores ( self.cur_mol, m, 1.0 )
                print ssQ
                fp = open ( "/Users/greg/Desktop/ssQ.txt", "a" )
                fp.write ( "%s\t%s\t%.3f\t%.3f\n" % (m.name, res, ssQ, ssQMax) )
                fp.close()

        qscores.sseQscores ( self.cur_mol, self.cur_dmap, 3.0 )


    def CalcSseQSel (self) :

        selRess = chimera.selection.currentResidues()
        if len ( selRess ) == 0 :
            return

        selRes = selRess[0]

        if self.cur_dmap == None :
            umsg ( "Select a map" )
            return

        print " - %d sel res" % len(selRess)

        minD, maxD = qscores.MinMaxD (self.cur_dmap)

        mol = selRes.molecule
        allRess = []
        for r in mol.residues :
            if r.id.chainId == selRes.id.chainId :
                allRess.append ( [r.id.position, r] )
        allRess.sort()
        allRess = [r for ri,r in allRess]

        #for r in allRess[:30] :
        #    print r.type, r.id.position

        sses = qscores.SSEs ( allRess )

        print " - %d sses, %d res" % ( len(sses), len(allRess) )

        for startRi, endRi, ss, ress in sses :
            if ss == "H" :
                if selRes.id.position >= startRi and selRes.id.position <= endRi :
                    print " - found sel sse, %d-%d, %d res" % (startRi, endRi, len(ress))
                    rr = qscores.sseQscore ( ress, "H", self.cur_dmap, 3.0, showRes=selRes, log=1, numPts=6, toRAD=3.0, dRAD=0.1, minD=minD, maxD=maxD, fitg=0 )
                    #rr = qscores.sseQscore ( ress, "H", self.cur_dmap, 3.0, showRes=None, log=1, numPts=6, toRAD=3.0, dRAD=0.1, minD=minD, maxD=maxD, fitg=1 )



    def Q_show (self) :

        # show sigma for a side chain

        selAts = chimera.selection.currentAtoms()
        if len ( selAts ) == 0 :
            return

        dmap = self.cur_dmap


        selAtom = selAts[0]
        r = selAtom.residue
        print ""
        print "Res: %s - %d.%s - %s - Atom: %s" % (r.type, r.id.position, r.id.chainId, r.molecule.name, selAtom.name)

        print " - in map: %s" % dmap.name

        if 1 or not hasattr ( r.molecule, 'bbats' ) :
            qscores.SetBBAts(r.molecule)
            r.molecule.bbats = True

        removeMods = []
        for m in chimera.openModels.list() :
            if "RAD points" in m.name :
                removeMods.append ( m )
        #chimera.openModels.remove ( removeMods )


        ats = [at for at in self.cur_mol.atoms if not at.element.name == "H"]
        if self.showH.get() :
            ats = self.cur_mol.atoms

        points = _multiscale.get_atom_coordinates ( ats, transformed = False )
        print " - search tree: %d/%d ats" % ( len(ats), len(r.molecule.atoms) )
        allAtTree = AdaptiveTree ( points.tolist(), ats, 1.0)

        sigma = float(self.sigma.get())

        #allAtTree = None
        #print "-"

        import time
        start = time.time()

        #sigma = RadAts ( r.scAtoms, dmap, allAtTree=allAtTree, show=1, log=1, numPts=30, toRAD=2, dRAD=0.5 )

        if 0 :
            print "_sigma____________________________"

            sigma = RadAts ( [selAtom], dmap, allAtTree=allAtTree, show=1, log=1, numPts=30, toRAD=2, dRAD=0.5 )
            res = sigma * numpy.pi * numpy.sqrt(2.0)

            end = time.time()
            print "%s - sigma: %.3f, res: %.3f, time: %f" % ( selAtom.name, sigma, res, (end - start) )

        elif 1 :
            print "_Q_score____________________________"

            minD, maxD = qscores.MinMaxD ( dmap )
            print " - mind: %.3f, maxd: %.3f" % (minD, maxD)


            if 0 :
                minD = numpy.min(M)
                print " - min before masking: %.4f" % minD
                points = _multiscale.get_atom_coordinates ( ats, transformed = False )
                _contour.affine_transform_vertices ( points, Matrix.xform_matrix( dmap.openState.xform.inverse() ) )
                mdata = VolumeData.zone_masked_grid_data ( dmap.data, points, 0.5 )
                M = mdata.full_matrix ()
                minD = numpy.min(M)
                print " - min after masking: %.4f" % minD
                M = numpy.where ( M == 0.0, numpy.ones_like(M)*(minD-0.2), M )
                import _volume
                points = _volume.high_indices(M, minD-0.1)
                fpoints = points.astype(numpy.single)
                fpoint_weights = M[points[:,2],points[:,1],points[:,0]]
                minD = numpy.min(fpoint_weights)
                print " - min of mask pts: %.4f" % minD


            #sigma = 2.0 / (numpy.pi * numpy.sqrt(2.0))
            #sigma = 0.4

            qs, yds, err = 0,0,0

            if 0 :
                rr = qscores.Qscore ( [selAtom], dmap, sigma, allAtTree=allAtTree, show=1, log=1, numPts=20, toRAD=2.0, dRAD=0.5, minD=minD, maxD=maxD, fitg=1 )
                qs, yds, err = rr

            elif 1 :
                rr = qscores.Qscore ( [selAtom], dmap, sigma, allAtTree=allAtTree, show=0, log=1, numPts=30, toRAD=2.0, dRAD=0.1, minD=minD, maxD=maxD, fitg=0 )
                qs, yds, err = rr

            else :
                qs = qscores.Qscore ( [selAtom], dmap, sigma, allAtTree=allAtTree, show=0, log=0, numPts=8, toRAD=2.0, dRAD=0.1, minD=minD, maxD=maxD, fitg=0 )

            end = time.time()
            print " - sigma: %.3f, Q-score: %.3f, time: %f" % ( sigma, qs, (end - start) )

            print "Atoms in %d.%s %s" % (selAtom.residue.id.position, selAtom.residue.id.chainId, selAtom.residue.type)
            #print "-"


            if 0 :
                avg, N = 0.0, 0.0
                #bbAts, scAts, baseAts, sugarAts = [], [], [], []

                for at in selAtom.residue.atoms :
                    at.Q = qscores.Qscore ( [at], dmap, sigma, allAtTree=allAtTree, show=0, log=0, numPts=8, toRAD=2.0, dRAD=0.1, minD=minD, maxD=maxD )
                    print " - %s : %.2f" % (at.name, at.Q)
                    if 1 or at.isSC :
                        avg += at.Q
                        N += 1.0

                    #if at.residue.isNA and at.isBB : bbAts.append ( at )
                    #if at.residue.isNA and at.isSugar : sugarAts.append ( at )
                    #if at.residue.isNA and at.isBase : baseAts.append ( at )

                    #if at.residue.isProt and at.isBB : bbAts.append ( at )
                    #if at.residue.isProt and at.isSC : scAts.append ( at )


                if selAtom.residue.isNA :
                    print "NA:"
                    print " - backbone Q: %.2f" % numpy.average ( [at.Q for at in selAtom.residue.bbAtoms] )
                    #print " - sugar Q: %.2f" % numpy.average ( [at.Q for at in sugarAts] )
                    print " - base Q: %.2f" % numpy.average ( [at.Q for at in selAtom.residue.scAtoms] )

                if selAtom.residue.isProt :
                    print "Protein:"
                    print " - backbone Q: %.2f" % numpy.average ( [at.Q for at in selAtom.residue.bbAtoms] )
                    print " - side chain Q: %.2f" % numpy.average ( [at.Q for at in selAtom.residue.scAtoms] )

                if N > 0 :
                    print "All:"
                    #print " - avg sc Q: %.2f" % (avg/N)
                    print " - avg Q: %.2f" % (avg/N)


    def CalcSelR (self) :

        atoms = chimera.selection.currentAtoms()
        if len ( atoms ) == 0 :
            umsg ( "No selected atoms found" )
            return

        dmap = self.cur_dmap
        mol = atoms[0].molecule

        #umsg ( "Calculating Q-scores of %d atoms..." % len(atoms) )

        print "R-score"
        print " -", dmap.name
        print " -", mol.name, " - %d atoms" % len(atoms)


        #qscores.RScore ( atoms, dmap )
        qscores.RScore ( atoms, dmap, minRes=1.0, maxRes=15.0, dRes = 0.2, xf=None )


    # .....
    def CalcSelQ (self) :

        # show sigma for a side chain

        atoms = chimera.selection.currentAtoms()
        if len ( atoms ) == 0 :
            umsg ( "No selected atoms found" )
            return

        dmap = self.cur_dmap
        mol = atoms[0].molecule

        umsg ( "Calculating Q-scores of %d atoms..." % len(atoms) )

        sigma = float(self.sigma.get())

        #sigma = 0.4
        print " - in map: %s" % self.cur_dmap.name
        print " - mol: %s" % mol.name
        print " - sigma: %.2f" % sigma

        if 1 or not hasattr ( mol.name, 'bbats' ) :
            qscores.SetBBAts(mol)
            mol.bbats = True

        ats = [at for at in mol.atoms if not at.element.name == "H"]
        if self.showH.get() :
            ats = mol.atoms

        points = _multiscale.get_atom_coordinates ( ats, transformed = False )
        print " - search tree: %d/%d ats" % ( len(ats), len(mol.atoms) )
        allAtTree = AdaptiveTree ( points.tolist(), ats, 1.0)

        minD, maxD = qscores.MinMaxD ( dmap )
        print " - minD %.3f, maxD %.3f" % (minD, maxD)

        import time
        start = time.time()

        from chimera import tasks, CancelOperation
        task = tasks.Task('Calculating Q-scores', modal = True)

        avg, avgBB, avgSC, numBB, numSC = 0.0, 0.0, 0.0, 0, 0
        avgSug, numSug = 0.0,0.0

        import traceback

        nonhatoms = []
        try :

            for ai, at in enumerate ( atoms ) :

                if at.element.name == "H" :
                    continue

                nonhatoms.append ( at )
                at.Q = qscores.Qscore ( [at], dmap, sigma, allAtTree=allAtTree, show=0, log=0, numPts=8, toRAD=2.0, dRAD=0.1, minD=minD, maxD=maxD, fitg=0 )
                at.bfactor = at.Q
                avg += at.Q

                if at.isBB :
                    avgBB += at.Q
                    numBB += 1
                if at.isSC :
                    avgSC += at.Q
                    numSC += 1
                if at.isSugar :
                    avgSug += at.Q
                    numSug += 1

                if (ai+1) % 10 == 0 :
                    leftTime = qscores.TimeLeftStr (ai, len(atoms), time.time() - start)
                    status ( "Calculating Q scores - atom %d/%d - eta: %s" % (ai+1, len(atoms), leftTime) )
                    print ".",
                    #task.updateStatus( "Calculating Q scores - atom %d/%d - %s in %s.%d.%s - eta: %s" % (ai+1, len(atoms), at.name, at.residue.type, at.residue.id.position, at.residue.id.chainId, leftTime) )
                    task.updateStatus( "Calculating Q scores - atom %d/%d - eta: %s" % (ai+1, len(atoms), leftTime) )

        except Exception, err:
            umsg ( "Something went wrong..." )
            print Exception, err
            traceback.print_exc()
            return

        finally :
            task.finished()

        for at in nonhatoms :
            print " - atom: %s %d.%s %s : %.3f" % (at.residue.type, at.residue.id.position, at.residue.id.chainId, at.name, at.Q)

        if numBB > 0 :
            print "%d BB/Pho atoms, Q=%.3f" % (numBB, avgBB/float(numBB))
        if numSC > 0 :
            print "%d SC/Base atoms, Q=%.3f" % (numSC, avgSC/float(numSC))
        if numSug > 0 :
            print "%d Sugar, Q=%.3f" % (numSug, avgSug/float(numSug))

        avgq = avg / float(len(nonhatoms))
        if len(atoms) > 1 :
            umsg ( "Q-score of %d atoms: %.3f" % (len(nonhatoms), avgq) )
        else :
            umsg ( "Q-score of %d atom: %.3f" % (len(nonhatoms), avgq) )

        res = float(self.mapRes.get())
        cc1, cc2 = ResCC ( mol, atoms, res, dmap )

        cc1_, cc2_ = ResCC ( mol, nonhatoms, res, dmap )

        numClash = 0.0
        if 1 and allAtTree :
            for at in nonhatoms :
                anear = allAtTree.searchTree ( at.coord().data(), 2.0 )
                for nat in anear :
                    if nat.residue != at.residue :
                        v = at.coord() - nat.coord()
                        if v.length < 1.8 :
                            numClash += 1.0
                            break

        clashScore = numClash / float(len(nonhatoms))

        res = float ( self.mapRes.get() )
        sigma = float ( self.sigma.get() )
        expQ, lowQ, highQ, eqn = qscores.ExpectedQScore ( res, sigma )
        print ( "For res %.2f, sigma %.1f, expected Q-score: %.3f (%.3f -- %.3f)" % (res, sigma, expQ, lowQ, highQ) )

        if 0 :
            # for some ligands stats...
            R = atoms[0].residue
            rid = "%s.%d.%s" % (R.type, R.id.position, R.id.chainId)

            #print mol.name, rid, avgq, cc1, cc2, cc1_, cc2_, res

            print "\nMol Name\tRes Id\tQ\tCC\tCCm\tCC(noh)\tCCm(noh)\tClash\tClashes"
            print "%s\t%s\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.2f\t%.0f (%.2f)\n" % (mol.name, rid, avgq, cc1, cc2, cc1_, cc2_, clashScore, numClash, res)

            if not os.path.isfile ("/Users/greg/Desktop/txt.txt") :
                fp = open ( "/Users/greg/Desktop/txt.txt", "a" )
                fp.write ( "Mol Name\tRes Id\tQ\tCC\tCCm\tCC(noh)\tCCm(noh)\tClash\tClashes\n" )
                fp.close()
            fp = open ( "/Users/greg/Desktop/txt.txt", "a" )
            fp.write ( "%s\t%s\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.2f\t%.0f\n" % (mol.name, rid, avgq, cc1, cc2, cc1_, cc2_, clashScore, numClash) )
            fp.close()








    def CalcSelQOpen ( self ) :

        for m in chimera.openModels.list() :

            if type(m) != chimera.Molecule :
                continue

            resN = None
            for r in m.residues :

                if r.type == "PTQ" :
                #if r.type == "PEE" or r.type == "ACB" :
                #if r.type == "F86" :
                    resN = r
                    break

            if resN :

                print "\n\n-------------- %s -------- %s.%d.%s" % (m.name, resN.type, resN.id.position, resN.id.chainId)

                chimera.selection.clearCurrent ()
                chimera.selection.addCurrent ( resN.atoms )
                self.CalcSelQ ()





    def AProfs (self) :

        mol = self.cur_mol
        if self.cur_mol == None :
            umsg ("Select a molecule first")
            return []

        chainId = self.chain.get()

        dmap = self.cur_dmap
        print " - in map: %s" % dmap.name

        if 1 or not hasattr ( mol, 'bbats' ) :
            qscores.SetBBAts(mol)
            mol.bbats = True

        ats = [at for at in self.cur_mol.atoms if not at.element.name == "H"]
        points = _multiscale.get_atom_coordinates ( ats, transformed = False )
        print " - search tree: %d/%d ats" % ( len(ats), len(mol.atoms) )
        allAtTree = AdaptiveTree ( points.tolist(), ats, 1.0)

        #sigma = 0.4
        minD, maxD = qscores.MinMaxD ( dmap )
        print " - mind: %.3f, maxd: %.3f" % (minD, maxD)

        sigma = float(self.sigma.get())

        def doAt (at, arr) :
            rr = qscores.Qscore ( [at], dmap, sigma, allAtTree=allAtTree, show=0, log=0, numPts=10, toRAD=3.0, dRAD=0.1, minD=minD, maxD=maxD, fitg=1 )
            Qscore, yds, err = rr
            #print len(yds)
            #if len(yds) == 31 :
            arr.append ( [Qscore,err] + yds.tolist() )
            if 0 :
                print "%.3f\t%.5f\t%s.%d.%s" % (at.Q, err, at.residue.type, at.residue.id.position, at.name),
                for y in yds : print "\t%f" % y,
                print ""
            else :
                print ".",


        bb_atn_q, sc_atn_q = {}, {}
        for at in mol.atoms :
            if at.residue.isProt and (at.name == 'C' or at.name == 'O' or at.name == 'N' or at.name == "CA") :
                if at.name in bb_atn_q :
                    bb_atn_q[at.name].append ( [at.Q, at] )
                else :
                    bb_atn_q[at.name] = [[at.Q, at]]

            atn = "%s(%s)" % (at.residue.type, at.name)
            if atn in sc_atn_q :
                sc_atn_q[atn].append ( [at.Q, at] )
            else :
                sc_atn_q[atn] = [[at.Q, at]]



        N = 60
        print "N = %d" % N

        # BB
        bb_c, bb_n, bb_o, bb_ca = [], [], [], []
        if 0 :
            for an, aa in [ ["C",bb_c], ["N",bb_n], ["O",bb_o], ["CA",bb_ca]] :
                print "___",an,"___";
                A = bb_atn_q[an];
                A.sort ( reverse=True, key=lambda x: x[0] )
                print "%d - " % len(A);
                i = 0
                for q, at in A[:N] :
                    if q > 0.8 :
                        doAt (at, aa)
                    i += 1;
                    print "%d" % i,
                print ""


        # SC
        asp_o, glu_o, arg_n, leu_c, val_c = [], [], [], [], []
        if 0 :
            for an, aa in [ ["ASP(OD1)",asp_o], ["ASP(OD2)",asp_o]] :
                print "___",an,"___"; A = sc_atn_q[an]; A.sort ( reverse=True, key=lambda x: x[0] )
                for q, at in A[0:N] :
                    if q > 0.8 : doAt (at, aa)
                print ""
            for an, aa in [ ["GLU(OE1)",glu_o], ["GLU(OE1)",glu_o]] :
                print "___",an,"___"; A = sc_atn_q[an]; A.sort ( reverse=True, key=lambda x: x[0] )
                for q, at in A[0:N] :
                    if q > 0.8 : doAt (at, aa)
                print ""
            for an, aa in [ ["ARG(NH1)",arg_n], ["ARG(NH2)",arg_n]] :
                print "___",an,"___"; A = sc_atn_q[an]; A.sort ( reverse=True, key=lambda x: x[0] )
                for q, at in A[0:N] :
                    if q > 0.8 : doAt (at, aa)
                print ""
            for an, aa in [ ["LEU(CD1)",leu_c], ["LEU(CD2)",leu_c]] :
                print "___",an,"___"; A = sc_atn_q[an]; A.sort ( reverse=True, key=lambda x: x[0] )
                for q, at in A[0:N] :
                    if q > 0.8 : doAt (at, aa)
                print ""
            for an, aa in [ ["VAL(CG1)",val_c], ["VAL(CG2)",val_c]] :
                print "___",an,"___"; A = sc_atn_q[an]; A.sort ( reverse=True, key=lambda x: x[0] )
                for q, at in A[0:N] :
                    if q > 0.8 : doAt (at, aa)
                print ""

        # HOH, ion
        hoh_o, i_i = [], []
        if 1 :
            for an, aa in [ ["HOH(O)",hoh_o], ["MG(MG)",i_i]] :
                print "___",an,"___"; A = sc_atn_q[an]; A.sort ( reverse=True, key=lambda x: x[0] )
                for q, at in A[0:N] :
                    if q > 0.8 : doAt (at, aa)
                print ""


        for r in mol.residues :
            if 0 :
                if r.type == "ASP" :
                    for at in [r.atomsMap["OD1"][0], r.atomsMap["OD2"][0]] :
                        if at.Q > 0.8 : doAt (at, asp_o)
                if r.type == "GLU" :
                    for at in [r.atomsMap["OE1"][0], r.atomsMap["OE2"][0]] :
                        if at.Q > 0.8 : doAt (at, glu_o)
            if 0 :
                if r.type == "VAL" :
                    for at in [r.atomsMap["CG1"][0], r.atomsMap["CG2"][0]] :
                        if at.Q > 0.8 : doAt (at, val_c)
                if r.type == "LEU" :
                    for at in [r.atomsMap["CD1"][0], r.atomsMap["CD2"][0]] :
                        if at.Q > 0.8 : doAt (at, leu_c)
            if 0 :
                if r.type == "ARG" :
                    for at in [r.atomsMap["NH1"][0], r.atomsMap["NH2"][0]] :
                        if at.Q > 0.8 : doAt (at, arg_n)

                #if r.type == "LEU" :
                #    for at in [r.atomsMap["CD1"][0], r.atomsMap["CD2"][0]] :
                #        if at.Q > 0.8 : doAt (at, leu_c)


        def outAt (arr, label, w="avg") :

            arr.sort ( reverse=True, key=lambda x: x[0] )

            #K = 10
            #aa = numpy.array ( arr[0:K] )
            aa = numpy.array ( arr )

            if w=="p" :
                print "Q\tAvgD - ", label
                for qa in aa :
                    for d in qa :
                        print "%f\t" % d,
                    print ""
                return

            s = numpy.std(aa,axis=0)
            m = numpy.mean(aa,axis=0)

            print label, "\t", aa.shape,

            #print label,
            if w == "avg" :
                for i in range(len(s)) :
                    print "\t%f" % m[i],
            else :
                for i in range(len(s)) :
                    print "\t%f" % s[i],

            print ""


        print ""
        print "Res\tQ\tErr",
        for yi in range(31) : print "\t%f" % (yi*.1),
        print ""

        if 0 :
            for w in ["p", "avg", "std"] :
                outAt ( val_c, "VAL(CG)", w )
                outAt ( leu_c, "LEU(CD)", w )
                outAt ( arg_n, "ARG(NH)", w )
                outAt ( asp_o, "ASP(OD)", w )
                outAt ( glu_o, "GLU(OE)", w )
                print ""

        if 0 :
            for w in ["p", "avg", "std"] :
                outAt ( bb_c, "C", w )
                outAt ( bb_ca, "CA", w )
                outAt ( bb_n, "N", w )
                outAt ( bb_o, "O", w )

        if 1 :
            for w in ["p", "avg", "std"] :
                outAt ( hoh_o, "Water(O)", w )
                outAt ( i_i, "Ion", w )



    def Ligs ( self ) :

        mol = self.cur_mol
        if self.cur_mol == None :
            umsg ("Select a molecule first")
            return []

        chainId = self.chain.get()

        dmap = self.cur_dmap
        print " - in map: %s" % dmap.name

        if 1 or not hasattr ( mol, 'bbats' ) :
            qscores.SetBBAts(mol)
            mol.bbats = True

        ats = [at for at in self.cur_mol.atoms if not at.element.name == "H"]
        points = _multiscale.get_atom_coordinates ( ats, transformed = False )
        print " - search tree: %d/%d ats" % ( len(ats), len(mol.atoms) )
        allAtTree = AdaptiveTree ( points.tolist(), ats, 1.0)

        #sigma = 0.6
        minD, maxD = qscores.MinMaxD ( dmap )
        print " - mind: %.3f, maxd: %.3f" % (minD, maxD)


        remAts = []
        showAts = []
        for r in mol.residues :
            if not r.isProt and not r.isNA :
                #print at.residue.id.position, at.residue.type,
                for at in r.atoms :
                    #print " - %s - %.2f" % (at.name, at.Q)

                    if hasattr ( at, 'Q1' ) and hasattr ( at, 'Q2' ) :
                        if at.Q > 0.8 and at.Q2 > 0.8 and at.Q1 > 0.8 :
                            print " -3- %s - %.2f, %.2f, %.2f" % (at.name, at.Q, at.Q1, at.Q2)
                            at.display = True
                            showAts.append ( at )
                        else :
                            at.display = False
                            remAts.append ( at )

                    else :
                        #print " -1- %s - %.2f" % (at.name, at.Q)
                        if at.Q > 0.8 :
                            at.display = True
                            showAts.append ( at )
                        else :
                            at.display = False
                            remAts.append ( at )

        print "Showing %d, Hiding %d" % ( len(showAts), len(remAts) )

        if 0 :
            for at in remAts :
                mol.deleteAtom ( at )




    def Scale ( self ) :

        mol = self.cur_mol
        if self.cur_mol == None :
            umsg ("Select a molecule first")
            return []

        chainId = self.chain.get()

        dmap = self.cur_dmap
        print " - scale map: %s" % dmap.name


        def SetStep ( D, S ) :
            # Preserve index origin.
            index_origin = D.data.xyz_to_ijk((0,0,0))
            #print " - origin 0 :", index_origin
            D.data.set_step ( S )
            xyz_origin = [x0-x for x0,x in zip(D.data.ijk_to_xyz((0,0,0)),D.data.ijk_to_xyz(index_origin))]
            D.data.set_origin(xyz_origin)
            #print " - origin 1 :", xyz_origin


        xf0 = mol.openState.xform
        s0 = dmap.data.step

        max_Xf, max_Avg, max_S = None, None, None

        RES = 2.0 # float(self.mapRes.get())
        print " - res: %.1f" % RES

        vals = []

        avg, cc, ccm = FitMolToMap ( mol, dmap, RES )
        MapUp ( dmap, showMesh = False, color=(.7,.7,.7,1) )

        print "Initial avg: %f, step: %.4f" % (avg, dmap.data.step[0])
        vals.append ( [avg, cc, ccm, dmap.data.step[0], dmap.data.step[1], dmap.data.step[2]] )


        D = 0.0001
        for i in range ( 100 ) :

            S = dmap.data.step
            S_ = ( S[0] - D, S[1] - D, S[2] - D )
            SetStep ( dmap, S_  )

            avg, cc, ccm = FitMolToMap ( mol, dmap, RES )
            MapUp ( dmap, showMesh = False, color=(.7,.7,.7,1) )

            vals.append ( [avg, cc, ccm, S_[0], S_[1], S_[2]] )

            #print " %f - %f" % (S_[0], cc)

            if i % 10 == 0 :
                status ( "Step: %f - %d/%d" % (S_[0], i+1, 100) )
                print ".",

            if max_Avg == None or max_Avg < avg :
                max_Avg = avg
                max_Xf = mol.openState.xform
                max_S = S_

        print ""

        vals.reverse ()

        SetStep ( dmap, s0  )
        MapUp ( dmap, showMesh = False, color=(.7,.7,.7,1) )
        mol.openState.xform = xf0

        for i in range ( 100 ) :

            S = dmap.data.step
            S_ = ( S[0] + D, S[1] + D, S[2] + D )
            SetStep ( dmap, S_ )

            avg, cc, ccm = FitMolToMap ( mol, dmap, RES )
            MapUp ( dmap, showMesh = False, color=(.7,.7,.7,1) )

            vals.append ( [avg, cc, ccm, S_[0], S_[1], S_[2]] )

            #print " %f - %f" % (S_[0], cc)

            if i % 10 == 0 :
                status ( "Step: %f + %d/%d" % (S_[0], i+1, 100) )
                print ".",

            if max_Avg == None or max_Avg < avg :
                max_Avg = avg
                max_Xf = mol.openState.xform
                max_S = S_


        print ""
        print "Max avg: %f, step: %.4f" % (max_Avg, max_S[0])
        SetStep ( dmap, max_S  )
        MapUp ( dmap, showMesh = False, color=(.7,.7,.7,1) )
        mol.openState.xform = max_Xf

        nout = os.path.splitext(dmap.data.path)[0] + "_scales.txt"
        fout = open ( nout, "w" )
        for avg, cc, ccm, s0, s1, s2 in vals :
            fout.write ( "%f\t%f\t%f\t%f\t%f\t%f\n" % (s0, s1, s2, avg, cc, ccm) )
        fout.close()
        print " - wrote ", nout





    def Asn ( self ) :

        print "ASN - show"

        mol = self.cur_mol
        if self.cur_mol == None :
            umsg ("Select a molecule first")
            return []

        chainId = self.chain.get()

        dmap = self.cur_dmap
        print " - scale map: %s" % dmap.name


        totAt, showAt = 0, 0

        tot = {}
        rids = []


        for r in self.cur_mol.residues :

            if r.id.chainId != chainId :
                continue

            rids.append ( [r.id.position, r] )

        rids.sort ()

        i = 0
        for ri, r in rids :
            if i > 2 :
                r2 = rids[i-2]
                if (r.type == "SER" or r.type == "THR") and r2[1].type == "ASN" :
                    print "%s - %d.%s" % (r2[1].type, r2[1].id.position, r2[1].id.chainId)

                    chimera.selection.addCurrent ( r2[1] )


            i += 1



    def DelSel ( self ) :

        mol = chimera.selection.currentMolecules()[0]

        for b in chimera.selection.currentBonds() :
            mol.deleteBond(b)

        for at in chimera.selection.currentAtoms() :
            mol.deleteAtom(at)


    def Take ( self ) :

        mols = []
        for m in chimera.openModels.list() :
            if m.display == True and type(m) == chimera.Molecule :
                mols.append ( m )

        print " - %s" % mols[0].name
        print " - %s" % mols[1].name

        m1, m2 = mols

        chainId = self.chain.get()

        rids = {}
        for r in m1.residues :
            if r.id.chainId == chainId :
                rids[r.id.position] = r


        for r in m2.residues :
            if not r.id.position in rids :
                #print " - %d %s %s" % (r.id.position, r.type, r.id.chainId)
                chimera.selection.addCurrent ( r )


        return


        aMap = dict()
        for ri, r in enumerate ( m2.residues ) :
            nres = m1.newResidue (r.type, chimera.MolResId(chainId, r.id.position))
            for at in r.atoms :
                nat = m1.newAtom (at.name, chimera.Element(at.element.number))
                aMap[at] = nat
                nres.addAtom( nat )
                p = chimera.Point ( at.coord().x, at.coord().y, at.coord().z )
                nat.setCoord ( p )

        for bond in m2.bonds :
            try :
                nb = m1.newBond ( aMap[bond.atoms[0]], aMap[bond.atoms[1]] )
                nb.display = nb.Smart
            except :
                pass



    def Domains ( self ) :

        print "dms"

        mol = self.cur_mol
        rmap = {}
        for r in mol.residues :
            rmap[r.id.position] = r
            r.dms = None

        dms = []
        fp = open ( "/Users/greg/Box Sync/20 Ribozyme - Zhaoming/L21RNA_DMS_0000.JustWT.txt" )
        for l in fp :
            #print l,
            s = l.split()
            #print s
            try :
                rid = int(s[0])
            except :
                continue
            if rid in rmap :
                r = rmap[rid]
                r.dms = float(s[1])
                #print "res %d - %g" % (rid, r.dms)
                dms.append ( r.dms )
            else :
                print "res %d - x" % rid

        print "min: %g" % numpy.min ( dms )
        print "max: %g" % numpy.max ( dms )
        print "avg: %g" % numpy.mean ( dms )
        print "std: %g" % numpy.std ( dms )

        dmin = numpy.min ( dms )
        dmax = numpy.max ( dms )
        dmean = numpy.mean ( dms )
        dstd = numpy.std ( dms )

        chimera.selection.clearCurrent ()

        for r in mol.residues :
            #rmap[r.id.position] = r
            #r.dms = None
            if r.dms == None :
                r.ribbonColor = chimera.MaterialColor ( .4, .4, .4, 1.0 )
            else :

                R = numpy.array ( [1,0,0] )
                G = numpy.array ( [0,1,0] )

                f = (r.dms - dmin) / ( dmax - dmin )
                col = R * f + G * (1.0-f)

                if r.dms > 0.05 :
                    col = R
                    chimera.selection.addCurrent ( r )
                    chimera.selection.addCurrent ( rmap[r.id.position-1] )
                    chimera.selection.addCurrent ( rmap[r.id.position+1] )

                else :
                    col = G

                r.ribbonColor = chimera.MaterialColor ( col[0], col[1], col[2], 1.0 )


    def ColorDomains ( self ) :

        print "color domains"

        mol = self.cur_mol

        mfrom = None
        if 0 :
            for m in chimera.openModels.list() :
                if type(m) == chimera.Molecule and m.display == True and m != mol :
                    mfrom = m
                    break


        if mfrom == None :
            rmap = {}
            for r in mol.residues :
                if r.id.position in rmap :
                    rmap[r.id.position].append ( r )
                else :
                    rmap[r.id.position] = [r]
                r.dms = None
                r.ribbonColor = chimera.MaterialColor ( .8, .8, .8, 1.0 )

            dms = []
            #fp = open ( "/Users/greg/GDriveS/_data/Ribozyme2/sec.txt" )
            print mol.openedAs
            mpath = os.path.split ( mol.openedAs[0] )[0] + "/sec.txt"
            print mpath
            fp = open ( mpath )
            for l in fp :
                #print l,
                s = l.split()
                #print s

                C = s[1].split(",")
                C = ( float(C[0]), float(C[1]), float(C[2]) )

                print s[0], C,

                r = s[2].split(",")
                for rs in r :
                    print rs,
                    be = rs.split("-")
                    for i in range ( int(be[0]), int(be[1])+1 ) :
                        if not i in rmap :
                            print " - res %d not in rmap" % i
                        else :
                            for r in rmap[i] :
                                r.ribbonColor = chimera.MaterialColor ( C[0], C[1], C[2], 1.0 )

                print ""

        else :

            print "from:", mfrom.name
            print "to:", mol.name

            rmap = {}
            for r in mfrom.residues :
                rmap[r.id.chainId + "%d"%r.id.position] = r

            for r in mol.residues :
                rid = r.id.chainId + "%d"%r.id.position
                if not rid in rmap :
                    print "r %s not found in %s" % (rid, mol.name)
                    continue
                rf = rmap[rid]
                for at in r.atoms :
                    if at.name in rf.atomsMap :
                        atf = rf.atomsMap[at.name][0]
                        at.setCoord ( atf.coord() )
                    else :
                        print " - at %s in %s.%d.%s not found in %s.%d.%s" % (at.name, r.type, r.id.position, r.id.chainId, rf.type, rf.id.position, rf.id.chainId)
                        break



    def AddH ( self ) :

        print "addh"

        selAt = chimera.selection.currentAtoms()[0]

        print selAt.name

        aN = selAt
        aC1 = selAt.residue.atomsMap["CE1"][0]
        aC2 = selAt.residue.atomsMap["CD2"][0]

        v1 = aC1.coord() - aN.coord(); v1.normalize()
        v2 = aC2.coord() - aN.coord(); v2.normalize()

        avgV = v1 + v2
        avgV.normalize()

        nat = selAt.molecule.newAtom ( "HNE2", chimera.Element(1))
        selAt.residue.addAtom( nat )
        nat.drawMode = nat.EndCap
        nat.setCoord ( aN.coord() - avgV * 1.0 )
        nat.display = True
        if nat.element.name.upper() in atomColors : nat.color = atomColors[nat.element.name.upper()]

        nb = selAt.molecule.newBond ( aN, nat )
        nb.display = nb.Smart
        nb.drawMode = nb.Stick








    def AddDiS ( self ) :

        print ""
        print "AddDiS"

        selAts = chimera.selection.currentAtoms()

        if len(selAts) != 2 :
            umsg ( "Select two atoms" )
            return

        at1, at2 = selAts
        if at1.name != "OG" or at1.residue.type != "SER" :
            umsg ( "Check atoms" )
            return

        if at2.name != "OG" or at2.residue.type != "SER" :
            umsg ( "Check atoms" )
            return

        mol1 = at1.molecule
        mol2 = at2.molecule

        if mol1 != mol2 :
            umsg ( "Not same molecule" )
            return

        nb = mol1.newBond ( at1, at2 )
        nb.display = nb.Smart
        nb.drawMode = nb.Stick




    def AddResProt ( self ) :

        #startI = self.seqRes [ max(self.seqSel[0],0) ].id.position
        #endI = self.seqRes [ min(self.seqSel[1],len(self.seqRes)-1) ].id.position

        print ""
        print "AddRes"

        mol = self.cur_mol
        if self.cur_mol == None :
            umsg ("Select a molecule first")
            return []

        chainId = self.chain.get()

        dmap = self.cur_dmap
        print " - in map: %s" % dmap.name

        startRi = self.seqRes[0].id.position + self.seqSel[0]
        endRi = self.seqRes[0].id.position + self.seqSel[1]
        numRes = endRi - startRi + 1

        print " - sel %d - %d" % (self.seqSel[0], self.seqSel[1])
        #print "res %d - %d" % (startI, endI)
        print " - res %d - %d, %d res" % (startRi, endRi, numRes)

        seq = self.addRess.get().upper().strip().replace(" ", "")
        print " - seq:", seq

        from chimera.resCode import protein1to3
        #from chimera.resCode import nucleic1t3
        for i in range ( len(seq) ) :
            if not seq[i] in protein1to3 :
                umsg ( "Sequence position %d '%s' not known" % (i+1, seq[i]) )
                return

        if len(seq) != numRes :
            umsg ( "%s is %d, need %d" % (seq, len(seq), numRes) )

        molbuild.BuildModLoop ( mol, startRi, endRi, seq, chainId )



    def Occ ( self ) :

        mol = self.cur_mol
        if self.cur_mol == None :
            umsg ("Select a molecule first")
            return []

        #chainId = self.chain.get()

        for at in mol.atoms :
            #if at.residue.id.chainId == chainId :
            ats = at.residue.atomsMap[at.name]

            if len(ats) > 1 :
                alts = {}
                for at in ats :
                    alts[at.altLoc] = at

                locs = alts.keys()
                locs.sort()
                #print at.name, at.residue.type, at.residue.id.position, locs

                sum = 0.0
                occ = 1.0 / float(len(locs))
                occ = round(occ * 100.0)/100.0
                for l in locs[:-1] :
                    alts[l].occupancy = occ
                    sum += occ
                alts[locs[-1]].occupancy = 1.0 - sum

            else :
                ats[0].occupancy = 1.0




    def RibD ( self ) :

        print self.cur_mol.name
        print self.cur_dmap.name

        minR, maxR = 1e9, -1e9

        L, H = numpy.array([1.0,0,0]), numpy.array([0,1.0,0])
        l, h = 2.0, 7.0

        for r in self.cur_mol.residues :
            points = _multiscale.get_atom_coordinates ( r.atoms, transformed = False )
            dvals = self.cur_dmap.interpolated_values ( points, self.cur_mol.openState.xform )
            dval = numpy.average(dvals)
            minR = min ( minR, dval )
            maxR = max ( maxR, dval )

            f = (dval - l) / (h-l)
            C = f * L + (1-f) * H
            r.ribbonColor = chimera.MaterialColor ( C[0], C[1], C[2], 1.0 )
            for at in r.atoms :
                at.color = r.ribbonColor

        print minR, maxR



    def AfColor ( self ) :

        import axes; reload (axes)
        newMod = _surface.SurfaceModel()
        newMod = axes.AddArrow4 ( chimera.Point(0,0,0), chimera.Vector(0,0,1), 4.0, (1,.2,.2,1), 0.2, newMod, 0.4, 0.4 )
        # ( pos, v, d, clr=(0,1,1,1), rad=0.2, mol=None, hrad=3.0, hlen=3.0 )
        chimera.openModels.add ( [newMod] )

        print self.cur_mol.name
        #print self.cur_dmap.name
        numAbove, num = 0.0, 0.0
        for r in self.cur_mol.residues :
            if 'CA' in r.atomsMap :
                bf = r.atomsMap['CA'][0].bfactor
                num += 1.0
                if bf > 80 :
                    r.ribbonColor = chimera.MaterialColor ( .3, .3, 1.0, 1.0 )
                    numAbove += 1.0
                else :
                    r.ribbonColor = chimera.MaterialColor ( 1.0, .5, .5, 1.0 )

        print " %.0f / %.0f pLDDT above 80 -> %.2f %%" % (numAbove, num, numAbove/num*100.0)


    def TunnelVis0 ( self ) :

        for m in chimera.openModels.list() :
            if m.display == True :
                print m.name, len(m.atoms)
                for at in m.atoms :
                    at.radius = at.bfactor
                    at.color = chimera.MaterialColor ( 0.8, 0.8, .2, 0.2 )

    def TunnelVis ( self ) :

        minR, maxR = 1e9, 0.0
        for m in chimera.openModels.list() :
            if m.display == True and type(m) == chimera.Molecule :
                print m.name, len(m.atoms)
                for at in m.atoms :
                    minR = min (at.bfactor, minR)
                    maxR = max (at.bfactor, maxR)

        print "minR %.2f maxR %.2f" % (minR, maxR)

        C1 = numpy.array ( [0.9, 0.2, 0.2] )
        C2 = numpy.array ( [0.9, 0.9, 0.2] )

        for m in chimera.openModels.list() :
            if m.display == True and type(m) == chimera.Molecule :
                print m.name, len(m.atoms)
                for at in m.atoms :
                    at.radius = at.bfactor
                    f = 1.0 - (at.bfactor - minR) / (maxR - minR)
                    C = C1 * f + C2 * (1.0-f)
                    at.color = chimera.MaterialColor ( C[0], C[1], C[2], 1.0 )


    def CompareChains ( self ) :

        print "-----"
        m1, c1 = self.cur_mol, self.chain.get ()
        print m1.name, c1
        m2, c2 = chimera.selection.currentMolecules()[0], chimera.selection.currentResidues()[0].id.chainId
        print m2.name, c2

        print "-----"

        rmap = {}
        for r1 in m1.residues :
            if r1.id.chainId == c1 :
                rmap[r1.id.position] = r1

        for r2 in m2.residues :
            if r2.id.chainId == c2 :
                if r2.id.position not in rmap :
                    print " - residues %d.%s not in chain %s" % (r2.id.position, c2, c1)
                    continue
                r1 = rmap[r2.id.position]
                if r1.type != r2.type :
                    print " - residue %d.%s not same as %d.%s" % (r2.id.position, c2, r1.id.position, c1)
                    continue
                for a2 in r2.atoms :
                    if a2.name not in r1.atomsMap :
                        print " - atom %s in res %d.%s not in res %d.%s" % (a2.name, r2.id.position, r2.id.chainId, r1.id.position, r1.id.chainId)
                        #break




    def RMSD ( self ) :

        mols = []
        for m in chimera.openModels.list() :
            if m.display == True and type(m) == chimera.Molecule :
                mols.append ( m )

        if len(mols) != 2 :
            umsg ( "Make at least two molecules visible" )
            return

        m1, m2 = mols

        qscores.SetBBAts ( m1 )
        qscores.SetBBAts ( m2 )

        print "\nRMSD"
        print "%s : %s" % (m1.name, m2.name)

        atids = {}
        rmap = {}
        for r in m1.residues :
            if r.isProt :
                rmap["%d.%s"%(r.id.position, r.id.chainId)] = r
                for at in r.atoms :
                    if len(r.atomsMap[at.name]) > 1 :
                        # ignore alt conformations...
                        continue
                    else :
                        #atId = "%d.%s.%s.%s" % (r.id.position,r.id.chainId,at.name,at.altLoc)
                        atId = "%d.%s.%s" % (r.id.position,r.id.chainId,at.name)
                        atids[atId] = at


        sums, N = {"All":0.0, "BB":0.0, "SC":0.0}, {"All":0.0, "BB":0.0, "SC":0.0}

        for r2 in m2.residues :

            if r2.isProt :

                rId = "%d.%s" % (r2.id.position, r2.id.chainId)
                if rId not in rmap :
                    #print " - res %s not in m1" % rId
                    continue

                r1 = rmap[rId]

                for at2 in r2.atoms :

                    if len(r2.atomsMap[at2.name]) > 1 :
                        # ignore alt conformations...
                        continue

                    if at2.element.name == "H" :
                        continue

                    at2Id = "%d.%s.%s" % (r2.id.position,r2.id.chainId,at2.name)

                    if at2Id not in atids :
                        #print " - atom %s not in m1" % (at2Id)
                        continue

                    at1 = atids[at2Id]

                    #atPos = m2.openState.xform.inverse().apply ( at1.xformCoord() )
                    v = at2.xformCoord() - at1.xformCoord()

                    v2 = v.length * v.length

                    sums["All"] += v2; N["All"] += 1.0
                    if at2.isBB :
                        sums["BB"] += v2; N["BB"] += 1.0
                    else :
                        sums["SC"] += v2; N["SC"] += 1.0


        for k in sums.keys() :
            #print "%s\t%.3f\t%.3f" % (k, sums[k], N[k] )
            print "%s\t%.3f" % (k, numpy.sqrt(sums[k] / N[k]) )
        print ""



    def Rotas ( self, res ) :


        ctrRes = chimera.Vector(0,0,0)
        for at in res.atoms :
            ctrRes += at.coord().toVector()

        ctrRes = ctrRes / float ( len(res.atoms) )
        ctrRes = chimera.Point ( ctrRes[0], ctrRes[1], ctrRes[2] )
        #print " - res %s %d.%s ctr " % ( res.type, res.id.position, res.id.chainId ), ctrRes

        #print " - in %s" % self.cur_dmap.name

        treeAts, treeAtsAll = [], []
        #print " - %d atoms in %s" % ( len(res.molecule.atoms), res.molecule.name )
        for at in res.molecule.atoms :
            d = (at.coord() - ctrRes).length
            if d < 40.0 :
                treeAtsAll.append ( at )
                if at.residue != res :
                    treeAts.append ( at )

        #print " - %d atoms within 40 - %d all" % ( len(treeAts), len(treeAtsAll) )

        points = _multiscale.get_atom_coordinates ( treeAts, transformed = False )
        atTree = AdaptiveTree ( points.tolist(), treeAts, 2.0)

        points = _multiscale.get_atom_coordinates ( treeAtsAll, transformed = False )
        atTreeAll = AdaptiveTree ( points.tolist(), treeAtsAll, 2.0)



        #print rmols

        rotas = []
        bbdep, rmols = getRotamers ( res, log=False )

        for ri, rmol in enumerate ( rmols ) :

            rotres = rmol.residues[0]
            rotres.rotamerProb = rmol.rotamerProb

            #print ri, rmol.rotamerProb

            to_ats = [ res.atomsMap['N'][0],res.atomsMap['CA'][0],res.atomsMap['CB'][0] ]
            rot_ats = [ rotres.atomsMap['N'][0],rotres.atomsMap['CA'][0],rotres.atomsMap['CB'][0] ]
            xf, rmsd = chimera.match.matchAtoms ( to_ats, rot_ats )


            clash = False
            rotres.clashes = False
            rotAts = []
            rotPos = []

            for ai, rat in enumerate ( rotres.atoms ) :

                atPos = xf.apply(rat.coord())
                rat.setCoord ( atPos )

                if rat.name != "C" and rat.name != "N" and rat.name != "CA" and rat.name != "O" :

                    rotAts.append ( rat )
                    rotPos.append ( atPos )
                    nearAts = self.AtsWithinPt ( atPos.data(), 2.0, atTree )
                    if len(nearAts) > 0 :
                        rotres.clashes = True
                        #for d, a in nearAts :
                        #    print " - at %s - %.2f - at %s in %d.%s" % (rat.name, d, a.name, a.residue.id.position, a.residue.id.chainId )
                        #    break
                        break


            if rotres.clashes :
                continue

            #dvals = dmap.interpolated_values ( apos, r.molecule.openState.xform )

            if self.cur_dmap == None :
                umsg ( "No map selected" )
                return

            #rotres.CC, ccm = ccAts ( rotAts, self.cur_dmap, resolution=3.0, mol=res.molecule )
            #rotres.AvgD = avgdAts ( rotAts, self.cur_dmap, mol=res.molecule )

            mol = res.molecule
            dmap = self.cur_dmap

            molg = MyMolMapX2 ( rotAts, 3.0, dmap.data.step[0], chimera.Xform.identity() )
            fpoints, fpoint_weights = fit_points_g ( molg, 1e-2 )
            map_values = dmap.interpolated_values ( fpoints, mol.openState.xform )
            #print map_values
            olap, rotres.CC, bbCCm = FitMap.overlap_and_correlation ( fpoint_weights, map_values )


            dvals = self.cur_dmap.interpolated_values ( rotPos, res.molecule.openState.xform )
            #print dvals
            rotres.AvgD = numpy.average(dvals)

            avgQ = 0
            minD, maxD = qscores.MinMaxD ( dmap )
            #print "%d | " % ri,
            for at in rotAts :
                Qs = qscores.Qscore ( [at], dmap, 0.6, allAtTree=atTreeAll, show=0, log=0, numPts=8, toRAD=2.0, dRAD=0.1, minD=minD, maxD=maxD, mol=mol )
                #print "%s:%.3f " % (at.name, ccm),
                avgQ += Qs


            rotres.Q = avgQ / float(len(rotAts))
            #print " | %.3f" % rotres.Q



            rotas.append ( rotres )
            #break

        return rotas



    def Rotas_ ( self, res ) :

        rotas = []

        bbdep, rmols = getRotamers ( res, log=False )

        #print rmols

        for ri, rmol in enumerate ( rmols ) :

            rotres = rmol.residues[0]
            rotres.rotamerProb = rmol.rotamerProb

            #print ri, rmol.rotamerProb

            to_ats = [ res.atomsMap['N'][0],res.atomsMap['CA'][0],res.atomsMap['CB'][0] ]
            rot_ats = [ rotres.atomsMap['N'][0],rotres.atomsMap['CA'][0],rotres.atomsMap['CB'][0] ]
            xf, rmsd = chimera.match.matchAtoms ( to_ats, rot_ats )

            for ai, rat in enumerate ( rotres.atoms ) :
                #if rat.name == "C" or rat.name == "N" or rat.name == "CA" or rat.name == "O" :
                #    continue
                rat.setCoord ( xf.apply(rat.coord()) )

            rotas.append ( rotres )

        return rotas


    def ApplyRota ( self, res, rota ) :
        for at in res.atoms :
            if at.name == "C" or at.name == "N" or at.name == "CA" or at.name == "O" :
                continue

            rotaAt = rota.atomsMap[at.name][0]
            at.setCoord ( rotaAt.coord() )


    def HohRota ( self ) :

        res = chimera.selection.currentResidues()[0]
        print "Res %d.%s %s" % (res.id.position, res.id.chainId, res.type)

        rotas = self.Rotas ( res )

        print "#\tProb\tCC\tAvg.D.\tQ"


        #rotas.sort ( reverse=True, key=lambda r: r.CC )
        rotas.sort ( reverse=True, key=lambda r: r.Q )


        for ri, r in enumerate ( rotas ) :
            #print " - %d, prob %.5f, cc " % (ri, r.rotamerProb),
            print "%d\t%f\t%f\t%f\t%f" % (ri+1, r.rotamerProb, r.CC, r.AvgD, r.Q),

            if r.clashes :
                print "--x--"
            else :
                print ""


        if len(rotas) > 0 :
            #ri = int ( numpy.floor ( ( random.random() * len(rotas) ) ) )
            print " - applying %d/%d" % (1, len(rotas))
            self.ApplyRota ( res, rotas[0] )

        self.rotas = rotas
        self.rotaAt = 0
        self.rotaRes = res


    def HohRotaL ( self ) :

        if not hasattr ( self, 'rotas' ) :
            return

        self.rotaAt = max ( self.rotaAt - 1, 0 )
        rota = self.rotas[self.rotaAt]
        print " - applying rota %d/%d - prob %f, cc %f" % (self.rotaAt+1, len(self.rotas), rota.rotamerProb, rota.CC)
        self.ApplyRota ( self.rotaRes, rota )


    def HohRotaR ( self ) :

        if not hasattr ( self, 'rotas' ) :
            return

        self.rotaAt = min ( self.rotaAt+1, len(self.rotas)-1 )
        rota = self.rotas[self.rotaAt]
        print " - applying rota %d/%d - prob %f, cc %f" % (self.rotaAt+1, len(self.rotas), rota.rotamerProb, rota.CC)
        self.ApplyRota ( self.rotaRes, rota )






    def ResMap ( self ) :


        print " - resmap - "

        mol = self.cur_mol
        if self.cur_mol == None :
            umsg ("Select a molecule first")
            return []

        chainId = self.chain.get()

        dmap = self.cur_dmap
        print " - scale map: %s" % dmap.name

        rmap = None
        for m in chimera.openModels.list() :
            if "resmap" in m.name :
                rmap = m

        print "mol:", self.cur_mol.name
        print "resmap:", rmap.name


        #points = _multiscale.get_atom_coordinates ( mol.atoms, transformed = False )

        molPath = os.path.splitext(mol.openedAs[0])[0]
        mapName = os.path.splitext(rmap.name)[0]

        nname = molPath + "__R__" + mapName + ".txt"
        print " - q vs resmap:", nname

        fp = open ( nname, "w" )


        for at in mol.atoms :
            res = rmap.interpolated_values ( [at.coord().data()], mol.openState.xform )
            fp.write ( "%f\t%f\n" % (at.Q, res)  )


        fp.close()
        print " - done"






    def BB_Sigma (self) :

        selAts = chimera.selection.currentAtoms()
        if len ( selAts ) == 0 :
            return

        dmap = self.cur_dmap


        a = selAts[0]
        r = a.residue
        print "Res: %s - %d.%s - %s - Atom: %s" % (r.type, r.id.position, r.id.chainId, r.molecule.name, a.name)

        if 1 or not hasattr ( r.molecule, 'bbats' ) :
            qscores.SetBBAts(r.molecule)
            r.molecule.bbats = True

        removeMods = []
        for m in chimera.openModels.list() :
            if "RAD points" in m.name :
                removeMods.append ( m )
        chimera.openModels.remove ( removeMods )


        ats = [at for at in self.cur_mol.atoms if not at.element.name == "H"]
        points = _multiscale.get_atom_coordinates ( ats, transformed = False )
        print " - search tree: %d/%d ats" % ( len(ats), len(r.molecule.atoms) )
        allAtTree = AdaptiveTree ( points.tolist(), ats, 1.0)

        #allAtTree = None
        #print "-"

        import time
        start = time.time()

        sigma = RadAts ( r.bbAtoms, dmap, allAtTree=allAtTree, show=0, log=1, numPts=10, toRAD=2, dRAD=0.25 )

        end = time.time()

        print "%s - rad: %.3f, time: %f" % ( a.name, sigma, (end - start) )




    def ZScoreSel (self) :

        selAts = chimera.selection.currentAtoms()
        if len ( selAts ) == 0 :
            return

        dmap = self.cur_dmap


        a = selAts[0]
        r = a.residue
        print "Res: %s - %d.%s - %s - Atom: %s" % (r.type, r.id.position, r.id.chainId, r.molecule.name, a.name)

        if not hasattr ( r.molecule, 'bbats' ) :
            qscores.SetBBAts(r.molecule)
            r.molecule.bbats = True

        removeMods = []
        for m in chimera.openModels.list() :
            if "RAD points" in m.name :
                removeMods.append ( m )
            if "SC " in m.name :
                removeMods.append ( m )
        chimera.openModels.remove ( removeMods )


        scZ, cc = zRotSideChain ( r.molecule, r, 3.0, dmap, show=True )
        print "- scZ %.3f, cc %.3f" % (scZ, cc)
        #print "%f\t%f\t%f" % (r.sigma,scZ,cc)




    def RotaZ1 (self) :

        selAts = chimera.selection.currentAtoms()
        if len ( selAts ) == 0 :
            return

        dmap = self.cur_dmap


        a = selAts[0]
        r = a.residue
        print "Res: %s - %d.%s - %s - Atom: %s" % (r.type, r.id.position, r.id.chainId, r.molecule.name, a.name)

        if not hasattr ( r.molecule, 'bbats' ) :
            qscores.SetBBAts(r.molecule)
            r.molecule.bbats = True

        removeMods = []
        for m in chimera.openModels.list() :
            if "RAD points" in m.name :
                removeMods.append ( m )
            if "SC " in m.name :
                removeMods.append ( m )
        chimera.openModels.remove ( removeMods )


        ats = [at for at in self.cur_mol.atoms if not at.element.name == "H"]
        points = _multiscale.get_atom_coordinates ( ats, transformed = False )
        print " - search tree: %d/%d ats" % ( len(ats), len(r.molecule.atoms) )
        allAtTree = AdaptiveTree ( points.tolist(), ats, 1.0)


        rZ = RadZ ( r.scAtoms, dmap, allAtTree=allAtTree, show=0, log=1, numPts=10, toRAD=2 )


        #scZ, cc = zRotSideChain ( r.molecule, r, 3.0, dmap, show=True )
        print "- radZ %.3f " % (rZ)
        #print "%f\t%f\t%f" % (r.sigma,scZ,cc)




    def R1 (self) :

        selAts = chimera.selection.currentAtoms()
        if len ( selAts ) == 0 :
            return

        dmap = self.cur_dmap


        a = selAts[0]
        r = a.residue
        print "Res: %s - %d.%s - %s - Atom: %s" % (r.type, r.id.position, r.id.chainId, r.molecule.name, a.name)

        if not hasattr ( r.molecule, 'bbats' ) :
            qscores.SetBBAts(r.molecule)
            r.molecule.bbats = True

        removeMods = []
        for m in chimera.openModels.list() :
            if "RAD points" in m.name :
                removeMods.append ( m )
            if "SC " in m.name :
                removeMods.append ( m )
        chimera.openModels.remove ( removeMods )



        ress = []
        bbAtoms = []
        allAtoms = []
        for r in a.molecule.residues :
            if r.id.chainId == a.residue.id.chainId :
                ress.append ( r )
                bbAtoms.extend ( r.bbAtoms )
                allAtoms.extend ( r.atoms )

        avgD = avgdAts ( allAtoms, dmap )
        bbAvgD = avgdAts ( bbAtoms, dmap )
        print " - avgd - all: %f, bb: %f" % (avgD, bbAvgD)

        r = a.residue
        if len(r.scAtoms) > 0 :
            scAvgD = avgdAts ( r.scAtoms, dmap )
            r.SCBBr = scAvgD / bbAvgD
            print " - residue %s.%d, %d side chain atoms, avgd: %.5f, r: %.5f" % ( r.type, r.id.position, len(r.scAtoms), scAvgD, r.SCBBr/bbAvgD )
        else :
            r.SCBBr = None
            print " - residue %s.%d - no side chain atoms" % ( r.type, r.id.position )







    def AlignRes1 ( self ) :

        if self.cur_mol == None :
            umsg ("Select a molecule first")
            return

        chainId = self.chain.get()
        if len(chainId) == 0 :
            umsg ("Select a chain first")
            return

        if self.cur_dmap == None :
            umsg ("Select a map first")
            return


        #SetBBAts ( self.cur_mol )
        last_x = 0.0
        last_y = 0.0

        r0, exR0, xtR0 = None, None, None

        alAts = []
        if self.exType == "ASP" : alAts = ["CG","OD1","OD2"]
        if self.exType == "LEU" : alAts = ["CG","CD1","CD2"]
        if self.exType == "GLU" : alAts = ["CD","OE1","OE2"]
        if self.exType == "TYR" : alAts = ["OH","CE1","CE2","CD1","CD2","CG","CB"]


        for r in self.cur_mol.residues :
            if r.id.chainId == chainId and r.type == self.exType :
                print " - res %s %d" % (r.type, r.id.position)

                if r0 == None :
                    r0 = r

                    r.exMaps[0].display = True
                    r.exMaps[1].display = False

                    #r.xtMaps[0].display = False
                    #r.xtMaps[1].display = False

                    for at in r.atoms :
                        at.display = at.name in alAts

                else :

                    exR0 = r0.exMol.residues[0]
                    exR = r.exMol.residues[0]
                    ats0, ats = [], []
                    for atName in alAts :
                        ats0.append ( exR0.atomsMap[atName][0] )
                        ats.append ( exR.atomsMap[atName][0] )

                    for at in r.atoms :
                        at.display = at.name in alAts

                    #aCG0, aOD10, aOD20 = exR0.atomsMap['CG'][0], exR0.atomsMap['OD1'][0], exR0.atomsMap['OD2'][0],
                    #aCG, aOD1, aOD2 = exR.atomsMap['CG'][0], exR.atomsMap['OD1'][0], exR.atomsMap['OD2'][0],

                    #xf, rmsd = chimera.match.matchPositions ( pts_o, pts_c )
                    #xf, rmsd = chimera.match.matchAtoms ( [aCG0, aOD10, aOD20], [aCG, aOD1, aOD2] )
                    xf, rmsd = chimera.match.matchAtoms ( ats0, ats )
                    print " - rmsd: ", rmsd

                    #from _multiscale import get_atom_coordinates
                    #points = get_atom_coordinates ( atoms, transformed = True )

                    #exR.xf0 = r.exMol.openState.xform

                    mxf = r0.exMol.openState.xform
                    mxf.multiply ( xf )
                    r.exMol.openState.xform = mxf
                    r.exMaps[0].openState.xform = mxf
                    r.exMaps[1].openState.xform = mxf
                    r.exMaps[0].display = True
                    r.exMaps[1].display = False

                    #r.xtMaps[0].display = False
                    #r.xtMaps[1].display = False


                    #break


    def AlignRes2 ( self ) :

        if self.cur_mol == None :
            umsg ("Select a molecule first")
            return

        chainId = self.chain.get()
        if len(chainId) == 0 :
            umsg ("Select a chain first")
            return

        if self.cur_dmap == None :
            umsg ("Select a map first")
            return


        #SetBBAts ( self.cur_mol )
        last_x = 0.0
        last_y = 0.0

        r0, exR0, xtR0 = None, None, None


        for r in self.cur_mol.residues :
            if r.id.chainId == chainId and r.type == "ASP" :
                print " - res %s %d" % (r.type, r.id.position)

                if r0 == None :
                    r0 = r

                    r.exMaps[0].display = False
                    r.exMaps[1].display = False

                    r.xtMaps[0].display = True
                    r.xtMaps[1].display = False


                else :

                    r.exMaps[0].display = False
                    r.exMaps[1].display = False

                    exR0 = r0.xtMol.residues[0]
                    aCB0, aCG0, aOD10, aOD20 = exR0.atomsMap['CB'][0], exR0.atomsMap['CG'][0], exR0.atomsMap['OD1'][0], exR0.atomsMap['OD2'][0],

                    exR = r.xtMol.residues[0]
                    aCB, aCG, aOD1, aOD2 = exR.atomsMap['CB'][0], exR.atomsMap['CG'][0], exR.atomsMap['OD1'][0], exR.atomsMap['OD2'][0],

                    #xf, rmsd = chimera.match.matchPositions ( pts_o, pts_c )
                    xf, rmsd = chimera.match.matchAtoms ( [aCB0, aCG0, aOD10, aOD20], [aCB, aCG, aOD1, aOD2] )
                    print " - rmsd: ", rmsd

                    #from _multiscale import get_atom_coordinates
                    #points = get_atom_coordinates ( atoms, transformed = True )

                    #exR.xf0 = r.exMol.openState.xform

                    mxf = r0.xtMol.openState.xform
                    mxf.multiply ( xf )
                    r.xtMol.openState.xform = mxf
                    r.xtMaps[0].openState.xform = mxf
                    r.xtMaps[1].openState.xform = mxf
                    r.xtMaps[0].display = True
                    r.xtMaps[1].display = False


                    #break



    def Avg ( self ) :

        print " -- finding base map --- "
        largestMap = None
        maxD = 0
        for m in OML(modelTypes = [VolumeViewer.volume.Volume]) :
            if m.display == True :
                d = numpy.sum ( m.data.size )
                if d > maxD :
                    maxD = d
                    largestMap = m

        print " - largest map: ", largestMap.name
        dmap = largestMap
        dmap.display = False


        fmap = None
        avgMat = dmap.data.full_matrix()
        N = 0.0

        print " ----------- Averaging... ---------------------"

        for m in OML(modelTypes = [VolumeViewer.volume.Volume]) :
            if m.display == True and m != dmap :
                print m.name

                df_mat = self.Map2Map ( m, dmap )
                m.display = False
                N = N + 1.0
                avgMat = avgMat + df_mat


        print " ----------- n=%f ---------------------" % N

        avgMat = avgMat / N
        df_data = VolumeData.Array_Grid_Data ( avgMat, dmap.data.origin, dmap.data.step, dmap.data.cell_angles, name="avg" )

        MapFromData ( df_data, "Avg", dmap, False )
        MapFromData ( df_data, "Avg", dmap, True )


        #df_v = VolumeViewer.volume.volume_from_grid_data ( df_data )
        #df_v.name = "Avg"
        #df_v.openState.xform = dmap.openState.xform

        #nv = self.ShrinkMap ( df_v, 1e-3 )



    def Map2Map ( self, densitiesFromMap, toGridOfMap, mask = False ) :

        fmap = toGridOfMap
        dmap = densitiesFromMap

        import _contour
        n1, n2, n3 = fmap.data.size[0], fmap.data.size[1], fmap.data.size[2]
        f_points = VolumeData.grid_indices( (n1,n2,n3), numpy.single )  # i,j,k indices
        _contour.affine_transform_vertices( f_points, fmap.data.ijk_to_xyz_transform )

        d_vals = dmap.interpolated_values ( f_points, fmap.openState.xform )
        df_mat = d_vals.reshape( (n3,n2,n1) )

        if mask :
            f_mat = fmap.data.full_matrix()
            f_mask = numpy.where ( f_mat > fmap.surface_levels[0], numpy.ones_like(f_mat), numpy.zeros_like(f_mat) )
            df_mat = df_mat * f_mask

        return df_mat




    def CloseExtracted ( self ) :

        if self.cur_mol == None :
            umsg ("Select a molecule first")
            return

        chainId = self.chain.get()
        if len(chainId) == 0 :
            umsg ("Select a chain first")
            return

        if self.cur_dmap == None :
            umsg ("Select a map first")
            return


        for r in self.cur_mol.residues :

            if hasattr ( r, "exMaps" ) :
                chimera.openModels.close ( r.exMaps ); del r.exMaps

            if hasattr ( r, "xtMaps" ) :
                chimera.openModels.close ( r.xtMaps ); del r.xtMaps

            if hasattr ( r, "exMol" ) :
                chimera.openModels.close ( [r.exMol] ); del r.exMol

            if hasattr ( r, "xtMol" ) :
                chimera.openModels.close ( [r.xtMol] ); del r.xtMol

        for m in chimera.openModels.list() :
            if m.name == "Avg" or m.name == "Avg_mesh" :
                chimera.openModels.close ( [m] )





    def Extract ( self ) :

        if self.cur_mol == None :
            umsg ("Select a molecule first")
            return

        chainId = self.chain.get()
        if len(chainId) == 0 :
            umsg ("Select a chain first")
            return

        if self.cur_dmap == None :
            umsg ("Select a map first")
            return


        #SetBBAts ( self.cur_mol )
        last_x = 0.0
        last_y = 0.0


        print "Extracting - %s - %s - %s" % (self.cur_dmap.name, self.cur_mol.name, chainId)

        #self.exType = "TYR"
        #self.exType = "GLU"
        #self.exType = "ASP"
        self.exType = "LEU"

        yzAts = { "ASP" : ["CB","CG","OD1"],
                  "GLU" : ["CG","CD","OE1"],
                  "TYR" : ["CB","CZ","CD1"],
                  "LEU" : ["CB","CG","CD1"]
        }

        for r in self.cur_mol.residues :

            if r.id.chainId == chainId and r.type == self.exType :

                print " - res %s %d" % (r.type, r.id.position)

                self.ExtractRes ( r, self.cur_mol, self.cur_dmap, last_x, last_y, yzAts[self.exType] )

                #self.ExtendRes ( r, self.cur_mol, self.cur_dmap, last_x, -8.0, thrF=0.8 )

                last_x += 7.0

                #break



    def ExtractRes ( self, r, mol, dmap, atX, atY, xyAts ) :

        nmol, nres = CopyRess ( [r] )
        nmol.name = mol.name + "_%s_%d" % (r.type, r.id.position)
        chimera.openModels.add ( [nmol] )
        nmol.openState.xform = mol.openState.xform

        for at in nmol.atoms :
            #at.drawMode = 3
            if at.element.name.upper() in atomColors : at.color = atomColors[at.element.name.upper()]
            #at.radius = at.radius * 0.8

        mname = dmap.name + "_%s_%d" % (r.type, r.id.position)

        #aCB, aCG, aOD1 = r.atomsMap['CB'][0], r.atomsMap['CG'][0], r.atomsMap['OD1'][0]
        aCB, aCG, aOD1 = r.atomsMap[xyAts[0]][0], r.atomsMap[xyAts[1]][0], r.atomsMap[xyAts[2]][0]

        dmap, mmap = ExtractDen ( r.atoms, dmap, mname, boundRad=2.0, showMesh=True )
        r.exMol = nmol
        r.exMaps = [dmap, mmap]

        X = aOD1.coord() - aCB.coord(); X.normalize()
        Y = aCG.coord() - aCB.coord(); Y.normalize()
        Z = chimera.cross ( X, Y ); Z.normalize()
        X = chimera.cross ( Y, Z ); Y.normalize()

        xf = chimera.Xform.coordFrame ( X, Y, Z, aCB.coord(), True ).inverse()
        xf.premultiply ( chimera.Xform.translation(atX, atY, 0) )

        nmol.openState.xform = xf
        dmap.openState.xform = xf
        if mmap : mmap.openState.xform = xf


    def ExtendRes ( self, r, mol, dmap, atX, atY, thrF=0.75 ) :

        nmol, nres = CopyRess ( [r] )
        nmol.name = mol.name + "_%s_%d_ext" % (r.type, r.id.position)
        chimera.openModels.add ( [nmol] )
        nmol.openState.xform = mol.openState.xform

        for at in nmol.atoms :
            at.drawMode = 3
            if at.element.name.upper() in atomColors : at.color = atomColors[at.element.name.upper()]
            at.radius = at.radius * 0.8

        mname = dmap.name + "_%s_%d_ext" % (r.type, r.id.position)


        R = nres[0]
        R.O, R.N, R.C, R.CA = R.atomsMap["O"][0], R.atomsMap["N"][0], R.atomsMap["C"][0], R.atomsMap["CA"][0]
        R.CB, R.CG, R.OD1, R.OD2 = R.atomsMap["CB"][0], R.atomsMap["CG"][0], R.atomsMap["OD1"][0], R.atomsMap["OD2"][0]

        bones = []
        bones.append ( Bone(R.CA, R.N, R.CB) )
        bones.append ( Bone(R.CA, R.C, R.CB) )
        bones.append ( Bone(R.C, R.O, R.CA) )

        bones.append ( Bone(R.CA, R.CB, R.N) )
        bones.append ( Bone(R.CG, R.CB, R.OD1) )
        bones.append ( Bone(R.CG, R.OD1, R.OD2) )
        bones.append ( Bone(R.CG, R.OD2, R.OD1) )

        for bi, bo in enumerate ( bones ) :
            if GetMod ( "bone_%d.mrc" % bi ) != None : chimera.openModels.close ( "bone_%d.mrc" % bi )
            if GetMod ( "bone_%d.mrc_mesh" % bi ) != None : chimera.openModels.close ( "bone_%d.mrc_mesh" % bi )
            bo.dmap = BoneMap ( bo, dmap, 1.0, "bone_%d.mrc" % bi, show = False, showMesh=True )

        v1 = R.CB.coord() - R.CA.coord(); v1.normalize()
        v2 = R.CB.coord() - R.CG.coord(); v2.normalize()
        ang = numpy.arccos ( v1*v2 ) * 180.0/numpy.pi
        ax = chimera.cross ( v1, v2 ); ax.normalize()

        print "CB-CG: %.2f" % (-ang + 180)

        T = chimera.Xform.translation ( R.CB.coord().toVector() )
        T.multiply ( chimera.Xform.rotation ( ax, -ang + 180 ) )
        T.multiply ( chimera.Xform.translation ( R.CB.coord().toVector()*-1.0 ) )

        for an in ["CG", "OD1", "OD2"] :
            at = R.atomsMap[an][0]
            at.setCoord ( T.apply (at.coord()) )

        #MoldMap2 ( bones, rmaps[0], rmaps[1] )

        d1 = diha ( R.N, R.CB, R.CG, R.OD1 )
        d2 = diha ( R.N, R.CB, R.CG, R.OD2 )
        ang = d1 if numpy.abs(d1) < numpy.abs(d2) else d2
        print "CG dihedral - ", d1, d2, " -> ", ang
        ax = R.CG.coord() - R.CB.coord(); ax.normalize()

        T = chimera.Xform.translation ( R.CG.coord().toVector() )
        T.multiply ( chimera.Xform.rotation ( ax, -ang ) )
        T.multiply ( chimera.Xform.translation ( R.CG.coord().toVector()*-1.0 ) )

        for an in ["OD1", "OD2"] :
            at = R.atomsMap[an][0]
            at.setCoord ( T.apply (at.coord()) )

        dmap, dmesh = MapForAtoms ( R.atoms, dmap, mname, showMesh=True, thrF=thrF )
        MoldMap2 ( bones, dmap, dmesh )
        r.xtMol = nmol
        r.xtMaps = [dmap, dmesh]


        X = R.OD1.coord() - R.CB.coord(); X.normalize()
        Y = R.CG.coord() - R.CB.coord(); Y.normalize()
        Z = chimera.cross ( X, Y ); Z.normalize()
        X = chimera.cross ( Y, Z ); Y.normalize()

        xf = chimera.Xform.coordFrame ( X, Y, Z, R.CB.coord(), True ).inverse()
        xf.premultiply ( chimera.Xform.translation(atX, atY, 0) )

        nmol.openState.xform = xf
        dmap.openState.xform = xf
        if dmesh : dmesh.openState.xform = xf




    def asp ( self ) :

        N = 1

        framei = 0
        mpath = "/Users/greg/Desktop/frames"
        for f in os.listdir ( mpath ) :
            if f.endswith(".png") :
                os.remove( mpath + "/" + f )

        dmap, mol = VisMapMod()
        resolution = 3.0 * dmap.data.step[0]

        print "Map: %s, mol: %s" % (dmap.name, mol.name)
        res = chimera.selection.currentResidues()[0]
        print " - res: %s %d.%s" % (res.type, res.id.position, res.id.chainId)
        z = None

        nname = "%s_%d" % ( res.type, res.id.position )

        #for na in ["ASP","molded.mrc","skinned.mrc"] :
        #    m = GetMod ( na )
        #    if m != None :
        #        chimera.openModels.close ( [m] )


        nmol = GetMod ( nname + ".pdb" )
        if nmol == None :
            nmol, nres = CopyRess ( [res] )
            nmol.name = nname + ".pdb"
            chimera.openModels.add ( [nmol] )
            nmol.openState.xform = mol.openState.xform

            xf = nmol.openState.xform
            #xf.multiply ( chimera.Xform.translation ( 0,0,5 ) )
            nmol.openState.xform = xf

            for at in nmol.atoms:
                at.drawMode = 3
                if at.element.name.upper() in atomColors :
                    at.color = atomColors[at.element.name.upper()]
                    at.radius = at.radius * 0.8

        nres = nmol.residues
        R = nres[0]

        R.O = R.atomsMap["O"][0]
        R.N = R.atomsMap["N"][0]
        R.C = R.atomsMap["C"][0]
        R.CA = R.atomsMap["CA"][0]
        R.CB = R.atomsMap["CB"][0]
        R.CG = R.atomsMap["CG"][0]
        R.OD1 = R.atomsMap["OD1"][0]
        R.OD2 = R.atomsMap["OD2"][0]


        bones = []
        bones.append ( Bone(R.CA, R.N, R.CB) )
        bones.append ( Bone(R.CA, R.C, R.CB) )
        bones.append ( Bone(R.C, R.O, R.CA) )

        bones.append ( Bone(R.CA, R.CB, R.N) )
        bones.append ( Bone(R.CG, R.CB, R.OD1) )
        bones.append ( Bone(R.CG, R.OD1, R.OD2) )
        bones.append ( Bone(R.CG, R.OD2, R.OD1) )

        for bi, bo in enumerate ( bones ) :
            if GetMod ( "bone_%d.mrc" % bi ) != None : chimera.openModels.close ( "bone_%d.mrc" % bi )
            if GetMod ( "bone_%d.mrc_mesh" % bi ) != None : chimera.openModels.close ( "bone_%d.mrc_mesh" % bi )
            bo.dmap = BoneMap ( bo, dmap, 1.0, "bone_%d.mrc" % bi, show = False, showMesh=True )


        v1 = R.CB.coord() - R.CA.coord(); v1.normalize()
        v2 = R.CB.coord() - R.CG.coord(); v2.normalize()
        ang = numpy.arccos ( v1*v2 ) * 180.0/numpy.pi
        print ang
        ax = chimera.cross ( v1, v2 ); ax.normalize()

        dmap.display = False
        mol.display = False

        NB = 2
        #N = 90
        toAng = -ang + 180
        dAng = toAng / float(N)

        print "CB-CG: %.2f/%.2f deg" % (toAng, dAng)

        rmaps = None

        for i in range ( N ) :

            print i,

            T = chimera.Xform.translation ( R.CB.coord().toVector() )
            #T.multiply ( chimera.Xform.rotation ( ax, -ang + 180 ) )
            T.multiply ( chimera.Xform.rotation ( ax, dAng ) )
            T.multiply ( chimera.Xform.translation ( R.CB.coord().toVector()*-1.0 ) )

            for an in ["CG", "OD1", "OD2"] :
                at = R.atomsMap[an][0]
                at.setCoord ( T.apply (at.coord()) )

            #SkinMap ( R.atoms, bones, NB, dmap, 2.0, "skinned.mrc", True)
            #MoldMap ( R.atoms, bones, dmap, "molded.mrc", showMesh=True )

            if rmaps == None :
                rmaps = MapForAtoms ( R.atoms, dmap, nname+".mrc", showMesh=True )
            #    for m in rmaps :
            #        if m != None :
            #            m.openState.xform = nmol.openState.xform

            MoldMap2 ( bones, rmaps[0], rmaps[1] )


            if N > 1 :
                chimera.viewer.postRedisplay()
                self.toplevel_widget.update_idletasks ()
                chimera.printer.saveImage ( mpath + "/%06d.png" % framei )
                framei += 1

        print ""

        if 1 :

            d1 = diha ( R.N, R.CB, R.CG, R.OD1 )
            d2 = diha ( R.N, R.CB, R.CG, R.OD2 )
            ang = d1 if numpy.abs(d1) < numpy.abs(d2) else d2
            print "CG dihedral - ", d1, d2, " -> ", ang
            ax = R.CG.coord() - R.CB.coord(); ax.normalize()

            toAng = -ang
            dAng = toAng / float( max(N/2,1) )
            print "CG dihedral -- %.2f/%.2f deg" % (toAng, dAng)

            for i in range ( max(N/2,1) ) :

                print i,

                T = chimera.Xform.translation ( R.CG.coord().toVector() )
                T.multiply ( chimera.Xform.rotation ( ax, dAng ) )
                T.multiply ( chimera.Xform.translation ( R.CG.coord().toVector()*-1.0 ) )

                for an in ["OD1", "OD2"] :
                    at = R.atomsMap[an][0]
                    at.setCoord ( T.apply (at.coord()) )

                #print "%d bones" % len(bones)
                #PtsToMapSkinD ( R.atoms, bones, NB, dmap, 2.0, "skinned.mrc", True)
                #MoldMap ( R.atoms, bones, dmap, "molded.mrc", showMesh=True )
                MoldMap2 ( bones, rmaps[0], rmaps[1] )

                if N > 1 :
                    chimera.viewer.postRedisplay()
                    self.toplevel_widget.update_idletasks ()
                    chimera.printer.saveImage ( mpath + "/%06d.png" % framei )
                    framei += 1




        if N > 1 :
            args = [ "/Users/greg/_mol/Chimera.app/Contents/Resources/bin/ffmpeg", "-r", "30",
                "-i", mpath + "/%06d.png", "-y", "-qscale", "1", "-b", "9000", "-vcodec", "mpeg4",  # mpeg4 libx264
                "-f", "mov", mpath+"/__ares.mov" ]

            print "- running: "
            for a in args : print a,
            print ""

            import subprocess
            subprocess.call ( args )
            print "done!\n"







def CalcSCBBr ( mol, cid, dmap ) :

    print "Calculating sc-bb ratios..."

    ress = []
    bbAtoms = []
    allAtoms = []
    scAtoms = []
    for r in mol.residues :
        if cid == None or r.id.chainId == cid :
            ress.append ( r )
            bbAtoms.extend ( r.bbAtoms )
            allAtoms.extend ( r.atoms )
            scAtoms.extend ( r.scAtoms )

    bbAvgD, scAvgD = avgdAts ( bbAtoms, dmap ),  avgdAts ( scAtoms, dmap )
    print " - avgd - bb: %.3f, sc: %.3f" % (bbAvgD, scAvgD)

    bbCC, bbCCm = ccAts ( bbAtoms, dmap, 2.0)
    print " - all bb cc: %.3f, ccm: %.3f" % (bbCC, bbCCm)

    cc, ccm = ccAts ( allAtoms, dmap, 2.0)
    print " - all cc: %.3f, ccm: %.3f" % (cc, ccm)


    dr, ccr, ccmr = [], [], []
    for r in ress :
        if len(r.scAtoms) > 0 :
            scAvgD = avgdAts ( r.scAtoms, dmap )
            #rbbAvgD = avgdAts ( r.bbAtoms, dmap )
            r.SCBBr = scAvgD / bbAvgD
            dr.append ( scAvgD / bbAvgD )

            scCC, scCCm = ccAts ( r.scAtoms, dmap, 2.0)
            ccr.append ( scCC/bbCC )
            ccmr.append ( scCCm/bbCCm )

            r.SCBBr = scCCm

        else :
            r.SCBBr = None

    print " - avg-r d:%.3f, cc:%.3f, ccm: %.3f" % ( numpy.average ( dr ), numpy.average ( ccr ), numpy.average ( ccmr ) )

    return cc, ccm, numpy.average ( dr ), numpy.average ( ccr ), numpy.average ( ccmr )


def ccAts ( atoms, dmap, resolution=3.0, mol=None ) :

    if mol == None :
        mol = atoms[0].molecule

    molg = MyMolMapX ( mol, atoms, resolution, dmap.data.step[0], chimera.Xform.identity() )
    fpoints, fpoint_weights = fit_points_g ( molg, 1e-3 )
    map_values = dmap.interpolated_values ( fpoints, mol.openState.xform )
    olap, bbCC, bbCCm = FitMap.overlap_and_correlation ( fpoint_weights, map_values )
    return bbCC, bbCCm


def avgdAts ( atoms, dmap, mol=None ) :

    if mol == None :
        mol = atoms[0].molecule

    if len(atoms) < 1 :
        #print " - no atoms" % len(atoms)
        return 0

    from _multiscale import get_atom_coordinates
    apos = get_atom_coordinates(atoms, transformed = False)
    dvals = dmap.interpolated_values ( apos, mol.openState.xform )
    #print dvals
    return numpy.average(dvals)





def RadZ ( atoms, dmap, allAtTree = None, show=0, log=0, numPts=10, toRAD=2.0 ) :

    if len(atoms) == 0 :
        #print " - no RAD atoms?"
        return None

    #pts = []
    #for at in atoms :
    #    p = at.coord()
    #    pts.append ( [p[0], p[1], p[2]] )

    from _multiscale import get_atom_coordinates
    pts = get_atom_coordinates(atoms, transformed = False)

    d_vals = dmap.interpolated_values ( pts, atoms[0].molecule.openState.xform )
    avg0 = numpy.average ( d_vals )


    #dRAD, toRAD, RAD = 0.2, 1.8, 0.1
    RAD = toRAD
    zscore = None

    outRad = RAD*0.9
    #outRad2 = outRad * outRad
    pts = []
    for at in atoms :
        npts = (numPts * RAD*RAD / (dRAD*dRAD)) if show else numPts
        npts = int ( npts )
        #print RAD, dRAD, numPts, " -> ", npts
        outPts = SpherePts ( at.coord(), RAD, npts )
        for pt in outPts :
            if allAtTree != None :
                vPt = numpy.array ( [pt[0], pt[1], pt[2]] )
                opointsNear = allAtTree.searchTree ( [pt[0], pt[1], pt[2]], outRad )
                if len(opointsNear) > 0 :
                    if 0 :
                        clash = False
                        for p in opointsNear :
                            v = vPt - p.coord().data()
                            sqSum = numpy.sum ( v * v )
                            if sqSum < outRad2 :
                                clash = True
                                break
                        if clash == False :
                            pts.append ( [pt[0], pt[1], pt[2]] )
                else :
                    pts.append ( [pt[0], pt[1], pt[2]] )
            else :
                pts.append ( [pt[0], pt[1], pt[2]] )

    if show :
        AddSpherePts ( pts, (.8,.2,.8,0.5), 0.1, "RAD points %.1f" % RAD )

    if len (pts) < 1 :
        if log :
            print " - no points for RAD %.1f - %d.%s - " % (RAD, atoms[0].residue.id.position, atoms[0].residue.type),
            print "SC" if atoms[0].isSC else "BB"

    else :
        d_vals = dmap.interpolated_values ( pts, atoms[0].molecule.openState.xform )
        avg = numpy.average ( d_vals )
        sdev = numpy.std ( d_vals )

        if sdev < 1e-4 : sdev = 1e-4
        zscore = (avg0 - avg) / sdev #(scores[0] - avg) / stdev
        #print " - scores: avg %.4f, std %.4f, z-score %.4f" % (avg, stdev, zscore )

        if log :
            print " - q at rad %.2f, avg0 %.3f, avg %.3f, stdev %.4f, z %.3f, %d pts" % (RAD, avg0, avg, sdev, zscore, len(pts))


    return zscore






def CurMolAndChain () :

    segModDialog = getdialog ()
    if segModDialog != None :

        if segModDialog.cur_mol == None :
            segModDialog.cur_mol = chimera.Molecule()
            segModDialog.cur_mol.name = "Model"
            #chimera.openModels.add ( [mol], noprefs = True )
            chimera.openModels.add ( [segModDialog.cur_mol] )
            segModDialog.struc.set ( segModDialog.cur_mol.name )

            try :
                segModDialog.cur_mol.openState.xform = chimera.openModels.list()[0].openState.xform
            except :
                pass

        chainId = segModDialog.chain.get()
        if len(chainId) == 0 :
            chainId = "A"
            segModDialog.chain.set ( chainId )

        return segModDialog.cur_mol, chainId

    return None, ""


def VisMapMod () :

    mol, map = None, None

    for m in OML(modelTypes = [chimera.Molecule]) :
        if m.display :
            mol = m

    for m in OML(modelTypes = [VolumeViewer.volume.Volume]) :
        if m.display :
            map = m

    return map, mol



def ZScoresVis ( ) :

    map, mol = VisMapMod()

    if mol != None and map != None :
        ZScores ( mol, map)
    else :
        print "Did not find visible mol and map"



def ZScores ( mol, map ) :

    resolution = 3.0 * map.data.step[0]
    print "Mol: %s, Map: %s -- res %.1f" % (mol.name, map.name, resolution)

    qscores.SetBBAts ( mol )


    cmap = {}
    for r in mol.residues :

        if r.id.chainId in cmap :
            cmap[r.id.chainId].append ( [r.id.position, r] )
        else :
            cmap[r.id.chainId] = [ [r.id.position, r] ]


    #ress = cmap['0']

    allBB, allSC = [], []

    for cid, ress in cmap.iteritems() :
        print " - chain %s" % cid

        ress.sort ()
        ares = [el[1] for el in ress]

        zscores = []
        if 0 :
            sses = SSEs ( ares )
            for el in sses :
                si, ei, ss, elRess = el
                zscore, ccs = zBB ( mol, elRess, resolution, map )
                #print ss, si, "-", ei, zscore
                if zscore != None :
                    zscores.append ( zscore )
                for r in elRess :
                    r.bbZ = zscore

        else :
            bbs = BBsegs ( self.seqRes )
            W = 3
            print " - %d BB segments" % len(bbs)
            for bb in bbs :
                print "  %d res, %d-%d" % (len(bb),bb[0].id.position,bb[-1].id.position)

                for ri, r in enumerate ( bb ) :
                    firstRi = max ( 0, ri-(W-1)/2 )
                    lastRi = min ( len(bb)-1, ri+(W-1)/2 )
                    ress = bb[firstRi:lastRi+1]
                    zscore, ccs = zBB ( self.cur_mol, ress, resolution, map )
                    if zscore != None :
                        zscores.append ( zscore )

        avgBB = 0
        if len(zscores) > 0 :
            avgBB = numpy.average ( zscores )
            allBB.extend ( zscores )
            #print " - BB - min %.2f max %.2f, avg %.2f" % (min(zscores), max(zscores), avgBB )
        #else :
        #    print " - BB - no zscores?"


        avgSC = 0
        zscores = CalcRotaZ ( map, mol, ares )
        if len(zscores) > 0 :
            avgSC = numpy.average(zscores)
            #print " - SC - min %.2f max %.2f, avg %.2f" % (min(zscores), max(zscores), numpy.average(zscores) )
            allSC.extend ( zscores )
        #else :
        #    print " - SC - no zscores?"

        print "Chain %s - %d res - avgBB %.2f, avgSC %.2f" % ( cid, len(ares), avgBB, avgSC )


    print ""

    avgBB = 0
    if len(avgBB) > 0 :
        avgBB = numpy.average(allBB)
        print "BB All - %d scores - min %.2f max %.2f, avg %.2f" % (len(allBB), min(allBB), max(allBB), avgBB )
    else :
        print "BB - no zscores?"

    avgSC = 0
    if len(allSC) > 0 :
        avgSC = numpy.average(allSC)
        print "SC All - %d scores - min %.2f max %.2f, avg %.2f" % (len(allSC), min(allSC), max(allSC), avgSC )
    else :
        print "SC - no zscores?"

    print ""







def CalcRotaZ ( dmap, mol, ress ) :

    A = []
    resolution = 3.0 * dmap.data.step[0]

    for ri, res in enumerate ( ress ) :

        if 1 :
            if res.isProt :
                res.scZ, cc = zRotSideChain ( mol, res, resolution, dmap, show=False )
            elif res.isNA :
                res.scZ = zRotBase ( mol, res, resolution, dmap, show=False )
            else :
                print "?_%d.%s_%s" % (res.id.position, res.id.chainId, res.type)
                res.scZ = 0
            res.qSC = res.scZ

        else :
            res.scZ = zShakeSC ( mol, res, resolution, dmap, show=False )


        if res.scZ != None :
            A.append ( res.scZ )


    #avgA, stdA = numpy.average ( A ), numpy.std ( A )
    #umsg ( "Avg side chain Z-score: %.3f" % ( avgA ) )
    return A



def MoveSC () :

    map, mol = VisMapMod()
    resolution = 3.0 * map.data.step[0]

    print "Map: %s, mol: %s" % (map.name, mol.name)
    res = chimera.selection.currentResidues()[0]
    print " - res: %s %d.%s" % (res.type, res.id.position, res.id.chainId)
    z = None

    if 1 :
        if res.isProt :
            z, cc = zRotSideChain ( mol, res, resolution, map, True )
        elif res.isNA :
            z = zRotBase ( mol, res, resolution, map, True )

    else :
        z = zShakeSC ( mol, res, resolution, map, True )

    print z


def score3 (R) :

    selAts = chimera.selection.currentAtoms()
    if len ( selAts ) == 0 :
        return

    dmap = getdialog ().cur_dmap

    a = selAts[0]
    r = a.residue
    print "Res: %s - %d.%s - %s - Atom: %s" % (r.type, r.id.position, r.id.chainId, r.molecule.name, a.name)

    if not hasattr ( r.molecule, 'bbats' ) :
        qscores.SetBBAts(r.molecule)
        r.molecule.bbats = True

    removeMods = []
    for m in chimera.openModels.list() :
        if "RAD points" in m.name :
            removeMods.append ( m )
    chimera.openModels.remove ( removeMods )


    ats = [at for at in self.cur_mol.atoms if not at.element.name == "H"]
    points = _multiscale.get_atom_coordinates ( ats, transformed = False )
    print " - search tree: %d/%d ats" % ( len(ats), len(r.molecule.atoms) )
    allAtTree = AdaptiveTree ( points.tolist(), ats, 1.0)


    #allAtTree = None
    #print "-"

    import time
    start = time.time()

    #r.sdev = RadAts ( selAts, dmap, allAtTree=allAtTree, show=1, log=0, numPts=40, toRAD=2, dRAD=0.5 )
    r.sigma = RadAts ( r.scAtoms, dmap, allAtTree=allAtTree, show=0, log=1, numPts=20, toRAD=2, dRAD=0.5 )

    end = time. time()

    print "%s - rad: %.3f, time: %f" % ( a.name, r.sigma, (end - start) )

    scZ, cc = zRotSideChain ( r.molecule, r, R, dmap, show=False )
    print " - cc %.3f, scZ %.3f " % (cc, scZ)
    print "%f\t%f\t%f" % (r.sigma, cc, scZ)




def zShakeSC ( mol, res, resolution, dmap, show=False ) :

    atoms = res.scAtoms

    if len(atoms) < 1 :
        #print " - no sc atoms" % len(atoms)
        return None

    score0 = 0
    scores, scorest = [], []
    T = 1
    trange = [-T*1.0, 0.0, T*1.0]
    #trange = [-T*2.0, -T, 0.0, T, T*2.0]

    fout = None
    if show :
        fout = open ("/Users/greg/Desktop/sc.txt", "w")

    moved = False

    for xx in trange :
        for yy in trange :
            for zz in trange :

                v = chimera.Vector(xx,yy,zz)
                xfT = chimera.Xform.translation ( chimera.Vector(xx,yy,zz) )

                molg = MyMolMapX ( mol, atoms, resolution, dmap.data.step[0], xfT )

                fpoints, fpoint_weights = fit_points_g ( molg )
                map_values = dmap.interpolated_values ( fpoints, mol.openState.xform )
                olap, corr1, corr2 = FitMap.overlap_and_correlation ( fpoint_weights, map_values )

                if numpy.fabs(xx) < .01 and numpy.fabs(yy) < .01 and numpy.fabs(zz) < .01 :
                    score0 = corr1
                else :
                    scores.append ( corr1 )
                    if fout :

                        #if not moved :
                        nmol, cress = CopyRess ( [res] )
                        for nr in cress :
                            for nat in nr.atoms :
                                try :
                                    nat.setCoord ( xfT.apply ( nat.coord() ) )
                                except :
                                    pass
                        #chimera.openModels.add ( [nmol] )
                        nmol.name = "S_%.0f_%.0f_%.0f" % (xx,yy,zz)
                        moved = True

                        scorest.append ( [corr1, [xx,yy,zz], nmol] )


    if fout :
        scorest.sort ()
        #scorest.reverse ()
        scorest = scorest[0:len(scorest)/2]
        if fout :
            fout.write ( "%.0f,%.0f,%.0f\t%f\n" % (0,0,0, score0) )
            for sc, t, nmol in scorest:
                fout.write ( "%.0f,%.0f,%.0f\t%f\n" % (t[0],t[1],t[2], sc) )
                chimera.openModels.add ( [nmol] )
                qscores.SetBBAts ( nmol )
                for at in nmol.atoms :
                    at.display = at.isSC


        fout.close()

    if 1 :
        scores.sort ()
        #scores.reverse ()
        scores = scores[0:len(scores)/2]

    #print ""
    avg = numpy.average ( scores )  #numpy.average ( scores[1:] )
    stdev = numpy.std ( scores ) #numpy.std ( scores[1:] )
    if stdev < 1e-8 :
        #print " - nostdev"
        return None
    zscore = (score0 - avg) / stdev #(scores[0] - avg) / stdev
    #print " - scores: avg %.4f, std %.4f, z-score %.4f" % (avg, stdev, zscore )
    #fout.close()

    return zscore




def zRotSideChain ( mol, r, resolution, dmap, show=False ) :

    r.CA, r.CB, r.CG = None, None, None
    try :
        r.CA = r.atomsMap["CA"][0]
        r.CB = r.atomsMap["CB"][0]
    except :
        pass

    if "CG" in r.atomsMap :
        r.CG = r.atomsMap["CG"][0]
    elif "CG1" in r.atomsMap :
        r.CG = r.atomsMap["CG1"][0]
    elif "CG2" in r.atomsMap :
        r.CG = r.atomsMap["CG2"][0]
    elif "OG" in r.atomsMap :
        r.CG = r.atomsMap["OG"][0]
    elif "SG" in r.atomsMap :
        r.CG = r.atomsMap["SG"][0]

    if r.CA == None or r.CB == None or r.CG == None :
        #print r.type, " - no ats"
        return None, None

    resolution = 3.0 * dmap.data.step[0]

    scores = []

    #molg = MyMolMap ( mol, r.atoms, resolution, dmap.data.step[0] )
    #fpoints, fpoint_weights = fit_points_g ( molg )
    #map_values = dmap.interpolated_values ( fpoints, mol.openState.xform )
    #olap_0, corr1_0, corr2_0 = FitMap.overlap_and_correlation ( fpoint_weights, map_values )

    rats = r.scAtoms
    nrats = []
    for at in rats :
        try :
            at.p0 = at.coord()
            nrats.append ( at )
        except :
            pass

    fout = None
    if show :
        fout = open ("/Users/greg/Desktop/sc.txt", "w")

    #for ri, rmol in enumerate ( rmols[0:10] ) :
    for deg in range (0, 360, 36) :

        RotAts ( nrats, r.CA, r.CB, deg )

        if fout :
            nmol, cress = CopyRess ( [r] )
            chimera.openModels.add ( [nmol] )
            nmol.name = "SC %d %.0f" % (r.id.position, deg)
            nr = nmol.residues[0]
            qscores.SetBBAts ( nmol )
            for at in nr.atoms :
                if at.isBB :
                    at.display = False
                else :
                    at.display = True

        corr = ResCC ( mol, nrats, resolution, dmap )
        scores.append ( corr )

        for at in nrats :
            at.setCoord ( at.p0 )

    if fout :
        for sci, sc in enumerate ( scores ):
            fout.write ( "%d\t%f\n" % (sci*36, sc) )

        fout.close()

    zscore1 = None
    if len(scores) > 3 :
        avg = numpy.average ( scores[1:] )
        stdev = numpy.std ( scores[1:] )
        zscore1 = ( (scores[0] - avg) / stdev ) if stdev > 1e-5 else 0
        #print " -0- avg %.4f, std %.4f, z-score %.4f" % (avg, stdev, zscore1 )
        #print scores
        #print " -1- avg %.4f, std %.4f, z-score %.4f" % (avg, stdev, zscore1 )


    return zscore1, scores[0]




def zRotBase ( mol, r, resolution, dmap, show=False ) :

    resolution = 3.0 * dmap.data.step[0]

    scores = []

    rats = r.scAtoms
    nrats = []
    for at in rats :
        try :
            if at.element.name == "H" :
                continue
            at.p0 = at.coord()
            nrats.append ( at )
        except :
            pass

    fout = None
    if show :
        fout = open ("/Users/greg/Desktop/sc.txt", "w")

    #for ri, rmol in enumerate ( rmols[0:10] ) :
    for deg in range (0, 360, 36) :

        RotAts ( nrats, r.atomsMap["C1'"][0], r.baseAt, deg )

        if fout :
            nmol, cress = CopyRess ( [r] )
            chimera.openModels.add ( [nmol] )
            nmol.name = "SC %d %.0f" % (r.id.position, deg)
            nr = nmol.residues[0]
            qscores.SetBBAts ( nmol )
            for at in nr.atoms :
                if at.isBB :
                    at.display = False
                else :
                    at.display = True

        corr = ResCC ( mol, nrats, resolution, dmap )
        scores.append ( corr )

        for at in nrats :
            at.setCoord ( at.p0 )

    if fout :
        for sci, sc in enumerate ( scores ):
            fout.write ( "%d\t%f\n" % (sci*36, sc) )

        fout.close()

    zscore1 = None
    if len(scores) > 3 :
        avg = numpy.average ( scores[1:] )
        stdev = numpy.std ( scores[1:] )
        zscore1 = ( (scores[0] - avg) / stdev ) if stdev > 1e-5 else 0
        #print " -1- avg %.4f, std %.4f, z-score %.4f" % (avg, stdev, zscore1 )


    return zscore1







def MoveBB () :

    map, mol = VisMapMod()
    resolution = 3.0 * map.data.step[0]

    print "Map: %s, mol: %s" % (map.name, mol.name)
    z, cc = zBB ( mol, chimera.selection.currentResidues(), resolution, map, True )
    print z



def zBB ( mol, ress, resolution, dmap, show=False ) :

    atoms = []
    for r in ress :
        #if 'C' in r.atomsMap : atoms.append ( r.atomsMap['C'][0] )
        #if 'N' in r.atomsMap : atoms.append ( r.atomsMap['N'][0] )
        #if 'CA' in r.atomsMap : atoms.append ( r.atomsMap['CA'][0] )
        #if 'O' in r.atomsMap : atoms.append ( r.atomsMap['O'][0] )
        atoms.extend ( r.bbAtoms )
        atoms.extend ( r.scAtoms )

    if len(atoms) < 1 :
        #print " - no atoms" % len(atoms)
        return [0,0]

    score0 = 0
    scores, scorest = [], []
    T = 2
    trange = [-T*1.0, 0.0, T*1.0]
    #trange = [-T*2.0, -T, 0.0, T, T*2.0]

    fout = None
    if show :
        fout = open ("/Users/greg/Desktop/sse.txt", "w")

    moved = False

    for xx in trange :
        for yy in trange :
            for zz in trange :

                v = chimera.Vector(xx,yy,zz)
                xfT = chimera.Xform.translation ( chimera.Vector(xx,yy,zz) )

                molg = MyMolMapX ( mol, atoms, resolution, dmap.data.step[0], xfT )

                fpoints, fpoint_weights = fit_points_g ( molg )
                map_values = dmap.interpolated_values ( fpoints, mol.openState.xform )
                olap, corr1, corr2 = FitMap.overlap_and_correlation ( fpoint_weights, map_values )

                if numpy.fabs(xx) < .01 and numpy.fabs(yy) < .01 and numpy.fabs(zz) < .01 :
                    score0 = corr2
                else :
                    scores.append ( corr2 )
                    if fout :
                        scorest.append ( [corr2, [xx,yy,zz]] )

                        if not moved :
                            nmol, cress = CopyRess ( ress )
                            for nr in cress :
                                for nat in nr.atoms :
                                    try :
                                        nat.setCoord ( xfT.apply ( nat.coord() ) )
                                    except :
                                        pass
                            chimera.openModels.add ( [nmol] )
                            nmol.name = "T_%.0f_%.0f_%.0f" % (xx,yy,zz)
                            moved = True


    if fout :
        scorest.sort ()
        scorest.reverse ()
        scorest = scorest[len(scorest)/2:]
        if fout :
            fout.write ( "%.0f,%.0f,%.0f\t%f\n" % (0,0,0, score0) )
            for sc, t in scorest:
                fout.write ( "%.0f,%.0f,%.0f\t%f\n" % (t[0],t[1],t[2], sc) )

        fout.close()

    if 0 :
        scores.sort ()
        scores.reverse ()
        scores = scores[len(scores)/2:]

    #print ""
    avg = numpy.average ( scores )  #numpy.average ( scores[1:] )
    stdev = numpy.std ( scores ) #numpy.std ( scores[1:] )
    if stdev < 1e-8 :
        #print " - nostdev"
        return [0,0]
    zscore = (score0 - avg) / stdev #(scores[0] - avg) / stdev
    #print " - scores: avg %.4f, std %.4f, z-score %.4f" % (avg, stdev, zscore )
    #fout.close()

    return [zscore, score0]







def CopyRess ( res ) :

    nmol = chimera.Molecule()
    ress = [None] * len ( res )

    aMap = dict()
    for ri, r in enumerate ( res ) :
        nres = nmol.newResidue (r.type, chimera.MolResId(r.id.chainId, r.id.position))
        ress[ri] = nres
        for at in r.atoms :
            nat = nmol.newAtom (at.name, chimera.Element(at.element.number))
            aMap[at] = nat
            nres.addAtom( nat )
            p = chimera.Point ( at.coord().x, at.coord().y, at.coord().z )
            nat.setCoord ( p )
            nat.coord0 = chimera.Point ( at.coord().x, at.coord().y, at.coord().z )
            #if at.name == "C" or at.name == 'CA' or at.name == 'O' or at.name == "N" :
            #    at.display = False



    for bond in res[0].molecule.bonds :
        try :
            nb = nmol.newBond ( aMap[bond.atoms[0]], aMap[bond.atoms[1]] )
            nb.display = nb.Smart
        except :
            pass

    for r in ress :
        r.CA, r.CB, r.CG = None, None, None
        try :
            r.CA = r.atomsMap["CA"][0]
            r.CB = r.atomsMap["CB"][0]
            r.CG = r.atomsMap["CG"][0]
        except :
            pass

    return nmol, ress



def RotAts (rats, a1, a2, deg) :

    # phi: N -> CA
    p1, p2 = a1.coord(), a2.coord()
    v = p2 - p1; v.normalize()

    xf = chimera.Xform.translation ( p1.toVector() )
    xf.multiply ( chimera.Xform.rotation ( v, deg ) )
    xf.multiply ( chimera.Xform.translation ( p1.toVector() * -1.0 ) )

    #for at in res.atoms :
    #    if at.name != 'C' and at.name != 'CA' and at.name != 'N' and at.name != 'CB' and at.name != 'O' :
    for at in rats :
            at.setCoord ( xf.apply (at.coord()) )



def MyMolMapX2 ( atoms, resolution, step=1.0, xf=None ) :

    from math import sqrt, pi

    pad = 3*resolution
    cutoff_range = 5 # in standard deviations
    sigma_factor = 1/(pi*sqrt(2)) # standard deviation / resolution

    from _multiscale import get_atom_coordinates
    xyz = get_atom_coordinates(atoms, transformed = False)

    # Transform coordinates to local coordinates of the molecule containing
    # the first atom.  This handles multiple unaligned molecules.
    # Or if on_grid is specified transform to grid coordinates.
    #m0 = atoms[0].molecule

    #xf = m0.openState.xform
    #import Matrix as M
    #M.transform_points(xyz, M.xform_matrix(xf.inverse()))

    anum = [a.element.number for a in atoms]

    grid = bounding_grid(xyz, step, pad, [])
    grid.name = ""

    sdev = resolution * sigma_factor
    add_gaussians(grid, xyz, anum, sdev, cutoff_range, [])

    #return grid, molecules
    return grid




def molecule_grid_dataX (m0, atoms, resolution, step, pad, xfT, cutoff_range, sigma_factor, transforms = [], csys = None):

    from _multiscale import get_atom_coordinates
    xyz = get_atom_coordinates(atoms, transformed = True)

    # Transform coordinates to local coordinates of the molecule containing
    # the first atom.  This handles multiple unaligned molecules.
    # Or if on_grid is specified transform to grid coordinates.
    #m0 = atoms[0].molecule
    xf = m0.openState.xform
    xf.multiply ( xfT )
    import Matrix as M
    M.transform_points(xyz, M.xform_matrix(xf.inverse()))
    if csys:
        xf.premultiply(csys.xform.inverse())
    tflist = M.coordinate_transform_list(transforms, M.xform_matrix(xf))

    anum = [a.element.number for a in atoms]

    molecules = set([a.molecule for a in atoms])
    if len(molecules) > 1:
        name = 'molmap res %.3g' % (resolution,)
    else:
        name = 'molmap %s res %.3g' % (m0.name, resolution)

    grid = bounding_grid(xyz, step, pad, tflist)
    grid.name = name

    sdev = resolution * sigma_factor
    add_gaussians(grid, xyz, anum, sdev, cutoff_range, tflist)

    #return grid, molecules
    return grid


def MyMolMapX ( m0, atoms, resolution, step, xf ) :

    #from MoleculeMap import molecule_grid_data
    from math import sqrt, pi
    #from chimera import openModels as om
    #from VolumeViewer import volume_from_grid_data

    atoms = tuple(atoms)

    pad = 3*resolution
    cutoff_range = 5 # in standard deviations
    sigma_factor = 1/(pi*sqrt(2)) # standard deviation / resolution
    transforms,csys = [], None
    display_threshold = 0.95

    return molecule_grid_dataX (m0, atoms, resolution, step, pad, xf, cutoff_range, sigma_factor, transforms, csys)



def MyMolMap ( m0, atoms, resolution, step ) :

    #from MoleculeMap import molecule_grid_data
    from math import sqrt, pi
    from chimera import openModels as om
    from VolumeViewer import volume_from_grid_data

    atoms = tuple(atoms)

    pad = 3*resolution
    cutoff_range = 5 # in standard deviations
    sigma_factor = 1/(pi*sqrt(2)) # standard deviation / resolution
    transforms,csys = [], None
    #display_threshold = 0.95

    return molecule_grid_data(m0, atoms, resolution, step, pad, None, cutoff_range, sigma_factor, transforms, csys)




def molecule_grid_data(m0, atoms, resolution, step, pad, on_grid,
                       cutoff_range, sigma_factor,
                       transforms = [], csys = None):



    from _multiscale import get_atom_coordinates
    xyz = get_atom_coordinates(atoms, transformed = True)

    # Transform coordinates to local coordinates of the molecule containing
    # the first atom.  This handles multiple unaligned molecules.
    # Or if on_grid is specified transform to grid coordinates.
    #m0 = atoms[0].molecule
    xf = on_grid.openState.xform if on_grid else m0.openState.xform
    import Matrix as M
    M.transform_points(xyz, M.xform_matrix(xf.inverse()))
    if csys:
        xf.premultiply(csys.xform.inverse())
    tflist = M.coordinate_transform_list(transforms, M.xform_matrix(xf))

    anum = [a.element.number for a in atoms]

    molecules = set([a.molecule for a in atoms])
    if len(molecules) > 1:
        name = 'molmap res %.3g' % (resolution,)
    else:
        name = 'molmap %s res %.3g' % (m0.name, resolution)

    if on_grid:
        from numpy import float32
        grid = on_grid.region_grid(on_grid.region, float32)
    else:
        grid = bounding_grid(xyz, step, pad, tflist)
    grid.name = name

    sdev = resolution * sigma_factor
    add_gaussians(grid, xyz, anum, sdev, cutoff_range, tflist)

    #return grid, molecules
    return grid




def ResCC ( mol, rats, resolution, dmap ) :

    molg = MyMolMap ( mol, rats, resolution, dmap.data.step[0] )

    #if 0 :
    #    fmap = VolumeViewer.volume.volume_from_grid_data ( molg )
    #    fmap.name = "res molmap!"
    #    fpoints, fpoint_weights = fit_points(fmap, False)
    #    map_values = dmap.interpolated_values ( fpoints, fmap.openState.xform )
    #    olap, corr1, corr2 = FitMap.overlap_and_correlation ( fpoint_weights, map_values )
    #    scores.append ( corr1 )
    #    chimera.openModels.close ( [fmap] )
    #else :

    fpoints, fpoint_weights = fit_points_g ( molg, 0.1 )
    map_values = dmap.interpolated_values ( fpoints, mol.openState.xform )
    olap, corr1, corr2 = FitMap.overlap_and_correlation ( fpoint_weights, map_values )
    return corr1, corr2



def fit_points_g (fdata, threshold = 0.3) :

    mat = fdata.full_matrix()

    import _volume
    points = _volume.high_indices(mat, threshold)
    fpoints = points.astype(numpy.single)
    fpoint_weights = mat[points[:,2],points[:,1],points[:,0]]

    nz = numpy.nonzero( fpoint_weights )[0]
    if len(nz) < len (fpoint_weights) :
        fpoints = numpy.take( fpoints, nz, axis=0 )
        fpoint_weights = numpy.take(fpoint_weights, nz, axis=0)

    transform_vertices ( fpoints, fdata.ijk_to_xyz_transform )

    if 0 : print "FitPoints from %s with threshold %.4f, %d nonzero" % (
        fmap.name, threshold, len(nz) )

    return fpoints, fpoint_weights


def fit_points (fmap, threshold = 1e-5) :

    mat = fmap.data.full_matrix()

    import _volume
    points = _volume.high_indices(mat, threshold)
    fpoints = points.astype(numpy.single)
    fpoint_weights = mat[points[:,2],points[:,1],points[:,0]]

    nz = numpy.nonzero( fpoint_weights )[0]
    if len(nz) < len (fpoint_weights) :
        fpoints = numpy.take( fpoints, nz, axis=0 )
        fpoint_weights = numpy.take(fpoint_weights, nz, axis=0)

    from _contour import affine_transform_vertices as transform_vertices
    transform_vertices ( fpoints, fmap.data.ijk_to_xyz_transform )
    #transform_vertices ( fpoints, Matrix.xform_matrix( fmap.openState.xform ) )

    if 0 : print "FitPoints from %s with threshold %.4f, %d nonzero" % (
        fmap.name, threshold, len(nz) )

    return fpoints, fpoint_weights





# -----------------------------------------------------------------------------
#
def bounding_grid(xyz, step, pad, transforms):

    xyz_min, xyz_max = point_bounds(xyz, transforms)
    origin = [x-pad for x in xyz_min]
    from math import ceil
    shape = [int(ceil((xyz_max[a] - xyz_min[a] + 2*pad) / step)) for a in (2,1,0)]
    from numpy import zeros, float32
    matrix = zeros(shape, float32)
    from VolumeData import Array_Grid_Data
    grid = Array_Grid_Data(matrix, origin, (step,step,step))
    return grid


# -----------------------------------------------------------------------------
#
def add_gaussians(grid, xyz, weights, sdev, cutoff_range, transforms = []):

    from numpy import zeros, float32, empty
    sdevs = zeros((len(xyz),3), float32)
    for a in (0,1,2):
        sdevs[:,a] = sdev / grid.step[a]

    import Matrix as M
    if len(transforms) == 0:
        transforms = [M.identity_matrix()]
    from _gaussian import sum_of_gaussians
    ijk = empty(xyz.shape, float32)
    matrix = grid.matrix()
    for tf in transforms:
        ijk[:] = xyz
        M.transform_points(ijk, M.multiply_matrices(grid.xyz_to_ijk_transform, tf))
        sum_of_gaussians(ijk, weights, sdevs, cutoff_range, matrix)

    from math import pow, pi
    normalization = pow(2*pi,-1.5)*pow(sdev,-3)
    matrix *= normalization



# -----------------------------------------------------------------------------
#
def point_bounds(xyz, transforms = []):

    from _multiscale import bounding_box
    if transforms :
        from numpy import empty, float32
        xyz0 = empty((len(transforms),3), float32)
        xyz1 = empty((len(transforms),3), float32)
        txyz = empty(xyz.shape, float32)
        import Matrix as M
        for i, tf in enumerate(transforms) :
            txyz[:] = xyz
            M.transform_points(txyz, tf)
            xyz0[i,:], xyz1[i,:] = bounding_box(txyz)
        xyz_min, xyz_max = xyz0.min(axis = 0), xyz1.max(axis = 0)
    else:
        xyz_min, xyz_max = bounding_box(xyz)

    return xyz_min, xyz_max





# ---------------------------------------------------------------------------------




def SkinMap ( atoms, bones, N, dmap, atomRad, nname, showMesh = False ) :

    from _multiscale import get_atom_coordinates
    points = get_atom_coordinates ( atoms, transformed = True )

    import _contour
    points0 = numpy.copy ( points )
    _contour.affine_transform_vertices ( points0, Matrix.xform_matrix( dmap.openState.xform.inverse() ) )


    nn3, nn2, nn1 = dmap.data.size

    npoints = VolumeData.grid_indices ( (int(nn1), int(nn2), int(nn3) ), numpy.single)  # i,j,k indices
    _contour.affine_transform_vertices ( npoints, dmap.data.ijk_to_xyz_transform )
    _contour.affine_transform_vertices ( npoints, Matrix.xform_matrix( dmap.openState.xform ) )
    _contour.affine_transform_vertices ( npoints, Matrix.xform_matrix( atoms[0].molecule.openState.xform.inverse() ) )

    for bo in bones :
        bo.MakeFrame ()

    if N == 1 :
        for pi, p in enumerate ( npoints ) :

            cbone, minDist = None, 1e9
            for bo in bones :
                d = bo.DistToPoint ( p )

                if d < minDist :
                    minDist = d
                    cbone = bo

            pt = cbone.SkinPoint ( p )
            npoints[pi] = pt

    else :

        for pi, p in enumerate ( npoints ) :

            dbos = []
            for bo in bones :
                dbos.append ( [bo.DistToPoint ( p ), bo] )

            dbos.sort()

            totD = 0.0
            sp = numpy.array ( [0,0,0] )
            for i in range ( N ) :
                d, bo = dbos[i]
                sp = sp + numpy.array ( bo.SkinPoint ( p ) ) * d
                totD += d

            npoints[pi] = sp / totD


    _contour.affine_transform_vertices ( npoints, Matrix.xform_matrix( atoms[0].molecule.openState.xform ) )
    _contour.affine_transform_vertices ( npoints, Matrix.xform_matrix( dmap.openState.xform.inverse() ) )

    dvals = dmap.interpolated_values ( npoints, dmap.openState.xform )
    nmat = dvals.reshape( (nn3,nn2,nn1) )
    ndata = VolumeData.Array_Grid_Data ( nmat, dmap.data.origin, dmap.data.step, dmap.data.cell_angles )

    mdata = VolumeData.zone_masked_grid_data ( ndata, points0, atomRad )

    MapFromData ( mdata, nname, dmap, False )
    if showMesh :
        MapFromData ( mdata, nname, dmap, True )



def ExtractDen ( atoms, dmap, nname, boundRad = 2.0, showMesh = False) :

    from _multiscale import get_atom_coordinates
    points1 = get_atom_coordinates ( atoms, transformed = False )
    #COM, U, S, V = prAxes ( points )

    bound = 4.0
    li,lj,lk = numpy.min ( points1, axis=0 ) - (bound, bound, bound)
    hi,hj,hk = numpy.max ( points1, axis=0 ) + (bound, bound, bound)

    nstep = (dmap.data.step[0], dmap.data.step[1], dmap.data.step[2] )

    n1 = int ( numpy.ceil ( (hi - li + 1) / nstep[0] ) )
    n2 = int ( numpy.ceil ( (hj - lj + 1) / nstep[1] ) )
    n3 = int ( numpy.ceil ( (hk - lk + 1) / nstep[2] ) )

    O = chimera.Point ( li, lj, lk  )
    #O = atoms[0].molecule.openState.xform.apply ( O )

    #print " - new map origin:", nO

    npoints = VolumeData.grid_indices ( (n1, n2, n3), numpy.single)  # i,j,k indices
    S = dmap.data.step

    _contour.affine_transform_vertices ( npoints, ((S[0], 0.0, 0.0, O[0]), (0.0, S[1], 0.0, O[1]), (0.0, 0.0, S[1], O[2])) )
    #_contour.affine_transform_vertices ( npoints, Matrix.xform_matrix( dmap.openState.xform.inverse() ) )

    dvals = dmap.interpolated_values ( npoints, atoms[0].molecule.openState.xform )
    nmat = dvals.reshape( (n3,n2,n1) )

    ndata = VolumeData.Array_Grid_Data ( nmat, O, nstep, dmap.data.cell_angles, name = nname )

    #_contour.affine_transform_vertices ( points1, Matrix.xform_matrix( atoms[0].molecule.openState.xform ) )
    #_contour.affine_transform_vertices ( points, Matrix.xform_matrix( dmap.openState.xform.inverse() ) )
    mdata = VolumeData.zone_masked_grid_data ( ndata, points1, boundRad )

    dmap = MapFromData ( mdata, nname, dmap, False, color=(.7,.7,.7,.2) )
    dmap.openState.xform = atoms[0].molecule.openState.xform
    dmesh = None

    if showMesh :
        dmesh = MapFromData ( mdata, nname, dmap, True )
        dmesh.openState.xform = atoms[0].molecule.openState.xform

    return [dmap, dmesh]






def BoneMap ( bone, dmap, atomRad, nname, show = False, showMesh = False ) :

    #_contour.affine_transform_vertices ( points, Matrix.xform_matrix( dmap.openState.xform.inverse() ) )
    #mdata = VolumeData.zone_masked_grid_data ( dmap.data, points, atomRad )

    from _multiscale import get_atom_coordinates
    atoms = [bone.a1, bone.a2]
    points = get_atom_coordinates ( atoms, transformed = True )

    import _contour
    points1 = numpy.copy ( points )
    _contour.affine_transform_vertices ( points1, Matrix.xform_matrix( dmap.openState.xform.inverse() ) )
    points0 = numpy.copy ( points1 )
    _contour.affine_transform_vertices ( points1, dmap.data.xyz_to_ijk_transform )

    bound = int ( numpy.ceil( atomRad / dmap.data.step[0] ) ) + 1
    li,lj,lk = numpy.min ( points1, axis=0 ) - (bound, bound, bound)
    hi,hj,hk = numpy.max ( points1, axis=0 ) + (bound, bound, bound)

    n1 = hi - li + 1
    n2 = hj - lj + 1
    n3 = hk - lk + 1

    #print " - bounds - %d %d %d --> %d %d %d --> %d %d %d" % ( li,lj,lk, hi,hj,hk, n1,n2,n3 )

    #nmat = numpy.zeros ( (n1,n2,n3), numpy.float32 )
    #dmat = dmap.full_matrix()

    nstep = (dmap.data.step[0], dmap.data.step[1], dmap.data.step[2] )
    #nstep = (fmap.data.step[0]/2.0, fmap.data.step[1]/2.0, fmap.data.step[2]/2.0 )

    nn1 = int ( round (dmap.data.step[0] * float(n1) / nstep[0]) )
    nn2 = int ( round (dmap.data.step[1] * float(n2) / nstep[1]) )
    nn3 = int ( round (dmap.data.step[2] * float(n3) / nstep[2]) )

    O = dmap.data.origin
    #print " - %s origin:" % dmap.name, O
    nO = ( O[0] + float(li) * dmap.data.step[0],
           O[1] + float(lj) * dmap.data.step[1],
           O[2] + float(lk) * dmap.data.step[2] )

    #print " - new map origin:", nO

    wmat = numpy.zeros ( (nn3,nn2,nn1), numpy.float32 )
    ndata = VolumeData.Array_Grid_Data ( wmat, nO, nstep, dmap.data.cell_angles )

    npoints = VolumeData.grid_indices ( (nn1, nn2, nn3), numpy.single)  # i,j,k indices
    npointsi = numpy.copy ( npoints )
    _contour.affine_transform_vertices ( npoints, ndata.ijk_to_xyz_transform )
    _contour.affine_transform_vertices ( npoints, Matrix.xform_matrix( dmap.openState.xform ) )
    _contour.affine_transform_vertices ( npoints, Matrix.xform_matrix( atoms[0].molecule.openState.xform.inverse() ) )

    for pi, p in enumerate ( npoints ) :

        i,j,k = npointsi[pi]
        d = bone.DistToPoint ( p )
        if d < atomRad :
            wmat[k,j,i] = 1.0
        else :
            wmat[k,j,i] = 1.0 / numpy.power (1+d-atomRad,8)

    _contour.affine_transform_vertices ( npoints, Matrix.xform_matrix( atoms[0].molecule.openState.xform ) )
    _contour.affine_transform_vertices ( npoints, Matrix.xform_matrix( dmap.openState.xform.inverse() ) )

    dvals = dmap.interpolated_values ( npoints, dmap.openState.xform )
    nmat = dvals.reshape( (nn3,nn2,nn1) )

    bone.ndata = VolumeData.Array_Grid_Data ( nmat*wmat, nO, nstep, dmap.data.cell_angles, name = nname )
    bone.xfmod = dmap

    if show :

        from random import random as rand
        clr = ( rand()*.5+.1, rand()*.5+.1, rand()*.5+.1 )

        bone.dmap = MapFromData ( bone.ndata, nname, dmap, showMesh, color = clr )
        bone.dmap.openState.xform = dmap.openState.xform





def MoldMap ( atoms, bones, dmap, nname, showMesh = False ) :


    ndata = dmap.data
    nn3, nn2, nn1 = dmap.data.size
    nO = dmap.data.origin
    nmat = numpy.zeros ( (nn3,nn2,nn1), numpy.float32 )
    nstep = (dmap.data.step[0], dmap.data.step[1], dmap.data.step[2] )

    if 1 :
        ndata = DataForAtoms ( atoms, dmap )

    npoints = VolumeData.grid_indices ( (nn1, nn2, nn3), numpy.single)  # i,j,k indices
    _contour.affine_transform_vertices ( npoints, ndata.ijk_to_xyz_transform )

    #_contour.affine_transform_vertices ( npoints, Matrix.xform_matrix( dmap.openState.xform ) )
    #_contour.affine_transform_vertices ( npoints, Matrix.xform_matrix( atoms[0].molecule.openState.xform.inverse() ) )

    for bone in bones :

        npointsc = numpy.copy ( npoints )

        _contour.affine_transform_vertices ( npointsc, Matrix.xform_matrix( bone.Xf().inverse() ) )
        _contour.affine_transform_vertices ( npointsc, Matrix.xform_matrix( bone.Xf0() ) )

        #_contour.affine_transform_vertices ( npoints, Matrix.xform_matrix( atoms[0].molecule.openState.xform ) )
        #_contour.affine_transform_vertices ( npoints, Matrix.xform_matrix( dmap.openState.xform.inverse() ) )

        _contour.affine_transform_vertices ( npointsc, bone.ndata.xyz_to_ijk_transform )

        p2mt = Matrix.xform_matrix ( chimera.Xform.identity() )
        #dvals, outvals = VolumeData.interpolate_volume_data ( npointsc, p2mt, bone.dmap.data.matrix(), method='linear' )
        dvals, outvals = VolumeData.interpolate_volume_data ( npointsc, p2mt, bone.ndata.matrix(), method='linear' )

        bmat = dvals.reshape( (nn3,nn2,nn1) )
        #nmat = nmat + bmat
        nmat = numpy.maximum ( nmat, bmat )

    #nmat = nmat / float ( len(bones) )

    ndata = VolumeData.Array_Grid_Data ( nmat, nO, nstep, dmap.data.cell_angles, name = nname )

    MapFromData ( ndata, nname, dmap, False )
    if showMesh :
        MapFromData ( ndata, nname, dmap, True )




def MoldMap2 ( bones, dmap, dmesh ) :


    ndata = dmap.data
    nn1, nn2, nn3 = dmap.data.size
    nO = dmap.data.origin
    nmat = numpy.zeros ( (nn3,nn2,nn1), numpy.float32 )
    nstep = (dmap.data.step[0], dmap.data.step[1], dmap.data.step[2] )

    npoints = VolumeData.grid_indices ( (nn1, nn2, nn3), numpy.single)  # i,j,k indices
    _contour.affine_transform_vertices ( npoints, ndata.ijk_to_xyz_transform )

    _contour.affine_transform_vertices ( npoints, Matrix.xform_matrix( dmap.openState.xform ) )
    _contour.affine_transform_vertices ( npoints, Matrix.xform_matrix (bones[0].a1.molecule.openState.xform.inverse()) )

    for bone in bones :

        npointsc = numpy.copy ( npoints )

        _contour.affine_transform_vertices ( npointsc, Matrix.xform_matrix( bone.Xf().inverse() ) )
        _contour.affine_transform_vertices ( npointsc, Matrix.xform_matrix( bone.Xf0() ) )

        _contour.affine_transform_vertices ( npointsc, Matrix.xform_matrix (bone.a1.molecule.openState.xform) )
        _contour.affine_transform_vertices ( npointsc, Matrix.xform_matrix( dmap.openState.xform.inverse() ) )

        _contour.affine_transform_vertices ( npointsc, bone.ndata.xyz_to_ijk_transform )

        p2mt = Matrix.xform_matrix ( chimera.Xform.identity() )
        #dvals, outvals = VolumeData.interpolate_volume_data ( npointsc, p2mt, bone.dmap.data.matrix(), method='linear' )
        dvals, outvals = VolumeData.interpolate_volume_data ( npointsc, p2mt, bone.ndata.matrix(), method='linear' )

        bmat = dvals.reshape( (nn3,nn2,nn1) )
        #nmat = nmat + bmat
        nmat = numpy.maximum ( nmat, bmat )

    #nmat = nmat / float ( len(bones) )

    #ndata = VolumeData.Array_Grid_Data ( nmat, nO, nstep, dmap.data.cell_angles, name = nname )
    dmap.data.full_matrix()[:,:,:] = nmat[:,:,:]
    dmap.data.values_changed()
    MapUp ( dmap, False )

    if dmesh != None :
        dmesh.data.full_matrix()[:,:,:] = nmat[:,:,:]
        dmesh.data.values_changed()
        MapUp ( dmesh, True )




def DataForAtoms ( atoms, dmap, nname = "data for atoms" ) :

    from _multiscale import get_atom_coordinates
    points = get_atom_coordinates ( atoms, transformed = True )

    points1 = numpy.copy ( points )
    _contour.affine_transform_vertices ( points1, Matrix.xform_matrix( dmap.openState.xform.inverse() ) )
    #points0 = numpy.copy ( points1 )
    _contour.affine_transform_vertices ( points1, dmap.data.xyz_to_ijk_transform )

    bound = 5
    li,lj,lk = numpy.min ( points1, axis=0 ) - (bound, bound, bound)
    hi,hj,hk = numpy.max ( points1, axis=0 ) + (bound, bound, bound)

    n1 = hi - li + 1
    n2 = hj - lj + 1
    n3 = hk - lk + 1

    nstep = (dmap.data.step[0], dmap.data.step[1], dmap.data.step[2] )

    nn1 = int ( round (dmap.data.step[0] * float(n1) / nstep[0]) )
    nn2 = int ( round (dmap.data.step[1] * float(n2) / nstep[1]) )
    nn3 = int ( round (dmap.data.step[2] * float(n3) / nstep[2]) )

    O = dmap.data.origin
    nO = ( O[0] + float(li) * dmap.data.step[0],
           O[1] + float(lj) * dmap.data.step[1],
           O[2] + float(lk) * dmap.data.step[2] )

    nmat = numpy.zeros ( (nn3,nn2,nn1), numpy.float32 )
    ndata = VolumeData.Array_Grid_Data ( nmat, nO, nstep, dmap.data.cell_angles, name = nname )
    return ndata



def MapForAtoms ( atoms, dmap, nname, showMesh=False, thrF = 1.0 ) :

    ndata = DataForAtoms ( atoms, dmap, nname )

    m1 = MapFromData ( ndata, nname, dmap, False, thrF=thrF )
    m2 = None

    if showMesh :
        m2 = MapFromData ( ndata, nname, dmap, True, thrF=thrF )

    return [m1,m2]


def MapUp (dmap, showMesh = False, color=(.7,.7,.7,1)) :

    ro = VolumeViewer.volume.Rendering_Options()
    ro.smoothing_factor = .3
    ro.smoothing_iterations = 2
    ro.surface_smoothing = False
    ro.square_mesh = True
    ro.line_thickness = 1

    dmap.update_surface ( False, ro )
    for sp in dmap.surfacePieces :
        v, t = sp.geometry
        if len(v) == 8 and len(t) == 12 :
            sp.display = False
        else :
            if showMesh :
                sp.color = (color[0]/2.0, color[1]/2.0, color[2]/2.0, 1.0)
                sp.displayStyle = sp.Mesh
            else :
                sp.color = (color[0], color[1], color[2], color[3])


def MapFromData ( ndata, nname, dmap, showMesh, thrF=1.0, color=(.7,.7,.7,1) ) :

    if showMesh :
        m = GetMod ( nname + "_mesh" )
        if m != None :
            chimera.openModels.close ( [m] )
    else :
        m = GetMod ( nname )
        if m != None :
            chimera.openModels.close ( [m] )


    nv = VolumeViewer.volume.volume_from_grid_data ( ndata )
    nv.openState.xform = dmap.openState.xform
    nv.name = nname
    if showMesh :
        nv.name = nname + "_mesh"
    nv.region = ( nv.region[0], nv.region[1], [1,1,1] )
    nv.surface_levels[0] = dmap.surface_levels[0] * thrF

    MapUp(nv, showMesh, color)
    return nv



def diha ( a1, a2, a3, a4 ) :
    #n1 = vnorm ( a1.coord(), a2.coord(), a3.coord() )
    #n2 = vnorm ( a2.coord(), a3.coord(), a4.coord() )
    #return numpy.arccos ( n2 * n1 * -1.0 ) * 180.0 / numpy.pi

    # http://math.stackexchange.com/questions/47059/how-do-i-calculate-a-dihedral-angle-given-cartesian-coordinates
    b1 = a2.coord() - a1.coord()
    b2 = a3.coord() - a2.coord()
    b3 = a4.coord() - a3.coord()

    n1 = chimera.cross ( b1, b2 ); n1.normalize()
    n2 = chimera.cross ( b2, b3 ); n2.normalize()
    m1 = chimera.cross ( n1, b2 ); m1.normalize()

    x = n1 * n2
    y = m1 * n2

    return -1.0 * numpy.arctan2 ( y, x) * 180.0 / numpy.pi


def angle ( a1, a2, a3 ) :
    n1 = a1.coord() - a2.coord()
    n2 = a3.coord() - a2.coord()
    return numpy.arccos ( (n2/n1.length) * (n1/n2.length) )  * 180.0 / numpy.pi


class Bone (object) :

    def __init__ (self, a1, a2, a3) :
        BoneInit ( self, a1, a2, a3 )

    def CS ( self ) :
        return CS ( a1.coord(), a2.coord(), a3.coord() )

    def CS0 ( self ) :
        return CS ( a1.coord0, a2.coord0, a3.coord0 )

    def Xf ( self ) :
        X,Y,Z = CS ( self.a1.coord(), self.a2.coord(), self.a3.coord() )
        return chimera.Xform.coordFrame ( X, Y, Z, self.a1.coord(), True )

    def Xf0 ( self ) :
        X,Y,Z = CS ( self.a1.coord0, self.a2.coord0, self.a3.coord0 )
        return chimera.Xform.coordFrame ( X, Y, Z, self.a1.coord0, True )

    def MakeFrame ( self ) :
        BoneMakeFrame ( self )

    def DistToPoint ( self, pt ) :
        return BoneDistToPoint ( self, pt )

    def SkinPoint ( self, pt ) :
        return BoneSkinPoint ( self, pt )


def BoneInit (bo, a1, a2, a3) :
    bo.a1, bo.a2, bo.a3 = a1, a2, a3
    bo.X0, bo.Y0, bo.Z0 = CS ( a1.coord0, a2.coord0, a3.coord0 )
    bo.F0 = chimera.Xform.coordFrame ( bo.X0, bo.Y0, bo.Z0, bo.a1.coord0, True )

def BoneMakeFrame ( bo ) :
    bo.X, bo.Y, bo.Z = CS ( bo.a1.coord(), bo.a2.coord(), bo.a3.coord() )
    bo.F = chimera.Xform.coordFrame ( bo.X, bo.Y, bo.Z, bo.a1.coord(), True )
    bo.F = bo.F.inverse()


def CS ( p1, p2, p3 ) :
    X = p2 - p1; X.normalize()
    Y = p3 - p1; Y.normalize()
    Z = chimera.cross ( X, Y ); Z.normalize()
    Y = chimera.cross ( Z, X ); Y.normalize()
    return X,Y,Z


def BoneDistToPoint ( bo, pt ) :

    pt = chimera.Point(pt[0], pt[1], pt[2])
    V = bo.a2.coord() - bo.a1.coord()
    v = pt - bo.a1.coord()
    t = V * v
    if t < 0.0 :
        return v.length
    elif t > 1.0 :
        return (pt-bo.a2.coord()).length
    else :
        lp = bo.a1.coord() + (V*t)
        return (pt-lp).length


def BoneSkinPoint ( bo, pt ) :

    #bo.X, bo.Y, bo.Z = CS ( bo.a1.coord(), bo.a2.coord(), bo.a3.coord() )
    #x = chimera.Xform.coordFrame ( bo.X, bo.Y, bo.Z, bo.a1.coord(), True )
    #x = x.inverse()
    #y = chimera.Xform.coordFrame ( bo.X0, bo.Y0, bo.Z0, bo.a1.coord0, True )

    pt = chimera.Point ( pt[0], pt[1], pt[2] )
    pt = bo.F.apply ( pt )
    pt = bo.F0.apply ( pt )
    return [pt[0], pt[1], pt[2]]





# -----------------------------------------------------------------------------



# -----------------------------------------------------------------------------
#


def FitMolToMap ( mol, dmap, RES, doTranslate = True, doRotate = True ) :

    import FitMap

    fpoints = _multiscale.get_atom_coordinates ( mol.atoms, transformed = True )
    fpoint_weights = numpy.ones ( len(mol.atoms), numpy.float32 )

    darray = dmap.data.matrix()

    xyz_to_ijk_tf = dmap.data.xyz_to_ijk_transform

    dmm = Matrix.invert_matrix ( Matrix.xform_matrix ( dmap.openState.xform ) )

    mm = Matrix.multiply_matrices ( dmap.data.xyz_to_ijk_transform, dmm )


    map_values, outside = VolumeData.interpolate_volume_data(fpoints, mm, darray)

    #olap0, cc0, other = FitMap.overlap_and_correlation ( fpoint_weights, map_values )
    avg0 = numpy.average ( map_values )
    #print " - 0 - ", avg0,

    move_tf, stats = FitMap.locate_maximum(fpoints, fpoint_weights,
                                    darray, mm,
                                    max_steps = 1000,
                                    ijk_step_size_min = 0.01,
                                    ijk_step_size_max = 0.5,
                                    optimize_translation = doTranslate,
                                    optimize_rotation = doRotate,
                                    metric = 'sum product',
                                    request_stop_cb = None)

    xf = chimera_xform ( move_tf )
    avg1 = stats['average map value']
    #print " - 1 - ", avg1

    xfm = mol.openState.xform
    xfm.premultiply ( xf )
    mol.openState.xform = xfm

    molg = MyMolMapX ( mol, mol.atoms, RES, dmap.data.step[0], chimera.Xform.identity() )
    fpoints, fpoint_weights = fit_points_g ( molg, 0.22 )
    map_values = dmap.interpolated_values ( fpoints, mol.openState.xform )
    import FitMap
    mmolap, cc, ccm = FitMap.overlap_and_correlation ( fpoint_weights, map_values )
    #print "Molmap - olap: %f, CC: %f, CCm: %f" % (mmolap, mmcorr1, mmcorr2)

    return avg1, cc, ccm




# ---------------------------------------------------

def getdialog ( create=False ) :

    from chimera import dialogs
    d = dialogs.find ( dlgName, create=False )
    return d



def close_dialog () :
    from chimera import dialogs


def setro (ro) :
    from chimera import dialogs
    d = dialogs.find ( "volume viewer", create=False )
    if d :
        d.surface_options_panel.set_gui_from_rendering_options (ro)
        #d.redisplay_needed_cb()


def vold () :
    from chimera import dialogs
    d = dialogs.find ( "volume viewer", create=False )
    d.surface_options_panel.line_thickness.set(2)
    d.redisplay_needed_cb()
    set_gui_from_rendering_options



def show_dialog () :

    from chimera import dialogs

    d = dialogs.find ( dlgName, create=False )
    if d :
        print " - found old diag"
        d.toplevel_widget.update_idletasks ()
        d.Close()
        d.toplevel_widget.update_idletasks ()

    dialogs.register (MapQ_Dialog.name, MapQ_Dialog, replace = True)

    d = dialogs.find ( dlgName, create=True )
    # Avoid transient dialog resizing when created and mapped for first time.
    d.toplevel_widget.update_idletasks ()
    d.enter()

    return d



def GetMod ( name ) :
    for m in chimera.openModels.list() :
        if m.name == name :
            return m
    return None




#def GetVisibleMol () :
#    for m in chimera.openModels.list() :
#        if m.display == True and type(m) == chimera.Molecule :
#            return m
#    return None

NA = {
    "A" : { "baseAtoms" : ["","",""] }
}


class NA ( object ):

    type
