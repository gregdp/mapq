
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
    import molref
    reload (molref)
except :
    pass

try :
    from segment_dialog import current_segmentation, segmentation_map
    import molbuild
    reload (molbuild)
except :
    pass

mapqVersion = ""
showDevTools = True

try :
    import Segger
    mapqVersion = Segger.mapqVersion
    showDevTools = Segger.showDevTools
    print "-segger"
except :
    mapqVersion = "1.6.4"
    showDevTools = False


import qscores
#reload (qscores)


gSigma = 0.6

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
        buttons = ( "SegMod", "Select", "Log", "Close" )
    else :
        buttons = ( "Log", "Close" )
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
        ff.grid(column=0, row=row, sticky='w', pady=0, padx=5)

        if 1 :
            ff = Tkinter.Frame(f)
            ff.grid(column=0, row=row, sticky='w', pady=5, padx=0)

            l = Tkinter.Label(ff, text=' Map:', anchor=Tkinter.W)
            l.grid(column=0, row=0, sticky='w')

            self.dmap = Tkinter.StringVar(parent)
            self.dmapMB  = Tkinter.Menubutton ( ff, textvariable=self.dmap, relief=Tkinter.RAISED, width=20 )
            self.dmapMB.grid (column=1, row=0, sticky='we', padx=1)
            self.dmapMB.menu  =  Tkinter.Menu ( self.dmapMB, tearoff=0, postcommand=self.MapMenu )
            self.dmapMB["menu"]  =  self.dmapMB.menu

            self.cur_dmap = None
            self.SetVisMap ()


            l = Tkinter.Label(ff, text='Model:', anchor=Tkinter.W)
            l.grid(column=2, row=0, sticky='w')

            self.struc = Tkinter.StringVar(parent)
            self.strucMB  = Tkinter.Menubutton ( ff, textvariable=self.struc, relief=Tkinter.RAISED, width=20 )
            self.strucMB.grid (column=3, row=0, sticky='we', padx=1)
            self.strucMB.menu  =  Tkinter.Menu ( self.strucMB, tearoff=0, postcommand=self.StrucMenu )
            self.strucMB["menu"]  =  self.strucMB.menu

            self.cur_mol = None
            self.cur_chains = []
            self.SetVisMol ()


            l = Tkinter.Label(ff, text=" Chain:" )
            l.grid(column=4, row=0, sticky='w')

            self.chain = Tkinter.StringVar(parent)
            self.chainMB  = Tkinter.Menubutton ( ff, textvariable=self.chain, relief=Tkinter.RAISED, width=4 )
            self.chainMB.grid (column=5, row=0, sticky='we', padx=1)
            self.chainMB.menu  =  Tkinter.Menu ( self.chainMB, tearoff=0, postcommand=self.ChainMenu )
            self.chainMB["menu"]  =  self.chainMB.menu

            l = Tkinter.Label(ff, text=" Show:" )
            l.grid(column=6, row=0, sticky='w')


            b = Tkinter.Button(ff, text="Chain", command=self.AllChain)
            b.grid (column=7, row=0, sticky='w', padx=1)

            b = Tkinter.Button(ff, text="All", command=self.AllChains)
            b.grid (column=8, row=0, sticky='w', padx=1)

            b = Tkinter.Button(ff, text="Sel.", command=self.ShowOnlySel)
            b.grid (column=9, row=0, sticky='w', padx=1)

            b = Tkinter.Button(ff, text="At.", command=self.SetSelAtoms)
            b.grid (column=10, row=0, sticky='w', padx=1)

            b = Tkinter.Button(ff, text="Rib.", command=self.SetSelRibbon)
            b.grid (column=11, row=0, sticky='w', padx=1)

            b = Tkinter.Button(ff, text="SCs", command=self.ShowSCs)
            b.grid (column=12, row=0, sticky='w', padx=1)

            b = Tkinter.Button(ff, text="~SCs", command=self.HideSCs)
            b.grid (column=13, row=0, sticky='w', padx=1)

            b = Tkinter.Button(ff, text="W", command=self.Wire)
            b.grid (column=14, row=0, sticky='w', padx=1)






        if 1 :

            l = Tkinter.Label(ff, text='   Zoom:', fg="#777")
            l.grid(column=35, row=0, sticky='e')

            b = Tkinter.Button(ff, text="-", command=self.ZoomMinus)
            b.grid (column=36, row=0, sticky='w', padx=0)

            b = Tkinter.Button(ff, text="+", command=self.ZoomPlus)
            b.grid (column=37, row=0, sticky='w', padx=0)

            b = Tkinter.Button(ff, text="<", command=self.ZoomBegin)
            b.grid (column=38, row=0, sticky='w', padx=0)

            b = Tkinter.Button(ff, text=">", command=self.ZoomEnd)
            b.grid (column=39, row=0, sticky='w', padx=0)




        self.mapRes = Tkinter.StringVar(f)
        self.mapRes.set ( "3" )


        if 1 :

            row += 1
            ff = Tkinter.Frame(f)
            ff.grid(column=0, row=row, sticky='w', pady=0, padx=0)

            b = Tkinter.Label(ff, text=" Res:")
            b.grid (column=0, row=0, sticky='w', padx=0, pady=5)

            e = Tkinter.Entry(ff, width=3, textvariable=self.mapRes)
            e.grid(column=1, row=0, sticky='w', padx=0, pady=5)

            #fff = Tkinter.Frame(ff, borderwidth=1, padx=2, pady=2, relief=Tkinter.GROOVE)
            #fff.grid(column=10, row=0, sticky='e', pady=0, padx=5)

            l = Tkinter.Label(ff, text=' Q-scores:', anchor=Tkinter.W, font = 'TkCaptionFont')
            l.grid(column=2, row=0, sticky='w')


            #b = Tkinter.Button(fff, text="Sigma", command=self.CalcAllSigma )
            #b.grid (column=2, row=0, sticky='w', padx=5)

            #b = Tkinter.Button(fff, text="RadZ", command=self.CalcAllRadZ )
            #b.grid (column=4, row=0, sticky='w', padx=5)

            if isModelZ :

                b = Tkinter.Button(ff, text="Z-scores", command=self.CalcZScores )
                b.grid (column=5, row=0, sticky='w', padx=5)

            else :

                b = Tkinter.Button(ff, text="Calc", command=self.CalcAllQ )
                b.grid (column=5, row=0, sticky='w', padx=5)

                b = Tkinter.Button(ff, text="Calc(P)", command=self.CalcAllQp )
                b.grid (column=6, row=0, sticky='w', padx=5)

                b = Tkinter.Button(ff, text="Load", command=self.GetQsFromFile )
                b.grid (column=7, row=0, sticky='w', padx=5)

                b = Tkinter.Button(ff, text="Sel", command=self.CalcSelQ )
                b.grid (column=8, row=0, sticky='w', padx=5)

                #b = Tkinter.Button(fff, text="R", command=self.CalcAllR )
                #b.grid (column=5, row=0, sticky='w', padx=5)

                #b = Tkinter.Button(fff, text="R", command=self.CalcAllR )
                #b.grid (column=5, row=0, sticky='w', padx=5)




            if 0 :

                self.colorMod = Tkinter.StringVar()
                self.colorMod.set ( 'sc' )

                b = Tkinter.Button(ff, text="Color:", command=self.DoColor)
                b.grid (column=20, row=0, sticky='w', padx=5)

                c = Tkinter.Radiobutton(ff, text="Bb", variable=self.colorMod, value = 'bb')
                c.grid (column=21, row=0, sticky='w')

                c = Tkinter.Radiobutton(ff, text="Sc", variable=self.colorMod, value = 'sc')
                c.grid (column=22, row=0, sticky='w')

                c = Tkinter.Radiobutton(ff, text="Rand", variable=self.colorMod, value = 'rand')
                c.grid (column=23, row=0, sticky='w')


            else :

                l = Tkinter.Label(ff, text=' On:', fg="#777")
                l.grid(column=20, row=0, sticky='e')


                b = Tkinter.Button(ff, text="Bb", command=self.DoColorBB)
                b.grid (column=21, row=0, sticky='w', padx=5)

                b = Tkinter.Button(ff, text="Sc", command=self.DoColorSC)
                b.grid (column=22, row=0, sticky='w', padx=5)

                b = Tkinter.Button(ff, text="Res", command=self.DoColorRes)
                b.grid (column=23, row=0, sticky='w', padx=5)

                if not isModelZ :
                    b = Tkinter.Button(ff, text="At", command=self.DoColorAtoms)
                    b.grid (column=24, row=0, sticky='w', padx=5)

                b = Tkinter.Button(ff, text="Rand", command=self.DoColorRandom)
                b.grid (column=25, row=0, sticky='w', padx=5)





            l = Tkinter.Label(ff, text='', fg="#000")
            l.grid(column=25, row=0, sticky='ens')


            #ff = Tkinter.Frame(ff, borderwidth=1, padx=2, pady=2, relief=Tkinter.GROOVE)
            #ff.grid(column=30, row=0, sticky='e', pady=0, padx=5)

            l = Tkinter.Label(ff, text='Select: ', fg="#000", font = 'TkCaptionFont')
            l.grid(column=35, row=0, sticky='ens')

            #oft = Hybrid.Checkbutton(ff, 'Ribbon', True)
            #oft.button.grid(column = 36, row = 0, sticky = 'w')
            #self.showRibbon = oft.variable
            #self.showRibbon.set ( 1 )

            #oft = Hybrid.Checkbutton(ff, 'Side Chains', True)
            #oft.button.grid(column = 37, row = 0, sticky = 'w')
            #self.showAtoms = oft.variable
            #self.showRibbon.set ( 1 )

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
            oft = Tkinter.Checkbutton( ff, text="Vol", variable=self.preserveVol, command=self.preserveVolCb)
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

            b = Tkinter.Button(ff, text="!", command=self.SelReLoad)
            b.grid (column=46, row=0, sticky='w', padx=0)

            if 0 and showDevTools :

                b = Tkinter.Button(ff, text="L", command=self.SelLoad)
                b.grid (column=47, row=0, sticky='w', padx=5)



            #b = Tkinter.Button(ff, text="Clear", command=self.ClearSel)
            #b.grid (column=40, row=0, sticky='w', padx=5)

            #self.keepExMap = Tkinter.IntVar()
            #self.keepExMap.set(0)
            #oft = Tkinter.Checkbutton( ff, text="Keep Extracted Maps", variable=self.keepExMap, command=self.keepExMapCb)
            #oft.grid(column = 40, row = 0, sticky = 'w')



        # ----------- select panel ----------------------------------

        if 1 :

            row += 1
            op = Hybrid.Popup_Panel(f)
            ff = op.frame
            ff.grid(row = row, column = 0, sticky = 'news')
            ff.grid_remove()
            #ff.columnconfigure(0, weight=1)
            self.selPanel = op.panel_shown_variable

            #ff = Tkinter.Frame(f)
            #ff.grid(column=0, row=row, sticky='w', pady=0, padx=5)

            l = Tkinter.Label(ff, text=' Sel:', font = 'TkCaptionFont')
            l.grid(column=1, row=0, sticky='w', pady=5)


            if 0 :
                #b = Tkinter.Button(ff, text="Asp", command=self.asp )
                #b.grid (column=1, row=0, sticky='w', padx=5)

                b = Tkinter.Button(ff, text="Extr", command=self.Extract )
                b.grid (column=2, row=0, sticky='w', padx=5)

                b = Tkinter.Button(ff, text="Al 1", command=self.AlignRes1 )
                b.grid (column=3, row=0, sticky='w', padx=5)

                b = Tkinter.Button(ff, text="Al 2", command=self.AlignRes2 )
                b.grid (column=4, row=0, sticky='w', padx=5)

                b = Tkinter.Button(ff, text="Avg", command=self.Avg )
                b.grid (column=5, row=0, sticky='w', padx=5)

                b = Tkinter.Button(ff, text="~Extr", command=self.CloseExtracted )
                b.grid (column=6, row=0, sticky='w', padx=5)


                #b = Tkinter.Button(ff, text="Sbb", command=self.BB_Sigma )
                #b.grid (column=8, row=0, sticky='w', padx=5)

                #b = Tkinter.Button(ff, text="Z", command=self.ZScoreSel )
                #b.grid (column=9, row=0, sticky='w', padx=5)

                #b = Tkinter.Button(ff, text="Zr", command=self.RotaZ1 )
                #b.grid (column=10, row=0, sticky='w', padx=5)

                #b = Tkinter.Button(ff, text="R1", command=self.R1 )
                #b.grid (column=11, row=0, sticky='w', padx=5)

                #b = Tkinter.Button(ff, text="ExA", command=self.ExCustA )
                #b.grid (column=12, row=0, sticky='w', padx=5)

                #b = Tkinter.Button(ff, text="ExB", command=self.ExCustB )
                #b.grid (column=13, row=0, sticky='w', padx=5)

                #b = Tkinter.Button(ff, text="ExC", command=self.ExCustC )
                #b.grid (column=14, row=0, sticky='w', padx=5)


            b = Tkinter.Button(ff, text="S-sel", command=self.S_sel )
            b.grid (column=20, row=0, sticky='w', padx=5)

            b = Tkinter.Button(ff, text="Q-sel", command=self.Q_sel )
            b.grid (column=21, row=0, sticky='w', padx=5)

            if 0 :
                b = Tkinter.Button(ff, text="Q-show", command=self.Q_show )
                b.grid (column=22, row=0, sticky='w', padx=5)

                b = Tkinter.Button(ff, text="SA-Q", command=self.SA_Q )
                b.grid (column=23, row=0, sticky='w', padx=5)


            #b = Tkinter.Button(ff, text="Ats", command=self.ShowAts)
            #b.grid (column=25, row=0, sticky='w', padx=10)

            if 1 :
                b = Tkinter.Button(ff, text="Alts", command=self.FindAlts)
                b.grid (column=28, row=0, sticky='w', padx=5)

                b = Tkinter.Button(ff, text="X-Alts", command=self.DelAlts)
                b.grid (column=29, row=0, sticky='w', padx=5)

            if 0 :
                b = Tkinter.Button(ff, text="APro", command=self.AProfs)
                b.grid (column=28, row=0, sticky='w', padx=5)

            #b = Tkinter.Button(ff, text="Ligs", command=self.Ligs)
            #b.grid (column=43, row=0, sticky='w', padx=5)

            #b = Tkinter.Button(ff, text="Scale", command=self.Scale)
            #b.grid (column=44, row=0, sticky='w', padx=5)



            b = Tkinter.Label(ff, text="   Str:")
            b.grid (column=30, row=0, sticky='w', padx=0, pady=5)

            self.selText = Tkinter.StringVar(f)
            self.selText.set ( "" )
            e = Tkinter.Entry(ff, width=20, textvariable=self.selText)
            e.grid(column=31, row=0, sticky='w', padx=5, pady=5)


            b = Tkinter.Button(ff, text="Sel", command=self.SelText)
            b.grid (column=32, row=0, sticky='w', padx=5)



            b = Tkinter.Label(ff, text="Rad:")
            b.grid (column=33, row=0, sticky='w', padx=0, pady=5)

            self.maskRad = Tkinter.StringVar(f)
            self.maskRad.set ( "2.5" )
            e = Tkinter.Entry(ff, width=3, textvariable=self.maskRad)
            e.grid(column=34, row=0, sticky='w', padx=5, pady=5)


            b = Tkinter.Button(ff, text="AddSel", command=self.AdSel)
            b.grid (column=35, row=0, sticky='w', padx=5)


            b = Tkinter.Button(ff, text="Nr", command=self.ShowNear)
            b.grid (column=40, row=0, sticky='w', padx=5)

            b = Tkinter.Button(ff, text="Ds", command=self.ShowDists)
            b.grid (column=41, row=0, sticky='w', padx=5)

            b = Tkinter.Button(ff, text="Inter", command=self.Inter)
            b.grid (column=42, row=0, sticky='w', padx=5)


            b = Tkinter.Button(ff, text="Occ", command=self.Occ)
            b.grid (column=43, row=0, sticky='w', padx=5)


            b = Tkinter.Button(ff, text="Rmsd", command=self.RMSD)
            b.grid (column=44, row=0, sticky='w', padx=5)



        if 0 :

            row += 1
            ff = Tkinter.Frame(f)
            ff.grid(column=0, row=row, sticky='w', pady=0, padx=5)


            b = Tkinter.Label(ff, text="   Atom:")
            b.grid (column=15, row=0, sticky='w', padx=0, pady=5)

            self.addText = Tkinter.StringVar(f)
            self.addText.set ( "Ca" )
            e = Tkinter.Entry(ff, width=10, textvariable=self.addText)
            e.grid(column=16, row=0, sticky='w', padx=5, pady=5)


            b = Tkinter.Button(ff, text="Add", command=self.AddAtom)
            b.grid (column=17, row=0, sticky='w', padx=5)





        if 1 :
            row += 1
            op = Hybrid.Popup_Panel(f)
            ff = op.frame
            ff.grid(row = row, column = 0, sticky = 'news')
            ff.grid_remove()
            #ff.columnconfigure(0, weight=1)
            self.modPanel = op.panel_shown_variable


            oft = Hybrid.Checkbutton(ff, 'Gaps', True)
            oft.button.grid(column = 2, row = 0, sticky = 'w')
            self.showGaps = oft.variable
            #self.showRibbon.set ( 1 )


            b = Tkinter.Label(ff, text="   Add:")
            b.grid (column=6, row=0, sticky='w', padx=0, pady=5)

            self.addRess = Tkinter.StringVar(f)
            #self.addRess.set ( "vsgtngtkrf" )
            self.addRess.set ( "NAG" )
            e = Tkinter.Entry(ff, width=30, textvariable=self.addRess)
            e.grid(column=7, row=0, sticky='w', padx=5, pady=5)


            b = Tkinter.Button(ff, text="Add", command=self.AddRes)
            b.grid (column=11, row=0, sticky='w', padx=5)

            b = Tkinter.Button(ff, text="N-", command=self.AddResN)
            b.grid (column=12, row=0, sticky='w', padx=5)

            b = Tkinter.Button(ff, text="C-", command=self.AddResC)
            b.grid (column=13, row=0, sticky='w', padx=5)

            b = Tkinter.Button(ff, text="Ref", command=self.Refine)
            b.grid (column=14, row=0, sticky='w', padx=5)

            b = Tkinter.Button(ff, text="Del", command=self.DelSel)
            b.grid (column=15, row=0, sticky='w', padx=5)

            b = Tkinter.Button(ff, text="Take", command=self.Take)
            b.grid (column=16, row=0, sticky='w', padx=5)

            b = Tkinter.Button(ff, text="D", command=self.DMS)
            b.grid (column=17, row=0, sticky='w', padx=5)

            b = Tkinter.Button(ff, text="S", command=self.SS)
            b.grid (column=18, row=0, sticky='w', padx=5)


        dummyFrame = Tkinter.Frame(parent, relief='groove', borderwidth=1)
        Tkinter.Frame(dummyFrame).pack()
        dummyFrame.grid(row=row,column=0,columnspan=7, pady=3, sticky='we')
        row += 1


        global msg
        msg = Tkinter.Label(parent, width = 60, anchor = 'w', justify = 'left', fg="red", pady=5, padx=10)
        msg.grid(column=0, row=row, sticky='ew')
        self.msg = msg

        self.showingAtoms = False


        if len ( self.cur_chains ) > 0 :
            self.chain.set ( self.cur_chains[0] )
            #self.ShowCh ( self.cur_chains[0] )
            self.GetSeq ()

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


        self.modPanel.set(showDevTools)
        self.selPanel.set(showDevTools)




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






    def Select ( self ) :
        self.selPanel.set (not self.selPanel.get())

    def SegMod ( self ) :
        self.modPanel.set (not self.modPanel.get())

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
            self.dmap.set ( dmap.name + " (%d)" % dmap.id )
            self.cur_dmap = dmap


    def MapMenu ( self ) :
        #print "Map menu..."
        self.dmapMB.menu.delete ( 0, 'end' )   # Clear menu
        self.cur_dmap = None
        self.dmap.set("")
        mlist = OML(modelTypes = [VolumeViewer.volume.Volume])
        for m in mlist :
            self.dmapMB.menu.add_radiobutton ( label=m.name+" (%d)"%m.id, variable=self.dmap,
                                command=lambda m=m: self.MapSelected(m) )


    def MapSelected ( self, dmap ) :

        self.cur_dmap = dmap
        print "Selected " + dmap.name

        self.GetSeq ()
        self.ZoomBegin ()


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
            self.struc.set ( mol.name + " (%d)" % mol.id )
            self.cur_mol = mol
            self.cur_chains = self.GetChains ( mol )
            SetBBAts ( mol )


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
                self.chain.set ( "" )
            elif self.chain.get() in self.cur_chains :
                print " - ch " + self.chain.get() + " already sel"
                #self.ShowCh ( self.chain.get() )
            else :
                self.chain.set ( self.cur_chains[0] )
                #self.ShowCh ( self.chain.get() )


        SetBBAts ( mol )
        self.parent.after(100, self.DoSeq)


    def DoSeq ( self ) :
        print "after 100"
        self.GetSeq ()
        self.ZoomBegin ()


    def ChainSelected ( self, ch ) :
        print " - sel chain: ", ch, self.chain.get()
        #self.ShowCh ( ch )
        self.parent.after(100, self.DoSeq)



    def StrucMenu ( self ) :
        self.strucMB.menu.delete ( 0, 'end' )   # Clear menu
        mlist = OML(modelTypes = [chimera.Molecule])
        for m in mlist :
            self.strucMB.menu.add_radiobutton ( label=m.name+" (%d)"%m.id, variable=self.struc,
                                           command=lambda m=m: self.StrucSelected(m) )

    def ChainMenu ( self ) :
        self.chainMB.menu.delete ( 0, 'end' )   # Clear menu
        print " - chain menu"
        print self.cur_chains
        for ch in self.cur_chains :
            self.chainMB.menu.add_radiobutton ( label=ch, variable=self.chain, command=lambda ch=ch: self.ChainSelected(ch) )

        self.chainMB.menu.add_radiobutton ( label="All", variable=self.chain, command=lambda ch="All": self.ChainSelected("All") )



    def DoColor ( self ) :

        print "color...", self.colorMod.get()

        #colSC = self.colorSC.get()
        #colRand = self.colorRand.get()

        if self.colorMod.get() == "rand" :
            self.RandColorChains()
        else :
            self.UpdateModColor ()

        #if self.colorMap.get() :
        #    self.UpdateSurfColor ()



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
            umsg ( "No molecule/chain selected?" )
            return

        if not hasattr (self, 'scores') :
            umsg ( "No scores... press Q, Qp, or Qf button first" )
            return

        foundScore = False
        for ri, r in enumerate ( self.seqRes ) :
            if r != None  and hasattr (r, 'Q') :
                foundScore = True

        if not foundScore :
            umsg ( "No scores... press Calc button first" )
            return


        minScore, maxScore = 0,0
        if colorMod == "sc" :
            minScore, maxScore = self.minScore1, self.maxScore1
        else :
            minScore, maxScore = self.minScore2, self.maxScore2

        cH = numpy.array( [0.0,1.0,0.0] )
        cL = numpy.array( [1.0,0.0,0.0] )

        for ri, r in enumerate ( self.seqRes ) :
            sc = None
            if r == None :
                continue
            #sc = self.scores[ri] if colorSC else self.scores2[ri]
            if colorMod == "sc" :
                if hasattr (r, 'scQ') :
                    sc = r.scQ
                else :
                    sc = 0
            elif colorMod == "bb" :
                sc = r.bbQ if hasattr (r, 'bbQ') else 0
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

        SetBBAts ( self.cur_mol )
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
            risel["%d.%s" % (r.id.position, r.id.chainId)] = 1

        for r in m.residues :
            rid = "%d.%s" % (r.id.position, r.id.chainId)
            if rid in risel :
                r.ribbonDisplay = not self.showingAtoms
                for at in r.atoms :
                    if at.element.name == "H" :
                        at.display = False
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

        atMap = {}
        for r in m.residues :

            altScores = {}
            for at in r.atoms :
                if at.isSC :
                    alt = "_" if at.altLoc == '' else at.altLoc
                    if alt in altScores :
                        altScores[alt].append ( at.Q )
                    else :
                        altScores[alt] = [at.Q]

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
        else :
            self.SetDrawMode ( self.GetCurRess(), showRibbon = False )

        self.showingAtoms = True




    def SetDrawMode ( self, ress, showRibbon = None ) :

        #if showRibbon == None :
        #    showRibbon = segmod_dialog().showRibbon.get()

        #showRibbon = True

        #SetBBAts ( ress[0].molecule )

        for m in chimera.openModels.list() :
            if type(m) == chimera.Molecule :
                if not hasattr ( m, 'bbats' ) :
                    SetBBAts(m)
                    m.bbats = True


        atMap = {}
        #atI = 0
        #c1 = (1.0,0.0,0.0,1)
        #c1 = (1.0,0.0,0.0,1)
        for res in ress :
            for at in res.atoms :

                if not hasattr (res, 'isProt') :
                    SetBBAts (res.molecule)

                if res.isProt or res.isNA :
                    at.drawMode = at.EndCap
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
            if not hasattr ( mol, 'bbats' ) :
                SetBBAts(mol)
                mol.bbats = True

        ress = chimera.selection.currentResidues()
        if len(ress) == 0 :
            ress = self.GetCurRess()

        for res in ress :
            #if res.id.position in rs and res.id.chainId == cid :
            for at in res.atoms :
                #at.drawMode = at.EndCap

                if not hasattr (at, 'isBB') :
                    SetBBAts (at.molecule)

                at.display = at.isBB
                if at.residue.isNA :
                    at.display = at.isBB and not at.isSugar

                #try :
                #    at.color = atomColors[at.element.name.upper()]
                #except :
                #    at.color = atomColors[" "]


    def ShowSCs ( self ) :

        for mol in chimera.selection.currentMolecules() :
            if not hasattr ( mol, 'bbats' ) :
                SetBBAts(mol)
                mol.bbats = True

        ress = chimera.selection.currentResidues()
        if len(ress) == 0 :
            ress = self.GetCurRess()

        for res in ress :
            for at in res.atoms :
                #at.drawMode = at.EndCap
                at.display = True
                try :
                    at.color = atomColors[at.element.name.upper()]
                except :
                    at.color = atomColors[" "]


    def ShowNear ( self ) :

        for mol in chimera.selection.currentMolecules() :
            if not hasattr ( mol, 'bbats' ) :
                SetBBAts(mol)
                mol.bbats = True

        ress = chimera.selection.currentResidues()
        if len(ress) == 0 :
            ress = self.GetCurRess()

        print "Near %d res:" % len(ress)
        for r in ress :
            print "%s.%d.%s - %d atoms" % (r.type, r.id.position, r.id.chainId, len(r.atoms))



        ats = [at for at in self.cur_mol.atoms if not at.element.name == "H"]
        points = _multiscale.get_atom_coordinates ( ats, transformed = False )
        print " - search tree: %d/%d ats" % ( len(ats), len(self.cur_mol.atoms) )
        allAtTree = AdaptiveTree ( points.tolist(), ats, 2.0)


        chimera.selection.clearCurrent ()

        nearRes = {}
        for r in ress :
            nats = self.AtsWithin ( r.atoms, 6.0, allAtTree )
            for at in nats :
                nearRes[at.residue] = 1

        for r in nearRes.keys() :
            print " -- %s.%d.%s - %d atoms" % (r.type, r.id.position, r.id.chainId, len(r.atoms))
            #chimera.selection.mergeCurrent ( chimera.selection.EXTEND, chimera.selection.OSLSelection ("") )
            if r in ress :
                continue
            chimera.selection.addCurrent ( r )
            for at in r.atoms :
                #at.drawMode = at.EndCap
                at.display = True
                if at.element.name.upper() in atomColors :
                    at.color = atomColors[at.element.name.upper()]




    def Inter ( self ) :

        for mol in chimera.selection.currentMolecules() :
            if not hasattr ( mol, 'bbats' ) :
                SetBBAts(mol)
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

        SetBBAts ( self.cur_mol )

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
                elif r.type in nucleic3to1 :
                    self.seq = self.seq + nucleic3to1[r.type]
                    self.conf = self.conf + "9"
                    self.predi = "C"
                    self.pred = self.pred + self.predi
                    self.seqRes.append ( r )
                    self.seqRi.append ( ri )
            else :
                if self.showGaps.get() :
                    self.seq = self.seq + "."
                    self.conf = self.conf + "9"
                    self.pred = self.pred + "C"
                    self.seqRes.append ( None )
                    self.seqRi.append ( ri )





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

        status ( "Getting secondary structure elements..." )

        ok = True
        try :
            print self.cur_dmap.name
        except :
            status ( "Selected map not found; please choose another map" )
            self.dmap.set ("")
            ok = False

        try :
            print self.cur_mol.name
        except :
            status ( "Selected model not found; please choose another model" )
            self.struc.set ("")
            self.chain.set ("")
            self.RemoveSeq ()
            ok = False

        if not ok :
            return


        resolution = 3.0 * self.cur_dmap.data.step[0]
        #resolution = 3.0
        umsg ( "Calculating backbone Z-scores..." )

        zscores2 = []

        if 0 : # old
            sses = SSEs ( self.seqRes )
            #print " - ",len(sses),"sse for ", len(ress), "res"

            atI = 1

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
                    r.bbQ = zscore
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


        ok = True
        try :
            print self.cur_dmap.name
        except :
            status ( "Selected map not found; please choose another map" )
            self.dmap.set ("")
            ok = False

        try :
            print self.cur_mol.name
        except :
            status ( "Selected model not found; please choose another model" )
            self.struc.set ("")
            self.chain.set ("")
            self.RemoveSeq ()
            ok = False

        if not ok :
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


        ok = True
        try :
            print self.cur_dmap.name
        except :
            status ( "Selected map not found; please choose another map" )
            self.dmap.set ("")
            ok = False

        try :
            print self.cur_mol.name
        except :
            status ( "Selected model not found; please choose another model" )
            self.struc.set ("")
            self.chain.set ("")
            self.RemoveSeq ()
            ok = False

        if not ok :
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


        ok = True
        try :
            print self.cur_dmap.name
        except :
            status ( "Selected map not found; please choose another map" )
            self.dmap.set ("")
            ok = False

        try :
            print self.cur_mol.name
        except :
            status ( "Selected model not found; please choose another model" )
            self.struc.set ("")
            self.chain.set ("")
            self.RemoveSeq ()
            ok = False

        if not ok :
            return


        cid = self.chain.get()

        ats = [at for at in self.cur_mol.atoms if not at.element.name == "H"]
        points = _multiscale.get_atom_coordinates ( ats, transformed = False )
        print " - search tree: %d/%d ats" % ( len(ats), len(self.cur_mol.atoms) )
        allAtTree = AdaptiveTree ( points.tolist(), ats, 1.0)

        if 0 :
            for r in self.cur_mol.residues :
                if hasattr ( r, 'Q' ) : del r.Q
                if hasattr ( r, 'scQ' ) : del r.scQ
                if hasattr ( r, 'bbQ' ) : del r.bbQ


        umsg ( "Calculating Q-scores - see bottom of main window for status or to cancel..." )

        Qavg = qscores.CalcQ (self.cur_mol, self.chain.get(), self.cur_dmap, gSigma, allAtTree=allAtTree, log=True )
        qscores.SaveQStats ( self.cur_mol, self.chain.get(), self.cur_dmap, gSigma, float(self.mapRes.get()) )
        self.ShowQScores ()

        #umsg ( "Average Q-score for %s: %.2f" % (self.cur_mol.name, Qavg) )
        umsg ( "Done Q-scores for %s" % (self.cur_mol.name) )



    def CalcAllQp (self) :

        ok = True
        try :
            print self.cur_dmap.name
        except :
            status ( "Selected map not found; please choose another map" )
            self.dmap.set ("")
            ok = False

        try :
            print self.cur_mol.name
        except :
            status ( "Selected model not found; please choose another model" )
            self.struc.set ("")
            self.chain.set ("")
            self.RemoveSeq ()
            ok = False

        if not ok :
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
                if hasattr ( r, 'scQ' ) : del r.scQ
                if hasattr ( r, 'bbQ' ) : del r.bbQ



        qscores.CalcQp (self.cur_mol, cid, self.cur_dmap, gSigma, allAtTree=None )
        qscores.SaveQStats ( self.cur_mol, self.chain.get(), self.cur_dmap, gSigma, float(self.mapRes.get()) )

        self.ShowQScores ()




    def ShowQScores (self) :

        cid = self.chain.get()

        scBB, scSC = [], []

        for r in self.cur_mol.residues :
            if cid == None or cid == "All" or r.id.chainId == cid :
                if r.isProt or r.isNA :
                    r.score1 = r.scQ
                    r.score2 = r.bbQ
                    if r.bbQ != None : scBB.append ( r.bbQ )
                    if r.scQ != None : scSC.append ( r.scQ )
                else :
                    r.score1 = r.Q
                    r.score2 = r.Q


        #bbRes = numpy.power ( numpy.e, (self.avgScore2 - 8.0334) / -4.128 ) # y = -4.128ln(x) + 8.0334
        #scRes = numpy.power ( numpy.e, (self.avgScore - 4.8261) / -3.097 ) # y = -3.097ln(x) + 4.8261
        #scRes = (self.avgScore2 - 3.507) / -0.721
        #bbRes = (self.avgScore - 6.1234) / -0.9191


        try :
            scMin, scMax, scAvg = min(scSC), max(scSC), numpy.average(scSC)
            bbMin, bbMax, bbAvg = min(scBB), max(scBB), numpy.average(scBB)


            print "Average Q sc : %.2f - %.2f, avg %.2f" % (scMin, scMax, scAvg )
            print "Average Q bb : %.2f - %.2f, avg %.2f" % (bbMin, bbMax, bbAvg )

            self.GetMaxScores()

        except :
            pass

        self.UpdateSeq ()



    def QuickQ (self) :

        ok = True
        try :
            print self.cur_dmap.name
        except :
            status ( "Selected map not found; please choose another map" )
            self.dmap.set ("")
            ok = False

        try :
            print self.cur_mol.name
        except :
            status ( "Selected model not found; please choose another model" )
            self.struc.set ("")
            self.chain.set ("")
            self.RemoveSeq ()
            ok = False

        if not ok :
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
                if hasattr ( r, 'scQ' ) : del r.scQ
                if hasattr ( r, 'bbQ' ) : del r.bbQ



        CalcQp (self.cur_mol, cid, self.cur_dmap, gSigma, allAtTree=None )



        scBB, scSC = [], []

        for r in self.cur_mol.residues :
            if cid == None or r.id.chainId == cid :
                if r.isProt or r.isNA :
                    r.score1 = r.scQ
                    r.score2 = r.bbQ
                    if r.bbQ != None : scBB.append ( r.bbQ )
                    if r.scQ != None : scSC.append ( r.scQ )
                else :
                    r.score1 = r.Q
                    r.score2 = r.Q

        scMin, scMax, scAvg = min(scSC), max(scSC), numpy.average(scSC)
        bbMin, bbMax, bbAvg = min(scBB), max(scBB), numpy.average(scBB)


        print " - Average Q sc : %.2f - %.2f, avg %.2f" % (scMin, scMax, scAvg )
        print " - Average Q bb : %.2f - %.2f, avg %.2f" % (bbMin, bbMax, bbAvg )


        self.minScore1, self.maxScore1 = 0.0,1
        self.minScore2, self.maxScore2 = 0.0,1


        self.UpdateSeq ()
        qscores.SaveQStats ( self.cur_mol, self.chain.get(), self.cur_dmap, gSigma, float(self.mapRes.get()) )






    def GetQsFromFile (self) :

        ok = True
        try :
            print self.cur_dmap.name
        except :
            status ( "Selected map not found; please choose another map" )
            self.dmap.set ("")
            ok = False

        try :
            print self.cur_mol.name
        except :
            status ( "Selected model not found; please choose another model" )
            self.struc.set ("")
            self.chain.set ("")
            self.RemoveSeq ()
            ok = False

        if not ok :
            return


        chainId = self.chain.get()


        if 0 :
            for r in self.cur_mol.residues :
                if hasattr ( r, 'Q' ) : del r.Q
                if hasattr ( r, 'scQ' ) : del r.scQ
                if hasattr ( r, 'bbQ' ) : del r.bbQ



        molPath = os.path.splitext(self.cur_mol.openedAs[0])[0]
        mapName = os.path.splitext(self.cur_dmap.name)[0]
        nname = molPath + "__Q__" + mapName + ".pdb"

        if not os.path.isfile ( nname ) :
            umsg ( "Q scores not found for this map and file - press Calc first" )
            return


        halfMap1, halfMap2 = "half_map_1" in mapName, "half_map_2" in mapName
        if halfMap1 : print " - half map 1"
        if halfMap2 : print " - half map 2"

        rids = {}
        for r in self.cur_mol.residues :
            rids["%d.%s" % (r.id.position,r.id.chainId)] = r


        if 0 :
            from _multiscale import get_atom_coordinates
            pts = get_atom_coordinates(atoms, transformed = False)
            A, B = maxD - minD, minD
            d_vals = dmap.interpolated_values ( pts, mol.openState.xform ).astype(numpy.float64, copy=False)

        minD, maxD = qscores.MinMaxD ( self.cur_dmap )

        M = self.cur_dmap.data.full_matrix()
        minD, maxD = numpy.min(M), numpy.max(M)

        maxM = numpy.max(M)
        minM = numpy.min(M)

        maxD = min ( numpy.average(M)+numpy.std(M)*10, maxM )
        minD = max ( numpy.average(M)-numpy.std(M)*1, minM )


        # http://www.wwpdb.org/documentation/file-format-content/format33/sect9.html#ATOM
        fin = open ( nname, "r" )
        for line in fin :
            if line[0:4] == "ATOM" or line[0:6] == "HETATM" :
                aname, aloc, cid, resi, occ, bfac = line[12:16].strip(), line[16:17].strip(), line[21], int(line[22:26]), float ( line[54:60] ), float ( line[60:66] )

                #if occ < 1.0 :
                rid = "%d.%s" % (resi,cid)

                if rid in rids :
                    r = rids[rid]

                    #if line[0:6] == "HETATM" :
                    #    print rid, r.id.position, r.type

                    if aname in r.atomsMap :
                        ats = r.atomsMap[aname]
                        found = False
                        for at in ats :
                            if at.altLoc == aloc :
                                at.Q = bfac
                                at.bfactor = 150.0 * (1.0 - at.Q)
                                #at.bfactor = 0

                                #at.occupancy = 1.0 # max(0,at.Q)

                                dval = self.cur_dmap.interpolated_values ( [ at.coord()  ], self.cur_mol.openState.xform ).astype(numpy.float64, copy=False)[0]
                                V = (dval - minD) / (maxD - minD)
                                #V = 0.5 * max ( min(V,1.0), 0.0 ) + 0.5
                                V = max ( min(V,1.0), 0.0 )
                                #at.occupancy = V
                                #at.occupancy = 1.0
                                #at.bfactor = at.bfactor / V

                                if halfMap1 : at.Q1 = at.Q
                                if halfMap2 : at.Q2 = at.Q
                                found = True
                        if not found :
                            #print " -xx- %d.%s - atom %s - loc %s" % (resi, cid, aname, aloc)
                            continue
                    else :
                        #print " -xx- %d.%s - atom %s" % (resi,cid, aname)
                        continue

                else :
                    #print " -xx- %d.%s " % (resi,cid)
                    continue


        fin.close ()

        qscores.QStats1 (self.cur_mol, chainId)
        qscores.SaveQStats ( self.cur_mol, self.chain.get(), self.cur_dmap, gSigma, float(self.mapRes.get()) )

        if 0 :
            self.SaveQsBfs ( self.cur_mol, 50.0 )
            self.SaveQsBfs ( self.cur_mol, 100.0 )
            self.SaveQsBfs ( self.cur_mol, 150.0 )
            self.SaveQsBfs ( self.cur_mol, 200.0 )
            self.SaveQsBfs ( self.cur_mol, 300.0 )

        scBB, scSC = [], []

        doRess = []
        for r in self.cur_mol.residues :
            if r.id.chainId == chainId :
                doRess.append ( r )

        print "Q for %d res..." % ( len(doRess) )
        for r in doRess :

            qscores.CalcResQ (r, None, None, useOld=True )

            if r.isProt or r.isNA :
                r.score1 = r.scQ
                r.score2 = r.scQ
                if r.bbQ != None : scBB.append ( r.bbQ )
                if r.scQ != None : scSC.append ( r.scQ )
            else :
                r.score1 = r.Q
                r.score2 = r.Q



        if len (scSC) > 0 :
            scMin, scMax, scAvg = min(scSC), max(scSC), numpy.average(scSC)
            bbMin, bbMax, bbAvg = min(scBB), max(scBB), numpy.average(scBB)

            print " - average Q sc : %.2f - %.2f, avg %.2f" % (scMin, scMax, scAvg )
            print " - average Q bb : %.2f - %.2f, avg %.2f" % (bbMin, bbMax, bbAvg )
            print ""

            if 1 :
                self.minScore1, self.maxScore1 = 0.0,1
                self.minScore2, self.maxScore2 = 0.0,1
            else :
                self.minScore1, self.maxScore1 = 0.0,max(scSC)
                self.minScore2, self.maxScore2 = 0.0,max(scBB)

            self.UpdateSeq ()


        self.ShowQScores ()


        #self.QStats ()
        #self.QStatsRNA()



    def SaveQsBfs ( self, mol, f ) :

        for at in mol.atoms :
            at.bfactor = f * (1.0 - at.Q)

        molPath = os.path.splitext(mol.openedAs[0])[0]

        nname = molPath + "__Bf%.0f__.pdb" % f
        print " - saving %s" % nname
        chimera.PDBio().writePDBfile ( [mol], nname )



    def QStats ( self ) :

        mol, dmap, chainId = self.cur_mol, self.cur_dmap, self.chain.get()

        ress = []
        for r in mol.residues :
            if r.id.chainId == chainId :
                ress.append ( r )


        sByType = {}
        rByType = {}

        def addType (tp, r, score) :
            if not tp in sByType :
                rByType[tp] = []
                sByType[tp] = []
            rByType[tp].append ( [score, r] )
            sByType[tp].append ( [score] )


        for r in ress :
            if r.isProt and r.type == "LEU" :
                avg = (r.atomsMap["CD1"][0].Q + r.atomsMap["CD2"][0].Q)/2.0
                addType ( "LEU(CD)", r, avg )
            if r.isProt and r.type == "LEU" and r.id.position==114 :
                avg = (r.atomsMap["CD1"][0].Q + r.atomsMap["CD2"][0].Q)/2.0
                addType ( "LEU_114(CD)", r, avg )

            if r.isProt and r.type == "VAL" :
                avg = (r.atomsMap["CG1"][0].Q + r.atomsMap["CG2"][0].Q)/2.0
                addType ( "VAL(CG)", r, avg )
            if r.isProt and r.type == "VAL" and r.id.position==33 :
                avg = (r.atomsMap["CG1"][0].Q + r.atomsMap["CG2"][0].Q)/2.0
                addType ( "VAL_33(CG)", r, avg )

            if r.isProt and r.type == "ARG" :
                avg = (r.atomsMap["NH1"][0].Q + r.atomsMap["NH2"][0].Q)/2.0
                addType ( "ARG(NH)", r, avg )
            if r.isProt and r.type == "ARG" and r.id.position==76 :
                avg = (r.atomsMap["NH1"][0].Q + r.atomsMap["NH2"][0].Q)/2.0
                addType ( "ARG_76(NH)", r, avg )
            if r.isProt and r.type == "ARG" and r.id.position==9 :
                avg = (r.atomsMap["NH1"][0].Q + r.atomsMap["NH2"][0].Q)/2.0
                addType ( "ARG_9(NH)", r, avg )

            if r.isProt and r.type == "LYS" :
                avg = r.atomsMap["NZ"][0].Q
                addType ( "LYS(NZ)", r, avg )

            if r.isProt and r.type == "ASP" :
                avg = (r.atomsMap["OD1"][0].Q + r.atomsMap["OD2"][0].Q)/2.0
                addType ( "ASP(OD)", r, avg )
            if r.isProt and r.type == "ASP" and r.id.position==42 :
                avg = (r.atomsMap["OD1"][0].Q + r.atomsMap["OD2"][0].Q)/2.0
                addType ( "ASP_42(OD)", r, avg )
            if r.isProt and r.type == "ASP" and r.id.position==131 :
                avg = (r.atomsMap["OD1"][0].Q + r.atomsMap["OD2"][0].Q)/2.0
                addType ( "ASP_131(OD)", r, avg )
            if r.isProt and r.type == "ASP" and r.id.position==171 :
                avg = (r.atomsMap["OD1"][0].Q + r.atomsMap["OD2"][0].Q)/2.0
                addType ( "ASP_171(OD)", r, avg )

            if r.isProt and r.type == "GLU" :
                avg = (r.atomsMap["OE1"][0].Q + r.atomsMap["OE2"][0].Q)/2.0
                addType ( "GLU(OE)", r, avg )
            if r.isProt and r.type == "GLU" and r.id.position==17 :
                avg = (r.atomsMap["OE1"][0].Q + r.atomsMap["OE2"][0].Q)/2.0
                addType ( "GLU_17(OE)", r, avg )
            if r.isProt and r.type == "GLU" and r.id.position==27 :
                avg = (r.atomsMap["OE1"][0].Q + r.atomsMap["OE2"][0].Q)/2.0
                addType ( "GLU_27(OE)", r, avg )
            if r.isProt and r.type == "GLU" and r.id.position==67 :
                avg = (r.atomsMap["OE1"][0].Q + r.atomsMap["OE2"][0].Q)/2.0
                addType ( "GLU_67(OE)", r, avg )
            if r.isProt and r.type == "GLU" and r.id.position==134 :
                avg = (r.atomsMap["OE1"][0].Q + r.atomsMap["OE2"][0].Q)/2.0
                addType ( "GLU_134(OE)", r, avg )


            if r.isProt or r.isNA :
                if r.scQ :
                    addType ( r.type, r, r.scQ )
            else :
                addType ( r.type, r, r.Q )

        avgs = []
        for rtype, ra in sByType.iteritems () :
            avgs.append ( [numpy.average (ra), rtype, numpy.std (ra)] )


        from chimera.resCode import protein3to1
        from chimera.resCode import nucleic3to1

        # sort by avg score
        #avgs.sort ( reverse=True, key=lambda x: x[0] )

        # sort by residue type
        avgs.sort ( reverse=False, key=lambda x: x[1] )


        mapName = os.path.splitext(dmap.name)[0]
        molName = os.path.splitext(mol.name)[0]
        mdir, mpfile = os.path.split(dmap.data.path)
        foname = mdir + "/" + mapName + "__" + molName + ".txt"


        print " - scores to: " + foname
        fp = open (foname,"w")

        for avgScore, rtype, sdev in avgs :

            rscores = rByType[rtype]
            if len(rscores) > 0 :
                rscores.sort ( reverse=True, key=lambda x: x[0] )
                hr = rscores[0]
                R = hr[1]
                highestScore = hr[0]
                numRes = len(rscores)

                rts = ""
                if R.isProt : rts = protein3to1[R.type]
                elif R.isNA : rts = nucleic3to1[R.type]

                print "%s\t%s\t%d\t%f\t%f\t%d\t.%s\t%f" % (rtype, rts, numRes, avgScore, sdev, R.id.position, R.id.chainId, highestScore)
                fp.write ( "%s\t%s\t%d\t%f\t%f\t%d\t.%s\t%f\n" % (rtype, rts, numRes, avgScore, sdev, R.id.position, R.id.chainId, highestScore) )

        fp.close()




    def QStatsRNA ( self ) :

        mol, dmap, chainId = self.cur_mol, self.cur_dmap, self.chain.get()

        #SetBBAts ( mol )

        print ""
        print "RNA stats for chain %s" % chainId
        print ""

        ress = []
        for r in mol.residues :
            if r.id.chainId == chainId :
                ress.append ( r )

        sByType = {}
        rByType = {}

        def addType (tp, r, score) :
            if not tp in sByType :
                rByType[tp] = []
                sByType[tp] = []
            rByType[tp].append ( [score, r] )
            sByType[tp].append ( [score] )


        scAts = []
        bbAts = []
        allAts = []

        for r in ress :
            if r.isNA :

                avg = numpy.average ( [at.Q for at in r.scAtoms] )
                #addType ( nucleic3to1[r.type] + "_SC", r, avg )
                addType ( r.type + "_SC", r, avg )

                avg = numpy.average ( [at.Q for at in r.bbAtoms] )
                #addType ( nucleic3to1[r.type] + "_BB", r, avg )
                addType ( r.type + "_BB", r, avg )

                scAts.extend ( r.scAtoms )
                bbAts.extend ( r.bbAtoms )
                allAts.extend ( [at for at in r.atoms if at.element.name != "H"] )


        avgQ = numpy.average ( [at.Q for at in allAts] )
        avgQbb = numpy.average ( [at.Q for at in bbAts] )
        avgQsc = numpy.average ( [at.Q for at in scAts] )

        sQ = numpy.std ( [at.Q for at in allAts] )
        sQbb = numpy.std ( [at.Q for at in bbAts] )
        sQsc = numpy.std ( [at.Q for at in scAts] )

        avgs = []
        for rtype, ra in sByType.iteritems () :
            avgs.append ( [numpy.average (ra), rtype, numpy.std (ra)] )


        from chimera.resCode import protein3to1
        from chimera.resCode import nucleic3to1

        # sort by avg score
        #avgs.sort ( reverse=True, key=lambda x: x[0] )

        # sort by residue type
        avgs.sort ( reverse=False, key=lambda x: x[1] )


        mapName = os.path.splitext(dmap.name)[0]
        molName = os.path.splitext(mol.name)[0]
        mdir, mpfile = os.path.split(dmap.data.path)
        foname = mdir + "/" + mapName + "__" + molName + "_rscores.txt"


        print " - scores to: " + foname
        fp = open (foname,"w")

        print ""
        print "Map\tModel\tQ_All\tQ_Backbone\tQ_SideChain\tStdQ_All\tStdQ_Backbone\tStdQ_SideChain"
        print "%s\t%s\t%f\t%f\t%f\t%f\t%f\t%f" % (mapName, molName, avgQ, avgQbb, avgQsc, sQ, sQbb, sQsc)
        print ""

        fp.write ( "%n" )
        fp.write ( "Map\tModel\tQ_All\tQ_Backbone\tQ_SideChain\tStdQ_All\tStdQ_Backbone\tStdQ_SideChain" )
        fp.write ( "%s\t%s\t%f\t%f\t%f\t%f\t%f\t%f" % (mapName, molName, avgQ, avgQbb, avgQsc, sQ, sQbb, sQsc) )
        fp.write ( "%n" )


        for avgScore, rtype, sdev in avgs :

            rscores = rByType[rtype]
            if len(rscores) > 0 :
                rscores.sort ( reverse=True, key=lambda x: x[0] )
                hr = rscores[0]
                R = hr[1]
                highestScore = hr[0]
                numRes = len(rscores)

                rts = ""
                if R.isProt : rts = protein3to1[R.type]
                elif R.isNA : rts = nucleic3to1[R.type]

                print "%s\t%s\t%d\t%f\t%f\t%d\t.%s\t%f" % (rtype, rts, numRes, avgScore, sdev, R.id.position, R.id.chainId, highestScore)
                fp.write ( "%s\t%s\t%d\t%f\t%f\t%d\t.%s\t%f\n" % (rtype, rts, numRes, avgScore, sdev, R.id.position, R.id.chainId, highestScore) )

        fp.close()




    def SA_Q (self) :

        ress = []
        try :
            ress = self.seqRes
        except :
            pass

        if len ( ress ) == 0 :
            umsg ( "No molecule/chain selected?" )
            return


        ok = True
        try :
            print self.cur_dmap.name
        except :
            status ( "Selected map not found; please choose another map" )
            self.dmap.set ("")
            ok = False

        try :
            print self.cur_mol.name
        except :
            status ( "Selected model not found; please choose another model" )
            self.struc.set ("")
            self.chain.set ("")
            self.RemoveSeq ()
            ok = False

        if not ok :
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

                if hasattr (r, 'scQ') and r.scQ != None :
                    fp.write ( "%s\t%d\t%f\t%f\n" % (r.type, r.id.position, r.scQ, r.SAArea) )
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


        ok = True
        try :
            print self.cur_dmap.name
        except :
            status ( "Selected map not found; please choose another map" )
            self.dmap.set ("")
            ok = False

        try :
            print self.cur_mol.name
        except :
            status ( "Selected model not found; please choose another model" )
            self.struc.set ("")
            self.chain.set ("")
            self.RemoveSeq ()
            ok = False

        if not ok :
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


        ok = True
        try :
            print self.cur_dmap.name
        except :
            status ( "Selected map not found; please choose another map" )
            self.dmap.set ("")
            ok = False

        try :
            print self.cur_mol.name
        except :
            status ( "Selected model not found; please choose another model" )
            self.struc.set ("")
            self.chain.set ("")
            self.RemoveSeq ()
            ok = False

        if not ok :
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

        RES = float(self.mapRes.get())

        avgQrna = -0.1574 * RES + 1.0673 # rna
        avgQprot = -0.1794 * RES + 1.1244 # protein
        avgQIon =  -0.1103 * RES + 1.0795 # ion
        avgQWater =  -0.0895 * RES + 1.0001 # water

        print " - res %.2f - exp Q-score: %.2f" % (RES, avgQprot)

        self.minScore1, self.maxScore1 = 0.0,avgQprot
        self.minScore2, self.maxScore2 = 0.0,avgQprot



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
                self.PtsToMap ( points, dmap, R, dmap.name + label, False, alpha=0.2 if self.showMesh.get() else 0.4 )
                if self.showMesh.get () :
                    self.PtsToMap ( points, dmap, R, dmap.name + label + "_mesh", True )




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
            SetBBAts(self.cur_mol)
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


    def PtsToMap ( self, points, dmap, atomRad, nname, showMesh = False, alpha=0.2 ) :

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
                            status ( "model not found" )
                            self.chain.set("")
                            self.struc.set("")
                            self.RemoveSeq ()
                            return
                    else :
                        try :
                            status ( "Sequence: ?/? %d/%d" % ( ri, resEnd.id.position ) )
                        except :
                            status ( "model not found" )
                            self.chain.set("")
                            self.struc.set("")
                            self.RemoveSeq ()
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
            SetBBAts(r.molecule)
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

        start = time.time()
        qq = qscores.Qscore  ( [selAtom], dmap, gSigma, allAtTree=allAtTree, show=0, log=1, numPts=8, toRAD=2.0, dRAD=0.1, minD=minD, maxD=maxD, fitg=0 )
        end = time.time()
        print " - time: %f" % ( (end - start) )

        start = time.time()
        qq = qscores.Qscore ( [selAtom], dmap, gSigma, allAtTree=allAtTree, show=0, log=1, numPts=50, toRAD=2.0, dRAD=0.1, minD=minD, maxD=maxD, fitg=0 )
        end = time.time()
        print " - time: %f" % ( (end - start) )

        #CC, CCm, yds, err = rr

        print r



    def Q_sel (self) :

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
            SetBBAts(r.molecule)
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

        #allAtTree = None
        #print "-"

        import time
        start = time.time()


        if 0 :
            print "_sigma____________________________"

            sigma = RadAts ( [selAtom], dmap, allAtTree=allAtTree, show=1, log=1, numPts=30, toRAD=2, dRAD=0.5 )
            res = sigma * numpy.pi * numpy.sqrt(2.0)

            end = time.time()
            print "%s - sigma: %.3f, res: %.3f, time: %f" % ( selAtom.name, sigma, res, (end - start) )


        elif 1 :
            print ""
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

            start = time.time()
            qs = qscores.Qscore ( [selAtom], dmap, gSigma, allAtTree=allAtTree, show=0, log=0, numPts=8, toRAD=2.0, dRAD=0.1, minD=minD, maxD=maxD, fitg=0 )
            #qs, yds, err = qs
            end = time.time()
            print " - sigma: %.3f, Q-score: %.3f, time: %f" % ( gSigma, qs, (end - start) )


            if 1 :
                #def QscorePt ( atPt, xfI, dmap, sigma, allAtTree = None, log=0, numPts=8, toRAD=2.0, dRAD=0.5, minD=None, maxD=None, fitg=0 ) :

                pt = selAtom.coord().data()
                xfI = selAtom.molecule.openState.xform

                #pt = selAtom.xformCoord().data()
                #xfI = chimera.Xform()

                start = time.time()
                qs = qscores.QscorePt ( pt, xfI, dmap, gSigma, allAtTree=allAtTree, log=0, numPts=8, toRAD=2.0, dRAD=0.1, minD=minD, maxD=maxD, fitg=0 )
                #qs, yds, err = qs
                end = time.time()
                print " - sigma: %.3f, Q-score Pt: %.3f, time: %f" % ( gSigma, qs, (end - start) )


            if 1 :
                print "Atoms in %d.%s %s" % (selAtom.residue.id.position, selAtom.residue.id.chainId, selAtom.residue.type)
                #print "-"

                avg, N = 0.0, 0.0
                #bbAts, scAts, baseAts, sugarAts = [], [], [], []

                for at in selAtom.residue.atoms :

                    at.Q = qscores.Qscore ( [at], dmap, gSigma, allAtTree=allAtTree, show=0, log=0, numPts=8, toRAD=2.0, dRAD=0.1, minD=minD, maxD=maxD )
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
            SetBBAts(r.molecule)
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
                rr = qscores.Qscore ( [selAtom], dmap, gSigma, allAtTree=allAtTree, show=1, log=1, numPts=20, toRAD=2.0, dRAD=0.5, minD=minD, maxD=maxD, fitg=1 )
                qs, yds, err = rr

            elif 1 :
                rr = qscores.Qscore ( [selAtom], dmap, gSigma, allAtTree=allAtTree, show=0, log=1, numPts=30, toRAD=2.0, dRAD=0.1, minD=minD, maxD=maxD, fitg=0 )
                qs, yds, err = rr

            else :
                qs = qscores.Qscore ( [selAtom], dmap, gSigma, allAtTree=allAtTree, show=0, log=0, numPts=8, toRAD=2.0, dRAD=0.1, minD=minD, maxD=maxD, fitg=0 )

            end = time.time()
            print " - sigma: %.3f, Q-score: %.3f, time: %f" % ( gSigma, qs, (end - start) )

            print "Atoms in %d.%s %s" % (selAtom.residue.id.position, selAtom.residue.id.chainId, selAtom.residue.type)
            #print "-"


            if 0 :
                avg, N = 0.0, 0.0
                #bbAts, scAts, baseAts, sugarAts = [], [], [], []

                for at in selAtom.residue.atoms :
                    at.Q = qscores.Qscore ( [at], dmap, gSigma, allAtTree=allAtTree, show=0, log=0, numPts=8, toRAD=2.0, dRAD=0.1, minD=minD, maxD=maxD )
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



    def CalcSelQ (self) :

        # show sigma for a side chain

        atoms = chimera.selection.currentAtoms()
        if len ( atoms ) == 0 :
            umsg ( "No selected atoms found" )
            return

        dmap = self.cur_dmap
        mol = atoms[0].molecule

        umsg ( "Calculating Q-scores of %d atoms..." % len(atoms) )


        #selAtom = selAts[0]
        #r = selAtom.residue
        #print "Res: %s - %d.%s - %s - Atom: %s" % (r.type, r.id.position, r.id.chainId, r.molecule.name, selAtom.name)

        #sigma = 0.4
        print " - in map: %s" % dmap.name
        print " - mol: %s" % mol.name
        print " - sigma: %.2f" % gSigma

        if 1 or not hasattr ( mol.name, 'bbats' ) :
            SetBBAts(mol)
            mol.bbats = True

        ats = [at for at in mol.atoms if not at.element.name == "H"]
        if self.showH.get() :
            ats = mol.atoms

        points = _multiscale.get_atom_coordinates ( ats, transformed = False )
        print " - search tree: %d/%d ats" % ( len(ats), len(mol.atoms) )
        allAtTree = AdaptiveTree ( points.tolist(), ats, 1.0)

        minD, maxD = qscores.MinMaxD ( dmap )

        import time
        start = time.time()


        from chimera import tasks, CancelOperation
        task = tasks.Task('Calculating Q-scores', modal = True)

        avg = 0.0

        import traceback


        try :

            for ai, at in enumerate ( atoms ) :

                at.Q = qscores.Qscore ( [at], dmap, gSigma, allAtTree=allAtTree, show=0, log=0, numPts=8, toRAD=2.0, dRAD=0.1, minD=minD, maxD=maxD, fitg=0 )
                at.bfactor = at.Q
                avg += at.Q

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


        cc = ResCC ( mol, atoms, 4.0, dmap )
        print " - CC: %.3f" % cc

        for at in atoms :
            print " - atom: %s %d.%s %s : %.2f" % (at.residue.type, at.residue.id.position, at.residue.id.chainId, at.name, at.Q)

        avgq = avg / float(len(atoms))
        if len(atoms) > 1 :
            umsg ( "Q-score of %d atoms: %.2f" % (len(atoms), avgq) )
        else :
            umsg ( "Q-score of %d atom: %.2f" % (len(atoms), avgq) )




    def AProfs (self) :

        mol = self.cur_mol
        if self.cur_mol == None :
            umsg ("Select a molecule first")
            return []

        chainId = self.chain.get()

        dmap = self.cur_dmap
        print " - in map: %s" % dmap.name

        if 1 or not hasattr ( mol, 'bbats' ) :
            SetBBAts(mol)
            mol.bbats = True

        ats = [at for at in self.cur_mol.atoms if not at.element.name == "H"]
        points = _multiscale.get_atom_coordinates ( ats, transformed = False )
        print " - search tree: %d/%d ats" % ( len(ats), len(mol.atoms) )
        allAtTree = AdaptiveTree ( points.tolist(), ats, 1.0)

        #sigma = 0.4
        minD, maxD = qscores.MinMaxD ( dmap )
        print " - mind: %.3f, maxd: %.3f" % (minD, maxD)



        def doAt (at, arr) :
            rr = qscores.Qscore ( [at], dmap, gSigma, allAtTree=allAtTree, show=0, log=0, numPts=10, toRAD=3.0, dRAD=0.1, minD=minD, maxD=maxD, fitg=1 )
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
            SetBBAts(mol)
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


    def AddAtom ( self ) :

        if 0 :
            print self.addText.get()

            ats = chimera.selection.currentAtoms()

            if len(ats) == 1 :
                at = ats[0]
                print "Placing near: ", at.name

            mol = self.cur_mol
            if self.cur_mol == None :
                umsg ("Select a molecule first")
                return []

            chainId = "D"

            nres = mol.newResidue ("CA", chimera.MolResId(chainId, 5))
            nat = mol.newAtom ('CA', chimera.Element(20))
            nres.addAtom( nat )
            nat.setCoord ( at.coord() )
            nat.drawMode = nat.Sphere
            nat.color = chimera.MaterialColor( 0.0, 1.0, 0.0, 1.0 )
            nat.display = True

        else :

            iAt = 1
            for at in chimera.selection.currentAtoms() :

                chainId = "B"

                nres = self.cur_mol.newResidue ("MG", chimera.MolResId(chainId, iAt))
                nat = self.cur_mol.newAtom ('MG', chimera.Element(12))
                nres.addAtom( nat )
                nat.setCoord ( at.coord() )
                nat.drawMode = nat.Sphere
                nat.color = chimera.MaterialColor( 1.0, 0.0, 0.0, 1.0 )
                nat.display = True

                print " - added atom %s, %d.%s -> %d.%s" % (at.name, at.residue.id.position, at.residue.id.chainId, iAt, chainId)
                iAt += 1


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



    def DMS ( self ) :

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
                r.ribbonColor = chimera.MaterialColor ( .7, .7, .7, 1.0 )
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


    def SS ( self ) :

        print "dms"

        mol = self.cur_mol
        rmap = {}
        for r in mol.residues :
            rmap[r.id.position] = r
            r.dms = None

            r.ribbonColor = chimera.MaterialColor ( .7, .7, .7, 1.0 )

        dms = []
        fp = open ( "/Users/greg/_data/Ribozyme/sec.txt" )
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
                        r = rmap[i]
                        r.ribbonColor = chimera.MaterialColor ( C[0], C[1], C[2], 1.0 )

            print ""


    def AddRes ( self ) :

        #startI = self.seqRes [ max(self.seqSel[0],0) ].id.position
        #endI = self.seqRes [ min(self.seqSel[1],len(self.seqRes)-1) ].id.position

        print ""
        print "AddRes"

        mol = self.cur_mol
        if self.cur_mol == None :
            umsg ("Select a molecule first")
            return []

        chainId = self.chain.get()

        resToAdd = self.addRess.get().upper().strip().replace(" ", "")
        print " - res to add:", resToAdd

        selAt, selReg = None, None
        import _surface
        import _molecule
        for c in chimera.selection.currentContents()[0] :
            if type(c) == _surface.SurfacePiece :
                print " - sp",
                selSp = c
                if hasattr ( selSp, 'region' ) :
                    selReg = selSp.region
                    print " - reg: %d" % selReg.rid
                else :
                    print "?"
            elif type(c) == _molecule.Atom :
                selAt = c
                print " - atom: %s" % selAt.name

        if resToAdd.lower() == "nag" :
            print " - adding nag"
            self.AddNAG ( selAt, selReg )



    def AddNAG ( self, selAt, selReg ) :

        nmol = chimera.PDBio().readPDBfile ( "/Users/greg/_data/NL63/NAG.pdb" )[0]
        print " - read %s - %d atoms" % ( nmol.name, len(nmol.atoms) )

        lastRid = 0
        cid = selAt.residue.id.chainId
        for r in selAt.molecule.residues :
            if r.id.chainId == cid :
                if r.id.position > lastRid :
                    lastRid = r.id.position

        if selAt.residue.type == "ASN" :
            atN = selAt.residue.atomsMap["ND2"][0]
            atO = selAt.residue.atomsMap["OD1"][0]
            atCG = selAt.residue.atomsMap["CG"][0]
            vN = atN.coord() - atCG.coord(); vN.normalize()
            vO = atO.coord() - atCG.coord(); vO.normalize()
            vA = chimera.cross ( vO, vN ); vA.normalize()

            vC1 = chimera.Xform.rotation (vA, 124.669) .apply (vN*-1.0)
            vC1.normalize()
            pC1 = atN.coord() + vC1 * 1.450

            pC1_ = nmol.residues[0].atomsMap["C1"][0].coord()
            pO1_ = nmol.residues[0].atomsMap["O1"][0].coord()
            vC1_ = pC1_ - pO1_; vC1_.normalize()

            rax = chimera.cross ( vC1_, vC1 ); rax.normalize()
            ang = numpy.arccos ( vC1_ * vC1 ) * 180.0 / numpy.pi

            xf = chimera.Xform.translation ( pC1_.toVector() * -1.0 )
            xf.premultiply ( chimera.Xform.rotation(rax, ang) )
            xf.premultiply ( chimera.Xform.translation ( pC1.toVector() ) )

            aMap = {}
            nres = selAt.molecule.newResidue ( nmol.residues[0].type, chimera.MolResId(cid, lastRid+1))

            for at in nmol.atoms :
                if at.element.name == "H" :
                    continue
                elif 1 and at.name == "O1" :
                    continue

                nat = selAt.molecule.newAtom (at.name, chimera.Element(at.element.number))
                aMap[at] = nat
                nres.addAtom( nat )
                nat.drawMode = nat.EndCap
                nat.setCoord ( xf.apply ( at.coord()) )
                nat.display = True
                if nat.element.name.upper() in atomColors : nat.color = atomColors[nat.element.name.upper()]

            for bond in nmol.bonds :
                if bond.atoms[0] in aMap and bond.atoms[1] in aMap :
                    nb = selAt.molecule.newBond ( aMap[bond.atoms[0]], aMap[bond.atoms[1]] )
                    nb.display = nb.Smart
                    nb.drawMode = nb.Stick

            nb = selAt.molecule.newBond ( atN, nres.atomsMap["C1"][0] )
            nb.display = nb.Smart
            nb.drawMode = nb.Stick

            if selReg == None :
                return

            segMap = selReg.segmentation.seg_map

            print " - map:", segMap.name
            zoneR = segMap.data.step[0]/2.0
            rpoints = numpy.concatenate ( [selReg.map_points() for r in [selReg]], axis=0 ).astype ( numpy.float32 )
            rdata = VolumeData.zone_masked_grid_data ( segMap.data, rpoints, zoneR )
            rmat = rdata.matrix()

            ##gdata = VolumeData.Array_Grid_Data ( ndata.full_matrix(), segMap.data.origin, segMap.data.step, segMap.data.cell_angles, name = "atom masked" )
            #nv = VolumeViewer.volume.volume_from_grid_data ( rdata )
            #nv.name = "helix mask 2"

            maxAng, maxD, angD = 0, -1e9, 1.0
            for ang in range ( 0, int(round(360.0/angD)), 1 ) :

                xf = chimera.Xform.translation ( pC1.toVector() * -1.0 )
                xf.premultiply ( chimera.Xform.rotation(vC1, angD) )
                xf.premultiply ( chimera.Xform.translation ( pC1.toVector() ) )

                for at in nres.atoms :
                    at.setCoord ( xf.apply(at.coord()) )

                points = _multiscale.get_atom_coordinates ( nres.atoms, transformed = True )
                _contour.affine_transform_vertices ( points, Matrix.xform_matrix(segMap.openState.xform.inverse()) )
                #_contour.affine_transform_vertices ( points, M )

                values, outside = VolumeData.interpolate_volume_data ( points, rdata.xyz_to_ijk_transform, rmat )

                #values = nv.interpolated_values ( points, selAt.molecule.openState.xform )
                #olap, corr, other = overlap_and_correlation ( rpoint_weights, rmap_values )
                avgD = numpy.average ( values )
                #print "%.1f\t%.4f" % (ang, avgD)

                if avgD > maxD :
                    maxD = avgD
                    maxAng = round(float(ang)/angD)

            print "Max ang: %.3f" % maxAng
            xf = chimera.Xform.translation ( pC1.toVector() * -1.0 )
            xf.premultiply ( chimera.Xform.rotation(vC1, maxAng) )
            xf.premultiply ( chimera.Xform.translation ( pC1.toVector() ) )

            for at in nres.atoms :
                at.setCoord ( xf.apply(at.coord()) )





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





    def AddRes_ ( self ) :

        seq = self.addRess.get()
        print "Seq:", seq
        if len(seq) > 0 :
            molref.AddRes ( seq[0] )


    def AddResN ( self ) :

        print "Adding"
        seq = self.addRess.get()
        print "Seq:", seq

        ress = chimera.selection.currentResidues()
        if len(ress) == 0 :
            umsg ( "Select a residue" )
            return
        elif len(ress) > 1 :
            umsg ( "Select only a residue to add before" )
            return
        print ress

        mol = ress[0].molecule

        if len(seq) > 0 :
            molref.AddResN ( seq[0], ress[0] )


        SetBBAts(mol)



    def AddResC ( self ) :

        print "Adding"


    def Refine ( self ) :

        ress = chimera.selection.currentResidues()
        if len(ress) == 0 :
            umsg ( "Select some residues..." )
            return

        dmap = self.cur_dmap

        molref.Refine ( ress, dmap )




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




    def RMSD ( self ) :

        mols = []
        for m in chimera.openModels.list() :
            if m.display == True and type(m) == chimera.Molecule :
                mols.append ( m )

        if len(mols) != 2 :
            umsg ( "Make at least two molecules visible" )
            return

        m1, m2 = mols

        SetBBAts ( m1 )
        SetBBAts ( m2 )

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
            SetBBAts(r.molecule)
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
            SetBBAts(r.molecule)
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
            SetBBAts(r.molecule)
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
            SetBBAts(r.molecule)
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

    SetBBAts ( mol )


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





def BBsegs ( ress ) :

    bbs = []

    firstRi, atRi = 0, 1
    for r in ress[1:] :
        if ress[atRi].id.position > ress[atRi-1].id.position + 1 or r.rtype == "?" :
            bbs.append ( ress[firstRi:atRi] )
            firstRi = atRi
        atRi += 1

    bbs.append ( ress[firstRi:atRi] )

    return bbs




def SSEs ( allRess ) :

    if len(allRess) < 1 :
        return []

    sses, ss = [], ""

    res, rStart = allRess[0], allRess[0]
    #print "  - at first res / pos: %d " % res.id.position
    if res.isHelix :
        ss = "H"
    elif res.isSheet or res.isStrand :
        ss = "E"
    else :
        ss = "_"

    ress = [ res ]
    lastRes = rStart
    for res in allRess [1:] :

        if res.id.position > lastRes.id.position + 1 :
            print " - gap at", res.id.position
            sses.append ( [rStart.id.position, lastRes.id.position, ss, ress] )
            ress = []
            rStart = res
            if res.isHelix :
                ss = "H"
            elif res.isSheet or res.isStrand :
                ss = "E"
            else :
                ss = "_"

        if res.isHelix :
            if ss != "H" :
                #print "%s -> H - at %d rid %d | %d->%d, %d res" % (ss, i, res.id.position, rStart.id.position, lastRes.id.position, len(ress))
                sses.append ( [rStart.id.position, lastRes.id.position, ss, ress] )
                ress = []
                rStart = res
                ss = "H"
        elif res.isSheet or res.isStrand :
            if ss != "E" :
                #print "%s -> E - at %d rid %d | %d->%d, %d res" % (ss, i, res.id.position, rStart.id.position, lastRes.id.position, len(ress))
                sses.append ( [rStart.id.position, lastRes.id.position, ss, ress] )
                ress = []
                rStart = res
                ss = "E"
        else :
            if ss == "H" or ss == "E" :
                #print "%s -> _ at %d rid %d | %d->%d, %d res" % (ss, i, res.id.position, rStart.id.position, lastRes.id.position, len(ress))
                sses.append ( [rStart.id.position, lastRes.id.position, ss, ress] )
                ress = []
                rStart = res
                ss = "_"

        ress.append ( res )
        lastRes = res

    #print "Done at rid %d - %s | %d->%d, %d res" % ( res.id.position, ss, rStart.id.position, res.id.position, len(ress))
    sses.append ( [rStart.id.position, res.id.position, ss, ress] )
    return sses




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
            res.scQ = res.scZ

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
        SetBBAts(r.molecule)
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
                SetBBAts ( nmol )
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
            SetBBAts ( nmol )
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
            SetBBAts ( nmol )
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
    display_threshold = 0.95

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

    fpoints, fpoint_weights = fit_points_g ( molg )
    map_values = dmap.interpolated_values ( fpoints, mol.openState.xform )
    olap, corr1, corr2 = FitMap.overlap_and_correlation ( fpoint_weights, map_values )
    return corr1



def fit_points_g (fdata, threshold = 1e-5) :

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



def SetBBAts ( mol ) :

    #if hasattr ( mol, "bbats" ) :
    #    return
    #mol.bbats = True

    print " - setting bbAts in %s" % mol.name
    for r in mol.residues :

        #r.isProt = "C" in r.atomsMap and "CA" in r.atomsMap and "N" in r.atomsMap
        #r.isProt = "CA" in r.atomsMap
        #r.isNA = "O3'" in r.atomsMap and "O5'" in r.atomsMap

        from chimera.resCode import nucleic3to1
        from chimera.resCode import protein3to1
        protein3to1['HSD'] = protein3to1['HIS']
        protein3to1['HSE'] = protein3to1['HIS']

        r.isProt = r.type in protein3to1
        r.isNA = r.type in nucleic3to1

        #r.score1 = None
        #r.score2 = None

        if r.isProt :
            r.rtype = "prot"
        elif r.isNA :
            r.rtype = "na"
        else :
            r.rtype = "?"


        if r.isNA :
            try :
                if nucleic3to1[r.type] == "G" :
                    r.baseAt = r.atomsMap["N9"][0]
                elif nucleic3to1[r.type] == "C" :
                    r.baseAt = r.atomsMap["N1"][0]
                elif nucleic3to1[r.type] == "A" :
                    r.baseAt = r.atomsMap["N9"][0]
                elif nucleic3to1[r.type] == "U" :
                    r.baseAt = r.atomsMap["N1"][0]
            except :
                #print " - baseAt not found - "
                pass


        r.bbAtoms = []
        r.scAtoms = []

        if r.isProt :
            for a in r.atoms :
                if a.element.name == "H" :
                    a.isBB, a.isSC = False, False
                    continue
                n = a.name
                a.isBB = n=="C" or n=="CA" or n=="O" or n=="N" or n=="OT1" or n=="OT2"
                a.isSC = not a.isBB
                if a.isBB :
                    r.bbAtoms.append ( a )
                else :
                    r.scAtoms.append ( a )

                a.isSugar, a.isBase = False, False

        elif r.isNA :
            for a in r.atoms :
                if a.element.name == "H" :
                    a.isBB, a.isSC = False, False
                    continue
                n = a.name

                a.isBB = n=="P" or n=="O1P" or n=="O2P" or n=="OP1" or n=="OP2" or n=="O5'" or n=="C5'" or n=="O3'"
                a.isSugar = n=="C1'" or n=="C2'" or n=="O4'" or n=="O2'" or n=="C3'" or n=="C4'"
                a.isBB = a.isBB or a.isSugar

                a.isBase = False

                if nucleic3to1[r.type] == "G" :
                    a.isBase = n=="N9" or n=="C8" or n=="N7" or n=="C5" or n=="C4" or n=="C6" or n=="O6" or n=="N1" or n=="C2" or n=="N2" or n=="N3"

                elif nucleic3to1[r.type] == "C" :
                    a.isBase = n=="N1" or n=="C2" or n=="O2" or n=="N3" or n=="C4" or n=="N4" or n=="C5" or n=="C6"

                elif nucleic3to1[r.type] == "A" :
                    a.isBase = n=="N9" or n=="C8" or n=="N7" or n=="C5" or n=="C4" or n=="N3" or n=="C2" or n=="N1" or n=="C6" or n=="N6"

                elif nucleic3to1[r.type] == "U" :
                    a.isBase = n=="N1" or n=="C2" or n=="O2" or n=="N3" or n=="C4" or n=="O4" or n=="C5" or n=="C6"

                else :
                    #print " -x- NA res %d.%s is ?" % (r.id.position, r.type)
                    break

                a.isSC = a.isBase

                #if nucleic3to1[r.type] == "G" :
                #    r.isBase = n=="" or n=="" or n=="" or n=="" or n=="" or n=="" or n=="" or n=="" or n="" or n="" or n=""
                #    r.baseAt = r.atomsMap["N9"][0]

                if a.isBB :
                    r.bbAtoms.append ( a )
                else :
                    r.scAtoms.append ( a )

        else :
            for a in r.atoms :
                a.isBB, a.isSC, a.isSugar, a.isBase = False, False, False, False



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
