
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
    from segment_dialog import current_segmentation, segmentation_map
except :
    pass


# reload ( mapq_chimera.mapq ); mapq_chimera.mapq.show_dialog()

OML = chimera.openModels.list

devMenu = False
isModelZ = False

mapqVersion = "1.5.3"


dlgName = "mapqdlg"
dlgTitle = "MapQ (v%s)" % mapqVersion
dlgHelp = 'https://cryoem.slac.stanford.edu/ncmi/resources/software/mapq'

if isModelZ :
    devMenu = False
    dlgName = "modelzdlg"
    dlgTitle = "ModelZ (v1.2)"
    dlgHelp = 'https://cryoem.slac.stanford.edu/ncmi/resources/software/modelz'


atomColors = {'C' : chimera.MaterialColor (0.565,0.565,0.565),
            'Cbb' : chimera.MaterialColor (0.2,0.6,0.2),
            'S' : chimera.MaterialColor (1.000,1.000,0.188),
            'O' : chimera.MaterialColor (1.000,0.051,0.051),
            'N' : chimera.MaterialColor (0.188,0.314,0.973),
            'P' : chimera.MaterialColor (1.0, 0.502, 0.0),
            'H' : chimera.MaterialColor (0.9,.9,.9),
            ' ' : chimera.MaterialColor (0.2,1,.2)
            }


atomColors2 = {'C' : (0.565,0.565,0.565,1),
            'Cbb' : (0.2,0.6,0.2,1),
            'S' : (1.000,1.000,0.188,1),
            'O' : (1.000,0.051,0.051,1),
            'N' : (0.188,0.314,0.973,1),
            'P' : (1.0, 0.502, 0.0,1),
            'H' : (0.9,.9,.9,1),
            ' ' : (0.7,.9,.7)
            }


ac = { 'O' : chimera.MaterialColor( .9, .2, .2, 1.0 ),
        'C' : chimera.MaterialColor( .7, .7, .7, 1.0 ),
        'N' : chimera.MaterialColor( .2, .2, .9, 1.0 ),
        'H' : chimera.MaterialColor( 1, 1, 1, 1.0 ),
        ' ' : chimera.MaterialColor( .2, .2, .2, 1.0 ),
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
    buttons = ( "Close" )
    title = dlgTitle
    help = dlgHelp


    def fillInUI(self, parent):

        self.group_mouse_mode = None

        tw = parent.winfo_toplevel()
        self.toplevel_widget = tw
        tw.withdraw()

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
            ff.grid(column=0, row=row, sticky='w', pady=5, padx=10)

            l = Tkinter.Label(ff, text='Map:', anchor=Tkinter.W)
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

            if len ( self.cur_chains ) > 0 :
                self.chain.set ( self.cur_chains[0] )
                #self.ShowCh ( self.cur_chains[0] )
                self.GetSeq ()


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


        if 1 :

            row += 1
            ff = Tkinter.Frame(f)
            ff.grid(column=0, row=row, sticky='w', pady=0, padx=5)


            fff = Tkinter.Frame(ff, borderwidth=1, padx=2, pady=2, relief=Tkinter.GROOVE)
            fff.grid(column=10, row=0, sticky='e', pady=0, padx=5)

            l = Tkinter.Label(fff, text='Q-scores:', anchor=Tkinter.W)
            l.grid(column=1, row=0, sticky='w')


            #b = Tkinter.Button(fff, text="Sigma", command=self.CalcAllSigma )
            #b.grid (column=2, row=0, sticky='w', padx=5)

            #b = Tkinter.Button(fff, text="RadZ", command=self.CalcAllRadZ )
            #b.grid (column=4, row=0, sticky='w', padx=5)

            if isModelZ :

                b = Tkinter.Button(fff, text="Z-scores", command=self.CalcZScores )
                b.grid (column=5, row=0, sticky='w', padx=5)

            else :

                b = Tkinter.Button(fff, text="Calc", command=self.CalcAllQ )
                b.grid (column=2, row=0, sticky='w', padx=5)

                if 1 or devMenu :
                    b = Tkinter.Button(fff, text="Calc(MP)", command=self.CalcAllQp )
                    b.grid (column=3, row=0, sticky='w', padx=5)

                b = Tkinter.Button(fff, text="Load", command=self.GetQsFromFile )
                b.grid (column=4, row=0, sticky='w', padx=5)

                b = Tkinter.Button(fff, text="Sel", command=self.Q_sel2 )
                b.grid (column=5, row=0, sticky='w', padx=5)

            #b = Tkinter.Button(fff, text="R", command=self.CalcAllR )
            #b.grid (column=5, row=0, sticky='w', padx=5)

            #b = Tkinter.Button(fff, text="R", command=self.CalcAllR )
            #b.grid (column=5, row=0, sticky='w', padx=5)



            if 0 :

                self.colorMod = Tkinter.StringVar()
                self.colorMod.set ( 'sc' )

                b = Tkinter.Button(ff, text="Color:", command=self.DoColor)
                b.grid (column=20, row=0, sticky='w', padx=5)

                c = Tkinter.Radiobutton(ff, text="BB", variable=self.colorMod, value = 'bb')
                c.grid (column=21, row=0, sticky='w')

                c = Tkinter.Radiobutton(ff, text="SC", variable=self.colorMod, value = 'sc')
                c.grid (column=22, row=0, sticky='w')

                c = Tkinter.Radiobutton(ff, text="Rand", variable=self.colorMod, value = 'rand')
                c.grid (column=23, row=0, sticky='w')


            else :

                l = Tkinter.Label(ff, text=' Color:', fg="#777")
                l.grid(column=20, row=0, sticky='e')


                b = Tkinter.Button(ff, text="Bb", command=self.DoColorBB)
                b.grid (column=21, row=0, sticky='w', padx=5)

                b = Tkinter.Button(ff, text="SC", command=self.DoColorSC)
                b.grid (column=22, row=0, sticky='w', padx=5)

                b = Tkinter.Button(ff, text="Res", command=self.DoColorRes)
                b.grid (column=23, row=0, sticky='w', padx=5)

                if not isModelZ :
                    b = Tkinter.Button(ff, text="Atoms", command=self.DoColorAtoms)
                    b.grid (column=24, row=0, sticky='w', padx=5)

                b = Tkinter.Button(ff, text="Random", command=self.DoColorRandom)
                b.grid (column=25, row=0, sticky='w', padx=5)





            l = Tkinter.Label(ff, text='', fg="#000")
            l.grid(column=25, row=0, sticky='ens')


            ff = Tkinter.Frame(ff, borderwidth=1, padx=2, pady=2, relief=Tkinter.GROOVE)
            ff.grid(column=30, row=0, sticky='e', pady=0, padx=5)

            l = Tkinter.Label(ff, text='Sequence select: ', fg="#000")
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

            b = Tkinter.Button(ff, text="<", command=self.KeepBack)
            b.grid (column=42, row=0, sticky='w', padx=5)

            if 0 and devMenu :
                b = Tkinter.Button(ff, text="R", command=self.SelReLoad)
                b.grid (column=42, row=0, sticky='w', padx=5)

                b = Tkinter.Button(ff, text="L", command=self.SelLoad)
                b.grid (column=43, row=0, sticky='w', padx=5)



            #b = Tkinter.Button(ff, text="Clear", command=self.ClearSel)
            #b.grid (column=40, row=0, sticky='w', padx=5)

            #self.keepExMap = Tkinter.IntVar()
            #self.keepExMap.set(0)
            #oft = Tkinter.Checkbutton( ff, text="Keep Extracted Maps", variable=self.keepExMap, command=self.keepExMapCb)
            #oft.grid(column = 40, row = 0, sticky = 'w')


        self.mapRes = Tkinter.StringVar(f)
        self.mapRes.set ( "3" )


        if devMenu :

            row += 1
            ff = Tkinter.Frame(f)
            ff.grid(column=0, row=row, sticky='w', pady=0, padx=5)


            b = Tkinter.Label(ff, text="Res:")
            b.grid (column=1, row=0, sticky='w', padx=0, pady=5)

            e = Tkinter.Entry(ff, width=3, textvariable=self.mapRes)
            e.grid(column=2, row=0, sticky='w', padx=5, pady=5)


            b = Tkinter.Label(ff, text="   Sel:")
            b.grid (column=6, row=0, sticky='w', padx=0, pady=5)

            self.selText = Tkinter.StringVar(f)
            self.selText.set ( "" )
            e = Tkinter.Entry(ff, width=60, textvariable=self.selText)
            e.grid(column=7, row=0, sticky='w', padx=5, pady=5)


            b = Tkinter.Button(ff, text="Sel", command=self.SelText)
            b.grid (column=8, row=0, sticky='w', padx=5)


            b = Tkinter.Label(ff, text="Rad:")
            b.grid (column=10, row=0, sticky='w', padx=0, pady=5)

            self.maskRad = Tkinter.StringVar(f)
            self.maskRad.set ( "2.5" )
            e = Tkinter.Entry(ff, width=3, textvariable=self.maskRad)
            e.grid(column=11, row=0, sticky='w', padx=5, pady=5)


            b = Tkinter.Button(ff, text="AddSel", command=self.AdSel)
            b.grid (column=12, row=0, sticky='w', padx=5)


            b = Tkinter.Label(ff, text="   Atom:")
            b.grid (column=15, row=0, sticky='w', padx=0, pady=5)

            self.addText = Tkinter.StringVar(f)
            self.addText.set ( "Ca" )
            e = Tkinter.Entry(ff, width=10, textvariable=self.addText)
            e.grid(column=16, row=0, sticky='w', padx=5, pady=5)

            b = Tkinter.Button(ff, text="Add", command=self.AddAtom)
            b.grid (column=17, row=0, sticky='w', padx=5)


            b = Tkinter.Button(ff, text="Nr", command=self.ShowNear)
            b.grid (column=40, row=0, sticky='w', padx=5)

            b = Tkinter.Button(ff, text="Ds", command=self.ShowDists)
            b.grid (column=41, row=0, sticky='w', padx=5)

            b = Tkinter.Button(ff, text="Inter", command=self.Inter)
            b.grid (column=42, row=0, sticky='w', padx=5)


        dummyFrame = Tkinter.Frame(parent, relief='groove', borderwidth=1)
        Tkinter.Frame(dummyFrame).pack()
        dummyFrame.grid(row=row,column=0,columnspan=7, pady=3, sticky='we')
        row += 1


        global msg
        msg = Tkinter.Label(parent, width = 60, anchor = 'w', justify = 'left', fg="red", pady=5, padx=10)
        msg.grid(column=0, row=row, sticky='ew')
        self.msg = msg

        self.showingAtoms = False

        #umsg ( 'Select one or more segmented regions then press "Place Points" to start' )


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
        self.dmapMB.menu.delete ( 0, 'end' )   # Clear menu
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

            mlist = OML(modelTypes = [chimera.Molecule])
            for m in mlist :
                m.display = False

            mol.display = True

            self.cur_chains = self.GetChains ( mol )

            if len(self.cur_chains) == 0 :
                self.chain.set ( "" )
            elif self.chain.get() in self.cur_chains :
                print " - ch " + self.chain.get() + " already sel"
                self.ShowCh ( self.chain.get() )
            else :
                self.chain.set ( self.cur_chains[0] )
                self.ShowCh ( self.chain.get() )

            self.GetSeq ()
            self.ZoomBegin ()
            SetBBAts ( mol )



    def ChainSelected ( self, ch ) :
        print " - sel chain: ", ch, self.chain.get()
        self.ShowCh ( ch )
        self.GetSeq ()
        self.ZoomBegin ()



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
            self.chainMB.menu.add_radiobutton ( label=ch, variable=self.chain,
                                            command=lambda ch=ch: self.ChainSelected(ch) )

        self.chainMB.menu.add_radiobutton ( label="All", variable=self.chain,
                                        command=lambda ch=ch: self.ChainSelected(ch) )



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
        for sc in self.scores :
            if sc != None :
                foundScore = True

        if not foundScore :
            umsg ( "No scores... press Q, Qp, or Qf button first" )
            return


        minScore, maxScore = 0,0
        if colorMod == "sc" :
            minScore, maxScore = self.minSCscore, self.maxSCscore
        else :
            minScore, maxScore = self.minBBscore, self.maxBBscore

        cH = numpy.array( [0.0,1.0,0.0] )
        cL = numpy.array( [1.0,0.0,0.0] )

        for ri, r in enumerate ( self.seqRes ) :
            sc = None
            #sc = self.scores[ri] if colorSC else self.scores2[ri]
            if colorMod == "sc" :
                sc = r.scQ
            elif colorMod == "bb" :
                sc = r.bbQ
            else :
                sc = r.Q

            if sc == None  :
                r.ribbonColor = chimera.MaterialColor ( .7, .7, .7, 1.0 )
                for at in r.atoms :
                    #at.color = r.ribbonColor
                    try :
                        at.color = atomColors[at.name[0]]
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
                        at.color = atomColors[at.name[0]]
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

        #ct = {}
        #for r in self.cur_mol.residues: ct[r.id.chainId] = 1
        #clist = ct.keys()
        #clist.sort()

        for r in self.cur_mol.residues :
            if r.id.chainId == chainId :
                if ("CA" in r.atomsMap and "N" in r.atomsMap and "C" in r.atomsMap) or ("O3'" in r.atomsMap and "O5'" in r.atomsMap)  :
                    r.ribbonDisplay = True
                    r.ribbonDrawMode = 2
                else :
                    r.ribbonDisplay = False
                    for at in r.atoms :
                        at.drawMode = at.Ball
                        at.display = True
            else :
                if ("CA" in r.atomsMap and "N" in r.atomsMap and "C" in r.atomsMap) or ("O3'" in r.atomsMap and "O5'" in r.atomsMap)  :
                    r.ribbonDisplay = False
                    r.ribbonDrawMode = 2
                else :
                    r.ribbonDisplay = False
                    for at in r.atoms :
                        at.drawMode = at.Ball
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
                if res.isProt or res.isNA :
                    at.drawMode = at.EndCap
                    at.display = True # not showRibbon
                    if at.element.name in atomColors :
                        at.color = atomColors[at.element.name]
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



    def ShowAts ( self ) :

        for mod in chimera.openModels.list() :
            if type(mod) == chimera.Molecule and mod.display == True :

                #cid = "1"
                #rs = [520, 521, 635, 575, 298, 550, 525, 639, 551, 303, 547, 305, 519]

                cid = "4"
                rs = [38, 42, 242, 244, 246, 181, 182, 135, 251, 94, 98, 91, 95, 284]

                #cid = "E"
                #rs = [128, 33, 136]

                for res in mod.residues :
                    #if res.id.position in rs and res.id.chainId == cid :
                    if res.id.position in rs :
                        for at in res.atoms :
                            at.drawMode = at.EndCap
                            at.display = True
                            try :
                                at.color = atomColors[at.name[0]]
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
                at.display = at.isBB
                if at.residue.isNA :
                    at.display = at.isBB and not at.isSugar

                #try :
                #    at.color = atomColors[at.name[0]]
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
                    at.color = atomColors[at.name[0]]
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


        nearRes = {}
        for r in ress :
            nats = self.AtsWithin ( r.atoms, 4.0, allAtTree )
            for at in nats :
                nearRes[at.residue] = 1

        for r in nearRes.keys() :
            print " -- %s.%d.%s - %d atoms" % (r.type, r.id.position, r.id.chainId, len(r.atoms))
            for at in r.atoms :
                #at.drawMode = at.EndCap
                at.display = True
                if at.name in atomColors :
                    at.color = atomColors[at.name[0]]









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

        m = self.cur_mol
        print " - cur mol:", m.name

        ct = {}
        for r in m.residues: ct[r.id.chainId] = 1
        clist = ct.keys()
        clist.sort()

        atsMap = {}
        for r in m.residues :
            show = True if r.id.chainId == ch else False
            if ("CA" in r.atomsMap and "N" in r.atomsMap and "C" in r.atomsMap) or ("O3'" in r.atomsMap and "O5'" in r.atomsMap) :
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
                        if at.element.name in atomColors :
                            at.color = atomColors[at.element.name]
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

        try :
            print self.cur_mol.name
        except :
            print " - mol may have been closed"
            return

        self.GetSeqFromStruc ( self.cur_mol, self.chain.get() )

        if len(self.seq) > 0 :

            print "-- seq from open mol -- %d res" % len(self.seq)
            print self.seq

            self.seqt = []
            self.seqSheetR = [None] * len(self.seq)
            self.seqHelixR = [None] * len(self.seq)
            self.seqScoreR = [None] * len(self.seq)
            self.seqScoreR2 = [None] * len(self.seq)
            self.scores2 = [None] * len(self.seq)
            self.scores = [None] * len(self.seq)

            self.UpdateSeqFont ()

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
            del self.seqText

        self.seqSel = None
        self.seq = ""
        self.UpdateSeqSel ()



    def GetSeqFromStruc ( self, mol, chainId ) :

        print "Getting seq from %s, %s" % (mol.name, chainId)

        self.conf = ""
        self.pred = ""
        self.seq = ""
        self.seqRes = []

        from chimera.resCode import protein3to1
        from chimera.resCode import nucleic3to1
        protein3to1['HSD'] = protein3to1['HIS']

        rids = {}
        for r in mol.residues :
            if r.id.chainId == chainId :
                if r.type in protein3to1 or r.type in nucleic3to1 :
                    rids[r.id.position] = r


        ris = rids.keys()
        ris.sort()

        for ri in ris :
            r = rids[ri]
            if r.type in protein3to1 :
                self.seq = self.seq + protein3to1[r.type]
                self.conf = self.conf + "9"
                self.predi = "C"
                if r.isSheet :
                    self.predi = "E"
                if r.isHelix :
                    self.predi = "H"
                self.pred = self.pred + self.predi
                self.seqRes.append ( r )
            elif r.type in nucleic3to1 :
                self.seq = self.seq + nucleic3to1[r.type]
                self.conf = self.conf + "9"
                self.predi = "C"
                self.pred = self.pred + self.predi
                self.seqRes.append ( r )




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

        self.minSCscore, self.maxSCscore = 0,2
        self.minBBscore, self.maxBBscore = 0,4

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


        print "Average Sigma sc : %.2f - %.2f, avg %.2f | %.2f - %.2f, avg %.2f" % (scMin, scMax, scAvg, 1.0/scMin, 1.0/scMax, 1.0/scAvg)
        print "Average Sigma bb : %.2f - %.2f, avg %.2f | %.2f - %.2f, avg %.2f" % (bbMin, bbMax, bbAvg, 1.0/bbMin, 1.0/bbMax, 1.0/bbAvg)


        self.minSCscore, self.maxSCscore = 0.0,0.5
        self.minBBscore, self.maxBBscore = 0.0,0.2

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



        CalcQ ( self.cur_mol, self.chain.get(), self.cur_dmap, allAtTree=allAtTree, log=True )
        SaveQ ( self.cur_mol, self.chain.get(), self.cur_dmap, float(self.mapRes.get()) )


        self.scores, self.scores2 = [], []
        scBB, scSC = [], []

        for r in self.cur_mol.residues :
            if cid == None or r.id.chainId == cid :
                if r.isProt or r.isNA :
                    self.scores2.append ( r.bbQ )
                    self.scores.append ( r.scQ )
                    if r.bbQ != None :
                        scBB.append ( r.bbQ )
                    if r.scQ != None :
                        scSC.append ( r.scQ )
                else :
                    self.scores2.append ( r.Q )
                    self.scores.append ( r.Q )


        #bbRes = numpy.power ( numpy.e, (self.avgScore2 - 8.0334) / -4.128 ) # y = -4.128ln(x) + 8.0334
        #scRes = numpy.power ( numpy.e, (self.avgScore - 4.8261) / -3.097 ) # y = -3.097ln(x) + 4.8261
        #scRes = (self.avgScore2 - 3.507) / -0.721
        #bbRes = (self.avgScore - 6.1234) / -0.9191


        try :
            scMin, scMax, scAvg = min(scSC), max(scSC), numpy.average(scSC)
            bbMin, bbMax, bbAvg = min(scBB), max(scBB), numpy.average(scBB)


            print "Average Q sc : %.2f - %.2f, avg %.2f" % (scMin, scMax, scAvg )
            print "Average Q bb : %.2f - %.2f, avg %.2f" % (bbMin, bbMax, bbAvg )


            self.minSCscore, self.maxSCscore = 0.0,1
            self.minBBscore, self.maxBBscore = 0.0,1

            self.UpdateSeq ()

        except :
            pass





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



        CalcQp (self.cur_mol, cid, self.cur_dmap, allAtTree=None )



        self.scores, self.scores2 = [], []
        scBB, scSC = [], []

        for r in self.cur_mol.residues :
            if cid == None or r.id.chainId == cid :
                if r.isProt or r.isNA :
                    self.scores2.append ( r.bbQ )
                    self.scores.append ( r.scQ )
                    if r.bbQ != None :
                        scBB.append ( r.bbQ )
                    if r.scQ != None :
                        scSC.append ( r.scQ )
                else :
                    self.scores2.append ( r.Q )
                    self.scores.append ( r.Q )

        scMin, scMax, scAvg = min(scSC), max(scSC), numpy.average(scSC)
        bbMin, bbMax, bbAvg = min(scBB), max(scBB), numpy.average(scBB)


        print " - Average Q sc : %.2f - %.2f, avg %.2f" % (scMin, scMax, scAvg )
        print " - Average Q bb : %.2f - %.2f, avg %.2f" % (bbMin, bbMax, bbAvg )


        self.minSCscore, self.maxSCscore = 0.0,1
        self.minBBscore, self.maxBBscore = 0.0,1


        self.UpdateSeq ()
        SaveQ ( self.cur_mol, self.chain.get(), self.cur_dmap, float(self.mapRes.get()) )




    def Q_sel2 (self) :

        # show sigma for a side chain

        atoms = chimera.selection.currentAtoms()
        if len ( atoms ) == 0 :
            umsg ( "No selected atoms found" )
            return

        dmap = self.cur_dmap
        mol = atoms[0].molecule


        #selAtom = selAts[0]
        #r = selAtom.residue
        #print "Res: %s - %d.%s - %s - Atom: %s" % (r.type, r.id.position, r.id.chainId, r.molecule.name, selAtom.name)

        print " - in map: %s" % dmap.name
        print " - mol: %s" % mol.name

        if 1 or not hasattr ( mol.name, 'bbats' ) :
            SetBBAts(mol)
            mol.bbats = True

        ats = [at for at in mol.atoms if not at.element.name == "H"]
        points = _multiscale.get_atom_coordinates ( ats, transformed = False )
        print " - search tree: %d/%d ats" % ( len(ats), len(mol.atoms) )
        allAtTree = AdaptiveTree ( points.tolist(), ats, 1.0)

        sigma = 0.6
        minD, maxD = MinMaxD ( dmap )

        import time
        start = time.time()


        from chimera import tasks, CancelOperation
        task = tasks.Task('Calculating Q-scores', modal = True)

        avg = 0.0

        import traceback


        try :

            for ai, at in enumerate ( atoms ) :

                cc, ccm = RadCC ( [at], dmap, sigma, allAtTree=allAtTree, show=0, log=0, numPts=8, toRAD=2.0, dRAD=0.1, minD=minD, maxD=maxD, fitg=0 )

                at.Q = ccm
                at.occupancy = ccm
                avg += ccm

                end = time.time()
                totSec = end - start

                leftTime = ""
                leftSec = 0.0
                iPerSec = float(ai) / totSec
                if iPerSec > 0 :
                    leftSec = float ( len(atoms) - ai ) / iPerSec
                    leftHour = numpy.floor ( leftSec / 60.0 / 60.0 )
                    leftSec = leftSec - leftHour * 60.0 * 60.0
                    leftMin = numpy.floor ( leftSec / 60.0 )
                    leftSec = leftSec - leftMin * 60.0
                    leftTime = "%.0f:%.0f:%.0f" % (leftHour, leftMin, leftSec)

                if (ai+1) % 10 == 0 :
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


        #cc = ResCC ( mol, atoms, 4.0, dmap )
        #print " - CC: %.3f" % cc


        avgq = avg / float(len(atoms))
        if len(atoms) > 1 :
            umsg ( "Q-score of %d atoms: %.2f" % (len(atoms), avgq) )
        else :
            umsg ( "Q-score of %d atom: %.2f" % (len(atoms), avgq) )

        for at in atoms :
            print "%s %d.%s %s : %.2f" % (at.residue.type, at.residue.id.position, at.residue.id.chainId, at.name, at.Q)




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
                                at.occupancy = at.Q = occ
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

        QStats1 (self.cur_mol, chainId)


        if 1 :
            SaveQ ( self.cur_mol, chainId, self.cur_dmap, float(self.mapRes.get()) )



        self.scores, self.scores2 = [], []
        scBB, scSC = [], []

        doRess = []
        for r in self.cur_mol.residues :
            if r.id.chainId == chainId :
                doRess.append ( r )

        print "Q for %d res..." % ( len(doRess) )
        for r in doRess :

            CalcResQ (r, None, None, useOld=True )

            if r.isProt or r.isNA :
                self.scores2.append ( r.bbQ )
                self.scores.append ( r.scQ )
                if r.bbQ != None : scBB.append ( r.bbQ )
                if r.scQ != None : scSC.append ( r.scQ )
            else :
                self.scores2.append ( r.Q )
                self.scores.append ( r.Q )


        if len (scSC) > 0 :
            scMin, scMax, scAvg = min(scSC), max(scSC), numpy.average(scSC)
            bbMin, bbMax, bbAvg = min(scBB), max(scBB), numpy.average(scBB)

            print " - average Q sc : %.2f - %.2f, avg %.2f" % (scMin, scMax, scAvg )
            print " - average Q bb : %.2f - %.2f, avg %.2f" % (bbMin, bbMax, bbAvg )
            print ""

            self.minSCscore, self.maxSCscore = 0.0,1
            self.minBBscore, self.maxBBscore = 0.0,1

            self.UpdateSeq ()


        #self.QStats ()
        #self.QStatsRNA()










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







    def UpdateSeqFont ( self ) :
        # http://stackoverflow.com/questions/4296249/how-do-i-convert-a-hex-triplet-to-an-rgb-tuple-and-back

        if not hasattr ( self, 'seq' ) :
            print " - update seq font - no seq"
            return

        print "seq len %d, text w %d" % ( len(self.seq), self.tw )

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

        if hasattr ( self, 'seqText' ) :
            self.Canvas.coords ( self.seqText, x_at, y_at )
            self.Canvas.itemconfigure ( self.seqText, font=self.font )
        else :
            self.seqText = self.Canvas.create_text( x_at, y_at, text=self.seq, font=self.font, anchor='w')


        #self.UpdateSeqSel ()




    def UpdateSeq ( self ) :

        if not hasattr ( self, 'seq' ) :
            print " - update seq - no seq"
            return

        x_at = self.seqX
        y_at = self.seqY + self.seqH/2

        if hasattr ( self, 'seqText' ) :
            self.Canvas.coords ( self.seqText, x_at, y_at )
        else :
            self.seqText = self.Canvas.create_text( x_at, y_at, text=self.seq, font=self.font, anchor='w')

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

                sc = None
                try :
                    sc = self.scores[si]
                except :
                    #continue
                    pass

                if sc == None :
                    if self.seqScoreR[si] != None :
                        self.Canvas.delete ( self.seqScoreR[si] )
                    self.seqScoreR[si] = None
                else :
                    xx0 = self.seqX + si * self.tw + 2
                    xx1 = xx0 + self.tw - 2
                    h = (sc - self.minSCscore) / (self.maxSCscore - self.minSCscore)
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

                bb = None
                try :
                    bb = self.scores2[si]
                except :
                    #continue
                    pass

                if bb == None :
                    if self.seqScoreR2[si] != None :
                        self.Canvas.delete ( self.seqScoreR2[si] )
                    self.seqScoreR2[si] = None
                else :
                    xx0 = self.seqX + si * self.tw + 2
                    xx1 = xx0 + self.tw - 2
                    h = (bb - self.minBBscore) / (self.maxBBscore - self.minBBscore)
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
            status ( "Start sequence sel at %d" % (seqI+1) )
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
        print "b1 up - ctrl - ", event.x, event.y
        self.B1_Up ( event )


    def B1_Up_Shift ( self, event ) :
        print "b1 up - shift - "
        self.B1_Up ( event )

    def B1_Up_Alt ( self, event ) :
        print "b1 up - alt - "
        self.B1_Up ( event )


    def B1_Up (self, event) :
        print "b1 up - ", event.x, event.y

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
                    startI = self.seqRes [ max(self.seqSel[0],0) ].id.position
                    endI = self.seqRes [ min(self.seqSel[1],len(self.seqRes)-1) ].id.position
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
        print "mod: ", self.modX, " seq:", self.seqX


    def KeepBack ( self ) :
        print " - keep - remove last..."

        if hasattr ( self, 'prevSel' ) and len(self.prevSel) > 0 :
            self.prevSel.pop()

            chimera.selection.clearCurrent()

            for s in self.prevSel :
                print " -s- adding to sel:", s
                chimera.selection.mergeCurrent ( chimera.selection.EXTEND, chimera.selection.OSLSelection (s) )

            if self.selExtract.get () :
                self.ShowSel ()



    def SelReLoad ( self ) :

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






    def ShowSel ( self ) :

        #showRibbon = self.showRibbon.get()
        showRibbon = not self.showingAtoms # self.showRibbon.get()
        showLigands = self.showLigands.get()
        showSC = True # self.showAtoms.get()

        atoms = []
        scores = []
        selResM = {}
        for r in chimera.selection.currentResidues () :
            rid = "%d.%s" % (r.id.position, r.id.chainId)
            selResM [rid] = 1

        if self.cur_mol == None :
            return

        if 1 or not hasattr ( self.cur_mol, 'bbats' ) :
            SetBBAts(self.cur_mol)
            self.cur_mol.bbats = True


        for r in self.cur_mol.residues :
            rid = "%d.%s" % (r.id.position, r.id.chainId)
            if rid in selResM :

                if hasattr (r, 'scZ') and r.scZ != None :
                    scores.append(r.scZ)

                r.ribbonDisplay = showRibbon

                for at in r.atoms :
                    if at.element.name == "H" :
                        at.display = False
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
                            at.color = atomColors[at.element.name]
                            #if at.element.name == "C" :
                            #    at.color = atomColors['Cbb']
                        else :
                            at.color = atomColors[at.element.name]

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
                            if at.element.name in atomColors :
                                at.color = atomColors[at.element.name]
                            atoms.append ( at )
                            ligAts.append ( at )
                    else :
                        for at in r.atoms :
                            at.display = False

            #chimera.selection.clearCurrent ()
            chimera.selection.addCurrent ( ligAts )


        #for bond in self.seqRes[0].molecule.bonds :
        #    bond.display = bond.Smart
            #if bond.atoms[0] in atMap and bond.atoms[1] in atMap :
            #    #bond.display = bond.Smart
            #    bond.display = bond.Smart
            #else :
            #    #bond.display = bond.Never
            #    bond.display = bond.Smart


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
                sp.color = (0.7, 0.7, 0.7, 0.2)


    def PtsToMap ( self, points, dmap, atomRad, nname, showMesh = False ) :

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
        nv.surface_levels[0] = dmap.surface_levels[0]
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
                    sp.color = (0.7, 0.7, 0.7, 0.1)


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
                status ( "Sequence selected %d - %d" % (self.seqSel[0]+1, self.seqSel[1]+1) )
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
                    resEnd = self.seqRes [ len(self.seqRes) - 1 ]

                    try :
                        status ( "Sequence: %s/%s %d/%d" % ( self.seq[si], res.type, res.id.position, resEnd.id.position ) )
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

            #self.seqX += event.delta * 10

            self.mag = self.mag + event.delta
            if self.mag > 15 : self.mag = 15
            if self.mag < 2 : self.mag = 2

            self.font = tkFont.Font(family='Courier', size=(self.mag), weight='normal')
            #self.boldFont = tkFont.Font(family='Courier', size=(self.mag+4), weight='bold')
            self.tw = self.font.measure ( "a" )

            #GetSegMod().seqX = self.seqX
            self.UpdateSeqFont ()
            self.UpdateSeq ()

            # ['__doc__', '__module__', 'char', 'delta', 'height', 'keycode', 'keysym', 'keysym_num', 'num', 'send_event', 'serial', 'state', 'time', 'type', 'widget', 'width', 'x', 'x_root', 'y', 'y_root']
            #print dir(event)
            #print event.delta
            status ( "Mag: %d" % self.mag )



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
        self.UpdateSeq ()

    def ZoomEnd ( self ) :
        self.seqX = - ( len(self.seq) - 50 ) * self.tw
        self.UpdateSeq ()



    def isSelected ( self, fmap ) :
        for sp in fmap.surfacePieces :
            if sp in Surface.selected_surface_pieces() :
                return True
        return False














def AddSpherePts ( pts, clr, rad, mname = "RAD points" ) :

    from chimera import elements, Coord, Atom, MolResId

    ptsMol = GetMod ( mname )

    res = None
    if ptsMol == None:
        from chimera import Molecule, openModels
        ptsMol = Molecule()
        ptsMol.name = mname
        ptsMol.isRealMolecule = False
        openModels.add ( [ptsMol], noprefs = True )
        res = ptsMol.newResidue('marker', chimera.MolResId('1', 1) )
    else :
        res = ptsMol.residues[0]

    for pt in pts :
        a = ptsMol.newAtom('', elements.H)
        res.addAtom(a)

        a.setCoord ( chimera.Point(*pt) )  # ( chimera.Point(*xyz) )
        a.radius = rad
        a.drawMode = Atom.Sphere
        a.color = chimera.MaterialColor ( *clr )
        a.surfaceCategory = 'markers'



def SpherePts ( ctr, rad, N ) :

    thetas, phis = [], []
    from math import acos, sin, cos, sqrt, pi
    for k in range ( 1, N+1 ) :
        h = -1.0 + ( 2.0*float(k-1)/float(N-1) )
        phis.append ( acos(h) )
        thetas.append ( 0 if k == 1 or k == N else
                        (thetas[k-2] + 3.6/sqrt(N*(1.0-h**2.0))) % (2*pi) )

    pts = [None] * N
    for i, theta, phi in zip ( range(N), thetas, phis ):
        v = chimera.Vector (sin(phi)*cos(theta), sin(phi)*sin(theta), cos(phi))
        #if numpy.abs ( v.length - 1.0 ) > 1e-3 :
        #    print "x"
        pt = ctr + v * rad
        pts[i] = pt

    return pts



import threading


def Calc_ ( label="", res=0.0, MP=True ) :

    print "Calc Q scores:", label

    from VolumeViewer import Volume
    vols = chimera.openModels.list(modelTypes = [Volume])
    if len(vols) == 0 :
        print " - no volumes loaded"
        return
    dmap = vols[0]
    print " - dmap: %s" % dmap.name
    print " - res: %s" % res

    #fp = open ( "/Users/greg/_data/_mapsq/scores.txt", "a" )
    #fp.write ( "%s...\n" % dmap.name.split("_")[0]  )
    #fp.close ()

    from chimera import Molecule
    mols = chimera.openModels.list(modelTypes = [Molecule])
    if len(mols) == 0 :
        print " - no molecules loaded"
        return
    mol = mols[0]
    print " - mol: %s" % mol.name
    SetBBAts ( mol )

    ats = [at for at in mol.atoms if not at.element.name == "H"]
    points = _multiscale.get_atom_coordinates ( ats, transformed = False )
    print " - search tree: %d/%d ats" % ( len(ats), len(mol.atoms) )
    #allAtTree = AdaptiveTree ( points.tolist(), ats, 1.0)
    allAtTree = None


    cc, ccm, dr, q, qcc, emr = 0,0,0,0,0,0
    #bbRadZ, scRadZ, scRotaZ = 0,0,0

    if 0 :

        #cc, ccm, dr, ccr, ccmr = CalcSCBBr ( mol, mol.residues[0].id.chainId, dmap )
        cc, ccm, dr, ccr, ccmr = CalcSCBBr ( mol, None, dmap )

    if 1 :
        #bbSig, scSig = CalcSigma ( mol, mol.residues[0].id.chainId, dmap, allAtTree, useOld=False, log=False )
        #bbRadZ, scRadZ = CalcRadZ ( mol, mol.residues[0].id.chainId, dmap, allAtTree, useOld=False, log=False )

        #q, qcc = CalcQp ( mol, mol.residues[0].id.chainId, dmap, allAtTree=allAtTree )
        q, qcc = CalcQp ( mol, None, dmap, allAtTree=allAtTree )

        #q, qcc = CalcQ ( mol, None, dmap, allAtTree=allAtTree )
        #q, qcc = CalcQp ( mol, None, dmap, allAtTree=allAtTree )

        print ""
        print "Avg. Q score: %.3f" % q
        print ""

    if 0 :
        bbRadZ, scRadZ = CalcRadZ ( mol, None, dmap, allAtTree, useOld=False, log=False )

    if 0 :
        print 'Side Chain Rota-Z for %d ress' % len(mol.residues)
        Zs = CalcRotaZ ( dmap, mol, mol.residues )
        scRotaZ = numpy.average ( Zs )

    if 0 :
        emr = emringer (dmap, mol)

    if 0 :
        at = 14
        fp = None
        if os.path.isdir("/Users/greg/Dropbox/_mapsq") :
            fp = open ( "/Users/greg/Dropbox/_mapsq/scores%d_Q_allc_%s.txt" % (at, label), "a" )
        elif os.path.isdir("/home/greg/Dropbox/_mapsq") :
            fp = open ( "/home/greg/Dropbox/_mapsq/scores%d_Q_allc_%s.txt" % (at, label), "a" )
        else :
            fp = open ( "scores%d_Q_allc_%s.txt" % (at, label), "a" )

        #fp.write ( "%s\t%s\t%f\t%f\t%f\t%f\t%f\t%f\n" % (dmap.name.split(".")[0], mol.name.split(".")[0], cc, ccm, dr, q, qcc, emr)  )
        #fp.write ( "%s\t%f\n" % (dmap.name.split("_")[0], scRotaZ)  )
        #fp.write ( "%s\t%s\t%f\t%f\t%f\t%f\t%f\t%f\n" % (dmap.name.split(".")[0], mol.name.split(".")[0], cc, ccm, dr, q, qcc, emr)  )
        #fp.write ( "%s\t%s\t%s\t%f\t%f\t%f\n" % (dmap.name, mol.name, res, q, cc, ccm)  )

        nProt = len ( [at for at in mol.atoms if at.residue.isProt == True] )
        nNA = len ( [at for at in mol.atoms if at.residue.isNA == True] )

        fp.write ( "%s\t%s\t%s\t%d\t%d\n" % (dmap.name, mol.name, res, nProt, nNA)  )

        fp.close ()




def emringer ( dmap, mol ) :

    print "----- %s ____________ EMRINGER ____________ %s -----" % (dmap.name, mol.name)

    cdir = os.getcwd()
    print " - now in: ", cdir

    #print " - splitting " + mol.openedAs[0]
    mpath, mname = os.path.split ( mol.openedAs[0] )
    dpath, dname = os.path.split ( dmap.data.path )

    bs = os.path.splitext ( mol.openedAs[0] )[0]


    print " - copying mol file... removes symmetry/connect stuff"
    fin = open ( mol.openedAs[0], "r" )
    fout = open ( bs + "_.pdb", "w" )
    for line in fin :
        if "ATOM" in line or "HETATM" in line :
            fout.write ( line )
    fin.close ()
    fout.close ()


    phPath = "/Users/greg/_mol/phenix-1.14-3260/build/bin/"
    #phPath = "/Users/greg/_mol/phenix-1.15rc3-3435/build/bin/"

    args = [phPath+'phenix.emringer', dmap.data.path, bs+"_.pdb" ]
    print "running: ",
    for arg in args : print arg,
    print ""

    outf = mpath + '/' + '_out.txt'
    errf = mpath + '/' + '_err.txt'
    fout = open ( outf, "w" )
    ferr = open ( errf, "w" )
    import subprocess
    p = subprocess.Popen(args, stdout=fout, stderr=ferr, cwd=mpath)
    p.wait()
    fout.close()
    ferr.close()

    print " - getting score from " + outf
    score = -100
    fin = open ( outf )
    for l in fin :
        if "EMRinger Score:" in l :
            s = l [ len("EMRinger Score:")+1 : ]
            print "Score: ", s
            score = float( s )
            print " - found score: %.3f" % score

    print " - removing ", bs + "_.pdb"
    import shutil
    try :
        os.remove ( bs + "_.pdb" )
        os.remove ( bs + "__emringer.pkl" )
        os.remove ( bs + "__emringer.csv" )
        shutil.rmtree ( bs + "__emringer_plots" )
        print " - done"
    except :
        print "  -- did not find"

    return score


def refine ( dmap, mol, res ) :

    print "----- %s ____________ REFINE ____________ %s -----" % (dmap.name, mol.name)

    cdir = os.getcwd()
    print " - now in: ", cdir

    #print " - splitting " + mol.openedAs[0]
    mpath, mname = os.path.split ( mol.openedAs[0] )
    dpath, dname = os.path.split ( dmap.data.path )

    bs = os.path.splitext ( mol.openedAs[0] )[0]


    print " - copying mol file... removes symmetry/connect stuff"
    fin = open ( mol.openedAs[0], "r" )
    fout = open ( bs + "_.pdb", "w" )
    for line in fin :
        if "ATOM" in line or "HETATM" in line :
            fout.write ( line )
    fin.close ()
    fout.close ()


    phPath = "/Users/greg/_mol/phenix-1.14-3260/build/bin/"
    phPath = "/Users/greg/_mol/phenix-1.15rc3-3435/build/bin/"

    args = [phPath+'phenix.real_space_refine', dmap.data.path, bs+"_.pdb", "resolution=%.1f"%res ]
    print "running: ",
    for arg in args : print arg,
    print ""

    outf = mpath + '/' + '_out.txt'
    errf = mpath + '/' + '_err.txt'
    fout = open ( outf, "w" )
    ferr = open ( errf, "w" )
    import subprocess
    p = subprocess.Popen(args, stdout=fout, stderr=ferr, cwd=mpath)
    p.wait()
    fout.close()
    ferr.close()

    print " - getting score from " + outf
    score = -100
    fin = open ( outf )
    for l in fin :
        if "EMRinger Score:" in l :
            s = l [ len("EMRinger Score:")+1 : ]
            print "Score: ", s
            score = float( s )
            print " - found score: %.3f" % score

    print " - removing ", bs + "_.pdb"
    import shutil
    try :
        os.remove ( bs + "_.pdb" )
        os.remove ( bs + "__emringer.pkl" )
        os.remove ( bs + "__emringer.csv" )
        shutil.rmtree ( bs + "__emringer_plots" )
        print " - done"
    except :
        print "  -- did not find"

    return score


def refdir ( rdir ) :

    print "Refining in", rdir



def CalcR_ ( label = "" ) :

    print "Calc all scores -", label

    from VolumeViewer import Volume
    dmap = chimera.openModels.list(modelTypes = [Volume])[0]
    print " - dmap: %s" % dmap.name

    #fp = open ( "/Users/greg/_data/_mapsq/scores.txt", "a" )
    #fp.write ( "%s...\n" % dmap.name.split("_")[0]  )
    #fp.close ()

    from chimera import Molecule
    mol = chimera.openModels.list(modelTypes = [Molecule])[0]
    print " - mol: %s" % mol.name
    SetBBAts ( mol )


    mapName = os.path.splitext(dmap.name)[0]
    molName = os.path.splitext(mol.name)[0]
    ddir, dfile = os.path.split(dmap.data.path)

    molFile = mol.openedAs[0]
    mdir, mfile = os.path.split(molFile)

    print "PhFmap -- " + molFile

    RES = 3.0
    print " -- res %.1f -- " % RES

    outFile = molFile + "_r%.0f" % RES + "_fmodel.ccp4"

    if not os.path.isfile ( outFile ) :

        phPath = "/usr/local/phenix-1.14-3260/build/bin/"

        args = [phPath+'phenix.fmodel', "high_resolution=%.1f"%RES, "scattering_table=electron", "generate_fake_p1_symmetry=True", molFile ]
        print "running: ",
        for arg in args : print arg,
        print ""

        fout = open ( mdir + '/' + '_0_fmodel.log', "w" )
        import subprocess
        p = subprocess.Popen(args, stdout=fout, cwd=mdir)
        p.wait()
        fout.close()

        print ""
        args = [phPath+'phenix.mtz2map', "high_resolution=%.1f"%RES, "include_fmodel=true", "scattering_table=electron", molFile, molFile + ".mtz" ]
        print "running: ",
        for arg in args : print arg,
        print ""

        fout = open ( mdir + '/' + '_1_mtz2map.log', "w" )
        p = subprocess.Popen(args, stdout=fout, cwd=mdir)
        p.wait()
        fout.close()

        print " - renaming to:", outFile
        os.rename( molFile + "_fmodel.ccp4", outFile )
        os.remove( molFile + ".mtz" )


    print " - loading map:", outFile
    dm = chimera.openModels.open ( outFile )[0]



    molg = MyMolMapX ( mol, mol.atoms, RES, dmap.data.step[0], chimera.Xform.identity() )
    fpoints, fpoint_weights = fit_points_g ( molg, 0.1 )
    map_values = dmap.interpolated_values ( fpoints, mol.openState.xform )

    mmolap, mmcorr1, mmcorr2 = FitMap.overlap_and_correlation ( fpoint_weights, map_values )
    print "Molmap - olap: %f, CC: %f, CCm: %f" % (mmolap, mmcorr1, mmcorr2)

    fpoints, fpoint_weights = fit_points_g ( dm.data, 5.0 )
    map_values = dmap.interpolated_values ( fpoints, dm.openState.xform )
    olap, phcorr1, phcorr2 = FitMap.overlap_and_correlation ( fpoint_weights, map_values )
    print "Phenix - olap: %f, CC: %f, CCm: %f" % (olap, phcorr1, phcorr2)

    #fpoints, fpoint_weights = fit_points_g ( dmap.data, -1e6 )
    #map_values = dm.interpolated_values ( fpoints, dmap.openState.xform )
    #olap, corr1, corr2 = FitMap.overlap_and_correlation ( fpoint_weights, map_values )
    #print "Phenix - olap: %f, CC: %f, CCm: %f" % (olap, corr1, corr2)


    print "%f\t%f\t%f\t%f" % (mmcorr1, mmcorr2, phcorr1, phcorr2)

    fp = open ( "/Users/greg/Dropbox/_mapsq/scores3_R_%s.txt" % label, "a" )
    fp.write ( "%s\t%f\t%f\t%f\t%f\n" % (dmap.name.split("_")[0], mmcorr1, mmcorr2, phcorr1, phcorr2)  )
    fp.close ()








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





def CalcRadZ ( mol, cid, dmap, allAtTree, useOld=False, log=False ) :


    print "Rad-Z Scores"
    print " - map: %s" % dmap.name
    print " - mol: %s, chain: %s" % (mol.name, cid if cid != None else "_all_")


    ress = []
    for r in mol.residues :
        if cid == None or r.id.chainId == cid :
            if not useOld :
                ress.append ( r )
            elif not hasattr (r, 'scS' ) :
                ress.append ( r )

    print " - residues to do: %d" % len(ress)


    for ri, r in enumerate ( ress ) :

        r.scZ = RadZ ( r.scAtoms, dmap, allAtTree=allAtTree, show=0, log=0, numPts=10, toRAD=2 )
        r.bbZ = RadZ ( r.bbAtoms, dmap, allAtTree=allAtTree, show=0, log=0, numPts=10, toRAD=2 )

        if log and ri % 10 == 0 :
            status ( "Calculating - res %d/%d" % (ri, len(ress)) )
            print ".",


    scoresBB, scoresSC = [], []
    for r in mol.residues :
        if cid == None or r.id.chainId == cid :
            if r.bbZ != None :
                scoresBB.append ( r.bbZ )
            if r.scZ != None :
                scoresSC.append ( r.scZ )

    print " - avg radz - side chain %.1f, backbone %.1f" % (numpy.average(scoresSC), numpy.average(scoresBB) )

    return numpy.average(scoresBB), numpy.average(scoresSC)




def qwork (num, ress, dmap, allAtTree, log):

    print 'qwork %d - %d res, %d - %d' % (num, len(ress), ress[0].id.position, ress[-1].id.position)

    for ri, r in enumerate ( ress ) :
        r.scZ = RadAts ( r.scAtoms, dmap, allAtTree=allAtTree, show=0, log=0, numPts=10, toRAD=2, dRAD=0.2 )
        r.bbZ = RadAts ( r.bbAtoms, dmap, allAtTree=allAtTree, show=0, log=0, numPts=10, toRAD=2, dRAD=0.2 )

        if num == 0 and log :
            status ( "Calculating Q scores - %d/%d" % (ri, len(ress)) )
            print ".",



def CalcSigma ( mol, cid, dmap, allAtTree, useOld=False, log=False ) :


    print "Sigma Scores"
    print " - map: %s" % dmap.name
    print " - mol: %s, chain: %s" % (mol.name, cid if cid != None else "_all_")

    ress = []
    for r in mol.residues :
        if cid == None or r.id.chainId == cid :
            if not useOld :
                ress.append ( r )
            elif not hasattr (r, 'scS' ) :
                ress.append ( r )

    print " - residues to do: %d" % len(ress)



    if 0 :

        import multiprocessing, threading
        N = 4 # multiprocessing.cpu_count()
        print " - cores: %d" % N
        dn = len(ress) / N

        threads = []
        for i in range(N):
            l = i * dn
            h = (i+1)*dn if i != N-1 else len(ress)
            #print "t %d, %d-%d" % (i, l, h)

            #t = threading.Thread(target=qwork, args=(i,ress[l:h], dmap, allAtTree))
            #threads.append(t)
            #t.start()

            #t = threading.Thread(name='d%d'%i, target=qwork, args=(i,ress[l:h], dmap, allAtTree, log))
            #t.setDaemon(True)
            #t.start()
            #threads.append(t)

            #print __name__
            if 1 or __name__ == '__main__':
                p = ctx.Process(target=qwork, args=(i,ress[l:h], dmap, allAtTree, log))
                p.start()
                threads.append(p)

        for i, t in enumerate(threads) :
            print "j %d" % (i)
            t.join()

    else :

        for ri, r in enumerate ( ress ) :

            r.bbZ = RadAts ( r.bbAtoms, dmap, allAtTree=allAtTree, show=0, log=0, numPts=10, toRAD=2, dRAD=0.2 )
            r.scZ = RadAts ( r.scAtoms, dmap, allAtTree=allAtTree, show=0, log=0, numPts=10, toRAD=2, dRAD=0.2 )

            if log and ri % 10 == 0 :
                status ( "Calculating - res %d/%d" % (ri, len(ress)) )
                print ".",



    scoresBB, scoresSC = [], []

    ress = []
    for r in mol.residues :
        if cid == None or r.id.chainId == cid :
            ress.append ( r )
            if r.bbZ != None : scoresBB.append ( r.bbZ )
            if r.scZ != None : scoresSC.append ( r.scZ )

    #sc = [x for x in scores if x is not None]
    #scSC = [1.0/x for x in scoresSC if x is not None]
    #scBB = [1.0/x for x in scoresBB if x is not None]

    #print " - %d res, SC min %.2f max %.2f, avg %.2f" % (len(ress), min(scSC), max(scSC), numpy.average(scSC))
    print " - avg sigma - side chain %.1f, backbone %.1f" % (numpy.average(scoresSC), numpy.average(scoresBB) )


    if 0 :

        sByType = {}
        rByType = {}
        for r in ress :
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


        mapName = os.path.splitext(dmap.name)[0]
        molName = os.path.splitext(mol.name)[0]
        mdir, mpfile = os.path.split(dmap.data.path)
        foname = mdir + "/" + mapName + "__" + molName + ".txt"


        print " - scores to: " + foname
        fp = open (foname,"w")

        for avgScore, rtype in avgs :

            rscores = rByType[rtype]
            rscores.sort ( reverse=False, key=lambda x: x[0] )
            hr = rscores[0]
            R = hr[1]
            highestScore = hr[0]
            numRes = len(rscores)

            rts = ""
            if R.isProt : rts = protein3to1[rtype]
            else : rts = nucleic3to1[rtype]

            print "%s\t%s\t%d\t%f\t%d\t.%s\t%f" % (rtype, rts, numRes, avgScore, R.id.position, R.id.chainId, highestScore)
            fp.write ( "%s\t%s\t%d\t%f\t%d\t.%s\t%f\n" % (rtype, rts, numRes, avgScore, R.id.position, R.id.chainId, highestScore) )

        fp.close()


    return numpy.average(scoresBB), numpy.average(scoresSC)


def CalcResQ (r, dmap, sigma, allAtTree=None, numPts=8, toRAD=2.0, dRAD=0.1, minD=0.0, maxD=1.0, useOld=False ) :

    scQ, bbQ, Q, numSC, numBB = 0.0, 0.0, 0.0, 0.0, 0.0
    for at in r.atoms :
        if at.element.name == "H" :
            continue

        if not hasattr ( at, 'Q' ) or not useOld :
            #cc, ccm = RadCC ( [at], dmap, sigma, allAtTree=allAtTree, show=0, log=0, numPts=numPts, toRAD=toRAD, dRAD=dRAD, minD=minD, maxD=maxD )
            #at.Q = ccm
            #at.CC = cc
            at.Q = 0
            at.CC = 0

        Q += at.Q
        if r.isProt or r.isNA :
            if at.isBB :
                bbQ += at.Q
                numBB += 1.0
            else :
                scQ += at.Q
                numSC += 1.0

    if r.isProt or r.isNA :
        if int(numSC) != len(r.scAtoms) :
            print " - res %d.%s.%s - %.0f/%d sc atoms" % (r.id.position,r.type,r.id.chainId, numSC, len(r.scAtoms))

        if numSC > 0 :
            r.scQ = scQ / numSC
        else :
            r.scQ = None

        if numBB > 0 :
            r.bbQ = bbQ / numBB
        else :
            r.bbQ = None

    r.Q = Q / float ( len(r.atoms) )



def CalcQ_ ( mol, cid, dmap, sigma=0.5, allAtTree=None, useOld=False, log=False ) :

    print "Q Scores - in parallel"
    print " - map: %s" % dmap.name
    print " - mol: %s, chain: %s" % (mol.name, cid if cid != None else "_all_")

    ress = []
    for r in mol.residues :
        if cid == None or r.id.chainId == cid :
            ress.append ( r )

    print " - residues to do: %d" % len(ress)


    import multiprocessing
    threads = multiprocessing.cpu_count() / 2
    print 'calc q using %d threads' % threads

    # Avoid periodic Python context switching.
    import sys
    original_check_interval = sys.getcheckinterval()
    sys.setcheckinterval(1000000000)

    # Define thread class for fitting.
    from threading import Thread
    class Q_Thread(Thread):
        def __init__(self, ress, ti):
            Thread.__init__(self)
            self.ress = ress
            self.ti = ti
        def run(self):
            print "run - %d - %d" % (self.ti, len(ress))
            for ri, r in enumerate ( self.ress ) :
                #CalcResQ (r, dmap, sigma, allAtTree=allAtTree, numPts=2, toRAD=2.0, dRAD=0.2 )
                #print "%d-%d/%d" % (ti,ri/len(self.ress)),
                for at in r.atoms :
                    if at.element.name != "H" :
                        cc, ccm = RadCC ( [at], dmap, sigma, allAtTree=allAtTree, show=0, log=0, numPts=8, toRAD=2.0, dRAD=0.5 )


    # Starts threads with each calculating an equal number of fits.
    n  = len(ress)
    g = [ress[(n*c)/threads:(n*(c+1))/threads] for c in range(threads)]
    threads = []
    for mi, ml in enumerate(g) :
        #print "%d - %d, %d-%d" % (mi,len(ml),ml[0].id.position,ml[-1].id.position)
        t = Q_Thread(ml,mi)
        threads.append(t)

    for t in threads:
        t.start()
    print ""

    # Wait for all threads to finish
    for t in threads:
        t.join()

    # Restore periodic context switching.
    sys.setcheckinterval(original_check_interval)

    # Collect fit results from all threads.
    #for t in threads:
    #    print "",





def CalcQ ( mol, cid, dmap, sigma=0.6, allAtTree=None, useOld=False, log=False ) :


    minD, maxD = MinMaxD ( dmap )

    print ""
    print "Q Scores"
    print " - map: %s" % dmap.name
    print " - mol: %s, chain: %s" % (mol.name, cid if cid != None else "_all_")
    print " - sigma: %.2f" % sigma
    print " - mind: %.3f, maxd: %.3f" % (minD, maxD)

    SetBBAts ( mol )


    atoms = []

    import time
    start = time.time()

    #ress = []
    for r in mol.residues :
        if cid == None or cid == "All" or r.id.chainId == cid :
            for at in r.atoms :
                if at.element.name == "H" :
                    continue
                atoms.append ( at )

    print " - atoms to do: %d" % len(atoms)

    #for ai, at in enumerate ( atoms[0:2] ) :
    #    cc, ccm = RadCC ( [at], dmap, sigma, allAtTree=allAtTree, show=0, log=0, numPts=8, toRAD=2.0, dRAD=0.1, minD=minD, maxD=maxD )


    from chimera import tasks, CancelOperation
    task = tasks.Task('Calculating Q-scores', modal = True)


    try :

        for ai, at in enumerate ( atoms ) :

            cc, ccm = RadCC ( [at], dmap, sigma, allAtTree=allAtTree, show=0, log=0, numPts=8, toRAD=2.0, dRAD=0.1, minD=minD, maxD=maxD )

            at.Q = ccm
            at.occupancy = ccm

            end = time.time()
            totSec = end - start

            leftTime = ""
            leftSec = 0.0
            iPerSec = float(ai) / totSec
            if iPerSec > 0 :
                leftSec = float ( len(atoms) - ai ) / iPerSec
                leftHour = numpy.floor ( leftSec / 60.0 / 60.0 )
                leftSec = leftSec - leftHour * 60.0 * 60.0
                leftMin = numpy.floor ( leftSec / 60.0 )
                leftSec = leftSec - leftMin * 60.0
                leftTime = "%.0f:%.0f:%.0f" % (leftHour, leftMin, leftSec)
                #hah


            if (ai+1) % 10 == 0 :
                if log :
                    status ( "Calculating Q scores - atom %d/%d - eta: %s" % (ai+1, len(atoms), leftTime) )
                    #print ".",

            if (ai+1) % 100 == 0 :
                print " - at atom %d/%d - eta: %s " % (ai+1, len(atoms), leftTime)

            #task.updateStatus( "Calculating Q scores - atom %d/%d - %s in %s.%d.%s - eta: %s" % (ai+1, len(atoms), at.name, at.residue.type, at.residue.id.position, at.residue.id.chainId, leftTime) )
            task.updateStatus( "Calculating Q scores - atom %d/%d - eta: %s" % (ai+1, len(atoms), leftTime) )

    #sc = [x for x in scores if x is not None]
    #scSC = [1.0/x for x in scoresSC if x is not None]
    #scBB = [1.0/x for x in scoresBB if x is not None]
    except :
        umsg ( "Something went wrong..." )
        return


    finally :
        task.finished()


    end = time.time()
    print ""
    print " - done, time: %f" % ( end-start )
    totSec = end - start
    totMin = numpy.floor ( totSec / 60.0 )
    totSec = totSec - totMin * 60.0
    print " - done, time: %.0f min, %.1f sec" % ( totMin, totSec )

    molPath = os.path.splitext(mol.openedAs[0])[0]
    mapName = os.path.splitext(dmap.name)[0]

    try :
        nname = molPath + "__Q__" + mapName + ".pdb"
        print "Saving pdb with Q-scores:", nname
        chimera.PDBio().writePDBfile ( [mol], nname )
        #umsg ( "Done Q-scores - saved %s with Q-scores in B-factor column" % nname )
    except :
        pass


    #umsg ( "Done Q-scores" )


    q, qcc = QStats1 ( mol, cid )

    print "Map: %s, Model: %s" % (dmap.name, mol.name)
    print " --> avg. Q score: %.3f" % q
    print ""

    return q, qcc







def CalcQForOpenModelsRess () :

    from VolumeViewer import Volume
    dmap = chimera.openModels.list(modelTypes = [Volume])[0]
    print " - dmap: %s" % dmap.name


    minD, maxD = MinMaxD ( dmap )
    print " - mind: %.3f, maxd: %.3f" % (minD, maxD)

    #fp = open ( "/Users/greg/_data/_mapsq/scores.txt", "a" )
    #fp.write ( "%s...\n" % dmap.name.split("_")[0]  )
    #fp.close ()

    from chimera import Molecule
    mol = chimera.openModels.list(modelTypes = [Molecule])[0]
    print " - mol: %s" % mol.name
    SetBBAts ( mol )


    #rids = {}
    #for r in mol.residues :
    #    rids["%d.%s" % (r.id.position,r.id.chainId)] = r

    atids = {}
    for r in mol.residues :
        for at in r.atoms :
            r = at.residue
            altLoc = '_' if at.altLoc == '' else at.altLoc
            atids["%d.%s.%s.%s" % (r.id.position,r.id.chainId,at.name,altLoc)] = at


    ats = [at for at in mol.atoms if not at.element.name == "H"]
    points = _multiscale.get_atom_coordinates ( ats, transformed = False )
    print " - search tree: %d/%d ats" % ( len(ats), len(mol.atoms) )
    allAtTree = AdaptiveTree ( points.tolist(), ats, 1.0)


    fin = open ( os.path.splitext ( dmap.data.path )[0] + ".txt" )
    fout = open ( os.path.splitext ( dmap.data.path )[0] + "_out.txt", "w" )
    foutn = os.path.splitext ( dmap.data.path )[0] + "_stat.txt"

    sig_at = []

    for l in fin :
        #print l,
        sigma, atIdStr = l.split()
        if not atIdStr in atids :
            print " - atid not found: ", atIdStr
        at = atids[atIdStr.strip()]
        sigma = float(sigma)
        sig_at.append ( [sigma, at, atIdStr] )

    fs = open ( foutn, "w" ); fs.write ( "%d/%d" % (0,len(sig_at) ) ); fs.close()

    import time
    start = time.time()

    i = 0
    for sigma, at, atId in sig_at :
        #print "%d.%s.%s" % (r.id.position,r.id.chainId,at.name),
        cc, ccm = RadCC ( [at], dmap, sigma, allAtTree=allAtTree, show=0, log=0, numPts=8, toRAD=2.0, dRAD=0.1, minD=minD, maxD=maxD )
        #print cc, ccm
        fout.write ( "%s %f %f\n" % (atId, cc, ccm) )

        if i%10 == 0 :

            end = time.time()
            totSec = end - start

            leftTime = ""
            leftSec = 0.0
            iPerSec = float(i) / totSec
            if iPerSec > 0 :
                leftSec = float ( len(sig_at) - i ) / iPerSec
                leftHour = numpy.floor ( leftSec / 60.0 / 60.0 )
                leftSec = leftSec - leftHour * 60.0 * 60.0
                leftMin = numpy.floor ( leftSec / 60.0 )
                leftSec = leftSec - leftMin * 60.0
                leftTime = "%.0f:%.0f:%.0f" % (leftHour, leftMin, leftSec)


            fs = open ( foutn, "w" ); fs.write ( "%d/%d - %s" % (i+1,len(sig_at),leftTime) ); fs.close()



        i += 1




    fin.close()
    fout.close()

    fs = open ( foutn, "w" ); fs.write ( "done" ); fs.close()



def QsFromFile ( mol, nname ) :

    rids = {}
    for r in mol.residues :
        rids["%d.%s" % (r.id.position,r.id.chainId)] = r


    # http://www.wwpdb.org/documentation/file-format-content/format33/sect9.html#ATOM
    try :
        fin = open ( nname, "r" )
    except :
        #print " - file not found"
        return False

    print " - Qs from file: %s" % nname

    for line in fin :
        if line[0:4] == "ATOM" or line[0:6] == "HETATM" :
            aname, aloc, cid, resi, occ, bfac = line[12:16].strip(), line[16:17].strip(), line[21], int(line[22:26]), float ( line[54:60] ), float ( line[60:66] )
            #if occ < 1.0 :
            rid = "%s.%s" % (resi,cid)
            if rid in rids :
                r = rids[rid]

                if aname in r.atomsMap :
                    ats = r.atomsMap[aname]
                    found = False
                    for at in ats :
                        if at.altLoc == aloc :
                            at.occupancy = at.Q = occ
                            found = True
                    if not found :
                        #print " -xx- %s.%s - atom %s - loc %s" % (resi, cid, aname, aloc)
                        continue
                else :
                    #print " -xx- %s.%s - atom %s" % (resi,cid, aname)
                    continue


    fin.close ()

    return True


def Calc ( chimeraPath, numProc, res=3.0 ) :

    print "Calc Q scores, v%s:" % mapqVersion
    print " - chimera path: ", chimeraPath
    print " - num processors: ", numProc
    print " - resolution: ", res


    from VolumeViewer import Volume
    vols = chimera.openModels.list(modelTypes = [Volume])
    if len(vols) == 0 :
        print " - no volumes loaded"
        return

    dmap = vols[0]
    print " - volume: %s" % dmap.name

    from chimera import Molecule

    mols = chimera.openModels.list(modelTypes = [Molecule])
    if len(mols) == 0 :
        print " - no molecules loaded"
        return

    for mi, mol in enumerate (mols) :

        print ""
        print "Model %d/%d: %s" % (mi+1, len(mols), mol.name)
        SetBBAts ( mol )

        ats = [at for at in mol.atoms if not at.element.name == "H"]
        points = _multiscale.get_atom_coordinates ( ats, transformed = False )
        print " - search tree: %d/%d ats" % ( len(ats), len(mol.atoms) )
        allAtTree = AdaptiveTree ( points.tolist(), ats, 1.0)
        #allAtTree = None

        if numProc == 1 :
            CalcQ ( mol, None, dmap, allAtTree=allAtTree )
        else :
            CalcQp ( mol, None, dmap, allAtTree=allAtTree, numProc=numProc )

        SaveQ ( mol, "All", dmap, res )



def CalcQp ( mol, cid, dmap, sigma=0.6, allAtTree=None, useOld=True, log=False, numProc=None ) :


    molPath = os.path.splitext(mol.openedAs[0])[0]
    mapName = os.path.splitext(dmap.name)[0]
    nname = molPath + "__Q__" + mapName + ".pdb"

    if 0 :
        if QsFromFile ( mol, nname ) :
            q, qcc = QStats1 ( mol, cid )
            return q, qcc


    import sys
    if numProc == None :
        try :
            numProc = int ( sys.argv[-1] )
        except :
            print " -- did not find # processors as last argument in the command"

            import multiprocessing
            numProc = multiprocessing.cpu_count() / 2
            print " - num processors detected:", numProc


    M = dmap.data.full_matrix()
    minD, maxD = numpy.min(M), numpy.max(M)

    print "Q Scores - p - %d" % numProc
    print " - map: %s" % dmap.name
    print " - mol: %s, chain: %s" % (mol.name, cid if cid != None else "_all_")
    print " - sigma: %.2f" % sigma
    print " - mind: %.3f, maxd: %.3f" % (minD, maxD)

    import time
    start = time.time()

    SetBBAts ( mol )

    ress = []
    atoms = []
    for r in mol.residues :
        if cid == None or cid == "All" or r.id.chainId == cid :
            if 1 or r.isNA :
                for at in r.atoms :
                    if 0 or not at.element.name == "H" :
                        atoms.append ( at )

    print " - atoms to do: %d" % len(atoms)

    import subprocess
    mapPath = os.path.split ( dmap.data.path )[0]
    mapBase = os.path.splitext ( dmap.data.path )[0]

    print "cmd:",
    #print sys.argv
    for arg in sys.argv :
        print arg,
    print ""

    # '/Users/greg/_mol/Chimera.app/Contents/Resources/share/__main__.py'
    #chiPath = os.path.split ( sys.argv[0] )[0]
    #mapQPPath = os.path.join ( chiPath, 'Segger' )
    #mapQPPath = os.path.join ( chiPath, 'mapqp.py' )
    #print " -- path to mapQ script:", mapQPPath


    # backtrack through path to chimera script to find executable
    print "Finding Chimera exec:"
    atPath = os.path.split ( sys.argv[0] )[0]
    chiPath = None
    for i in range (100) :
        atPath = os.path.split ( atPath )[0]
        print " --- at %s" % atPath
        for f1, f2 in ( ('bin', 'chimera'), ('MacOS', 'chimera') ) :
            tryPath = os.path.join(os.path.join(atPath,f1),f2)
            if os.path.isfile(tryPath) :
                print " - found path: %s" % tryPath
                chiPath = tryPath
                break

        if chiPath != None or len(atPath) <= 1 :
            break

    if chiPath == None :
        print " - did not find Chimera exec path"
        return

    #print " - path to Chimera:", chiPath

    dir_path = os.path.dirname(os.path.realpath(__file__))
    inDir = os.path.split(dir_path)[0]
    print "Working dir:", inDir
    mapQPPath = os.path.join ( inDir, 'mapq' )
    mapQPPath = os.path.join ( mapQPPath, 'mapqp.py' )
    print " -- path to mapQ script:", mapQPPath

    mapBase = mapBase + "_Q-score-mp"

    n = len(atoms)
    g = [atoms[(n*c)/numProc:(n*(c+1))/numProc] for c in range(numProc)]
    procs = []
    for mi, atoms1 in enumerate(g) :

        ress1 = atoms1[0].residue
        ressN = atoms1[-1].residue
        print " - %d/%d, %d-%d" % (mi+1, numProc, ress1.id.position, ressN.id.position)

        fout = open ( mapBase + "_%d.txt" % mi, "w" )
        for at in atoms1 :
            r = at.residue
            altLoc = '_' if at.altLoc == '' else at.altLoc
            fout.write ( "%.3f %d.%s.%s.%s\n" % (sigma, r.id.position,r.id.chainId,at.name,altLoc) )
        fout.close()

        nmap_path = mapBase + "_%d.mrc" % mi
        #print " -> ", nmap_path
        nmap = MaskMapResize ( atoms1, 4.0, dmap, nmap_path )
        #nmap.write_file ( nmap_path , "mrc" )

        args = [chiPath, '--nogui', '--silent', '--nostatus', mol.openedAs[0], nmap_path, mapQPPath]
        if mi == 0 :
            print "running proc:",
            for arg in args :
                print arg,
            print ""

        fout = open ( mapBase + "_%d.log" % mi, "w" )
        foute = open ( mapBase + "_%d_err.log" % mi, "w" )
        p = subprocess.Popen(args, stdout=fout, stderr=foute, cwd=inDir)
        procs.append ( [mi, p, fout, foute] )

    print ""
    print "Waiting...",
    for mi, p, fout, foute in procs :
        p.wait()
        fout.close()
        foute.close()
        print "%d" % mi,
    print ""

    atids = {}
    for r in mol.residues :
        for at in r.atoms :
            r = at.residue
            altLoc = '_' if at.altLoc == '' else at.altLoc
            atids["%d.%s.%s.%s" % (r.id.position,r.id.chainId,at.name,altLoc)] = at

    print ""
    print "Getting...",
    for mi, p, fout, foute in procs :
        fin = mapBase + "_%d_out.txt" % mi
        #print " - getting from: ", fin
        fp = open ( fin )
        for l in fp :
            #print " - ", l
            atId, cc, ccm = l.split()
            at = atids[atId.strip()]
            #at = r.atomsMap[atName][0]
            at.Q = float(ccm)
            at.CC = float(cc)
            at.occupancy = at.Q

        fp.close()

        if mi == 0 :
            print ""
            print ""
            print "__Out %d__" % mi
            foute = open ( mapBase + "_%d.log" % mi, "r" )
            for l in foute :
                print l,
            print ""
            foute.close()


            print "__Err %d__" % mi
            foute = open ( mapBase + "_%d_err.log" % mi, "r" )
            for l in foute :
                print l,
            print ""
            foute.close()

        #print " - removing..."
        os.remove ( mapBase + "_%d_out.txt" % mi )
        os.remove ( mapBase + "_%d_stat.txt" % mi )
        os.remove ( mapBase + "_%d.txt" % mi )
        os.remove ( mapBase + "_%d.mrc" % mi )
        os.remove ( mapBase + "_%d.log" % mi )
        os.remove ( mapBase + "_%d_err.log" % mi )
        print "%d" % mi,

    print ""


    end = time.time()
    print ""
    print " - done, time: %f" % ( end-start )
    totSec = end - start
    totMin = numpy.floor ( totSec / 60.0 )
    totSec = totSec - totMin * 60.0
    print " - done, time: %.0f min, %.1f sec" % ( totMin, totSec )


    molPath = os.path.splitext(mol.openedAs[0])[0]
    mapName = os.path.splitext(dmap.name)[0]

    nname = molPath + "__Q__" + mapName + ".pdb"
    print "Saving pdb with Q-scores:", nname
    chimera.PDBio().writePDBfile ( [mol], nname )

    q, qcc = QStats1 ( mol, cid )

    print "Map: %s, Model: %s" % (dmap.name, mol.name)
    print " --> avg. Q score: %.3f" % q
    print ""

    return q, qcc



def QStats1 ( mol, chainId ) :

    totQ, totCC, totN = 0.0, 0.0, 0.0
    #QT, QN = { "Protein":0.0, "Nucleic":0.0, "Other":0.0 }, { "Protein":0.0, "Nucleic":0.0, "Other":0.0}
    QT, QN = {}, {}
    QT_, QN_ = {}, {}

    doRess = []

    for r in mol.residues :
        #if r.id.chainId == chainId or chainId == None :
        doRess.append ( r )

    print ""
    print "Q for %d res..." % ( len(doRess) )
    for r in doRess :


        #if not r.isNA : continue
        #if not r.isProt : continue

        CalcResQ (r, None, None, useOld=True )

        for at in r.atoms :
            if at.element.name == "H" :
                continue
            if hasattr ( at, "Q") :
                #if 0 or at.residue.isProt :
                if 1 or r.isNA :
                    totQ += at.Q
                    totN += 1.0
                    if hasattr(at, "CC") :
                        totCC += at.CC

                tp = "Other"
                if at.residue.isProt : tp = "Protein"
                elif at.residue.isNA : tp = "Nucleic"
                else : tp = at.residue.type

                if tp in QT :
                    QT[tp] += at.Q; QN[tp] += 1.0
                else :
                    QT[tp] = at.Q; QN[tp] = 1.0

                tps = "chain " + r.id.chainId + " (" + tp.lower() + ")"
                if tps in QT_ :
                    QT_[tps] += at.Q; QN_[tps] += 1.0
                else :
                    QT_[tps] = at.Q; QN_[tps] = 1.0

    try :
        umsg ( "Model Q-score: %.3f" % (totQ/totN) )
    except :
        pass

    #for tp in ["Other", "Protein", "Nucleic"] :
    print ""
    print "Chain\tAvg.Q-score\tEst.Res.(A)"
    tpk = QT_.keys()
    tpk.sort()
    for tp in tpk :
        if QN_[tp] > 0 :
            avgQ = QT_[tp]/QN_[tp]
            avgR = 0
            if "nucleic" in tp.lower() :
                avgR = (avgQ-1.0673)/-0.1574
            else :
                avgR = (avgQ-1.1244)/-0.1794
            print " %s\t%.3f\t%.2f" % (tp, avgQ, avgR )
        else :
            print " %s\tn/a" % (tp)

    #for tp in ["Other", "Protein", "Nucleic"] :
    print ""
    print "Type\tAvg.Q-score\tEst.Res.(A)"
    for tp in QT.keys() :
        if QN[tp] > 0 :
            avgQ = QT[tp]/QN[tp]
            avgR = 0
            if "nucleic" in tp.lower() :
                avgR = (avgQ-1.0673)/-0.1574
            else :
                avgR = (avgQ-1.1244)/-0.1794
            print " %s\t%.3f\t%.2f" % (tp, avgQ, avgR )
        else :
            print " %s\tn/a" % (tp)

    print ""

    return totQ/totN, totCC/totN




def SaveQ ( mol, chainId, dmap, RES=3.0 ) :

    avgQrna = -0.1574 * RES + 1.0673 # rna
    avgQprot = -0.1794 * RES + 1.1244 # protein

    cres = {}
    for r in mol.residues :
        try :
            cres[r.id.chainId].append ( [r.id.position, r] )
        except :
            cres[r.id.chainId] = [ [r.id.position, r] ]


    molPath = os.path.splitext(mol.openedAs[0])[0]
    mapName = os.path.splitext(dmap.name)[0]
    nname = molPath + "__Q__" + mapName + "_" + chainId + ".txt"
    #nname = molPath + "__Q__" + mapName + "_" + cid + ".txt"
    fp = open (nname, "w")

    print "Saving per-chain & per-residue Q-scores:"
    print " -> res=", RES
    print " -> file:", nname

    fp.write ( "Chain\tQ_chain\tEst.Res.\tExpectedQ@%.2f\n" % RES )

    chains = cres.keys()
    chains.sort()

    for cid in chains :
        if 0 or cid == chainId or chainId == "All" :

            tps = {}
            resAtoms = []
            rs = cres[cid]
            for ri, r in rs :
                resAtoms.extend ( r.atoms )
                tp = "Other"
                if r.isProt : tp = "Protein"
                if r.isNA : tp = "Nucleic"
                tps[tp] = 1

            ctypes = ""
            for tp in tps.keys() :
                ctypes = (ctypes + tp) if len(ctypes) == 0 else (ctypes + "," + tp)

            cQ = numpy.average ( [at.Q for at in resAtoms if at.element.name != "H"] )

            formula = "=-0.1775 * %.2f + 1.1192" % RES
            estRes = (cQ - 1.1192) / -0.1775
            if "Nucleic" in ctypes :
                formula = "= -0.1377 * %.2f + 0.9973" % RES
                estRes = (cQ - 0.9973) / -0.1377

            fp.write ( "%s\t%.2f\t%.2f\t%s\t(%s)\n" % (cid, cQ, estRes, formula, ctypes) )

            #print " - cid: %s - %s - %.2f" % (cid, ctypes, cQ)


    fp.write ( "\n" )
    fp.write ( "Protein: avgQ = -0.1775 * RES + 1.1192\n" )
    fp.write ( "Nucleic: avgQ = -0.1377 * RES + 0.9973\n" )
    fp.write ( "\n" )

    for cid in cres.keys () :

        if cid == chainId or chainId == "All" :

            fp.write ( "Chain %s\t\t\t\t\t\t\t\tAverage over 1 residue\t\t\t\t\tAverage over 2 residues\t\t\t\t\tAverage over 3 residues\t\t\t\t\tAverage over 5 residues\n\n" % cid )

            fp.write ( "Chain\tRes\tRes #\tQ_backBone\tQ_sideChain\tQ_residue\tExpectedQ@%.2f\t\t" % RES )
            fp.write ( "Q_backBone(avg-1)\tQ_sideChain(avg-1)\tQ_residue(avg-1)\tExpectedQ@%.2f\t\t" % RES )
            fp.write ( "Q_backBone(avg-2)\tQ_sideChain(avg-2)\tQ_residue(avg-2)\tExpectedQ@%.2f\t\t" % RES )
            fp.write ( "Q_backBone(avg-3)\tQ_sideChain(avg-3)\tQ_residue(avg-3)\tExpectedQ@%.2f\t\t" % RES )
            fp.write ( "Q_backBone(avg-5)\tQ_sideChain(avg-5)\tQ_residue(avg-5)\tExpectedQ@%.2f\t\n" % RES )

            #cid = 'A'
            rs = cres[cid]

            #print " - cid: %s - " % (cid)

            rs.sort()
            #for i in range (10) :
            #    print rs[i]

            ress = []
            Qs, AV, CC = [], [], []
            for ri, r in rs :

                #if not r.isProt and not r.isNA :
                #    print " - cid: %s - r %d - not prot or RNA" % (cid, r.id.position)
                #    continue

                ress.append (r)

                r.Q = numpy.average ( [at.Q for at in r.atoms if at.element.name != "H"] )

                r.qBB, r.qSC = 0, 0
                if len(r.bbAtoms) > 0 :
                    r.qBB = numpy.average ( [at.Q for at in r.bbAtoms if at.element.name != "H"] )
                if len(r.scAtoms) > 0 :
                    r.qSC = numpy.average ( [at.Q for at in r.scAtoms if at.element.name != "H"] )
                Qs.append ( [r.qBB, r.qSC, r.Q] )

                if 0 :
                    ad = avgdAts ( r.atoms, dmap )
                    aSC, aBB = 0, 0
                    if len(r.scAtoms) > 0 :
                        aSC = avgdAts ( r.scAtoms, dmap )
                    if len(r.bbAtoms) > 0 :
                        aBB = avgdAts ( r.bbAtoms, dmap )
                    AV.append ( [ad, aBB, aSC] )

                if 0 :
                    cc, ccm = ccAts ( r.atoms, dmap, RES )
                    ccSC, ccmSC = ccAts ( r.scAtoms, dmap, RES )
                    ccBB, ccmBB = ccAts ( r.bbAtoms, dmap, RES )
                    CC.append ( [cc, ccBB, ccSC] )
                    #CC.append ( [ccm, ccmBB, ccmSugar, ccmBase] )



            def N ( A, i, ind, N ) :
                #for i, a in enumerate ( A ) :
                sum, n = 0, 0
                for j in range ( i-N, i+N+1 ) :
                    if j >= 0 and j < len(A) :
                        sum += A[j][ind]
                        n += 1.0
                return sum/n


            last_i = None
            for i, r in enumerate ( ress ) :

                if (r.isNA or r.isProt) and last_i != None :
                    ii = last_i+1
                    while ii < r.id.position :
                        avgQ = avgQrna if r.isNA else avgQprot
                        fp.write ( "%s\t%s\t%d\t\t\t\t%f\t\t" % (r.id.chainId, "", ii, avgQ ) )
                        fp.write ( "\t\t\t%f\t\t" % (avgQ) )
                        fp.write ( "\t\t\t%f\t\t" % (avgQ) )
                        fp.write ( "\t\t\t%f\t\t" % (avgQ) )
                        fp.write ( "\t\t\t%f\n" % (avgQ) )
                        ii += 1

                if r.isNA :
                    fp.write ( "%s\t%s\t%d\t%f\t%f\t%f\t%f\t\t" % (r.id.chainId, r.type, r.id.position, r.qBB, r.qSC, r.Q, avgQrna ) )
                    fp.write ( "%f\t%f\t%f\t%f\t\t" % (N(Qs,i,0,1), N(Qs,i,1,1), N(Qs,i,2,1), avgQrna ) )
                    fp.write ( "%f\t%f\t%f\t%f\t\t" % (N(Qs,i,0,2), N(Qs,i,1,2), N(Qs,i,2,2), avgQrna ) )
                    fp.write ( "%f\t%f\t%f\t%f\t\t" % (N(Qs,i,0,3), N(Qs,i,1,3), N(Qs,i,2,3), avgQrna ) )
                    fp.write ( "%f\t%f\t%f\t%f\n" % (N(Qs,i,0,5), N(Qs,i,1,5), N(Qs,i,2,5), avgQrna ) )
                elif r.isProt :
                    fp.write ( "%s\t%s\t%d\t%f\t%f\t%f\t%f\t\t" % (r.id.chainId, r.type, r.id.position, r.qBB, r.qSC, r.Q, avgQprot ) )
                    fp.write ( "%f\t%f\t%f\t%f\t\t" % (N(Qs,i,0,1), N(Qs,i,1,1), N(Qs,i,2,1), avgQprot ) )
                    fp.write ( "%f\t%f\t%f\t%f\t\t" % (N(Qs,i,0,2), N(Qs,i,1,2), N(Qs,i,2,2), avgQprot ) )
                    fp.write ( "%f\t%f\t%f\t%f\t\t" % (N(Qs,i,0,3), N(Qs,i,1,3), N(Qs,i,2,3), avgQprot ) )
                    fp.write ( "%f\t%f\t%f\t%f\n" % (N(Qs,i,0,5), N(Qs,i,1,5), N(Qs,i,2,5), avgQprot ) )
                else :
                    fp.write ( "%s\t%s\t%d\t\t\t%f\t%f\t\t" % (r.id.chainId, r.type, r.id.position, r.Q, avgQprot ) )
                    fp.write ( "\t\t%f\t%f\t\t" % (N(Qs,i,2,1), avgQprot ) )
                    fp.write ( "\t\t%f\t%f\t\t" % (N(Qs,i,2,2), avgQprot ) )
                    fp.write ( "\t\t%f\t%f\t\t" % (N(Qs,i,2,3), avgQprot ) )
                    fp.write ( "\t\t%f\t%f\n" % (N(Qs,i,2,5), avgQprot ) )

                last_i = r.id.position


            fp.write ( "\n\n" )


    fp.close()
    print ""




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


def MinMaxD ( dmap ) :
    M = dmap.data.full_matrix()
    maxM = numpy.max(M)
    minM = numpy.min(M)

    maxD = min ( numpy.average(M)+numpy.std(M)*10, maxM )
    minD = max ( numpy.average(M)-numpy.std(M)*1, minM )

    # xray
    #maxD = min ( numpy.average(M)+numpy.std(M)*3.5, maxM )
    #minD = max ( numpy.average(M)-numpy.std(M)*0.77, minM )

    #print "%s - %.2f->%.2f, %.2f->%.2f" % (dmap.name, minD, maxD, minM, maxM )
    #minD = numpy.min(M)
    #minD, maxD = numpy.min(M), numpy.max(M)
    return minD, maxD


def RadCC ( atoms, dmap, sigma, allAtTree = None, show=0, log=0, numPts=8, toRAD=2.0, dRAD=0.5, minD=None, maxD=None, fitg=0, mol=None ) :

    if minD == None or maxD == None :
        minD, maxD = MinMaxD (dmap)

    #sigma = 1.0

    if len(atoms) == 0 :
        #print " - no RAD atoms?"
        return None

    from _multiscale import get_atom_coordinates
    pts = get_atom_coordinates(atoms, transformed = False)
    #print " __%s__ " % (atoms[0].name), pts[0]


    A, B = maxD - minD, minD
    refG = A * numpy.exp ( -0.5 * numpy.power(0.0/sigma,2) ) + B
    #print " - refg: ", refG

    # g_vals should have the reference gaussian...
    g_vals = (numpy.ones ( [len(pts)*numPts,1] ) * refG).astype(numpy.float64, copy=False)

    if mol == None :
        mol = atoms[0].molecule


    # r_avg holds the average values and number of points at each radial distance
    d_vals = dmap.interpolated_values ( pts, mol.openState.xform ).astype(numpy.float64, copy=False)

    d_vals = numpy.repeat ( d_vals, numPts )

    avgV = numpy.average ( d_vals )
    r_avg = [ [0,avgV,len(pts)*numPts] ]



    # make smaller atom tree...
    if 1 and allAtTree != None :
        ats_near = []
        for at in atoms :
            anear = allAtTree.searchTree ( at.coord().data(), toRAD*2.0 )
            ats_near.extend ( anear )

        points = _multiscale.get_atom_coordinates ( ats_near, transformed = False )
        if log :
            print " - new search tree: %d pts" % ( len(ats_near) )
        allAtTree = AdaptiveTree ( points.tolist(), ats_near, 1.0)



    #olap, corr1, corr2 = FitMap.overlap_and_correlation ( fpoint_weights, map_values )
    #dRAD, toRAD, RAD = 0.2, 1.8, 0.1
    RAD = dRAD
    i = 1.0
    while RAD < toRAD + 0.01 :
        outRad = RAD*0.9
        outRad2 = outRad * outRad
        #outRad2 = outRad * outRad
        pts = []
        for at in atoms :
            #npts = (numPts * RAD*RAD / (dRAD*dRAD)) if show else numPts
            #npts = numPts * (RAD*RAD / (dRAD*dRAD))
            npts = numPts # 8 # int ( npts )
            #print RAD, dRAD, numPts, " -> ", npts
            for i in range (0, 100) :
                outPts = SpherePts ( at.coord(), RAD, npts+i*2 )
                at_pts, at_pts_i = [None]*len(outPts), 0
                for pt in outPts :
                    vPt = [pt[0], pt[1], pt[2]]
                    apt = numpy.array ( vPt )
                    if allAtTree != None :
                        opointsNear = allAtTree.searchTree ( vPt, outRad )

                        if 1 :
                            foundNearPt = False
                            for npt in opointsNear :
                                v = apt - npt.coord().data()
                                r2 = numpy.sum ( v * v )
                                if r2 < outRad2 :
                                    foundNearPt = True
                                    break
                            if not foundNearPt :
                                at_pts[at_pts_i] = vPt
                                at_pts_i += 1

                        else :
                            if len(opointsNear) == 0 :
                                at_pts[at_pts_i] = vPt
                                at_pts_i += 1
                    else :
                        at_pts[at_pts_i] = vPt
                        at_pts_i += 1
                #if log :
                #    print " - %d, %d pts" % (i, len(at_pts))
                if at_pts_i >= npts or i >= 95 :
                    pts.extend ( at_pts[0:at_pts_i] )
                    break

        if show :
            AddSpherePts ( pts, (.6,.6,.6,0.4), 0.1, "RAD points %.1f" % RAD )

        if len (pts) < 1 :
            if log :
                print " - no points for RAD %.1f - %d.%s - " % (RAD, atoms[0].residue.id.position, atoms[0].residue.type),
                print "SC" if atoms[0].isSC else "BB"

            r_avg.append ( [RAD,0,0] )


        else :
            d_vals_n = dmap.interpolated_values ( pts, mol.openState.xform )
            d_vals = numpy.append ( d_vals, d_vals_n )
            avg = numpy.average ( d_vals_n )

            #gv = A * numpy.exp ( -0.5 * numpy.power(x/sdev,2) ) + B
            #A, B = GV, 0
            #A, B = GV - minD, minD
            A,B = maxD - minD, minD
            gv = A * numpy.exp ( -0.5 * numpy.power(RAD/sigma,2) ) + B
            g_vals = numpy.append ( g_vals, numpy.ones([len(pts),1]) * gv )

            r_avg.append ( [RAD,avg,len(pts)] )

            #if log :
            #    print "%.1f\t%f\t%f\t%d" % (RAD, avg, gv, len(pts))

        RAD += dRAD
        i+=1

    if log :
        min, max = r_avg[0][1], r_avg[0][1]
        for RAD, avg, numPts in r_avg :
            if avg < min : min = avg
            if avg > max : max = avg
        A,B = max-min, min
        A,B = maxD - minD, minD
        #A,B = GV - minD, minD
        for RAD, avg, numPts in r_avg :
            gv = A * numpy.exp ( -0.5 * numpy.power(RAD/sigma,2) ) + B
            #print "%.1f\t%f\t%f\t%d" % (RAD, avg+0.02, gv+0.02, numPts)
            print "%.1f\t%f\t%f\t%d" % (RAD, avg, gv, numPts)

    #d_vals = d_vals + 0.02
    #g_vals = g_vals + 0.02

    olap, CC, CCm = FitMap.overlap_and_correlation ( d_vals, g_vals )
    if log :
        print "olap: %.3f cc: %.3f, ccm: %.3f" % (olap, CC, CCm)
        print "%f\t%f\t%f" % (olap, CC, CCm)

    if fitg :
        if log : print "fitting gaussian : "
        #V, N = [ [x[0],x[1]] for x in r_avg ], float(len(r_avg))
        V, N = [ [x[0],x[1]] for x in r_avg[0:15] ], float(15)

        sdev, A, B = optSGD ( V, 5000, 1.0 )
        sdev, A, B = optSGD ( V, 5000, 0.1, sdev, A, B )
        err = numpy.sqrt(err3(V,sdev,A,B)/N)
        if log : print " sgd - sdev: %.4f, A %.4f, B %.4f, err: %f" % (sdev, A, B, err)
        sdev2, A2, B2 = optGN ( V, 0.0001, sdev, A, B )
        if sdev2 != None :
            sdev, A, B = sdev2, A2, B2
            err = numpy.sqrt(err3(V,sdev,A,B)/N)
            print "max:", r_avg[0][1]
            errp = err / r_avg[0][1] * 100.0
            if log : print "  gn - sdev: %.4f, A %.4f, B %.4f, err: %f (%.1f%%)" % (sdev, A, B, err, errp)

        yds, i = numpy.zeros ( len(r_avg) ), 0
        for x, y, n in r_avg:
            gv = A * numpy.exp ( -0.5 * numpy.power(x/sdev,2) ) + B
            #yds[i] = y - gv
            yds[i] = y
            if log : print "%.1f\t%f\t%f" % (x, y, gv)
            i += 1

        return CC, CCm, yds, err

    else :
        return CC, CCm



def RadAts ( atoms, dmap, allAtTree = None, show=0, log=0, numPts=20, toRAD=2.0, dRAD=0.1 ) :

    if len(atoms) == 0 :
        #print " - no RAD atoms?"
        return None

    #pts = []
    #for at in atoms :
    #    p = at.coord()
    #    pts.append ( [p[0], p[1], p[2]] )

    from _multiscale import get_atom_coordinates
    pts = get_atom_coordinates(atoms, transformed = False)

    RD_, X, Y = [], [], []
    d_vals = dmap.interpolated_values ( pts, atoms[0].molecule.openState.xform )
    avg = numpy.average ( d_vals )

    RD_.append ( [0,avg] ); X.append (0); Y.append (avg)


    #dRAD, toRAD, RAD = 0.2, 1.8, 0.1
    RAD = dRAD
    i = 1.0
    while RAD < toRAD + 0.01 :
        outRad = RAD*0.9
        outRad2 = outRad * outRad
        pts = []
        for at in atoms :
            npts = (numPts * RAD*RAD / (dRAD*dRAD)) if show else numPts
            npts = int ( npts )
            #print RAD, dRAD, numPts, " -> ", npts
            outPts = SpherePts ( at.coord(), RAD, npts )
            for pt in outPts :
                ppt = [pt[0], pt[1], pt[2]]
                if allAtTree != None :
                    vPt = numpy.array ( ppt )
                    opointsNear = allAtTree.searchTree ( ppt, outRad )
                    if 1 :
                        clash = False
                        for p in opointsNear :
                            v = vPt - p.coord().data()
                            sqSum = numpy.sum ( v * v )
                            if sqSum < outRad2 :
                                clash = True
                                break
                        if clash == False :
                            pts.append ( ppt )

                    else :
                        if len(opointsNear) == 0 :
                            pts.append ( ppt )
                else :
                    pts.append ( ppt )

        if show :
            AddSpherePts ( pts, (.6,.6,.6,0.4), 0.1, "RAD points %.1f" % RAD )

        if len (pts) < 1 :
            if log :
                print " - no points for RAD %.1f - %d.%s - " % (RAD, atoms[0].residue.id.position, atoms[0].residue.type),
                print "SC" if atoms[0].isSC else "BB"

        else :
            d_vals = dmap.interpolated_values ( pts, atoms[0].molecule.openState.xform )
            avg = numpy.average ( d_vals )
            RD_.append ( [RAD,avg] );
            if log :
                print RAD, avg, len(pts)
                X.append (RAD); Y.append (avg)

        RAD += dRAD

    #minSd = opt0 ( RD_, 0.1 )
    #if minSd != None :
    #    if show :
    #        print " SD0: %.1f" % minSd

    sdev = toRAD
    slope = 0

    if RD_[0][1] <=  RD_[-1][1] :
        sdev = 10.0

    else :

        #for i in range ( len(RD_) ) :
        #    RD_[i][1] = RD_[i][1] - RD_[-1][1]
        #    if log :
        #        Y[i] = Y[i] - Y[-1]


        #import time
        #start = time.time()
        sdev, A, B = optSGD ( RD_, 9000, 0.2 )
        sdev, A, B = optSGD ( RD_, 9000, 0.02, sdev, A, B )
        sdev, A, B = optSGD ( RD_, 9000, 0.002, sdev, A, B )
        #end = time.time()
        #if log : print " sgd - sdev: %.4f, A %.4f, B %.4f -- %f" % (sdev, A, B, (end - start))
        sdev = sdev
        if log : print " sgd - sdev: %.4f, A %.4f, B %.4f" % (sdev, A, B)

        #start = time.time()
        #sdev, A, B = optGN ( RD_, 0.0001 )
        #print " gn - sdev: %.4f, A %.4f, B %.4f -- %f" % (sdev, A, B, (end - start))
        #end = time.time()

        if 1 :
            if 0 and sdev != None :

                if log :
                    print " gn1 - sdev: %.4f, A %.4f, B %.4f" % (sdev, A, B)

            else :
                sdev, A, B = optSGD ( RD_, 10000, 0.01 )

                if log :
                    print " sgd - sdev: %.4f, A %.4f, B %.4f" % (sdev, A, B)

                sdev2, A2, B2 = optGN ( RD_, 0.0001, sdev, A, B )
                if sdev2 != None :
                    sdev, A, B = sdev2, A2, B2
                    if log :
                        print " gn2 - sdev: %.4f, A %.4f, B %.4f" % (sdev, A, B)
                #else :
                #    return 10.0


        if log :
            r = numpy.polyfit ( X, Y, 1, rcond=None, full=False, w=None, cov=False)
            print " sdev: %.4f, A %.4f, B %.4f // slope: %.4f y %.4f" % (sdev, A, B, r[0], r[1])

            #A, B = 0.26+0.08, -0.08
            lastX = 0
            for i in range ( len(RD_) ) :
                x, y = RD_[i]
                gv = A * numpy.exp ( -0.5 * numpy.power(x/sdev,2) ) + B
                gvRef = A * numpy.exp ( -0.5 * numpy.power(x/0.5,2) ) + B
                lv = x * r[0] + r[1]
                print "%.1f\t%f\t%f\t%f" % (x, y, gv, gvRef)
                lastX = x

            if 1 :
                x = lastX + dRAD
                #while x < min(4 * sdev,50.0) :
                while x < min(10.0,50.0) :
                    gv = A * numpy.exp ( -0.5 * numpy.power(x/sdev,2) ) + B
                    gvRef = A * numpy.exp ( -0.5 * numpy.power(x/0.5,2) ) + B
                    lv = x * r[0] + r[1]
                    print "%.1f\t\t%f\t%f" % (x, gv, gvRef)
                    x += dRAD


    #return abs(sdev), abs(slope)
    return abs(sdev)


def optGN ( V, err, S=None, A=None, B=None ) :

    y0 = V[0][1]
    yN = V[-1][1]

    if S == None :
        S = 0.5
        A = y0+yN
        B = yN

    an = numpy.array ( [A,B,S] )
    #print " _ -- A %.3f B %.3f s %.3f" % (A, B, S)

    reg = 1.0
    badMatCount = 0

    for i in range ( 1000 ) :

        J = numpy.zeros ( [len(V),3] )
        e = numpy.zeros ( [len(V),1] )

        err0 = 0
        j = 0
        for x,y in V :
            expv = numpy.exp ( -0.5 * numpy.power(x/S,2) )
            v = A * expv + B
            yd = v - y
            err0 += yd * yd
            #print "%.2f,%.2f/%.2f(%.2f)" % (x, y, v, yd),

            dA = expv
            dB = 1
            dS = A*x*x*numpy.power(S,-3) * expv
            J[j,:] = [dA, dB, dS]
            e[j,0] = yd
            j += 1

        Jt = numpy.transpose(J)

        try :
            J_ = numpy.dot ( numpy.linalg.inv ( numpy.dot(Jt,J) ), Jt )
        except :
            #print " - bad matrix?"
            #print numpy.dot(Jt,J)
            badMatCount += 1

            if badMatCount > 3 :
                return None, None, None

            from numpy import random as R
            an = numpy.array ( [R.random()*(y0+yN),R.random()*yN,R.random()*10.0] )
            A,B,S = an[0], an[1], an[2]
            #print " ? -- A %.3f B %.3f s %.3f" % (A, B, S)
            reg = 1.0

            continue

        ad = numpy.dot ( J_, e )
        ann = an - ( ad[:,0] * reg )
        A,B,S = ann[0], ann[1], ann[2]

        err1 = err3 ( V, S, A, B )
        #if err1 > err0 :
        #    reg = reg * 0.1
        #    if reg < err :
        #        break
        #else :
        an = ann
        #print " %d -- A %.3f B %.3f s %.3f - err %.3f, reg %.5f" % (i, A, B, S, err1, reg)

        if abs(err0 - err1) < err :
            #print " - done"
            break

        i += 1

    return S,A,B



def optSGD ( V, N, err, S=None, A=None, B=None ) :

    if S == None :
        y0 = V[0][1]
        yN = V[-1][1]
        S = 0.5
        A = y0+yN
        B = yN

    from numpy import random

    lastE = err3 ( V, S, A, B )
    #while True :
    for i in range(N) :

        S_ = S + random.normal ( 0, err ) # mean, sigma
        A_ = A + random.normal ( 0, err ) # mean, sigma
        B_ = B + random.normal ( 0, err ) # mean, sigma

        e = err3 ( V, S_, A_, B_ )
        #print "%d %.2f %f %f %.4f" % (i, sdAt, e, numpy.log(e), dd)
        if e < lastE :
            S, A, B = S_, A_, B_
            lastE = e

    return S,A,B


def err3 ( XYz, sd, A, B ) :

    y0 = XYz[0][1]
    err = 0
    #for x,y in XYz[1:] :
    for x,y in XYz :
        yd = y - A * numpy.exp ( -0.5 * numpy.power(x/sd,2) ) - B
        err += yd * yd
    #err /= float(len(XYz))
    return err



def err ( XYz, sd ) :

    y0 = XYz[0][1]
    err = 0
    for x,y in XYz[1:] :
        yd = y - y0 * numpy.exp ( -0.5 * numpy.power(x/sd,2) )
        err += yd * yd
    #err /= float(len(XYz))
    return err


def opt0 ( RD_, dStep ) :

    sd = 0.1
    y0 = RD_[0][1]
    minSd, minErr, N = None, 1e99, float ( len(RD_)-1 )
    while sd < 10.0 :

        err = 0
        for x,y in RD_[1:] :
            yd = y - y0 * numpy.exp ( -0.5 * numpy.power(x/sd,2) )
            err += yd * yd
        err /= N

        #print err

        if err < minErr :
            minErr = err
            minSd = sd

        sd += dStep


def opt ( V, maxErr ) :

    dd = 1.0
    sdAt = 0.1
    lastE = err ( V, sdAt )
    #while True :
    for i in range(10000) :
        sdAt += dd
        e = err ( V, sdAt )
        #print "%d %.2f %f %f %.4f" % (i, sdAt, e, numpy.log(e), dd)
        if e >= lastE :
            dd *= -0.75
            if abs(dd) < maxErr :
                return sdAt
        lastE = e
    return sdAt






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

    dmap = MapFromData ( mdata, nname, dmap, False )
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





def MaskMapResize ( atoms, R, dmap, fout=None ) :


    import _multiscale
    import _contour
    import _volume
    from _contour import affine_transform_vertices as transform_vertices
    from VolumeData import grid_indices, zone_masked_grid_data, interpolate_volume_data

    points = _multiscale.get_atom_coordinates ( atoms, transformed = True )
    #print " %d points" % len(points)
    fpoints = points


    if 0 :
        _contour.affine_transform_vertices ( points, Matrix.xform_matrix( dmap.openState.xform.inverse() ) )
        mdata = VolumeData.zone_masked_grid_data ( dmap.data, points, R )

        #mdata = VolumeData.Array_Grid_Data ( mdata.full_matrix(), dmap.data.origin, dmap.data.step, dmap.data.cell_angles, name = "atom masked" )

        mat = mdata.full_matrix()
        threshold = 1e-3

        points = _volume.high_indices(mat, threshold)
        fpoints = points.astype(numpy.single)
        fpoint_weights = mat[points[:,2],points[:,1],points[:,0]]

        #print " %d points" % len(points)

        nz = numpy.nonzero( fpoint_weights )[0]
        #print " %d pts nonzero" % len(nz)
        if len(nz) > 0 and len(nz) < len (fpoint_weights) :
            fpoints = numpy.take( fpoints, nz, axis=0 )
            fpoint_weights = numpy.take(fpoint_weights, nz, axis=0)

    else :
        _contour.affine_transform_vertices ( fpoints, Matrix.xform_matrix( dmap.openState.xform.inverse() ) )
        #transform_vertices( fpoints, dmap.data.ijk_to_xyz_transform )
        transform_vertices( fpoints, dmap.data.xyz_to_ijk_transform )

    #print " - %s mask %d atoms, %d nonzero points" % ( dmap.name, len(atoms), len(nz) )

    #transform_vertices( fpoints,  Matrix.xform_matrix( fmap.openState.xform ) )
    #transform_vertices( fpoints,  Matrix.xform_matrix( dmap.openState.xform.inverse() ) )
    #transform_vertices ( fpoints, dmap.data.xyz_to_ijk_transform )

    bound = 6
    li,lj,lk = numpy.min ( fpoints, axis=0 ) - (bound, bound, bound)
    hi,hj,hk = numpy.max ( fpoints, axis=0 ) + (bound, bound, bound)

    n1 = hi - li + 1
    n2 = hj - lj + 1
    n3 = hk - lk + 1

    #print " - bounds - %d %d %d --> %d %d %d --> %d %d %d" % ( li, lj, lk, hi, hj, hk, n1,n2,n3 )

    #nmat = numpy.zeros ( (n1,n2,n3), numpy.float32 )
    #dmat = dmap.full_matrix()

    nstep = (dmap.data.step[0], dmap.data.step[1], dmap.data.step[2] )

    nn1 = int ( round (dmap.data.step[0] * float(n1) / nstep[0]) )
    nn2 = int ( round (dmap.data.step[1] * float(n2) / nstep[1]) )
    nn3 = int ( round (dmap.data.step[2] * float(n3) / nstep[2]) )

    O = dmap.data.origin
    #print " - %s origin:" % dmap.name, O
    nO = ( O[0] + float(li) * dmap.data.step[0],
           O[1] + float(lj) * dmap.data.step[1],
           O[2] + float(lk) * dmap.data.step[2] )

    #print " - new map origin:", nO

    ox = round ( nO[0]/dmap.data.step[0] ) * dmap.data.step[0]
    oy = round ( nO[1]/dmap.data.step[1] ) * dmap.data.step[1]
    oz = round ( nO[2]/dmap.data.step[2] ) * dmap.data.step[2]

    nO = ( ox, oy, oz )

    #print " - new map origin:", nO


    nmat = numpy.zeros ( (nn1,nn2,nn3), numpy.float32 )
    ndata = VolumeData.Array_Grid_Data ( nmat, nO, nstep, dmap.data.cell_angles )

    npoints = grid_indices ( (nn1, nn2, nn3), numpy.single)  # i,j,k indices
    transform_vertices ( npoints, ndata.ijk_to_xyz_transform )

    dvals = dmap.interpolated_values ( npoints, dmap.openState.xform )
    #dvals = numpy.where ( dvals > threshold, dvals, numpy.zeros_like(dvals) )
    #nze = numpy.nonzero ( dvals )

    nmat = dvals.reshape( (nn3,nn2,nn1) )

    ndata = VolumeData.Array_Grid_Data ( nmat, nO, nstep, dmap.data.cell_angles )

    if fout == None :
        try : nv = VolumeViewer.volume.add_data_set ( ndata, None )
        except : nv = VolumeViewer.volume.volume_from_grid_data ( ndata )
        dmap_base = os.path.splitext(dmap.name)[0]
        dmap_path = os.path.splitext (dmap.data.path)[0]
        nv.name = dmap_base + "_masked"
        nv.openState.xform = dmap.openState.xform
        return nv

    else :

        from VolumeData import save_grid_data
        #d = self.grid_data()
        format = save_grid_data(ndata, fout, None, {}, False)
        #print " - saved data"





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

        r.isProt = r.type in protein3to1
        r.isNA = r.type in nucleic3to1

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
