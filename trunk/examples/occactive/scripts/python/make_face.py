##Copyright 2009-2015 Thomas Paviot (tpaviot@gmail.com)
##
##This file is part of pythonOCC.
##
##pythonOCC is free software: you can redistribute it and/or modify
##it under the terms of the GNU Lesser General Public License as published by
##the Free Software Foundation, either version 3 of the License, or
##(at your option) any later version.
##
##pythonOCC is distributed in the hope that it will be useful,
##but WITHOUT ANY WARRANTY; without even the implied warranty of
##MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
##GNU Lesser General Public License for more details.
##
##You should have received a copy of the GNU Lesser General Public License
##along with pythonOCC.  If not, see <http://www.gnu.org/licenses/>.

import math

from OCC.gp import (gp_Pnt, gp_Sphere, gp_Ax3, gp_Dir, gp_Circ, gp_Ax2,
                    gp_Pnt2d, gp_Dir2d)
from OCC.BRepBuilderAPI import (BRepBuilderAPI_MakeEdge,
                                BRepBuilderAPI_MakeFace,
                                BRepBuilderAPI_MakeWire)
from OCC.TColgp import TColgp_Array2OfPnt
from OCC.GeomAPI import GeomAPI_PointsToBSplineSurface
from OCC.GeomAbs import GeomAbs_C2
from OCC.Geom2d import Geom2d_Line
from OCC.BRepLib import breplib_BuildCurves3d
from OCC.Quantity import Quantity_Color, Quantity_NOC_PINK

from OCC.Display.SimpleGui import init_display
display, start_display, add_menu, add_function_to_menu = init_display()


def face():

    #The brown face
    circle = gp_Circ(gp_Ax2(gp_Pnt(0, 0, 0), gp_Dir(1, 0, 0)), 80)
    Edge1 = BRepBuilderAPI_MakeEdge(circle, 0, math.pi)
    Edge2 = BRepBuilderAPI_MakeEdge(gp_Pnt(0, 0, -80), gp_Pnt(0, -10, 40))
    Edge3 = BRepBuilderAPI_MakeEdge(gp_Pnt(0, -10.00001, 40), gp_Pnt(0, 0, 80))

    ##TopoDS_Wire YellowWire
    MW1 = BRepBuilderAPI_MakeWire(Edge1.Edge(), Edge2.Edge(), Edge3.Edge())

    print MW1.IsDone()
    
    if MW1.IsDone():
        yellow_wire = MW1.Wire()
    brown_face = BRepBuilderAPI_MakeFace(yellow_wire)
	
    display.DisplayColoredShape(brown_face.Face(), 'BLUE')

if __name__ == '__main__':
    face()
    start_display()