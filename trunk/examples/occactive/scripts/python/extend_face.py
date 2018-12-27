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

from OCC.gp import *
from OCC.BRepBuilderAPI import *
from OCC.BRepOffsetAPI import *
from OCC.TColgp import *
from OCC.GeomAPI import *
from OCC.GeomAbs import *
from OCC.Geom2d import *
from OCC.BRepLib import *
from OCC.Quantity import *
from OCC.BRepPrimAPI import *
from OCC.BRepAlgoAPI import *
from OCC.TopOpeBRepTool import *
from OCC.TopAbs import *
from OCC.TopExp import *
from OCC.TopoDS import *
        
from OCC.Display.SimpleGui import init_display

display, start_display, add_menu, add_function_to_menu = init_display()

def build():

    Edge1 = BRepBuilderAPI_MakeEdge(gp_Pnt(0, 0, 0), gp_Pnt(100, 0, 0))
    Edge2 = BRepBuilderAPI_MakeEdge(gp_Pnt(100, 0, 0), gp_Pnt(100, 100, 0))
    Edge3 = BRepBuilderAPI_MakeEdge(gp_Pnt(100, 100, 0), gp_Pnt(0, 100, 0))
    Edge4 = BRepBuilderAPI_MakeEdge(gp_Pnt(0, 100, 0), gp_Pnt(0, 0, 0))

    wire1 = BRepBuilderAPI_MakeWire(Edge1.Edge(), Edge2.Edge(), Edge3.Edge(), Edge4.Edge())
    face1 = BRepBuilderAPI_MakeFace(wire1.Wire())

    display.DisplayColoredShape(face1.Face(), 'BLUE')
        
    v = gp_Vec(50, 0, 0);    
    face2 = BRepPrimAPI_MakeSweep(Edge2.Edge(), v);

    display.DisplayColoredShape(face2.Shape(), 'RED')
    
    
if __name__ == '__main__':
    build()
    start_display()