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

    #The blue face
    Edge1 = BRepBuilderAPI_MakeEdge(gp_Pnt(0, 0, 0), gp_Pnt(100, 0, 0))
    Edge2 = BRepBuilderAPI_MakeEdge(gp_Pnt(100, 0, 0), gp_Pnt(100, 100, 0))
    Edge3 = BRepBuilderAPI_MakeEdge(gp_Pnt(100, 100, 0), gp_Pnt(0, 100, 0))
    Edge4 = BRepBuilderAPI_MakeEdge(gp_Pnt(0, 100, 0), gp_Pnt(0, 0, 0))

    wire1 = BRepBuilderAPI_MakeWire(Edge1.Edge(), Edge2.Edge(), Edge3.Edge(), Edge4.Edge())
    face1 = BRepBuilderAPI_MakeFace(wire1.Wire())
	
    #display.DisplayColoredShape(face1.Face(), 'BLUE')
    
    #The red face
    Edge5 = BRepBuilderAPI_MakeEdge(gp_Pnt(50, -50, -50), gp_Pnt(50, -50, 100))
    Edge6 = BRepBuilderAPI_MakeEdge(gp_Pnt(50, -50, 100), gp_Pnt(50, 50, 100))
    Edge7 = BRepBuilderAPI_MakeEdge(gp_Pnt(50, 50, 100), gp_Pnt(50, 50, -50))
    Edge8 = BRepBuilderAPI_MakeEdge(gp_Pnt(50, 50, -50), gp_Pnt(50, -50, -50))
    
    wire2 = BRepBuilderAPI_MakeWire(Edge5.Edge(), Edge6.Edge(), Edge7.Edge(), Edge8.Edge())
    face2 = BRepBuilderAPI_MakeFace(wire2.Wire())
	
    #display.DisplayColoredShape(face2.Face(), 'RED')

    Edge9 = BRepBuilderAPI_MakeEdge(gp_Pnt(70, 20, -50), gp_Pnt(70, 20, 100))
    Edge10 = BRepBuilderAPI_MakeEdge(gp_Pnt(70, 20, 100), gp_Pnt(70, 50, 100))
    Edge11 = BRepBuilderAPI_MakeEdge(gp_Pnt(70, 50, 100), gp_Pnt(70, 50, -50))
    Edge12 = BRepBuilderAPI_MakeEdge(gp_Pnt(70, 50, -50), gp_Pnt(70, 20, -50))
    
    wire3 = BRepBuilderAPI_MakeWire(Edge9.Edge(), Edge10.Edge(), Edge11.Edge(), Edge12.Edge())
    face3 = BRepBuilderAPI_MakeFace(wire3.Wire())
    
    shape1 = BRepAlgoAPI_Fuse(face1.Face(), face2.Face())
    shape2 = BRepAlgoAPI_Fuse(shape1.Shape(), face3.Face())

    #display.DisplayColoredShape(shape2.Shape(), 'YELLOW')
    
    ex1 = TopExp_Explorer(shape2.Shape(), TopAbs_FACE)
    #ex1.Next()
    face = topods_Face(ex1.Current())        
    display.DisplayColoredShape(face, 'BLUE')

    ex2 = TopExp_Explorer(face, TopAbs_WIRE)
    while ex2.More():
        wire = topods_Wire(ex2.Current())
        ex3 = TopExp_Explorer(wire, TopAbs_EDGE)
        while ex3.More():
            edge = topods_Edge(ex3.Current())
            
            if (edge.Orientation() == TopAbs_INTERNAL):
                print "Internal"
            else:
                print "External"
                
            ex3.Next()
        ex2.Next()
    
if __name__ == '__main__':
    build()
    start_display()