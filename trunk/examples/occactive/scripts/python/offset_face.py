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

    box = BRepPrimAPI_MakeBox(10., 20., 30.)
    display.DisplayColoredShape(box.Shape(), 'BLUE')
    
    ex1 = TopExp_Explorer(box.Shape(), TopAbs_FACE)
    while ex1.More():
        face = topods_Face(ex1.Current())
        
        b = BRepOffsetAPI_MakeOffsetShape(face, 2, 0);
        b.Build();
    
        display.DisplayColoredShape(b.Shape(), 'RED')
        ex1.Next()
    
    
if __name__ == '__main__':
    build()
    start_display()