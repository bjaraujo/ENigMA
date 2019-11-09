from OCC.Display.SimpleGui import init_display
from OCC.BRepPrimAPI import BRepPrimAPI_MakeBox
 
display, start_display, add_menu, add_function_to_menu = init_display()

box = BRepPrimAPI_MakeBox(10., 20., 30.).Shape() 
display.DisplayShape(box, update=True)
start_display()