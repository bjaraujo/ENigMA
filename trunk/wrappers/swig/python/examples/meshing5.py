from ENigMA import ENigMA

import tkinter

edgeMesh = ENigMA.CMshMeshDouble()

node1 = ENigMA.CMshNodeDouble(0.0, 0.0, 0.0)
node2 = ENigMA.CMshNodeDouble(1.0, 0.0, 0.0)
node3 = ENigMA.CMshNodeDouble(1.0, 1.0, 0.0)
node4 = ENigMA.CMshNodeDouble(0.0, 1.0, 0.0)

edgeMesh.addNode(1, node1)
edgeMesh.addNode(2, node2)
edgeMesh.addNode(3, node3)
edgeMesh.addNode(4, node4)

element1 = ENigMA.CMshElementDouble(ENigMA.ET_BEAM)
element1.addNodeId(1)
element1.addNodeId(2)

element2 = ENigMA.CMshElementDouble(ENigMA.ET_BEAM)
element2.addNodeId(2)
element2.addNodeId(3)

element3 = ENigMA.CMshElementDouble(ENigMA.ET_BEAM)
element3.addNodeId(3)
element3.addNodeId(4)

element4 = ENigMA.CMshElementDouble(ENigMA.ET_BEAM)
element4.addNodeId(4)
element4.addNodeId(1)

edgeMesh.addElement(1, element1)
edgeMesh.addElement(2, element2)
edgeMesh.addElement(3, element3)
edgeMesh.addElement(4, element4)

triangleMesher = ENigMA.CMshTriangleMesherDouble()

meshSize = 0.026

edgeMesh.generateFaces(1E-3)
triangleMesher.remesh(edgeMesh, meshSize)
interiorPoints = ENigMA.StdVectorCGeoCoordinateDouble()
triangleMesher.generate(edgeMesh, 9999, interiorPoints, meshSize, meshSize, meshSize, 1E-6)

surfaceMesh = triangleMesher.mesh()
triangleMesher.flipEdges(surfaceMesh)
triangleMesher.relaxNodes(surfaceMesh)

print(surfaceMesh.nbNodes())
print(surfaceMesh.nbElements())

master = tkinter.Tk()

scale = 400

w = tkinter.Canvas(master, width=scale, height=scale)
w.pack()

for i in range(0, surfaceMesh.nbElements()):
    elementId = surfaceMesh.elementId(i)
    element = surfaceMesh.element(elementId)

    nodeId0 = element.nodeId(0)
    nodeId1 = element.nodeId(1)
    nodeId2 = element.nodeId(2)

    point0 = surfaceMesh.node(nodeId0)*scale
    point1 = surfaceMesh.node(nodeId1)*scale
    point2 = surfaceMesh.node(nodeId2)*scale
    
    w.create_polygon(point0.x(), point0.y(), point1.x(), point1.y(), point2.x(), point2.y(), fill="black", outline="red")
    
tkinter.mainloop()