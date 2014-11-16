class Vertex(object):
    def __init__(self, id):
        self.atmID = id
        self.edges = []
    def getAtmID(self):
        return self.atmID
    def getEdges(self):
        return self.edges
    def setAtmId(self, id):
        self.atmID = id
    def addEdge(self, ed):
        self.edges.append(ed)
        
class GraphRep(object):
    def __init__(self):
        self.vertices = {}
    def getVertices(self):
        return self.vertices
    def addVertex(self, id):
        self.vertices[id] = Vertex(id)
    def addEdge(self, a, b):
        self.vertices[a].addEdge(b)
        self.vertices[b].addEdge(a)
    def getVertexID(self, id):
        return self.vertices[id]
    def find_path(self, graph, start, end, path=[]):
        path = path + [start]
        if start == end: return [path]
        if start not in graph: return []
        paths = []
        for node in graph.getVertexID(start).getEdges():
            if node not in path:
                newpaths = self.find_path(self.getVertices(), start, end, path)
                for newpath in newpaths: paths.append(newpath)
        return paths
        
        
def find_path(graph, start, end, path=[]):
    path = path + [start]
    if start == end:
        return [path]
    if start not in graph:
        return []
    paths = []
    for node in graph[start]:
        if node not in path:
            newpaths = find_path(graph, node, end, path)
            for newpath in newpaths:
                paths.append(newpath)
    return paths

def find_loops(graph):
    loops = []
    sortedloops = []
    for i in graph:
        for j in graph[i]:
            paths = find_path(graph, i, j)
            if len(paths) <= 1: continue
            for k in paths:
                if len(k) <= 2: continue
                if sorted(k) not in sortedloops:
                    loops.append(k)
                    sortedloops.append(sorted(k))
    if len(loops) > 0:
        return loops
    return None