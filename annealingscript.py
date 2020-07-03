import dwavebinarycsp
from dwave.system.samplers import DWaveSampler
from dwave.system.composites import EmbeddingComposite
import networkx as nx
from dwave.cloud import Client
import neal
import time
import numpy as np
import pandas as pd

def read_graph(name):
    num_vertices = 0
    num_edges = 0
    list_edges = []
    file = open("data/" + name , "r")
    for line in file:
        tokens = line.split()
        if(tokens[0] == "c"):
            #this line is a comment
            continue
        elif(tokens[0] == "p"and tokens[1] == "edge"):
            # save number of vertices and edges
            num_vertices = tokens[2]
            num_edges = tokens[3]
        else:
            # this is a line describing an edge
            temp = (tokens[1], tokens[2])
            list_edges.append(temp)
    file.close()
    return int(num_vertices), int(num_edges), list_edges

def isSampleInSamples(s, count):
    for i in range(len(count)):
        if (count[i]['sample']== s.sample):
            count[i]['num_occurrences']+= s.num_occurrences
            return True
    return False

def try_k_coloring(k_colors, graph_name, is_simulated = False):
    print("processing ", graph_name, "... ")
    num_vertices, num_edges, list_edges = read_graph(graph_name)
    vertices = [str(i+1) for i in range(num_vertices)]
    csp = dwavebinarycsp.ConstraintSatisfactionProblem(dwavebinarycsp.BINARY)
    one_color_configurations = set()

    def not_same_color(v1, v2):
        #constraint: not to adyacent nodes share same color
        return not (v1 and v2)

    for i in range(k_colors):
        one_color_configurations.add(tuple(1 if i == j else 0 for j in range(k_colors)))

    #constraint: just one color per vertex
    for vertex in vertices:
        variables = [vertex+"c"+str(i) for i in range(k_colors)]
        csp.add_constraint(one_color_configurations, variables)

    for edge in list_edges:
        v1, v2 = edge
        for i in range(k_colors):
            variables = [str(v1)+"c"+ str(i), str(v2)+"c" + str(i)]
            csp.add_constraint(not_same_color, variables)

    def plot_map(self):
        G = nx.Graph()
        G.add_nodes_from(vertices)
        G.add_edges_from(list_edges)
        # Translate from binary to integer color representation
        color_map = {}
        for province in vertices:
              for i in range(k_colors):
                if sample[province+"c"+str(i)]:
                    color_map[province] = i
        # Plot the sample with color-coded nodes
        node_colors = [color_map.get(node) for node in G.nodes()]
        nx.draw_circular(G, with_labels=True, node_color=node_colors, node_size=3000, cmap=plt.cm.rainbow)
        plt.show()

    bqm = dwavebinarycsp.stitch(csp)
    counts = []
    idSample = 0
    if not is_simulated:
        client_qpu = Client.from_config()
        client_cpu = Client.from_config(profile='prod')

        # Set up a solver using the local system’s default D-Wave Cloud Client configuration file
        # and sample 50 times
        sampler = EmbeddingComposite(DWaveSampler())         # doctest: +SKIP
        start_time = time.time()
        response = sampler.sample(bqm, num_reads=1024)         # doctest: +SKIP
        elapsed_time = time.time() - start_time
        for s in response.data():
            if not isSampleInSamples(s, counts):
                counts.append({'id': idSample,'sample':s.sample, 'num_occurrences':s.num_occurrences, 'energy':s.energy})
                idSample+=1

        # Plot the lowest-energy sample if it meets the constraints

        sample = next(response.samples())      # doctest: +SKIP
        if not csp.check(sample):              # doctest: +SKIP
            print("Failed to color map")
        print("execution time: ", elapsed_time)
        return counts,response
    else:
        sampler = neal.SimulatedAnnealingSampler()
        start_time = time.time()
        response = sampler.sample(bqm, num_reads=1024)
        elapsed_time = time.time() - start_time
        sample = next(response.samples())      # doctest: +SKIP

        for s in response.data():

            if not isSampleInSamples(s, counts):
                counts.append({'id': idSample,'sample':s.sample, 'num_occurrences':s.num_occurrences, 'energy':s.energy})
                idSample+=1

        if not csp.check(sample):              # doctest: +SKIP
            print("Failed to color map")
        print("execution time: ", elapsed_time)
        return counts, response

try_k_coloring(1, "graph1.txt", True)
try_k_coloring(2, "graph2.txt", True)
try_k_coloring(1, "graph3.txt", True)
try_k_coloring(2, "graph4.txt", True)
try_k_coloring(2, "graph5.txt", True)
try_k_coloring(2, "graph6.txt", True)
try_k_coloring(4, "graph.txt", True)
