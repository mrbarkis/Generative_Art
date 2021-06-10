import networkx as nx
import numpy as np
import random
import matplotlib.pyplot as plt
from numpy.random import RandomState
import copy
from shapely.geometry import LineString
from itertools import combinations


def node_distance(Graph, n1, n2):
    """Calculate the distance between nodes n1 and n2 of Graph"""
    a = Graph.nodes[n1]["pos"]
    b = Graph.nodes[n2]["pos"]
    return np.linalg.norm(a - b)

def node_pos_dict(Graph):
    """Make a dictionary {node: [x_coordinate, y_coordinate]}."""
    return {k: v.tolist() for k, v in Graph.nodes.data("pos")}

def memoize(func):
    cache = dict()

    def memoized_func(*args):
        if args in cache:
            return cache[args]
        result = func(*args)
        cache[args] = result
        return result

    return memoized_func


# Cache distance calculations
node_distance = memoize(node_distance)


def get_length(Graph):
    """Calculate the total length of the graph """
    return sum([node_distance(Graph, n1, n2) for (n1, n2) in Graph.edges])
    
def get_nearby(Graph,node, reach=1):
    return list(Graph.graph['hash_grid'].get_nearby_cities(node, reach))


def get_rnd_graph(N, size=(640, 360), margin=(0, 0, 0, 0), seed=123):
    """ Make a random graph with N nodes"""
    #np.random.seed(seed)
    random.seed(seed)
    top, bottom, left, right = margin
    width, height = size
    G = nx.Graph()
    G.add_nodes_from(
        [
            (
                n,
                {
                    "pos": np.asarray(
                        [
                            random.randint(left, width - right),
                            random.randint(top, height - bottom),
                        ]
                    )
                },
            )
            for n in range(N)
        ]
    )

    G.add_edges_from([(n, (n + 1) % N) for n in range(N)])
    G.graph["length"] = get_length(G)
    G.graph["N"] = N
    return G


def get_rnd_segment(Graph, N_seg, reach, prng):
    """Select a random segment from graph

    Args:
        N_seg (int): Number of nodes in the segment.

    Returns:
        seg (list): List of nodes in the segment.
        ends (list): The end nodes of the segment.
        nbrs (list): Current neighbors of the segment.
        new_nbrs (list): New neighbors for the segment.
    """
    # Sample random segment from the Graph:
    seg = set()
    n = prng.randint(low=0, high=Graph.graph["N"])
    while len(seg) < N_seg:
        seg.add(n)
        n = prng.choice(list(Graph.neighbors(n)))

    # Find the ends and neighbours of the segment:
    nbrs = []
    ends = []
    for n in seg:
        for nbr in Graph.neighbors(n):
            if nbr not in seg:
                nbrs.append(nbr)
                ends.append(n)

    # Sample new neighbours for the segment:
    #n1 = np.random.randint(low=0, high=Graph.graph["N"])
    rnd_end = prng.choice(ends)
    nearby = set(get_nearby(Graph, rnd_end, reach)).difference(seg)
   
    while len(nearby) == 0: #no nodes nearby
        #print('reach increased')
        reach = 2 * reach
        nearby = set(get_nearby(Graph, rnd_end, reach)).difference(seg)
    
    n1 = prng.choice(list(nearby))
    #    n1 = np.random.randint(low=0, high=Graph.graph["N"])

    n2 = prng.choice(list(Graph.neighbors(n1)))
    
    while n2 in seg:
        n2 = prng.choice(list(Graph.neighbors(n1)))

    new_nbrs = [n1, n2]

    return (list(seg), ends, nbrs, new_nbrs)


def remove_segment(Graph, seg, ends, nbrs, new_nbrs):
    """ Remove segment from Graph."""
    edge1, edge2 = list(zip(nbrs, ends))
    # Remove the segment:
    Graph.remove_edge(*edge1)
    Graph.remove_edge(*edge2)
    # Connect the resulting loose ends:
    Graph.add_edge(*nbrs)


def add_segment(Graph, seg, ends, nbrs, new_nbrs):
    """ Connect segment to new_nbrs."""
    edge1, edge2 = list(zip(new_nbrs, ends))
    # Break the Graph between new_nbrs:
    Graph.remove_edge(*new_nbrs)
    # Connect new_nbrs with the segment:
    Graph.add_edge(*edge1)
    Graph.add_edge(*edge2)


def move_segment(Graph, seg, ends, nbrs, new_nbrs):
    """ Move segment between nbrs to new_nbrs."""
    remove_segment(Graph, seg, ends, nbrs, new_nbrs)
    add_segment(Graph, seg, ends, nbrs, new_nbrs)


def diff_E(Graph, seg, ends, nbrs, new_nbrs):
    """ Calculate the change in total length, if the segment is moved."""
    nbr1_to_nbr2, end1_to_nbr1, end2_to_nbr2 = (
        node_distance(Graph, nbrs[0], nbrs[1]),
        node_distance(Graph, nbrs[0], ends[0]),
        node_distance(Graph, nbrs[1], ends[1]),
    )
    
    new_nbr1_to_nbr2, new_end1_to_nbr1, new_end2_to_nbr2 = (
        node_distance(Graph, new_nbrs[0], new_nbrs[1]),
        node_distance(Graph, new_nbrs[0], ends[0]),
        node_distance(Graph, new_nbrs[1], ends[1]),
    )
    
    diff_E = (
        -new_nbr1_to_nbr2
        + nbr1_to_nbr2
        + new_end1_to_nbr1
        - end1_to_nbr1
        + new_end2_to_nbr2
        - end2_to_nbr2
    )

    return diff_E
    
    
def plot_Graph(ax, Graph, kT, kTs, N_walk=0):
    kT_min, kT_max = kTs[-1], kTs[0]
    ax.cla()
    #ax.set_title(
    #    "T = %.f, Walk = %.f, Length = %.2f" % (T, N_walk, Graph.graph["length"])
    #)
    nx.draw(
        Graph,
        pos=node_pos_dict(Graph),
        with_labels=False,
        node_size=10,
        ax=ax,
        node_color=[kT] * Graph.graph["N"],
        cmap=plt.cm.plasma,
        vmin=kT_min,
        vmax=kT_max,
        edge_color=[kT] * Graph.graph["N"],
        edge_cmap=plt.cm.plasma,
        edge_vmin=kT_min,
        edge_vmax=kT_max
    )


def plot_log(ax, track_kT, track_E, kTs, N_cities=233):
    kT_min, kT_max = kTs[-1] / 1000, kTs[0] / 1000
    ax.cla()
    #ax.set_ylim(top=60)
    ax.set(xlim=(kT_max, kT_min), title= 'avg. distance between cities vs. available "thermal energy"', ylabel="E / N [km]", xlabel="kT [km]")
    ax.tick_params(
        axis="both",
        which="both",
        bottom=True,
        labelbottom=True,
        left=True,
        labelleft=True,
    )
    ax.scatter(np.asarray(track_kT) / 1000, np.asarray(track_E) / N_cities, cmap=plt.cm.plasma, c=np.asarray(track_kT) / 1000, vmin=kT_min, vmax=kT_max)


def path_to_list(G):
    node_coordinates = node_pos_dict(G)
    visited = set()
    q = []
    current_node = 0
    while current_node not in visited:
        visited.add(current_node)
        q.append(node_coordinates[current_node])

        nbrs = list(G.neighbors(current_node))
        if nbrs[0] not in visited:
            current_node = nbrs[0]
        elif nbrs[1] not in visited:
            current_node = nbrs[1]

    
    return q
    
def reset(G_start, G_best, seed):
    return (copy.deepcopy(G_start),
            copy.deepcopy(G_best),
            RandomState(seed),
            [],
            [])


def get_coordinates(Graph, node):
    return tuple(Graph.nodes[node]['pos'])


def count_intersections(Graph):
    count = 0
    for edge_1, edge_2 in combinations(Graph.edges(), 2):
        segment_1 = LineString([get_coordinates(Graph, edge_1[n]) for n in range(2)])
        segment_2 = LineString([get_coordinates(Graph, edge_2[n]) for n in range(2)])
        count += segment_1.intersects(segment_2)
    
    return count - Graph.graph['N'] # each segment intersects with its neighbour