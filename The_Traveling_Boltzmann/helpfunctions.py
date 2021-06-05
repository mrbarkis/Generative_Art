import networkx as nx
import numpy as np
import random
import matplotlib.pyplot as plt



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


def get_rnd_segment(Graph, N_seg, reach):
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
    n = np.random.randint(low=0, high=Graph.graph["N"])
    while len(seg) < N_seg:
        seg.add(n)
        n = np.random.choice(list(Graph.neighbors(n)))

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
    rnd_end = np.random.choice(ends)
    nearby = set(get_nearby(Graph, rnd_end, reach)).difference(seg)
   
    while len(nearby) == 0: #no nodes nearby
        #print('reach increased')
        reach = 2 * reach
        nearby = set(get_nearby(Graph, rnd_end, reach)).difference(seg)
    
    n1 = np.random.choice(list(nearby))
    #    n1 = np.random.randint(low=0, high=Graph.graph["N"])

    n2 = np.random.choice(list(Graph.neighbors(n1)))
    
    while n2 in seg:
        n2 = np.random.choice(list(Graph.neighbors(n1)))

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
    
    
def plot_Graph(ax, Graph, T, T_range=(0,300), N_walk=0):
    Tmin, Tmax = T_range
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
        node_color=[T] * Graph.graph["N"],
        cmap=plt.cm.plasma,
        vmin=Tmin,
        vmax=Tmax,
        edge_color=[T] * Graph.graph["N"],
        edge_cmap=plt.cm.plasma,
        edge_vmin=Tmin,
        edge_vmax=Tmax,
    )


def plot_log(axs, temps, lengths, T_range =(0,300), k_b=2, N_cities=233):
    ax, ax_twinx, ax_twiny = axs
    Tmin, Tmax = T_range
    ax.cla()
    ax.set(xlim=(Tmax, Tmin - 1), ylabel="Total Distance (Energy) [km]", xlabel="Temperature [K]")
    ax.tick_params(
        axis="both",
        which="both",
        bottom=True,
        labelbottom=True,
        left=True,
        labelleft=True,
    )
    ax.scatter(temps, lengths, cmap=plt.cm.plasma, c=temps, vmin=Tmin, vmax=Tmax)
   
    # set up twiny axis
    ax_twiny.cla()
    ax_twiny.set(xlim=(Tmax, Tmin - 1), xlabel='Available "Thermal Energy" kT [km]')
    ax_twiny.set_xticks(ax.get_xticks())
    ax_twiny.set_xbound(ax.get_xbound())
    ax_twiny.set_xticklabels([int(k_b * T / 1000) for T in ax.get_xticks()])
    
     # set up twiny axis
    ax_twinx.cla()
    ax_twinx.set(ylabel='Average Distance Between Stops [km]')
    ax_twinx.set_yticks(ax.get_yticks())
    ax_twinx.set_ybound(ax.get_ybound())
    ax_twinx.set_yticklabels([int(e / N_cities) for e in ax.get_yticks()])


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