import matplotlib.pyplot as plt
import networkx as nx


def get_colors(G: nx.classes.graph.Graph) -> list:
    # getting colors

    color_map = {6: "cyan", 8: "orange", 7: "magenta"}

    colors = []
    for atom_id in G.nodes():
        atom = G.nodes[atom_id]
        atomic_num = atom["atomic_num"]

        if atomic_num in color_map:
            colors.append(color_map[atomic_num])
        else:
            colors.append("gray")

    # getting colors

    return colors


def get_labels(G: nx.classes.graph.Graph) -> dict:
    """
    Input:
        Graph
    """
    atomic_nums = nx.get_node_attributes(G, "atomic_num")
    isotopes = nx.get_node_attributes(G, "isotope")
    num_element = {
        6: "C",
        7: "N",
        8: "O",
        0: "*",
        16: "S",
    }

    labels = {}

    for atom_id in atomic_nums:
        atomic_num = atomic_nums[atom_id]

        if atomic_num == 0:
            if isotopes[atom_id] != 0:
                label = "%d*" % (isotopes[atom_id])
            else:
                label = "*"
        else:
            label = num_element.get(atomic_num)
            if label is None:
                label = atomic_nums[atom_id]
        labels[atom_id] = label

    return labels


def draw_mol_graph(G: nx.classes.graph.Graph, out: str = "mol_graph.png") -> None:
    """

    Input:
        G - graph representing a molecule

    Process:
        We get labels for the atom nodes in the Graph
        based on their element (or dummyAtom status)
        We plot an image for a graph G


    """
    labels = get_labels(G)

    colors = get_colors(G)

    nx.draw(G, labels=labels, with_labels=True, node_color=colors, node_size=800)

    plt.show()
    plt.savefig(out)
    return
