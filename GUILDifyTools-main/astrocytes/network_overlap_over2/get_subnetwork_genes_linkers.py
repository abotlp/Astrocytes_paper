#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import sys
import networkx as nx
import numpy as np

########################
# HELPER FUNCTIONS
########################

def parse_seeds_file(seeds_file):
    """
    Reads `seeds.txt` which has 6 columns:
      0: BIANA ID
      1: UniProt ID
      2: Gene Symbol
      3: Description
      4: Gene ID
      5: Equivalent Entries

    We do NOT have a score. We'll assign a dummy score=1.0 to every seed.
    We'll treat BIANA ID (column 0) as the unique node key in our subnetwork.

    Return a dict:
      {
        user_entity_id : (uniprot, gene_symbol, score, rank, "seed")
      }
    """
    seeds_dict = {}
    rank_counter = 1

    with open(seeds_file, 'r') as f:
        header = f.readline()  # skip or store the header line
        for line in f:
            line = line.strip()
            if not line:
                continue

            parts = line.split('\t')
            if len(parts) < 6:
                # skip lines that don't have at least 6 columns
                continue

            user_entity_id = parts[0].strip()  # e.g. "71751"
            uniprot_id     = parts[1].strip()  # e.g. "O43707"
            gene_symbol    = parts[2].strip()  # e.g. "ACTN4"
            # columns 3,4,5 are description, gene ID, eq entries, not used for logic here

            # Assign a dummy score=1.0 for bridging logic
            score = 1.0
            # We'll mark all as "seed"
            seed_status = "seed"

            seeds_dict[user_entity_id] = (uniprot_id, gene_symbol, score, rank_counter, seed_status)
            rank_counter += 1

    return seeds_dict


def create_graph_instance(networks_dir, diamond=False):
    """
    Same as your original code, reading edge_scores.txt (or edge_diamond_file.txt).
    """
    edge_diamond_file = os.path.join(networks_dir, "edge_diamond_file.txt")
    edge_scores_file  = os.path.join(networks_dir, "edge_scores.txt")
    G = nx.Graph()

    if diamond and os.path.exists(edge_diamond_file):
        with open(edge_diamond_file, 'r') as f:
            for line in f:
                node1, node2 = line.strip().split('\t')
                G.add_node(node1)
                G.add_node(node2)
                G.add_edge(node1, node2)
    else:
        with open(edge_scores_file, 'r') as f:
            for line in f:
                line = line.strip()
                if not line:
                    continue
                parts = line.split()
                if len(parts) == 3:
                    node1, foo, node2 = parts
                    G.add_node(node1)
                    G.add_node(node2)
                    G.add_edge(node1, node2)

    return G


def get_shortest_paths_nodes(all_shortest_paths_file, nodes):
    """
    Reads precomputed 'all_shortest_paths.txt' for only the 'nodes' we care about.
    Return { nodeA : { nodeB: "nodeA,nodeX,nodeB" }, ... }
    """
    all_nodes_to_path = {}
    if not os.path.exists(all_shortest_paths_file):
        return all_nodes_to_path  # empty => we'll compute on the fly if needed

    with open(all_shortest_paths_file, 'r') as inp_fd:
        # first line: '-' plus node list
        line = inp_fd.readline().strip()
        if not line:
            return all_nodes_to_path
        nodes_network = line.split('\t')[1:]
        for line in inp_fd:
            fields = line.strip().split('\t')
            if not fields:
                continue
            curr_node = fields[0]
            if curr_node in nodes:
                paths = fields[1:]
                dictionary = dict(zip(nodes_network, paths))
                all_nodes_to_path[curr_node] = dictionary
            if len(all_nodes_to_path) == len(nodes):
                break
    return all_nodes_to_path


def decide_longest_connected_component(components, node_to_score):
    """
    Same logic: pick the connected component with the highest mean node score.
    """
    if not components:
        return set(), []
    max_len = 0
    lccs = []
    for comp in components:
        if len(comp) > max_len:
            max_len = len(comp)
            lccs.append(comp)

    if len(lccs) == 1:
        other_components = []
        for c in components:
            if c != lccs[0]:
                other_components.append(c)
        return set(lccs[0]), other_components
    else:
        best_mean_score = 0
        best_comp = []
        for comp in lccs:
            mean_score = np.mean([node_to_score[n] for n in comp if n in node_to_score])
            if mean_score > best_mean_score:
                best_mean_score = mean_score
                best_comp = comp
        others = []
        for c in components:
            if c != best_comp:
                others.append(c)
        return set(best_comp), others


def find_shortest_paths_two_components(cc, lcc, network, node_to_path=None, decide_by_score=True, node_to_score=None):
    """
    Finds shortest paths bridging cc -> lcc.
    If decide_by_score is True, pick path w/ best (max or mean) internal node score.
    """
    shortest_paths = []
    min_len = None
    for n1 in cc:
        for n2 in lcc:
            if node_to_path and (n1 in node_to_path) and (n2 in node_to_path[n1]):
                sp = node_to_path[n1][n2].split(',')
            else:
                sp = nx.shortest_path(network, source=n1, target=n2)
                if node_to_path is not None:
                    node_to_path.setdefault(n1, {})
                    node_to_path[n1][n2] = ",".join(sp)
            plen = len(sp)
            if min_len is None:
                min_len = plen
                shortest_paths = [sp]
            elif plen < min_len:
                min_len = plen
                shortest_paths = [sp]
            elif plen == min_len:
                shortest_paths.append(sp)

    if decide_by_score:
        best_max_score  = None
        best_mean_score = None
        final_paths = []
        for sp in shortest_paths:
            internals = sp[1:-1]
            if not internals:
                final_paths.append(sp)
                continue
            scores = [node_to_score[n] for n in internals if n in node_to_score]
            if len(scores) != len(internals):
                # some node not in node_to_score => skip or partially handle
                continue
            mx = max(scores)
            me = np.mean(scores)
            if best_max_score is None:
                best_max_score  = mx
                best_mean_score = me
                final_paths = [sp]
            else:
                if mx > best_max_score:
                    best_max_score  = mx
                    best_mean_score = me
                    final_paths = [sp]
                elif mx == best_max_score:
                    if me > best_mean_score:
                        best_mean_score = me
                        final_paths = [sp]
                    elif me == best_mean_score:
                        final_paths.append(sp)
        return final_paths, node_to_path
    else:
        return shortest_paths, node_to_path


#############################
# MAIN FUNCTION
#############################

def create_subnetwork_from_seeds_no_scores(
    network_json_file,
    seeds_file,
    networks_dir,
    output_dir,
    network_instance=None,
    species="9606",
    diamond=False
):
    """
    Builds a subnetwork from a seeds file that has *no* score column.
    We parse the seeds, assigning *score=1.0* to each. Then do the usual bridging logic.
    Finally, output a SIF, a JSON, and a PNG.
    """

    import matplotlib
    matplotlib.use('Agg')
    from matplotlib import pyplot as plt

    # 1) Parse seeds => user_entity_id -> (uniprot, gene_symbol, 1.0, rank, "seed")
    seeds_dict = parse_seeds_file(seeds_file)

    # Build node->score dictionary for bridging logic
    node_to_score = {}
    for ue_id, vals in seeds_dict.items():
        node_to_score[ue_id] = vals[2]  # the dummy 1.0

    # 2) Build or use existing network
    if network_instance is None:
        network_instance = create_graph_instance(networks_dir, diamond)

    # 3) Possibly read all_shortest_paths.txt
    all_shortest_paths_file = os.path.join(networks_dir, "all_shortest_paths.txt")
    if not os.path.exists(all_shortest_paths_file):
        print("all_shortest_paths.txt not found; bridging will compute on the fly.")
    else:
        print("Using existing shortest paths file:", all_shortest_paths_file)

    # 4) Build the subnetwork of seeds
    subG = nx.Graph()
    edges = []

    # We'll parse edge_scores.txt => keep edges if both ends are in seeds_dict
    edge_scores_file = os.path.join(networks_dir, "edge_scores.txt")
    if not os.path.exists(edge_scores_file):
        print("Error: edge_scores.txt not found in", networks_dir)
        return

    with open(edge_scores_file, 'r') as ef:
        for line in ef:
            line = line.strip()
            if not line:
                continue
            parts = line.split()
            if len(parts) != 3:
                continue
            n1, w, n2 = parts
            if n1 in seeds_dict:
                subG.add_node(n1)
                if n2 in seeds_dict:
                    subG.add_edge(n1, n2)
                    edges.append((n1, n2))
            if n2 in seeds_dict:
                subG.add_node(n2)

    # 5) Identify disconnected components among seeds => connect with bridging
    components = list(nx.connected_components(subG))
    if len(components) > 1:
        print("We have {} disconnected seed components. Attempting bridging...".format(len(components)))
        all_nodes_to_path = get_shortest_paths_nodes(all_shortest_paths_file, subG.nodes())
        lcc, others = decide_longest_connected_component(components, node_to_score)
        connectors = set()
        for cc in others:
            sps, all_nodes_to_path = find_shortest_paths_two_components(
                cc, lcc, network_instance,
                node_to_path=all_nodes_to_path,
                decide_by_score=True,
                node_to_score=node_to_score
            )
            for sp in sps:
                if len(sp) < 2:
                    continue
                internals = sp[1:-1]
                for i in range(len(sp)-1):
                    a = sp[i]
                    b = sp[i+1]
                    if not subG.has_node(a):
                        subG.add_node(a)
                    if not subG.has_node(b):
                        subG.add_node(b)
                    if not subG.has_edge(a, b):
                        subG.add_edge(a, b)
                        edges.append((a, b))
                for mid in internals:
                    connectors.add(mid)
            lcc.update(cc)
    else:
        connectors = set()

    # 6) Write subnetwork.sif
    subnetwork_file = os.path.join(output_dir, "subnetwork_genes_linkers.sif")
    with open(subnetwork_file, 'w') as sf:
        for (a, b) in edges:
            sf.write("%s interaction %s\n" % (a, b))
    

    # 7) Build JSON
    # seeds_dict => user_entity_id -> (uniprot, gene_symbol, 1.0, rank, "seed")
    rows = []
    for ue_id, vals in seeds_dict.items():
        uniprot = vals[0]
        gene    = vals[1]
        score   = vals[2]
        rank    = vals[3]
        tnode   = vals[4]  # "seed"

        size = score * 100.0
        if gene == '-':
            gene = uniprot

        row = ("{\n data: { id: \"%s\", size: %f, xref: \"%s\", label: \"%s\", "
               "score: %f, rank: \"%s\", type: \"%s\" } \n}" %
               (ue_id, size, uniprot, gene, score, rank, tnode))
        rows.append(row)

        if tnode == 'seed':
            subG.add_node(ue_id, color='#B0FF00', shape='8', size=size)
        else:
            subG.add_node(ue_id, color='yellow', shape='o', size=size)

    # bridging nodes
    for cnode in connectors:
        if cnode not in seeds_dict:
            # bridging node => treat as "linker"
            size = 100.0
            row = ("{\n data: { id: \"%s\", size: %f, xref: \"%s\", label: \"%s\", "
                   "score: 1.0, rank: \"9999\", type: \"linker\" } \n}" %
                   (cnode, size, cnode, cnode))
            rows.append(row)
            subG.add_node(cnode, color='gray', shape='o', size=size)

    node_str = ", ".join(rows)
    edge_list = []
    for (a, b) in edges:
        edge_txt = ("{\n data: { id: \"%s-%s\", source: \"%s\", target: \"%s\" } \n}" %
                    (a, b, a, b))
        edge_list.append(edge_txt)
    edge_str = ", ".join(edge_list)

    network_json_file = os.path.join(output_dir, network_json_file)
    with open(network_json_file, 'w') as jf:
        jf.write("[\n%s,\n%s\n]" % (node_str, edge_str))

    # 8) Plot
    import matplotlib
    matplotlib.use('Agg')
    from matplotlib import pyplot as plt

    nodePos = nx.fruchterman_reingold_layout(subG)
    nodeShapes = set((d["shape"] for (_, d) in subG.nodes(data=True)))

    for aShape in nodeShapes:
        shape_nodes = [n for (n, dd) in subG.nodes(data=True) if dd["shape"] == aShape]
        shape_colors= [subG.nodes[n]["color"] for n in shape_nodes]
        shape_sizes = [subG.nodes[n]["size"]  for n in shape_nodes]
        nodes_obj = nx.draw_networkx_nodes(subG,
                                           nodePos,
                                           node_shape=aShape,
                                           nodelist=shape_nodes,
                                           node_color=shape_colors,
                                           node_size=shape_sizes,
                                           linewidths=0.5)
        nodes_obj.set_edgecolor('black')
    nx.draw_networkx_edges(subG, nodePos, ax=None)
    plt.axis('off')
    network_png = os.path.join(output_dir, "network_plot_seeds.png")
    plt.savefig(network_png, format='PNG')

    print("Subnetwork SIF written to:", subnetwork_file)
    print("Cytoscape.js JSON written to:", network_json_file)
    print("PNG saved to:", network_png)


if __name__ == "__main__":
    # Example usage:
    seeds_file         = "seeds.txt"  # 6 columns (BIANA ID, UniProt, Gene Symbol, ...)
    network_json_file  = "subnetwork_genes_linkers.json"
    networks_dir       = "/home/patricia/BIANA_phy"
    output_dir         = "."
    species            = "9606"
    diamond            = False

    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    create_subnetwork_from_seeds_no_scores(
        network_json_file=network_json_file,
        seeds_file=seeds_file,
        networks_dir=networks_dir,
        output_dir=output_dir,
        network_instance=None,  # let it create from edge_scores
        species=species,
        diamond=diamond
    )
