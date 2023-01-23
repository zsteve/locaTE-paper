import numpy as np
import networkx as nx
import itertools
import matplotlib

def draw(g, gene_names, thresh = None, lst = None, cmap_dict = {0 : "MyGrey", 1 : "MyGrey"}, layout_args = ""):
    # get weights from g
    edges,weights = zip(*nx.get_edge_attributes(g,'weight').items())
    weights = np.array(list(weights))
    weights /= np.max(weights)
    # node and edgelist
    nl = list(g.nodes())
    el = list(g.edges())
    if thresh is not None:
        edge_idx = np.where(weights > np.quantile(weights, thresh))[0]
        edgelist = [el[i] for i in edge_idx]
    elif lst is not None:
        edgelist = [tuple(i) for i in lst]
    node_idx = list(set(itertools.chain(*edgelist)))
    nodelist = [nl[i] for i in node_idx]
    # take subgraph
    g_sub = g.subgraph(nodelist)
    # now plot
    edges,weights = zip(*nx.get_edge_attributes(g_sub,'weight').items())
    _,ref = zip(*nx.get_edge_attributes(g_sub,'ref').items())
    weights = np.array(list(weights))
    weights /= np.max(weights)
    pos = nx.nx_agraph.graphviz_layout(g_sub, prog = 'fdp', args = layout_args)
    for (k, v) in cmap_dict.items():
        edgelist_color = [e for e in edgelist if g_sub.edges[e]['ref'] == k]
        edge_idx_color = np.where([e in edgelist_color for e in g_sub.edges()])[0]
        arrows = nx.draw_networkx_edges(g_sub, pos, 
                                        edgelist = edgelist_color, edge_color = weights[edge_idx_color], alpha = weights[edge_idx_color], 
                                        width = 2.5*weights[edge_idx_color], edge_cmap = matplotlib.colormaps[v], edge_vmin = -0.1, edge_vmax = 1.1, node_size = 600)
        for a, w in zip(arrows, weights[edge_idx_color]):
            # from https://stackoverflow.com/questions/67251763/how-to-obtain-non-rounded-arrows-on-fat-lines
            a.set_mutation_scale(20 + w)
            a.set_joinstyle('miter')
            a.set_capstyle('butt')
    nodes, centrality = zip(*nx.get_node_attributes(g_sub,'centrality').items())
    centrality /= np.max(centrality)
    nx.draw_networkx_nodes(g, pos, nodelist = nodes, node_color = centrality, cmap = matplotlib.colormaps["viridis"], alpha = 0.5, node_size = 600)
    nx.draw_networkx_labels(g, pos, labels = {x : gene_names[x] for (c, x) in zip(centrality, nodes)}, font_size = 14);
    return g_sub