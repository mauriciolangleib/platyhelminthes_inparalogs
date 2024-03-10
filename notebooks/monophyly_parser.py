def create_edges_monophyly_assemblies_at_most_many_to_one(tree, leaves_graph, target_assemblies):
    for node in tree.get_monophyletic(values = target_assemblies, target_attr = 'assembly'):
            # create only if the monophyletic group contains: a. only sequences of desired assemblies (and for all assemblies), and at most in realtionship 1-to-many
            # count for assemblies
            assemblies_in = [leaf.split('.')[0] for leaf in node.get_leaf_names()]
            # create dictionary with counts
            assembly_count_dict = {}
            for assembly in target_assemblies:
                assembly_count_dict.update({assembly: assemblies_in.count(assembly)})
            # first verifying that only target assemblies are in this monophyletic group
            assemblies_in_set = set(assemblies_in)
            target_assemblies_set = set(target_assemblies)
            if (len(target_assemblies_set.difference(assemblies_in_set)) == 0) and (len(assemblies_in_set.difference(target_assemblies_set)) == 0):
                # now checking that, at most, one value for assembly is greater than one
                total_count = len(assembly_count_dict.keys())
                counts_upper_one = len([count for count in assembly_count_dict.values() if count > 1])
                if (total_count - counts_upper_one) >= (total_count - 1):
                    # adding edges to target graph
                    node_leafs = node.get_leaf_names()
                    # connect allt the pairwise combinations of them
                    for pairwise_nodes in itertools.combinations(node_leafs, 2):
                        # get them
                        node_a, node_b = pairwise_nodes
                        # add edge
                        leaves_graph.add_edge(node_a, node_b)
                        
# refine the function function
def getting_assembly_monophyletic_groups(tree, target_assemblies):
    # create the graph
    tree_graph = nx.Graph()
    # add nodes
    [tree_graph.add_node(leaf) for leaf in tree]
    # create the assembly dictionary
    node2assembly_dict = {}
    [node2assembly_dict.update({leaf.name: leaf.name.split('.')[0]}) for leaf in tree]
    # setting this features into the tree
    for leaf in tree:
        leaf.add_features(assembly = node2assembly_dict.get(leaf.name, 'none'))
    # traverse the tree in a 'postorder' modality
    for node in tree.traverse('postorder'):
        if node != tree:
            tree.set_outgroup(node)
        # performing the algorithm
        create_edges_monophyly_assemblies_at_most_many_to_one(tree = tree, leaves_graph = tree_graph, target_assemblies = target_assemblies)
    # get connected components
    connected_comps = [list(group) for group in list(nx.connected_components(tree_graph)) if len(list(group)) > 1]
    # return result
    return connected_comps
