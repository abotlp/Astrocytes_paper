# randomize_network.py

import os
import network_utilities

# Set the path to your network file in SIF format
network_file_in_sif = '/home/patricia/BIANA_phy/edge_scores.txt'

# Set the delimiter used in your SIF file (e.g., ' ' for space or '\t' for tab)
delim = '\t'  # Adjust as per your file's delimiter

# Specify the number of random networks you want to generate
n_sample = 100  # For example, generate 100 random networks

# Specify the output directory for random networks
random_networks_folder = 'sampled_graphs'
if not os.path.exists(random_networks_folder):
    os.makedirs(random_networks_folder)  # Python 2.7 compatible, no exist_ok argument

# Create random network files
print("Creating random networks")
sampling_prefix = os.path.join(random_networks_folder, "sampled_graph.sif.")

# Check if sampled networks already exist
if os.path.exists("{}{}".format(sampling_prefix, n_sample)):
    print("\tSampled networks exist, skipping this step!")
else:
    # Create the network from the SIF file
    g = network_utilities.create_network_from_sif_file(
        network_file_in_sif=network_file_in_sif,
        use_edge_data=True,
        delim=delim
    )

    # Generate random networks
    for i in range(1, n_sample + 1):
        g_sampled = network_utilities.randomize_graph(
            graph=g,
            randomization_type="preserve_topology_and_node_degree"
        )
        output_file = "{}{}".format(sampling_prefix, i)
        network_utilities.output_network_in_sif(g_sampled, output_file)
        print("\tGenerated random network: {}".format(output_file))
