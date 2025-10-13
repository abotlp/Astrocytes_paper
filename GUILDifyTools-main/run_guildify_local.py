import sys, os, re
import argparse
import math
import networkx as nx
import numpy as np
import operator
import shutil
from scipy.stats import fisher_exact
from scipy.stats import hypergeom
from statsmodels.sandbox.stats.multicomp import multipletests
import os, sys, re



def main():

    options = parse_user_arguments()
    run_guildify_local(options)

def parse_user_arguments(*args, **kwds):
    """
    Parses the arguments of the program.
    Examples:
    python /home/patricia/GUILDifyTools/scripts/run_guildify_local.py -i "TP53; BRCA1; BRCA2" -n /home/patricia/guildify2/data/BIANA_phy/9606/all -o /home/quim/data/guildify_local_outputs/example -s netscore -re 1 -it 1 -t 9606 -ns BIANA_phy -ti all -dr /var/www/html/sbi/guildp/data/drug_info.txt
    python /home/patricia/GUILDifyTools/scripts/run_guildify_local.py -i "TP53; BRCA1; BRCA2" -n /var/www/html/sbi/guildify2/data/BIANA/9606/all -o /home/quim/data/guildify_local_outputs/example_molt_simple_guildify2_quim -s netscore -re 1 -it 1 -t 9606 -ns BIANA -ti all -dr /var/www/html/sbi/guildify2/data/drug_info.txt
    python /home/patricia/GUILDifyTools/scripts/run_guildify_local.py -i "TP53; BRCA1; BRCA2" -n /var/www/html/sbi/guildify2/data/BIANA/9606/all -o /home/quim/data/guildify_local_outputs/example_simple_guildify2_diamond -s diamond -re 1 -it 1 -t 9606 -ns BIANA -ti all -dr /var/www/html/sbi/guildify2/data/drug_info.txt
    python /home/patricia/GUILDifyTools/scripts/run_guildify_local.py -i BRCA1 -n /var/www/html/sbi/guildp/data/BIANA_phy/9606/all -o /home/quim/data/guildify_local_outputs/example_BRCA1 -s netscore -re 1 -it 1 -t 9606 -ns BIANA_phy -ti all -dr /var/www/html/sbi/guildp/data/drug_info.txt
    """

    parser = argparse.ArgumentParser(
        description = "Generate the profiles of GUILDify",
        epilog      = "@oliva's lab 2023")
    parser.add_argument('-i','--input',dest='input',action = 'store',
                        help = """Input seeds (in gene symbol notation). If there is more than one seed, separate them using semicolons (e.g. "TP53; BRCA1; BRCA2")""")
    parser.add_argument('-n','--networks_dir',dest='networks_dir',action = 'store',
                        help = """ Folder containing all the networks for a specific database (e.g. BIANA_phy). """)
    parser.add_argument('-o','--output_dir',dest='output_dir',action = 'store',
                        help = """ Output directory. """)
    parser.add_argument('-s','--scoring',dest='scoring',action = 'store', default="netscore",
                        help = """ Scoring function name: netscore / netzcore / netshort / netcombo / diamond. """)
    parser.add_argument('-re','--repetitions',dest='repetitions',action = 'store', default=3,
                        help = """ Number of repetitions (for netscore and netcombo) """)
    parser.add_argument('-it','--iterations',dest='iterations',action = 'store', default=2,
                        help = """ Number of iterations (for netscore, netzcore and netcombo) """)
    parser.add_argument('-t','--taxid',dest='taxid',action = 'store', default='9606',
                        help = """ Taxonomy ID of the species of the network: 9606 (human) / 10090 (mouse) / 3702 (plant) / 6239 (worm) /
                                   7227 (fly) / 10116 (rat) / 4932 (yeast) """)
    parser.add_argument('-ns','--network_source',dest='network_source',action = 'store', default='BIANA_phy',
                        help = """ Name of the type of network considered (e.g. 'BIANA_phy' or 'BIANA_func') """)
    parser.add_argument('-ti','--tissue',dest='tissue',action = 'store', default='all',
                        help = """ Type of tissue. If all of them considered, enter 'all' """)
    parser.add_argument('-dr','--drug_info_file',dest='drug_info_file',action = 'store', default='/home/patricia/BIANA_phy/drug_info.txt',
                        help = """ File containing the drugs-target interactions considered by GUILDify """)
    parser.add_argument('-gu','--guild_executable_path',dest='guild_executable_path',action = 'store', default='/home/patricia/guild/guild',
                        help = """ Path to the executable program of GUILD """)
    parser.add_argument('-di','--diamond_executable_path',dest='diamond_executable_path',action = 'store', default='/home/patricia/DIAMOnD/DIAMOnD.py',
                        help = """ Path to the executable program of DIAMOND """)

    options=parser.parse_args()

    return options


#################
#################
# MAIN FUNCTION #
#################
#################

def run_guildify_local(options):
    """
    Generates GUILDify profiles
    """
    # Example of outputs
    #ls /var/www/html/sbi/guildify2/data/sessions/8075b9ad-3cec-4138-83c5-8a03d3a2fab8/

    #-------------------#
    #  PROCESS INPUTS   #
    #-------------------#

    # Create a directory for GUILD results
    output_dir = options.output_dir
    create_directory(output_dir)

    # Get associated gene_symbols
    words = options.input.split(";")
    if len(words) > 1:
        # Gene symbol list if query contains ";" 
        keywords = [ word.strip() for word in words if word != "" ]
    else:
        keywords = words
    print(keywords)

    # Read node info file
    networks_dir = options.networks_dir
    node_info_file = os.path.join(networks_dir, "node_info.txt")
    user_entity_to_values = get_user_entities_info(node_info_file)
    gene_symbol_to_user_entity = get_gene_symbol_to_user_entity_dict(user_entity_to_values)

    # Create seed info file
    seed_info_file = os.path.join(output_dir, "seeds.txt")
    seed_user_entity_ids = create_seed_info_file(keywords, user_entity_to_values, gene_symbol_to_user_entity, seed_info_file)

    # Save the species
    species = str(options.taxid).lower()
    species_file = os.path.join(output_dir, "species.txt")
    open(species_file, "w").write("{}\n".format(species))

    # Save the tissue
    tissue = str(options.tissue).lower()
    tissue_file = os.path.join(output_dir, "tissue.txt")
    open(tissue_file, "w").write("{}\n".format(tissue))

    # Save the network source
    network_source = str(options.network_source)
    network_source_file = os.path.join(output_dir, "network_source.txt")
    open(network_source_file, "w").write("{}\n".format(network_source))
    
    # Save scoring file
    scoring_functions = ['netcombo', 'netscore', 'netzcore', 'netshort', 'diamond', 'functionflow', 'netrank', 'netwalk', 'netprop']
    scoring_function_to_short = {'netscore': 'ns', 'netzcore': 'nz', 'netshort': 'nd', 'diamond': 'di', 'functionflow': 'ff', 'netrank': 'nr', 'netwalk': 'nw','netprop': 'np'}
    if options.scoring.lower() in scoring_functions:
        scoring = [options.scoring.lower()]
        if "netcombo" in scoring:
            scoring = ['netscore', 'netzcore', 'netshort', 'diamond', 'functionflow', 'netrank', 'netwalk', 'netprop']
            
           

    else:
        print('Scoring function {} does not exist. Introduce any of the following scoring functions: {}'.format(options.scoring.lower(), ', '.join(scoring_functions)))
        sys.exit(10)
    
    iterations = int(options.iterations)
    repetitions = int(options.repetitions)
    scoring_parameters = {}
    if 'netscore' in scoring:
        scoring_parameters['ns'] = [repetitions, iterations]
    if 'netzcore' in scoring:
        scoring_parameters['nz'] = [5]
    if 'netshort' in scoring:
        scoring_parameters['nd'] = None
    if 'diamond' in scoring:
        scoring_parameters['di'] = None
    if 'functionflow' in scoring:
        scoring_parameters['ff'] = [5, 1.0]
    if 'netrank' in scoring:
        scoring_parameters['nr'] = None
    if 'netwalk' in scoring:
        scoring_parameters['nw'] = None
    if 'netprop' in scoring:
        scoring_parameters['np'] = None

    print(scoring_parameters)
    scoring_file = os.path.join(output_dir, "scoring.txt")
    with open(scoring_file, "w") as scoring_fd:
        for scoring_type, parameters_type in scoring_parameters.iteritems():
            if parameters_type is not None:
                scoring_fd.write("{}\t{}\n".format(scoring_type, "\t".join(map(lambda x: str(x), parameters_type))))
            else:
                scoring_fd.write("{}\n".format(scoring_type))


    #-------------------------------------#
    #  SCORE NODES IN THE NETWORK (GUILD) #
    #-------------------------------------#

    # Prepare files before running GUILD
    prepare(seed_user_entity_ids, networks_dir, output_dir, default_seed_score=1.0, default_non_seed_score=0.01, n_sample_graph=100)

    # Run GUILD
    
    for scoring_type, parameters_type in scoring_parameters.iteritems():
        output_scores_file = os.path.join(output_dir, "output_scores.txt.{}".format(scoring_type))
        if not fileExist(output_scores_file):
            score_command = decide_scoring_commands(scoring_parameters=scoring_parameters, scoring_short=scoring_type, networks_dir=networks_dir, output_dir=output_dir, guild_executable_path=options.guild_executable_path, diamond_executable_path=options.diamond_executable_path)
            print('GUILDIFY INFO:\tExecuting GUILD with the following parameters: Score={} / Parameters={}.\n'.format(scoring_type, parameters_type))
            os.system(score_command)
            print('GUILDIFY INFO:\tGUILD has finished.\n')
        else:
            print('GUILDIFY INFO:\tThe scoring of the network with GUILD was already done for: Score={} / Parameters={}.\n'.format(scoring_type, parameters_type))


    #--------------------------------------------------------#
    #  PROCESS GUILD OUTPUT AND CALCULATE ADDITIONAL RESULTS #
    #--------------------------------------------------------#

    # Create guild scores file
    n_node, n_node_diamond, network_instance = create_output_file(networks_dir, output_dir, [scoring_function_to_short[sc] for sc in scoring], user_entity_to_values, seed_user_entity_ids, diamond=(scoring[0]=="diamond"))

    # Calculate the biological validation
    position_values, enrichment_values, seed_values, validated_nodes, non_validated_nodes, rank_range, enrichment_cutoff = calculate_biological_validation_using_Carlota(networks_dir, output_dir, user_entity_to_values, seed_user_entity_ids, network_instance, species, criteria_list=['fdr_bh', 'bonferroni'])
    # Modified code to include 5% and 10%
    top_enrichment = 'enrich_' + str(enrichment_cutoff)
    top_percentages = [1, 2, 5, top_enrichment]
    calculate_enrichment_of_top_nodes(networks_dir, output_dir, species, top_percentages, criteria_list=['fdr_bh', 'bonferroni'], diamond=(scoring[0]=="diamond"))

    top_percentages_drug = [1, 2, 5, top_enrichment]
    # Create network json and drug files
    for top in top_percentages_drug:
        get_network_json(str(top), networks_dir, output_dir, network_instance, options.drug_info_file, species, diamond=(scoring[0]=="diamond"))

    return

#######################
#######################
# SECONDARY FUNCTIONS #
#######################
#######################


def fileExist(file):
    """
    Checks if a file exists AND is a file
    """
    return os.path.exists(file) and os.path.isfile(file)


def create_directory(directory):
    """
    Checks if a directory exists and if not, creates it.
    """
    try:
        os.stat(directory)
    except:
        os.mkdir(directory)
    return


def get_user_entities_info(file_name, user_entity_ids = None):
    """
    Get the information of the user entities from 'file_name'.
    """
    user_entity_to_values = {}
    # Get information of these user entities
    if user_entity_ids is not None:
        user_entity_ids = set(user_entity_ids)
    header = True
    for line in open(file_name):
        if header:
            header = False
            continue
        #print(line)
        fields = line.strip().split('\t')
        #print(len(fields))
        user_entity_id, entry_ids, gene_symbols, descriptions, gene_ids, other_ids = line.strip().split("\t")

        # Remove isoforms
        entry_ids = remove_isoforms(entry_ids)

        if user_entity_ids is not None:
            if user_entity_id not in user_entity_ids:
                continue
        user_entity_to_values[user_entity_id] = [entry_ids, gene_symbols.strip(), descriptions, gene_ids, other_ids]
    return user_entity_to_values 


def remove_isoforms(entry_ids):
    """
    Remove the isoforms from a list of uniprot accessions.
    """
    # Check how many non-isoforms are there
    non_isoforms = [ entry for entry in entry_ids.split('; ') if len(entry.split('-')) == 1 ]

    # Remove the isoforms
    old_entry_ids = entry_ids.split('; ')
    new_entry_ids = set()
    for entry in old_entry_ids:
        if len(non_isoforms) > 0:
            if len(entry.split('-')) == 1:
                new_entry_ids.add(entry)
        else:
            new_entry_ids.add(entry)
    entry_ids = '; '.join(list(new_entry_ids))
    return entry_ids


def get_gene_symbol_to_user_entity_dict(user_entity_to_values):
    """
    Get a dictionary that maps gene symbols to BIANA user entities.
    """
    gene_symbol_to_user_entity = {}
    for user_entity_id in user_entity_to_values:
        [entry_ids, gene_symbols, descriptions, gene_ids, other_ids] = user_entity_to_values[user_entity_id]
        for gene_symbol in gene_symbols.split('; '):
            gene_symbol_to_user_entity.setdefault(gene_symbol, set()).add(user_entity_id)
    return gene_symbol_to_user_entity


def create_seed_info_file(keywords, user_entity_to_values, gene_symbol_to_user_entity, output_file):
    """
    Create file containing the seeds and the information about the seeds.
    """
    # Get the user entities of the seeds
    seed_user_entity_ids = set()
    for keyword in keywords:
        if keyword in gene_symbol_to_user_entity:
            user_entities_keyword = gene_symbol_to_user_entity[keyword]
            for user_entity in user_entities_keyword:
                seed_user_entity_ids.add(user_entity)
    # Write the output file containing information about each seed
    with open(output_file, "w") as out_fd:
        out_fd.write("BIANA ID\tUniProt ID\tGene Symbol\tDescription\tGene ID\tEquivalent Entries\n") # Write header
        for user_entity_id in seed_user_entity_ids:
            entry_ids, gene_symbols, descriptions, gene_ids, equivalent_entries = user_entity_to_values[user_entity_id]
            if entry_ids == "-":
                continue
            entry_ids = entry_ids.strip("; ")
            inner_values = [user_entity_id, entry_ids, gene_symbols, descriptions, gene_ids, equivalent_entries]
            out_fd.write("{}\n".format("\t".join(inner_values)))
    return seed_user_entity_ids


def create_edge_scores_as_node_scores_file(edges, node_to_score, edge_to_score, edge_scores_file, ignored_nodes = None, default_score = 0):
    """
    Creates edge score file from node association scores, intended comparing netshort with other algorithms without using other edge reliability/relevance score.
    """
    with open(edge_scores_file, 'w') as edge_scores_fd:
        for u,v in edges:
            if ignored_nodes is not None and u in ignored_nodes:
                score_u = default_score
            else:
                if node_to_score.has_key(u):
                    score_u = node_to_score[u]
                else:
                    score_u = default_score
            if ignored_nodes is not None and v in ignored_nodes:
                score_v = default_score
            else:
                if node_to_score.has_key(v):
                    score_v = node_to_score[v]
                else:
                    score_v = default_score
            weight = 1 # before it was default_score but makes no change since all included in edge_to_score has 1 as minimum anyways 
            if (u,v) in edge_to_score:
                weight = edge_to_score[(u,v)]
            elif (v,u) in edge_to_score:
                weight = edge_to_score[(v,u)]
            edge_scores_fd.write("%s %f %s\n" % (u, weight*(score_u + score_v) / 2, v))
    return


def prepare(seed_user_entity_ids, networks_dir, output_dir, default_seed_score=1.0, default_non_seed_score=0.01, n_sample_graph=100):
    """
    Prepare files needed to use GUILD.
    """

    # Get all user entities in the network
    user_entity_ids_all = set()
    edges = []
    seed_to_score = {}
    edge_to_score = {}
    edge_scores_file = os.path.join(networks_dir, "edge_scores.txt")
    for line in open(edge_scores_file, 'r'):
        id1, score, id2 = line.strip().split()
        user_entity_ids_all.add(id1)
        user_entity_ids_all.add(id2)
        edges.append((id1, id2))
        edge_to_score[(id1, id2)] = float(score)
    seed_user_entity_ids = set(seed_user_entity_ids)

    # Create initial node scores file
    node_scores_file = os.path.join(output_dir, "node_scores.txt")
    with open(node_scores_file, 'w') as node_scores_fd:
        for user_entity_id in user_entity_ids_all:
            if user_entity_id in seed_user_entity_ids:
                score = default_seed_score
                seed_to_score[user_entity_id] = default_seed_score
            else:
                score = default_non_seed_score
            node_scores_fd.write("{} {:.6f}\n".format(user_entity_id, score))



    # Create initial edge scores as node scores file for NetShort
    edge_scores_as_node_scores_file = os.path.join(output_dir, "edge_scores_as_node_scores.txt")
    create_edge_scores_as_node_scores_file(edges = edges, node_to_score = seed_to_score, edge_to_score = edge_to_score, edge_scores_file = edge_scores_as_node_scores_file, ignored_nodes = None, default_score = default_non_seed_score)

    # Create sampled graphs for NetZcore
    sampled_graph_dir = os.path.join(networks_dir, "sampled_graphs")
    sampled_graph_prefix = os.path.join(sampled_graph_dir, "sampled_graph.sif.")


def decide_scoring_commands(scoring_parameters, scoring_short, networks_dir, output_dir, guild_executable_path="/home/patricia/guild/guild", diamond_executable_path="/home/patricia/DIAMOnD/DIAMOnD.py"):
    """
    Define the scoring command.
    """
    edge_scores_file = os.path.join(networks_dir, "edge_scores.txt")
    edge_scores_spaces_file = os.path.join(networks_dir, "edge_scores_spaces.txt")
    edge_diamond_file = os.path.join(networks_dir, "edge_diamond_file.txt")
    sampled_graph_dir = os.path.join(networks_dir, "sampled_graphs")
    sampled_graph_prefix = os.path.join(sampled_graph_dir, "sampled_graph.sif.")
    node_scores_file = os.path.join(output_dir, "node_scores.txt")
    edge_scores_as_node_scores_file = os.path.join(output_dir, "edge_scores_as_node_scores.txt")
    output_scores_file = os.path.join(output_dir, "output_scores.txt")
    score_log_file = os.path.join(output_dir, "log.txt")
    seed_info_file = os.path.join(output_dir, "seeds.txt")
    random_walk_script_path = ("/home/patricia/guild/src/random_walk.r")
    if scoring_short == "ns":
        score_command = "{} -s s -n \"{}\" -e \"{}\" -o \"{}\" -r {} -i {} &> \"{}\"".format(guild_executable_path, node_scores_file, edge_scores_file, output_scores_file + ".ns", scoring_parameters["ns"][0], scoring_parameters["ns"][1], score_log_file + ".ns")
    elif scoring_short == "nz":
        score_command = "{} -s z -n \"{}\" -e \"{}\" -o \"{}\" -i {} -x {} -d \"{}\" &> \"{}\"".format(guild_executable_path, node_scores_file, edge_scores_file, output_scores_file+".nz", scoring_parameters["nz"][0], 100, sampled_graph_prefix, score_log_file+".nz")
    elif scoring_short == "nd":
        score_command = "{} -s d -n \"{}\" -e \"{}\" -o \"{}\" &> \"{}\"".format(guild_executable_path, node_scores_file, edge_scores_as_node_scores_file, output_scores_file + ".nd", score_log_file + ".nd")
    elif scoring_short == "di":
        score_command = 'python {} {} {} {} {} {} &> {}'.format(diamond_executable_path, edge_diamond_file, seed_info_file, '500', '1', output_scores_file+".di", score_log_file+".di")
    elif scoring_short == "ff":
        score_command = "{} -s f -n \"{}\" -e \"{}\" -o \"{}\" -i {} -t {} &> \"{}\"".format(guild_executable_path, node_scores_file, edge_scores_file, output_scores_file + ".ff", scoring_parameters["ff"][0], scoring_parameters["ff"][1], score_log_file + ".ff")
    elif scoring_short == "nr":
        score_command = "{} -s r -n \"{}\" -e \"{}\" -o \"{}\" &> \"{}\"".format(guild_executable_path, node_scores_file, edge_scores_file, output_scores_file + ".nr", score_log_file + ".nr")
    elif scoring_short == "nw":
        score_command = "R --slave --args {} {} {} < {} &> {}".format(node_scores_file, edge_scores_spaces_file, output_scores_file + ".nw", random_walk_script_path, score_log_file + ".nw")
    elif scoring_short == "np":
        score_command = "R --slave --args {} {} {} 1 < {} &> {}".format(node_scores_file, edge_scores_spaces_file, output_scores_file + ".np", random_walk_script_path, score_log_file + ".np")

    print(score_command)
    return score_command



def create_output_file(networks_dir, output_dir, scoring_methods_short, user_entity_to_values, seed_user_entity_ids, diamond=False):
    """
    Creates an output file for the GUILD calculation.
    """
    edge_scores_file = os.path.join(networks_dir, "edge_scores.txt")
    output_file = os.path.join(output_dir, "guild_scores.txt")
    output_scores_file = os.path.join(output_dir, "output_scores.txt")
    n_node = None
    n_node_diamond = None
    network_instance = None

    # Check if the file is already created and if it is a diamond file
    # If so, count all the nodes in the diamond file and return True
    if os.path.exists(output_file):
        n_node = -1 # header
        for line in open(output_file):
            n_node += 1
        if diamond:
            if os.path.exists(edge_scores_file):
                # Read edge scores file and 
                # get the number of nodes in the network
                f = open(edge_scores_file)
                nodes = set()
                for line in f:
                    node1, interaction, node2 = line.strip().split()
                    nodes.add(node1)
                    nodes.add(node2)
                n_node_diamond = len(nodes)
                f.close()
            else:
                n_node_diamond = n_node
        return n_node, n_node_diamond, network_instance

    # Check if the scoring files are created
    # If all of them exist, then we create the output file
    score_files = [output_scores_file + "." + x for x in scoring_methods_short]
    if all(os.path.exists(x) for x in score_files):

        # Always create DIAMOnD-specific output if DIAMOnD is part of the scoring methods
        if 'di' in scoring_methods_short:
            print("GUILDIFY INFO:\tCreating DIAMOnD-specific output file.\n")
            n_node_diamond, network_instance = create_output_file_for_diamond(networks_dir, output_dir, user_entity_to_values, seed_user_entity_ids)

        # Check if there is only one score or multiple
        # If multiple, we combine them!
        if len(score_files) == 1:
            output_scores_file = score_files[0]
        else:
            score_combined(score_files, output_scores_file)

        # Create a NetworkX graph instance
        network_instance = create_graph_instance(networks_dir, diamond)

        # Get the scores from the scores file
        values = []
        for line in open(output_scores_file):
            user_entity_id, score = line.strip().split()
            values.append((user_entity_id, float(score)))
        values.sort(lambda x, y: cmp(y[1], x[1]))
        user_entity_ids = map(lambda x: x[0], values)

        # Write the output file
        f = open(output_file, 'w')
        f.write("BIANA ID\tUniProt ID\tGene Symbol\tSeed\tDescription\tGene ID\tEquivalent Entries\tGUILD Score\tDegree\n") # % "\t".join(BIANAQuery.NODE_ATTRIBUTES[2:]))
        i = 0
        for user_entity_id, score in values:
            if user_entity_id not in user_entity_to_values:
                continue
            entry_ids, gene_symbols, descriptions, gene_ids, equivalent_entries = user_entity_to_values[user_entity_id]
            if entry_ids == "-":
                continue
            entry_ids = entry_ids.strip("; ")
            inner_values = [user_entity_id]
            is_seed = "0"
            if user_entity_id in seed_user_entity_ids:
                is_seed = "1"
            inner_values.extend([entry_ids, gene_symbols, is_seed, descriptions, gene_ids, equivalent_entries])
            inner_values.append(str(score))

            # Add the degree
            if user_entity_id in network_instance.nodes():
                inner_values.append(str(network_instance.degree([user_entity_id])[user_entity_id]))
            else:
                inner_values.append('0') # If it is not in the network, the degree is 0

            f.write("{}\n".format("\t".join(inner_values)))
            # skip ambigious entries
            i += 1
        n_node = i
        return n_node, n_node_diamond, network_instance
    return n_node, n_node_diamond, network_instance

def create_output_file_for_diamond(networks_dir, output_dir, user_entity_to_values, seed_user_entity_ids):
    """
    Creates an output file specifically for DIAMOnD calculations in a simplified format: user_entity_id and score.
    Preserves the original DIAMOnD output by copying it directly to the rank file before processing.
    """

    # Step 1: Copy the DIAMOnD output to the rank file (preserve original format and header)
    output_scores_file = os.path.join(output_dir, "output_scores.txt") + '.di'
    output_rank_file = os.path.join(output_dir, "output_rank.txt.di")
    shutil.copy(output_scores_file, output_rank_file)  # Copies the file directly

    # Step 2: Create a graph instance
    network_instance = create_graph_instance(networks_dir, diamond=True)

    # Step 3: Read the DIAMOnD output and gather the ranking
    values = []
    with open(output_scores_file, 'r') as output_scores_fd:
        for line in output_scores_fd:
            if line[0] == '#':  # Skip header or comment lines
                continue
            rank, user_entity_id = line.strip().split()
            values.append((user_entity_id, int(rank)))

    # Step 4: Write the simplified output file
    output_file = os.path.join(output_dir, "output_scores.txt.di")
    with open(output_file, 'w') as output_file_fd:

        # Step 4.1: Add seed nodes with score of 1.0
        for seed in seed_user_entity_ids:
            score = 1.0  # Assign the highest score for seed nodes
            output_file_fd.write("{}\t{:.6f}\n".format(seed, score))

        # Step 4.2: Add DIAMOnD ranked nodes with calculated score
        total_ranks = len(values)
        for user_entity_id, rank in values:
            # Calculate score based on rank
            score = (total_ranks - rank) / float(total_ranks)
            output_file_fd.write("{}\t{:.6f}\n".format(user_entity_id, score))

        # Step 4.3: Add nodes that were not included in DIAMOnD ranking
        user_entities_with_rank = set([user_entity_id for user_entity_id, _ in values])
        user_entities_without_rank = [
            user_entity_id for user_entity_id in network_instance.nodes()
            if user_entity_id not in user_entities_with_rank and user_entity_id not in seed_user_entity_ids
        ]

        for user_entity_id in user_entities_without_rank:
            score = 0.0  # Nodes without rank are assigned a score of 0.0
            output_file_fd.write("{}\t{:.6f}\n".format(user_entity_id, score))

    # Return the total number of nodes and the network instance
    return len(network_instance.nodes()), network_instance


def calculate_diamond_score(rank, number_of_ranks):
    """
    Calculates a score for DIAMOnD based on the rank:
    (number_of_ranks - rank) / (number_of_ranks)
    """
    return float((number_of_ranks - rank)) / float(number_of_ranks)

def create_graph_instance(networks_dir, diamond=False):
    """
    Parses the edge score file and adds the nodes and edges
    into a NetworkX Graph instance.
    """
    edge_diamond_file = os.path.join(networks_dir, "edge_diamond_file.txt")
    edge_scores_file = os.path.join(networks_dir, "edge_scores.txt")
    G=nx.Graph()
    if diamond:
        with open(edge_diamond_file, 'r') as file_fd:
            for line in file_fd:
                node1, node2 = line.strip().split('\t')
                G.add_node(node1)
                G.add_node(node2)
                G.add_edge(node1, node2)
    else:
        with open(edge_scores_file, 'r') as file_fd:
            for line in file_fd:
                node1, foo, node2 = line.strip().split()
                G.add_node(node1)
                G.add_node(node2)
                G.add_edge(node1, node2)
    network_instance = G
    return network_instance


def create_graph_instance_with_score(networks_dir, user_entity_id_to_score, diamond=False):
    """
    Parses the edge score file and adds the nodes and edges
    into a NetworkX Graph instance.
    """
    edge_diamond_file = os.path.join(networks_dir, "edge_diamond_file.txt")
    edge_scores_file = os.path.join(networks_dir, "edge_scores.txt")
    G=nx.Graph()
    if diamond:
        with open(edge_diamond_file, 'r') as file_fd:
            for line in file_fd:
                node1, node2 = line.strip().split('\t')
                G.add_node(node1, score=user_entity_id_to_score[node1])
                G.add_node(node2, score=user_entity_id_to_score[node2])
                G.add_edge(node1, node2)
    else:
        with open(edge_scores_file, 'r') as file_fd:
            for line in file_fd:
                node1, foo, node2 = line.strip().split()
                G.add_node(node1, score=user_entity_id_to_score[node1])
                G.add_node(node2, score=user_entity_id_to_score[node2])
                G.add_edge(node1, node2)
    network_instance = G
    return


def score_combined(scores_file_list, output_scores_file):
    """
    Calculates a combined score based on normalized scores of each scoring method (NetCombo)
    """
    node_to_scores = {}
    for scores_file in scores_file_list:
        node_to_score_inner = {}
        for line in open(scores_file):
            node, score = line.strip().split() 
            node_to_score_inner[node] = float(score)
        mean, sigma = calc_mean_and_sigma(node_to_score_inner.values())
        for node, score in node_to_score_inner.iteritems():
            node_to_scores.setdefault(node, []).append((score-mean)/sigma)
    values = []
    for node, scores in node_to_scores.iteritems():
        score = sum(scores) / len(scores)
        values.append((score, node))
    values.sort()
    min_v, max_v = min(values)[0], max(values)[0]
    f = open(output_scores_file, 'w')
    for score, node in values:
        score = (score-min_v) / (max_v-min_v)
        f.write("%s\t%f\n" % (node, score))
    f.close()
    return


def calc_mean_and_sigma(alist):
    def sq(x): return x*x
    def mean(list):
        n=len(list)
        return reduce(operator.add, list)/float(n)
    def sigma(list):
        n=len(list)
        m = mean(list)
        mSq=reduce(operator.add, map(sq, list))/float(n)
        s = mSq-sq(m)
        try:
            s = math.sqrt(s)
        except:
            if s<0.0000000001:
                s = 0
            else:
                raise
        return s
    return mean(alist), sigma(alist)


def get_user_entity_info_and_scores(guild_scores_file, n_start, n_end, tax_id):
    """
    Parse the guild_scores.txt file and get the information from user entities and scores.
    """
    values = []
    ueids = []
    i = 0
    f = open(guild_scores_file)
    f.readline() # skip header
    for line in f:
        i += 1
        if i < n_start:
            continue
        if i > n_end:
            break
        fields = line.strip().split("\t")
        if len(fields) == 8:
            [user_entity_id, entry_ids, gene_symbols, is_seed, descriptions, gene_ids, equivalent_entries, score] = fields
        elif len(fields) == 9:
            [user_entity_id, entry_ids, gene_symbols, is_seed, descriptions, gene_ids, equivalent_entries, score, degree] = fields

        new_entry_ids = [ ]
        last_entry = ''
        # Check how many non-isoforms are there
        for entry_id in entry_ids.split("; "):

            # Skip isoforms
            if '-' in entry_id:
                last_entry = entry_id
                continue

            if entry_id.startswith("KEGG:"):
                entry = '<a href="http://www.genome.jp/dbget-bin/www_bget?%s:%s" title="%s">%s</a>' % (BIANAQuery.TAX_ID_TO_KEGG_PREFIX[tax_id], entry_id[len("KEGG:"):], "\n".join(descriptions.split("; ")), entry_id)
            else:
                entry = '<a href="http://www.uniprot.org/uniprot/%s" title="%s">%s</a>' % (entry_id, "\n".join(descriptions.split("; ")), entry_id)
            new_entry_ids.append(entry)

        # If only isoforms, introduce the last one
        if len(new_entry_ids) < 1 and last_entry != '':
            if last_entry.startswith("KEGG:"):
                entry = '<a href="http://www.genome.jp/dbget-bin/www_bget?%s:%s" title="%s">%s</a>' % (BIANAQuery.TAX_ID_TO_KEGG_PREFIX[tax_id], entry_id[len("KEGG:"):], "\n".join(descriptions.split("; ")), entry_id)
            else:
                entry = '<a href="http://www.uniprot.org/uniprot/%s" title="%s">%s</a>' % (entry_id, "\n".join(descriptions.split("; ")), entry_id)
            new_entry_ids.append(entry)


        if gene_ids == "-":
            new_gene_ids = gene_ids
        else:
            #new_gene_ids = [ '<a href="http://www.ncbi.nlm.nih.gov/gene/%s">%s</a>' % (gene_id, gene_id) for gene_id in gene_ids.split("; ") ]
            new_gene_ids = [ '<a title="Equivalent Entries matched by BIANA:\n%s" href="http://www.ncbi.nlm.nih.gov/gene/%s">%s</a>' % (", ".join(equivalent_entries.split("; ")), gene_id, gene_id) for gene_id in gene_ids.split("; ") ]
        if int(is_seed) ==1:
            gene_symbol_text = "<span class='success'>seed</span> " + ",<br/>".join(gene_symbols.split("; ")) #+ "&nbsp;&nbsp;<span class='success'>seed</span>"
        else:
            gene_symbol_text = ",<br/>".join(gene_symbols.split("; "))
        #print(i, "<br/>".join(new_gene_ids), "<br/>".join(new_entry_ids), gene_symbol_text, "%.4f" % float(score))
        values.append((i, ",<br/>".join(new_gene_ids), ",<br/>".join(new_entry_ids), gene_symbol_text, "%.4f" % float(score))) # user_entity_id, "<br/>".join(descriptions.split("; "))
        #values.append((i, "<br/>".join(new_gene_ids), "<br/>".join(new_entry_ids), gene_symbol_text, "%.4f" % float(score), '<input type=checkbox name="%s" id ="%s" class="network_selection" title="Select node in the network"/>' %(user_entity_id, user_entity_id) )) # user_entity_id, "<br/>".join(descriptions.split("; "))
        ueids.append(user_entity_id)
    f.close()
    return values, ueids


def calculate_biological_validation_using_Carlota(networks_dir, output_dir, user_entity_to_values, seed_user_entity_ids, network_instance, species, criteria_list=['fdr_bh', 'bonferroni']):
    """
    Calculate the functional enrichment of the seeds, obtaining GO terms significantly enriched.
    Calculate how many of the nodes in the ranking are part of the enriched GO terms.
    """

    selected_method_to_file = {
        'GObp' : os.path.join(networks_dir, 'GObp2gene.{}.txt'.format(species)), 
        'GOmf' : os.path.join(networks_dir, 'GOmf2gene.{}.txt'.format(species)), 
        'Reactome' : os.path.join(networks_dir, 'Reactome2gene.{}.txt'.format(species))
    }

    # This is for Plotly.js
    javascript_output_file = os.path.join(output_dir, 'plotly_structures.txt')

    # Check if it exists
    # If it does not exist, calculate the biological validation
    # If it does not exist, calculate the biological validation
    if not os.path.exists(javascript_output_file):

        bv_instance = BiologicalValidation(None, None)

        # Use get_geneIDs_from_info_file to standardize seed IDs:
        node_info_file = os.path.join(networks_dir, 'node_info.txt')
        seed_to_values, _ = get_geneIDs_from_info_file(node_info_file, seed_user_entity_ids)

        # Extract standardized gene IDs from seeds
        seed_geneids = set()
        for seed in seed_to_values:
            for geneID in seed_to_values[seed]:
                seed_geneids.add(geneID)

        num_seeds = len(seed_geneids)

        #-------------------------------------------------#
        #  Calculate biological validation using Carlota  #
        #-------------------------------------------------#

        genes_in_terms_enriched = set()
        terms_enriched = set()

        for method in selected_method_to_file:

            # Check if the associations file is created. If not, create it!
            associations_file = selected_method_to_file[method]
            if not os.path.exists(associations_file):
                create_association_file(method, associations_file, network_instance)

            # Parse the association file
            dicbp, dicbp2, term2namebp = load_functional_terms(method, networks_dir, species)
            Nbp = len(dicbp.keys())
            NBbp = len(dicbp2.keys())


            for criteria in criteria_list:

                # Functional enrichment of seeds
                enrichment_file = os.path.join(output_dir, 'enrichment.{}.{}.seeds.txt'.format(method, criteria))
                test_passed_f1_bp, terms_l1_bp, term_to_values_bp = functional_enrichment(dicbp2, Nbp, list(seed_geneids), criteria)
                # Write the output file
                if len(seed_geneids) > 0:
                    print(seed_geneids)
                    write_enrichment_file_Carlota(term_to_values_bp, term2namebp, enrichment_file)
                else:
                    with open(enrichment_file, 'w') as out_fd:
                        out_fd.write('# Term ID\tTerm name\tnum of genes\tnum of total genes\tP-value\tP-value corrected')

                if criteria == 'fdr_bh' and method in ['GObp', 'GOmf']:
                    for term in term_to_values_bp:
                        pval, pval_corrected, x, m = term_to_values_bp[term]
                        if float(pval_corrected) < 0.5:
                            terms_enriched.add(term)
                            if term in dicbp2:
                                for gene in dicbp2[term]:
                                    genes_in_terms_enriched.add(gene)

        # Get the enriched genes in the ranking
        output_file = os.path.join(output_dir, "guild_scores.txt")
        user_entity_ranking, user_entity_to_geneids = bv_instance.parse_guild_scores_file(output_file)
        enriched_positions, non_enriched_positions, seed_ranking = bv_instance.get_enriched_genes_in_ranking(output_dir, user_entity_ranking, user_entity_to_geneids, seed_user_entity_ids, genes_in_terms_enriched)

        # Calculate enrichment in sliding window
        maximum_rank = 500 # We define the maximum rank to calculate as 500
        if num_seeds*2 > maximum_rank: # but if the number of seeds multiplied by 2 exceeds 500, then we define it as num_seedsx2 (but maximum 2000)
            maximum_rank = num_seeds*2
            if num_seeds > 2000:
                maximum_rank = 2000
        positions, pvalues, cutoff, cutoff_right_interval = bv_instance.calculate_enrichment_in_all_sliding_window_positions(output_dir, num_seeds, enriched_positions, non_enriched_positions, seed_ranking, maximum_rank=maximum_rank)

        # Output data structures (for Plotly.js)
        position_values, enrichment_values, seed_values, validated_nodes, non_validated_nodes, rank_range = bv_instance.output_plotly_information(positions, pvalues, seed_ranking, maximum_rank, cutoff_right_interval, output_dir)
        return position_values, enrichment_values, seed_values, validated_nodes, non_validated_nodes, rank_range, cutoff_right_interval

    # If the output file exists, recover it!
    else:
        # For Plotly.js
        counter = 0
        javascript_structures = [[],[],[],[],[],[],[]]
        with open(javascript_output_file, 'r') as javascript_output_fd:
            for line in javascript_output_fd:
                if line[0] == '#':
                    counter+=1
                    continue
                content = line.strip().split(',')
                javascript_structures[counter-1] += content
        position_values, enrichment_values, seed_values, validated_nodes, non_validated_nodes, cutoff_right_interval, rank_range = javascript_structures
        return position_values, enrichment_values, seed_values, validated_nodes, non_validated_nodes, rank_range, cutoff_right_interval[0]



def create_association_file(selected_method, network_instance, associations_file, node_info_file, data_dir="/home/patricia/BIANA_phy"):

    # Parse the node info file
    node_to_geneIDs, _ = get_geneIDs_from_info_file(node_info_file)

    # Get the Gene IDs of the network
    geneIDs_in_network = set()
    for node in network_instance.nodes():
        if node in node_to_geneIDs:
            for geneID in node_to_geneIDs[node]:
                geneIDs_in_network.add(geneID)

    # Parse the raw association files
    go_dir = os.path.join(data_dir, 'go')
    reactome_file = os.path.join(data_dir, 'go/NCBI2Reactome.txt')
    evidence_codes_reactome = ['IEA', 'TAS']
    evidence_codes_go = ['EXP', 'IDA', 'IMP', 'IGI', 'IEP', 'ISS', 'ISA', 'ISM', 'ISO']
    if selected_method == 'Reactome':
        # Evidence codes in Reactome: 'IEA', 'TAS'
        term2genes, term2name = parse_reactome(reactome_file, evidence_codes_reactome, species)
    elif selected_method == 'GObp':
        term2genes, term2name = parse_gene2go(go_dir, 'Process', evidence_codes_go, species)
    elif selected_method == 'GOmf':
        term2genes, term2name = parse_gene2go(go_dir, 'Function', evidence_codes_go, species)
    else:
        print('Incorrect method.')
        sys.exit(10)

    with open(associations_file, 'w') as out_fd:
        for termID in term2genes:
            term_name = term2name[termID]
            for geneID in term2genes[termID]:
                if geneID in geneIDs_in_network:
                    out_fd.write('{}\t{}\t{}\n'.format(termID, term_name, geneID))


def get_geneIDs_from_info_file(info_file, seeds=None):
    """
    Parse the functional enrichment file created by GOAtools
    and provide every row of interest.
    """

    node_to_geneIDs = {}
    geneID_to_nodes = {}
    if seeds:
        seeds = set([str(seed) for seed in seeds])

    with open(info_file, 'r') as f:
        for line in f:
            fields = line.strip().split('\t')
            if fields[0] == 'BIANA ID':
                continue
            node = fields[0]
            geneIDs = fields[4].split(';')

            if seeds:
                if node in seeds:
                    node_to_geneIDs[node] = geneIDs
                    for geneID in geneIDs:
                        geneID_to_nodes.setdefault(geneID, [])
                        geneID_to_nodes[geneID].append(node)
            else:
                node_to_geneIDs[node] = geneIDs
                for geneID in geneIDs:
                    geneID_to_nodes.setdefault(geneID, [])
                    geneID_to_nodes[geneID].append(node)

    return node_to_geneIDs, geneID_to_nodes


def parse_gene2go(go_dir, type_go, evidence_codes_used, taxID):

    if int(taxID) != 4932:
        go2name = {}
        go2geneids = {}
        gene2go_file = os.path.join(go_dir, 'gene2go')
        with open(gene2go_file, 'r') as inp_fd:
            first_line = inp_fd.readline()
            for line in inp_fd:
                #tax_id GeneID  GO_ID   Evidence    Qualifier   GO_term PubMed  Category
                tax_id, geneID, goID, evidence, qualifier, go_name, pubmed, category  = line.strip().split('\t')
                if category != type_go or tax_id != str(taxID):
                    continue
                if evidence in evidence_codes_used:
                    go2geneids.setdefault(goID, set()).add(geneID)
                    go2name[goID] = go_name
    else:
        taxID_to_gaf_file = {
            9606 : 'goa_human.gaf',
            10090 : 'mgi.gaf',
            10116 : 'rgd.gaf',
            6239 : 'wb.gaf',
            7227 : 'fb.gaf',
            3702 : 'tair.gaf',
            4932 : 'sgd.gaf',
        } 
        taxID_to_geneinfo_file = {
            9606 : 'Homo_sapiens.gene_info',
            10090 : 'Mus_musculus.gene_info',
            10116 : 'Rattus_norvegicus.gene_info',
            6239 : 'Caenorhabditis_elegans.gene_info',
            7227 : 'Drosophila_melanogaster.gene_info',
            3702 : 'Arabidopsis_thaliana.gene_info',
            4932 : 'Saccharomyces_cerevisiae.gene_info',
        } 
        from goatools.obo_parser import GODag
        from goatools.associations import read_gaf
        go_obo_file = os.path.join(go_dir, 'go-basic.obo')
        obodag = GODag(go_obo_file)
        gaf_file = os.path.join(go_dir, 'goa_files/{}'.format(taxID_to_gaf_file[int(taxID)]))
        geneinfo_file = os.path.join(go_dir, 'goa_files/{}'.format(taxID_to_geneinfo_file[int(taxID)]))
        gene2gos = read_gaf(gaf_file, taxids=[int(taxID)], evidence_set=evidence_codes_used)
        gos = set()
        for gene in gene2gos:
            for go in gene2gos[gene]:
                gos.add(go)
        print('Genes: {}, GOs: {}'.format(len(gene2gos), len(gos)))

        geneid2gos = {}
        go2name = {}
        with open(geneinfo_file, 'r') as inp_fd:
            first_line = inp_fd.readline()
            for line in inp_fd:
                fields = line.strip().split('\t')
                geneID = fields[1]
                gene = fields[5]
                #print(geneID, gene)
                if ':' in gene:
                    gene = gene.split(':')[1]
                if gene != '-' and geneID != '-':
                    if gene in gene2gos:
                        for go in gene2gos[gene]:
                            geneid2gos.setdefault(geneID, set()).add(go)
                            go_name = obodag[go].name
                            go2name[go] = go_name
                    else:
                        #print('Gene {} not found'.format(gene))
                        pass
        gos = set()
        for gene in geneid2gos:
            for go in geneid2gos[gene]:
                gos.add(go)
        print('GeneIDs: {}, GOs: {}'.format(len(geneid2gos), len(gos)))
        go2geneids = associations.get_b2aset(geneid2gos)
    return go2geneids, go2name


def parse_reactome(reactome_file, evidence_codes, taxID):
    """
    Parse Reactome to NCBI gene file.
    """
    taxID2species = {
        9606 : 'Homo sapiens',
        10090 : 'Mus musculus',
        10116 : 'Rattus norvegicus',
        6239 : 'Caenorhabditis elegans',
        7227 : 'Drosophila melanogaster',
        3702 : 'Arabidopsis thaliana',
        4932 : 'Saccharomyces cerevisiae'
    }
    selected_species = taxID2species[int(taxID)]
    reactome2genes = {}
    reactome2name = {}
    with open(reactome_file, 'r') as inp_fd:
        for line in inp_fd:
            fields = line.strip().split('\t')
            # 1 R-HSA-114608    https://reactome.org/PathwayBrowser/#/R-HSA-114608  Platelet degranulation  TAS Homo sapiens
            geneID, reactomeID, reactome_url, reactome_name, evidence_code, species_name = fields
            if evidence_code not in evidence_codes:
                continue
            if selected_species == species_name:
                reactome2genes.setdefault(reactomeID, set()).add(geneID)
                reactome2name[reactomeID] = reactome_name
    return reactome2genes, reactome2name


def load_functional_terms(selected_method, associations_dir, taxID):
    from collections import defaultdict
    file_method = {'GOmf': 'GOmf2gene.{}.txt'.format(taxID), 'GObp': 'GObp2gene.{}.txt'.format(taxID), 'Reactome': 'Reactome2gene.{}.txt'.format(taxID)}
    data_file = os.path.join(associations_dir, '{}'.format(file_method[selected_method]))
    dicbp = defaultdict(list)
    dicbp2 = defaultdict(list)
    term2name = {}
    fi3 = open(data_file)
    for line in fi3.readlines():
        gene = line.rstrip().split("\t")
        if len(gene) == 2:
            gobp = gene[0]
            geneid = gene[1]
        elif len(gene) == 3:
            gobp = gene[0]
            name = gene[1]
            geneid = gene[2]
            term2name[gobp] = name

        if gobp not in dicbp[geneid]:
            dicbp[geneid].append(gobp)

        if geneid not in dicbp2[gobp]:
            dicbp2[gobp].append(geneid)

    fi3.close()
    #print(selected_method, 'data loaded')
    return dicbp, dicbp2, term2name


def functional_enrichment(dicbp2, N, geneset_toenrich, criteria):
    """
    Calculate the functional enrichment using the method described by Carlota.
    """
    pvals = {}
    genes_enriched = {}
    go_enriched = False
    k = len(geneset_toenrich)
    terms_l = []
    test_passed_f = False
    term_to_values = {}
    for term in dicbp2.keys():
        m = len(dicbp2[term])
        xl = [y for y in dicbp2[term] if y in geneset_toenrich]
        x = len(xl)

        if x != 0:
            go_enriched = True
            xlist = []
            for i in range(x, m + 1):
                xlist.append(i)

            # calculation of the hypervalue
            # dhyper = robjects.r['dhyper']
            # xboh = robjects.IntVector(xlist)
            # dhypervalue = dhyper(xboh, m, (N - m), k, log=False)
            dhypervalue = hypergeom.pmf(xlist, N, m, k)
            # threshold of enrichment
            pvals[term] = sum(dhypervalue)
            genes_enriched[term] = [x, m] # Quim: addition to get genes enriched
    if go_enriched:
        pvals_values = list(pvals.values())
        terms = list(pvals.keys())
        pvals_corrected = multipletests(pvals_values, alpha=0.05, method=criteria, is_sorted=False,
                                        returnsorted=False)

        for i in range(0, len(terms)):
            if list(pvals_corrected[1])[i] < 0.05:
                # Quim: addition to get the p-values and genes enriched
                pval = pvals_values[i]
                pval_corrected = list(pvals_corrected[1])[i]
                term = terms[i]
                x, m = genes_enriched[term]
                term_to_values[term] = [pval, pval_corrected, x, m]
                ####################################
                test_passed_f = True
                terms_l.append(terms[i])
    return test_passed_f, terms_l, term_to_values


def write_enrichment_file_Carlota(term_to_values, term2name, enrichment_file):
    with open(enrichment_file, 'w') as out_fd:
        out_fd.write('# Term ID\tTerm name\tnum of genes\tnum of total genes\tP-value\tP-value corrected\n')
        for term, values in sorted(term_to_values.items(),key=lambda x: x[1][1],reverse=False):
            pval, pval_corrected, x, m = values
            if term in term2name:
                term_name = term2name[term]
            else:
                term_name = '-'
            out_fd.write('{}\t{}\t{}\t{}\t{}\t{}\n'.format(term, term_name, x, m, pval, pval_corrected))
    return


def calculate_enrichment_of_top_nodes(networks_dir, output_dir, species, top_percentages, criteria_list=['fdr_bh', 'bonferroni'], diamond=False):
    """
    Calculate the functional enrichment of the seeds, obtaining GO terms significantly enriched.
    Calculate how many of the nodes in the ranking are part of the enriched GO terms.
    """

    selected_method_to_file = {
        'GObp' : os.path.join(networks_dir, 'GObp2gene.{}.txt'.format(species)), 
        'GOmf' : os.path.join(networks_dir, 'GOmf2gene.{}.txt'.format(species)), 
        'Reactome' : os.path.join(networks_dir, 'Reactome2gene.{}.txt'.format(species))
    }

    # Check if the enrichment files exist
    # If one of them does not exist, we will calculate the functional enrichemnt
    calculate_enrichment = False
    for top_percentage in top_percentages:
        for selected_method in selected_method_to_file:
            for criteria in criteria_list:
                enrichment_file = os.path.join(output_dir, 'enrichment.{}.{}.{}.txt'.format(selected_method, criteria, str(top_percentage)))
                if not os.path.exists(enrichment_file):
                    calculate_enrichment = True

    if calculate_enrichment:

        # Parse GUILD scores file
        output_file = os.path.join(output_dir, "guild_scores.txt")
        user_entity_ranking, user_entity_to_geneids = parse_guild_scores_file(output_file)

        # Create the functional enrichment for the three methods
        for selected_method in selected_method_to_file:

            # Check if the associations file is created
            # If not, create it!
            associations_file = selected_method_to_file[selected_method]
            if not os.path.exists(associations_file):
                create_association_file(selected_method, associations_file)

            # Read associations file
            dicbp, dicbp2, term2name = load_functional_terms(selected_method, networks_dir, species)
            N = len(dicbp.keys())
            NB = len(dicbp2.keys())

            # For each percentage, calculate the functional enrichment
            for percentage in top_percentages:

                # Get the GeneIDs of the top nodes
                top_geneIDs = set()
                top_user_entity_id_to_values = get_top_scoring_ueids_and_values(networks_dir, output_dir, percentage, diamond=diamond)
                for user_entity in top_user_entity_id_to_values:
                    if user_entity in user_entity_to_geneids:
                        for geneID in user_entity_to_geneids[user_entity]:
                            if geneID == '-':
                                continue
                            top_geneIDs.add(geneID)

                for criteria in criteria_list:

                    # Calculate the functional enrichment
                    test_passed_f1, terms_l1, term_to_values = functional_enrichment(dicbp2, N, list(top_geneIDs), criteria)

                    # Write the output file
                    enrichment_file = os.path.join(output_dir, 'enrichment.{}.{}.{}.txt'.format(selected_method, criteria, str(percentage)))
                    if len(top_geneIDs) > 0:
                        write_enrichment_file_Carlota(term_to_values, term2name, enrichment_file)
                    else:
                        with open(enrichment_file, 'w') as out_fd:
                            out_fd.write('# Term ID\tTerm name\tnum of genes\tnum of total genes\tP-value\tP-value corrected')


        # Write the number of GO terms used in a file.
        # This will be used in the calculation of overlap of functions
        num_gos_file = os.path.join(networks_dir, "number_of_GOs.txt")
        if not os.path.exists(num_gos_file):
            with open(num_gos_file, 'w') as num_gos_fd:
                for method in selected_method_to_file:
                    dicbp, dicbp2, term2name = load_functional_terms(method, networks_dir, species)
                    gos = dicbp2.keys()
                    num_gos_fd.write('{}\t{}\n'.format(method, len(gos)))

    return


def node_top_scoring(node_to_vals, threshold, output_file):
    """
    Creates a profile with the most relevant nodes using a percentage threshold of the score provided by the user.

    Parameters:
        @node_to_vals:          Dictionary of all the nodes to a tuple containing (score, pvalue)
        @threshold:             Top percentage of the pvalue_file in which we will cut to obtain the most relevant nodes
        @output_file:           Resulting file which will contain the most relevant nodes
    """

    # Get top scoring, i.e. nodes that are in a given percentage of top score
    ntop=float(threshold)*len([x for x in node_to_vals.iteritems()])/100.0

    top_nodes = set()
    last_score = ''
    ii=0

    f = open(output_file, 'w')
    f.write('#Id Score P-value\n')

    # Now, write the profile with the top scoring nodes
    for node, vals in sorted(node_to_vals.items(),key=lambda x: x[1][0],reverse=True):
        score, pval = vals
        if ii < ntop:
            top_nodes.add(node)
            last_score = score # The last score is saved, so that if there is a score which is above the top threshold but has the same value than the last node, it is included as well                    
            f.write('{} {} {}\n'.format(node, score, pval))
            ii=ii+1
        else:
            if score == last_score: # Here, if a score above the threshold has the same score as the last, it is also recorded
                top_nodes.add(node)
                f.write('{} {} {}\n'.format(node, score, pval))
            else:
                break
    f.close()

    print('GUILDIFY INFO:\t{} file created.\n'.format(output_file))

    return


def edge_top_scoring(network_file, node_to_vals, threshold, output_file):
    """
    Creates profiles with the most relevant edges using the thresholds provided by the user.
    The profile is created selecting the top most relevant nodes, and creating the subsequent scored subnetwork.

    Parameters:
        @network_file:          File containing the scored network
        @node_to_vals:          Dictionary of all the nodes to a tuple containing (score, pvalue)
        @threshold:             Top percentage of the network in which we will cut to obtain the most relevant edges
        @output_file:           Resulting file which will contain the most relevant edges
    """

    # Create a network from the network file
    g = create_network_from_sif_file(network_file, use_edge_data=True)

    # Get top scoring, i.e. nodes that are in a given percentage of top score
    ntop=float(threshold)*len([x for x in node_to_vals.iteritems()])/100.0

    top_nodes = set()
    ii=0
    last_score = ''

    # Filter the nodes by top scoring nodes
    for node, vals in sorted(node_to_vals.items(),key=lambda x: x[1][0],reverse=True):
        score, pval = vals
        if ii < ntop:
            top_nodes.add(node)
            last_score = score # The last score is saved, so that if there is a score which is above the top threshold but has the same value than the last node, it is included as well
            ii=ii+1
        else:
            if score == last_score: # Here, if a score above the threshold has the same score as the last, it is also recorded
                top_nodes.add(node)
            else:
                break

    # Get subnetwork induced by top scoring nodes
    g_sub = get_subgraph(g, top_nodes)

    f = open(output_file, 'w')

    # The subnetwork induced by the top scoring nodes is scored and recorded in a file
    for u, v in g_sub.edges():
        score_u = node_to_vals[u][0]
        score_v = node_to_vals[v][0]
        score = (score_u + score_v) / 2
        f.write('{}\t{:f}\t{}\n'.format(u, score, v))

    f.close()

    print('GUILDIFY INFO:\t{} file created.\n'.format(output_file))

    return


def functional_top_scoring(obodag, geneid2gos, all_nodes, top_nodes, output_file, temp_file):
    """
    Creates profiles with the most relevant functions using the thresholds provided by the user.
    The profile is created selecting the top most relevant nodes, and computing a functional enrichment analysis of them.

    Parameters:
        @obodag:                Dictionary containing the Gene Ontology
        @geneid2gos:            Dictionary containing the equivalences from geneID to GOs
        @all_nodes:             A list containing all the nodes of the network
        @top_nodes:             A list containing the selected nodes to do the enrichment analysis
        @output_file:           Resulting file which will contain the most relevant edges
        @temp_file:             A file where the functional enrichment will be temporary calculated
    """

    # Get integer of background genes (because the geneIDs in gene2go are integers)
    new_all_nodes = []
    for x in all_nodes:
        if x == '-' or x == '':
            continue
        new_all_nodes.append(int(x))

    # Get integer of top genes  (because the geneIDs in gene2go are integers)
    new_top_nodes = []
    for x in top_nodes:
        if x == '-' or x == '':
            continue
        new_top_nodes.append(int(x))

    calculate_functional_enrichment_profile(obodag, geneid2gos, new_top_nodes, new_all_nodes, temp_file, output_file)
    #GOA.calculate_functional_enrichment_profile(obodag, geneid2gos, top_nodes, all_nodes, temp_file, output_file)

    print("GUILDIFY INFO:\t%s file created.\n" %(output_file))

    return


def score_network(network_file, node_to_vals, output_file):
    """
    Scores a network
    """

    # Create a network from the network file
    g = create_network_from_sif_file(network_file, use_edge_data=True)

    # Get all nodes
    all_nodes = node_to_vals.keys()

    f  = open(output_file, 'w')
    # Here, the network is scored and stored in a file
    for u, v in g.edges():
        score_u = node_to_vals[u][0]
        score_v = node_to_vals[v][0]
        score = (float(score_u) + float(score_v)) / 2
        f.write('{}\t{:f}\t{}\n'.format(u, score, v))
    f.close()
 

    return


def create_network_from_sif_file(network_file_in_sif, use_edge_data = False, delim = None, include_unconnected=True):
    """
    Creates a NetworkX graph object from a sif file
    """
    setNode, setEdge, dictDummy, dictEdge = get_nodes_and_edges_from_sif_file(network_file_in_sif, store_edge_type = use_edge_data, delim = delim)
    g=nx.Graph()
    if include_unconnected:
        g.add_nodes_from(setNode)
    if use_edge_data:
        for e,w in dictEdge.iteritems():
            u,v = e
            g.add_edge(u,v,{'w':w})
    else:
        g.add_edges_from(setEdge)
    return g


def get_nodes_and_edges_from_sif_file(file_name, store_edge_type = False, delim=None, data_to_float=True):
    """
    Parse sif file into node and edge sets and dictionaries
    returns setNode, setEdge, dictNode, dictEdge
    store_edge_type: if True, dictEdge[(u,v)] = edge_value
    delim: delimiter between elements in sif file, if None all whitespaces between letters are considered as delim
    """
    setNode = set()
    setEdge = set()
    dictNode = {}
    dictEdge = {}
    f=open(file_name)
    for line in f:
        if delim is None:
            words = line[:-1].split()
        else:
            words = line[:-1].split(delim)
        id1 = words[0]
        setNode.add(id1)
        if len(words) == 2:
            if data_to_float:
                score = float(words[1])
            else:
                score = words[1]
            dictNode[id1] = score
        elif len(words) == 3: 
            id2 = words[2]
            setNode.add(id2)
            setEdge.add((id1, id2))
            if store_edge_type:
                if data_to_float:
                    dictEdge[(id1, id2)] = float(words[1])
                else:
                    dictEdge[(id1, id2)] = words[1]
    f.close()
    if len(setEdge) == 0:
        setEdge = None
    if len(dictNode) == 0:
        dictNode = None
    if len(dictEdge) == 0:
        dictEdge = None
    return setNode, setEdge, dictNode, dictEdge


def get_subgraph(G, nodes):
    """
    NetworkX subgraph method wrapper
    """
    return G.subgraph(nodes)


def calculate_functional_enrichment_profile(obodag, geneid2go, top_nodes, all_nodes, temp_file, output_file):

    print("{N:,} annotated genes".format(N=len(geneid2go)))

    # Get a dictionary of GOs to their geneids (it will be used for the Log of Odds calculus)
    #from goatools import associations
    #go2geneids = associations.get_b2aset(geneid2go)

    # Define the GOEnrichmentStudy object
    from goatools.go_enrichment import GOEnrichmentStudy
    goeaobj = GOEnrichmentStudy(
            all_nodes, # List of background genes
            geneid2go, # geneid/GO associations
            obodag, # Ontologies
            propagate_counts = False,
            alpha = 0.05, # default significance cut-off
            methods = ['fdr_bh']) # we can use bonferroni, fdr_bh...

    # Writing the temorary file (without the Log of odds ratio)
    goea_results_all = goeaobj.run_study(top_nodes)
    goeaobj.wr_tsv(temp_file, goea_results_all)

    results = []
    with open(temp_file, 'r') as ft:
        first_line = ft.readline()
        for line in ft:
            fields = line.strip().split('\t')
            go = fields[0]
            type_func = fields[1]
            name = fields[3]
            ratio_in_study = fields[4]
            ratio_in_pop = fields[5]
            p_uncorrected = float(fields[6])
            depth = fields[7]
            study_count = fields[8]
            p_corrected = float(fields[9])

            # Calculating the Log of odds ratio
            n_in_study, all_in_study = ratio_in_study.split('/')
            n_in_pop, all_in_pop = ratio_in_pop.split('/')
            q = int(n_in_study) # num genes assoc with function in top
            m = int(n_in_pop) # num genes assoc with func in population
            k = int(all_in_study) # num genes in top
            t = int(all_in_pop) # num genes in population
            log_of_odds = calculate_log_of_odds_ratio(q, k, m, t)
            #print('q: {}, m: {}, k: {}, t: {}'.format(q, m, k, t))
            results.append([q, m, log_of_odds, p_uncorrected, p_corrected, go, name, type_func])

    # Writing the output file
    with open(output_file, 'w') as fo:
        fo.write('# num of genes\tnum of total genes\tLog of odds ratio\tP-value\tAdjusted p-value\tGO term ID\tGO term name\tGO type\n')

        for line in sorted(results, key=lambda x: x[2], reverse=True):
            q, m, log_of_odds, p_uncorrected, p_corrected, go, name, type_func = line
            if q <= 0: # Skip the functions that do not have at least 1 top gene associated (the log of odds ratio is -infinite!!)
                continue
            if log_of_odds < 0: # Skip the log of odds ratio that are negative (we will only use functions with log of odds ratio from 0 to inf)
                continue
            if type_func.upper() == 'BP' or type_func.upper() == 'MF':
                new_line = '{}\t{}\t{:.3f}\t{:.2E}\t{:.2E}\t{}\t{}\t{}\n'.format(q, m, log_of_odds, p_uncorrected, p_corrected, go, name, type_func)
                fo.write(new_line)

    # Removing the temporary file
    command = "rm {}".format(temp_file)
    os.system(command)

    return


def calculate_log_of_odds_ratio(q, k, m, t):
    """Calculates the log of the odds ratio"""
    # print("q: {}".format(q))
    # print("k: {}".format(k))
    # print("m: {}".format(m))
    # print("t: {}".format(t))
    odds_ratio = ( (float(q)/float(k)) / (float(m)/float(t)) )
    #odds_ratio = ( (float(q)/(float(k)-float(q))) / (float(m)/((float(t)-float(m)))) )
    if odds_ratio == 0:
        return -float('inf')
    else:
        return float(math.log(odds_ratio, 2))


def parse_guild_scores_file(guild_scores_file):
    """
    Parse the output file (guild scores file) and get the values of interest
    """
    user_entity_ranking = []
    user_entity_to_geneids = {}
    with open(guild_scores_file, 'r') as output_file_fd:
        output_file_fd.readline() # skip header
        for line in output_file_fd:
            fields = line.strip().split("\t")
            if len(fields) == 8:
                [user_entity_id, entry_ids, gene_symbols, is_seed, descriptions, gene_ids, equivalent_entries, score] = fields
            elif len(fields) == 9:
                [user_entity_id, entry_ids, gene_symbols, is_seed, descriptions, gene_ids, equivalent_entries, score, degree] = fields
            user_entity_ranking.append(user_entity_id)
            user_entity_to_geneids.setdefault(user_entity_id, set())
            gene_ids = gene_ids.split(';')
            for geneid in gene_ids:
                user_entity_to_geneids[user_entity_id].add(geneid)
    return user_entity_ranking, user_entity_to_geneids


def get_top_scoring_ueids_and_values(networks_dir, output_dir, percentage, diamond=False):
    """ 
    Parse top scoring nodes and their values.
    """

    edge_scores_file = os.path.join(networks_dir, "edge_scores.txt")
    output_file = os.path.join(output_dir, "guild_scores.txt")

    # Differenciate between top percentage and top enrichment
    if str(percentage)[0:7] == 'enrich_':
        enrich = True
        percentage = int(percentage.split('_')[1])
    else:
        enrich = False
        percentage = float(percentage)

    values = {}
    i = 0
    if diamond:
        # Read edge scores file and 
        # get the number of nodes in the network
        f = open(edge_scores_file)
        nodes = set()
        for line in f:
            node1, interaction, node2 = line.strip().split()
            nodes.add(node1)
            nodes.add(node2)
        network_nodes = len(nodes)
        f.close()
        # Read guild scores file and 
        # get the number of scored nodes
        f = open(output_file)
        f.readline() # skip header
        for line in f:
            i += 1
        scored_nodes = i
        f.close()
        i = 0
        n_total = network_nodes
    else:
        f = open(output_file)
        f.readline() # skip header
        for line in f:
            i += 1
        n_total = i
        f.close()
        i = 0

    # Calculate the last position that we will include in the subnetwork
    if enrich:
        n = percentage
    else:
        n = round(float(percentage) * n_total / 100) # Number of nodes in top percentage
        if diamond:
            if n > scored_nodes:
                n = scored_nodes

    f = open(output_file)
    f.readline() # skip header
    for line in f:
        i += 1
        if i>n:
            break
        fields = line.strip().split("\t")
        if len(fields) == 8:
            [user_entity_id, entry_ids, gene_symbols, is_seed, descriptions, gene_ids, equivalent_entries, score] = line.strip().split("\t")
        elif len(fields) == 9:
            [user_entity_id, entry_ids, gene_symbols, is_seed, descriptions, gene_ids, equivalent_entries, score, degree] = line.strip().split("\t")
        # Somehow the html tags get escaped
        #new_entry_ids = []
        #for entry_id in entry_ids.split("; "):
        #    entry = "<a href='http://www.uniprot.org/uniprot/%s' title='%s'>%s</a>" % (entry_id, "\n".join(descriptions.split("; ")), entry_id)
        #    new_entry_ids.append(entry)
        #if gene_ids == "-":
        #    new_gene_ids = gene_ids
        #else:
        #    new_gene_ids = [ '<a title=\"Equivalent Entries matched by BIANA:\n%s\" href=\"http://www.ncbi.nlm.nih.gov/gene/%s\">%s\</a>' % (", ".join(equivalent_entries.split("; ")), gene_id, gene_id) for gene_id in gene_ids.split("; ") ]
        genes = gene_symbols.split("; ")
        min_idx = 0
        if len(genes) > 1:
            min_len = 100
            for k, gene in enumerate(genes):
                if len(gene) < min_len:
                    min_idx = k
        gene = genes[min_idx]
        val = "non-seed"
        if is_seed == "1":
            val = "seed"
        #values[user_entity_id] = (i, gene_ids, entry_ids, gene, float(score), val)
        values[user_entity_id] = (entry_ids, gene, float(score), i, val)
    f.close()
    return values


def get_network_json(percentage, networks_dir, output_dir, network_instance, drug_info_file, species, diamond=False):
    """
    Create the JSON files to visualize the top-ranking networks in Cytoscape.
    """
    # Create Cytoscape.js network files
    network_json_file = os.path.join(output_dir, "network.json") + percentage
    if not os.path.exists(network_json_file):
        create_network_json_file_for_cytoscapejs(network_json_file, percentage, networks_dir, output_dir, network_instance, drug_info_file, species, diamond=False)
    network_json_new = "".join(open(network_json_file).readlines())
    return None, network_json_new


def create_network_json_file_for_cytoscapejs(network_json_file, percentage, networks_dir, output_dir, network_instance, drug_info_file, species, diamond=False):
    """
    Create a JSON file of the network for Cytoscape.js
    """

    def truncate_text(txt):
        if len(txt) > 30:
            txt = txt[:25]+"..."
        return txt

    f = open(network_json_file, 'w')

    # Parse the GUILD scores file
    output_file = os.path.join(output_dir, "guild_scores.txt")
    top_user_entity_id_to_values = get_top_scoring_ueids_and_values(networks_dir, output_dir, percentage, diamond=diamond)
    all_nodes_to_values = get_top_scoring_ueids_and_values(networks_dir, output_dir, 100, diamond=diamond)

    # Get the score of the user entities
    all_nodes_to_score = {}
    for ueid in all_nodes_to_values:
        score = all_nodes_to_values[ueid][2]
        all_nodes_to_score[ueid] = score

    # Get network as a networkX graph, with combined of score from both sessions
    if not network_instance:
        create_graph_instance(networks_dir, diamond)

    # Get all shortest paths between all nodes
    all_shortest_paths_file = os.path.join(networks_dir, "all_shortest_paths.txt")
    if not os.path.exists(all_shortest_paths_file):
        print('calculating shortest paths...')
        n=1
        network_nodes = network_instance.nodes()
        all_shortest_paths = nx.all_pairs_shortest_path(network_instance)
        #print(nx.__version__)
        all_nodes_to_path = {}
        for node1, node_to_path in all_shortest_paths:
            for node2 in node_to_path:
                path = ','.join(node_to_path[node2])
                all_nodes_to_path.setdefault(node1, {})
                all_nodes_to_path[node1][node2] = path
            #print(n, node1)
            n+=1
        print('...shortest paths calculated.')
        print('writing shortest paths...')
        with open(all_shortest_paths_file, 'w') as out_fd:
            out_fd.write('-\t{}\n'.format('\t'.join(network_nodes)))
            for node1 in network_nodes:
                out_fd.write('{}'.format(node1))
                for node2 in network_nodes:
                    if node1 in all_nodes_to_path and node2 in all_nodes_to_path[node1]:
                        path = all_nodes_to_path[node1][node2]
                        out_fd.write('\t{}'.format(path))
                    else:
                        out_fd.write('\t-')
                out_fd.write('\n')
        print('...shortest paths written.')


    node_str = ""
    edge_str = ""
    subnetwork_file = os.path.join(output_dir, "subnetwork.sif")
    f2 = open(subnetwork_file + "." + percentage, 'w')
    subnetwork_drugs = nx.Graph()

    # Add drug info
    if species == "9606": 
        # Get drug and target info from drug info file
        col_to_idx, drug_to_values, uniprot_to_drugs = get_drug_and_target_info(drug_info_file)
        # Calculate drug scores
        drug_to_score, drug_to_rank = get_drug_to_score(col_to_idx, drug_to_values, top_user_entity_id_to_values)
        drug_to_score = normalize_drug_to_score(drug_to_score)
        drugs = set()
        for ueid, values in top_user_entity_id_to_values.iteritems():
            uniprots = values[0].split("; ")
            for uniprot in uniprots:
                if uniprot in uniprot_to_drugs:
                    for drug in uniprot_to_drugs[uniprot]:
                        if drug in drug_to_score and drug in drug_to_rank:
                            drugs.add((drug, drug_to_score[drug], drug_to_rank[drug]))
                            edge_str += "{\n data: { id: \"%s-%s\", source: \"%s\", target: \"%s\" } \n}" % (ueid, drug, ueid, drug)
                            edge_str += ",\n"
                            f2.write("%s interaction %s\n" % (ueid, drug))
                            subnetwork_drugs.add_edge(ueid, drug)
                        else:
                            #print('Drug without score: {}'.format(drug))
                            pass
        if len(drugs) > 0:
            node_str += ", ".join(map(lambda x: "{\n data: { id: \"%s\", size: %f, xref: \"%s\", label: \"%s\", score: %f, rank: \"%s\", type: \"drug\" } \n}" % (x[0], 100*x[1], x[0], truncate_text(drug_to_values[x[0]][col_to_idx["name"]]), x[1], x[2]), drugs))
            node_str += ",\n"

        # Output drug info file for this subnetwork
        drug_file = os.path.join(output_dir, "drugs.txt")
        f3 = open(drug_file + "." + percentage, 'w')
        f3.write("drugbank id\tname\tgroups\tpubchem id\tdescription\tindication\ttarget(s)\tscore\trank\n")
        drugs = list(drugs)
        drugs.sort(lambda x, y: cmp(y[1], x[1]))
        for drug, score, rank in drugs:
            values = drug_to_values[drug]
            f3.write("%s\t%s\t%s\t%s\n" % (drug, "\t".join(values), score, rank))
        f3.close()

    # Write the subnetwork file
    nodes = set()
    edges = []
    subnetwork = nx.Graph()
    edge_scores_file = os.path.join(networks_dir, "edge_scores.txt")
    for line in open(edge_scores_file):
        id1, w, id2 = line.strip().split()
        if id1 in top_user_entity_id_to_values:
            nodes.add(id1)
            subnetwork.add_node(id1)
            if id2 in top_user_entity_id_to_values:
                edges.append((id1, id2))
                subnetwork.add_edge(id1, id2)
                f2.write("%s interaction %s\n" % (id1, id2))
        if id2 in top_user_entity_id_to_values:
            nodes.add(id2)
            subnetwork.add_node(id2)
    f2.close()


    #### Calculate shortest paths for unconnected nodes ####

    # Define the connected components
    components = list(nx.connected_components(subnetwork))
    if not components or len(components) == 0:
        if len(subnetwork.nodes()) > 0:
            components = [set([node]) for node in subnetwork.nodes()]

    # Get shortest paths
    all_nodes_to_path = get_shortest_paths_nodes(all_shortest_paths_file, subnetwork.nodes())

    # Compare the longest connected component with the rest of components
    # and find the shortest paths between them
    connectors = set()
    if len(components) > 1:
        lcc, other_components = decide_longest_connected_component(components, all_nodes_to_score)
        for cc in other_components:
            sps, all_nodes_to_path = find_shortest_paths_two_components(cc, lcc, network_instance, node_to_path=all_nodes_to_path, decide_by_score=True, node_to_score=all_nodes_to_score)
            for sp in sps:
                new_nodes = sp[1:-1]
                for x in xrange(len(sp)):
                    if x == len(sp)-1:
                        break
                    node1 = str(sp[x])
                    node2 = str(sp[x+1])
                    lcc.add(node1) # Add the new nodes in the lcc
                    lcc.add(node2)
                    if not sorted([node1, node2]) in edges:
                        edges.append(sorted([node1, node2]))
                        subnetwork.add_edge(node1, node2)
                    if node1 in new_nodes:
                        connectors.add(node1)
                        # Add edges with other nodes of subnetwork
                        neighbors1 = network_instance.neighbors(node1)
                        for ngbr in neighbors1:
                            if ngbr in subnetwork.nodes():
                                if not sorted([node1, ngbr]) in edges:
                                    edges.append(sorted([node1, ngbr]))
                                    subnetwork.add_edge(node1, ngbr)
                    if node2 in new_nodes:
                        connectors.add(node2)
                        # Add edges with other nodes of subnetwork
                        neighbors2 = network_instance.neighbors(node2)
                        for ngbr in neighbors2:
                            if ngbr in subnetwork.nodes():
                                if not sorted([node2, ngbr]) in edges:
                                    edges.append(sorted([node2, ngbr]))
                                    subnetwork.add_edge(node2, ngbr)
            # Add the nodes of the component in the lcc
            for node in cc:
                lcc.add(node)

    # Get the composition of the subnetwork of proteins and the subnetwork of drugs
    final_subnetwork = nx.compose(subnetwork, subnetwork_drugs)

    # Write the JSON file
    rows = []
    for x in top_user_entity_id_to_values:
        size = top_user_entity_id_to_values[x][2]*100 # Score multiplied by 100
        xref = top_user_entity_id_to_values[x][0]
        label = top_user_entity_id_to_values[x][1]
        score = top_user_entity_id_to_values[x][2]
        rank = top_user_entity_id_to_values[x][3]
        type_node = top_user_entity_id_to_values[x][4]
        if label == '-':
            label = xref
        row = "{\n data: { id: \"%s\", size: %f, xref: \"%s\", label: \"%s\", score: %f, rank: \"%s\", type: \"%s\" } \n}" % (x, size, xref, label, score, rank, type_node)
        rows.append(row)
        if type_node=='seed':
            final_subnetwork.add_node(x, color='#B0FF00', shape='8', size=size)
        else:
            final_subnetwork.add_node(x, color='yellow', shape='o', size=size)
    for x in connectors:
        size = all_nodes_to_values[x][2]*100 # Score multiplied by 100
        xref = all_nodes_to_values[x][0]
        label = all_nodes_to_values[x][1]
        score = all_nodes_to_values[x][2]
        rank = all_nodes_to_values[x][3]
        type_node = 'linker'
        if label == '-':
            label = xref
        row = "{\n data: { id: \"%s\", size: %f, xref: \"%s\", label: \"%s\", score: %f, rank: \"%s\", type: \"%s\" } \n}" % (x, size, xref, label, score, rank, type_node)
        rows.append(row)
        final_subnetwork.add_node(x, color='gray', shape='o', size=size)
    node_str += ", ".join(rows)
    edge_str += ", ".join(map(lambda x: "{\n data: { id: \"%s-%s\", source: \"%s\", target: \"%s\" } \n}" % (x[0],x[1],x[0],x[1]), edges))

    f.write("[\n%s,\n%s\n]"%(node_str, edge_str))
    f.close()

    if species == "9606": 
        for node in final_subnetwork.nodes():
            if node in drug_to_values:
                score = drug_to_score[node]
                final_subnetwork.add_node(node, color='#5093FF', shape='d', size=score*100)

    import matplotlib
    matplotlib.use('Agg')
    from matplotlib import pyplot as plt

    #Get all distinct node classes according to the node shape attribute
    nodePos = nx.fruchterman_reingold_layout(final_subnetwork)
    nodeShapes = set((aShape[1]["shape"] for aShape in final_subnetwork.nodes(data = True)))
    for aShape in nodeShapes:
        #...filter and draw the subset of nodes with the same symbol in the positions that are now known through the use of the layout.
        nodelist=[sNode[0] for sNode in filter(lambda x: x[1]["shape"]==aShape,final_subnetwork.nodes(data = True))]
        node_color=[sNode[1]["color"] for node in final_subnetwork.nodes(data = True) if node[1]["shape"]==aShape]
        node_size=[sNode[1]["size"] for node in final_subnetwork.nodes(data = True) if node[1]["shape"]==aShape]
        nodes = nx.draw_networkx_nodes(final_subnetwork,nodePos,node_shape=aShape, nodelist=nodelist, node_color=node_color, node_size=node_size, linewidths=0.5)
        nodes.set_edgecolor('black')
    #Finally, draw the edges between the nodes
    nx.draw_networkx_edges(final_subnetwork,nodePos, ax=None)
    #nx.draw(final_subnetwork, node_color=node_color, node_size=node_size, node_shape=node_shape, linewidths=1.0)

    network_plot = os.path.join(output_dir, 'network_plot.png')
    plt.axis('off')
    plt.savefig(os.path.join(network_plot + percentage), format="PNG")

    return 


def get_drug_and_target_info(drug_info_file):
    """
    Get the drug and target info from the drug_info.txt file.
    """
    uniprot_to_drugs = {}
    col_to_idx, drug_to_values = get_drug_info(drug_info_file)
    for drug, values in drug_to_values.iteritems():
        #print drug, values
        uniprots = values[col_to_idx["targets"]].split("; ")
        for uniprot in uniprots:
            uniprot_to_drugs.setdefault(uniprot, set()).add(drug)
    return col_to_idx, drug_to_values, uniprot_to_drugs



def get_drug_info(drug_info_file):
    drug_to_values = {}
    f = open(drug_info_file)
    header = f.readline().strip().split("\t")
    col_to_idx = dict((k, i) for i, k in enumerate(header[1:])) 
    for line in f:
        words = line.strip("\n").split("\t")
        drug_to_values[words[0]] = words[1:]
    return col_to_idx, drug_to_values


def calculate_drug_score_from_targets(values):
    val = 0.0
    for value in values:
        val += value * value
    return math.sqrt(val)


def get_drug_to_score(col_to_idx, drug_to_values, user_entity_id_to_values):
    """
    Calculate the score and ranking for the drugs.
    """
    drug_to_score = {}
    uniprot_to_score = {}
    # Get the scores for the proteins
    for ueid, values in user_entity_id_to_values.iteritems():
        for uniprot in values[0].split(";"):
            uniprot_to_score[uniprot] = values[2]

    # Assign a score to all the drugs associated to proteins
    for drug, values in drug_to_values.iteritems():
        scores = []
        targets = values[col_to_idx["targets"]].split("; ")
        for uniprot in targets:
            if uniprot in uniprot_to_score:
                scores.append(uniprot_to_score[uniprot])
        if len(scores) == 0:
            val = 0
            continue
        else:
            val = calculate_drug_score_from_targets(scores)
        drug_to_score[drug] = val #/ len(scores) # len(targets)

    # Calculate the ranking for all the drugs
    drug_to_rank = calculate_drug_to_rank(drug_to_score)

    return drug_to_score, drug_to_rank


def normalize_drug_to_score(drug_to_score):
    """
    Normalize the scores of the drugs from 1 to 0.
    """
    new_drug_to_score = {}
    i = 0
    highest_score = 1
    for drug,score in sorted(drug_to_score.items(), key=lambda x:x[1], reverse=True):
        if i == 0:
            highest_score = score
        new_score = float(score) / float(highest_score)
        new_drug_to_score[drug] = new_score
        i+=1
    return new_drug_to_score


def calculate_drug_to_rank(drug_to_score):
    """
    Calculate the ranks of the drugs
    """
    drug_to_rank = {}
    rank = 1
    initial_rank = 1
    previous_drug = None
    stored_drugs = []
    for drug,score in sorted(drug_to_score.items(), key=lambda x:x[1], reverse=True):
        drug_to_rank[drug] = rank
        if previous_drug and drug_to_score[previous_drug] == score:
            rank_output = '{}-{}'.format(initial_rank, rank)
            stored_drugs.append(drug)
            for stored_drug in stored_drugs:
                drug_to_rank[stored_drug] = rank_output
        else:
            stored_drugs = []
            initial_rank = rank
            stored_drugs.append(drug)
        rank+=1
        previous_drug = drug
    return drug_to_rank


def get_shortest_paths_nodes(all_shortest_paths_file, nodes):
    all_nodes_to_path = {}
    with open(all_shortest_paths_file, 'r') as inp_fd:
        nodes_network = inp_fd.readline().strip().split('\t')[1:]
        for line in inp_fd:
            fields = line.strip().split('\t')
            curr_node = fields[0]
            if curr_node in nodes:
                paths = fields[1:]
                dictionary = dict(zip(nodes_network, paths))
                all_nodes_to_path[curr_node] = dictionary
            if len(nodes) == len(all_nodes_to_path):
                break
    return all_nodes_to_path


def find_shortest_paths_two_components(cc, lcc, network, node_to_path=None, decide_by_score=True, node_to_score=None):
    """
    Find the shortest paths between two components.
    If decide_by_score is true, between all the shortest paths, we look for the one containing the node of maximum score
    """
    shortest_paths = []
    min_len = None
    for n1 in cc:
        for n2 in lcc:
            if node_to_path and n1 in node_to_path and n2 in node_to_path[n1]:
                sp = node_to_path[n1][n2].split(',')
            else:
                sp = list(nx.shortest_path(network, source=n1, target=n2))
                node_to_path.setdefault(n1, {})
                node_to_path[n1][n2] = sp
            path_length = len(sp)
            if not min_len:
                min_len = path_length
                shortest_paths.append(sp)
            elif path_length < min_len:
                min_len = path_length
                shortest_paths = []
                shortest_paths.append(sp)
            elif path_length == min_len:
                shortest_paths.append(sp)
    if decide_by_score:
        best_max_score = None
        best_mean_score = None
        final_shortest_paths=[]
        for sp in shortest_paths:
            new_nodes = sp[1:-1]
            if len(new_nodes)>0:
                scores = [node_to_score[node] for node in new_nodes if node in node_to_score]
                if len(scores) != len(new_nodes):
                    for node in new_nodes:
                        # In case that there are nodes in the path that are not
                        # in the node_info_file and therefore not in the node_to_score
                        # dictionary, we skip them.
                        if node not in node_to_score:
                            print(node)
                            pass
                    continue
                max_score = max(scores)
                mean_score = np.mean(scores)
                if not best_max_score:
                    best_max_score=max_score
                    best_mean_score=mean_score
                    final_shortest_paths.append(sp)
                elif max_score > best_max_score:
                    best_max_score=max_score
                    best_mean_score=mean_score
                    final_shortest_paths=[]
                    final_shortest_paths.append(sp)
                elif max_score == best_max_score:
                    if mean_score > best_mean_score:
                        best_mean_score=mean_score
                        final_shortest_paths=[]
                        final_shortest_paths.append(sp)
                    elif mean_score == best_mean_score:
                        final_shortest_paths.append(sp)
            else:
                final_shortest_paths.append(sp)
        return final_shortest_paths, node_to_path
    else:
        return shortest_paths, node_to_path


def decide_longest_connected_component(components, node_to_score):
    """
    Check the length of the components. If there are more than one lcc,
    decide which is the best one to consider as lcc.
    We get the component with highest mean score.
    """
    max_len=0
    lccs=[]
    for component in components:
        if len(component)>max_len:
            max_len=len(component)
            lccs.append(component)

    if len(lccs)==1:
        other_components = []
        for component in components:
            if component != lccs[0]:
                other_components.append(component)
        return set(lccs[0]), other_components
    else:
        best_mean_score=0
        lcc = []
        for component in lccs:
            mean_score = np.mean([node_to_score[node] for node in component])
            if mean_score>best_mean_score:
                best_mean_score=mean_score
                lcc = component
        other_components = []
        for component in components:
            if component != lcc:
                other_components.append(component)
        return set(lcc), other_components




#################
#### CLASSES ####
#################

class BiologicalValidation(object):
    """ 
    Class defining a biological validation object 
    """

    def __init__(self, gene_to_go_file, go_to_genes_file):
        """ 
        @param:    gene_to_go_file
        @pdef:     Path to the file containing the dictionary geneID : set([GOIDs])
        @ptype:    {String}

        @param:    go_to_genes_file
        @pdef:     Path to the file containing the dictionary GOID : set([geneID])
        @ptype:    {String}

        """

        self.gene_to_go_file = gene_to_go_file
        self.go_to_genes_file = go_to_genes_file

        #self.gene_to_go = cPickle.load(open(self.gene_to_go_file))
        #self.go_to_genes = cPickle.load(open(self.go_to_genes_file))


    ###########
    # METHODS #
    ###########


    def parse_seed_file(self, seed_file, seed_info=False):
        """
        Parse a file containing the seeds separated by new lines.
        """
        seeds = set()
        with open(seed_file, 'r') as seed_file_fd:
            if seed_info:
                first_line = seed_file_fd.readline() #Skip first line
            for line in seed_file_fd:
                fields = line.strip().split('\t')
                seed = fields[0]
                if seed_info:
                    geneids = fields[4].split(';')
                    for geneid in geneids:
                        seeds.add(geneid)
        return seeds


    def parse_background_genes_file(self, background_genes_file, node_info=False):
        """
        Parse a file containing the seeds separated by new lines.
        """
        background_genes = set()
        with open(background_genes_file, 'r') as background_genes_fd:
            if node_info:
                first_line = background_genes_fd.readline() #Skip first line
            for line in background_genes_fd:
                fields = line.strip().split('\t')
                background_gene = fields[0]
                if node_info:
                    geneids = fields[4].split(';')
                    for geneid in geneids:
                        background_genes.add(geneid)
        return background_genes


    def parse_guild_scores_file(self, guild_scores_file):
        """
        Parse the output file (guild scores file) and get the values of interest
        """
        user_entity_ranking = []
        user_entity_to_geneids = {}
        with open(guild_scores_file, 'r') as output_file_fd:
            output_file_fd.readline() # skip header
            for line in output_file_fd:
                fields = line.strip().split("\t")
                if len(fields) == 8:
                    [user_entity_id, entry_ids, gene_symbols, is_seed, descriptions, gene_ids, equivalent_entries, score] = fields
                elif len(fields) == 9:
                    [user_entity_id, entry_ids, gene_symbols, is_seed, descriptions, gene_ids, equivalent_entries, score, degree] = fields
                else:
                    print(line)
                    sys.exit(10)
                user_entity_ranking.append(user_entity_id)
                user_entity_to_geneids.setdefault(user_entity_id, set())
                gene_ids = gene_ids.split(';')
                for geneid in gene_ids:
                    user_entity_to_geneids[user_entity_id].add(geneid)
        return user_entity_ranking, user_entity_to_geneids


    def calculate_functional_enrichment(self, test_geneids, background_geneids, enrichment_file, correction = 'bonferroni', max_num_of_genes_per_GO = None):
        """
        Calculate the functions that are enriched in a set of 'test' genes
        with respect to a set of 'background' genes.
        correction = Method to correct the p-value: Bonferroni = 'bonferroni', Benjamini/Hochberg  (non-negative) = 'fdr_bh'
        max_num_of_genes_per_GO = The maximum number of genes that a GO has to have to be considered.
        """

        #------------------------------------------------------------#
        # Get the GOs associated to the GeneIDs of the tested genes  #
        #------------------------------------------------------------#

        test_geneids_with_GOs = set()
        GOs = set()

        # Get the test genes that have GOs
        # Get the GOs of the test genes
        for geneid in test_geneids:
            if geneid in self.gene_to_go:
                if max_num_of_genes_per_GO:
                    # Only accept genes that have GOs not exceeding 
                    # the maximum number of genes per GO
                    GO_accepted = False
                    for go_id in self.gene_to_go[geneid]:
                        if go_id in self.go_to_genes:
                            if len(self.go_to_genes[go_id]) <= max_num_of_genes_per_GO:
                                GO_accepted = True
                                GOs.add(go_id) # Add all the GOs from test genes in one set
                    if GO_accepted:
                        test_geneids_with_GOs.add(geneid) # Add the test genes containing GOs in one set
                else:
                    test_geneids_with_GOs.add(geneid) # Add the test genes containing GOs in one set
                    for go_id in self.gene_to_go[geneid]:
                        GOs.add(go_id) # Add all the GOs from test genes in one set

        # Get the background genes that have GOs
        background_geneids_with_GOs = set()
        for geneid in background_geneids:
            if geneid in self.gene_to_go:
                if max_num_of_genes_per_GO:
                    # Only accept genes that have GOs not exceeding 
                    # the maximum number of genes per GO
                    GO_accepted = False
                    for go_id in self.gene_to_go[geneid]:
                        if go_id in self.go_to_genes:
                            if len(self.go_to_genes[go_id]) <= max_num_of_genes_per_GO:
                                GO_accepted = True
                                break
                    if GO_accepted:
                        background_geneids_with_GOs.add(geneid) # Add the background genes containing GOs in one set
                else:
                    background_geneids_with_GOs.add(geneid) # Add the background genes containing GOs in one set

        non_test = background_geneids_with_GOs - test_geneids_with_GOs # Obtain all the non test genes


        #---------------------------------------------------------------------------#
        # Calculate if the GOs associated to test genes are significantly enriched  #
        #---------------------------------------------------------------------------#

        GOs_enriched = set()
        genes_in_GOs_enriched = set()
        GOs_analyzed = []
        pvalues = []
        oddsratios = []

        for go_id in GOs:

            test_annotated = set()
            non_test_annotated = set()

            if go_id in self.go_to_genes:
                for gene in self.go_to_genes[go_id]:
                    if gene in test_geneids_with_GOs:
                        test_annotated.add(gene) # Obtain the test genes annotated
                    else:
                        non_test_annotated.add(gene) # Obtain the non-test genes annotated
            test_non_annotated = test_geneids_with_GOs - test_annotated # Obtain the test genes non annotated
            non_test_non_annotated = non_test - non_test_annotated
            contingency_table = [[len(test_annotated), len(test_non_annotated)], [len(non_test_annotated), len(non_test_non_annotated)]]
            oddsratio, pvalue = fisher_exact(contingency_table)
            GOs_analyzed.append(go_id)
            pvalues.append(pvalue)
            oddsratios.append(oddsratio)

        # Calculate the Bonferroni correction
        hypotheses, pvalues_adjusted, sidak, bonf = multipletests(pvalues, method='bonferroni')

        # Write the results
        with open(enrichment_file, 'w') as enrichment_fd:

            enrichment_fd.write('#GOid\tOddsRatio\tP-value\tAdjP-value\n')

            for x in xrange(len(GOs_analyzed)):
                enrichment_fd.write('{}\t{:.3f}\t{:.3f}\t{:.3f}\n'.format(GOs_analyzed[x], float(oddsratios[x]), float(pvalues[x]), float(pvalues_adjusted[x])))
                if float(pvalues_adjusted[x]) < 0.5:
                    GOs_enriched.add(GOs_analyzed[x])
                    for gene in self.go_to_genes[GOs_analyzed[x]]:
                        genes_in_GOs_enriched.add(gene)

        print('Calculation of the GOs associated to test genes which are significantly enriched done!')

        return GOs_enriched, genes_in_GOs_enriched

    def recover_enrichment_of_seeds(self, enrichment_file):
        """
        Recover the GOs enriched from the seeds and the genes associated
        to these GOs.
        """
        GOs_enriched = set()
        genes_in_GOs_enriched = set()

        with open(enrichment_file, 'r') as enrichment_fd:
            for line in enrichment_fd:
                if line[0] == '#':
                    continue # Skip comments
                go_id, oddsratio, pvalue, p_adjusted = line.strip().split('\t')
                if float(p_adjusted) < 0.5:
                    GOs_enriched.add(go_id)
                    for gene in self.go_to_genes[go_id]:
                        genes_in_GOs_enriched.add(gene)

        return GOs_enriched, genes_in_GOs_enriched


    def get_enriched_genes_in_ranking(self, output_dir, user_entity_ranking, user_entity_to_geneids, seeds, genes_in_GOs_enriched):
        """
        Get the genes in the ranking that are associated to significantly
        enriched functions of the seeds.
        """

        #---------------------------------------------------------------------------#
        # Get the genes of the ranking that are part of significantly enriched GOs  #
        #---------------------------------------------------------------------------#

        enriched_positions = set()
        non_enriched_positions = set()
        seed_ranking = set()
        rank = 0
        seed_ranking_file = os.path.join(output_dir, 'seed_ranking.txt')

        with open(seed_ranking_file, 'w') as seed_ranking_fd:
            seed_ranking_fd.write('#Position\tSeedUserEntity\n')
            for user_entity in user_entity_ranking:
                rank+=1
                tp = False
                if user_entity in user_entity_to_geneids:
                    for geneid in user_entity_to_geneids[user_entity]:
                        if geneid in genes_in_GOs_enriched:
                            #print('Gene {} enriched!'.format(geneid))
                            tp = True
                            break
                if tp:
                    enriched_positions.add(rank)
                else:
                    non_enriched_positions.add(rank)
                if user_entity in seeds:
                    seed_ranking.add(rank)
                    seed_ranking_fd.write('{}\t{}\n'.format(rank, user_entity))

        return enriched_positions, non_enriched_positions, seed_ranking


    def calculate_initial_position_of_sliding_window(self, num_seeds):
        """
        Calculate the initial position where the sliding position
        will start its calculations.
        """
        # Calculate the initial position of the sliding window:
        # We know that the left part of the first interval in the sliding window is:
        # i - num_seeds / 2 = 1
        # So, if we want to know i (the initial position) we have to calculate:
        # i = 1 + num_seeds / 2
        initial_position = 1 + int( float(num_seeds) / float(2) )
        print('Initial position: {}'.format(initial_position))
        return initial_position


    def calculate_final_position_of_sliding_window(self, num_seeds, maximum_rank=500):
        """
        Calculate the initial position where the sliding position
        will start its calculations.
        - maximum_rank = The maximum number of positions that we want to include
                         in the calculations
        """
        # Calculate the final position of the sliding window:
        # We know that the right part of the last interval will be:
        # i + num_seeds / 2 = last_rank
        # So the final position will be:
        # i = last_rank - num_seeds / 2
        last_rank = maximum_rank
        final_position = last_rank - int( float(num_seeds) / float(2) )
        print('Final position: {}'.format(final_position))
        return final_position

    def calculate_enrichment_in_all_sliding_window_positions(self, output_dir, num_seeds, enriched_positions, non_enriched_positions, seed_ranking, maximum_rank=500):
        """
        Calculate the enrichment in all the positions of the sliding window.
        """
        #----------------------------------------------------------------------#
        # Calculate the enrichment in all the positions of the sliding window  #
        #----------------------------------------------------------------------#

        all_positions = range(self.calculate_initial_position_of_sliding_window(num_seeds), self.calculate_final_position_of_sliding_window(num_seeds, maximum_rank)+1)

        positions = [] # All the positions in a list
        pvalues = [] # Pvalues obtained
        # When below_cutoff = 0, it has not found any p-value under 0.05
        # When below_cutoff = 1, it is finding for the first time values under 0.05
        # When below_cutoff = 2, it has gone over 0.05
        below_cutoff = 0
        cutoff = 0
        cutoff_right_interval = 0

        enrichment_positions_file = os.path.join(output_dir, 'enrichment_positions.txt')

        with open(enrichment_positions_file, 'w') as enrichment_positions_fd:
            enrichment_positions_fd.write('#Position\tOddsRatio\tP-value\n')
            for position in all_positions:
                left_interval = position - int( float(num_seeds) / float(2) ) # Define left interval position
                right_interval = position + int( float(num_seeds) / float(2) ) # Define right interval position
                window = set(range(left_interval, right_interval+1)) # Define window positions
                non_window = set(all_positions) - window # Define non-window positions
                enriched_positions_window = enriched_positions & window
                non_enriched_positions_window = non_enriched_positions & window
                enriched_positions_non_window = enriched_positions & non_window
                non_enriched_positions_non_window = non_enriched_positions & non_window
                # Calculate if the enrichment is significant
                contingency_table = [[len(enriched_positions_window), len(non_enriched_positions_window)], [len(enriched_positions_non_window), len(non_enriched_positions_non_window)]]
                oddsratio, pvalue = fisher_exact(contingency_table)
                enrichment_positions_fd.write('{}\t{}\t{}\n'.format(position, oddsratio, pvalue))
                positions.append(position)
                pvalues.append(pvalue)
                if position in seed_ranking:
                    continue
                if below_cutoff == 0 and pvalue <= 0.05:
                    below_cutoff = 1
                if below_cutoff == 1 and pvalue <= 0.05:
                    cutoff = position # Get the last position where the cut-off is below 0.05
                    cutoff_right_interval = right_interval # Get the last node which is included in the interval
                elif below_cutoff == 1 and pvalue > 0.05:
                    below_cutoff = 2 # When pvalue is above 0.05, we stop getting the cutoff

        return positions, pvalues, cutoff, cutoff_right_interval

    def recover_enrichment_in_sliding_positions(self, enrichment_positions_file, seed_ranking_file, num_seeds):
        """
        Recover the calculation of the enrichment in all the sliding positions
        from the enrichment_positions.txt file.
        """
        positions = [] # All the positions in a list
        pvalues = [] # Pvalues obtained
        seed_ranking = set() # Ranking positions of seeds
        below_cutoff = 0
        cutoff = 0
        cutoff_right_interval = 0

        with open(seed_ranking_file, 'r') as seed_ranking_fd:
            for line in seed_ranking_fd:
                if line[0] == '#':
                    continue # Skip comments
                position, seed = line.strip().split('\t')
                seed_ranking.add(int(position))


        with open(enrichment_positions_file, 'r') as enrichment_positions_fd:
            for line in enrichment_positions_fd:
                if line[0] == '#':
                    continue # Skip comments
                position, oddsratio, pvalue = line.strip().split('\t')
                position = int(position)
                pvalue = float(pvalue)
                positions.append(position)
                pvalues.append(pvalue)
                right_interval = position + int( float(num_seeds) / float(2) ) # Define right interval position
                if below_cutoff == 0 and pvalue <= 0.05:
                    below_cutoff = 1
                if below_cutoff == 1 and pvalue <= 0.05:
                    cutoff = position # Get the last position where the cut-off is below 0.05
                    cutoff_right_interval = right_interval # Get the last node which is included in the interval
                elif below_cutoff == 1 and pvalue > 0.05:
                    below_cutoff = 2 # When pvalue is above 0.05, we stop getting the cutoff

        return positions, pvalues, cutoff, cutoff_right_interval, seed_ranking

    def output_plotly_information(self, positions, pvalues, seed_ranking, maximum_rank, cutoff_right_interval, output_dir):
        """
        Save the data introduced in Plotly as a file.
        """
        position_values = []
        enrichment_values = []
        seed_values = []
        validated_nodes = []
        non_validated_nodes = []
        rank_range = [1, maximum_rank]
        pos_index = 0
        for x in xrange(maximum_rank):
            rank = x+1
            if rank in positions:
                position_values.append(positions[pos_index])
                enrichment_values.append(pvalues[pos_index])
                if rank in seed_ranking:
                    seed_values.append(positions[pos_index])
                else:
                    seed_values.append('')
                if rank <= cutoff_right_interval:
                    validated_nodes.append(positions[pos_index])
                    non_validated_nodes.append('')
                elif rank > cutoff_right_interval:
                    non_validated_nodes.append(positions[pos_index])
                    validated_nodes.append('')
                pos_index+=1

        # Define the output file
        javascript_output_file = os.path.join(output_dir, 'plotly_structures.txt')
        # Write the output file
        with open(javascript_output_file, 'w') as javascript_output_fd:
            javascript_output_fd.write('#position_values\n')
            javascript_output_fd.write('{}\n'.format(','.join([str(x) for x in position_values])))
            javascript_output_fd.write('#enrichment_values\n')
            javascript_output_fd.write('{}\n'.format(','.join([str(x) for x in enrichment_values])))
            javascript_output_fd.write('#seed_values\n')
            javascript_output_fd.write('{}\n'.format(','.join([str(x) for x in seed_values])))
            javascript_output_fd.write('#validated_nodes\n')
            javascript_output_fd.write('{}\n'.format(','.join([str(x) for x in validated_nodes])))
            javascript_output_fd.write('#non_validated_nodes\n')
            javascript_output_fd.write('{}\n'.format(','.join([str(x) for x in non_validated_nodes])))
            javascript_output_fd.write('#cut-off right interval\n')
            javascript_output_fd.write('{}\n'.format(cutoff_right_interval))
            javascript_output_fd.write('#range\n')
            javascript_output_fd.write('{},{}\n'.format(rank_range[0], rank_range[1]))
        return position_values, enrichment_values, seed_values, validated_nodes, non_validated_nodes, rank_range


    def output_javascript_dictionaries(self, positions, pvalues, seed_ranking, maximum_rank, cutoff, cutoff_right_interval, output_dir):
        """
        Output the javascript dictionaries used by GUILDify
        to plot the final graphic.
        """
        #------------------------------------------#
        # Obtain the output javascript dictionary  #
        #------------------------------------------#        

        # Obtain a javascript dictionary of positions and p-values for the graphic
        #plot_values = '['      # It was used when ScatterChart, but now we do not use it with LineChart
        position_values = '['   # Data for the x axis, containing the values of the positions
        enrichment_values = '[' # Data for the y axis, containing the p-values of the enrichment
        line_values = '['       # Data for the y axis of the p-value=0.05 line
        point_sizes = '['       # Data of the sizes of the points
        point_colors = '['      # Data of the colors of the points
        pos_index = 0
        for x in xrange(maximum_rank):
            rank = x+1
            line_values += '0.05,' 
            if rank not in positions:
                position_values+='%s,'%(rank)
                enrichment_values+=','
                point_sizes += '0,' 
                point_colors+=','
            else:
                #plot_values+='{x:%s, y:%s},'%(positions[pos_index], pvalues[pos_index])
                position_values+='%s,'%(positions[pos_index])
                enrichment_values+='%s,'%(pvalues[pos_index])
                pos_index+=1
                # if rank in seed_ranking and rank != cutoff:
                #     point_colors+='%s,'%('window.chartColors.orange')
                #     point_sizes += '3,' 
                # elif rank == cutoff or rank == cutoff_right_interval:
                #     point_colors+='%s,'%('\'red\'')
                #     point_sizes += '6,' 
                # else:
                #     point_colors+='%s,'%('window.chartColors.green')
                #     point_sizes += '3,' 
                if rank in seed_ranking and rank <= cutoff_right_interval:
                    point_colors+='%s,'%('\'red\'')
                    point_sizes += '5,' 
                elif rank in seed_ranking and rank > cutoff_right_interval:
                    point_colors+='%s,'%('window.chartColors.orange')
                    point_sizes += '3,' 
                elif rank not in seed_ranking and rank <= cutoff_right_interval:
                    point_colors+='%s,'%('\'#2471a3\'')
                    point_sizes += '5,' 
                elif rank not in seed_ranking and rank > cutoff_right_interval:
                    point_colors+='%s,'%('window.chartColors.green')
                    point_sizes += '3,' 
        #plot_values += ']'
        position_values += ']'
        enrichment_values += ']'
        line_values += ']'
        point_sizes += ']'
        point_colors += ']'

        # Define the output file
        javascript_output_file = os.path.join(output_dir, 'javascript_structures.txt')
        # Write the output file
        with open(javascript_output_file, 'w') as javascript_output_fd:
            javascript_output_fd.write('#position_values\n')
            javascript_output_fd.write('{}\n'.format(position_values))
            javascript_output_fd.write('#enrichment_values\n')
            javascript_output_fd.write('{}\n'.format(enrichment_values))
            javascript_output_fd.write('#line_values\n')
            javascript_output_fd.write('{}\n'.format(line_values))
            javascript_output_fd.write('#point_sizes\n')
            javascript_output_fd.write('{}\n'.format(point_sizes))
            javascript_output_fd.write('#point_colors\n')
            javascript_output_fd.write('{}\n'.format(point_colors))
            javascript_output_fd.write('#cutoff_right_interval\n')
            javascript_output_fd.write('{}\n'.format(cutoff_right_interval))

        return position_values, enrichment_values, line_values, point_sizes, point_colors

    def parse_functional_profile_goatools(self, functional_profile):
        """
        Parses a functional profile obtained from an functional enrichment analysis of a list of nodes.
        Returns a dictionary with:
        - Key = GO term id
        - Value = Set --> (Log of odds ratio, Adjusted p-value)
        """
        go_id_to_values = {}
        with open(functional_profile, 'r') as f:
            for line in f:
                if line[0] == '#':
                    continue
                num_genes, total_genes, log_odds_ratio, pval, adj_pval, go_id, go_name, go_type = line.strip().split('\t')
                go_id_to_values[go_id] = (float(log_odds_ratio), float(adj_pval))
        return go_id_to_values



if  __name__ == "__main__":
    main()