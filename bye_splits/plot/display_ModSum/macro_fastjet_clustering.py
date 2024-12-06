import fastjet

def perform_clustering():
    # Define some dummy input data for a single event
    pseudojet_data = [
        fastjet.PseudoJet(1.1, 1.2, 1.3, 11.4),
        fastjet.PseudoJet(2.1, 2.2, 2.3, 12.4),
        fastjet.PseudoJet(3.1, 3.2, 3.3, 13.4),
    ]
    
    # Create a JetDefinition using the Anti-kt algorithm with radius parameter 0.6
    jet_definition = fastjet.JetDefinition(fastjet.antikt_algorithm, 0.6)
    
    # Perform the clustering using the ClusterSequence class
    cluster_sequence = fastjet.ClusterSequence(pseudojet_data, jet_definition)
    
    # Extract the inclusive jets
    inclusive_jets = cluster_sequence.inclusive_jets()
    
    # Print the clustered jets information
    print("Clustered jets using Anti-kt algorithm:")
    for jet in inclusive_jets:
        print(f"px: {jet.px():.2f}, py: {jet.py():.2f}, pz: {jet.pz():.2f}, E: {jet.E():.2f}, phi: {jet.phi():.2f}, eta: {jet.eta():.2f} ")

# Run the function
perform_clustering()
