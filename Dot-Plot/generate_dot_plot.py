import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import ListedColormap
from util import get_diagonals

def plot_protein_sequence_alignment(protein1, protein2):
    """
    Generate a dot plot for pairwise protein sequence alignment with identified diagonals.
    
    This function creates a dot plot that visualizes the pairwise alignment of two protein sequences.
    Matches between the sequences are indicated by red squares, and diagonals (both main and reverse)
    are marked with black and green lines, respectively.

    Args:
        protein1 (str): The first protein sequence (amino acids).
        protein2 (str): The second protein sequence (amino acids).
    
    Returns:
        None: Displays the dot plot using matplotlib.
    """
    # Initialize an empty alignment matrix for matches
    alignment_matrix = np.zeros((len(protein1), len(protein2), 3))

    # Populate the alignment matrix: 1 indicates a match between amino acids
    for i in range(len(protein1)):
        for j in range(len(protein2)):
            if protein1[i] == protein2[j]:
                alignment_matrix[i, j, 0] = 1

    # Extract main and reverse diagonals using get_diagonals
    main_diagonals, reverse_diagonals = get_diagonals(alignment_matrix)

    # Extract only the match (0 or 1) part for plotting
    plot_matrix = alignment_matrix[:, :, 0]

    # Define a custom colormap (red for matches, white for mismatches)
    colors = [(1, 1, 1), (1, 0, 0)]  # Red for match, white for mismatch
    cmap = ListedColormap(colors)

    # Create the dot plot using pcolor
    plt.pcolor(plot_matrix[::-1], cmap=cmap, edgecolors='k', linewidths=1)

    # Plot main and reverse diagonals
    for diagonal in main_diagonals:
        plt.plot(diagonal[0], diagonal[1], ".-", color='black', linewidth=3, markersize=10)
    for diagonal in reverse_diagonals:
        plt.plot(diagonal[0], diagonal[1], ".-", color='green', linewidth=3, markersize=10)

    plt.xticks(np.arange(len(protein2)) + 0.5, list(protein2))
    plt.yticks(np.arange(len(protein1)) + 0.5, list(protein1[::-1]))

    plt.xlabel("Protein Sequence 2")
    plt.ylabel("Protein Sequence 1")
    
    # Move the x-axis to the top
    plt.gca().xaxis.tick_top()
    plt.gca().xaxis.set_label_position('top')

    # Create the legends for matches and diagonals
    match_legend = plt.Rectangle((0, 0), 1, 1, color='red', label='Amino Acid Match')
    main_diagonal_legend = plt.Line2D([0], [0], marker='.', linestyle="-", color='black', markersize=10, label='Main Diagonal')
    reverse_diagonal_legend = plt.Line2D([0], [0], marker='.', linestyle="-", color='green', markersize=10, label='Reverse Diagonal')
    
    plt.legend(handles=[match_legend, main_diagonal_legend, reverse_diagonal_legend],
               loc='upper left', bbox_to_anchor=(1.05, 0.6), borderaxespad=0., fontsize=10)

    plt.tight_layout() 
    plt.show()

# Test Case for protein sequences
protein1 = "MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQANNLR"
protein2 = "MKTAYIAKQRQISFVKSHFSRQLEER"
plot_protein_sequence_alignment(protein1, protein2)
