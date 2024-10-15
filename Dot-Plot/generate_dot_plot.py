import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import ListedColormap
from util import get_diagonals

def plot_sequence_alignment(seq1, seq2):
    """
    Generate a dot plot for pairwise sequence alignment with identified diagonals.

    This function creates a dot plot that visualizes the pairwise alignment of two sequences.
    Matches between the sequences are indicated by red squares, and diagonals (both main and reverse)
    are marked with black and green lines, respectively.

    Args:
        seq1 (str or list): The first sequence (could be any type of sequence).
        seq2 (str or list): The second sequence (could be any type of sequence).
    
    Returns:
        None: Displays the dot plot using matplotlib.
    """
    # Initialize an empty alignment matrix for matches
    alignment_matrix = np.zeros((len(seq1), len(seq2), 3))

    # Populate the alignment matrix: 1 indicates a match between elements
    for i in range(len(seq1)):
        for j in range(len(seq2)):
            if seq1[i] == seq2[j]:
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

    plt.xticks(np.arange(len(seq2)) + 0.5, list(seq2))
    plt.yticks(np.arange(len(seq1)) + 0.5, list(seq1[::-1]))

    plt.xlabel("Sequence 2")
    plt.ylabel("Sequence 1")
    
    # Move the x-axis to the top
    plt.gca().xaxis.tick_top()
    plt.gca().xaxis.set_label_position('top')

    # Create the legends for matches and diagonals
    match_legend = plt.Rectangle((0, 0), 1, 1, color='red', label='Element Match')
    main_diagonal_legend = plt.Line2D([0], [0], marker='.', linestyle="-", color='black', markersize=10, label='Main Diagonal')
    reverse_diagonal_legend = plt.Line2D([0], [0], marker='.', linestyle="-", color='green', markersize=10, label='Reverse Diagonal')
    
    plt.legend(handles=[match_legend, main_diagonal_legend, reverse_diagonal_legend],
               loc='upper left', bbox_to_anchor=(1.05, 0.6), borderaxespad=0., fontsize=10)

    plt.tight_layout() 
    plt.show()

# DNA Test Case  
seq1 = "CTATTGACGTA"
seq2 = "CTATGAA"
plot_sequence_alignment(seq1, seq2)

