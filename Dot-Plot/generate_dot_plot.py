import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import ListedColormap
from util import get_diagonals

def plot_sequence_alignment(sequence1, sequence2):
    """
    Generate a dot plot for pairwise sequence alignment with identified diagonals.
    
    This function creates a dot plot that visualizes the pairwise alignment of two sequences
    (DNA, RNA, or protein). Matches between the sequences are indicated by red squares, 
    and diagonals (both main and reverse) are marked with black and green lines, respectively.

    Args:
        sequence1 (str): The first sequence (DNA, RNA, or protein).
        sequence2 (str): The second sequence (DNA, RNA, or protein).
    
    Returns:
        None: Displays the dot plot using matplotlib.
    """
    # Initialize an empty alignment matrix for matches
    alignment_matrix = np.zeros((len(sequence1), len(sequence2), 3))

    # Populate the alignment matrix: 1 indicates a match
    for i in range(len(sequence1)):
        for j in range(len(sequence2)):
            if sequence1[i] == sequence2[j]:
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

    plt.xticks(np.arange(len(sequence2)) + 0.5, list(sequence2))
    plt.yticks(np.arange(len(sequence1)) + 0.5, list(sequence1[::-1]))

    plt.xlabel("$Sequence 2$")
    plt.ylabel("$Sequence 1$")
    
    # Move the x-axis to the top
    plt.gca().xaxis.tick_top()
    plt.gca().xaxis.set_label_position('top')

    match_legend = plt.Rectangle((0, 0), 1, 1, color='red', label='Match')
    main_diagonal_legend = plt.Line2D([0], [0], marker='.', linestyle="-", color='black', markersize=10, label='Main Diagonal')
    reverse_diagonal_legend = plt.Line2D([0], [0], marker='.', linestyle="-", color='green', markersize=10, label='Reverse Diagonal')
    
    plt.legend(handles=[match_legend, main_diagonal_legend, reverse_diagonal_legend],
               loc='upper left', bbox_to_anchor=(1.05, 0.6), borderaxespad=0., fontsize=10)

    plt.tight_layout() 
    plt.show()

# Test Case
sequence1 = "MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQANNLR"
sequence2 = "MKTAYIAKQRQISFVKSHFSRQLEER"
plot_sequence_alignment(sequence1, sequence2)
