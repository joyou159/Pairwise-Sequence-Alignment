import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import ListedColormap
from util import get_diagonals

def plot_sequence_alignment(sequence1, sequence2):
    """
    Generate a dot plot for pairwise sequence alignment with identified diagonals.
    
    This function creates a dot plot that visualizes the pairwise alignment of two sequences.
    Matches between the sequences are indicated by red squares, and diagonals (both main and reverse)
    are marked with black and green lines, respectively.

    Args:
        sequence1 (str): The first nucleotide or protein sequence.
        sequence2 (str): The second nucleotide or protein sequence.
    
    Returns:
        None: Displays the dot plot using matplotlib.
    """
    alignment_matrix = np.zeros((len(sequence1), len(sequence2), 3)) 

    # Populate the alignment matrix (1 indicates a match between the sequences)
    for i in range(len(sequence1)):
        for j in range(len(sequence2)):
            if sequence1[i] == sequence2[j]:
                alignment_matrix[i, j, 0] = 1

    # Extract main and reverse diagonals
    main_diagonals, reverse_diagonals = get_diagonals(alignment_matrix)  

    # Extract only the match (0 or 1) part for plotting
    plot_matrix = alignment_matrix[:, :, 0]

    # Custom colormap for the plot
    colors = [(1, 1, 1), (1, 0, 0)]  # Red for matches, white for mismatches
    cmap = ListedColormap(colors)

    # Create the plot
    plt.pcolor(plot_matrix[::-1], cmap=cmap, edgecolors='k', linewidths=1)

    # Plot main and reverse diagonals
    for diagonal in main_diagonals:
        plt.plot(diagonal[0], diagonal[1], ".-", color='black', linewidth=3, markersize=10)
    for diagonal in reverse_diagonals:
        plt.plot(diagonal[0], diagonal[1], ".-", color='green', linewidth=3, markersize=10)

    # Set axis labels and ticks
    plt.xticks(np.arange(len(sequence2)) + 0.5, list(sequence2))
    plt.yticks(np.arange(len(sequence1)) + 0.5, list(sequence1[::-1]))

    plt.xlabel("$Sequence 2$")
    plt.ylabel("$Sequence 1$")
    
    # Move the x-axis to the top
    plt.gca().xaxis.tick_top()
    plt.gca().xaxis.set_label_position('top')

    # Create legend
    match_legend = plt.Rectangle((0, 0), 1, 1, color='red', label='Match')
    main_diagonal_legend = plt.Line2D([0], [0], marker='.', linestyle="-", color='black', markersize=10, label='Main Diagonal')
    reverse_diagonal_legend = plt.Line2D([0], [0], marker='.', linestyle="-", color='green', markersize=10, label='Reverse Diagonal')
    plt.legend(handles=[match_legend, main_diagonal_legend, reverse_diagonal_legend], bbox_to_anchor=(1.4, 1))

    plt.show()



# Test Case   
sequence1 = "CTATTGACGTAACAT"
sequence2 = "CTATTGAACAT"
plot_sequence_alignment(sequence1, sequence2)