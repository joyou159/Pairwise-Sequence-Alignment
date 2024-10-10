import matplotlib.pyplot as plt
from NW import *

def plot_needleman_wunsch(protein1, protein2, scoring_scheme):
    """
    Plot the Needleman-Wunsch scoring matrix and save the alignment result if specified.

    Parameters:
    protein1 (str): The first protein sequence.
    protein2 (str): The second protein sequence.
    scoring_scheme (dict): The scoring parameters used for alignment.
    file_name (str, optional): The name of the file to save the alignment. Defaults to None.
    """
    # Generate the scores matrix and the optimal alignment paths
    scores = needleman_wunsch(protein1, protein2, scoring_scheme)
    all_paths = find_all_paths(scores, protein1, protein2, scoring_scheme)

    # Create the plot
    plt.figure(figsize=(7, 6))
    im = plt.imshow(scores, cmap="coolwarm", aspect='auto')
    
    # Adjust colorbar size
    plt.colorbar(im, label='Alignment Score', shrink=0.75)

    # Display the score values in each cell
    for i in range(scores.shape[0]):
        for j in range(scores.shape[1]):
            plt.text(j, i, f'{int(scores[i, j])}', ha='center', va='center', color='black', fontsize=10)

    # Generate colors for each path
    colors = plt.cm.viridis(np.linspace(0, 1, len(all_paths)))

    # Plot each optimal path in a different color
    for idx, path in enumerate(all_paths):
        path_y, path_x = zip(*path)
        plt.plot(path_x, path_y, color=colors[idx], marker='o', markersize=20, label=f'Path {idx + 1}', alpha=0.5, linewidth=2)

    plt.xlabel("Protein 2", fontsize=14)
    plt.ylabel("Protein 1", fontsize=14)

    # Move Protein 2 labels to the top
    plt.gca().xaxis.set_ticks_position('top')
    plt.gca().xaxis.set_label_position('top')
    plt.xticks(ticks=np.arange(1, len(protein2) + 1), labels=list(protein2), fontsize=12)
    plt.yticks(ticks=np.arange(1, len(protein1) + 1), labels=list(protein1), fontsize=12)

    # Move legend outside of the plot
    plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0, ncol=2, borderpad=1)

    # Show plot with layout adjustment
    plt.tight_layout()
    plt.show()


protein1 = "CTATGAA"
protein2 = "CTATTGACGTA"
scoring_scheme = {'match_score': 5, 'mismatch_penalty': -2, 'gap_penalty': -4}
plot_needleman_wunsch(protein1, protein2, scoring_scheme)
