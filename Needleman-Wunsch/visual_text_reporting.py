import matplotlib.pyplot as plt
import os
from NW import *


def save_alignment_to_file(protein1, protein2, paths, optimal_score, file_name):
    """
    Save the alignment result to a file.

    Parameters:
    protein1 (str): The first protein sequence.
    protein2 (str): The second protein sequence.
    path (list): The optimal alignment path.
    optimal_score (int): The alignment score. 
    file_name (str): The name of the file to save the alignment.
    """
    if os.path.exists(file_name):
        os.remove(file_name)

    for ind, path in enumerate(paths):
        aligned_protein1 = []
        aligned_protein2 = []
        i_cache, j_cache = path[-1]
        for (i, j) in list(reversed(path))[1:]: # neglect for caching(0, 0)
            if i == i_cache and j != j_cache: # moved horizontally
                aligned_protein1.append('-')
                aligned_protein2.append(protein2[j - 1])  
            elif i != i_cache and j == j_cache: # moved vertically
                aligned_protein1.append(protein1[i - 1])
                aligned_protein2.append('-')
            elif i != i_cache and j != j_cache: # moved diagonally 
                aligned_protein1.append(protein1[i - 1])
                aligned_protein2.append(protein2[j - 1])
            i_cache, j_cache = i, j
    
        with open(file_name, 'a') as file:
            file.write(f"Alignment Result {ind+1} with Score ({optimal_score}):\n" + ''.join(aligned_protein1) + '\n' + ''.join(aligned_protein2) + '\n'+ "------------------" + '\n')

def plot_needleman_wunsch(protein1, protein2, scoring_scheme, file_name=None):
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

    plt.figure(figsize=(7, 6))
    im = plt.imshow(scores, cmap="coolwarm", aspect='auto')
    
    plt.colorbar(im, label='Alignment Score', shrink=0.75)

    for i in range(scores.shape[0]):
        for j in range(scores.shape[1]):
            plt.text(j, i, f'{int(scores[i, j])}', ha='center', va='center', color='black', fontsize=10)

    colors = plt.cm.viridis(np.linspace(0, 1, len(all_paths)))

    # Plot each optimal path in a different color
    for idx, path in enumerate(all_paths):
        path_y, path_x = zip(*path)
        plt.plot(path_x, path_y, color=colors[idx], marker='o', markersize=20, label=f'Path {idx + 1}', alpha=0.5, linewidth=2)

    plt.xlabel("Protein 2", fontsize=14)
    plt.ylabel("Protein 1", fontsize=14)

    plt.gca().xaxis.set_ticks_position('top')
    plt.gca().xaxis.set_label_position('top')
    plt.xticks(ticks=np.arange(1, len(protein2) + 1), labels=list(protein2), fontsize=12)
    plt.yticks(ticks=np.arange(1, len(protein1) + 1), labels=list(protein1), fontsize=12)

    plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0, labelspacing=1.2)

    plt.tight_layout()
    plt.show()
    print(all_paths)

    # report the alignment results.
    if file_name:
        save_alignment_to_file(protein1, protein2, all_paths, scores[-1,-1],file_name)



protein1 = "CTATTGACGTA"
protein2 = "CTATGAA"
scoring_scheme = {'match_score': 5, 'mismatch_penalty': -2, 'gap_penalty': -4}
plot_needleman_wunsch(protein1, protein2, scoring_scheme, "Needleman-Wunsch/alignment_result.txt")
