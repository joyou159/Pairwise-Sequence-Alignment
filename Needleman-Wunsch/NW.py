import numpy as np 

def needleman_wunsch(protein1, protein2, scoring_scheme):
    """
    Implements the Needleman-Wunsch algorithm for global sequence alignment.

    Args:
        protein1 (str): The first protein sequence (or any sequence of characters).
        protein2 (str): The second protein sequence (or any sequence of characters).
        scoring_scheme (dict): A dictionary containing the following keys:
            - 'match_score': Score for matching characters.
            - 'mismatch_penalty': Penalty for mismatched characters.
            - 'gap_penalty': Penalty for introducing gaps.

    Returns:
        numpy.ndarray: A 2D array (scoring matrix) where each element represents the optimal alignment score
                       between the prefixes of the two sequences up to that point.
    """
    match_score = scoring_scheme['match_score']
    mismatch_penalty = scoring_scheme['mismatch_penalty']
    gap_penalty = scoring_scheme['gap_penalty']
    
    n, m = len(protein1), len(protein2)
    scores = np.zeros((n + 1, m + 1))
    
    # Initialize the first row and column with gap penalties
    for i in range(1, n + 1):
        scores[i, 0] = gap_penalty * i
    for j in range(1, m + 1):
        scores[0, j] = gap_penalty * j
    
    # Fill the scoring matrix
    for i in range(1, n + 1):
        for j in range(1, m + 1):
            match = scores[i - 1, j - 1] + (match_score if protein1[i - 1] == protein2[j - 1] else mismatch_penalty)
            delete = scores[i - 1, j] + gap_penalty
            insert = scores[i, j - 1] + gap_penalty
            scores[i, j] = max(match, delete, insert)
    
    return scores

def traceback_path(scores, protein1, protein2, scoring_scheme):
    """
    Reconstructs the optimal alignment path from the scoring matrix.

    Args:
        scores (numpy.ndarray): The scoring matrix from the Needleman-Wunsch algorithm.
        protein1 (str): The first protein sequence (or any sequence of characters).
        protein2 (str): The second protein sequence (or any sequence of characters).
        scoring_scheme (dict): A dictionary containing the following keys:
            - 'match_score': Score for matching characters.
            - 'mismatch_penalty': Penalty for mismatched characters.
            - 'gap_penalty': Penalty for introducing gaps.

    Returns:
        list: A list of tuples representing the path through the scoring matrix (traceback path) from the
              bottom-right to the top-left, corresponding to the optimal alignment.
    """
    match_score = scoring_scheme['match_score']
    mismatch_penalty = scoring_scheme['mismatch_penalty']
    gap_penalty = scoring_scheme['gap_penalty']
    
    i, j = len(protein1), len(protein2)
    path = [(i, j)]
    
    while i > 0 or j > 0:
        if i > 0 and j > 0 and (scores[i, j] == scores[i - 1, j - 1] + 
                                (match_score if protein1[i - 1] == protein2[j - 1] else mismatch_penalty)):
            i, j = i - 1, j - 1
        elif i > 0 and scores[i, j] == scores[i - 1, j] + gap_penalty:
            i -= 1
        else:
            j -= 1
        path.append((i, j))
    
    return path
