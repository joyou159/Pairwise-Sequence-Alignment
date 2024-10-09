def diagonal_traversing(i, j, alignment_matrix, direction='main'):
    """
    Traverse a diagonal (main or reverse) of the alignment matrix starting from position (i, j).
    
    The function looks for continuous matches along the specified diagonal (main or reverse) and
    returns the x and y coordinates of the matching sequence if found. If only a single match is found,
    the diagonal is rejected.
    
    Args:
        i (int): Starting row index of the traversal.
        j (int): Starting column index of the traversal.
        alignment_matrix (ndarray): A 3D numpy array representing the alignment matrix.
            The first dimension indicates whether two elements match (0 or 1),
            the second dimension tracks if the cell has been checked for a main diagonal,
            and the third dimension tracks if it has been checked for a reverse diagonal.
        direction (str): The direction of diagonal traversal, either 'main' (top-left to bottom-right)
                         or 'reverse' (top-right to bottom-left).
    
    Returns:
        list: A list containing two lists [x_values, y_values] with coordinates of the diagonal
              if it contains more than one matching cell. Otherwise, an empty list is returned.
    """
    x_values = []
    y_values = []
    row = i
    column = j
    
    while alignment_matrix[row, column, 0]:
        if direction == 'main':
            alignment_matrix[row, column, 1] = 1  # Mark the cell as checked for the main diagonal
        else:
            alignment_matrix[row, column, 2] = 1  # Mark the cell as checked for the reverse diagonal
        
        # Center ticks and reverse y-axis
        x_values.append((column + 0.5))  
        y_values.append((alignment_matrix.shape[0] - row - 0.5))
        
        # Move in the specified direction
        row += 1
        if direction == 'main':
            column += 1
        else:
            column -= 1
        
        # Avoid indexing errors
        if row == alignment_matrix.shape[0] or column < 0 or column == alignment_matrix.shape[1]:
            break
    
    # Reject if only a single match is found
    if len(x_values) == 1:
        return []
    else:
        return [x_values, y_values]


def get_diagonals(alignment_matrix):
    """
    Identify and extract both main and reverse diagonals from the alignment matrix.
    
    This function iterates through the alignment matrix and calls the appropriate diagonal traversing
    functions to extract the main and reverse diagonals containing continuous matches. Each cell in the matrix
    can only contribute to one main diagonal and one reverse diagonal.
    
    Args:
        alignment_matrix (ndarray): A 3D numpy array representing the alignment matrix. The first dimension
            indicates whether two elements match (0 or 1), the second tracks if the cell has been checked
            for a main diagonal, and the third tracks if it has been checked for a reverse diagonal.

    Returns:
        tuple: A tuple (main_diagonals, reverse_diagonals) where each is a list of diagonal coordinates.
    """
    main_diagonals = []
    reverse_diagonals = []

    for i in range(alignment_matrix.shape[0]):
        for j in range(alignment_matrix.shape[1]):
            # Check for unvisited main diagonal
            if alignment_matrix[i, j, 0] == 1 and alignment_matrix[i, j, 1] == 0:
                curr_main = diagonal_traversing(i, j, alignment_matrix, direction='main')
                if curr_main:
                    main_diagonals.append(curr_main)
            
            # Check for unvisited reverse diagonal
            if alignment_matrix[i, j, 0] == 1 and alignment_matrix[i, j, 2] == 0:
                curr_reverse = diagonal_traversing(i, j, alignment_matrix, direction='reverse')
                if curr_reverse:
                    reverse_diagonals.append(curr_reverse)
    
    return main_diagonals, reverse_diagonals
