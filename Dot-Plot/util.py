def main_diagonal_traversing(i, j, alignment_matrix):
    """
    Traverse the main diagonal of the alignment matrix starting from position (i, j).
    
    The function looks for continuous matches along the main diagonal (i.e., from top-left to bottom-right)
    and returns the x and y coordinates of the matching sequence if found. If only a single match is found,
    the diagonal is rejected.
    
    Args:
        i (int): Starting row index of the traversal.
        j (int): Starting column index of the traversal.
        alignment_matrix (ndarray): A 3D numpy array representing the alignment matrix.
            The first dimension indicates whether two nucleotides match (0 or 1),
            the second dimension tracks if the cell has been checked for a main diagonal,
            and the third dimension tracks if it has been checked for a reverse diagonal.

    Returns:
        list: A list containing two lists [x_values, y_values] with coordinates of the diagonal
              if it contains more than one matching cell. Otherwise, an empty list is returned.
    """
    x_values = []
    y_values = []
    row = i 
    column = j
    while alignment_matrix[row, column, 0]:
        alignment_matrix[row, column, 1] = 1  # Mark the cell as checked for the main diagonal
        x_values.append((column + 0.5))  # Center ticks
        y_values.append((alignment_matrix.shape[0] - row - 0.5))  # Center ticks and reverse y-axis
        row += 1 
        column += 1
        if row == alignment_matrix.shape[0] or column == alignment_matrix.shape[1]:  # Avoid indexing errors
            break
    
    if len(x_values) == 1:  # Reject if only a single match is found
        return []
    else:
        return [x_values, y_values]
    


def reverse_diagonal_traversing(i, j, alignment_matrix):
    """
    Traverse the reverse diagonal of the alignment matrix starting from position (i, j).
    
    The function looks for continuous matches along the reverse diagonal (i.e., from top-right to bottom-left)
    and returns the x and y coordinates of the matching sequence if found. If only a single match is found,
    the diagonal is rejected.
    
    Args:
        i (int): Starting row index of the traversal.
        j (int): Starting column index of the traversal.
        alignment_matrix (ndarray): A 3D numpy array representing the alignment matrix.
            The first dimension indicates whether two nucleotides match (0 or 1),
            the second dimension tracks if the cell has been checked for a main diagonal,
            and the third dimension tracks if it has been checked for a reverse diagonal.

    Returns:
        list: A list containing two lists [x_values, y_values] with coordinates of the diagonal
              if it contains more than one matching cell. Otherwise, an empty list is returned.
    """
    x_values = []
    y_values = []
    row = i
    column = j
    while alignment_matrix[row, column, 0]:
        alignment_matrix[row, column, 2] = 1  # Mark the cell as checked for the reverse diagonal
        x_values.append((column + 0.5))  # Center ticks
        y_values.append((alignment_matrix.shape[0] - row - 0.5))  # Center ticks and reverse y-axis
        row += 1  
        column -= 1  # Traverse in reverse direction (decrement column)
        if row == alignment_matrix.shape[0] or column == -1:  # Avoid indexing errors
            break
    
    if len(x_values) == 1:  # Reject if only a single match is found
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
            indicates whether two nucleotides match (0 or 1), the second tracks if the cell has been checked
            for a main diagonal, and the third tracks if it has been checked for a reverse diagonal.

    Returns:
        tuple: A tuple (main_diagonals, reverse_diagonals) where each is a list of diagonal coordinates.
    """
    main_diagonals = []
    reverse_diagonals = []

    for i in range(alignment_matrix.shape[0]):
        for j in range(alignment_matrix.shape[1]):
            if alignment_matrix[i, j, 0] == 1 and alignment_matrix[i, j, 1] == 0:  # Check for unvisited main diagonal
                curr_main = main_diagonal_traversing(i, j, alignment_matrix)
                if curr_main:
                    main_diagonals.append(curr_main)
            if alignment_matrix[i, j, 0] == 1 and alignment_matrix[i, j, 2] == 0:  # Check for unvisited reverse diagonal
                curr_reverse = reverse_diagonal_traversing(i, j, alignment_matrix)
                if curr_reverse:
                    reverse_diagonals.append(curr_reverse)
    
    return main_diagonals, reverse_diagonals