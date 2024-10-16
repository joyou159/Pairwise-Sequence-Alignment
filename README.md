# Pairwise Sequence Alignment

This repository provides implementations for basic sequence alignment techniques, focusing on two popular methods: Dot Plot and Needleman-Wunsch algorithm. These techniques are widely used in bioinformatics to compare biological sequences, such as DNA, RNA, or protein sequences.

# Alignment Techniques

## Needleman-Wunsch Algorithm
### Overview
The Needleman-Wunsch algorithm is a global alignment technique used to align entire sequences from end to end. It uses a dynamic programming approach to find the optimal alignment based on a scoring scheme.

  <p align="center">
  <img src="READMD-Assets\NW.png" alt="Global alignment using NW algorithm" title="Global alignment using NW algorithm" width="450" />
  </p>

### How It Works
1. **Initialization**: A scoring matrix is created with the sequences along the horizontal and vertical axes. The first row and column are initialized with gap penalties.
2. **Matrix Filling**: The matrix is filled using a scoring scheme (match, mismatch, and gap penalties) to compute the optimal alignment scores.
3. **Traceback**: Starting from the bottom-right corner of the matrix, the algorithm traces back to the top-left corner to determine the alignment path.

### Scoring Scheme
The scoring scheme for the Needleman-Wunsch algorithm includes:
- **Match Score**: Positive score when two characters are identical.
- **Mismatch Penalty**: Negative score when characters do not match.
- **Gap Penalty**: Negative score for introducing a gap in the alignment.

The optimal alignment is found by maximizing the alignment score.

### Applications
- Comparing two complete sequences (e.g., aligning entire protein sequences).
- Studying evolutionary relationships by finding the best global alignment.
- Serving as a foundation for other alignment algorithms, such as Smith-Waterman.


## Dot Plot
### Overview
The Dot Plot method is a simple graphical approach used to compare two sequences. It displays similarities between the sequences in a matrix form, where each axis represents one of the sequences. 

  <p align="center">
  <img src="READMD-Assets\Dot_plot.png" alt="Highlighing the similarities between two sequences using dot plot" title="Highlighing the similarities between two sequences using dot plot" width="450" />
  </p>

### How It Works
1. The sequences are placed along the horizontal and vertical axes of a matrix.
2. A dot is placed in the matrix at positions where the corresponding elements of the sequences are identical (or similar based on a threshold).
3. The resulting pattern shows regions of similarity, such as diagonals indicating consecutive matches or repeating patterns.

### Applications
- Visual identification of repeating sequences.
- Locating regions of high similarity between sequences.
- Detecting inversions or translocations in genomic data.


## Usage

### Installation

1. **Clone the Repository:**
    ```bash
    git clone git@github.com:joyou159/Pairwise-Sequence-Alignment.git
    cd Pairwise-Sequence-Alignment
    ```
2. **Install Dependencies:**
    ```bash
    pip install -r requirements.txt
    ```

### Running the Code

-  **Dot Plot**
   - Use the provided function to generate a dot plot for two sequences:
     ```python
     from generate_dot_plot import plot_dot_plot
     sequence1 = "CTATTGACGTA"
     sequence2 = "CTATGAA"
     plot_dot_plot(sequence1, sequence2)
     ```

-  **Needleman-Wunsch Algorithm**
   - Plot the alignment scoring matrix and optionally save the alignment result:
     ```python
     from NW_alignment import plot_alignment
     sequence1 = "CTATTGACGTA"
     sequence2 = "CTATGAA"
     scoring_scheme = {'match_score': 5, 'mismatch_penalty': -2, 'gap_penalty': -4}
     plot_alignment(sequence1, sequence2, scoring_scheme, "./alignment_result.txt")
     ```


## References
- Needleman, S. B., & Wunsch, C. D. (1970). A general method applicable to the search for similarities in the amino acid sequence of two proteins. *Journal of Molecular Biology*.
- Mount, D. W. (2004). *Bioinformatics: Sequence and Genome Analysis*. Cold Spring Harbor Laboratory Press.

