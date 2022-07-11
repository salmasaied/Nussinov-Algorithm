import numpy as np
complementer = {'A':'U','U':'A','G':'C','C':'G'}


def intialize_matrix(sequence):
    '''
    Creates an intial LxL matrix; L = length of provided sequence
    INPUT:
        sequence: string. The given RNA sequence to predict its secondary structure.
    OUTPUT:
        matrix: 2D Array. The intial matrix.
    NOTE:
        We decided to fill null-values with 0 so its easier to trace by the eye.
    '''
    L = len(sequence)
    matrix = [[0 for i in range(L)] for j in range(L)]
    return matrix


def fill_matrix(sequence,matrix):
    '''
    Fills the inital matrix according to Nussinov's algorithm.
    INPUT:
        sequence: string. The given RNA sequence to predict its secondary structure.
        matrix: 2D Array. The intial matrix.
    OUTPUT:
        matrix: 2D Array. The solution matrix.
    '''
    L = len(matrix)
    '''Loop in diagonal stripes'''
    for y in range(2,L+1):
        for x in range(y,L+1):
            i = x - y
            j = x - 1
            Lower_D = matrix[i+1][j-1] + 1 if sequence[i] == complementer[sequence[j]] else matrix[i+1][j-1]
            Left = matrix[i][j-1]
            Bottom = matrix[i+1][j]
            fourth_case = max([matrix[i][k] + matrix[k+1][j] for k in range(i,j)])
            matrix[i][j] = max([Lower_D,Left,Bottom,fourth_case])
    return matrix


def traceback(sequence,matrix):
    '''
    Performs the traceback step over the matrix
    INPUT:
        sequence: string. The given RNA sequence to predict its secondary structure.
        matrix: 2D Array. The solution matrix.
    OUTPUT:
        .: string. RNA secondary structure.        
    '''
    L = len(matrix)
    i = 0
    j = L - 1
    return_string = list('.'*L)
    while matrix[i][j] != 0:
        if matrix[i][j] == matrix[i+1][j-1] + 1 and sequence[i] == complementer[sequence[j]] :
            return_string[i] = '('
            return_string[j] = ')'
            i = i + 1
            j = j - 1
        elif matrix[i][j] == matrix[i][j-1]:
            j = j - 1
        elif matrix[i][j] == matrix[i+1][j]:
            i = i + 1
        else:
            fourth_case_list = [matrix[i][k] + matrix[k+1][j] for k in range(i,j)]
            k = np.argmax(fourth_case_list) + i + 1
            max_val = np.max(fourth_case_list)
            if matrix[i][j] == max_val:
                j = k
    return ''.join(return_string)


def nussinov(sequence):
    '''
    Our implementation to the Nussinov Alogrithm,assembels both of the matrix filling
    & traceback steps together.
    INPUT:
        sequence: string. The given RNA sequence to predict its secondary structure.
    OUTPUT:
        solution: string. RNA secondary structure.
        filled_matrix: 2D Array. The solution matrix.
    '''
    matrix = intialize_matrix(sequence)
    filled_matrix = fill_matrix(sequence,matrix)
    solution = traceback(sequence,filled_matrix)
    print(' ',*sequence,sep = '  ')
    for i in range(len(sequence)):
        print(sequence[i],filled_matrix[i])
    print("Secondary Structure Prediction : ",solution)
    return solution,filled_matrix


nussinov('GGGAAAUCC') #First Example


    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    