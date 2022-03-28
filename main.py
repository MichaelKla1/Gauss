
matrix = [[-1.41, 2, 0],
[1, -1.41, 1],
      [0, 2,-1.41]
      ]


resvec = [1,1,1]


def buildIdentityMatrix(n):
    """
    Builds and returns Identity matrix nxn. if n is not integer returns None
    :param n: size of matrix (nxn)
    :return: Identity matrix nxn
    """
    if not isinstance(n, int):
        return None
    resmat = []
    i = 0

    while i < n:
        resrow = []
        j = 0
        while j < n:
            if i == j:
                resrow.append(1)
            else:
                resrow.append(0)
            j = j + 1
        resmat.append(resrow)
        i = i + 1
    return resmat


def isMatrix(m):
    """
    Returns True if m is matrix and False otherwise
    :param m: parameter to be checked if it is matrix
    :return: True if m is matrix and False otherwise
    """
    try:
        # check if input are matrixes
        iter(m)
        if len(m) < 1:
            return False  # matrix is empty
        matcolcount = len(m[0])
        for i in m:
            iter(i)
            if len(i) != matcolcount:
                return False  # columns of matrix have different number of elements
    except TypeError as e:
        return False  # input is not matrixes
    return True


def multiplyMatrix(matrix1, matrix2):
    """
    Return the result of multiplication matrix1 and matrix2
    :param matrix1: first matrix
    :param matrix2: second matrix
    :return: the result of multiplication matrix1 and matrix2
    """
    # check if input is correct
    check = any([True for item in matrix2 if type(item) is list])  # check if the second matrix is b vector in Ax=b
    try:
        # check if input are matrixes
        #if not isMatrix(matrix1) or not isMatrix(matrix2):
        #    return None  # one or both inputs are not matrixes
        iter(matrix1)
        iter(matrix2)
        mat1colcount = len(matrix1[0])
        mat1rowcount = len(matrix1)
        if check:
            mat2colcount = len(matrix2[0])
            mat2rowcount = len(matrix2)
        else:
            mat2colcount = 1;
            mat2rowcount = len(matrix2)
        if (mat1colcount != mat2rowcount):
            return None  # can't multiply matrixes

        # check if matrixes can be multiplied
    except TypeError as e:
        return None  # input is not matrixes
    k = 0
    resmat = []
    while k < mat1rowcount:
        i = 0
        resrow = []
        while i < mat2colcount:
            sum = 0
            j = 0
            while j < mat1colcount:
                if check:
                    sum = sum + matrix1[k][j] * matrix2[j][i]
                else:
                    sum = sum + matrix1[k][j] * matrix2[j]
                j = j + 1
            resrow.append(sum)
            i = i + 1
        if check:
            resmat.append(resrow)
        else:
            resmat.append(sum)
        k = k + 1
    return resmat
def print_mat(ls):#print the matrix
    for i in ls:
        print(f'{i}')
    print('\n')
def print_mul(m1, m2):#print the matrix in multiplication
    length = len(m1)
    for i in range(length):
        if length % 2 == 0 and length / 2 == i:
            print(f'({m1[i]}) x ({m2[i]}) = ')
            continue
        elif length % 2 == 1 and int(length / 2) == i:
            print(f'({m1[i]}) x ({m2[i]}) = ')
            continue
        print(f'({m1[i]})   ({m2[i]})')
    print('\n')
def solutionOfMatrixByGaussElimination(matrix, resvec):
    """
    Return the solution vector of square matrix that has result vector resvec and also print the elementary matrixes that was using during Gauss elimination
    :param matrix: matrix
    :param resvec: result vector
    :return: solution vector of square matrix that has result vector resvec
    """
    if not isMatrix(matrix):
        return []
    try:
        iter(resvec)
    except TypeError as e:
        return []  # resvec is not vector
    rowscount = len(matrix)
    if rowscount != len(resvec):
        return []  # number of rows in matrix doesn't match number of rows in result vector
    #matrix = copy.deepcopy(matrix)
    #resvec = copy.deepcopy(resvec)
    elementaryMatrixesList = []
    i = 0
    while i < rowscount:
        colscount = len(matrix[i])
        if colscount != rowscount:
            return []  # matrix is not squared
        rowdivider = matrix[i][i]
        s=i+1
        while s < rowscount:
            if abs(matrix[s][i]) > abs(rowdivider):
                # swap lines
                tempmat = buildIdentityMatrix(rowscount)
                temp = tempmat[s]
                tempmat[s] = tempmat[i]
                tempmat[i] = temp
                print_mul(tempmat, matrix)
                matrix = multiplyMatrix(tempmat,matrix)
                print_mat(matrix)

                flag = 1
                rowdivider = matrix[i][i]
                print_mul(tempmat, resvec)
                resvec = multiplyMatrix(tempmat,resvec)
                print_mat(resvec)
                #temp = resvec[s]
                #resvec[s] = resvec[i]
                #resvec[i] = temp

            s += 1
        #if flag == 0:
         #   return []  # same lines in matrix
        elementaryMatrix = buildIdentityMatrix(rowscount)
        elementaryMatrix[i][i] /= rowdivider
        print_mul(elementaryMatrix, matrix)
        matrix = multiplyMatrix(elementaryMatrix, matrix)
        print_mat(matrix)
        elementaryMatrixesList.insert(0, elementaryMatrix)
        print_mul(elementaryMatrix, resvec)
        resvec = multiplyMatrix(elementaryMatrix,resvec)
        print_mat(resvec)
        #resvec[i] = resvec[i] / rowdivider
        k = i
        while k < rowscount:
            if k != i:
                n = matrix[k][i]
                elementaryMatrix = buildIdentityMatrix(rowscount)
                elementaryMatrix[k][i] = -n * matrix[i][i]
                print_mul(elementaryMatrix, matrix)
                matrix = multiplyMatrix(elementaryMatrix, matrix)
                print_mat(matrix)
                elementaryMatrixesList.insert(0, elementaryMatrix)
                print_mul(elementaryMatrix, matrix)
                resvec = multiplyMatrix(elementaryMatrix,resvec)
                print_mat(resvec)

                #resvec[k] = resvec[k] - n * resvec[i]
            k += 1
        i += 1

    i=0
    while i < rowscount:
        colscount = len(matrix[i])
        if colscount != rowscount:
            return []  # matrix is not squared
        rowdivider = matrix[i][i]

        #if flag == 0:
         #   return []  # same lines in matrix
        elementaryMatrix = buildIdentityMatrix(rowscount)
        elementaryMatrix[i][i] /= rowdivider
        print_mul(elementaryMatrix, matrix)
        matrix = multiplyMatrix(elementaryMatrix, matrix)
        print_mat(matrix)
        elementaryMatrixesList.insert(0, elementaryMatrix)
        print_mul(elementaryMatrix, resvec)
        resvec = multiplyMatrix(elementaryMatrix,resvec)
        print_mat(resvec)
        #resvec[i] = resvec[i] / rowdivider
        k = 0
        while k < i:
            if k != i:
                n = matrix[k][i]
                elementaryMatrix = buildIdentityMatrix(rowscount)
                elementaryMatrix[k][i] = -n * matrix[i][i]
                print_mul(elementaryMatrix, matrix)
                matrix = multiplyMatrix(elementaryMatrix, matrix)
                print_mat(matrix)
                elementaryMatrixesList.insert(0, elementaryMatrix)
                print_mul(elementaryMatrix, matrix)
                resvec = multiplyMatrix(elementaryMatrix,resvec)
                print_mat(resvec)

                #resvec[k] = resvec[k] - n * resvec[i]
            k += 1
        i += 1

    counter = len(elementaryMatrixesList)
    """for i in elementaryMatrixesList: #print the elementary matrixes list
        print("E" + str(counter) + ": " + str(i))
        counter -= 1"""
    return resvec

def print_ans(mat):
    [print(f'x{i+1} = {mat[i]}') for i in range(len(mat))]
print_ans(solutionOfMatrixByGaussElimination(matrix, resvec))