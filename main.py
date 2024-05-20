from typing import List, TypeVar
from math import sqrt
import copy, os

T = TypeVar('T')
Matrix = List[List[T]]

def calculoDeterminante(my_order: int, my_matrix: Matrix) -> int:
    os.system('clear')
    sum = 0

    order = my_order
    matrix = copy.deepcopy(my_matrix)

    if order <= 0:
        print("Ordem invalida! Insira outra ordem")
        return -1

    if order == 1:
        return matrix[0][0]

    if order == 2:
        return matrix[0][0] * matrix[1][1] - matrix[1][0] * matrix[0][1]

    for col in list(range(len(matrix))):
        matrix_copy = copy.deepcopy(matrix[1:])

        for i in range(len(matrix_copy)):
            matrix_copy[i] = matrix_copy[i][0:col] + matrix_copy[i][col+1:]

        sum += (-1) ** (col % 2) * matrix[0][col] * calculoDeterminante(order - 1, matrix_copy)

    return int(sum)

def isLowerTriangular(matrix):
    rows = len(matrix)
    for i in range(rows):
        for j in range(i+1, rows):
            if matrix[i][j] != 0:
                return False
    return True

def sistemaTriangularInferior(my_order: int, my_matrix: Matrix, my_vector: List) -> List:
    os.system('clear')
    order = my_order
    matrix = copy.deepcopy(my_matrix)
    vector = copy.deepcopy(my_vector)

    if order == 0:
        print("Ordem invalida")
        return []
    
    if calculoDeterminante(order, matrix) == 0:
        print("Determinante igual a zero! Tente inserir outra matriz")
        return []

    if not isLowerTriangular(matrix):
        print("Matriz fornecida nao e triangular inferior")
        return []

    if order == 1:
        if matrix[0][0] == 0 and vector[0] == 0:
            print("Divisao por zero! Tente outra matriz e/ou vetor")
            return []

        if matrix[0][0] == 0:
            return vector[0]

        return vector[0] / matrix[0][0]

    solution_vector = [0] * len(vector)

    for i in range(len(matrix)):
        for j in range(i):
            solution_vector[i] -= matrix[i][j] * solution_vector[j]

        solution_vector[i] += vector[i]

        solution_vector[i] /= matrix[i][i]

    solution_vector = [round(x, 4) for x in solution_vector]

    return solution_vector

def isUpperTriangular(matrix):
    rows = len(matrix)
    for i in range(rows):
        for j in range(i):
            if matrix[i][j] != 0:
                return False
    return True

def sistemaTriangularSuperior(my_order: int, my_matrix: Matrix, my_vector: List) -> List:
    os.system('clear')
    order = my_order
    matrix = copy.deepcopy(my_matrix)
    vector = copy.deepcopy(my_vector)

    if order == 0:
        print("Ordem invalida")

    if calculoDeterminante(order, matrix) == 0:
        print("Determinante igual a zero! Tente inserir outra matriz")
        return []

    if not isUpperTriangular(matrix):
        print("Matriz fornecida nao e triangular superior")
        return []

    if order == 1:
        if matrix[-1][-1] == 0 and vector[-1] == 0:
            print("Divisao por zero! Tente outra matriz e/ou vetor")
            return []

        if matrix[-1][-1] == 0:
            return vector[0]

        return vector[-1] / matrix[-1][-1]

    solution_vector = [0] * len(vector)

    for i in reversed(range(len(matrix))):
        solution_vector[i] += vector[i]

        for j in range(i + 1, order):
            solution_vector[i] -= matrix[i][j] * solution_vector[j]

        solution_vector[i] /= matrix[i][i]

    solution_vector = [round(x, 4) for x in solution_vector]

    return solution_vector

def mult_matrix(my_matrix_m: Matrix, my_matrix_n: Matrix):
    matrix_m = copy.deepcopy(my_matrix_m)
    matrix_n = copy.deepcopy(my_matrix_n)

    transpose_n = list(zip(*matrix_n))
    return [[sum(el_m * el_n for el_m, el_n in zip(row_m, col_n)) for col_n in transpose_n] for row_m in matrix_m]

def pivot_matrix(my_order: int):
    order = my_order

    identity_matrix = [[float(i == j) for i in range(order)] for j in range(order)]

    return identity_matrix

def leadingPrincipalSubmatrix(my_matrix, my_order):
    matrix = copy.deepcopy(my_matrix)
    order = my_order

    return [row[:order] for row in matrix[:order]]

def isNonSingular(my_matrix, order):
    matrix = copy.deepcopy(my_matrix)

    for k in range(1, order + 1):
        submatrix = leadingPrincipalSubmatrix(matrix, k)
        if calculoDeterminante(k, submatrix) == 0:
            return False
    return True

def decomposicaoLU(my_order: int, my_matrix: Matrix, my_vector: List) -> List:
    os.system('clear')
    order = my_order
    matrix = copy.deepcopy(my_matrix)
    vector = copy.deepcopy(my_vector)

    if not isNonSingular(matrix, order):
        print("A matriz possui uma submatriz singular")
        return [] 
    
    L = [[0.0] * order for i in range(order)]
    U = [[0.0] * order for i in range(order)]

    P = pivot_matrix(order)
    PM = mult_matrix(P, matrix)

    for j in range(order):
        L[j][j] = 1.0

        for i in range(j + 1):
            sum_u = sum(U[k][j] * L[i][k] for k in range(i))
            U[i][j] = PM[i][j] - sum_u

        for i in range(j, order):
            sum_l = sum(U[k][j] * L[i][k] for k in range(j))
            L[i][j] = (PM[i][j] - sum_l) / U[j][j]

    y: List = sistemaTriangularInferior(len(L), L, vector) 
    x: List = sistemaTriangularSuperior(len(U), U, y)

    x = [round(z, 4) for z in x]

    return x

def isPositiveDefinite(my_matrix: Matrix, my_order: int):
    matrix = copy.deepcopy(my_matrix)
    order = my_order 

    for i in range(order):
        for j in range(i + 1, order):
            if my_matrix[i][j] != my_matrix[j][i]:
                print(f"Matriz nao e simetrica em ({i},{j})")
                return False

    for k in range(1, order + 1):
        submatrix = leadingPrincipalSubmatrix(matrix, k)
        det = calculoDeterminante(k, submatrix)
        if det <= 0:
            return False
    return True

def cholesky(my_order: int, my_matrix: Matrix, my_vector: List) -> List:
    os.system('clear')
    order = my_order
    matrix = copy.deepcopy(my_matrix)
    vector = copy.deepcopy(my_vector)

    if not isPositiveDefinite(matrix, order):
        print("Matriz nao e definida positiva")
        return []

    L = [[0.0] * order for i in range(order)]
    LT = [[0.0] * order for i in range(order)]

    for i in range(order):
        for k in range(i+1):
            sum_l = sum(L[i][j] * L[k][j] for j in range(k))

            if i == k:
                L[i][k] = sqrt(matrix[i][i] - sum_l)
            else:
                L[i][k] = (1.0 / L[k][k] * (matrix[i][k] - sum_l))

            
    LT = [list(row) for row in zip(*L)]

    y = sistemaTriangularInferior(len(L), L, vector) 
    x = sistemaTriangularSuperior(len(LT), LT, y)

    x = [round(z, 4) for z in x]

    return x

def gaussCompacto(my_order: int, my_matrix: Matrix, my_vector: List) -> List:
    os.system('clear')
    order = my_order
    matrix = copy.deepcopy(my_matrix)
    vector = copy.deepcopy(my_vector)

    if order <= 0:
        print("Ordem invalida")
        return []

    if not isNonSingular(matrix, order):
        print("A matriz possui uma submatriz singular! Tente outra matriz")
        return []

    solution_vector = [0] * len(vector)

    for k in range(order - 1):
        for i in range(k + 1, order):
            scalar = matrix[i][k] / matrix[k][k]
            for j in range(k, order):
                matrix[i][j] = matrix[i][j] - scalar * matrix[k][j]
            vector[i] = vector[i] - scalar * vector[k]

    solution_vector[order - 1] = vector[order - 1] / matrix[order - 1][order - 1]

    for i in range(order - 2, -1, -1):
        sum_vec = vector[i]
        for j in range(i + 1, order):
            sum_vec = sum_vec - matrix[i][j] * solution_vector[j]
        solution_vector[i] = sum_vec / matrix[i][i]

    solution_vector = [round(x, 4) for x in solution_vector]

    return solution_vector

def gaussJordan(my_order: int, my_matrix: Matrix, my_vector: List) -> List:
    os.system('clear')
    order = my_order
    matrix = copy.deepcopy(my_matrix)
    vector = copy.deepcopy(my_vector)

    if order <= 0:
        print("Ordem invalida")

    if not isNonSingular(matrix, order):
        print("A matriz possui uma submatriz singular")
        return []

    for k in range(order):
        for i in range(k + 1, order):
            if max(abs(matrix[i][k]) for i in range(len(matrix))) > abs(matrix[k][k]):
                for j in range(k, order):
                    matrix[k][j], matrix[i][j] = matrix[i][j], matrix[k][j]
                vector[k], vector[i] = vector[i], vector[k]

        pivot = matrix[k][k]
        for j in range(k, order):
            matrix[k][j] /= pivot
        vector[k] /= pivot

        for i in range(order):
            if i == k or matrix[i][k] == 0: continue
            factor = matrix[i][k]
            for j in range(k, order):
                matrix[i][j] -= factor * matrix[k][j]
            vector[i] -= factor * vector[k]

    vector = [round(x, 4) for x in vector]
        
    return vector

def jacobi_converges_row_criteria(matrix: Matrix) -> bool:
    n = len(matrix)
    iteration_matrix = [[0.0] * n for _ in range(n)]

    for i in range(n):
        if matrix[i][i] == 0:
            return False  
        for j in range(n):
            if i != j:
                iteration_matrix[i][j] = -matrix[i][j] / matrix[i][i]

    converges = max(sum(abs(iteration_matrix[i][j]) for j in range(n)) / abs(matrix[i][i]) for i in range(n))

    return converges < 1

def jacobi_converges_column_criteria(matrix: Matrix) -> bool:
    n = len(matrix)
    iteration_matrix = [[0.0] * n for _ in range(n)]

    for i in range(n):
        if matrix[i][i] == 0:
            return False  
        for j in range(n):
            if i != j:
                iteration_matrix[i][j] = -matrix[i][j] / matrix[i][i]

    converges = max(sum(abs(iteration_matrix[j][i]) for j in range(n)) / abs(matrix[i][i]) for i in range(n))

    return converges < 1

def uniform_norm(vector: List[float]) -> float:
    return max(abs(x) for x in vector)

def jacobi(my_order: int, my_matrix: Matrix, my_vector: List, my_vectorAprox: List, my_precision: float, my_max_it: int):
    os.system('clear')
    order = my_order
    matrix = copy.deepcopy(my_matrix)
    vector = copy.deepcopy(my_vector)
    vectorAprox = copy.deepcopy(my_vectorAprox)
    precision = my_precision
    max_it =  my_max_it

    if not jacobi_converges_row_criteria(matrix) and not jacobi_converges_column_criteria(matrix):
        print("Matriz não converge! Tente outra matriz")
        return

    newVec = vectorAprox[:]

    for k in range(max_it):
        for i in range(order):        
            d = vector[i]                  
             
            for j in range(order):     
                if(i != j):
                    d -= matrix[i][j] * vectorAprox[j]

            newVec[i] = d / matrix[i][i]

        solution_norm = uniform_norm(newVec)

        if k > 0:
            relative_change = uniform_norm([newVec[i] - vectorAprox[i] for i in range(order)]) / solution_norm
            if relative_change < precision:
                newVec = [round(x, 4) for x in newVec]
                return newVec

        vectorAprox = newVec[:]

    vectorAprox = [round(x, 4) for x in vectorAprox]

    return vectorAprox

def gauss_seidel_converges(matrix: Matrix) -> bool:

    n = len(matrix)
    iteration_matrix = [[0.0] * n for _ in range(n)]

    for i in range(n):
        if matrix[i][i] == 0:
            return False  
        for j in range(n):
            if i != j:
                iteration_matrix[i][j] = -matrix[i][j] / matrix[i][i]

    converges = max(sum(abs(iteration_matrix[i][j]) for j in range(n)) / abs(matrix[i][i]) for i in range(n))

    return converges < 1 
    
def gauss_seidel_converges_sassenfeld(matrix: Matrix) -> bool:
    n = len(matrix)
    beta_values = []

    for i in range(n):
        if matrix[i][i] == 0:
            return False  
        beta_i = sum(abs(matrix[i][j]) for j in range(n) if j != i) / abs(matrix[i][i])
        beta_values.append(beta_i)

    max_beta = max(beta_values)
    return max_beta < 1

def gaussSeidel(my_order: int, my_matrix: Matrix, my_vector: List, my_vectorAprox: List, my_precision: float, my_max_it: int) -> List:
    os.system('clear')
    order = my_order
    matrix = copy.deepcopy(my_matrix)
    vector = copy.deepcopy(my_vector)
    vectorAprox = copy.deepcopy(my_vectorAprox)
    precision = my_precision
    max_it =  my_max_it

    if not gauss_seidel_converges(matrix) and not gauss_seidel_converges_sassenfeld(matrix):
        print("Matriz não converge! Tente outra matriz")
        return []

    for k in range(max_it):
        for i in range(order):        
            d = vector[i]                  
             
            for j in range(order):     
                if(i != j):
                    d -= matrix[i][j] * vectorAprox[j]

            vectorAprox[i] = d / matrix[i][i]

        solution_norm = uniform_norm(vectorAprox)

        if k > 0:
            relative_change = uniform_norm([vectorAprox[i] - vector[i] for i in range(order)]) / solution_norm
            if relative_change < precision:
                vectorAprox = [round(x, 4) for x in vectorAprox]
                return vectorAprox


    vectorAprox = [round(x, 4) for x in vectorAprox]

    return vectorAprox

def matrizInversa(my_order:int, my_matrix: Matrix):
    os.system('clear')
    order = my_order
    matrix = copy.deepcopy(my_matrix)

    if order < 0:
        print("Ordem invalida")
        return [[]]

    if order == 0:
        print("Matriz nao inversivel ou singular")

    identity_matrix = [[float(i == j) for i in range(order)] for j in range(order)]
    inverse_matrix = [[0] * order for i in range(order)]

    while True:
        try:
            option = int(input("Aperte 1 para inverter usando decomposicao LU ou dois para Gauss compacto.\n"))
            if option in [1, 2]:
                break
            else:
                print("Entrada invalida, tente novamente.\n")
        except ValueError:
            print("ValueError, tente novamente.\n")

    if option == 1:
        L = [[0.0] * order for i in range(order)]
        U = [[0.0] * order for i in range(order)]
    
        P = pivot_matrix(order)
        PM = mult_matrix(P, matrix)
    
        for j in range(order):
            L[j][j] = 1.0
    
            for i in range(j + 1):
                sum_u = sum(U[k][j] * L[i][k] for k in range(i))
                U[i][j] = PM[i][j] - sum_u
    
            for i in range(j, order):
                sum_l = sum(U[k][j] * L[i][k] for k in range(j))
                L[i][j] = (PM[i][j] - sum_l) / U[j][j]
    
        for i in range(order):
            e_i = [identity_matrix[j][i] for j in range(order)]
            y = sistemaTriangularInferior(order, L, e_i)
            x = sistemaTriangularSuperior(order, U, y)
            for j in range(order):
                inverse_matrix[j][i] = x[j]
    
    elif option == 2:
        for i in range(order):
            e_i = [identity_matrix[j][i] for j in range(order)]
            x = gaussCompacto(order, matrix, e_i)
            for j in range(order):
                inverse_matrix[j][i] = x[j]

    for row in inverse_matrix:
        for i in range(len(row)):
            row[i] = round(row[i], 4)
    
    for row in inverse_matrix:
        print('[{}]'.format(', '.join(map(str, row))))

def show_matrix(vector: List, matrix: Matrix):
    os.system('clear')
    print("Vetor dos termos independentes: ")
    print(vector)

    print("Matriz: ")
    for row in matrix:
        print('[{}]'.format(', '.join(map(str, row))))

def get_values():
    os.system('clear')

    while True:
        try:
            order_input = input("Insira a ordem do sistema: ")
            order = int(order_input[0])   
            break
        except ValueError:
            print("Entrada inválida. Insira um número inteiro válido.")

    matrix: Matrix[float] = []
    vector: List[float] = []

    for i in range(order):
        while True:
            try:
                row_input = input(f"Insira os valores para a linha {i+1} (separado por espaco, {order} valores): ")
                row = [float(x) for x in row_input.split()]
                if len(row) != order:
                    print(f"Entrada invalida. Insira valores que sejam igual a ordem {order}.")
                    continue
                matrix.append(row)
                break
            except ValueError:
                print("Entrada invalida. Insira valores numericos")

    while True:
        try:
            vector_input = input(f"Insira os termos independentes (separado por espaco, {order} valores): ")
            vector = [float(x) for x in vector_input.split()]
            if len(vector) != order:
                print(f"Entrada invalida. Insira {order} valores.")
                continue
            break
        except ValueError:
            print("Entrada invalida. Insira valores numericos.")

    return order, matrix, vector

def get_vecProx_prec_maxIt(order):
    os.system('clear')
    while True:
        try:
            vec_input = input(f"Insira o vetor aproximado inicial (separado por espaco, {order} valores): ")
            vec = [float(x) for x in vec_input.split()]
            if len(vec) != order:
                print(f"Entrada invalida. Insira {order} valores.")
                continue
            break
        except ValueError:
            print("Entrada invalida. Insira valores numericos.")

    while True:
        try:
            precision = float(input("Insira a precisao desejada (numero decimal): "))
            break
        except ValueError:
            print("Entrada invalida. Insira um numero decimal.")

    while True:
        try:
            max_iterations = int(input("Insira o numero maximo de iteracoes: "))
            break
        except ValueError:
            print("Entrada invalida. Insira um numero inteiro.")

    return vec, precision, max_iterations

def get_option():
    while True:
        try:
            choice = int(input("Escolha uma opcao (1-13): "))
            if 1 <= choice <= 13:
                return choice
            else:
                print("Opcao invalida. Escolha um valor entre 1 e 12")
        except ValueError:
            print("Opcao invalida. Escolha um valor entre 1 e 12")

def draw_menu(): 
    print("Sistemas Lineares e Matriz Inversa")
    print("Criado por Eduardo Rodrigues Teixeira")
    print("Escolha uma opcao (1-13)")
    print("1  - Calcular determinante")
    print("2  - Sistema Triangular Inferior")
    print("3  - Sistema Triangular Superior")
    print("4  - Decomposicao LU")
    print("5  - Cholesky")
    print("6  - Gauss Compacto")
    print("7  - Gauss Jordan")
    print("8  - Jacobi")
    print("9  - Gauss Seidel")
    print("10 - Matriz Inversa")
    print("11 - Reinserir valores")
    print("12 - Apresentar matriz")
    print("13 - Sair")

def main():
    my_order, my_matrix, my_vector = get_values()

    menu_options = {
        1: calculoDeterminante,
        2: sistemaTriangularInferior,
        3: sistemaTriangularSuperior,
        4: decomposicaoLU,
        5: cholesky,
        6: gaussCompacto,
        7: gaussJordan,
        8: jacobi,
        9: gaussSeidel,
        10: matrizInversa,
        11: get_values,
        12: show_matrix
    }

    while True:
        draw_menu()
        user_choice = get_option()
        
        if user_choice == 13:
            print("Saindo do programa.")
            break
        elif user_choice == 11:
            my_order, my_matrix, my_vector = menu_options[user_choice]()
            input("Aperte qualquer tecla para continuar")
            os.system('clear')
        elif user_choice == 12:
            show_matrix(my_vector, my_matrix)
            input("Aperte qualquer tecla para continuar")
            os.system('clear')
        elif user_choice == 1:
            print(menu_options[user_choice](my_order, my_matrix))
            input("Aperte qualquer tecla para continuar")
            os.system('clear')
        elif user_choice == 10:
            menu_options[user_choice](my_order, my_matrix)
            input("Aperte qualquer tecla para continuar")
            os.system('clear')
        elif user_choice in [2, 3, 4, 5, 6, 7]:
            print(menu_options[user_choice](my_order, my_matrix, my_vector))
            input("Aperte qualquer tecla para continuar")
            os.system('clear')
        elif user_choice in [8, 9]:
            my_vecProx, my_precision, my_maxIt = get_vecProx_prec_maxIt(my_order)
            print(menu_options[user_choice](my_order, my_matrix, my_vector, my_vecProx, my_precision, my_maxIt))
            input("Aperte qualquer tecla para continuar")
            os.system('clear')

if __name__ == "__main__":
    main()
