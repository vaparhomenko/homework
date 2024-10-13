import time

def fib_recursive(n):
    if n <= 0:
        return 0
    elif n == 1:
        return 1
    else:
        return fib_recursive(n - 1) + fib_recursive(n - 2)

n = 15
start_time = time.time()
print(f"The {n}-th Fibonacci number by recursive method: {fib_recursive(n)}")
print(f"Recursive method execution time: {time.time() - start_time:.8f} seconds")


def fib_iterative(n):
    if n <= 0:
        return 0
    elif n == 1:
        return 1

    a, b = 0, 1
    for _ in range(2, n + 1):
        a, b = b, a + b

    return b


n = 500
start_time = time.time()
print(f"The {n}-th Fibonacci number by iterative method: {fib_iterative(n)}")
print(f"Iterative method execution time: {time.time() - start_time:.8f} seconds")

def multiply_matrices(A, B):
    return [
        [A[0][0] * B[0][0] + A[0][1] * B[1][0], A[0][0] * B[0][1] + A[0][1] * B[1][1]],
        [A[1][0] * B[0][0] + A[1][1] * B[1][0], A[1][0] * B[0][1] + A[1][1] * B[1][1]]
    ]

def matrix_power(M, n):
    if n == 1:
        return M
    elif n % 2 == 0:
        half_power = matrix_power(M, n // 2)
        return multiply_matrices(half_power, half_power)
    else:
        return multiply_matrices(M, matrix_power(M, n - 1))


def fib_matrix(n):
    if n == 0:
        return 0
    elif n == 1:
        return 1

    F = [[1, 1],
         [1, 0]]

    result = matrix_power(F, n)

    return result[0][1]  # F_n находится в позиции [0][1]

n = 500
start_time = time.time()
print(f"The {n}-th Fibonacci number by matrix method: {fib_matrix(n)}")
print(f"Matrix method execution time: {time.time() - start_time:.8f} seconds")


