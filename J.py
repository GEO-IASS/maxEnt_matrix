import default
import kernel

def J(alpha, A, omega_n, omega, C_real_inv, C_imag_inv):
    import numpy as np
    import numpy.linalg

    matrix = kernel.K_matrix_real(omega_n, omega)
    result = -np.transpose(matrix).dot(C_real_inv).dot(matrix)

    matrix = kernel.K_matrix_imag(omega_n, omega)
    result = result - np.transpose(matrix).dot(C_imag_inv).dot(matrix)

    Nomega = len(omega)
    for i in range(Nomega):
        result[i, i] = result[i, i] - alpha/A[i]
    return result
