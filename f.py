import default
import kernel

def f(alpha, G_real, G_imag, A, omega_n, omega, C_real_inv, C_imag_inv):
    import numpy as np
    import numpy.linalg 

    domega = omega[1] - omega[0]
    Nomega = len(omega)
    result = np.zeros(Nomega)
    for nw in range(Nomega):
        result[nw] = -alpha*domega*(1 + np.log(A[nw]/default.D(omega[nw])))
    
    matrix = kernel.K_matrix_real(omega_n, omega)
    vector_right = G_real - matrix.dot(A)
    result = result + np.transpose(matrix).dot(C_real_inv).dot(vector_right)

    matrix = kernel.K_matrix_imag(omega_n, omega)
    vector_right = G_imag - matrix.dot(A)
    result = result + np.transpose(matrix).dot(C_imag_inv).dot(vector_right)

    return result
