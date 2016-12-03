import kernel

def chi(G_real, G_imag, A, omega_n, omega, C_real_inv, C_imag_inv):
    import numpy as np
    import numpy.linalg

    Niom = len(G_real)
    Nomega = len(omega)

    vector = G_real - kernel.K_matrix_real(omega_n, omega).dot(A)
    result = np.transpose(vector).dot(C_real_inv).dot(vector)

    vector = G_imag - kernel.K_matrix_imag(omega_n, omega).dot(A)
    result = result + np.transpose(vector).dot(C_imag_inv).dot(vector)
    return result
