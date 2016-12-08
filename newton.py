import kernel
import J
import f
import nan

def diff(a, b):
    import numpy as np
    s = 0.0
    for i in range(len(a)):
        s = s + (a[i] - b[i])**2
    return np.sqrt(s)

def newton(alpha, G_real, G_imag, omega_n, omega, A_initial, C_real_inv, C_imag_inv):
    import numpy as np
    import numpy.linalg
    
    Nomega = len(omega)
    Niom = len(omega_n)
    
    iterationMax = 20
    counter = 0
    
    eps = 0.00001
    while(True):
        counter = counter + 1
        if (counter > iterationMax):
            break
        function = f.f(alpha, G_real, G_imag, A_initial, omega_n, omega, C_real_inv, C_imag_inv)
        Jacobian = J.J(alpha, A_initial, omega_n, omega, C_real_inv, C_imag_inv)
        A_updated = A_initial - numpy.linalg.inv(Jacobian).dot(function)
        if (nan.array_isnan(A_updated)):
            A_initial = A_updated
            break
        error = diff(A_initial, A_updated)
        if (error < eps):
            break
        print "counter = ", counter, ", diff = ", error
        A_initial = A_updated
    
    return A_initial
