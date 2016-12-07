import default
import kernel
import f
import J
import chi
import entropy
import newton
import printFile
import nan

def readFiles(Greal, Gimag):
    import os
    import sys
    import numpy as np
    import sys

    omega_n = []
    G_real = []
    G_imag = []

    try:
        ifile = open(Greal, "r")
    except:
        sys.exit(Greal + " does not exist. ")
    for index, string in enumerate(ifile):
        a = string.split()
        omega_n.append(float(a[0]))
        G_real.append(float(a[1]))
    ifile.close()

    try:
        ifile = open(Gimag, "r")
    except:
        sys.exit(Gimag + " does not exist. ")
    for index, string in enumerate(ifile):
        a = string.split()
        G_imag.append(float(a[1]))
    ifile.close()

    G_real = np.asarray(G_real)
    G_imag = np.asarray(G_imag)

    return omega_n, G_real, G_imag


def main():
    import os
    import sys
    import numpy as np
    import numpy.linalg

    Greal = "G_cc_real.txt"
    Gimag = "G_cc_imag.txt"
    omega_n, G_real, G_imag = readFiles(Greal, Gimag)
    Niom = len(omega_n)

    N = 90
    omega_lower = -5
    omega_upper = -omega_lower
    domega = (omega_upper - omega_lower)/float(N)
    omega = np.zeros(N+1)
    for i in range(N+1):
        omega[i] = omega_lower + i*domega
    Nomega = len(omega)
    A_initial = np.zeros(Nomega)
    if (not os.path.exists("A_initial.txt")):
        for i in range(len(A_initial)):
            A_initial[i] = default.D(omega[i])
        printFile.printFile(omega, A_initial, "A_initial.txt")
    else:
        omega = []
        A_initial = []
        ifile = open("A_initial.txt", "r")
        for index, string in enumerate(ifile):
            a = string.split()
            omega.append(float(a[0]))
            A_initial.append(float(a[1]))
        ifile.close()
        Nomega = len(omega)
        A_initial = np.asarray(A_initial)

    C_real = np.zeros((Niom, Niom))
    C_imag = np.zeros((Niom, Niom))

    ifile = open("CM_cc_real.txt", "r")
    for (index, string) in enumerate(ifile):
        a = string.split()
        rowIndex = int(a[0])-1
        colIndex = int(a[1])-1
        if (True):
            C_real[rowIndex, colIndex] = float(a[2])
    ifile.close()
    ifile = open("CM_cc_imag.txt", "r")
    for (index, string) in enumerate(ifile):
        a = string.split()
        rowIndex = int(a[0])-1
        colIndex = int(a[1])-1
        if (True):
            C_imag[rowIndex, colIndex] = float(a[2])
    ifile.close()           
    for i in range(Niom):
        C_real[i, i] = 0.01**2
        C_imag[i, i] = 0.01**2
    C_real_inv = numpy.linalg.inv(C_real)
    C_imag_inv = numpy.linalg.inv(C_imag)
    printFile.printMatrix(C_real_inv, "C_real_inv.txt")

    if (False):
        alpha = 1.0
        function = f.f(alpha, G_real, G_imag, A_initial, omega_n, omega, C_real_inv, C_imag_inv)
        Jacobian = J.J(alpha, A_initial, omega_n, omega, C_real_inv, C_imag_inv)
        print function
        printFile.printMatrix(Jacobian, "J.txt")
        return 0

    if (True):
        alpha = []
        for i in range(80):
            alpha.append(0.012*np.exp(-i*0.02))

        ofile = open("alpha.txt", "a")
        for i in range(len(alpha)):
            print "alpha = ", alpha[i]
            A_updated = newton.newton(alpha[i], G_real, G_imag, omega_n, omega, A_initial, C_real_inv, C_imag_inv)
            if (nan.array_isnan(A_updated)):
                break
            output = "A_updated_alpha_" + str(alpha[i]) + ".txt"
            printFile.printFile(omega, A_updated, output)
            os.system("cp " +  output + " A_initial.txt")
            ofile.write(str(alpha[i]) + "\n")
        ofile.close()
    else:
        alpha = 0.01
        print "alpha = ", alpha
        A_updated = newton(alpha, G, omega_n, omega, A_initial, C_real_inv, C_imag_inv)
        if (nan.array_isnan(A_updated)):
            printFile.printFile(omega, A_updated, "A_updated_alpha_" + str(alpha) + ".txt")
            os.system("cp A_updated_alpha_" + str(alpha) + ".txt" + "A_initial.txt")
    return 0

main()
