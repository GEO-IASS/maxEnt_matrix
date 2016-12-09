#!/usr/bin/env python

import default
import kernel
import f
import J
import chi
import entropy
import newton
import printFile
import nan

def read_spectral(fileName):
    import numpy as np
    omega = []
    A = []
    ifile = open(fileName, "r")
    for index, string in enumerate(ifile):
        a = string.split()
        omega.append(float(a[0]))
        A.append(float(a[1]))
    ifile.close()
    return omega, np.asarray(A)

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

    if (len(sys.argv) == 1):
        print "a0 = sys.argv[1], b0 = sys.argv[2]. alpha = a0*exp(-i*b0). "
        return -1
    if (len(sys.argv) == 2):
        a0 = float(sys.argv[1])
        b0 = 0.05
    if (len(sys.argv) == 3):
        a0 = float(sys.argv[1])
        b0 = float(sys.argv[2])
    
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
    model = np.zeros(Nomega)
    for i in range(Nomega):
        model[i] = default.D(omega[i])
    printFile.printFile(omega, model, "model.txt")
    if (not os.path.exists("A_initial.txt")):
        for i in range(len(A_initial)):
            A_initial[i] = default.D(omega[i])
        printFile.printFile(omega, A_initial, "A_initial.txt")
    else:
        omega, A_initial = read_spectral("A_initial.txt")
        Nomega = len(omega)
    
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
    eig = 0.005
    for i in range(Niom):
        C_real[i, i] = eig**2
        C_imag[i, i] = eig**2
    C_real_inv = numpy.linalg.inv(C_real)
    C_imag_inv = numpy.linalg.inv(C_imag)
    printFile.printMatrix(C_real_inv, "C_real_inv.txt")

    if (True):
        alpha = []
        for i in range(30):
            alpha.append(a0*np.exp(-i*b0))
        
        ofile = open("alpha.txt", "a")
        for i in range(len(alpha)):
            A_updated = newton.newton(alpha[i], G_real, G_imag, omega_n, omega, A_initial, C_real_inv, C_imag_inv)
            if (nan.array_isnan(A_updated)):
                omega, A_initial = read_spectral("A_initial.txt")
                continue
            output = "A_updated_alpha_" + str(alpha[i]) + ".txt"
            printFile.printFile(omega, A_updated, output)
            os.system("cp " +  output + " A_initial.txt")
            ofile.write(str(alpha[i]) + "\n")
            print "alpha = ", alpha[i]
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
