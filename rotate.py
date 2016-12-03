import default
import kernel
import f
import J
import chi
import entropy

def readFiles(Greal, Gimag):
    import os
    import sys
    import numpy as np
    import sys

    omega_n = []
    G = []
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

    for i in range(len(G_real)):
        G.append(G_real[i] + G_imag[i]*1j)

    return omega_n, G


def printFile(x, y, fileName):
    ofile = open(fileName, "w")
    for i in range(len(x)):
        ofile.write(str(x[i]) + "    " + str(y[i]) + "\n")
    ofile.close()

def isOrthogonal(matrix):
    import sys
    import numpy as np
    t = np.shape(matrix)
    if (t[0] != t[1]):
        sys.exit("Matrix is not square. ")
    dimension = t[0]
    result = np.transpose(matrix).dot(matrix)
    identity = np.zeros((dimension, dimension))
    for i in range(dimension):
        identity[i,i] = 1.0
    diff = result - identity
    norm = np.sqrt(np.trace(np.transpose(diff).dot(diff)))
    return norm == 0
    

def main():
    import os
    import sys
    import numpy as np
    import numpy.linalg

    Greal = "G_cc_real.txt"
    Gimag = "G_cc_imag.txt"
    omega_n, G = readFiles(Greal, Gimag)
    Niom = len(omega_n)


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

    eigenvalues_real, eigenvectors_real = np.linalg.eig(C_real)
    eigenvalues_imag, eigenvectors_imag = np.linalg.eig(C_imag)

    O_real = numpy.zeros((Niom, Niom))
    O_imag = numpy.zeros((Niom, Niom))
    for i in range(Niom):
        for j in range(Niom):
            O_real[i, j] = eigenvectors_real[i][j]
            O_imag[i, j] = eigenvectors_imag[i][j]

    if (not isOrthogonal(O_real)):
        print "O_real is not an orthogonal matrix. "
        return -1
    if (not isOrthogonal(O_imag)):
        print "O_imag is not an orthogonal matrix. "
        return -1



    C_real_inv = numpy.linalg.inv(C_real)
    C_imag_inv = numpy.linalg.inv(C_imag)



    return 0

main()
