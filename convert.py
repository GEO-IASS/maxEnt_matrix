def main():
    import sys
    import os
    import numpy as np

    if (len(sys.argv) != 2):
        print "fileName = sys.argv[1]. "
        return -1

    fileName = sys.argv[1]

    x = []
    y = []

    ifile = open(fileName, "r")
    for i, string in enumerate(ifile):
        a = string.split()
        x.append(float(a[0]))
        y.append(float(a[1]))
    ifile.close()

    ofile = open(fileName.replace(".txt", "_2pi.txt"), "w")
    for i in range(len(x)):
        ofile.write(str(x[i]) + "    " + str(y[i]/(2*np.pi)) + "\n")
    ofile.close()

    return 0

main()
