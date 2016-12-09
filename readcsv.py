def main():
    import os
    import sys

    if (len(sys.argv) != 2):
        print "fileName = sys.argv[1]. "
        return -1

    fileName = sys.argv[1]
    contents = []
    ifile = open(fileName, "r")
    for index, string in enumerate(ifile):
        string = string.replace(",", " ")
        contents.append(string.strip())
    ifile.close()

    ofile = open(fileName, "w")
    for i in range(len(contents)):
        ofile.write(contents[i] + "\n")
    ofile.close()
    return 

main()
