#!/usr/bin/env python

##############
# Line count #
##############

def blocks(files, size=65536):
    while True:
        b = files.read(size)
        if not b: break
        yield b

def getFileLineCount(filename, extensionType = "uncompressed"):

    import gzip, bz2
    
    if extensionType == "gz":
        f = gzip.open(filename, "r")
        return sum(bl.count(b"\n") for bl in blocks(f))
    elif extensionType == "bz2":
        f = bz2.open(filename, "r")
        return sum(bl.count(b"\n") for bl in blocks(f))
    else:
        f = open(filename, "r")
        return sum(bl.count("\n") for bl in blocks(f))
