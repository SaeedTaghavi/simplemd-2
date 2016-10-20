#!/usr/bin/env python

from __future__ import print_function
import os

headers = []
sources = []
for (dirpath, dirnames, filenames) in os.walk("."):
    for file in filenames:
        ext = file[len(file)-4:len(file)]
        if ext == ".hpp":
            headers.append( (dirpath, file ) )
        elif ext == ".cpp":
            sources.append( (dirpath, file ) )


depend = dict()
for path,file in headers+sources:
    target = "/".join((path,file))
    for line in open(target):
        if len(line) > 8:
            if line[0:8] == "#include":
                h = line[8:].split('"')
                if len(h) > 1:
                    #print(target,"uses",h[1])
                    if target not in depend:
                        depend[target] = set()
                    depend[target].add(h[1])

for target in depend:
    print("{0}: ".format(target), " ".join(depend[target]))

