#!/usr/bin/env python
# coding: utf-8

import sys
from rotation import *


def usage():
    print "usage: %s -m|-3|-4" % sys.argv[0]
    print "  -m\tmdview."
    print "  -3\tnx3a."
    print "  -4\tnx4a."
    sys.exit(1)

if len(sys.argv) == 2:
    mode = sys.argv[1]
else:
    usage()


def output_mdview(atoms):
    print "-center 0 0 0"
    print "-fold"
    print len(atoms)
    for label,x,y,z in atoms:
        print label,x,y,z



defr = dict()
while True:
    #read a line, anyway
    line = sys.stdin.readline()
    if len(line) == 0:
        break
    columns = line.split()
    #look up tags
    tag = columns[0]
    if tag in  ('@ID08',):
        line = sys.stdin.readline()
        id08 = line[0:8]
        if mode in ("-3","-4"):
            print "@ID08"
            print id08
    elif tag in ('@BOX3',):
        line = sys.stdin.readline()
        columns = line.split()
        box = map(float,columns[0:3])
        if mode == "-m":
            print "-length '(%s, %s, %s)'" % (box[0],box[1],box[2])
        else:
            print "@BOX3"
            print box[0],box[1],box[2]
    elif tag in  ('@DEFR',):
        line = sys.stdin.readline()
        id = line[0:8]
        line = sys.stdin.readline()
        nsite = int(line)
        sites = []
        intr  = []
        for site in range(nsite):
            line = sys.stdin.readline()
            columns = line.split()  # x,y,z,mass,label
            columns[0:4] = map(float,columns[0:4])
            sites.append(columns)
        for site in range(nsite):
            line = sys.stdin.readline()
            columns = line.split() #eps, sig, charge
            columns = map(float,columns[0:3])
            intr.append(columns)
        defr[id] = (sites,intr)
    elif tag in  ('@WTG6', '@WTG3', '@NX4A', '@NX3A'):
        #get the first line == number of molecules
        line = sys.stdin.readline()
        columns = line.split()
        nmol = int(columns[0])

        if mode == "-3":
            print "@NX3A"
        elif mode == "-4":
            print "@NX4A"
        if mode in ("-3", "-4"):
            print nmol
        #read molecular info
        atoms = []
        for i in range(nmol):
            line = sys.stdin.readline()
            columns = line.split()
            cx,cy,cz = map(float,columns[0:3])
            #last data, i-th molecule, 6 velocities
            if tag == '@WTG6':
                rotmat = map(float,columns[3:12])
            elif tag in ('@WTG3', '@NX4A'):
                quat = map(float,columns[3:7])
                rotmat = quat2rotmat(quat)
            elif tag in ('@NX3A'):
                euler = map(float,columns[3:7])
                quat = euler2quat(euler)
                rotmat = quat2rotmat(quat)
            if mode == "-m":
                mol = defr[id08]
                for site in mol[0]:
                    x,y,z,mass,label = site
                    xx = cx + x*rotmat[0] + y*rotmat[1] + z*rotmat[2]
                    yy = cy + x*rotmat[3] + y*rotmat[4] + z*rotmat[5]
                    zz = cz + x*rotmat[6] + y*rotmat[7] + z*rotmat[8]
                    atoms.append((label, xx,yy,zz))
            elif mode == "-3":
                euler = quat2euler(rotmat2quat(rotmat))
                print cx,cy,cz,euler[0],euler[1],euler[2]
            elif mode == "-4":
                quat = rotmat2quat(rotmat)
                print cx,cy,cz,quat[0],quat[1],quat[2],quat[3]
        if mode == "-m":
            output_mdview(atoms)
