#!/bin/python
"""
Script to shrink tajectories by removing frames
"""

from os import sys

num_particles = 0
step_size = 1
inp = "trajectory.xyz"
out = "trjconv.xyz"
skip = -1
move = 0
max_z = 0
for i, val in enumerate(sys.argv):
    if val == "-dt":
        step_size = int(sys.argv[i+1])
    if val == "-f":
        inp = sys.argv[i+1]
    if val == "-o":
        out = sys.argv[i+1]
    if val == "-s" or val == "--skip":
        skip = int(sys.argv[i+1])
    if val == "-m" or val == "--move":
        move = float(sys.argv[i+1])

print "Extracting every %s step from %s" % (step_size, sys.argv[1])

step_counter = 0
in_step = False
outfile = open(out, "w")

if skip != -1:
    print "Skipping %s frames" % skip


max_z = 0
with open(inp, "r") as infile:
    for line in infile:

        # get num particles from first line
        if not line.startswith("atom") and num_particles == 0:
            num_particles = int(line.split()[0])
            max_z = float(line.split()[5])
            print "Particles: %s" % num_particles
            print ""

        if not line.startswith("atom"):
            step_counter += 1

        if step_counter > skip and (step_counter == 0 or (step_counter % step_size) == 0):
             if move > 0 and line.startswith("atom"):
                split = line.split()
                split[3] = float(split[3]) + move
                if(split[3] > max_z):
                    split[3] -= max_z
                if(split[3] < 0):
                    split[3] += max_z
                outfile.write("%s %s %s %s\n" % (split[0],split[1],split[2],split[3]))
             else:
                outfile.write(line)

        if step_counter % 100 == 0:
            print(step_counter)
            sys.stdout.write("\033[F")

outfile.close()
