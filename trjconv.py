#!/bin/python
"""
Script to shrink tajectories by removing frames
"""

from os import sys

num_particles = 0
step_size = 1
inp = "trajectory.xyz"
out = "trjconv.xyz"
for i, val in enumerate(sys.argv):
    if val == "-dt":
        step_size = int(sys.argv[i+1])
    if val == "-f":
        inp = sys.argv[i+1]
    if val == "-o":
        out = sys.argv[i+1]

print "Extracting every %s step from %s" % (step_size, sys.argv[1])

step_counter = 0
in_step = False
outfile = open(out, "w")
with open(inp, "r") as infile:
    for line in infile:

        # get num particles from first line
        if not line.startswith("atom") and num_particles == 0:
            num_particles = int(line.split()[0])
            print "Particles: %s" % num_particles
            print ""

        if not line.startswith("atom"):
            step_counter += 1

        if step_counter == 0 or (step_counter % step_size) == 0:
             outfile.write(line)

        if step_counter % 100 == 0:
            print(step_counter)
            sys.stdout.write("\033[F")

outfile.close()
