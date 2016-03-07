#!/usr/bin/python

from numpy import *
from random import gauss

import sys, getopt

if len(sys.argv) != 5 and len(sys.argv) != 2:
  print 'Usage:'
  print '%s -h' % sys.argv[0]
  print '%s -b <type_of_beam> -n <number_particles>' % sys.argv[0]
  sys.exit("Incorrect number of arguments!")

try:
  opts, args = getopt.getopt(sys.argv[1:],"hb:n:",["help","beam=","number="])
except getopt.GetoptError:
  print 'Usage: %s -b <type_of_beam> -n <number_particles>' % sys.argv[0]
  sys.exit("Incorrect input")

beam = ''
Npart = 0

for opt, arg in opts:
  if opt in ("-h", "--help"):
    print 'Generate_InRays_IniBeam_gaussian.py -b <type_of_beam> -n <number_particles>'
    sys.exit()
  elif opt in ("-b", "--beam"):
    if arg != 'Core' and arg != 'Halo' :
      print 'Please type either "Core" or "Halo"!'
      sys.exit("Unkown input")
    else: 
      beam = arg
  elif opt in ("-n", "--number"):
    Npart += int(arg)
    if Npart == 0:
      print 'Please give the number of particles in the beam!'
      sys.exit("Missing input")


print 'Generating inital particle file for ', Npart, ' particles in ', beam
  
gamma=2544.0313111546
ex=5.08e-06
ey=3e-8
dp=0.0008

dx=0
dpx=0
betx=1.348753532
alfx=-0.7862388903
alfy=-3.285248127
bety=19.79462726

#file=open("inrays.madx","w") #Different output for madx files!
myfile=open("iniparticles.dat","w")

sx, spx, sy, spy, se = sqrt(ex*betx/gamma), sqrt(ex/betx/gamma), sqrt(ey*bety/gamma), sqrt(ey/bety/gamma), dp

Rx = 0;
Ry = 0;

j=1;
while (j < Npart+1):
  if beam == 'Core':
    Rx=gauss(0,1) 
    Ry=gauss(0,1)
  elif beam == 'Halo':
    Rx=gauss(0,5)
    Ry=gauss(0,10)
    if sqrt(abs(Rx)*abs(Rx)+abs(Ry)*abs(Ry)) <= 3 :
      continue 
  x1 = random.random() 
  y1 = random.random() 
  x = sx*Rx
  y = sy*Ry
  vpx = random.random()*spx - (alfx*sqrt(ex/betx/gamma))*x1
  vpy = random.random()*spy - (alfy*sqrt(ey/bety/gamma))*y1
  e = 1.3
  z = gauss(0,se)
  var_ener = z*10**3+e*10**3
  j=j+1

  print >>myfile, x*10**6, vpx*10**3, y*10**6, vpy*10**3, z, var_ener;
  #For different output file variant:
  #print >>file, "ptc_start, x=", x, ", px=", vpx, ", y=", y, ", py=", vpy, ", pt=", z, ";"

myfile.close()
#file.close()
