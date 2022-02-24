#!/usr/bin/python3
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
#                                                                   #
#   Calculate phosphorescent lifetimes according to                 #
#   K. Mori,  T. P. M. Goumans, E. van Lentheb  and  F. Wang        #
#   Phys. Chem. Chem. Phys., 2014,16, 14523-14530                   #
#                                                                   #
#   Using ORCA SOC-calculation as explained at                      #
#   https://www.orcasoftware.de/tutorials_orca/spec/SOC.html        #
#                                                                   #
#   Thanks to Bernardo De Souza @ The ORCA forum                    #
#                                                                   #
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
#

import sys, getopt, math

epsilon =  8.8541878128e-12
planck = 1.054571817e-34
me = 9.1093837015e-31
e = 1.602176634e-19
alpha = 7.2973525693e-3
cm2Eh = 4.55634e-6
kb = 1.380649e-23
Eh2Joule = 4.35974e-18
t0 = ((4*math.pi*epsilon)**2)*(planck**3)/(me*(e**4))
T = 300
states = 3


E = []
t = []
k = []
w = []

def printHelp():
    print ('test.py -i <inputfile> -t <300> -s <3>')
    sys.exit(2)

def getK(E,t):
    return 4.0/3.0/t0 *alpha**3*E**3*t

def printArray(array):
    for i in array:
        print(i)

def calculateTime(k, E):
    global T
    w = []
    tau_A = 0
    tau_B = 0
    t_A = 0
    t_B = 0
    for curr_k, curr_E in zip(k, E):
        w = math.exp(-1*(curr_E-E[0])*Eh2Joule/kb/T)
        tau_A = tau_A + w
        tau_B = tau_B + w*curr_k
        t_A = t_A + 1
        t_B = t_B + curr_k
#        print(tau_A, tau_B, t_A, t_B)

    if(tau_B != 0):
        tau = tau_A/tau_B
        print(len(k),'\t',"%.4e" %curr_k,'\t',  "%.2e" %(t_A/t_B*1e6),'\t',"%.2e" %(tau_A/tau_B*1e6))
        #print(tau* 1e6)
    else:
        print("Division by zero in case of state", len(k), " -- skipping this --")

def main(argv):
   inputfile = ''
   calculate = 0
   count = 0
   state = 0
   global T
   global states
   try:
      opts, args = getopt.getopt(argv,"hi:T:s:")
   except getopt.GetoptError:
      printHelp()

   for opt, arg in opts:
      if opt == '-h':
        printHelp()
      elif opt in ("-i"):
         inputfile = arg
      elif opt in ("-T"):
         T = float(arg)
      elif opt in ("-s"):
         states = float(arg)

   if(len(inputfile) == 0):
        printHelp()

   print('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
   print('~ Calculation of phosphorescent lifetimes according to Phys. Chem. Chem. Phys., 2014,16, 14523-14530 from SOC calculations using ORCA ~')
   print('~                        Please cite the above article and the appropriate ORCA methods                                               ~')
   print('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n')
   print('Using ',inputfile, ' and ',"%.0f"%states, ' states at a tempature of ', T, ' K\n\n')


   with open(inputfile, 'r') as file_in:
       lines = []
       for line in file_in:
        elements = line.split();

        if(len(elements) == 10):
            if((elements[6] == "TRANSITION") and (elements[7] == "ELECTRIC") and (elements[8] == "DIPOLE")):
                calculate = 1
                print('State\t     k\t\t𝜏(average)\t𝜏(Boltzman)')
                print('     \t   [1/s]\t   [μs]\t\t    [μs]\n')

        if(calculate == 1):
            count += 1

        if(count > 5 and  state < states):
                curr_E = float(elements[2])*cm2Eh
                curr_t = float(elements[5])
                E.append(curr_E)
                t.append(curr_t)
                k.append(getK(curr_E,curr_t))
                calculateTime(k, E)
                state += 1

if __name__ == "__main__":
   main(sys.argv[1:])
