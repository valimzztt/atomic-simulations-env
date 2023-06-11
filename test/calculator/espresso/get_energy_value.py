import re
import os 
curr_directory  = os.getcwd()
filename = curr_directory + "/espresso_nonhubbard.pwo"
energies = []
with open(filename,"r") as file:
    for line in file:
        pattern="!"
        if re.search(pattern, line):
            energies.append(re.findall(r'-?\d+', line))
    # Get the last calculated energy value before the calculation converges
    a = energies[0][0]
    b = energies[0][1]
    d = float(f'{a}.{b}')  
    print("The final energies is ", d)
