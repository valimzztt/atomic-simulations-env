# import required module
import os
import json
# assign directory
directory = r"C:\Users\Utente\Desktop\QMI Internship Summer 2023\HEO_calc\TiO2-hubbard-refinement"

data_dict = {}
temperatures =[]
info = []
for filename in os.listdir(directory):
    f = os.path.join(directory, filename)
    if os.path.isfile(f):
        if(filename.startswith("MC")):
            with open(f) as f:
                read = f.read()
                content = read.replace("'", "\"")
                start_index = content.find("{")
                end_index = content.find("}") + 1
                json_string = content[start_index:end_index] 
                # Convert the string to a JSON object
                mc_data = json.loads(json_string)
                print(mc_data)
                temperature = mc_data["temperature"]
                energy = mc_data["energy"]
                heat_capacity = mc_data["heat_capacity"]
                temperatures.append(temperature)
                info.append((energy, heat_capacity))
                
data = dict(zip(temperatures, info))
print(data)

import matplotlib.pyplot as plt

# Extract keys and values from the dictionary
temperatures = list(data.keys())
info = list(data.values())
energies = []
heat_capacities = []
for datapoint in info:
    energies.append(datapoint[0])
    heat_capacities.append(datapoint[1])


fig, ax = plt.subplots(2, 2)
fig.tight_layout()
plt.subplot(2, 1, 1) # row 1, col 2 index 1
plt.plot(temperatures, energies, "bo",  marker='.')
plt.xlabel('Temperature (K)')
plt.ylabel('Energy (eV)')
plt.title('Energy and Heat Capacity versus Temperature for MC simulation on (Ti, Zr)O2')

plt.subplot(2, 1, 2) # index 2
plt.plot(temperatures, heat_capacities, "bo",  marker='.')
plt.xlabel('Temperature (K)')
plt.ylabel('Heat capacity (J*K)')



# Set l

# Show the plot
plt.show()
              

                


                