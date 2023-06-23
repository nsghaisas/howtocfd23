## This program will plot charge density

import matplotlib.pyplot as plt
import numpy as np

charge=np.loadtxt('output/output0.txt',unpack=True,skiprows=1)
charge1=np.loadtxt('output/output40.txt',unpack=True,skiprows=1)
charge2=np.loadtxt('output/output80.txt',unpack=True,skiprows=1)


plt.title('Charge density',fontweight='bold',size=15)
plt.xlabel('Grid position',fontweight='bold',size=15)
plt.ylabel('Charge density',fontweight='bold',size=15)
#plt.yticks(y1)
plt.plot(charge[0], charge[2],'--or',label='Charge density')
plt.legend()
#plt.show()
#plt.plot(charge[0], charge[3],'--or',label='E')
#plt.show()
plt.plot(charge1[0], charge1[2],'--ob',label='Charge density')
plt.legend()
#plt.show()
#plt.plot(charge1[0], charge1[3],'--or',label='E')
#plt.show()
plt.plot(charge2[0], charge2[2],'--og',label='Charge density')
plt.legend()
#plt.show()
#plt.plot(charge2[0], charge2[3],'--or',label='E')
plt.show()


