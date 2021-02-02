# Imports 
from ROOT import TFile
import numpy as np
import matplotlib.pyplot as plt
from scipy import optimize

path="/home/aa/EPPT/run/"

load_file = TFile(path+"B5.root", "READ")
tree = load_file.Get('B5')

#Drift chamber 1
#Event 0 to initiate numpy array with the good shape and just loop on all other events after
tree.GetEntry(0)
dc1_xpos = np.array( getattr(tree, 'Dc1HitsVector_x') )
dc1_ypos = np.array( getattr(tree, 'Dc1HitsVector_y') )
dc1_zpos = np.array( getattr(tree, 'Dc1HitsVector_z') )

for entry in range(1, tree.GetEntries()):
    tree.GetEntry(entry)
    if (np.array( getattr(tree, 'Dc1HitsVector_x') ).size >= 5):
       dc1_xpos=np.vstack( ( dc1_xpos , np.array( getattr( tree , 'Dc1HitsVector_x' ) )[:5] ) )
    else : 
        print("Too less x hits. Reject event #", entry)
        
    if (np.array( getattr(tree, 'Dc1HitsVector_y') ).size >= 5):
       dc1_ypos=np.vstack( ( dc1_ypos , np.array( getattr( tree , 'Dc1HitsVector_y' ) )[:5] ) )
    else : 
        print("Too less y hits. Reject event #", entry)
    
    if (np.array( getattr(tree, 'Dc1HitsVector_z') ).size >= 5):
       dc1_zpos=np.vstack( ( dc1_zpos , np.array( getattr( tree , 'Dc1HitsVector_z' ) )[:5] ) )
    else : 
        print("Too less z hits. Reject event #", entry)

#Not sure how to deal with events with > 5 hits so for now just take the first 5 hits. Ignore events with less than 5 hits. 

#Same for DC2

dc2_xpos = np.array( getattr(tree, 'Dc2HitsVector_x') )
dc2_ypos = np.array( getattr(tree, 'Dc2HitsVector_y') )
dc2_zpos = np.array( getattr(tree, 'Dc2HitsVector_z') )

for entry in range(1, tree.GetEntries()):
    tree.GetEntry(entry)
    if (np.array( getattr(tree, 'Dc2HitsVector_x') ).size >= 5):
       dc2_xpos=np.vstack( ( dc2_xpos , np.array( getattr( tree , 'Dc2HitsVector_x' ) )[:5] ) )
    else : 
        print("Too less x hits. Reject event #", entry)
        
    if (np.array( getattr(tree, 'Dc2HitsVector_y') ).size >= 5):
       dc2_ypos=np.vstack( ( dc2_ypos , np.array( getattr( tree , 'Dc2HitsVector_y' ) )[:5] ) )
    else : 
        print("Too less y hits. Reject event #", entry)
    
    if (np.array( getattr(tree, 'Dc2HitsVector_z') ).size >= 5):
       dc2_zpos=np.vstack( ( dc2_zpos , np.array( getattr( tree , 'Dc2HitsVector_z' ) )[:5] ) )
    else : 
        print("Too less z hits. Reject event #", entry)

save = False

def linear(x,a,b):
    return a*x+b

zPosDC1 = np.array([-6.25, -5.75, -5.25, -4.75, -4.25])
zPosDC2 = np.array([2.25, 2.75, 3.25, 3.75, 4.25])

#Set x values to the good units
dc1_xpos_metres=dc1_xpos/1000
dc2_xpos_metres=dc2_xpos/1000


popt1,pcov1 = optimize.curve_fit(linear,zPosDC1,dc1_xpos_metres[event])
z1=np.linspace(-6.25,0)
x1=linear(x=z1,a=popt1[0],b=popt1[1])

popt2,pcov2 = optimize.curve_fit(linear,zPosDC2,dc2_xpos_metres[event])
z2=np.linspace(-0.5,4.25)
x2=linear(x=z2,a=popt2[0],b=popt2[1])


#PLOT
fig = plt.figure(figsize=(8,8))

axes= fig.add_axes([0.1,0.1,0.8,0.8])

axes.scatter(zPosDC1,dc1_xpos_metres[event],label='DC1')
axes.scatter(zPosDC2, dc2_xpos_metres[event] , label='DC2')
axes.plot(z1,x1)
axes.plot(z2,x2)
axes.set_ylim(-1.5E-2,1E-3)

square=plt.Rectangle((-1,-1),2,2,fill=False)
plt.gca().add_patch(square)

plt.legend()

plt.grid()
if save : plt.savefig("event"+str(event)+".png")
    
plt.title("Event #"+str(event))

plt.show()
