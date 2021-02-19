#Imports
from ROOT import TFile
import numpy as np
import matplotlib.pyplot as plt
from scipy import optimize , stats
import sys

if str(sys.argv[2]) == "momentumMode":
    momentumMode = True
    calorimetryMode = False
elif str(sys.argv[2]) == "calorimetryMode":
    momentumMode = False
    calorimetryMode = True
#Booleans to monitor 
debug = False


#Useful functions 
def linear(x,a,b):
    return a*x+b

def compute_p(event,B,L,plot=False,save=False):

    popt1,pcov1 = optimize.curve_fit(linear,zPosDC1,dc1_xpos_metres[event])
    z1=np.linspace(-6.25,0,626)
    x1=linear(x=z1,a=popt1[0],b=popt1[1])

    popt2,pcov2 = optimize.curve_fit(linear,zPosDC2,dc2_xpos_metres[event])
    z2=np.linspace(-0.5,4.25,476)
    x2=linear(x=z2,a=popt2[0],b=popt2[1])

    #Intersection point 
    z_intersect = (popt2[1] - popt1[1]) / (popt1[0] - popt2[0])
    x_intersect = popt1[0] * z_intersect + popt1[1]

    h = abs(zPosDC2[0] - z_intersect)
    X = abs(dc2_xpos_metres[event][0] - x_intersect)
    theta = np.arcsin(X/h)
    p = (0.3*B*L) / (2 * np.sin(theta/2) )

    if plot : 
    
        fig = plt.figure(figsize=(8,8))

        axes= fig.add_axes([0.1,0.1,0.8,0.8])

        axes.scatter(zPosDC1,dc1_xpos_metres[event],label='DC1')
        axes.scatter(zPosDC2, dc2_xpos_metres[event] , label='DC2')
        axes.plot(z1,x1)
        axes.plot(z2,x2)
        axes.scatter(z_intersect,x_intersect,color='black')
        axes.set_ylim(-1.5E-2,1E-3)

        square=plt.Rectangle((-1,-1),2,2,fill=False)
        plt.gca().add_patch(square)



        plt.legend()

        plt.grid()
        if save : plt.savefig("fig/event"+str(event)+".png")
    
        plt.title("Event #"+str(event))

        plt.show()
    
    return(p)


#File loading and set some constants
path="/home/aa/EPPT/run/"

case = str(sys.argv[1])

print("Case : ", case)

if (case == 'case1') or (case == 'case1_PbBefore') or (case == 'case1_PbAfter') or (case == 'case1_positrons') or (case == 'case1_protons'):
    P_TRUE = 100
    B_TRUE = 0.5
elif (case == 'case2') or (case == 'case2_PbBefore') or (case == 'case2_PbAfter') :
    P_TRUE = 100
    B_TRUE = 0.25
elif (case == 'case3') or (case == 'case3_PbBefore') or (case == 'case3_PbAfter') :
    P_TRUE = 100
    B_TRUE = 1
elif (case == 'case4') or (case == 'case4_PbBefore') or (case == 'case4_PbAfter') :
    P_TRUE = 50
    B_TRUE = 0.5
elif (case == 'case5') or (case == 'case5_PbBefore') or (case == 'case5_PbAfter') :
    P_TRUE = 200
    B_TRUE = 0.5

if (case == 'case1') or (case == 'case2') or (case == 'case3') or (case == 'case4') or (case == 'case5') :
    WINDOW = 10
else :
    WINDOW = 50

if debug :
    WINDOW = 10000

#Load file and tree 
load_file = TFile(path+"B5_"+case+".root", "READ")
tree = load_file.Get('B5')

if momentumMode :

    ###Obtain all positions for DC1 and DC2

    #Event 0 to initiate numpy array and just vstack after
    tree.GetEntry(0)
    dc1_xpos = np.array( getattr(tree, 'Dc1HitsVector_x') ) [:5]
    dc1_ypos = np.array( getattr(tree, 'Dc1HitsVector_y') ) [:5]
    dc1_zpos = np.array( getattr(tree, 'Dc1HitsVector_z') ) [:5]

    dc2_xpos = np.array( getattr(tree, 'Dc2HitsVector_x') ) [:5]
    dc2_ypos = np.array( getattr(tree, 'Dc2HitsVector_y') ) [:5]
    dc2_zpos = np.array( getattr(tree, 'Dc2HitsVector_z') ) [:5]

    for entry in range(1, tree.GetEntries()):
        tree.GetEntry(entry)
        if (np.array( getattr(tree, 'Dc1HitsVector_x') ).size >= 5) and (np.array( getattr(tree, 'Dc2HitsVector_x') ).size >= 5) :
            dc1_xpos=np.vstack( ( dc1_xpos , np.array( getattr( tree , 'Dc1HitsVector_x' ) )[:5] ) )
            dc2_xpos=np.vstack( ( dc2_xpos , np.array( getattr( tree , 'Dc2HitsVector_x' ) )[:5] ) )
        else : 
            print("Too less x hits. Reject event #", entry)
            
        if (np.array( getattr(tree, 'Dc1HitsVector_y') ).size >= 5) and (np.array( getattr(tree, 'Dc2HitsVector_y') ).size >= 5) :
            dc1_ypos=np.vstack( ( dc1_ypos , np.array( getattr( tree , 'Dc1HitsVector_y' ) )[:5] ) )
            dc2_ypos=np.vstack( ( dc2_ypos , np.array( getattr( tree , 'Dc2HitsVector_y' ) )[:5] ) )
        else : 
            print("Too less y hits. Reject event #", entry)
        
        if (np.array( getattr(tree, 'Dc1HitsVector_z') ).size >= 5) and (np.array( getattr(tree, 'Dc2HitsVector_z') ).size >= 5) :
            dc1_zpos=np.vstack( ( dc1_zpos , np.array( getattr( tree , 'Dc1HitsVector_z' ) )[:5] ) )
            dc2_zpos=np.vstack( ( dc2_zpos , np.array( getattr( tree , 'Dc2HitsVector_z' ) )[:5] ) )
        else : 
            print("Too less z hits. Reject event #", entry)


    #Set z values
    zPosDC1 = np.array([-6.25, -5.75, -5.25, -4.75, -4.25])
    zPosDC2 = np.array([2.25, 2.75, 3.25, 3.75, 4.25])

    #Set x values to the good units
    dc1_xpos_metres=dc1_xpos/1000
    dc2_xpos_metres=dc2_xpos/1000

    p_values = []
    for event in range(dc1_xpos_metres.shape[0]):
        p = compute_p(event,B_TRUE,2) 
        if p > P_TRUE-WINDOW and p < P_TRUE+WINDOW : 
            p_values.append(p)

    m, s = stats.norm.fit(p_values) # get mean and standard deviation 

    resol = s/P_TRUE

    fig = plt.figure(figsize=(8,8))

    if debug :
        y, x, _ = plt.hist(p_values,bins=50,range=(min(p_values),max(p_values)),density=False);
    else :
        y, x, _ = plt.hist(p_values,bins=50,range=(m-WINDOW,m+WINDOW),density=False);

    plt.vlines(m,0,y.max()+2,linestyles='dashed')
    plt.xlabel('p [GeV]',fontsize=25)
    plt.ylabel('$N_{events}$',fontsize=25)
    plt.xticks(fontsize=15) 
    plt.yticks(fontsize=15)
    plt.title('$\mu=$'+str(round(m,2)) + '\t' + '$\sigma=$'+str(round(s,2)) + '\t' + 'Resolution = ' + str(round(resol,4)) ,fontsize=20 )

    plt.savefig('fig/pDistribution_'+case+'.png')
    plt.close()


    print("mu = ",m)
    print("sigma = ",s)

if calorimetryMode :
    ECEnergy = []
    for entry in range(0, tree.GetEntries()):
        tree.GetEntry(entry)
        if (np.array( getattr(tree, 'ECEnergy')) /1000 > 0): #Just to be sure there isn't any negative values
            ECEnergy.append( np.array( getattr(tree, 'ECEnergy')) /1000 )
        
    HCEnergy = []
    for entry in range(0, tree.GetEntries()):
        tree.GetEntry(entry)
        if (np.array( getattr(tree, 'HCEnergy')) /1000 > 0):
            HCEnergy.append( np.array( getattr(tree, 'HCEnergy')) /1000 )

    TotalEnergy = ECEnergy + HCEnergy

    #ECAL deposit plot 
    fig = plt.figure(figsize=(8,8))

    y, x, _ = plt.hist(ECEnergy,bins=100,range=(min(ECEnergy),max(ECEnergy)),density=False)

    plt.xlabel('ECAL energy deposit [GeV]',fontsize=25)
    plt.ylabel('$N_{events}$',fontsize=25)
    plt.xticks(fontsize=15) 
    plt.yticks(fontsize=15)

    plt.savefig('fig/ECALDeposit_'+case+'.png')
    plt.close()

    #HCAL deposit plot (basicaly the same code)
    fig = plt.figure(figsize=(8,8))

    y, x, _ = plt.hist(HCEnergy,bins=100,range=(min(HCEnergy),max(HCEnergy)),density=False)

    plt.xlabel('HCAL energy deposit [GeV]',fontsize=25)
    plt.ylabel('$N_{events}$',fontsize=25)
    plt.xticks(fontsize=15) 
    plt.yticks(fontsize=15);

    plt.savefig('fig/HCALDeposit_'+case+'.png')
    plt.close()

    #Total energy deposit
    fig = plt.figure(figsize=(8,8))

    y, x, _ = plt.hist(TotalEnergy,bins=100,range=(min(TotalEnergy),max(TotalEnergy)),density=False)

    plt.xlabel('Total energy deposit [GeV]',fontsize=25)
    plt.ylabel('$N_{events}$',fontsize=25)
    plt.xticks(fontsize=15) 
    plt.yticks(fontsize=15);

    plt.savefig('fig/TotalEnergyDeposit_'+case+'.png')
    plt.close()

print("End of execution \n")