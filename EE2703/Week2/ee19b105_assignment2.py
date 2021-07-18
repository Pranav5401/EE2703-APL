"""
        EE2703 Applied Programming Lab - 2021
            Assignment 2: EE19B105
"""
from sys import argv, exit
import os
from pathlib import Path
import numpy as np
import math
import cmath

class node():
    def __init__(self,name,number):
        self.name = name
        self.number = number                #GND will be 0 and as we parse each line we will give new number to each new node encountered
        self.voltage = 0                    #Assign 0 initially, will get values after solving matrix which can be inserted here


def matrix_formation():                     #Function to fill all values into the matrices M,b using the elemenets' values            
    if w!=0:
        M = np.zeros((len(nodes_list)+len(voltage_sources),len(nodes_list)+len(voltage_sources)),dtype=np.complex)       #Size of matrix is going to no. of nodes + no. of Voltage sources
        b = np.zeros(len(nodes_list)+ len(voltage_sources),dtype=np.complex)
    else:
        M = np.zeros((len(nodes_list)+len(voltage_sources),len(nodes_list)+len(voltage_sources)))
        b = np.zeros(len(nodes_list)+ len(voltage_sources))

    if w==0:
        for element in passive_elements:
            M[(element[-2].number),(element[-2].number)] += 1/element[0]                                  #Add 1/Z at the from-node,from-nodeth entry according to Modified Nodal Analysis    
            M[(element[-2].number),(element[-3].number)] += -1/element[0]                                 #Add -1/Z at the from-node,to-nodeth entry according to Modified Nodal Analysis
            M[(element[-3].number),(element[-3].number)] += 1/element[0]                                  #Add 1/Z at the to-node,to-nodeth entry according to Modified Nodal Analysis  
            M[(element[-3].number),(element[-2].number)] += -1/element[0]                                 #Add -1/Z at the to-node,from-nodeth entry according to Modified Nodal Analysis
            
        for source in voltage_sources:
            M[len(nodes_list)+voltage_sources.index(source),source[-2].number] += +1                      #Since we have the relation V(from-node) - V(to-node) = value of source   
            M[len(nodes_list)+voltage_sources.index(source),source[-3].number] += -1         #Following the notation that value of voltage source is potential rise from the from-node to the to-node  
            b[len(nodes_list)+voltage_sources.index(source)] = source[0] 
            
            M[source[-2].number,len(nodes_list)+voltage_sources.index(source)] += 1                       #Since we are assuming some current to flow through the voltage source
            M[source[-3].number,len(nodes_list)+voltage_sources.index(source)] += -1

        for source in current_sources:                                                      #Following the notation that value of current source is current flowing from the from-node to the to-node
            b[source[-2].number] += -source[0]                                                            #Current source's value being subtracted from the from-node 
            b[source[-3].number] += source[0]                                                             #Current source's value being added to the to-node                
            
        return M,b            
    else:
        for element in passive_elements:
            M[(element[-2].number),(element[-2].number)] += 1/element[0]                                  #Add 1/Z at the from-node,from-nodeth entry according to Modified Nodal Analysis    
            M[(element[-2].number),(element[-3].number)] += -1/element[0]                                 #Add -1/Z at the from-node,to-nodeth entry according to Modified Nodal Analysis
            M[(element[-3].number),(element[-3].number)] += 1/element[0]                                  #Add 1/Z at the to-node,to-nodeth entry according to Modified Nodal Analysis  
            M[(element[-3].number),(element[-2].number)] += -1/element[0]                                 #Add -1/Z at the to-node,from-nodeth entry according to Modified Nodal Analysis
            
        for source in voltage_sources:
            M[len(nodes_list)+voltage_sources.index(source),source[-2].number] += +1                      #Since we have the relation V(to-node) - V(from-node) = value of source   
            M[len(nodes_list)+voltage_sources.index(source),source[-3].number] += -1         #Following the notation that value of voltage source is potential rise from the from-node to the to-node  
            b[len(nodes_list)+voltage_sources.index(source)] = source[1] 
            
            M[source[-2].number,len(nodes_list)+voltage_sources.index(source)] += 1                       #Since we are assuming some current to flow through the voltage source
            M[source[-3].number,len(nodes_list)+voltage_sources.index(source)] += -1

        for source in current_sources:                                                      #Following the notation that value of current source is current flowing from the from-node to the to-node
            b[source[-2].number] += -source[1]                                              #Current source's value being subtracted from the from-node ( phase in radians ) 
            b[source[-3].number] += source[1]                                               #Current source's value being added to the to-node ( phase in radians )               
            
        return M,b            
        


def myspice(lines):                                             #Function which checks if all tokens are of valid form and stores them, also adds nodes to a distinct set            
    nodes_set = {'GND'}
    Output = []                                                 #Output list which has all the tokens stored in reverse order    
    line_count = len(lines)+1                                   #Since we are traversing in reversed order setting count at maximum + 1
    for line in reversed(lines):
        line_count = line_count-1
        tokens = line.split()
        l = len(tokens)
        
        if line[0] == "#":
            continue
            
        if l==4 or (l>4 and tokens[4][0] =='#'):               #check if length of tokens is 4 i.e type of component is R/L/C ( check for # as 5th element since there may be comments )
            element = tokens[0]
            n1 = tokens[1]
            n2 = tokens[2]
            value = tokens[3]
            
            
            if n1 == n2:
                print("The from-node and to-node cannot be same in line number ",line_count)
                exit()
            
            if not (n1.isalnum() and n2.isalnum()):            #Check if node names are alphanumeric, if not tell user netlist file has wrong inputs
                print("Node names should be alphanumeric")
                exit()
            
            if n1 in nodes_set:
                pass
            else:
                nodes_set.add(n1)        
            if n2 in nodes_set:
                pass
            else:
                nodes_set.add(n2)        
            
            token = [value,tokens[2],tokens[1],element]
            Output.append(token)
        
        elif l == 6 or (l>6 and tokens[6][0] =='#'):           #check if length of tokens is 6 i.e type of component is ac source ( check for # as 7th element since there may be comments )
            element = tokens[0]
            n1 = tokens[1]
            n2 = tokens[2]
            ac = tokens[3]
            value = tokens[4]
            phase = tokens[5]
            
            if n1 == n2:
                print("The from-node and to-node cannot be same in line number ",line_count)
                exit() 
            
            if not (n1.isalnum() and n2.isalnum()):      #Check if node names are alphanumeric, if not tell user netlist file has wrong inputs
                print("Node names should be alphanumeric")
                exit()
            
            if n1 in nodes_set:
                pass
            else:
                nodes_set.add(n1)       
            if n2 in nodes_set:
                pass
            else:
                nodes_set.add(n2)
            
            token = [phase,value,ac,tokens[2],tokens[1],element]
            Output.append(token)
            
        elif l==5 or (l>5 and tokens[5][0] =='#'):             #check if length of tokens is 5 i.e type of component is dc source ( check for # as 6th element since there may be comments )
            element = tokens[0]
            n1 = tokens[1]
            n2 = tokens[2]
            dc = tokens[3]
            value = tokens[4]
            
            if n1 == n2:
                print("The from-node and to-node cannot be same in line number ",line_count)
                exit()
            
            if not (n1.isalnum() and n2.isalnum()):            #Check if node names are alphanumeric, if not tell user netlist file has wrong inputs
                print("Node names should be alphanumeric")
                exit()
            
            if n1 in nodes_set:
                pass
            else:
                nodes_set.add(n1)       
            if n2 in nodes_set:
                pass
            else:
                nodes_set.add(n2)
            
            token = [value,dc,tokens[2],tokens[1],element]
            Output.append(token)
        
        else:
            print('The given .netlist file has entries of wrong format between .circuit and .end lines')        
    
    return Output,nodes_set



CIRCUIT = '.circuit'
END = '.end'
AC = '.ac'

if len(argv) != 2:                                           #Making sure correct number of arguments are given as input       
    print('\nUsage: %s <inputfile>' % argv[0])
    exit()

if os.path.isfile(Path(argv[1])) == False:      #Checking if such a file exists, if not show error ( instead of doing try except if loop can also take care of error where user gives wrong input )  
    print('\nPlease give input file path as an existing one')
    exit()
    
    
with open(argv[1]) as f:
    lines = f.readlines()
    start = -1; end = -2; w=0
    for line in lines:              # extracting circuit definition start and end lines
        if AC == line[:len(AC)]:           #Check if circuit has AC element
            ac_tokens = line.split()    
            try:
                w = float(ac_tokens[-1])     #Assign and store its angular frequency if AC element is there
            except:
                print("The frequency of an element should be numeric")
                exit()       
        if CIRCUIT == line[:len(CIRCUIT)]:
            start = lines.index(line)
        elif END == line[:len(END)]:
            end = lines.index(line)
            
    if start >= end:                # validating circuit block
        print('Invalid circuit definition')
        exit(0)

Output,nodes_set = myspice(lines[start+1:end])
passive_elements = []               #Contains all the R,L,C
voltage_sources = []                #Contains all the independant voltage sources
current_sources = []                #Contains all the independant current sources
count = 0                           #To check for ground connection in circuit

for token in Output:
    try:
        token[0] = float(token[0])
    except:
        print("Value of each element should be numeric")
        exit()    
        
for token in Output:
    if token[-2] == "GND" or token[-3] == "GND":
        count = count + 1
if count == 0:
    print("There is no ground in the circuit, error!")
    exit()

nodes_list = [node('GND',0)]                                         #list to which we will add entries of the class node by creating them in the next for loop
node_number = 1                                         #0 will be directly assigned to GND node
for x in nodes_set:
    if x == 'GND':
        continue
    else:
        nodes_list.append(node(x,node_number))
        node_number += 1
                
for token in Output:
    for x in nodes_list:
        if x.name == token[-2]:
            token[-2] = x
        elif x.name == token[-3]:
            token[-3] = x    


for token in Output:
    if token[-1][0] == 'R' or token[-1][0] == 'L' or token[-1][0] == 'C':
        passive_elements.append(token)
    elif token[-1][0] == 'V':
        voltage_sources.append(token)
    elif token[-1][0] == 'I':
        current_sources.append(token)
        
print("The frequency is ",w)
w = 2*math.pi*w                                          #Converting Hz frequency to radians/sec
for element in passive_elements:
    if element[-1][0] == "C" and w != 0:
        element[0] = complex(0,-1/(w*element[0]))
    elif element[-1][0] == "C" and w == 0:              #For DC sources C is like a open circuit
        element[0] = 1e20
    elif element[-1][0] == "L" and w != 0:
        element[0] = complex(0,(w*element[0]))
    elif element[-1][0] == "L" and w == 0:
        element[0] = 1e-20

for source in voltage_sources:
    if (w!=0 and source[-4] == 'ac'):
        source[1] = (float(source[1])/2)*complex(math.cos(source[0]),math.sin(source[0]))     #Since p-p values are given halve it, store the complex value   
    if (w!=0 and source[-4] == 'dc'):
        source[1] = float(source[0])                    #While forming matrix I have done different cases for w = 0 and w!=0 so when w!=0 value of source I am using to feed matrix is from source[1]
    if (w==0):   
        source[0] = float(source[0])    
for source in current_sources:
    if (w!=0 and source[-4] == 'ac'):
        source[1] = (float(source[1])/2)*complex(math.cos(source[0]),math.sin(source[0]))     #Since p-p values are given halve it, store the complex value   
    if (w!=0 and source[-4] == 'dc'):
        source[1] = float(source[0])                    #While forming matrix I have done different cases for w = 0 and w!=0 so when w!=0 value of source I am using to feed matrix is from source[1]
    if (w==0):   
        source[0] = float(source[0])    

        
M,b = matrix_formation()                                #Calling the matrix formation function    
print(M)
print(b)
       
try:
    X = np.linalg.solve(M,b)
except:
    print("Unsolvable Matrix")
    exit()    

#Since this X has final values for all voltages but we did not use the criterion that GND is at 0 voltage till now, we should now shift all nodal voltages so that GND is at 0 
temp = X[0]
for i in range(len(nodes_list)):
    X[i] = X[i] - temp

num_node = 0                 #Will use it to assign voltages to each node in the following for loop
for n in nodes_list:
    n.voltage = X[num_node] 
    num_node += 1
    print("The voltage at node " +n.name+ " is ",n.voltage) 

for source in voltage_sources:
    print("The current through voltage source " +source[-1]+ " is ",X[len(nodes_list)+voltage_sources.index(source)])              
    
    
    
    
