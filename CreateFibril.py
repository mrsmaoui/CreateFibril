#!/usr/bin/env python
import os
import math
import string
import fileinput
import tempfile
import sys
import matrix
from subprocess import call

#-------------------------- START PARSING INPUT COMMAND --------------------------#

error = 'ERROR: ENTER COMMAND AS FOLLOWS:\n./CreateFibril.py -f filename -c chains -r repetitions\n'
error = error + '\t-f filename: pdb file ex. "-f 2RNM.pdb" (required)\n'
error = error + '\t-c chains: chains to repeat ex. "-c 2,4" (required)\n'
error = error + '\t-r repetitions: number of times to repeat chain(s) ex. "-r 4" (required integer)\n'
error = error + '\t-d distance: max distance in Angstroms between the new added chain(s) and the last chain in filename ex. "-d 4" (required int)\n'
error = error + '\t-a angle: the rotation angle in degrees between the new added chains and the last chain in filename ex "-a 30.2" (float) [0, 360] (required)\n'
error = error + '\t-Polygon n: Create a fibril with n protofilaments in a polygon conformation\n'
error = error + '\t-Ring n: Create a fibril with n protofilaments in a ring conformation\n'
error = error + '\t-Radius r: Create a fibril with n protofilaments centered r Angstroms away from axis of rotation (required if you select Polygon or Ring)\n'
error = error + '\t-fp x,y,z,residue : Fibril axis initial point + residue number it is on ex. "-fp 2.452,3.343,4.433,37" (required)\n'
error = error + '\t-nr n : Number of residues ex. "-nr 37" (required)\n'

#verifying input variables are correct
for i in range(len(sys.argv)):
    if( i % 2 == 1):	#check odd values
        #exit if command doesn't exist
        if not (sys.argv[i] == '-f' or sys.argv[i] == '-c' or sys.argv[i] == '-r' or sys.argv[i] == '-d' or sys.argv[i] == '-a' or sys.argv[i] == '-Polygon' or sys.argv[i] == '-Ring' or sys.argv[i] == '-Radius' or sys.argv[i] == '-fp' or sys.argv[i] == '-nr'):
            print error
            sys.exit(0)

#get input parameters
filename = ''
chains = ''
repititions = ''
distance = ''
iterate_distance = ''
angle = ''
iterate_angle = ''
polygon = ''
ring = ''
radius = ''
fp = ''
nr = ''

try:
    for i in range(len(sys.argv)):
        if(sys.argv[i] == '-f'): filename = sys.argv[i+1]
        elif (sys.argv[i] == '-c'): chains = sys.argv[i+1]
        elif (sys.argv[i] == '-r'): repititions = sys.argv[i+1]
        elif (sys.argv[i] == '-d'): distance = sys.argv[i+1]
        elif (sys.argv[i] == '-a'): angle = sys.argv[i+1]
        elif (sys.argv[i] == '-Polygon'): polygon = sys.argv[i+1]
        elif (sys.argv[i] == '-Ring'): ring = sys.argv[i+1]
        elif (sys.argv[i] == '-Radius'): radius = sys.argv[i+1]
        elif (sys.argv[i] == '-fp'): fp = sys.argv[i+1]
        elif (sys.argv[i] == '-nr'): nr = int(sys.argv[i+1])

except Exception:
    print error
    sys.exit(0)

#check that input file is of pdb format
if not filename[len(filename)-4:len(filename)] == '.pdb':
    print error
    sys.exit(1)

filename = filename[0:(len(filename)-4)]

#parse chains, repititions, distance, and angle and check if bounds are integers or floats
chain = [0, 0]

try:
    chain[0] = int(chains[0:1])
    chain[1] = int(chains[2:3])
    repititions = int(repititions)
    if(distance == ''):
        distance = 4.8
        iterate_distance = 'yes'
    else:
        distance = float(distance)
        iterate_distance = 'no'

    if(angle == ''):
        angle = 35
        iterate_angle = 'no'
    else:
        angle = float(angle)
        iterate_angle = 'no'


except ValueError:
    print "here"
    print error
    sys.exit(0)

if(chain[1] < chain[0] or chain[0] <= 0 or chain[1] == 1 or repititions <= 0): 
	print error
	sys.exit(0)

if(angle < -360 or angle > 360 or distance < 0):
	print error
	sys.exit(0)

if(nr == ''):
    print error
    sys.exit(0)

LJ = 'false'
Coulomb = 'false'
Solvation = 'false'
temperature = ''


if( polygon == '0' ):
    polygon = ''
if( ring == '0'):
    ring = ''

if( not( polygon == '') or not( ring == '') ):
    if( radius == '' ):
        print error
        exit(0)
    else:
        radius = int(radius) 
                
       

#there are 2 locations for ext to edit. This one and another.
ext = 'r%sc%s%sd%sa%s' % (repititions, chain[0], chain[1], distance, angle)

		#-------------------------- END PARSING INPUT COMMAND --------------------------#

def main():
    global distance
    global radius
    global angle
    global ext
    fixed_angle = angle

    # Create 1 protofilament
    if(polygon == '' and ring == ''):
        createSingleStructure()

    # RINGS OR POLYGONS
    if not (polygon == ''):
        createPolygon(polygon)
    if not (ring == ''):
        createRing(ring)
    

def createSingleStructure():

####################### AUTOMATE THIS: DONE ##########################
        
    point = fp.split(',')           #point for fibril axis direction
    p = [float(point[0]), float(point[1]), float(point[2])]  

####################### AUTOMATE THIS ##########################

    filament = createProtofilament(getMonomerInfo(1, 'single'), distance, angle, p, 0)
    createFile(filament)

        
        

#---------- Appends all the filaments to an output file -------------# 
def createFile(filaments):
        file_pdb_transform = open(filename +ext+ '_transform.pdb', 'w')
        for line in filaments:
                file_pdb_transform.write(line)
        file_pdb_transform.close()
        call('cat %s%s_init.pdb > %s%s.pdb' % (filename, ext, filename, ext), shell=True)
        call('cat %s%s_transform.pdb >> %s%s.pdb' % (filename, ext, filename, ext), shell=True)
        call('rm -f %s%s_r_chain.pdb %s%s_repeat_chain.pdb %s%s_init.pdb %s%s_transform.pdb' % (filename, ext, filename, ext, filename, ext, filename, ext), shell=True)
	call('chmod a+r %s%s.pdb' % (filename, ext), shell=True)



def createProtofilament(monomerInfo, d, ang, rotation_point, extra_rotation):
        global ext                              		
        ext = 'r%sc%s%sd%sa%s' % (repititions, chain[0], chain[1], d, ang)
        filament = []    
	# monomerInfo = [x, y, z, last atom, last chain, last residue, base_distance, rotation_point]
        axis_direction = [monomerInfo[0], monomerInfo[1], monomerInfo[2]]
        number_of_atoms = int(monomerInfo[3])
        number_of_chains = monomerInfo[4]
        base_distance = monomerInfo[5]
        number = monomerInfo[6]

		
	# for loop on number of repititions
	i = 0
	while( i < int(repititions)):
		i = i + 1
			
		# translate: for each atom
                t = [ [1,0,0,(axis_direction[0]*(base_distance+d)*(i-1))], [0,1,0,(axis_direction[1]*(base_distance+d)*(i-1))], [0,0,1,(axis_direction[2]*(base_distance+d)*(i-1))], [0,0,0,1] ]
			
		#translate rotation point
		rp = matrix.zero(4,1)
		rp[0][0] = rotation_point[0]
		rp[1][0] = rotation_point[1]
		rp[2][0] = rotation_point[2]
		rp[3][0] = 1
		range_rotation_point = matrix.mult(t, rp)
	
		theta = ang*(i-1) + extra_rotation
		r = computeRotationMatrix(range_rotation_point, axis_direction, theta)

		#current_chain = last_chain

                #open file'
		file_pdb = open(filename+ext+'_repeat_chain.pdb', 'r')
		fileLines = file_pdb.readlines()
	        for one_line in fileLines:
	        	ones = string.split(one_line)
			if(len(ones) > 1):
	       		        if (ones[0] == 'ATOM'):
                                        ones[4] = chr( ord(ones[4]) + number_of_chains*(i-1) + number_of_chains*repititions*(number-1) ) 
                                        ones[1] = int(ones[1]) + int(number_of_atoms)*(i-1) + (int(number_of_atoms)*(repititions))*(number-1)
                                        

####################### AUTOMATE THIS: DONE ##########################

                                        ones[5] = int(ones[5]) + (ord(ones[4]) - ord('A') )*(nr)    #number of residues
                                        #ones[5] = int(ones[5]) + (ord(ones[4]) - ord('A') )*(37)    #number of residues is 37 for amylin

####################### AUTOMATE THIS ##########################

                                        
                                        #special fix for Capital and minuscale same chains
                                        if( ord(ones[4]) - ord('a') >= 0):
                                                ones[4] = chr( (ord(ones[4]) + 26) )
                                        
                                        p = matrix.zero(4,1)
                                        p[0][0] = float(ones[6])
                                        p[1][0] = float(ones[7])
                                        p[2][0] = float(ones[8])
                                        p[3][0] = 1
							
                                        new_positions_3d = matrix.mult( r, matrix.mult(t, p) )
                                        ones[0] = addSpaces(ones[0], 6, 'after')
                                        ones[1] = addSpaces(ones[1], 5, 'before')
                                        ones[2] = addSpaces(ones[2], 4, 'after')
                                        ones[3] = addSpaces(ones[3], 3, 'before')
                                        ones[5] = addSpaces(ones[5], 4, 'before')
                                        new_positions_3d[0][0] = '%.3f' % new_positions_3d[0][0]
                                        new_positions_3d[0][0] = addSpaces(new_positions_3d[0][0], 8, 'before')
                                        new_positions_3d[1][0] = '%.3f' % new_positions_3d[1][0]
                                        new_positions_3d[1][0] = addSpaces(new_positions_3d[1][0], 8, 'before')
                                        new_positions_3d[2][0] = '%.3f' % new_positions_3d[2][0]
                                        new_positions_3d[2][0] = addSpaces(new_positions_3d[2][0], 8, 'before')
                                        new_line = '%s%s %s %s %s%s    %s%s%s\n' % (ones[0], ones[1], ones[2], ones[3], ones[4], ones[5], new_positions_3d[0][0], new_positions_3d[1][0], new_positions_3d[2][0])
					filament.append(new_line)
										
                file_pdb.close()
		# end repititions for loop

	return filament	
        

       
# ------------ Creates a polygon based fibril --------------- #	
def createPolygon(number):
    if(number == ''):
        return
    else:
        number = int(number)

        ######################## AUTOMATE THIS ###########################
        # strand vector
        #ATOM    215  CA  -0.679  -4.351  29.876
        #ATOM    230  CA  5.499  -7.115  30.352
        #mag = math.sqrt(6.178*6.178 + 2.764*2.764 + 0.476*0.476)
        #d1 = [-6.178/mag, 2.764/mag, -0.476/mag]

        print "\n*******************************\nCreateFibril will build a Polygon structure. Please specify the following parameters:\n"
        ep3 = raw_input("Enter mid point of polygon structure (Point on the center of the facing strand) x,y,z: ") # p = [-23.881,-0.306,28.842]
        ep1 = raw_input("Enter point1 of strand vector x,y,z: ")  #-6.970,-2.332,29.758
        ep2 = raw_input("Enter point2 of strand vector x,y,z: ")  #14.878,-11.615,29.695
        print "\n"

        mag = math.sqrt( math.pow( (float(ep2.split(',')[0])-float(ep1.split(',')[0])),2) + math.pow( (float(ep2.split(',')[1])-float(ep1.split(',')[1])),2) +math.pow( (float(ep2.split(',')[2])-float(ep1.split(',')[2])),2) )
        d1 = [ (float(ep1.split(',')[0])-float(ep2.split(',')[0]))/mag, (float(ep1.split(',')[1])-float(ep2.split(',')[1]))/mag, (float(ep1.split(',')[2])-float(ep2.split(',')[2]))/mag]

        # Fibril axis direction
        #ATOM    222  CA  3.121  -4.207  29.774
        #ATOM    222  CA  2.897  -4.365  34.574
        #result = [0.224, 0.158, -4.8]
        #x1 = 0.224
        #y1 = 0.158
        #z1 = -4.8	

        #mag2 = math.sqrt(0.224*0.224 + 0.158*0.158 + 4.8*4.8)
        #x1 = x1/mag2
        #y1 = y1/mag2
        #z1 = z1/mag2
        res = getMonomerInfo(1, 'polygon')

        x1 = float(res[0])
        y1 = float(res[1])
        z1 = float(res[2])

        x2 = -d1[0]
        y2 = -d1[1]
        z2 = -d1[2]
        # cross product vector
        d2 = [ - y1*z2 + y2*z1, - z1*x2 + z2*x1, - x1*y2 + x2*y1]

        #p = [3.121, -4.207, 29.774] #ATOM    222  CA  ASN A 201       3.121  -4.207  29.774, point on beta sheet starting of fibril
        p = [float(ep3.split(',')[0]), float(ep3.split(',')[1]), float(ep3.split(',')[2])]  # head point of ring structure. Point on the curve of the amyloid

        #########################################################################
	

        if(number == 2):
                
                # range of rotation points a in [-7, 4] b in [5, 10]
                a = 0
                b = radius
                rotation_point = [0, 0, 0]
                rotation_point[0] = p[0] + ( a*d1[0] + b*d2[0] )
                rotation_point[1] = p[1] + ( a*d1[1] + b*d2[1] )
                rotation_point[2] = p[2] + ( a*d1[2] + b*d2[2] )
                monomerInfo = getMonomerInfo(1, 'polygon')
                filament1 = createProtofilament(monomerInfo, distance, angle, rotation_point, 0)
                monomerInfo = getMonomerInfo(2, 'polygon')
                filament2 = createProtofilament(monomerInfo, distance, angle, rotation_point, 180) 
                createFile(filament1 + filament2)
                
        if(number == 3):
                
               
                # range of rotation points a in [,]
                a = radius
                rotation_point = [0, 0, 0]
                rotation_point[0] = p[0] + ( a*d2[0] )
                rotation_point[1] = p[1] + ( a*d2[1] )
                rotation_point[2] = p[2] + ( a*d2[2] )

                
                monomerInfo = getMonomerInfo(1, 'polygon')
                filament1 = createProtofilament(monomerInfo, distance, angle, rotation_point, 0)
                monomerInfo = getMonomerInfo(2, 'polygon')
                filament2 = createProtofilament(monomerInfo, distance, angle, rotation_point, 120) 
                monomerInfo = getMonomerInfo(3, 'polygon')
                filament3 = createProtofilament(monomerInfo, distance, angle, rotation_point, 240)
                createFile(filament1 + filament2 + filament3)
                
                
        if(number == 4):
                
                # range of rotation points a in [,]
                a = radius
                rotation_point = [0, 0, 0]
                rotation_point[0] = p[0] + ( a*d2[0] )
                rotation_point[1] = p[1] + ( a*d2[1] )
                rotation_point[2] = p[2] + ( a*d2[2] )
                monomerInfo = getMonomerInfo(1, 'polygon')
                filament1 = createProtofilament(monomerInfo, distance, angle, rotation_point, 0)
                monomerInfo = getMonomerInfo(2, 'polygon')
                filament2 = createProtofilament(monomerInfo, distance, angle, rotation_point, 90) 
                monomerInfo = getMonomerInfo(3, 'polygon')
                filament3 = createProtofilament(monomerInfo, distance, angle, rotation_point, 180)
                monomerInfo = getMonomerInfo(4, 'polygon')
                filament4 = createProtofilament(monomerInfo, distance, angle, rotation_point, 270)
                createFile(filament1 + filament2 + filament3 + filament4)

        if(number == 5):
                
                # range of rotation points a in [,]
                a = radius
                rotation_point = [0, 0, 0]
                rotation_point[0] = p[0] + ( a*d2[0] )
                rotation_point[1] = p[1] + ( a*d2[1] )
                rotation_point[2] = p[2] + ( a*d2[2] )
                monomerInfo = getMonomerInfo(1, 'polygon')
                filament1 = createProtofilament(monomerInfo, distance, angle, rotation_point, 0)
                monomerInfo = getMonomerInfo(2, 'polygon')
                filament2 = createProtofilament(monomerInfo, distance, angle, rotation_point, 72) 
                monomerInfo = getMonomerInfo(3, 'polygon')
                filament3 = createProtofilament(monomerInfo, distance, angle, rotation_point, 144)
                monomerInfo = getMonomerInfo(4, 'polygon')
                filament4 = createProtofilament(monomerInfo, distance, angle, rotation_point, 216)
                monomerInfo = getMonomerInfo(5, 'polygon')
                filament5 = createProtofilament(monomerInfo, distance, angle, rotation_point, 288)
                createFile(filament1 + filament2 + filament3 + filament4 + filament5)

        if(number == 6):
                
                # range of rotation points a in [,]
                a = radius
                rotation_point = [0, 0, 0]
                rotation_point[0] = p[0] + ( a*d2[0] )
                rotation_point[1] = p[1] + ( a*d2[1] )
                rotation_point[2] = p[2] + ( a*d2[2] )
                monomerInfo = getMonomerInfo(1, 'polygon')
                filament1 = createProtofilament(monomerInfo, distance, angle, rotation_point, 0)
                monomerInfo = getMonomerInfo(2, 'polygon')
                filament2 = createProtofilament(monomerInfo, distance, angle, rotation_point, 60) 
                monomerInfo = getMonomerInfo(3, 'polygon')
                filament3 = createProtofilament(monomerInfo, distance, angle, rotation_point, 120)
                monomerInfo = getMonomerInfo(4, 'polygon')
                filament4 = createProtofilament(monomerInfo, distance, angle, rotation_point, 180)
                monomerInfo = getMonomerInfo(5, 'polygon')
                filament5 = createProtofilament(monomerInfo, distance, angle, rotation_point, 240)
                monomerInfo = getMonomerInfo(6, 'polygon')
                filament6 = createProtofilament(monomerInfo, distance, angle, rotation_point, 300)
                
                createFile(filament1 + filament2 + filament3 + filament4 + filament5 +filament6)

    

# ------------ Creates a Ring based fibril --------------- #	
def createRing(number):
        if(number == ''):
                return
        else:
                number = int(number)

####################### AUTOMATE THIS ##########################                
        # strand vector, Radius direction

        print "\n*******************************\nCreateFibril will build a Ring structure. Please specify the following parameters:\n"
        ep3 = raw_input("Enter head point of ring structure (Point on the curve of the amyloid) x,y,z: ") # p = [-23.881,-0.306,28.842]

        ep1 = raw_input("Enter point1 of strand vector x,y,z: ")  #-6.970,-2.332,29.758
        ep2 = raw_input("Enter point2 of strand vector x,y,z: ")  #14.878,-11.615,29.695
        print "\n"
        
        mag = math.sqrt( math.pow( (float(ep2.split(',')[0])-float(ep1.split(',')[0])),2) + math.pow( (float(ep2.split(',')[1])-float(ep1.split(',')[1])),2) +math.pow( (float(ep2.split(',')[2])-float(ep1.split(',')[2])),2) )
        d1 = [ (float(ep1.split(',')[0])-float(ep2.split(',')[0]))/mag, (float(ep1.split(',')[1])-float(ep2.split(',')[1]))/mag, (float(ep1.split(',')[2])-float(ep2.split(',')[2]))/mag]

        #mag2 = math.sqrt(21.848*21.848 + 9.283*9.283 + 0.063*0.063)
        #d12 = [-21.848/mag2, 9.283/mag2, 0.063/mag2]
        #p2 = [-23.881, -0.306, 28.842]     

        p = [float(ep3.split(',')[0]), float(ep3.split(',')[1]), float(ep3.split(',')[2])]  # head point of ring structure. Point on the curve of the amyloid
        

####################### AUTOMATE THIS ##########################
 
        if(number == 2):
                # range of rotation points a in [,]
                a = radius

                rotation_point = [0, 0, 0]
                rotation_point[0] = p[0] + ( a*d1[0] )
                rotation_point[1] = p[1] + ( a*d1[1] )
                rotation_point[2] = p[2] + ( a*d1[2] )
                
                monomerInfo = getMonomerInfo(1, 'ring')
                filament1 = createProtofilament(monomerInfo, distance, angle, rotation_point, 0)
                monomerInfo = getMonomerInfo(2, 'ring')
                filament2 = createProtofilament(monomerInfo, distance, angle, rotation_point, 180) 
                createFile(filament1 + filament2)
                             
        
        if(number == 3):
                # range of rotation points a in [,]
                a = radius
                rotation_point = [0, 0, 0]
                rotation_point[0] = p[0] + ( a*d1[0] )
                rotation_point[1] = p[1] + ( a*d1[1] )
                rotation_point[2] = p[2] + ( a*d1[2] )
                monomerInfo = getMonomerInfo(1, 'ring')
                filament1 = createProtofilament(monomerInfo, distance, angle, rotation_point, 0)
                monomerInfo = getMonomerInfo(2, 'ring')
                filament2 = createProtofilament(monomerInfo, distance, angle, rotation_point, 120) 
                monomerInfo = getMonomerInfo(3, 'ring')
                filament3 = createProtofilament(monomerInfo, distance, angle, rotation_point, 240)
                #createFile(filament1 + filament2 + filament3)
                createFile(filament1 + filament2 + filament3)

                
        if(number == 4):
                # range of rotation points a in [,]
                a = radius
                rotation_point = [0, 0, 0]
                rotation_point[0] = p[0] + ( a*d1[0] )
                rotation_point[1] = p[1] + ( a*d1[1] )
                rotation_point[2] = p[2] + ( a*d1[2] )
                monomerInfo = getMonomerInfo(1, 'ring')
                filament1 = createProtofilament(monomerInfo, distance, angle, rotation_point, 0)
                monomerInfo = getMonomerInfo(2, 'ring')
                filament2 = createProtofilament(monomerInfo, distance, angle, rotation_point, 90) 
                monomerInfo = getMonomerInfo(3, 'ring')
                filament3 = createProtofilament(monomerInfo, distance, angle, rotation_point, 180)
                monomerInfo = getMonomerInfo(4, 'ring')
                filament4 = createProtofilament(monomerInfo, distance, angle, rotation_point, 270)
                createFile(filament1 + filament2 + filament3 + filament4)
                

        if(number == 5):
                # range of rotation points a in [,]
                a = radius
                rotation_point = [0, 0, 0]
                rotation_point[0] = p[0] + ( a*d1[0] )
                rotation_point[1] = p[1] + ( a*d1[1] )
                rotation_point[2] = p[2] + ( a*d1[2] )
                monomerInfo = getMonomerInfo(1, 'ring')
                filament1 = createProtofilament(monomerInfo, distance, angle, rotation_point, 0)
                monomerInfo = getMonomerInfo(2, 'ring')
                filament2 = createProtofilament(monomerInfo, distance, angle, rotation_point, 72) 
                monomerInfo = getMonomerInfo(3, 'ring')
                filament3 = createProtofilament(monomerInfo, distance, angle, rotation_point, 144)
                monomerInfo = getMonomerInfo(4, 'ring')
                filament4 = createProtofilament(monomerInfo, distance, angle, rotation_point, 216)
                monomerInfo = getMonomerInfo(5, 'ring')
                filament5 = createProtofilament(monomerInfo, distance, angle, rotation_point, 288)
                createFile(filament1 + filament2 + filament3 + filament4 + filament5)

        if(number == 6):
                # range of rotation points a in [,]
                a = radius
                rotation_point = [0, 0, 0]
                rotation_point[0] = p[0] + ( a*d1[0] )
                rotation_point[1] = p[1] + ( a*d1[1] )
                rotation_point[2] = p[2] + ( a*d1[2] )
                monomerInfo = getMonomerInfo(1, 'ring')
                filament1 = createProtofilament(monomerInfo, distance, angle, rotation_point, 0)
                monomerInfo = getMonomerInfo(2, 'ring')
                filament2 = createProtofilament(monomerInfo, distance, angle, rotation_point, 60) 
                monomerInfo = getMonomerInfo(3, 'ring')
                filament3 = createProtofilament(monomerInfo, distance, angle, rotation_point, 120)
                monomerInfo = getMonomerInfo(4, 'ring')
                filament4 = createProtofilament(monomerInfo, distance, angle, rotation_point, 180)
                monomerInfo = getMonomerInfo(5, 'ring')
                filament5 = createProtofilament(monomerInfo, distance, angle, rotation_point, 240)
                monomerInfo = getMonomerInfo(6, 'ring')
                filament6 = createProtofilament(monomerInfo, distance, angle, rotation_point, 300)
                
                createFile(filament1 + filament2 + filament3 + filament4 + filament5 +filament6)


#----------- CHECKS FOR CORRECT RANGE -----------#
def checkCorrectRange(chain_start, chain_end, last_chain_file):
	last = ord(last_chain_file) - ord('A') + 1
	if last >= int(chain_end):
		return True
	else:
		return False
#----------- END OF CORRECT RANGE CHECKING -----------#	


def addSpaces(v, vlength, insert):
	v = str(v)
	spaces = vlength - len(v)

	for i in range(spaces):
		if insert == 'before':
			v = ' ' + v
		else:
			v = v + ' '

	return v

def getMonomerInfo(number, conformation):
	
	#1. creates upto the desired chains to build
	#2. create the file that has only the range to replicate

	file_pdb = open(filename+'.pdb', 'r')
        fileLines = file_pdb.readlines()
	first_line = ''
	last_line = ''

	
	# output file
        file_pdb_transform = open(filename + ext+'_init.pdb', 'w')

	# file containing chains to repeat
	file_repeat_chain = open(filename + ext+'_repeat_chain.pdb', 'w')

	normal = []
	CA_Fpoints = []
	CA_Ipoints = []

	first = 0
	for one_line in fileLines:
        	ones = string.split(one_line)
		if(len(ones) > 1):
                	if (ones[0] == 'ATOM' and ( ord(ones[4]) - ord('A') + 1) <= int(chain[1])):
				
				#format file correctly
				ones[0] = addSpaces(ones[0], 6, 'after')
				ones[1] = addSpaces(ones[1], 5, 'before')
				ones[2] = addSpaces(ones[2], 4, 'after')
                                ones[3] = addSpaces(ones[3], 3, 'before')
                                ones[5] = addSpaces(ones[5], 4, 'before')
				ones[6] = addSpaces(ones[6], 8, 'before')
				ones[7] = addSpaces(ones[7], 8, 'before')
				ones[8] = addSpaces(ones[8], 8, 'before')


####################### AUTOMATE THIS ##########################
				
                                residueFibrilAxisPassesThrough = fp.split(',')[3]
                                if(conformation == 'single'): # 2 points determine a vector Amylin uses CA 30
                                        # populate CA_Fpoints from all chains, will be used for direction vector
                                        
                                        residueFibrilAxisPassesThrough = fp.split(',')[3]
                                        if(ones[5] == addSpaces(residueFibrilAxisPassesThrough, 4, 'before') and ones[2] == addSpaces('CA', 4, 'after')):
                                        	CA_Fpoints.append([ float(ones[6]), float(ones[7]), float(ones[8]) ])
                                                			
                                        # will be used to calculate the distance from CA 200 on Chain[1] to CA 200 on Chain[0]
                                        if(ones[5] == addSpaces(residueFibrilAxisPassesThrough, 4, 'before') and ones[2] == addSpaces('CA', 4, 'after')):
                                                CA_Ipoints.append([ float(ones[6]), float(ones[7]), float(ones[8]) ])
                                                                 

                                if(conformation == 'polygon'):
                                        # populate CA_Fpoints from all chains, will be used for direction vector
                                        if(ones[5] == addSpaces(residueFibrilAxisPassesThrough, 4, 'before') and ones[2] == addSpaces('CA', 4, 'after')):
                                        	CA_Fpoints.append([ float(ones[6]), float(ones[7]), float(ones[8]) ])
			
                                        # will be used to calculate the distance from CA 268 on Chain[1] to CA 232 on Chain[0]
                                        if(ones[5] == addSpaces(residueFibrilAxisPassesThrough, 4, 'before') and ones[2] == addSpaces('CA', 4, 'after')):
                                                CA_Ipoints.append([ float(ones[6]), float(ones[7]), float(ones[8]) ])

                                if(conformation == 'ring'):
                                        # populate CA_Fpoints from all chains, will be used for direction vector
                                        if(ones[5] == addSpaces(residueFibrilAxisPassesThrough, 4, 'before') and ones[2] == addSpaces('CA', 4, 'after')):
                                        	CA_Fpoints.append([ float(ones[6]), float(ones[7]), float(ones[8]) ])
			
                                        # will be used to calculate the distance from CA 211
                                        if(ones[5] == addSpaces(residueFibrilAxisPassesThrough, 4, 'before') and ones[2] == addSpaces('CA', 4, 'after')):
                                                CA_Ipoints.append([ float(ones[6]), float(ones[7]), float(ones[8]) ])
                                                
####################### AUTOMATE THIS ##########################
                                                
				new_line = '%s%s %s %s %s%s    %s%s%s\n' % (ones[0], ones[1], ones[2], ones[3], ones[4], ones[5], ones[6], ones[7], ones[8])
                                
				if(ord(ones[4]) - ord('A') + 1 < int(chain[0]) ):
					file_pdb_transform.write(new_line)
					
	
				# write to chain file
				if( ord(ones[4]) - ord('A') + 1 >= int(chain[0]) ):
					file_repeat_chain.write(new_line)
					if first == 0:
                                                first = 1
                                                first_line = ones
                                        last_line = ones
        
	file_pdb.close()
	file_pdb_transform.close()
	file_repeat_chain.close()	

	
	# find direction of fibril axis as (x, y, z)
	# get 2 points, same N atom on same residue number of two different chains
	# subtract two points to get direction of axis (smalled chain number from bigger)
	# for 2RNM we use the CA atom on residue 266	
	normal_len = len(normal)
	
	x = CA_Fpoints[chain[1] - 1][0] - CA_Fpoints[chain[1] - 2][0]
	y = CA_Fpoints[chain[1] - 1][1]	- CA_Fpoints[chain[1] - 2][1]
	z = CA_Fpoints[chain[1] - 1][2] - CA_Fpoints[chain[1] - 2][2]
	
	#distance between repeated chains
	dx = CA_Fpoints[chain[1] - 1][0] - CA_Ipoints[chain[0] - 1][0]
	dy = CA_Fpoints[chain[1] - 1][1] - CA_Ipoints[chain[0] - 1][1]
	dz = CA_Fpoints[chain[1] - 1][2] - CA_Ipoints[chain[0] - 1][2]
	base_distance = math.sqrt(dx*dx + dy*dy + dz*dz)
	
	#normalize direction vector
	magnitude = math.sqrt(x*x + y*y + z*z)
	x = x/magnitude
	y = y/magnitude
	z = z/magnitude

        number_of_atoms = int(last_line[1]) - int(first_line[1]) + 1 + 1 # second 1 compensating for ter
        
        number_of_chains = ord(last_line[4]) - ord(first_line[4]) + 1
        
	monomerInfo = [x, y, z, number_of_atoms, number_of_chains, base_distance, number]
	#fixFirstChain(rotation_point, x, y, z, filename)
	return monomerInfo

        
        

def computeRotationMatrix(r_point, rotation_axis, theta):
	# To rotate with respect to an axis and a point we do the following:
	
	rotation_point = [r_point[0][0], r_point[1][0], r_point[2][0]]
	#1. Translate point to origin
	
	t = [ [1,0,0,-1*rotation_point[0]], [0,1,0,-1*rotation_point[1]], [0,0,1,-1*rotation_point[2]], [0,0,0,1] ]
	
	#2. Align rotation axis to x,y plane
	xy_vector = [ rotation_axis[0], 0, 0] # project to xy plane
	xy_length = math.sqrt( xy_vector[0] * xy_vector[0] + xy_vector[1] * xy_vector[1] + xy_vector[2] * xy_vector[2] )
	axis_vector = [rotation_axis[0], 0, rotation_axis[2]]
	axis_vector_length = math.sqrt(axis_vector[0]*axis_vector[0] + axis_vector[1]*axis_vector[1] + axis_vector[2]*axis_vector[2])
	xy_angle = math.acos( ( xy_vector[0]*axis_vector[0] + xy_vector[1]*axis_vector[1] + xy_vector[2]*axis_vector[2] ) / (xy_length*axis_vector_length) ) #angle = arccos ( (v1 dot v2) / magnitudes)
	if(xy_angle < 0):
		xy_angle = -1*xy_angle
		
	#if axis is in quadrant 1 or 3, use reverse angle
	if( (rotation_axis[0] > 0 and rotation_axis[2] < 0) or (rotation_axis[0] < 0 and rotation_axis[2] > 0) ):
		xy_angle = -1 * xy_angle
	
	r_xy = [ [math.cos(xy_angle),0,math.sin(xy_angle),0], [0,1,0,0], [-1*math.sin(xy_angle),0,math.cos(xy_angle),0], [0,0,0,1] ]
	
	#3. Align rotation axis to x axis
	x_vector = [ rotation_axis[0], 0, 0 ]
	x_vector_length = math.sqrt( x_vector[0]*x_vector[0] + x_vector[1]*x_vector[1] + x_vector[2]*x_vector[2] )
	if(rotation_axis[0] < 0):
		axis_vector_length = -1*axis_vector_length
	axis_vector = [axis_vector_length, rotation_axis[1], 0] #tricky, eh? the new vector has same x length as the size of the projected vector onto the xz plane of the rotation_vector
	axis_length = math.sqrt( axis_vector[0]*axis_vector[0] + axis_vector[1]*axis_vector[1] + axis_vector[2]*axis_vector[2] )
	x_angle = math.acos( (axis_vector[0]*x_vector[0] + axis_vector[1]*x_vector[1] + axis_vector[2]*x_vector[2]) / (x_vector_length*axis_length) )
	if(x_angle < 0):
		x_angle = -1*x_angle
	#if axis is in quadrant 1 or 3, use reverse angle
	if( (rotation_axis[0] > 0 and rotation_axis[1] > 0) or (rotation_axis[0] < 0 and rotation_axis[1] < 0) ):
                x_angle = -1 * x_angle
	r_x_angle = [ [math.cos(x_angle),-1*math.sin(x_angle),0,0], [math.sin(x_angle),math.cos(x_angle),0,0], [0,0,1,0], [0,0,0,1] ]

	#4. Perform angle rotation on x axis
	r_x = [ [1,0,0,0], [0,math.cos(theta*math.pi/180),-1*math.sin(theta*math.pi/180),0], [0,math.sin(theta*math.pi/180),math.cos(theta*math.pi/180),0], [0,0,0,1]]
	
	#5. inverse of 3
	r_x_angle_inv = [ [math.cos(-1*x_angle),-1*math.sin(-1*x_angle),0,0], [math.sin(-1*x_angle),math.cos(-1*x_angle),0,0], [0,0,1,0], [0,0,0,1] ]

	#6. inverse of 2
	r_xy_inv = [ [math.cos(-1*xy_angle),0,math.sin(-1*xy_angle),0], [0,1,0,0], [-1*math.sin(-1*xy_angle),0,math.cos(-1*xy_angle),0], [0,0,0,1] ]

	#7. inverse of 1
	t_inv = [ [1,0,0,rotation_point[0]], [0,1,0,rotation_point[1]], [0,0,1,rotation_point[2]], [0,0,0,1] ]

        #compute Rotation matrix
	Rotation = matrix.mult(t_inv, matrix.mult(r_xy_inv, matrix.mult(r_x_angle_inv, matrix.mult(r_x, matrix.mult(r_x_angle, matrix.mult(r_xy, t))))))
	return Rotation



main()
