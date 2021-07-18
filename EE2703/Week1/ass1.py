import sys
CIRCUIT='.circuit'
END='.end'

def spice1(line):
	for l in line:

		tokens=l.split()
		p=len(tokens)
		if (p==4 or (p>4 and tokens[4][0]=='#')):
			element_name=tokens[0]
			n1=tokens[1]
			n2=tokens[2]
			value=tokens[3]
			if n1.isalnum() and n2.isalnum():
				print('node names are alphanumeric')
				for h in reversed(tokens):
					print(h,end=' ')
			else:
				return
		elif p==6 or (p>6 and tokens[6][0]=='#'):
			element_name=tokens[0]
			n1=tokens[1]
			n2=tokens[2]
			n3=tokens[3]
			n4=tokens[4]
			value=tokens[5]
			if n1.isalnum() and n2.isalnum() and n3.isalnum() and n4.isalnum():
				print('node names are alphanumeric')
				for i in reversed(tokens):
					print(i,end=' ')
			else:
				return
		elif l==5 or (l>5 and tokens[5][0]=='#'):
			element_name=tokens[0]
			n1=tokens[1]
			n2=tokens[2]
			V=tokens[3]
			value = tokens[4]
			if n1.isalnum() and n2.isalnum() and n3.isalnum() and n4.isalnum():
				print('node names are alphanumeric')
				for i in reversed(tokens):
					print(i,end=' ')
			else:
				return
		return
    
if len(sys.argv) != 2:
    print('\nUsage: %s <inputfile>' % argv[0])
    exit()




try:
	with open(sys.argv[1]) as f:
		
		line=f.readlines()
		print(line)
		start=-1; end=-2
		for l in line:
			if CIRCUIT == l[:len(CIRCUIT)]:
				start = line.index(line)
			elif END == l[:len(END)]:
				end = line.index(line)
				break
		if start >= end:                
			print('Invalid circuit definition')
			exit(0)
		spice1(line[start+1:end])
except IOError:
	print('invalid file name')
	exit()
