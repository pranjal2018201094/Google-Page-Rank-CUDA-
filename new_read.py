f= open("filename.txt","r")
f2 = open("input2.txt","w+")
# for i in range(10):
fl =f.readlines()
for x in fl:
	temp = ' '.join(x.split())
	f2.write(temp+'\n')

f2.close()
f.close() 