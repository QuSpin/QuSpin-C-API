import numpy as np

s = 3
m=np.array(2,dtype=np.uint8)
mask = np.array(3,dtype=np.uint8)
zero = np.array(0,dtype=np.uint8)
ints = np.arange(16).astype(np.uint8)

for i in ints:
    r = i ^ (-(m<<s) ^ i) & (mask << s)
    print("{:04b} {:04b} {:04b}".format(i,r,m))
