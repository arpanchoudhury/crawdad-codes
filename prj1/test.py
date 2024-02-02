import numpy as np

a = np.zeros((3,3))

b = [(0,1), (2,3), (4,5)]
print(b)
for i in range(len(b)):
	if b[i] == (0,1):
		print(b[i])
