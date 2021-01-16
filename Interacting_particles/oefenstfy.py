import matplotlib.pyplot as plt
import math
import numpy as np 

lijst = [1,2,3,4,5,6]

for i in range(0, len(lijst)):
    lijst[i] = 0

## Special periodic conditions, making sure particle stays in the box
# Particle bounces from the side 
def Bouncing_from_sides(nieuw, max, min, box):
    while nieuw > min or nieuw < max:
        if nieuw > max: 
            nieuw = 2 * max - nieuw
        if nieuw < min: 
            nieuw = 2 * min - nieuw
        print(nieuw)
    return nieuw

nieuw = Bouncing_from_sides(4,3,-3, 6)
print (nieuw)
