import sys
with open('opacity_AGSS09_z0.02.txt', 'r') as f:
    for line in f:
        if '3.750     |' in line:
            print(line.strip())
            break
