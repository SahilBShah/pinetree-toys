from os_orig_genome import call_pt
import os
import subprocess

i = 1
#x = os.system('~/pinetree_test.py')

for i in range(100000):
    call_pt.pt_call()
    print("i =", i)
