import os
import subprocess

class call_pt:
    def pt_call():
        subprocess.check_call(["python3.7", "/home/sahil/pinetree-toys/three_genes_evolution/memory_leak_tests/pinetree_test.py"])
        os.system('~/pinetree-toys/three_genes_evolution/memory_leak_tests/pinetree_test.py')
