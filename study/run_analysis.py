import numpy as np
import pandas as pd

def readData(filename):
    data = pd.read_csv(filename)
    print(data)
    return 

if __name__ == "__main__":
    
    filename = "cylinder_399K_deg4_8cpu.csv"
    readData(filename)
