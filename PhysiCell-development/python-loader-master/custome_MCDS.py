import numpy as np
import pandas as pd
import scipy
import matplotlib.pyplot as plt 
from pyMCDS import pyMCDS
from scipy.integrate import odeint
import xml.etree.ElementTree as ET
import numpy as np
import pandas as pd
import scipy.io as sio
import sys
import warnings
from pathlib import Path


class cMCDS:
    def __init__(self):
        
    def count_live_cells(self):
        all_type= mcds.data['discrete_cells']['cell_type']
        cycle_models = mcds.data['discrete_cells']['cycle_model']
        number_of_cells = all_type.size
        live1_count = 0
        live2_count = 0
        for n in range(number_of_cells):
            cell_type = np.int( all_type[n] )
            cycle_model = np.int( cycle_models[n] )
    
    # cylce_model = 100 is dead cell
            if(cell_type==1 and cycle_model!=100):
                live1_count +=1
            if(cell_type==2 and cycle_model!=100):
                live2_count +=1
        dead_cell = number_of_cells-live1_count-live2_count
        return live1_count,live2_count
