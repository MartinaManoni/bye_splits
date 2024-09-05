_all_ = [ ]

import os
from pathlib import Path
import sys

parent_dir = os.path.abspath(__file__ + 3 * '/..')
sys.path.insert(0, parent_dir)

#from dash import dcc, html
import json
#import dash_daq as daq
import h5py
import re
import plotly.express.colors as px
import random
import logging
import time
import multiprocessing
from shapely.strtree import STRtree


log = logging.getLogger(__name__)

#from bye_splits.plot.display_plotly import yaml, np, pd, go, dbc
from utils import params, common 
from data_handle.data_process import *
from data_handle.geometry import GeometryData
import pandas as pd
import numpy as np
import yaml
import processingMS
import plotMS
import algosMS
import matplotlib
from scipy.spatial.distance import cdist
from shapely import geometry as geom
from shapely.geometry import Polygon, mapping, shape, MultiPolygon, Point
from shapely.ops import unary_union
from matplotlib.patches import Polygon as matPoly
from matplotlib.collections import PatchCollection
from concurrent.futures import ProcessPoolExecutor
from multiprocessing import cpu_count
#matplotlib.use("TkAgg")

class Helper():
    def __init__(self):
        pass

    def read_hdf5_structure(self,file_path, group_name='df', indent=0):

        def _print_group_structure(group, indent):
            for name in group:
                print(" " * indent + name)
                if isinstance(group[name], h5py.Group):
                    _print_group_structure(group[name], indent + 2)

        with h5py.File(file_path, 'r') as f:
            if group_name in f:
                group = f[group_name]
                _print_group_structure(group, indent)

    def read_all_block0_values(self,file_path, group_name='df'):
        with h5py.File(file_path, 'r') as f:
            if group_name in f:
                group = f[group_name]
                if 'block0_values' in group:
                    block0_values = group['block0_values']
                    print("Values in 'block0_values':")
                    for row in block0_values:
                        print(row)
                else:
                    print("Error: 'block0_values' not found in group {}.".format(group_name))
            else:
                print("Error: Group '{}' not found in HDF5 file.".format(group_name))

    def store_event(self, path, data):
        if isinstance(data, dict):
            for key, item in data.items():
                if isinstance(item, pd.DataFrame):
                    item.to_hdf(self.filename, path + str(key))
                else:
                    raise ValueError('Cannot save %s type'%type(item))
        else:
            data.to_hdf(self.filename, path)