import os
import sys

base_dir = os.path.dirname(__file__)
data_dir = os.path.join(base_dir,"..", "data")
sys.path.extend([os.path.join(base_dir, '../..')])

from enricher.features import *


def main():