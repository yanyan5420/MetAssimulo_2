import pandas as pd
import numpy as np


dd = {["1", "replicate_1"]: [1,2,3],
      ["2", "replicate_2"]: [1,2,3],
      ["3", "replicate_3"]: [1,2,3],}


name_list = [j for i, j in dd.keys()]
print(name_list)