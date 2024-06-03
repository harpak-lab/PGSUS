import numpy as np
import pandas as pd 

class component_error(object):

	def __init__(self, simulation_df):

		self.simulation_df = simulation_df

	def run(self):
		print(self.simulation_df.keys())

