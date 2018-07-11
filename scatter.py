import matplotlib
import matplotlib.pyplot as plt

import pandas as pd
import numpy as np
import datetime as dt

import csv
import sys
import os


def accumulate(data):
	sum = 0
	for i in range(0, len(data)):
		sum = sum + data[i]
		data[i] = sum

if __name__ == "__main__":
	args = sys.argv
	name_csv = "f_nn_input.txt"
	plt.style.use('ggplot') 
	#font = {'family' : 'meiryo'}
	#matplotlib.rc('font', **font)

	# data import
	# format as C strftime()
	d = pd.read_csv(name_csv, header=None)
	#accumulate(d[index_column])
	#accumulate(d[177])
	#accumulate(d[178])
	#accumulate(d[179])

	plot=d.plot.scatter(x=[0], y=[1], figsize=(16,6), alpha=0.5, rot=20)
	#plot=d.plot(y=[177, 178, 179], figsize=(16,6), alpha=0.5, rot=20)
	#plot.legend(prop={"size": 26})
	#plot_range_rect(plot)
	plot.get_figure().savefig("out.png", bbox_inches="tight")
