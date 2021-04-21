# encoding=utf-8

import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt


def plot_informationCriterion(aic, bic, n_components):

	plt.plot(n_components, aic, '-b', label='AIC')
	plt.plot(n_components, bic, 'r--', label='BIC')
	# plt.title("Akaike Information Criterion (AIC) vs components")
	plt.ylabel('Information Criterion')
	plt.xlabel('n_components')
	plt.legend()
	plt.savefig('Information Criterion plot for GC6 data.png')
	# plt.tight_layout()
	plt.show()

	pass


def plot_Qualscore_profile():
	pass


def plot_regressionforProfile():
	pass

