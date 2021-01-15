# encoding=utf-8

import matplotlib.pyplot as plt
import seaborn as sns


def plot_FLmodel_output(true_data, fitted_data):

	sns.set_style('darkgrid')
	sns.kdeplot(true_data, color='blue', shade=True, label='True Density', alpha=0.5)
	sns.kdeplot(fitted_data, color='red',shade=True, label='Estimated Density', alpha=0.5)
	plt.xlabel('Fragment-size (log-scale)')
	plt.ylabel('Probability Density')
	plt.legend()
	plt.savefig('Raw GC6 dataset.png')
	plt.show()


def plot_aic_and_bic(n_components, aic, bic):

	plt.plot(n_components, aic, '-b', label='AIC')
	plt.plot(n_components, bic, 'r--', label='BIC')
	plt.title("Akaike Information Criterion (AIC) vs components")
	plt.ylabel('Information Criterion')
	plt.xlabel('n_components')
	plt.savefig('AIC and BIC.png')
	plt.show()


def plot_regression():
	pass


def plot_pca():
	pass


