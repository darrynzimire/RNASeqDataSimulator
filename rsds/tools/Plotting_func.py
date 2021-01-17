# encoding=utf-8

import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
from sklearn.linear_model import LinearRegression
from sklearn.metrics import mean_squared_error
from sklearn.metrics import mean_squared_error, r2_score


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


def plot_regression(dfref, dfsim):
	merged_set = pd.merge(dfref, dfsim, how='inner', on='transcript_id')
	merged_set.head(5)
	merged_set = merged_set.iloc[:, [0, 3, 4, 10, 11]]
	merged_set.head(5)
	x = merged_set.iloc[:, [1]].round()
	y = merged_set.iloc[:, [4]].round()

	# plt.scatter(x, y)
	# data points
	# Model initialization
	regression_model = LinearRegression()
	# # Fit the data(train the model)
	regression_model.fit(x, y)
	# # Predict
	y_predicted = regression_model.predict(x)

	# # model evaluation
	rmse = mean_squared_error(y, y_predicted)
	r2 = r2_score(y, y_predicted)

	# # printing values
	print('Slope:', regression_model.coef_)
	print('Intercept:', regression_model.intercept_)
	print('Root mean squared error: ', rmse)
	print('R2 score: ', r2)

	# plotting values
	plt.scatter(x, y, s=10)
	plt.xlabel('real data')
	plt.ylabel('simulated data')

	# # predicted values
	plt.plot(x, y_predicted, color='r')
	plt.savefig('real vs simulated expression profile.png')
	plt.show()


# plt.show()


def plot_pca():
	pass


