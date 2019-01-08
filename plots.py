import numpy as np
from sklearn.linear_model import LogisticRegression #Logistic Regression.                     
import matplotlib.pyplot as plt #Matplotlib
from matplotlib.colors import ListedColormap


"""
Notes
----------
This function plots cumulative and principal explained variances.
"""
def plot_explained_variance(var_exp, var_exp_scaled, cum_var_exp, scale, numfig, saveplot):
	plt.figure(numfig)	
	plt.bar( range(1, len(var_exp)+1), var_exp_scaled, color='r', alpha=0.5, align='center',
		   label=str(scale) + r'$\times$' + 'individual explained variance', zorder=5 )
	plt.step( range(1, len(var_exp)+1), cum_var_exp, where='mid', 
		   label='cumulative explained variance', zorder=6 )
	plt.xlabel('Principal Components', fontsize=20)
	plt.ylabel('Explained variance ratio', fontsize=20)	
	plt.legend( loc='best' )
	plt.tight_layout()

	if(saveplot): plt.savefig( 'images/variance_ratio.png' )

	#plt.show()
	plt.close			


"""
Notes
----------
> This function plots class labels information for the purpose of illustration.
> Keep in mind that PCA is an unsupervised technique that doesn't use class label information.
"""	
def plot_labels(XX, yy, yy_pred, numfig, saveplot):
	misclassification = (yy_pred != yy)

	num_celltype = 10
	colors = ['r', 'b', 'g', 'c', 'gold', 'firebrick', 'm', 'y', '#7F3FBF', '#BF7F3F', '#654C4F']
	markers = ['s', 'x', 'o', '*', 'P', '>', '^', '<', 'v', 'd']
	ax = [(0,1), (0,2), (0,3), (0,4), (0,5), (1,2), (1,3), (1,4), (1,5), (1,6)] 

	for index_ax in range(len(ax)):
		plt.figure(numfig)
		for index in range(num_celltype):
			plt.plot( XX[ yy == np.unique(yy)[index], ax[index_ax][0] ], 
					XX[ yy == np.unique(yy)[index], ax[index_ax][1] ], 
					c=colors[index], label=np.unique(yy)[index], marker=markers[index],
					markersize=3.5, ls='', zorder=1 )

			plt.plot( XX[ (yy==np.unique(yy)[index])*misclassification, ax[index_ax][0] ], 
					XX[ (yy==np.unique(yy)[index])*misclassification, ax[index_ax][1] ], 
					c='k', label=np.unique(yy)[index]+'misclassified',
					marker=markers[index], markersize=3.5, ls='', zorder=2 )			

		plt.xlabel('Principal Component ' + str(ax[index_ax][0]), fontsize=20)
		plt.ylabel('Principal Component ' + str(ax[index_ax][1]), fontsize=20)
	
		plt.title('Number of Principal Components = ' + str(XX.shape[1]))	
		plt.legend(loc='center left', bbox_to_anchor=(1, 0.5), fontsize=7 )

		square_plot = True

		if(square_plot): plt.axes().set_aspect('equal', adjustable='datalim')
		else: plt.axis('scaled') 

		plt.tight_layout()
		if(saveplot): 
			plt.savefig( 'images/Labels_PCA_PC' + str(ax[index_ax][0]) + '-PC' + str(ax[index_ax][1]) \
			          + '_numPC_' + str(XX.shape[1]) + '.png' )
		#plt.show()
		plt.close()
		numfig+=1