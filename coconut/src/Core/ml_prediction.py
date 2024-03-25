"""Summary
"""
################################################################################
#   Copyright (C) 2016-2024 Neelesh Soni <neeleshsoni03@gmail.com>,
#   <neelesh.soni@alumni.iiserpune.ac.in>
#
#   This library is free software: you can redistribute it and/or modify
#   it under the terms of the GNU Lesser General Public License as published by
#   the Free Software Foundation, either version 3 of the License, or
#   (at your option) any later version.
#
#   This library is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#   GNU Lesser General Public License for more details.
#
#   You should have received a copy of the GNU Lesser General Public License
#   along with this library.  If not, see <http://www.gnu.org/licenses/>.
################################################################################

from numpy import array
import numpy as np
from sklearn.model_selection import StratifiedKFold, KFold
from sklearn.ensemble import RandomForestClassifier
from sklearn import preprocessing
from sklearn.preprocessing import label_binarize
from sklearn.feature_selection import VarianceThreshold
from sklearn.metrics import confusion_matrix
from sklearn.metrics import auc
import joblib
from sklearn.metrics import matthews_corrcoef
from sklearn.metrics import classification_report
from sklearn.metrics import accuracy_score
import matplotlib.pyplot as plt
from sklearn.metrics import RocCurveDisplay
from sklearn import metrics
from collections import OrderedDict
from random import shuffle

def load_dataset(fname):
	"""Summary
	
	Args:
	    fname (TYPE): Description
	
	Returns:
	    TYPE: Description
	"""
	inf =open(fname,'r')
	lines = inf.readlines()
	inf.close()

	dataX,dataY,keys=[],[],[]

	for l in lines:
		toks = l.strip().split(' ');
		dataX.append(toks[0:-2])
		keys.append(toks[-2])
		dataY.append(toks[-1])

	dataX,dataY = array(dataX).astype(np.float64),array(dataY).astype(np.int64)
	keys = array(keys)

	return dataX,dataY,keys

def pred_round(X):
	"""Summary
	
	Args:
	    X (TYPE): Description
	
	Returns:
	    TYPE: Description
	"""
	if X<0:
		return -1
	else:
		return 1

def train_RF_models(CCDDIR, CCNETDIR,CCPLOTDIR, DIMERS, KFold_CV=10):
	"""Summary
	
	Args:
	    CCDDIR (TYPE): Description
	    CCNETDIR (TYPE): Description
	    CCPLOTDIR (TYPE): Description
	    DIMERS (TYPE): Description
	    KFold_CV (TYPE): Description
	
	Returns:
	    TYPE: Description
	"""
	DimerLengths=[]
	rs = 1
	X,y,keys = load_dataset(CCDDIR+'/CC_dataset_dimers_RF_scores.txt')
	for ky in keys:
		#Get the hexpairs lengths
		DimerLengths.append(len(DIMERS[ky][3]))
	DimerLengths=np.array(DimerLengths)

	kfoldk = KFold_CV

	# Run classifier with cross-validation and plot ROC curves
	cv = StratifiedKFold(n_splits=kfoldk,shuffle=True,random_state=rs)
	
	AvgMCC=[];AvgACC=[]
	tprs = []; aucs = [];
	mean_fpr = np.linspace(0, 1, 100)

	fig, ax = plt.subplots()

	incorect_pred=OrderedDict()
	class_weight = dict({-1:1.0, 1:1.0})

	i=0;
	for j,(train, test) in enumerate(cv.split(X, y, DimerLengths)):

		clf = RandomForestClassifier(n_estimators=100,bootstrap=True,class_weight=class_weight,n_jobs=8,criterion='entropy',random_state=rs)
		clf = clf.fit(X[train], y[train])
	
		X_predictions=clf.predict(X[test])		
		
		print("Training Data Shape:",X[train].shape, "Testing Data Shape:",X[test].shape)

		print ("Confusion Matrix\n",confusion_matrix(y[test], X_predictions))
		
		mcc=matthews_corrcoef(y[test], X_predictions)
		acc=accuracy_score(y[test], X_predictions)
		AvgMCC.append(mcc);AvgACC.append(acc);
		
		
		print ("MCC:",round(mcc,2),'\t',round(acc,2))
		target_names = ['AP', 'P']
		
		#print(classification_report( y[test], X_predictions, target_names=target_names))

		joblib.dump(clf,CCNETDIR+'/trained_RF_network_'+str(i)+'.pkl')
		i += 1

		ax = plt.gca()
		
		viz = RocCurveDisplay.from_estimator(clf, X[test], y[test],
                         name='ROC fold {}'.format(j),
                         alpha=0.3, lw=1, ax=ax)

		interp_tpr = np.interp(mean_fpr, viz.fpr, viz.tpr)
		interp_tpr[0] = 0.0
		tprs.append(interp_tpr)
		aucs.append(viz.roc_auc)
	
	ax.plot([0, 1], [0, 1], linestyle='--', lw=2, color='r',
        label='Chance', alpha=.8)

	mean_tpr = np.mean(tprs, axis=0)
	mean_tpr[-1] = 1.0
	mean_auc = auc(mean_fpr, mean_tpr)
	
	std_auc = np.std(aucs)
	ax.plot(mean_fpr, mean_tpr, color='b',
        label=r'Mean ROC (AUC = %0.2f $\pm$ %0.2f)' % (mean_auc, std_auc),
        lw=2, alpha=.8)

	std_tpr = np.std(tprs, axis=0)
	tprs_upper = np.minimum(mean_tpr + std_tpr, 1)
	tprs_lower = np.maximum(mean_tpr - std_tpr, 0)
	ax.fill_between(mean_fpr, tprs_lower, tprs_upper, color='grey', alpha=.2,
                label=r'$\pm$ 1 std. dev.')

	ax.set(xlim=[-0.05, 1.05], ylim=[-0.05, 1.05],
       title="Receiver operating characteristic (ROC) Plot")
	ax.legend(loc="lower right")
	plt.savefig(CCPLOTDIR+"/ROC_Curve_"+str(kfoldk)+"_fold"+".pdf")
	#plt.show()
	plt.tight_layout()
	print(np.mean(AvgMCC),np.mean(AvgACC))
	
	return

def train_RF_models_plot(CCDDIR, CCNETDIR,CCPLOTDIR, DIMERS,model_obj):
	"""Summary
	
	Args:
	    CCDDIR (TYPE): Description
	    CCNETDIR (TYPE): Description
	    CCPLOTDIR (TYPE): Description
	    DIMERS (TYPE): Description
	    model_obj (TYPE): Description
	
	Returns:
	    TYPE: Description
	"""
	DimerLengths=[]
	rs = 1
	X,y,keys = load_dataset(CCDDIR+'/CC_dataset_dimers_RF_scores.txt')
	for ky in keys:
		#Get the hexpairs lengths
		DimerLengths.append(len(DIMERS[ky][3]))
	DimerLengths=np.array(DimerLengths)

	model_obj.load_RF_models()
	
	class_weight = dict({-1:1.0, 1:1.0})

	kfoldk = 10
	# Run classifier with cross-validation and plot ROC curves
	cv = StratifiedKFold(n_splits=kfoldk,shuffle=True,random_state=rs)

	idxi=0;
	PredAccuracy=[]
	PredAccuracy_X=[]
	PredAccuracy_S=[]
	Pred_accu_dimerlengthwise = [];dimerlengthsdict={}
	Problist=[]
	for j,(train, test) in enumerate(cv.split(X, y, DimerLengths)):

		clf = RandomForestClassifier(n_estimators=100,bootstrap=True,class_weight=class_weight,n_jobs=8,criterion='entropy',random_state=rs)
		clf = clf.fit(X[train], y[train])
	
		X_predictions=clf.predict(X[test])
		X_predictions_prob=clf.predict_proba(X[test])

		############PLOTTING###############################
		Correct_idxs = np.where(y[test]==X_predictions)
		Incorrect_idxs = np.where(y[test]!=X_predictions)
		idxs_P = np.where(y[test]==1)
		idxs_AP = np.where(y[test]==-1)		

		Correct_Prob_AP = X_predictions_prob[np.intersect1d(Correct_idxs,idxs_AP),0]
		Correct_Prob_P = X_predictions_prob[np.intersect1d(Correct_idxs,idxs_P),1]

		
		Incorrect_Prob_AP = X_predictions_prob[np.intersect1d(Incorrect_idxs,idxs_AP),0]
		Incorrect_Prob_P = X_predictions_prob[np.intersect1d(Incorrect_idxs,idxs_P),1]

		Problist.append([Correct_Prob_AP,Correct_Prob_P, Incorrect_Prob_AP, Incorrect_Prob_P ])
		
		for dlen in sorted(list(set(DimerLengths[test]))):
			Dimers_dlen = np.where(DimerLengths[test]==dlen)
			num_correct_preds = len(np.intersect1d(Dimers_dlen,Correct_idxs))
			num_incorrect_preds =  len(np.intersect1d(Dimers_dlen,Incorrect_idxs))
			pred_acc = round(100*num_correct_preds/(num_correct_preds+num_incorrect_preds),1)
			
			if dlen not in dimerlengthsdict.keys():
				dimerlengthsdict[dlen]=[]	
			dimerlengthsdict[dlen].append([pred_acc,num_correct_preds + num_incorrect_preds])
			

	dimerlengthsdict_means=OrderedDict()
	for k, v in dimerlengthsdict.items():

		dacc = 0;dtot_dimers = 0;
		for pred in v:
			dacc+=pred[0]*pred[1]
			dtot_dimers+=pred[1]

		dimerlengthsdict_means[k]=[dacc/dtot_dimers,dtot_dimers]


	dimerlengthsdict_means = dict(sorted(dimerlengthsdict_means.items(), key=lambda x: x[0]))
	sizecutoff=1
	PredAccuracy=[];PredAccuracy_X=[];PredAccuracy_S=[]
	PredAccuracy1=[];PredAccuracy_X1=[];PredAccuracy_S1=[]
	for k,v in dimerlengthsdict_means.items():

		if v[1]<=sizecutoff:
			continue

		PredAccuracy_X.append(k)
		PredAccuracy.append(v[0])
		PredAccuracy_S.append(v[1])

	
	fig = plt.figure(figsize=(6,4))
	plt.plot(PredAccuracy_X,PredAccuracy)
	plt.scatter(PredAccuracy_X,PredAccuracy,c="red",s=10*np.sqrt(PredAccuracy_S) )
	
	for i, txt in enumerate(PredAccuracy_S):
		plt.annotate(int(txt), (PredAccuracy_X[i]+0.2,PredAccuracy[i]))
	
	plt.xlabel("Number of Hex-Pairs in a Dimer")
	plt.ylabel("Prediction Accuracy (%)")
	plt.title("Prediction accuracy across dimer-length")
	lngd = plt.legend(("Prediction Accuracy","Number of Dimers"))
	lngd.legendHandles[1].set_sizes([20.0])
	plt.savefig(CCPLOTDIR+'/Prediction_accuracy_dimer_length.pdf')


	Correct_Prob_AP1,Correct_Prob_P1, Incorrect_Prob_AP1, Incorrect_Prob_P1 = Problist[0]
	for l in range(1,len(Problist)):
		Correct_Prob_AP1 = np.hstack((Correct_Prob_AP1,Problist[l][0]))
		Correct_Prob_P1 = np.hstack((Correct_Prob_P1,Problist[l][1]))
		Incorrect_Prob_AP1 = np.hstack((Incorrect_Prob_AP1,Problist[l][2]))
		Incorrect_Prob_P1 = np.hstack((Incorrect_Prob_P1,Problist[l][3]))


	fig, ax = plt.subplots(2,1, figsize=(12,6))
	APPreds = np.hstack((Correct_Prob_AP1,Incorrect_Prob_AP1))
	shuffle(APPreds)
	PPreds = np.hstack((Correct_Prob_P1,Incorrect_Prob_P1))
	shuffle(PPreds)

	ax[0].hlines(0.5,-20,len(APPreds),colors='red')
	sc1 = ax[0].scatter(range(0,len(APPreds)),APPreds,s=20,c=APPreds)
	ax[0].set_xlabel("Anti-Parallel Dimers")
	ax[0].set_ylabel("Probability")
	
	ax[0].annotate(len(Correct_Prob_AP1),(len(APPreds),0.7))
	ax[0].annotate(len(Incorrect_Prob_AP1),(len(APPreds),0.3))

	ax[1].hlines(0.5,-20,len(PPreds),colors='red')
	sc2 = ax[1].scatter(range(0,len(PPreds)),PPreds,s=20,c=PPreds)
	ax[1].set_xlabel("Parallel Dimers")
	ax[1].set_ylabel("Probability")
	
	ax[1].annotate(len(Correct_Prob_P1),(len(PPreds),0.7))
	ax[1].annotate(len(Incorrect_Prob_P1),(len(PPreds),0.3))

	plt.colorbar(sc1,ax=ax[0])
	plt.colorbar(sc2,ax=ax[1])
	plt.title("Dimer Probabilities")
	plt.tight_layout()
	plt.savefig(CCPLOTDIR+'/Incorrect_prediction_Probabiities.pdf')
	
	return

def train_bayesian_models(CCDDIR, CCNETDIR, CCPLOTDIR, DIMERS):
	"""Summary
	
	Args:
	    CCDDIR (TYPE): Description
	    CCNETDIR (TYPE): Description
	    CCPLOTDIR (TYPE): Description
	    DIMERS (TYPE): Description
	
	Returns:
	    TYPE: Description
	"""
	DimerLengths=[]
	rs = 1
	X,y,keys = load_dataset(CCDDIR+'/CC_dataset_dimers_Bayes_probabilities.txt')
	for ky in keys:
		#Get the hexpairs lengths
		DimerLengths.append(len(DIMERS[ky][3]))
	DimerLengths=np.array(DimerLengths)

	X_predictions=np.zeros(X.shape[0])
	for i,val in enumerate(X):
		if val[0]>=val[1]:
			X_predictions[i]=-1
		else:
			X_predictions[i]=1


	print ("Confusion Matrix\n",confusion_matrix(y, X_predictions))
	mcc=matthews_corrcoef(y, X_predictions)
	acc=accuracy_score(y, X_predictions)
	print ("MCC:",round(mcc,4),'\t',round(acc,4))
	target_names = ['AP', 'P']
	#print(classification_report( y[test], X_predictions, target_names=target_names))

	############PLOTTING###############################
	Correct_idxs = np.where(y==X_predictions)
	Incorrect_idxs = np.where(y!=X_predictions)
	idxs_P = np.where(y==1)
	idxs_AP = np.where(y==-1)

	Correct_Prob_AP = X[np.intersect1d(Correct_idxs,idxs_AP),0]
	Correct_Prob_P = X[np.intersect1d(Correct_idxs,idxs_P),1]
	
	Incorrect_Prob_AP = X[np.intersect1d(Incorrect_idxs,idxs_AP),0]
	Incorrect_Prob_P = X[np.intersect1d(Incorrect_idxs,idxs_P),1]

	Pred_accu_dimerlengthwise = []
	for dlen in sorted(list(set(DimerLengths))):
		Dimers_dlen = np.where(DimerLengths==dlen)
		num_correct_preds = len(np.intersect1d(Dimers_dlen,Correct_idxs))
		num_incorrect_preds =  len(np.intersect1d(Dimers_dlen,Incorrect_idxs))
		pred_acc = round(100*num_correct_preds/(num_correct_preds+num_incorrect_preds),1)
		
		#Consider only when total dimers are atleast 5	
		if num_correct_preds + num_incorrect_preds >10:
			Pred_accu_dimerlengthwise.append([dlen,pred_acc,num_correct_preds + num_incorrect_preds])
	Pred_accu_dimerlengthwise = np.array(Pred_accu_dimerlengthwise)
	
	fig = plt.figure(figsize=(6,4))
	plt.plot(Pred_accu_dimerlengthwise[:,0],Pred_accu_dimerlengthwise[:,1])
	plt.scatter(Pred_accu_dimerlengthwise[:,0],Pred_accu_dimerlengthwise[:,1],c="red",s=10*np.sqrt(Pred_accu_dimerlengthwise[:,2]) )
	for i, txt in enumerate(Pred_accu_dimerlengthwise[:,2]):
		plt.annotate(int(txt), (Pred_accu_dimerlengthwise[i,0]+0.2,Pred_accu_dimerlengthwise[i,1]))
	
	plt.xlabel("Number of Hex-Pairs in a Dimer")
	plt.ylabel("Prediction Accuracy (%)")
	plt.title("Prediction accuracy across dimer-length")
	lngd = plt.legend(("Prediction Accuracy","Number of Dimers"))
	lngd.legendHandles[1].set_sizes([20.0])
	plt.savefig(CCPLOTDIR+'/Prediction_accuracy_dimer_length.pdf')


	fig, ax = plt.subplots(2,1, figsize=(12,6))
	APPreds = np.hstack((Correct_Prob_AP,Incorrect_Prob_AP))
	shuffle(APPreds)
	PPreds = np.hstack((Correct_Prob_P,Incorrect_Prob_P))
	shuffle(PPreds)

	ax[0].hlines(0.5,-20,len(APPreds),colors='red')
	sc1 = ax[0].scatter(range(0,len(APPreds)),APPreds,s=20,c=APPreds)
	ax[0].set_xlabel("Anti-Parallel Dimers")
	ax[0].set_ylabel("Probability")
	
	ax[0].annotate(len(Correct_Prob_AP),(len(APPreds),0.7))
	ax[0].annotate(len(Incorrect_Prob_AP),(len(APPreds),0.3))

	ax[1].hlines(0.5,-20,len(PPreds),colors='red')
	sc2 = ax[1].scatter(range(0,len(PPreds)),PPreds,s=20,c=PPreds)
	ax[1].set_xlabel("Parallel Dimers")
	ax[1].set_ylabel("Probability")
	
	ax[1].annotate(len(Correct_Prob_P),(len(PPreds),0.7))
	ax[1].annotate(len(Incorrect_Prob_P),(len(PPreds),0.3))

	plt.colorbar(sc1,ax=ax[0])
	plt.colorbar(sc2,ax=ax[1])
	plt.title("Dimer Probabilities")
	plt.tight_layout()
	plt.savefig(CCPLOTDIR+'/Incorrect_prediction_Probabiities.pdf')

	plt.show()

	

	return

if __name__ == '__main__':
	train_RF_models()
	train_RF_models_plot()
	train_bayesian_models()