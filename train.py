#!/usr/bin/env python
import os
import sys
import numpy as np
# import pandas as pd
import xgboost as xgb
import pickle as pk
from xgboost import XGBClassifier, plot_importance, DMatrix
from root_numpy import root2array
from sklearn import model_selection#, preprocessing
from sklearn.metrics import accuracy_score, classification_report, roc_auc_score, roc_curve, auc
import matplotlib.pyplot as plt
from collections import OrderedDict

def get_barename(filestring):
  if '/' in filestring:
    return filestring.split('/')[-1].split('.')[0]
  else:
    return filestring.split('.')[0]

listfile = sys.argv[1]

savename = ''
for filename in open(listfile):
  namelist = get_barename(filename).split('_')
  namelist.remove('anTGCtree')
  savename += '_'.join(namelist)+'_'

signal = 'beamhalo'

save_directory = os.getcwd() + '/' + signal + '/'

# create directory in which to save model
try:  
	os.mkdir(save_directory)
except OSError:  
	print ("Creation of the directory %s failed" % save_directory)
else:  
	print ("Successfully created the directory %s " % save_directory)

# get root file paths for signal and backgrounds
signal_list = [line.rstrip('\n') for line in open(listfile) if signal in line]
print('Signal list: ', signal_list)
bg_list = [line.rstrip('\n') for line in open(listfile) if 'prompt' in line]
print('Background list: ', bg_list)

feature_names = ['phoR9','phoSigmaIEtaIEtaFull5x5','phoSigmaIEtaIPhiFull5x5','phoSigmaIPhiIPhiFull5x5','phoHoverE','phoSeedTime','phoNumHERHzside','phoNumESRHzside','phoSCEtaWidth','phoSCPhiWidth']
feature_symbols = ['R9',r'$\sigma_{i\eta i\eta}$',r'$\sigma_{i\eta i\phi}$',r'$\sigma_{i\phi i\phi}$','H/E','Seed Time','num HERH zside','num ESRH zside',r'SC $\eta$ width',r'SC $\phi$ width']
#feature_names = ['phoR9','phoSigmaIEtaIEtaFull5x5','phoSigmaIEtaIPhiFull5x5','phoSigmaIPhiIPhiFull5x5','phoHoverE','phoNumHERH','phoNumESRH','phoSCEtaWidth','phoSCPhiWidth']
#feature_symbols = ['R9',r'$\sigma_{i\eta i\eta}$',r'$\sigma_{i\eta i\phi}$',r'$\sigma_{i\phi i\phi}$','H/E','num HERH','num ESRH',r'SC $\eta$ width',r'SC $\phi$ width']
#feature_names = ['phoR9','phoSigmaIEtaIEtaFull5x5','phoSigmaIEtaIPhiFull5x5','phoSigmaIPhiIPhiFull5x5','phoHoverE','phoSCEtaWidth','phoSCPhiWidth']
#feature_symbols = ['R9',r'$\sigma_{i\eta i\eta}$',r'$\sigma_{i\eta i\phi}$',r'$\sigma_{i\phi i\phi}$','H/E',r'SC $\eta$ width',r'SC $\phi$ width']
#feature_names = ['phoR9','phoSigmaIEtaIEtaFull5x5','phoSigmaIEtaIPhiFull5x5','phoSigmaIPhiIPhiFull5x5','phoNumHERH','phoSCEtaWidth','phoSCPhiWidth','phoSeedTime']
#feature_symbols = ['R9',r'$\sigma_{i\eta i\eta}$',r'$\sigma_{i\eta i\phi}$',r'$\sigma_{i\phi i\phi}$','num HERH',r'SC $\eta$ width',r'SC $\phi$ width','Seed Time']
feature_dict = {}
for i in range(len(feature_names)):
  feature_dict[feature_names[i]] = feature_symbols[i]

strnumfeat = str(len(feature_names))
savedir = 'BDT_figures_'+savename+strnumfeat+'vars'
if not os.path.exists(savedir):
  os.makedirs(savedir)

# load signal from root files
_sigFeat = root2array(signal_list, 'EventTree', feature_names)
sigFeat = _sigFeat.reshape(_sigFeat.shape + (-1,))
print('# of signal events = ') + str(np.shape(sigFeat)[0])
# load backgrounds from root files
_bkgFeat = root2array(bg_list, 'EventTree', feature_names)
bkgFeat = _bkgFeat.reshape(_bkgFeat.shape + (-1,))
print '# of background events = ' + str(np.shape(bkgFeat)[0])

sigList = []
sigFeat = sigFeat.flatten()
for event in sigFeat:
  if len(event[0]) == 1:
    new_event = []
    for feature in event:
      new_event.append(feature[0])#feature[0]
    sigList.append(tuple(new_event))
  else:
    temp_list = []
    for feature in event:
      feature_list = np.split(feature,len(feature))
      temp_list.append(feature_list)
    for i in range(len(event[0])):
      new_event = []
      for j in range(len(temp_list)):
        new_event.append(temp_list[j][i][0])
      sigList.append(tuple(new_event))
sigFeat = np.array(sigList,dtype='float32')

bkgFeat = bkgFeat.flatten()
bkgList = []
for event in bkgFeat:
  if len(event[0]) == 1:
    new_event = []
    for feature in event:
      new_event.append(feature[0])
    bkgList.append(tuple(new_event))
bkgFeat = np.array(bkgList,dtype='float32')


print 'bkgFeat.shape: '+str(bkgFeat.shape[0])
print 'sigFeat.shape: '+str(sigFeat.shape[0])

if bkgFeat.shape[0] < sigFeat.shape[0]:
  sigFeat = sigFeat[:bkgFeat.shape[0],:]
else:
  bkgFeat = bkgFeat[:sigFeat.shape[0],:]

print 'bkgFeat.shape: '+str(bkgFeat.shape)
print 'sigFeat.shape: '+str(sigFeat.shape)
#plt.hist(bkgFeat[:,0])
#plt.show()
#raw_input("wait")

# data & background labels
sigY = np.full((np.shape(sigFeat)[0], 1), 1, dtype=int)
bkgY = np.full((np.shape(bkgFeat)[0], 1), 0, dtype=int)

# merge signal and background arrays
_Feat = np.vstack((sigFeat, bkgFeat))
_Y = np.vstack((sigY, bkgY))

# split into train & test sets
_test_size = 0.30
X_train, X_test, Y_train, Y_test = model_selection.train_test_split(_Feat, _Y, test_size=_test_size)

#print "X_train: "
#print X_train
#dtrain = DMatrix(X_train, label=Y_train, feature_names=feature_names)
#
#print "Y_train: "
#print Y_train
#
#print "X_test: "
#print X_test
#
#print "Y_test: "
#print Y_test


# fit model to training data
#model = XGBClassifier(base_score=0.5, colsample_bylevel=1, colsample_bytree=1,
#	gamma=0, learning_rate=0.1, max_delta_step=0, max_depth=5,
#	min_child_weight=1, missing=None, n_estimators=500, nthread=-1,
#	objective='binary:logistic', reg_alpha=0, reg_lambda=1,
#	scale_pos_weight=1, seed=0, silent=True, subsample=1)
model = XGBClassifier(max_depth=5, learning_rate=0.1, n_estimators=100,
    objective='binary:logistic', base_score=0.5, colsample_bylevel=1,
    gamma=0, colsample_bytree=1, max_delta_step=0, min_child_weight=1, 
    missing=None, reg_alpha=0, reg_lambda=1, scale_pos_weight=1, seed=0, subsample=1)

model.feature_names = feature_names

print(model)
model.fit(X_train, np.ravel(Y_train))

#save model
model_save_name = save_directory + '/antgc_' + signal + '_bdt'
model._Booster.dump_model(model_save_name+'.xgb')
model._Booster.save_model(model_save_name+'_bin.xgb')
pk.dump(model, open(model_save_name+'.pickle', 'wb'))
print 'Saved model ' + model_save_name + '(*.xgb, *.pickle)'

# save test and test sets
train_save_file = save_directory + '/train_set' + signal + '.txt'
test_save_file = save_directory + '/test_set' + signal + '.txt'
train_save = np.append(X_train, Y_train, axis=1)
test_save = np.append(X_test, Y_test, axis=1)
np.savetxt(train_save_file, train_save, delimiter=",")
np.savetxt(test_save_file, test_save, delimiter=",")
print('Saved train set ',train_save_file)
print('Saved test set ',test_save_file)

# # make predictions for test data
y_test_pred = model.predict(X_test)#,validate_features=True)
y_train_pred = model.predict(X_train)#,validate_features=True)
#predictions = [round(value) for value in y_test_pred]

# accuracy = accuracy_score(Y_test, predictions)
accuracy = model.score(X_test, Y_test)
print("Accuracy: %.2f%%" % (accuracy * 100.0))

plt.figure()
# reorder features by importance with labels
impnums = list(model.feature_importances_)
impnames = model.feature_names
implist = []
for i in range(len(impnames)):
  impdict = {}
  impdict['name'] = impnames[i]
  impdict['value'] = impnums[i]
  implist.append(impdict)
impodict = sorted(implist, key=lambda x: x['value'])
imponames = []
imponums = []
for entry in impodict:
  imponames.append(feature_dict[entry['name']])
  imponums.append(entry['value'])
y_pos = np.arange(len(imponames))
plt.barh(y_pos,imponums,height=0.5,tick_label=imponames)
offsets = list(0.25*np.ones(len(imponames)))
plt.yticks(y_pos+offsets,imponames)
plt.tight_layout()
plt.gcf().savefig(savedir+'/FeatureImportance_'+savename+strnumfeat+'.png')
plt.figure()
# BDT probability score plot
test_ones = np.ones((len(Y_test),1),dtype=int)
train_ones = np.ones((len(Y_train),1),dtype=int)
Y_test_inv = test_ones-Y_test
Y_train_inv = train_ones-Y_train
y_test_proba = model.predict_proba(X_test)[:,1]
y_test_proba_pos = np.compress(Y_test.flatten(),y_test_proba)
y_test_proba_neg = np.compress(Y_test_inv.flatten(),y_test_proba)
y_train_proba = model.predict_proba(X_train)[:,1]
y_train_proba_pos = np.compress(Y_train.flatten(),y_train_proba)
y_train_proba_neg = np.compress(Y_train_inv.flatten(),y_train_proba)
plt.hist(y_test_proba_pos,color='r',histtype='step',label='beam halo test')
plt.hist(y_train_proba_pos,color='b',histtype='step',label='beam halo train')
plt.hist(y_test_proba_neg,color='m',histtype='step',label='prompt test')
plt.hist(y_train_proba_neg,color='c',histtype='step',label='prompt train')
plt.legend(loc='upper center')
plt.title('BDT Score')
plt.xlabel('BDT Score')
plt.gcf().savefig(savedir+'/BDTscore_'+savename+strnumfeat+'.png')
Y_test_flat = Y_test.flatten()
Y_train_flat = Y_train.flatten()
fpr_test, tpr_test, threshold_test = roc_curve(Y_test,y_test_proba)
fpr_train, tpr_train, threshold_train = roc_curve(Y_train,y_train_proba)
roc_auc_test = auc(fpr_test,tpr_test)
roc_auc_train = auc(fpr_train,tpr_train)
plt.figure()
plt.title('Receiver Operator Characteristic')
plt.plot(fpr_test,tpr_test,'r',label='AUC test = {0}'.format(roc_auc_test))
plt.plot(fpr_train,tpr_train,'b',label='AUC train = {0}'.format(roc_auc_train))
plt.legend(loc='lower right')
plt.plot([0,1],[0,1],'k--')
#plt.plot([0.296],[0.983],'o')
#plt.text(0.306,0.953,'AN2010-115-v4-1')
plt.plot([0.00019],[0.044],'o')
plt.text(0.01019,0.044,'beamhalostatus.pdf')
plt.ylabel('True Positive Rate')
plt.xlabel('False Positive Rate')
plt.gcf().savefig(savedir+'/ROCplot_'+savename+strnumfeat+'.png')

plt.figure()
fig, ax = plt.subplots()
plt.title('Receiver Operator Characteristic')
ax.plot(fpr_test,tpr_test,'r',label='AUC test = {0}'.format(roc_auc_test))
rocfile = open(savedir+'/txtROC_'+savename+strnumfeat+'.txt','w')
#for i in range(len(tpr_test)):
#  if abs(fpr_test[i]-0.00019) < 0.0009 and fpr_test[i] > 0:
#    print fpr_test[i]
#    print tpr_test[i]
rocpairs = np.stack((fpr_test,tpr_test),axis=-1)
for pair in rocpairs:
  if pair[0] <= 0.01:
    print pair
rocfile.write('?')
ax.plot(fpr_train,tpr_train,'b',label='AUC train = {0}'.format(roc_auc_train))
ax.legend(loc='lower right')
ax.plot([0,1],[0,1],'k--')
ax.plot([0.00019],[0.044],'o',markersize=10)
ax.text(0.01019,0.044,'beamhalostatus.pdf')
plt.ylabel('True Positive Rate')
plt.xlabel('False Positive Rate')
ax.set_xlim(0,0.005)
ax.set_ylim(0,1)
plt.gcf().savefig(savedir+'/miniROCplot_'+savename+strnumfeat+'.png')
