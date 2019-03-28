import os
import csv
import pandas as pd
from itertools import izip

import numpy as np
import matplotlib.pyplot as plt
from itertools import cycle

from sklearn import svm, datasets
from sklearn import metrics

from sklearn.model_selection import train_test_split
from sklearn.preprocessing import label_binarize
from sklearn.multiclass import OneVsRestClassifier
from scipy import interp

class DataAnalysis(object):

    def __init__(self, volume, pdb_li, tract_labels):

        self.volume = volume
        self.pdb_li = pdb_li
        self.tract_labels =tract_labels
        self._all_scores_df()
        self._calculate_medians()


    def _all_scores_df(self):
        for pdb, tract in izip(self.pdb_li, self.tract_labels):
            try:
                target_df = pd.read_csv('pdb_files/{}/{}/tractability_{}/scores.csv'.format(pdb[1:3], pdb, self.volume))
            except IOError:
                print pdb, "FAILED"
                continue

            target_df['pdb'] = pdb
            target_df['tractability'] = tract
            all_li.append(target_df)

        self.all_df = pd.concat(all_li)

        self.all_df.to_csv('scores_validation_{}.csv'.format(self.volume), index=False,)


    def _calculate_medians(self):

        grouped = self.all_df.groupby('pdb', as_index=False).median()
        grouped = grouped.rename(index=str, columns={'scores': 'median_scores'})
        df = self.all_df.merge(grouped, on='pdb')
        df = df.sort_values(by=['median_scores'], ascending=[False])
        df['colour'] = '#D33C60'
        df.loc[(df['tractability'] == 'd'), 'colour'] = '#75D8D5'
        df['prot_id'] = df['tractability'] + '_' + df['pdb']
        self.df = df
        self.df2 = df.groupby('prot_id',sort=False)
        df3 = self.df2.first()
        n_classes = 2
        y_score = df3['median_scores']
        y_test = df3['tractability']
        self.y_test = y_test
        self.y_score = y_score

    def roc(self, return_fig = False):

        fpr, tpr, thresholds = metrics.roc_curve(self.y_test, self.y_score, pos_label='d')
        roc_auc = metrics.auc(fpr, tpr)


        if return_fig:
            plt.figure()
            lw = 2
            plt.plot(fpr, tpr, color='darkorange',
                     lw=lw, label='ROC curve (area = %0.2f)' % roc_auc)
            plt.plot([0, 1], [0, 1], color='navy', lw=lw, linestyle='--')
            plt.xlim([0.0, 1.0])
            plt.ylim([0.0, 1.05])
            plt.xlabel('False Positive Rate')
            plt.ylabel('True Positive Rate')
            plt.title('Tractability from Hotspot Score Distribution')
            plt.legend(loc="lower right")
            return plt
        return roc_auc

    def mcc(self, threshold =14):

        y_pred = []
        if self.y_score is None:
            self.roc()

        for s in self.y_score:
            if s > threshold:
                y_pred.append('d')
            else:
                y_pred.append('n')
        return metrics.matthews_corrcoef(self.y_test, y_pred)

    def acc(self, threshold =14):

        y_pred = []
        if self.y_score is None:
            self.roc()

        for s in self.y_score:
            if s > threshold:
                y_pred.append('d')
            else:
                y_pred.append('n')
        return metrics.accuracy_score(self.y_test, y_pred)


retry_li = ['1moq','4cox', '1lpz', '1vbm','1rv1','1jak']

pdb_li = []
tract_labels = []
with open('info.csv') as csv_file:
    pdb_info = csv.reader(csv_file, delimiter=',')
    all_li = []
    for line in pdb_info:
        pdb = line[0].lower()
        if line[2] == 'v':
            pdb_li.append(pdb)
            tract_labels.append(line[1])

volumes = [350,400,450,500,550]

for v in volumes:
    da = DataAnalysis(v,pdb_li, tract_labels)
    for threshold in xrange(10,20):

        print v,threshold, da.roc(),da.mcc(threshold=threshold), da.acc(threshold)

