# This file is used to generate a confusion matrix 

# Input data:
# 1) A file containing truth labels and predictions

from traceback import print_tb
from sklearn.metrics import roc_auc_score
from sklearn.metrics import average_precision_score
from sklearn.metrics import f1_score
from sklearn.metrics import matthews_corrcoef
from sklearn.metrics.cluster import adjusted_rand_score
from sklearn.metrics import confusion_matrix

from sklearn import metrics
import matplotlib.pyplot as plt

import tensorflow as tf
import pandas as pd

path_file_predictions = "/home/thatchau/demultiplexing_project/artificial_tissue_blood_experiment/prediction_artificial_data/results/TB_annotation_1/prediction_TB_combined_TBannotation1.csv"
#path_file_predictions = "/home/thatchau/demultiplexing_project/artificial_tissue_blood_experiment/prediction_artificial_data/results/prediction_TB_combined_ALL5executions.csv"

df_file_predictions = pd.read_csv(path_file_predictions)

print(df_file_predictions)


list_truth_labels = (df_file_predictions.iloc[:,6]).tolist()

list_predictions = (df_file_predictions.iloc[:,5]).tolist()

# get confusion matrix in the format (tn, fp, fn, tp)
tn, fp, fn, tp = confusion_matrix(list_truth_labels, list_predictions).ravel()

print("#### Confusion matrix ###")
print("True negative: " + str(tn))
print("False positive: " + str(fp))
print("False negative: " + str(fn))
print("True positive: " + str(tp))

confusion_matrix = metrics.confusion_matrix(list_truth_labels, list_predictions)

cm_display = metrics.ConfusionMatrixDisplay(confusion_matrix = confusion_matrix, display_labels = [False, True])

cm_display.plot()
plt.show()

