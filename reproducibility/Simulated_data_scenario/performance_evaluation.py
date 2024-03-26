# This file is used to evaluation the performance of the model using the metrics including:
# 1) AUCROC (macro), 2) AUCROC (micro), 3) AUCPR, 4) F1 score, 5) MCC, and 6) ARI
# Input data:
# 1) A file containing truth labels (1st column), and predictions (2nd) columns

from traceback import print_tb
from sklearn.metrics import roc_auc_score
from sklearn.metrics import average_precision_score
from sklearn.metrics import f1_score
from sklearn.metrics import matthews_corrcoef
from sklearn.metrics.cluster import adjusted_rand_score
from sklearn.metrics import confusion_matrix
from sklearn.metrics import classification_report

from sklearn import metrics
import matplotlib.pyplot as plt

import tensorflow as tf
import pandas as pd

#path_file_predictions = "/home/thatchau/demultiplexing_project/artificial_tissue_blood_experiment/prediction_artificial_data/results/TB_annotation_5/prediction_TB_monocyte_TBannotation5.csv"

path_file_predictions = "/home/thatchau/demultiplexing_project/artificial_tissue_blood_experiment/prediction_artificial_data/results/test_originator_function/meta_data_merged_tissue_blood_residents_DEMO_TB_TuN_annotation_usingWholeblood_commonImmune_integerClasses.csv"


df_file_predictions = pd.read_csv(path_file_predictions)

print(df_file_predictions)


list_truth_labels = (df_file_predictions.iloc[:,15]).tolist()

list_predictions = (df_file_predictions.iloc[:,16]).tolist()

list_truth_labels_one_hot = tf.keras.utils.to_categorical(list_truth_labels, num_classes=2)
list_predictions_one_hot = tf.keras.utils.to_categorical(list_predictions, num_classes=2)

# print(list_truth_labels)
# print(list_predictions)
# print()

# print(len(list_truth_labels))
# print(len(list_predictions))

# print(list_truth_labels_one_hot)
# print()
# print(list_predictions_one_hot)
# print(list_predictions_one_hot[0])

# get AUC (macro)
auc_macro = roc_auc_score(list_truth_labels_one_hot, list_predictions_one_hot, average = 'macro')

# get AUC (micro)
auc_micro = roc_auc_score(list_truth_labels_one_hot, list_predictions_one_hot, average = 'micro')

# get AUCPR
auc_pr = average_precision_score(list_truth_labels_one_hot, list_predictions_one_hot, average = 'micro')

# get F1 score
f1 = f1_score(list_truth_labels_one_hot, list_predictions_one_hot, average = 'micro')

# get MCC
mcc = matthews_corrcoef(list_truth_labels, list_predictions)

# get ARI
ari = adjusted_rand_score(list_truth_labels, list_predictions)

print("Working on: ", path_file_predictions)
print()
print("Result:")
print("auc_macro: ", auc_macro)
print("auc_micro: ", auc_micro)
print("auc_pr: ", auc_pr)
print("F1 score: ", f1)
print("MCC: ", mcc)
print("ARI: ", ari)

# Get classification report
print(classification_report(list_truth_labels, list_predictions))


# get confusion matrix in the format (tn, fp, fn, tp)
#tn, fp, fn, tp = confusion_matrix(list_truth_labels, list_predictions).ravel()

#print("#### Confusion matrix ###")
#print("True negative: " + str(tn))
#print("False positive: " + str(fp))
#print("False negative: " + str(fn))
#print("True positive: " + str(tp))

#confusion_matrix = metrics.confusion_matrix(list_truth_labels, list_predictions)

#cm_display = metrics.ConfusionMatrixDisplay(confusion_matrix = confusion_matrix, display_labels = [True, False])

#cm_display.plot()
#plt.show()
