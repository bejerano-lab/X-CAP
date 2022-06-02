import joblib
import lightgbm as lgb
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from collections import defaultdict
from sklearn.metrics import auc, precision_recall_curve, roc_curve
from annovar_utils import get_value_from_line, get_value_from_info
from variant_utils import Variant

############# X-CAP
def get_xcap_feature_names(filepath):
    header = open(filepath).readline().strip().split('\t')
    feature_names = header[: -1]
    return feature_names

def get_xcap_data(filepath):
    feature_names = get_xcap_feature_names(filepath)
    data = np.genfromtxt(filepath, delimiter='\t', skip_header=1)
    features = data[:, :-1]
    labels = data[:, -1].astype(int)
    return feature_names, features, labels

def get_xcap_predictions(clf, filepath, lightgbm=True):
    _, features, labels = get_xcap_data(filepath)
    if lightgbm:
        predictions = clf.predict(features)
    else: # sklearn GradientBoostingClassifier
        predictions = clf.predict_proba(features)[:, 1]
    return labels, predictions

xcap_test_scores = None
def load_xcap_percentiles(xcap_clf):
    global xcap_test_scores
    features_path = "/cluster/u/rrastogi/ECNN/results/d_original/xcap/featurize_0520/total_test.features"
    _, xcap_test_scores = get_xcap_predictions(xcap_clf, features_path)

def get_xcap_percentile(score):
    assert xcap_test_scores is not None
    num_less_than = (xcap_test_scores < score).sum()
    return 100 * (num_less_than / xcap_test_scores.size)

############## ALoFT
def get_zygosity_from_aloft_row(row):
    chrom = row[0].replace('chr', '')

    details = row[7]
    if "GT=" in details:
        gt = get_value_from_info(details, "GT", str)
        assert(gt in ["1/1", "0/1", "1/0"])
        return gt == "1/1"

    if chrom in ['X', 'Y']:
        return True

    return get_value_from_info(details, "HOMOZYGOUS", lambda x: bool(int(x)))

def get_prediction_from_aloft_row(row):
    homozygous = get_zygosity_from_aloft_row(row)
    score = float(row[13]) + float(row[11]) if homozygous else float(row[11])
    return score, homozygous

def get_aloft_predictions(aloft_results_filepath, stopgains_vcf):
    predictions_map = defaultdict(list)
    labels_map = {}
    zygosities_map = {}

    for line in open(aloft_results_filepath):
        if line.startswith("#"):
            continue
        row = line.strip().split('\t')
        preliftover_index = get_value_from_info(row[7], "PRE_LIFTOVER_INDEX", int)
        score, _ = get_prediction_from_aloft_row(row)
        predictions_map[preliftover_index].append(score)

    preliftover_index = 0
    for line in open(stopgains_vcf):
        if line.startswith("#"):
            continue
        row = line.strip().split('\t')
        label = get_value_from_info(row[7], "LABEL", int)
        zygosity = get_zygosity_from_aloft_row(row)
        labels_map[preliftover_index] = label
        zygosities_map[preliftover_index] = zygosity
        preliftover_index += 1

    predictions, labels, zygosities = [], [], []
    for i in range(len(labels_map)):
        score = np.mean(predictions_map[i]) if i in predictions_map else 0.0
        predictions.append(score)
        labels.append(labels_map[i])
        zygosities.append(zygosities_map[i])

    return np.asarray(labels), np.asarray(predictions), np.asarray(zygosities)

aloft_test_scores = None
def load_aloft_percentiles():
    global aloft_test_scores
    results_path = "/cluster/u/rrastogi/ECNN/results/d_original/aloft/total/aloft_output/hg19.aloft.lof.results"
    hg38_vcf_path = "/cluster/u/rrastogi/ECNN/results/d_original/aloft/total/hg38.vcf"
    _, aloft_test_scores, _ = get_aloft_predictions(results_path, hg38_vcf_path)

def get_aloft_percentile(score):
    num_less_than = (aloft_test_scores < score).sum()
    return 100 * (num_less_than / aloft_test_scores.size)

############## MutPred-LOF
def get_mutpred_predictions(input_filepath, results_filepath):
    labels = [get_value_from_line(line, "LABEL", int) for line in open(input_filepath)]

    predictions = []
    for line in open(results_filepath):
        row = line.strip().split('|')
        assert(len(row) == 3)
        predictions.append(float(row[1]))
    
    return np.asarray(labels), np.asarray(predictions)

mutpred_test_scores = None
def load_mutpred_percentiles():
    global mutpred_test_scores
    input_path = "/cluster/u/rrastogi/ECNN/results/d_original/mutpred/total/subset.tsv"
    results_path = "/cluster/u/rrastogi/ECNN/results/d_original/mutpred/total/results_output.txt"
    _, mutpred_test_scores = get_mutpred_predictions(input_path, results_path)

def get_mutpred_percentile(score):
    num_less_than = (mutpred_test_scores < score).sum()
    return 100 * (num_less_than / mutpred_test_scores.size)

############## Whole-Genome Predictors
def get_wgp_predictions(input_fname, results_fname, column_name):
    predictions_df = pd.read_csv(results_fname, sep='\t', header=0)
    assert column_name in predictions_df.columns
    
    default_value = None if column_name == "Eigen-raw_coding" else 0.0

    predictions = []
    labels = []
    for (_, row), line in zip(predictions_df.iterrows(), open(input_fname)):
        pred = row[column_name]
        pred = float(pred) if pred != "." else default_value
        predictions.append(pred)

        label = get_value_from_line(line, "LABEL", int)
        labels.append(label)
    
    if column_name == "Eigen-raw_coding":
        min_pred = min([p for p in predictions if p != default_value])
        predictions = [pred if pred != default_value else min_pred for pred in predictions]

    return np.asarray(labels), np.asarray(predictions)

def get_cadd_predictions(input_fname, results_fname):
    return get_wgp_predictions(input_fname, results_fname, "CADD_phred")

def get_dann_predictions(input_fname, results_fname):
    return get_wgp_predictions(input_fname, results_fname, "DANN_score")

def get_eigen_predictions(input_fname, results_fname):
    return get_wgp_predictions(input_fname, results_fname, "Eigen-raw_coding")

wgp_input_path = "/cluster/u/rrastogi/ECNN/results/d_original/preprocess/collate/total_test.tsv"
wgp_results_path = "/cluster/u/rrastogi/ECNN/results/d_original/wgp/total/avinput.hg38_multianno.txt"

cadd_test_scores = None
dann_test_scores = None
eigen_test_scores = None

def load_cadd_percentiles():
    global cadd_test_scores
    _, cadd_test_scores = get_cadd_predictions(wgp_input_path, wgp_results_path)

def load_dann_percentiles():
    global dann_test_scores
    _, dann_test_scores = get_dann_predictions(wgp_input_path, wgp_results_path)

def load_eigen_percentiles():
    global eigen_test_scores
    _, eigen_test_scores = get_eigen_predictions(wgp_input_path, wgp_results_path)

def get_cadd_percentile(score):
    num_less_than = (cadd_test_scores < score).sum()
    return 100 * (num_less_than / cadd_test_scores.size)

def get_dann_percentile(score):
    num_less_than = (dann_test_scores < score).sum()
    return 100 * (num_less_than / dann_test_scores.size)

def get_eigen_percentile(score):
    num_less_than = (eigen_test_scores < score).sum()
    return 100 * (num_less_than / eigen_test_scores.size)

############## PLOTS
def display(plt, image_path):
    plt.tight_layout()
    if image_path is not None:
        plt.savefig(image_path)
    else:
        plt.show()

##### PRECISION RECALL CURVE
def getPrcAucFromPredictions(y_true, y_pred):
    precision, recall, _ = precision_recall_curve(y_true, y_pred)
    return auc(recall, precision)

def plotPrCurve(y_true, y_pred, title, image_path=None):
    precision, recall, _ = precision_recall_curve(y_true, y_pred)
    auprc = auc(recall, precision)

    plt.plot(recall, precision)
    plt.xlabel('recall')
    plt.ylabel('precision')
    plt.xlim([0.0, 1.0])
    plt.ylim([0.0, 1.0])
    plt.title(f"{title}: AUC={round(auprc, 3)}")
    display(plt, image_path)

# Categories is a hashmap from category name -> (y_true, y_pred)
def plotPrCurves(categories, title, image_path=None):
    legend = []
    for category in categories:
        y_true, y_pred = categories[category]
        precision, recall, _ = precision_recall_curve(y_true, y_pred)
        auprc = auc(recall, precision)

        plt.plot(recall, precision)
        legend.append("{}: {}".format(category, round(auprc, 3)))

    plt.xlabel('recall')
    plt.ylabel('precision')
    plt.xlim([0.0, 1.0])
    plt.ylim([0.0, 1.0])
    plt.legend(legend)
    plt.title(title)
    display(plt, image_path)

##### RECEIVER OPERATING CHRACTERISTIC CURVE
def getRocAucFromPredictions(y_true, y_pred):
    fpr, tpr, _ = roc_curve(y_true, y_pred)
    return auc(fpr, tpr)

def plotRocCurve(y_true, y_pred, title, image_path=None):
    fpr, tpr, _ = roc_curve(y_true, y_pred)
    auroc = auc(fpr, tpr)

    plt.plot(fpr, tpr)
    plt.xlabel('false positive rate')
    plt.ylabel('true positive rate')
    plt.xlim([0.0, 1.0])
    plt.ylim([0.0, 1.0])
    plt.title(f"{title}: AUC={round(auroc, 3)}")
    display(plt, image_path)

# Categories is a hashmap from category name -> (y_true, y_pred)
def plotRocCurves(categories, title, image_path=None):
    legend = []
    for category in categories:
        y_true, y_pred = categories[category]
        fpr, tpr, _ = roc_curve(y_true, y_pred, drop_intermediate=False)
        auroc = auc(fpr, tpr)

        plt.plot(fpr, tpr)
        legend.append("{}: {}".format(category, round(auroc, 3)))

    plt.xlabel('FPR')
    plt.ylabel('TPR')
    plt.xlim([0.0, 1.0])
    plt.ylim([0.0, 1.0])
    plt.legend(legend)
    plt.title(title)
    display(plt, image_path)

##### HIGH-SENSITIVITY REGION RECEIVER OPERATING CHARACTERISTIC CURVE
def getHsrAucFromPredictions(y_true, y_pred):
    fpr, tpr, _ = roc_curve(y_true, y_pred)
    return getHsrAuc(fpr, tpr)

def getHsrAuc(fpr, tpr):
    # Need to special case if there's no tpr greater than 0.95 aside from 1. Specifically, we need to inject
    # the fpr @ tpr=0.95.
    if tpr[-2] < 0.95:
        m = (tpr[-1] - tpr[-2]) / (fpr[-1] - fpr[-2])
        tpr_des = 0.95
        fpr_des = fpr[-2] + (tpr_des - tpr[-2]) / m
            
        tpr = np.insert(tpr, -1, tpr_des)
        fpr = np.insert(fpr, -1, fpr_des)

    clamped_tpr = [min(0.95, x) for x in tpr]
    auc_diff = auc(fpr, tpr) - auc(fpr, clamped_tpr)
    return 20 * auc_diff

def plotHsrRocCurve(y_true, y_pred, title, image_path=None):
    fpr, tpr, _ = roc_curve(y_true, y_pred)
    hsr_auc = getHsrAuc(fpr, tpr)

    plt.plot(fpr, tpr)
    plt.xlabel('FPR')
    plt.ylabel('TPR')
    plt.xlim([0.0, 1.0])
    plt.ylim([0.95, 1.0])
    plt.title(f"{title}: high-sensitive region AUC={round(hsr_auc, 3)}")
    display(plt, image_path)

# Categories is a hashmap from category name -> (y_true, y_pred)
def plotHsrRocCurves(categories, title, image_path=None):
    legend = []
    for category in categories:
        y_true, y_pred = categories[category]
        fpr, tpr, _ = roc_curve(y_true, y_pred)
        hsr_auc = getHsrAuc(fpr, tpr)

        plt.plot(fpr, tpr)
        legend.append("{}: {}".format(category, round(hsr_auc, 3)))

    plt.xlabel('FPR')
    plt.ylabel('TPR')
    plt.xlim([0.0, 1.0])
    plt.ylim([0.95, 1.0])
    plt.legend(legend)
    plt.title(title)
    display(plt, image_path)

def getHsThreshold(y_true, y_pred):
    fpr, tpr, thresholds = roc_curve(y_true, y_pred)
    for rate, threshold in zip(tpr, thresholds):
        if rate >= 0.95:
            return threshold
    return None

def getHsTnr(y_true, y_pred):
    fpr, tpr, _ = roc_curve(y_true, y_pred)
    for f, t in zip(fpr, tpr):
        if t >= 0.95:
            return 1 - f 
    return None

##### FEATURE IMPORTANCES
def plotFeatureImportances(importances, features, title, image_path=None):
    sorted_indices = np.argsort(importances)[::-1]
    sorted_features = [features[i] for i in sorted_indices]
    sorted_importances = [importances[i] for i in sorted_indices]

    plt.figure()
    y = range(2 * len(sorted_importances), 0, -2)
    plt.barh(y, sorted_importances)
    plt.yticks(y, sorted_features)
    plt.xlabel("relative importance")
    plt.ylabel("features")
    plt.title(title)
    display(plt, image_path)

##### CLASSIFIER DEBUGGING
def plotPredictionDistribution(y_pred, title, bins=20, image_path=None):
    plt.hist(y_pred, bins=bins)
    plt.xlabel("prediction score")
    plt.ylabel("density")
    plt.title(title)
    display(plt, image_path)
