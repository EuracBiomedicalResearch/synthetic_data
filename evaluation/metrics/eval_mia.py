import argparse
import os
import numpy as np
from pysnptools.snpreader import Bed
from sklearn.decomposition import PCA
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import (
    roc_auc_score, roc_curve, precision_recall_curve,
    precision_score, recall_score, f1_score, confusion_matrix
)
from sklearn.model_selection import train_test_split
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from imblearn.under_sampling import RandomUnderSampler
from imblearn.over_sampling import RandomOverSampler
from imblearn.over_sampling import SMOTE
import seaborn as sns


def _str2bool(v):
    if isinstance(v, bool):
        return v
    if v.lower() in ("yes", "true", "t", "1", "y"):  # noqa: E741
        return True
    if v.lower() in ("no", "false", "f", "0", "n"):
        return False
    raise argparse.ArgumentTypeError("Boolean value expected.")


parser = argparse.ArgumentParser(description="Membership Attack Inference (MAI) evaluation")
parser.add_argument("--orig_prefix", required=True, help="Prefix to original PLINK files (without extension)")
parser.add_argument("--gen_prefix", required=True, help="Prefix to generated PLINK files (without extension)")
parser.add_argument("--out_dir", required=True, help="Directory to save output figures")
parser.add_argument("--balancing", choices=["none", "undersample", "oversample"], default="undersample")
parser.add_argument("--pca_components", type=int, default=20)
parser.add_argument("--train_fraction", type=float, default=0.7)
parser.add_argument("--shuffle", type=_str2bool, default=True)
parser.add_argument("--random_state", type=int, default=42)
args = parser.parse_args()

# Ensure output directory exists
out_dir = args.out_dir
os.makedirs(out_dir, exist_ok=True)

# ------------------------------
# Load PLINK data
# ------------------------------
orig_prefix = args.orig_prefix  # without .bed extension
gen_prefix = args.gen_prefix  # without .bed extension

orig_data = Bed(orig_prefix, count_A1=True).read()
gen_data = Bed(gen_prefix, count_A1=True).read()

X_orig = orig_data.val
X_gen = gen_data.val

print(f"Original shape: {X_orig.shape}, Generated shape: {X_gen.shape}")


# ------------------------------
# Impute missing values
# ------------------------------
def impute_missing(X):
    col_means = np.nanmean(X, axis=0)
    inds = np.where(np.isnan(X))

    X[inds] = np.take(col_means, inds[1])
    return X


X_orig = impute_missing(X_orig)

X_gen = impute_missing(X_gen)

# ------------------------------
# Build dataset (imbalanced)
# ------------------------------
X = np.vstack([X_orig, X_gen])
y = np.array([1] * X_orig.shape[0] + [0] * X_gen.shape[0])  # 1=in-sample, 0=out-sample


# ------------------------------
# PCA for dimensionality reduction
# ------------------------------
pca = PCA(n_components=args.pca_components, random_state=args.random_state)
X_pca = pca.fit_transform(X)


# ------------------------------
# Choose balancing strategy
# Options: "none", "undersample", "oversample"
balancing = args.balancing

if balancing == "undersample":
    rus = RandomUnderSampler(random_state=args.random_state)
    X_res, y_res = rus.fit_resample(X_pca, y)
elif balancing == "oversample":
    ros = SMOTE(random_state=args.random_state)
    X_res, y_res = ros.fit_resample(X_pca, y)
else:
    X_res, y_res = X_pca, y

print(f"Resampled shapes: {X_res.shape}, {y_res.shape}")


# ------------------------------
# Train adversary
# ------------------------------
test_size = 1.0 - float(args.train_fraction)
stratify_param = y_res if args.shuffle else None
X_train, X_test, y_train, y_test = train_test_split(
    X_res,
    y_res,
    test_size=test_size,
    shuffle=args.shuffle,
    stratify=stratify_param,
    random_state=args.random_state,
)

clf = LogisticRegression(max_iter=1000)
clf.fit(X_train, y_train)

y_pred_proba = clf.predict_proba(X_test)[:, 1]
y_pred = (y_pred_proba >= 0.5).astype(int)

# ------------------------------
# Metrics
# ------------------------------
auc = roc_auc_score(y_test, y_pred_proba)
prec = precision_score(y_test, y_pred)
rec = recall_score(y_test, y_pred)
f1 = f1_score(y_test, y_pred)
cm = confusion_matrix(y_test, y_pred)

print(f"AUC: {auc:.3f}")
print(f"Precision: {prec:.3f}")
print(f"Recall: {rec:.3f}")
print(f"F1-score: {f1:.3f}")
print("Confusion Matrix:\n", cm)

# ------------------------------
# Figures for paper
# ------------------------------
# ROC Curve
fpr, tpr, _ = roc_curve(y_test, y_pred_proba)
plt.figure()
plt.plot(fpr, tpr, label=f"AUC = {auc:.2f}")
plt.plot([0, 1], [0, 1], "k--")
plt.xlabel("False Positive Rate")
plt.ylabel("True Positive Rate")
plt.title(f"ROC Curve ({balancing})")
plt.legend()
plt.savefig(os.path.join(out_dir, f"roc_curve_{balancing}.png"), dpi=300)
plt.close()

# Precision-Recall Curve
prec_curve, rec_curve, _ = precision_recall_curve(y_test, y_pred_proba)
plt.figure()
plt.plot(rec_curve, prec_curve)
plt.xlabel("Recall")
plt.ylabel("Precision")
plt.title(f"Precision-Recall Curve ({balancing})")
plt.savefig(os.path.join(out_dir, f"pr_curve_{balancing}.png"), dpi=300)
plt.close()


# Plotting the confusion matrix
labels = ["Generated", "Original"]
plt.figure()
sns.heatmap(cm, annot=True, fmt='d', cmap='Blues',
            xticklabels=labels, yticklabels=labels)
plt.xlabel('Predicted Labels')
plt.ylabel('True Labels')
plt.title('Confusion Matrix')
plt.savefig(os.path.join(out_dir, f"test_confusion_matrix_{balancing}.png"), dpi=300)


