
import os
from sklearn import svm
import sklearn.metrics as metrics
from sklearn.preprocessing import label_binarize
import pandas as pd

# What datasets are there?
csvs = filter(lambda f: f.endswith(".csv"), os.listdir("data"))

for csv in csvs:
    # Get dataset
    print("Processing", csv)
    df = pd.read_csv(f"data/{csv}")

    # Get relevant data
    X = df.iloc[:, 2:-1].values
    y = df.iloc[:, -1].values

    # Create model
    model = svm.SVC(kernel="rbf", gamma="scale", class_weight={1: 8}, probability=True)

    # Fit model
    a = model.fit(X, y)

    # Get predictions
    predictions = model.predict(X)

    # Get probabilities
    probabilities = model.decision_function(X)
    roc_x, roc_y, _ = metrics.roc_curve(y, probabilities)

    # Get metrics
    print("Recall", metrics.recall_score(y, predictions))
    print("Precision", metrics.precision_score(y, predictions))
    print("AUC", metrics.auc(roc_x, roc_y))
