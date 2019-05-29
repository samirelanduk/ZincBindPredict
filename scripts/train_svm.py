import os
from sklearn import svm
import sklearn.metrics as metrics
import matplotlib.pyplot as plt
import pandas as pd

# What datasets are there?
csvs = filter(lambda f: f.endswith(".csv"), os.listdir("data"))

for csv in csvs:
    # Get dataset
    print("Processing", csv)
    df = pd.read_csv(f"data/{csv}")
    name = csv[:-4]

    # Get relevant data
    X = df.iloc[:, 2:-1].values
    y = df.iloc[:, -1].values

    # Define hyperparameters
    hyperparameters = {
     "kernel": "rbf", "gamma": "scale", "class_weight": {1: 5}, "probability": True
    }

    # Create model
    model = svm.SVC(**hyperparameters)

    # Fit model
    a = model.fit(X, y)

    # Get predictions
    predictions = model.predict(X)

    # Get probabilities
    probabilities = model.decision_function(X)
    roc_x, roc_y, _ = metrics.roc_curve(y, probabilities)

    # Save chart
    plt.plot(roc_x, roc_y, label="ROC curve")
    plt.plot([0, 1], [0, 1], linestyle="--", label="Random model")
    plt.xlabel("False Positive Rate")
    plt.ylabel("True Positive Rate")
    plt.legend(loc="lower right")
    plt.title(f"{name} ROC")
    plt.savefig(name + ".png")
    plt.clf()

    # Get metrics
    for key, value in hyperparameters.items():
        print(f"{key}: {value}")
    print()
    print("Recall", metrics.recall_score(y, predictions))
    print("Precision", metrics.precision_score(y, predictions))
    print("F1", metrics.f1_score(y, predictions))
    print("MCC", metrics.matthews_corrcoef(y, predictions))
    print("Precision", metrics.precision_score(y, predictions))
    print("AUC", metrics.auc(roc_x, roc_y))
    print("\n")
