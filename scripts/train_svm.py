import os
from sklearn import svm
from sklearn.metrics import recall_score, precision_score
import pandas as pd

# What datasets are there?
csvs = filter(lambda f: f.endswith("csv"), os.listdir("data"))

for csv in csvs:
    # Get dataset
    print("Processing", csv)
    df = pd.read_csv(f"data/{csv}")

    # Get relevant data
    X = df.iloc[:, 2:-1].values
    y = df.iloc[:, -1].values

    # Create model
    model = svm.SVC(gamma="scale", class_weight={1: 50})

    # Fit model
    model.fit(X, y)

    # Get predictions
    predicitions = model.predict(X)

    print("Recall", recall_score(y, predicitions))
    print("Precision", precision_score(y, predicitions))