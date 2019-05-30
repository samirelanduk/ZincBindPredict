import os
from datetime import datetime
from sklearn import neighbors
import sklearn.metrics as metrics
from sklearn.model_selection import train_test_split
import matplotlib.pyplot as plt
import pandas as pd

MODEL = "KNN"

def print_log(string="", log=f"logs/{MODEL}_{datetime.now().strftime('%Y-%m-%d_%H-%M-%S')}.log"):
    with open(log, "a") as f:
        f.write(string + "\n")
        print(string)


# Define hyperparameters
print_log("HYPERPARAMETERS")
hyperparameters = {
 "n_neighbors": 5
}
for key, value in hyperparameters.items():
    print_log(f"{key}: {value}")
print_log("\n")


# What datasets are there?
csvs = filter(lambda f: f.endswith(".csv"), os.listdir("data"))
for csv in csvs:
    # Get dataset
    print_log("Processing " + csv)
    df = pd.read_csv(f"data/{csv}")
    name = csv[:-4]

    # Get relevant data
    X = df.iloc[:, 2:-1].values
    y = df.iloc[:, -1].values

    # Test/train split
    X_train, X_test, y_train, y_test = train_test_split(
     X, y, test_size=0.2, random_state=42
    )

    # Create model
    model = neighbors.KNeighborsClassifier(**hyperparameters)

    # Fit model
    a = model.fit(X_train, y_train)

    for label, X, y in (("Train", X_train, y_train), ("Test", X_test, y_test)):
        print_log(" " + label.upper())
        # Get predictions
        predictions = model.predict(X)

        # Get probabilities
        probabilities = model.predict_proba(X)
        roc_x, roc_y, _ = metrics.roc_curve(y, probabilities[:, 1])

        # Save chart
        plt.plot(roc_x, roc_y, label="ROC curve")
        plt.plot([0, 1], [0, 1], linestyle="--", label="Random model")
        plt.xlabel("False Positive Rate")
        plt.ylabel("True Positive Rate")
        plt.legend(loc="lower right")
        plt.title(f"{MODEL} {name} ROC")
        plt.savefig(f"{MODEL}-{name}.png")
        plt.clf()

        # Get metrics
        print_log(f"  Recall {metrics.recall_score(y, predictions)}")
        print_log(f"  Precision {metrics.precision_score(y, predictions)}")
        print_log(f"  F1 {metrics.f1_score(y, predictions)}")
        print_log(f"  MCC {metrics.matthews_corrcoef(y, predictions)}")
        print_log(f"  AUC {metrics.auc(roc_x, roc_y)}")
        print_log("")
    print_log("")
