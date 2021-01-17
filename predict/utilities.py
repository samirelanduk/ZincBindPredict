import os
import numpy as np
from sklearn.model_selection import GridSearchCV
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import recall_score, precision_score, f1_score, matthews_corrcoef

def parse_data_args(args):
    """Processes the command line arguments and returns the appropriate
    settings."""

    similarities = [None]
    categories = ["structure", "sequence"]
    with open(os.path.join("data", "families.dat")) as f:
        families = f.read().splitlines()
    sizes = [1]
    for arg in args:
        if arg.startswith("--categories="):
            categories = [c for c in categories if c in arg[13:].split(",")]
        if arg.startswith("--limit="):
            families = [f for f in families if f in arg[8:].split(",")]
        if arg.startswith("--exclude="):
            families = [f for f in families if f not in arg[10:].split(",")]
        if arg.startswith("--clustering="):
            similarities = [None if s == "n" else float(s) for s in arg[13:].split(",")]
        if arg.startswith("--sizes="):
            sizes = [float(s) for s in arg[8:].split(",")]
    return (categories, families, similarities, sizes)


def get_data_subset(df, size):
    """Gets a random subset of a dataframe, keeping positive and negative
    proportions the same."""
    
    positives = df.loc[df["positive"] == 1]
    negatives = df.loc[df["positive"] != 1]
    positives_subset = positives.sample(frac=size)
    negatives_subset = negatives.sample(frac=size)
    return np.concatenate((positives_subset.values, negatives_subset.values))


def train_model(X, y):
    """Trains a random forest model from the training data given, after doing
    cross-validation to get the best hyperparameters."""

    hyperparams = get_ideal_hyperparameters(X, y)
    model = RandomForestClassifier(n_jobs=-1, **hyperparams)
    model.hyperparams_ = hyperparams # For later inspection
    model.fit(X, y)
    return model


def get_ideal_hyperparameters(X, y):
    """Uses cross-validation to get the ideal hyperparameters based on 
    training data given."""

    grid = {
        "n_estimators": [10, 100, 1000], "max_depth": [4, 6, 8, None],
        "criterion": ["gini", "entropy"], "max_features": ["sqrt", "log2"],
    }
    grid_search = GridSearchCV(RandomForestClassifier(), grid, n_jobs=-1, cv=5, scoring="f1")
    grid_search.fit(X, y)
    return grid_search.best_params_


def evaluate_model(model, X, y):
    """Takes a model, test input data, and test output data, and returns various
    metrics of the model using this test data."""

    y_pred = model.predict(X)
    return {
        "tp": len([y1 for y1, y2 in zip(y, y_pred) if y1 == 1 and y2 == 1]),
        "tn": len([y1 for y1, y2 in zip(y, y_pred) if y1 == 0 and y2 == 0]),
        "fp": len([y1 for y1, y2 in zip(y, y_pred) if y1 == 0 and y2 == 1]),
        "fn": len([y1 for y1, y2 in zip(y, y_pred) if y1 == 1 and y2 == 0]),
        "recall": round(recall_score(y, y_pred), 3),
        "precision": round(precision_score(y, y_pred), 3),
        "f1": round(f1_score(y, y_pred), 3),
        "mcc": round(matthews_corrcoef(y, y_pred), 3),
    }