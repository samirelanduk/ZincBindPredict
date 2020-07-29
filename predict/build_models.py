#! /usr/bin/env python3

import os
import sys

categories = {c: {} for c in os.listdir(
    os.path.join("data", "csv")
) if "." not in c}
for category in categories:
    categories[category] = {d[:-4]: {} for d in os.listdir(
        os.path.join("data", "csv", category)
    ) if d.endswith("csv")}

for arg in sys.argv:
    if arg.startswith("--categories="):
        allowed_categories = arg[13:].split(",")
        categories = {
            k: v for k, v in categories.items() if k in allowed_categories
        }
    if arg.startswith("--datasets="):
        for category in categories:
            allowed_datasets = arg[11:].split(",")
            categories[category] = {k: v for k, v in
             categories[category].items() if k in allowed_datasets}


for category in categories:
    print(category)
    for d in categories[category]:
        print(" ", d)



'''import sys
import os
import json
import pandas as pd
import numpy as np
import joblib
from itertools import product
import matplotlib.pyplot as plt
from sklearn.model_selection import train_test_split, cross_validate, learning_curve
from sklearn.metrics import recall_score, precision_score, f1_score, roc_auc_score, roc_curve
from collections import Counter
from models import MODELS

categories = ["structure", "sequence"]
families = open("data/families.dat").read().splitlines()
proportion = 1
loadhp = False
for arg in sys.argv:
    if arg.startswith("--categories="):
        categories = [c for c in arg[13:].split(",") if c in categories]
        continue
    if arg.startswith("--families="):
        families = [f for f in arg[11:].split(",") if f in families]
        continue
    if arg.startswith("--proportion="):
        proportion = float(arg[13:])
        assert 0 < proportion <= 1
        continue
    if arg.startswith("--loadhp"):
        loadhp = True
        continue

results = {}
try:
    with open(os.path.join("predict", "models", "results.json")) as f:
        old_results = json.load(f)
except: old_results = {}

for category in categories:
    print(f"{category} models")
    results[category] = {}

    for family in families:
        print(f"    {family}")
        results[category][family] = {}

        # Get the CSV data as a numpy array
        data = pd.read_csv(
            f"data/csv/{category}/{family}.csv",
            usecols=lambda c: c not in ["id", "site"]
        ).values
        np.random.shuffle(data)
        if proportion != 1: data = data[:int(len(data) * proportion)]

        # Get input and label objects
        X, y = data[:, :-1], data[:, -1]

        # Create test/train split
        X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=23)

        # Go through models
        trained_models = []
        for Model in MODELS:
            model_name = MODELS[Model]["name"]
            print("       ", model_name)
            results[category][family][model_name] = {}

            # What are the potential hyperparemeters?
            hp_space = MODELS[Model]["hyperparameters"]

            # Which hyperparameters give the best validation F1 score?
            hyperparameters = None
            if loadhp:
                try:
                    hyperparameters = old_results[category][family][model_name]["hyperparameters"]
                except: pass
            if hyperparameters:
                print("            Cross-validating with precalculated hyperparameters...")
                model = Model(**hyperparameters)
                scores = cross_validate(model, X_train, y_train, cv=5, scoring=["recall", "precision", "f1"])
                validation_recall = scores["test_recall"].mean()
                validation_precision = scores["test_precision"].mean()
                results[category][family][model_name]["hyperparameters"] = hyperparameters
            else:
                print("            Searching hyperparameter space...")
                hp_results = []
                for combo in product(*hp_space.values()):
                    hyperparameters = dict(zip(hp_space.keys(), combo))
                    model = Model(**hyperparameters)
                    scores = cross_validate(model, X_train, y_train, cv=5, scoring=["recall", "precision", "f1"])
                    hp_results.append((scores["test_f1"].mean(), hyperparameters))
                hyperparameters = max(hp_results, key=lambda h: h[0])[1]
                results[category][family][model_name]["hyperparameters"] = hyperparameters
                validation_recall = scores["test_recall"].mean()
                validation_precision = scores["test_precision"].mean()

            # Create learning curve
            print("            Calculating learning curve...")
            model = Model(**hyperparameters)
            train_sizes = [x for x in [
                0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.08, 0.1, 0.2, 0.3, 0.4, 0.6, 0.8, 1
            ] if 1 / (40 * x) <= proportion]
            scores = learning_curve(model, X_train, y_train, train_sizes=train_sizes, scoring="f1", cv=5, n_jobs=-1)
            plt.plot(scores[0], [v.mean() for v in scores[1]], label="Train F1")
            plt.plot(scores[0], [v.mean() for v in scores[2]], label="Test F1")
            plt.legend()
            plt.xlabel("Training set size")
            plt.ylabel("F1 Score")
            plt.title(f"{family}, {model_name.replace(' ', '-')}")
            plt.savefig(f"predict/models/{category}/{family}-{model_name.replace(' ', '-')}-learning-curve.svg")
            plt.clf()

            # Create AUC plot for these hyperparameters
            print("            ROC/AUC...")
            if "SVM" in model_name: hyperparameters["probability"] = True
            model = Model(**hyperparameters)
            model.fit(X_train, y_train)
            y_pred = model.predict_proba(X_test)[:,1]
            fpr, tpr, thresholds = roc_curve(y_test, y_pred, drop_intermediate=False)
            plt.plot(fpr, tpr, color="#643b78", label=f"ROC Curve (AUC: {round(roc_auc_score(y_test, y_pred), 3)})")
            plt.plot([0, 1], [0, 1], linestyle="--")
            plt.xlabel("False Positive Rate")
            plt.ylabel("True Positive Rate")
            plt.title(f"{family}, {model_name.replace(' ', '-')}")
            plt.legend()
            plt.savefig(f"predict/models/{category}/{family}-{model_name.replace(' ', '-')}-auc.svg")
            plt.clf()

            # Final test with best hyperparameters
            print("            Final tests...")        
            y_pred = model.predict(X_test)
            test_recall = recall_score(y_test, y_pred)
            test_precision = precision_score(y_test, y_pred)
            trained_models.append(model)
            
            # Final train and save
            print("            Saving final model...")
            model = Model(**hyperparameters)
            model.fit(X, y)
            model._validation_recall = validation_recall
            model._validation_precision = validation_precision
            model._test_recall = results[category][family][model_name]["recall"] = test_recall
            model._test_precision = results[category][family][model_name]["precision"] = test_precision
            joblib.dump(model, f"predict/models/{category}/{family}-{model_name.replace(' ', '-')}.joblib")
            if "SVM" in model_name:
                del results[category][family][model_name]["hyperparameters"]["probability"]
            with open(os.path.join("predict", "models", "results.json"), "w") as f:
                json.dump(results, f, indent=4)
        
        # Get ensemble of methods for this family
        print("        Ensemble")
        results[category][family]["ensemble"] = {}
        predictions = [model.predict(X_test) for model in trained_models]
        y_pred = [Counter(l).most_common(1)[0][0] for l in zip(*predictions)]
        results[category][family]["ensemble"]["recall"] = recall_score(y_test, y_pred)
        results[category][family]["ensemble"]["precision"] = precision_score(y_test, y_pred)


print("Final Metrics")
for category in results:
    print(" ", category)
    for family in results[category]:
        print("   ", family)
        for model in results[category][family]:
            print("     ", model)
            for metric in results[category][family][model]:
                if metric != "hyperparameters":
                    print(f"        {metric}: {round(results[category][family][model][metric], 3)}")'''