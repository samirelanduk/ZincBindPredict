#! /usr/bin/env python3

import sys
import os
import pandas as pd
import warnings
import joblib
import matplotlib.pyplot as plt
from sklearn.model_selection import train_test_split, cross_validate
from sklearn.metrics import recall_score, precision_score, roc_auc_score
from sklearn.ensemble import RandomForestClassifier

warnings.filterwarnings("ignore", category=FutureWarning)

categories = ["structure", "sequence"]
families = open("data/families.dat").read().splitlines()
for arg in sys.argv:
    if arg.startswith("--categories="):
        categories = [c for c in arg[13:].split(",") if c in categories]
        continue
    if arg.startswith("--families="):
        families = [f for f in arg[11:].split(",") if f in families]
        continue

for category in categories:
    print(f"{category} models")

    for family in families:
        print(f"    {family}")

        # Get data and divide it into training and test data
        all_data = pd.read_csv(
            f"data/csv/{category}/{family}.csv",
            usecols=lambda c: c not in ["id", "site"]
        ).values
        data_train, data_test = train_test_split(all_data, random_state=23)
        
        # Train model
        model = RandomForestClassifier()
        X_train, y_train = data_train[:, :-1], data_train[:, -1].astype("int")
        model.fit(X_train, y_train)

        # Evaluate model
        scores = cross_validate(
            model, X_train, y_train, cv=5, scoring=["recall", "precision", "roc_auc"]
        )
        model.validation_recall_ = scores["test_recall"].mean()
        model.validation_precision_ = scores["test_precision"].mean()
        model.validation_roc_auc_ = scores["test_roc_auc"].mean()

        # Test model
        X_test, y_test = data_test[:, :-1], data_test[:, -1].astype("int")
        y_pred = model.predict(X_test)
        model.test_recall_ = recall_score(y_test, y_pred)
        model.test_precision_ = precision_score(y_test, y_pred)
        model.test_roc_auc_ = roc_auc_score(y_test, y_pred)

        # Save model
        joblib.dump(model, f"predict/models/{category}/{family}-randomforest.joblib")

        

