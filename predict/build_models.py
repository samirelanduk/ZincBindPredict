#! /usr/bin/env python3

import sys
import os
import pandas as pd
import warnings
import joblib
from sklearn.model_selection import train_test_split
from sklearn.ensemble import RandomForestClassifier

warnings.filterwarnings("ignore", category=FutureWarning)

categories = ["structure", "sequence"]
families = open("data/families.dat").read().splitlines()
for arg in sys.argv:
    if arg.startswith("--categories="):
        types = [c for c in arg[8:].split(",") if c in categories]
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

        # Save model
        joblib.dump(model, f"predict/models/{category}/{family}-randomforest.joblib")

        

