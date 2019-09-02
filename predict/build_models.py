import os
import sys
import pandas as pd
import numpy as np
from sklearn.model_selection import train_test_split
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import recall_score, precision_score
from sklearn.model_selection import cross_validate
import joblib
from config import STRUCTURE_CONFIG, SEQUENCE_CONFIG

# Process arguments to filter models
arguments = {a.split("=")[0][2:]: a.split("=")[1] for a in sys.argv[1:]
 if a.startswith("--") and "=" in a}
categories = [c for c in arguments.get(
 "categories", "structure,sequence"
).split(",") if c in ["structure", "sequence"]]
families = arguments.get("families")
configs = {"structure": STRUCTURE_CONFIG, "sequence": SEQUENCE_CONFIG}
models = arguments.get("models")

# Go through each category
for category in categories:
    print(f"Building {category} models")
    # What families are there?
    csvs = [csv[:-4] for csv in os.listdir(f"data/csv/{category}")
     if csv.endswith("csv") and (not families or csv[:-4] in families)]
    
    # Go through each family
    for csv in csvs:
        print(f"    Processing {csv}...")
        # Load the dataset and save it as a numpy array
        data = pd.read_csv(f"data/csv/{category}/{csv}.csv").values
        np.save(f"predict/npy/{category}/{csv}", data)

        # Split into train and test data, and save
        data_train, data_test = train_test_split(data)
        data_train, data_test = data_train[:, 1:], data_test[:, 1:]
        np.save(f"predict/npy/{category}/{csv}_train", data_train)
        np.save(f"predict/npy/{category}/{csv}_test", data_test)

        # Get the models to be used
        for Model, hyperparameters in configs[category].items():
            if not models or Model.__name__ in models:
                print("       ", Model.__name__)

                # Create the model with hyperparameters
                model = Pipeline([
                 ("scaler", StandardScaler()),
                 ("model", Model(**hyperparameters))
                ])

                # Train model
                X_train, y_train = data_train[:, :-1], data_train[:, -1].astype("int")
                model.fit(X_train, y_train)

                # Validate model
                scores = cross_validate(
                 model, X_train, y_train, cv=5, scoring=["recall", "precision"]
                )
                model.validation_recall_ = scores["test_recall"].mean()
                model.validation_precision_ = scores["test_precision"].mean()

                # Test model
                X_test, y_test = data_test[:, :-1], data_test[:, -1].astype("int")
                y_pred = model.predict(X_test)
                model.test_recall_ = recall_score(y_test, y_pred)
                model.test_precision_ = precision_score(y_test, y_pred)

                # Save model
                joblib.dump(
                 model,
                 f"predict/models/{category}/{csv}-{Model.__name__}.joblib"
                )