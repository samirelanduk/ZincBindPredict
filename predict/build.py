import os
import pandas as pd
import numpy as np
from sklearn.model_selection import train_test_split
from sklearn.metrics import recall_score, precision_score
import joblib
from config import CONFIG

# What families are there?
csvs = [csv[:-4] for csv in os.listdir("data/csv") if csv.endswith("csv")]

# Go through each family
for csv in csvs:
    print(f"Processing {csv}...")
    # Load the dataset and save it as a numpy array
    data = pd.read_csv(f"data/csv/{csv}.csv").values
    np.save(f"predict/npy/{csv}", data)

    # Split into train and test data, and save
    data_train, data_test = train_test_split(data)
    data_train, data_test = data_train[:, 1:], data_test[:, 1:]
    np.save(f"predict/npy/{csv}_train", data_train)
    np.save(f"predict/npy/{csv}_test", data_test)

    # Get the models to be used
    for Model, hyperparameters in CONFIG.items():
        print("   ", Model.__name__)

        # Create the model with hyperparameters
        model = Model(**hyperparameters)

        # Train model
        X_train, y_train = data_train[:, :-1], data_train[:, -1].astype("int")
        model.fit(X_train, y_train)

        # Test model
        X_test, y_test = data_test[:, :-1], data_test[:, -1].astype("int")
        y_pred = model.predict(X_test)
        model.test_recall_ = recall_score(y_test, y_pred)
        model.test_precision_ = precision_score(y_test, y_pred)

        # Save model
        joblib.dump(model, f"predict/models/{csv}-{Model.__name__}.joblib")