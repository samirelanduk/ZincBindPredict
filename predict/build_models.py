#! /usr/bin/env python3

import warnings
warnings.warn = lambda *args, **kwargs: None
import sys
sys.path.append("../zincbindpredict")
import pandas as pd
import joblib
from tqdm import tqdm
from sklearn.model_selection import train_test_split
from predict.utilities import *

# Get options
categories, families, similarities, sizes = parse_data_args(sys.argv)

# Create progress bar
model_count = len(families) * len(sizes) * len(similarities)
if categories == ["structure"]: model_count /= len(similarities)
if len(categories) == 2: model_count += len(families) * len(sizes)
with tqdm(total=model_count) as pbar:

    # Go through each category
    for category in categories:

        # Go through each family
        for family in families:

            # Go through each similarity
            for sim in similarities if category == "sequence" else [None]:
                
                # Load data
                df = pd.read_csv(os.path.join(
                    "data", "csv", category,
                    f"{family}_{int(sim * 100)}.csv" if sim else f"{family}.csv"
                ))

                # Go through each dataset size
                for size in sizes:

                    # Get subset of data
                    data = get_data_subset(df, size)

                    # Get four datasets
                    X, y = data[:, :-1], data[:, -1]
                    X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2)

                    # Train model
                    model = train_model(X_train, y_train)
                    
                    # Evaluate model
                    model.results = evaluate_model(model, X_test, y_test)
                    
                    # Save model
                    model_name = f"{family}_{int(size * 100)}"
                    if sim is not None: model_name += f"_{int(sim * 100)}"
                    joblib.dump(model, os.path.join(
                        "predict", "models", category, f"{model_name}.joblib"
                    ))
                    pbar.update()