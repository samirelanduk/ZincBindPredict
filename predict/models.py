from sklearn.ensemble import RandomForestClassifier
from sklearn.neighbors import KNeighborsClassifier

MODELS = {
    RandomForestClassifier: {"name": "Random Forest", "hyperparameters": { 
        "criterion": ["gini"], "n_estimators": [100], "max_depth": [None]
    }},
    KNeighborsClassifier: {"name": "K Nearest Neighbor", "hyperparameters": {
        "n_neighbors": [1, 3, 5]
    }}
}