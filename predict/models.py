from sklearn.ensemble import RandomForestClassifier
from sklearn.neighbors import KNeighborsClassifier
from sklearn.svm import SVC

MODELS = {
    RandomForestClassifier: {"name": "Random Forest", "hyperparameters": { 
        "criterion": ["gini"], "n_estimators": [100], "max_depth": [None]
    }},
    KNeighborsClassifier: {"name": "K Nearest Neighbor", "hyperparameters": {
        "n_neighbors": [1, 3, 5]
    }},
    #SVC: {"name": "SVM", "hyperparameters": {
    #    "C": [0.5, 1, 2], "kernel": ["rbf", "linear", "sigmoid"],
    #    "gamma": ["scale", "auto"]
    #}}
    SVC: {"name": "SVM", "hyperparameters": {
        "C": [0.5], "kernel": ["rbf"],
        "gamma": ["auto"]
    }}
}