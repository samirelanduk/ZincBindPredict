from sklearn.neighbors import KNeighborsClassifier
from sklearn.ensemble import RandomForestClassifier
from sklearn.svm import SVC

CONFIG = {
    KNeighborsClassifier: {

    },
    RandomForestClassifier: {
        "n_estimators": 10
    },
    SVC: {
        "gamma": "scale",
        "probability": True
    }
}