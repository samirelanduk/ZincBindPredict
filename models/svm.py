from sklearn import svm
MODEL = svm.SVC

MODEL_NAME = "SVM"

HYPERPARAMETERS = {
 "kernel": "rbf", "gamma": "scale", "class_weight": {1: 5}, "probability": True
}
