from django.urls import path
from .views import *

urlpatterns = [
 path("predict/", predict),
 path("predict/structure/", predict_structure),
 path("predict/sequence/", predict_sequence),
 path("", root)
]
