from django.urls import path
from .views import *

urlpatterns = [
 path("structure/", predict_structure),
 path("structure/<int:id>/", job),
 path("sequence/", predict_sequence),
 path("", root)
]
