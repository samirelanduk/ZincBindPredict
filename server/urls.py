from graphene_django.views import GraphQLView
from graphene_file_upload.django import FileUploadGraphQLView
from django.urls import path

urlpatterns = [
 path("graphql/", FileUploadGraphQLView.as_view(graphiql=True)),
 path("", FileUploadGraphQLView.as_view(graphiql=True)),
]
