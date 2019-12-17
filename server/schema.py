import graphene

class Query(graphene.ObjectType):
    field = graphene.String()
    

schema = graphene.Schema(query=Query)