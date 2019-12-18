import time
import os
import json
from datetime import datetime, timedelta
from subprocess import Popen
import graphene
from graphene.relay import Connection, ConnectionField
from graphql import GraphQLError
from django.conf import settings
from .utilities import save_pdb_code

results = None

class ResidueType(graphene.ObjectType):

    id = graphene.String()
    name = graphene.String()
    


class SiteType(graphene.ObjectType):

    probability = graphene.Float()
    family = graphene.String()
    residues = graphene.List(ResidueType)
    model = graphene.Field(lambda: ModelType)

    def resolve_residues(self, info, **kwargs):
        return [ResidueType(**r) for r in self.residues]
    

    def resolve_model(self, info, **kwargs):
        return ModelType(**results[self.model])



class SiteConnection(Connection):

    class Meta:
        node = SiteType
    
    count = graphene.Int()

    def resolve_count(self, info, **kwargs):
        return len(self.edges)



class ModelType(graphene.ObjectType):

    name = graphene.String()
    sites = graphene.ConnectionField(SiteConnection)
    validation_recall = graphene.Float()
    validation_precision = graphene.Float()
    test_recall = graphene.Float()
    test_precision = graphene.Float()

    def resolve_sites(self, info, **kwargs):
        return [SiteType(**s) for s in self.sites]



class ModelConnection(Connection):

    class Meta:
        node = ModelType
    
    count = graphene.Int()

    def resolve_count(self, info, **kwargs):
        return len(self.edges)



class JobType(graphene.ObjectType):

    id = graphene.String()
    submitted = graphene.String()
    expires = graphene.String()
    status = graphene.String()
    models = graphene.ConnectionField(ModelConnection)
    sites = graphene.ConnectionField(SiteConnection)

    def resolve_submitted(self, info, **kwargs):
        return str(datetime.utcfromtimestamp(int(self.id) // 1000)) + " UTC"
    

    def resolve_expires(self, info, **kwargs):
        return str(datetime.utcfromtimestamp(
            int(self.id) // 1000
        ) + timedelta(days=settings.JOB_EXPIRATION)) + " UTC"
    

    def resolve_status(self, info, **kwargs):
        try:
            with open(f"server/jobs/{self.id}/status.txt") as f:
                return f.read()
        except FileNotFoundError:
            return "starting"
    

    def resolve_models(self, info, **kwargs):
        global results
        try:
            with open(f"server/jobs/{self.id}/results.json") as f:
                results = json.load(f)
            return [ModelType(**r) for r in results.values()]
        except FileNotFoundError:
            return []
    

    def resolve_sites(self, info, **kwargs):
        global results
        try:
            with open(f"server/jobs/{self.id}/results.json") as f:
                results = json.load(f)
            return [SiteType(**site) for sites in [m["sites"] for m in results.values()] for site in sites]
            
        except FileNotFoundError:
            return []



class SubmitStructure(graphene.Mutation):

    class Arguments:
        code = graphene.String()

    job = graphene.Field(JobType)

    def mutate(self, info, **kwargs):
        job_id = int(time.time() * 1000)
        os.mkdir(f"server/jobs/{job_id}")
        if not save_pdb_code(kwargs["code"], job_id):
            raise GraphQLError("That does not seem to be a valid PDB code")
        Popen(["server/structure_job.py", str(job_id)])
        return SubmitStructure(job=JobType(id=str(job_id)))



class Query(graphene.ObjectType):
    job = graphene.Field(JobType, id=graphene.String(required=True))

    def resolve_job(self, info, **kwargs):
        if os.path.exists(f"server/jobs/{kwargs['id']}"):
            return JobType(id=kwargs["id"])



class Mutations(graphene.ObjectType):
    submit_structure = SubmitStructure.Field()


schema = graphene.Schema(query=Query, mutation=Mutations)