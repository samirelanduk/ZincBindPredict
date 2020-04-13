import time
import os
import json
from datetime import datetime, timedelta
from subprocess import Popen
import graphene
import joblib
from graphene.relay import Connection, ConnectionField
from graphene_file_upload.scalars import Upload
from graphql import GraphQLError
from django.conf import settings
from .utilities import *

class ResidueType(graphene.ObjectType):

    identifier = graphene.String()
    name = graphene.String()



class SequenceSiteType(graphene.ObjectType):

    probability = graphene.Float()
    family = graphene.String()
    residues = graphene.String()



class SequenceJobType(graphene.ObjectType):
    
    id = graphene.String()
    time = graphene.String()
    status = graphene.String()
    type = graphene.String()
    protein = graphene.String()
    sites = graphene.List(SequenceSiteType)
    rejected = graphene.List(SequenceSiteType)

    def resolve_time(self, info, **kwargs):
        return datetime.utcfromtimestamp(
            int(self.id) / 1000
        ).strftime("%Y-%m-%d %H:%M:%S UTC")
    

    def resolve_sites(self, info, **kwargs):
        return [SequenceSiteType(**site) for site in self.sites]
    

    def resolve_rejected(self, info, **kwargs):
        return [SequenceSiteType(**site) for site in self.rejected]
    


class Query(graphene.ObjectType):
    
    sequence_job = graphene.Field(SequenceJobType, id=graphene.String(required=True))
    #structure_job = graphene.Field(StructureJobType, id=graphene.String(required=True))

    def resolve_sequence_job(self, info, **kwargs):
        with open(get_job_location(kwargs["id"])) as f:
            return SequenceJobType(**json.load(f))
    

    def resolve_structure_job(self, info, **kwargs):
        with open(get_job_location(kwargs["id"])) as f:
            return SequenceJobType(**json.load(f))



class SearchStructure(graphene.Mutation):

    class Arguments:
        structure = Upload(required=True)
        families = graphene.List(graphene.String, description="Site families to limit to")
        probability = graphene.Float(description="Probability threshold")
        find_half = graphene.Boolean(description="Look for half sites instead")
        use_families_models = graphene.Boolean(description="Use family based models")
        use_location_models = graphene.Boolean(description="Use location based models")

    job_id = graphene.String()

    def mutate(self, info, **kwargs):
        kwargs["job_id"] = str(int(time.time() * 1000))
        file_extension = str(kwargs["structure"]).split(".")[-1]
        file_location = get_job_location(kwargs["job_id"])[:-4] + file_extension
        with open(file_location, "wb") as f:
            f.write(kwargs["structure"].read())
        kwargs["structure"] = kwargs["job_id"] + "." + file_extension
        job = initialise_job(kwargs["job_id"], "structure", kwargs["structure"])
        save_job(job)
        Popen(["server/structure_job.py", json.dumps(kwargs)])
        return SearchStructure(job_id=kwargs["job_id"])



class SearchSequence(graphene.Mutation):

    class Arguments:
        sequence = graphene.String(required=True)
        families = graphene.List(graphene.String, description="Site families to limit to")
        probability = graphene.Float(description="Probability threshold")
    
    job_id = graphene.String()

    def mutate(self, info, **kwargs):
        kwargs["job_id"] = str(int(time.time() * 1000))
        job = initialise_job(kwargs["job_id"], "sequence", kwargs["sequence"])
        save_job(job)
        Popen(["server/sequence_job.py", json.dumps(kwargs)])
        return SearchSequence(job_id=kwargs["job_id"])



class Mutations(graphene.ObjectType):
    search_structure = SearchStructure.Field()
    search_sequence = SearchSequence.Field()


schema = graphene.Schema(query=Query, mutation=Mutations)