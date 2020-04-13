import json
from datetime import datetime
from subprocess import Popen
import graphene
import joblib
from graphene_file_upload.scalars import Upload
from .utilities import *

class ResidueType(graphene.ObjectType):

    identifier = graphene.String()
    name = graphene.String()



class SequenceSiteType(graphene.ObjectType):

    probability = graphene.Float()
    family = graphene.String()
    residues = graphene.String()



class StructureSiteType(graphene.ObjectType):

    probability = graphene.Float()
    family = graphene.String()
    residues = graphene.List(ResidueType)

    def resolve_residues(self, info, **kwargs):
        return [ResidueType(**res) for res in self.residues]



class AbstractJobType:

    id = graphene.String()
    status = graphene.String()
    protein = graphene.String()
    time = graphene.String()

    def resolve_time(self, info, **kwargs):
        return datetime.utcfromtimestamp(
            int(self.id) / 1000
        ).strftime("%Y-%m-%d %H:%M:%S UTC")



class SequenceJobType(AbstractJobType, graphene.ObjectType):
    
    
    sites = graphene.List(SequenceSiteType)
    rejected = graphene.List(SequenceSiteType)

    def resolve_sites(self, info, **kwargs):
        return [SequenceSiteType(**site) for site in self.sites]
    

    def resolve_rejected(self, info, **kwargs):
        return [SequenceSiteType(**site) for site in self.rejected]



class StructureJobType(AbstractJobType, graphene.ObjectType):

    sites = graphene.List(StructureSiteType)
    rejected = graphene.List(StructureSiteType)

    def resolve_sites(self, info, **kwargs):
        return [StructureSiteType(**site) for site in self.sites]
    

    def resolve_rejected(self, info, **kwargs):
        return [StructureSiteType(**site) for site in self.rejected]
    


class Query(graphene.ObjectType):
    
    sequence_job = graphene.Field(SequenceJobType, id=graphene.String(required=True))
    structure_job = graphene.Field(StructureJobType, id=graphene.String(required=True))

    def resolve_sequence_job(self, info, **kwargs):
        with open(get_job_location(kwargs["id"])) as f:
            return SequenceJobType(**json.load(f))
    

    def resolve_structure_job(self, info, **kwargs):
        with open(get_job_location(kwargs["id"])) as f:
            return StructureJobType(**json.load(f))



class SearchSequence(graphene.Mutation):

    class Arguments:
        sequence = graphene.String(required=True)
        families = graphene.List(graphene.String, description="Site families to limit to")
        probability = graphene.Float(description="Probability threshold")
    
    job_id = graphene.String()

    def mutate(self, info, **kwargs):
        job = initialise_job(kwargs["sequence"])
        save_job(job)
        kwargs["job_id"] = job["id"]
        Popen(["server/sequence_job.py", json.dumps(kwargs)])
        return SearchSequence(job_id=job["id"])



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
        job = initialise_job("")
        job["protein"] = kwargs["structure"] = save_structure_file(
            kwargs["structure"], job["id"])
        save_job(job)
        kwargs["job_id"] = job["id"]
        Popen(["server/structure_job.py", json.dumps(kwargs)])
        return SearchStructure(job_id=job["id"])
        


class Mutations(graphene.ObjectType):
    search_structure = SearchStructure.Field()
    search_sequence = SearchSequence.Field()


schema = graphene.Schema(query=Query, mutation=Mutations)