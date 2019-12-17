import time
import os
from datetime import datetime, timedelta
from subprocess import Popen
import graphene
from graphql import GraphQLError
from django.conf import settings
from .utilities import save_pdb_code

class JobType(graphene.ObjectType):

    id = graphene.String()
    submitted = graphene.String()
    expires = graphene.String()
    status = graphene.String()

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