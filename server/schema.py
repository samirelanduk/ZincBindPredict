import time
import os
import json
from datetime import datetime, timedelta
from subprocess import Popen
import graphene
import joblib
from graphene.relay import Connection, ConnectionField
from graphql import GraphQLError
from django.conf import settings
from .utilities import get_job_location, initialise_job

results = None

'''class ResidueType(graphene.ObjectType):

    id = graphene.String()
    name = graphene.String()
    


class SiteType(graphene.ObjectType):

    probability = graphene.Float()
    family = graphene.String()
    residues = graphene.List(ResidueType)
    vector = graphene.List(graphene.Float)
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

    method = graphene.String()
    family = graphene.String()
    category = graphene.String()
    recall = graphene.Float()
    precision = graphene.Float()
    hyperparameters = graphene.String()



class JobModelType(graphene.ObjectType):
    method = graphene.String()
    family = graphene.String()
    category = graphene.String()
    sites = graphene.ConnectionField(SiteConnection)
    rejected = graphene.ConnectionField(SiteConnection)
    recall = graphene.Float()
    precision = graphene.Float()
    test_recall = graphene.Float()
    test_precision = graphene.Float()
    test_roc_auc = graphene.Float()

    def resolve_sites(self, info, **kwargs):
        return [SiteType(**s) for s in self.sites]
    

    def resolve_rejected(self, info, **kwargs):
        return [SiteType(**s) for s in self.rejected]



class ModelConnection(Connection):

    class Meta:
        node = ModelType
    
    count = graphene.Int()

    def resolve_count(self, info, **kwargs):
        return len(self.edges)



class JobModelConnection(Connection):

    class Meta:
        node = JobModelType
    
    count = graphene.Int()

    def resolve_count(self, info, **kwargs):
        return len(self.edges)



class JobType(graphene.ObjectType):

    id = graphene.String()
    submitted = graphene.String()
    expires = graphene.String()
    status = graphene.String()
    models = graphene.ConnectionField(JobModelConnection)
    sites = graphene.ConnectionField(SiteConnection)
    rejected = graphene.ConnectionField(SiteConnection)

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
            for key, value in results.items():
                if key != "time":
                    value["family"] = value["name"].split("-")[0]
                    value["method"] = "".join(value["name"].split("-")[1:])
                    del value["name"]
            return [JobModelType(**r) for k, r in results.items() if k != "time"]
        except FileNotFoundError:
            return []
    

    def resolve_sites(self, info, **kwargs):
        global results
        try:
            with open(f"server/jobs/{self.id}/results.json") as f:
                results = json.load(f)
            return [SiteType(**site) for sites in [m["sites"] for k, m in results.items() if k != "time"] for site in sites]
        except FileNotFoundError:
            return []
    

    def resolve_rejected(self, info, **kwargs):
        global results
        try:
            with open(f"server/jobs/{self.id}/results.json") as f:
                results = json.load(f)
            return [SiteType(**site) for sites in [m["rejected"] for m in results.values()] for site in sites]
        except FileNotFoundError:
            return []



class Query(graphene.ObjectType):
    job = graphene.Field(JobType, id=graphene.String(required=True))
    model = graphene.Field(
        ModelType, category=graphene.String(required=True),
        family=graphene.String(required=True), method=graphene.String(required=True)
    )
    models = graphene.ConnectionField(
        ModelConnection, category=graphene.String(),
        family=graphene.String(), method=graphene.String()
    )

    def resolve_job(self, info, **kwargs):
        if os.path.exists(f"server/jobs/{kwargs['id']}"):
            return JobType(id=kwargs["id"])
    

    def resolve_model(self, info, **kwargs):
        with open("predict/models/results.json") as f:
            model_results = json.load(f)
        model = None
        for category in model_results:
            if kwargs["category"] == category: 
                for family in model_results[category]:
                    if kwargs["family"] == family: 
                        for model_name in model_results[category][family]:
                            method = model_name.replace(" ", "")
                            if kwargs["method"] == method: 
                                model_results[category][family][model_name]["category"] = category
                                model_results[category][family][model_name]["family"] = family
                                model_results[category][family][model_name]["method"] = method
                                model_results[category][family][model_name]["hyperparameters"] = ",".join(
                                    [f"{k}:{str(v)}" for k, v in model_results[category][family][model_name]["hyperparameters"].items()]
                                )
                                model = model_results[category][family][model_name]
        return ModelType(
            **{k: model[k] for k in ["method", "family", "category", "recall", "precision", "hyperparameters"]}
        ) if model else None
    

    def resolve_models(self, info, **kwargs):
        with open("predict/models/results.json") as f:
            model_results = json.load(f)
        models = []
        for category in model_results:
            if "category" not in kwargs or kwargs["category"] == category: 
                for family in model_results[category]:
                    if "family" not in kwargs or kwargs["family"] == family: 
                        for model_name in model_results[category][family]:
                            if model_name == "ensemble": continue
                            method = model_name.replace(" ", "")
                            if "method" not in kwargs or kwargs["method"] == method:
                                model_results[category][family][model_name]["category"] = category
                                model_results[category][family][model_name]["family"] = family
                                model_results[category][family][model_name]["method"] = method
                                model_results[category][family][model_name]["hyperparameters"] = ",".join(
                                    [f"{k}:{str(v)}" for k, v in model_results[category][family][model_name]["hyperparameters"].items()]
                                )
                                models.append(model_results[category][family][model_name])
        return [ModelType(
            **{k: model[k] for k in ["method", "family", "category", "recall", "precision", "hyperparameters"]}
        ) for model in models]

 
class SubmitStructure(graphene.Mutation):

    class Arguments:
        code = graphene.String()
        assembly = graphene.Int(required=False)

    job = graphene.Field(JobType)

    def mutate(self, info, **kwargs):
        job_id = int(time.time() * 1000)
        os.mkdir(f"server/jobs/{job_id}")
        if not save_pdb_code(kwargs["code"], job_id):
            raise GraphQLError("That does not seem to be a valid PDB code")
        Popen(["server/structure_job.py", str(job_id), str(kwargs.get("assembly", 0))])
        return SubmitStructure(job=JobType(id=str(job_id)))'''



class ResidueType(graphene.ObjectType):

    id = graphene.String()
    name = graphene.String()



class SiteType(graphene.ObjectType):

    probability = graphene.Float()
    family = graphene.String()
    residues = graphene.List(ResidueType)

    def resolve_residues(self, info, **kwargs):
        return [ResidueType(**res) for res in self.residues]



class JobType(graphene.ObjectType):
    
    id = graphene.String()
    time = graphene.String()
    status = graphene.String()
    type = graphene.String()
    protein = graphene.String()
    sites = graphene.List(SiteType)
    rejected = graphene.List(SiteType)

    def resolve_time(self, info, **kwargs):
        return datetime.utcfromtimestamp(
            int(self.id) / 1000
        ).strftime("%Y-%m-%d %H:%M:%S UTC")
    

    def resolve_sites(self, info, **kwargs):
        return [SiteType(**site) for site in self.sites]
    

    def resolve_rejected(self, info, **kwargs):
        return [SiteType(**site) for site in self.rejected]
    


class Query(graphene.ObjectType):
    
    job = graphene.Field(JobType, id=graphene.String(required=True))

    def resolve_job(self, info, **kwargs):
        with open(get_job_location(kwargs["id"])) as f:
            return JobType(**json.load(f))



class SearchStructure(graphene.Mutation):

    class Arguments:
        structure = graphene.String(required=True)
        families = graphene.List(graphene.String, description="Site families to limit to")
        probability = graphene.Float(description="Probability threshold")
        find_half = graphene.Boolean(description="Look for half sites instead")
        use_families_models = graphene.Boolean(description="Use family based models")
        use_location_models = graphene.Boolean(description="Use location based models")

    job_id = graphene.String()

    def mutate(self, info, **kwargs):
        kwargs["job_id"] = str(int(time.time() * 1000))
        job = initialise_job(kwargs["job_id"], "structure", kwargs["structure"])
        with open(get_job_location(kwargs["job_id"]), "w") as f:
            json.dump(job, f)
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
        with open(get_job_location(kwargs["job_id"]), "w") as f:
            json.dump(job, f)
        Popen(["server/sequence_job.py", json.dumps(kwargs)])
        return SearchSequence(job_id=kwargs["job_id"])



class Mutations(graphene.ObjectType):
    search_structure = SearchStructure.Field()
    search_sequence = SearchSequence.Field()


schema = graphene.Schema(query=Query, mutation=Mutations)