import os

def get_categories_from_arguments_and_filesystem(arguments):
    """Gets a dict of dicts for the available categories and datasets in the
    data folder, filtered by whatever the user specifies in arguments."""
    
    categories = {
        category: {
            dataset[:-4]: {} for dataset in os.listdir(
                os.path.join("data", "csv", category)
            ) if dataset.endswith("csv")
        } for category in os.listdir(
            os.path.join("data", "csv")
        ) if "." not in category
    }
    for arg in arguments:
        if arg.startswith("--categories="):
            allowed_categories = arg[13:].split(",")
            categories = {
                k: v for k, v in categories.items() if k in allowed_categories
            }
        if arg.startswith("--datasets="):
            for category in categories:
                allowed_datasets = arg[11:].split(",")
                categories[category] = {k: v for k, v in
                categories[category].items() if k in allowed_datasets}
    return categories