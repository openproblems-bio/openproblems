from . import utils


def get_image(task_name, function_type, function_name):
    """Get the name of the Docker image associated with a function."""
    fun = utils.get_function(task_name, function_type, function_name)
    # print docker image name
    return fun.metadata["image"]


def main(args):
    """Run the ``image`` subcommand."""
    image = get_image(args.task, args.function_type, args.name)
    return image
