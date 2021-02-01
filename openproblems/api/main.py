from .. import data
from . import parser
from . import utils
from . import tasks, list, image, load, run, evaluate  # noqa

SUBCOMMANDS = {
    utils.module_to_str(module): module
    for module in [tasks, list, image, load, run, evaluate]
}


def main():
    """Run the command-line interface."""
    argparser = parser.create_parser()
    args = argparser.parse_args()

    if args.parallel:
        data.no_cleanup()

    if args.subcommand is None:
        argparser.print_help()
    elif args.subcommand in SUBCOMMANDS:
        return SUBCOMMANDS[args.subcommand].main(args)
    else:
        raise NotImplementedError
