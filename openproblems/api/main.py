from .. import data
from . import parser
from . import utils

from . import tasks, list, image, load, run, evaluate, hash  # noqa

SUBCOMMANDS = {
    utils.module_to_str(module): module
    for module in [tasks, list, image, load, run, evaluate, hash]
}


def main(args=None, print=True):
    """Run the command-line interface."""
    argparser = parser.create_parser()
    args = argparser.parse_args(args=args)

    if args.parallel:
        data.no_cleanup()

    if args.subcommand is None:
        argparser.print_help()
    elif args.subcommand in SUBCOMMANDS:
        output = SUBCOMMANDS[args.subcommand].main(args)
        if print:
            utils.print_output(output)
            return 0
        else:
            return output
    else:
        raise NotImplementedError
