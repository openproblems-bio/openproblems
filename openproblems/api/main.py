from .. import __version__
from .. import data
from . import evaluate
from . import hash
from . import image
from . import list
from . import load
from . import parser
from . import run
from . import tasks
from . import utils

SUBCOMMANDS = {
    utils.module_to_str(module): module
    for module in [tasks, list, image, load, run, evaluate, hash]
}


def _main(args=None):
    argparser = parser.create_parser()
    args = argparser.parse_args(args=args)

    if args.version:
        return __version__

    if args.parallel:
        data.no_cleanup()

    if args.subcommand is None:
        argparser.print_help()
    elif args.subcommand in SUBCOMMANDS:
        return SUBCOMMANDS[args.subcommand].main(args)
    else:
        raise NotImplementedError


def main(args=None, print=True):
    """Run the command-line interface."""
    output = _main(args)
    if print and output:
        utils.print_output(output)
        return 0
    else:
        return output
