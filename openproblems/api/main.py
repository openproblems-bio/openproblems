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
    elif args.test_hash:
        hash.git_hash(__file__)
        return

    if args.parallel:
        data.no_cleanup()

    if args.subcommand is None:
        argparser.print_help()
    elif args.subcommand in SUBCOMMANDS:
        return SUBCOMMANDS[args.subcommand].main(args)
    else:
        raise NotImplementedError


def main(args=None, do_print=True):
    """Run the command-line interface.

    Since printing the output to stdout is important here,
    we redirect all other stdout to stderr.
    """
    with utils.RedirectStdout():
        output = _main(args)
    if do_print:
        utils.print_output(output)
        return 0
    else:
        return output
