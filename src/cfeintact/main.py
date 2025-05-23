
import os
import click
from importlib.metadata import version, PackageNotFoundError

import cfeintact.intact as it
import cfeintact.subtypes as st
from cfeintact.user_error import UserError
from cfeintact.log import logger


def get_working_folder(path):
    """
    Returns a default working folder.
    """

    if not os.path.exists(path):
        os.makedirs(path)

    return path


@click.group()
def cli():
    pass


@cli.command('version')
def get_version() -> None:
    try:
        print(version(__package__))
        exit(0)
    except PackageNotFoundError:
        print("CFEIntact is not installed.")
        exit(1)


@cli.command('check')
@click.argument(
    'input_file',
    type=click.Path(exists=True, readable=True, resolve_path=True)
)
@click.option(
    '--subtype',
    required=True,
    type=click.Choice(st.subtypes())
)
@click.option(
    '--check-packaging-signal/--ignore-packaging-signal', default=True
)
@click.option(
    '--check-rre/--ignore-rre', default=True
)
@click.option(
    '--check-major-splice-donor-site/--ignore-major-splice-donor-site', default=True
)
@click.option(
    '--check-hypermut/--ignore-hypermut', default=True
)
@click.option(
    '--check-long-deletion/--ignore-long-deletion', default=True
)
@click.option(
    '--check-nonhiv/--ignore-nonhiv', default=True
)
@click.option(
    '--check-scramble/--ignore-scramble', default=True
)
@click.option(
    '--check-internal-inversion/--ignore-internal-inversion', default=True
)
@click.option(
    '--check-unknown-nucleotides/--ignore-unknown-nucleotides', default=True
)
@click.option(
    '--check-small-orfs/--ignore-small-orfs', default=True
)
@click.option(
    '--check-distance/--ignore-distance', default=True
)
@click.option(
    '--output-csv/--output-json', default=True
)
@click.option(
    '--output',
    default=os.getcwd()
)
def intact(input_file: str, subtype: str, check_packaging_signal: bool,
           check_rre: bool, check_major_splice_donor_site: bool, check_hypermut: bool,
           check_long_deletion: bool, check_nonhiv: bool, check_scramble: bool, check_internal_inversion: bool,
           check_unknown_nucleotides: bool, check_small_orfs: bool, check_distance: bool, output_csv: bool,
           output: str) -> None:
    """
    Check consensus sequences for intactness.
    """

    logger.info('Intactness called.')
    folder = get_working_folder(output)

    try:
        it.check(
            folder, input_file, subtype, check_packaging_signal, check_rre,
            check_major_splice_donor_site, check_hypermut,
            check_long_deletion, check_nonhiv, check_scramble, check_internal_inversion,
            check_unknown_nucleotides, check_small_orfs, check_distance, output_csv,
        )
    except UserError as e:
        logger.fatal(e.fmt, *e.fmt_args)
        exit(e.code)

    exit(0)


if __name__ == "__main__": cli()
