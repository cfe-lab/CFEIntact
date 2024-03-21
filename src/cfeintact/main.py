
import os
import click
from importlib.metadata import version

import cfeintact.intact as it
import cfeintact.log as log
import cfeintact.subtypes as st


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
    if __package__ is None:
        print("CFEIntact is not installed.")
        exit(1)
    else:
        print(version(__package__))
        exit(0)


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
    '--check-packaging-signal/--ignore-packaging-signal', default=True)
@click.option(
    '--check-rre/--ignore-rre', default=True)
@click.option(
    '--check-major-splice-donor-site/--ignore-major-splice-donor-site', default=True)
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
    '--check-small-orfs/--ignore-small-orfs', default=True)
@click.option(
    '--output-csv/--output-json', default=True)
@click.option(
    '--working-folder',
    default=os.getcwd()
)
def intact(input_file: str, subtype: str, check_packaging_signal: str,
           check_rre: str, check_major_splice_donor_site: str, check_hypermut: str,
           check_long_deletion: str, check_nonhiv: str, check_scramble: str, check_internal_inversion: str,
           check_unknown_nucleotides: str, check_small_orfs: str, output_csv: str, working_folder: str) -> None:
    """
    Check consensus sequences for intactness.
    """

    log.info('Intactness called.')
    folder = get_working_folder(working_folder)

    it.intact(
        folder, input_file, subtype, check_packaging_signal, check_rre,
        check_major_splice_donor_site, check_hypermut,
        check_long_deletion, check_nonhiv, check_scramble, check_internal_inversion,
        check_unknown_nucleotides, check_small_orfs, output_csv
    )

    exit(0)


if __name__ == "__main__": cli()
