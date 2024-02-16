
import os
import click
import shutil
import uuid

import intact.intact as it
import util.log as log
import util.subtypes as st

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

@cli.command('intact')
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
    '--include-packaging-signal/--exclude-packaging-signal', default=True)
@click.option(
    '--include-rre/--exclude-rre', default=True)
@click.option(
    '--check-major-splice-donor-site/--ignore-major-splice-donor-site', default=True)
@click.option(
    '--run-hypermut/--no-hypermut', default=False
)
@click.option(
    '--check-long-deletion/--ignore-long-deletion', default=False
)
@click.option(
    '--check-nonhiv/--ignore-nonhiv', default=False
)
@click.option(
    '--check-scramble/--ignore-scramble', default=False
)
@click.option(
    '--check-internal-inversion/--ignore-internal-inversion', default=False
)
@click.option(
    '--check-unknown-nucleotides/--ignore-unknown-nucleotides', default=True
)
@click.option(
    '--include-small-orfs/--exclude-small-orfs', default=False)
@click.option(
    '--output-csv/--output-json', default=False)
@click.option(
    '--working-folder',
    default=os.getcwd()
)

def intact(input_file, subtype, include_packaging_signal,
           include_rre, check_major_splice_donor_site, run_hypermut,
           check_long_deletion, check_nonhiv, check_scramble, check_internal_inversion,
           check_unknown_nucleotides, include_small_orfs, output_csv, working_folder):
    """
    Check consensus sequences for intactness.
    """

    log.info('Intactness called.')
    folder = get_working_folder(working_folder)

    it.intact(
        folder, input_file, subtype, include_packaging_signal, include_rre,
        check_major_splice_donor_site, run_hypermut,
        check_long_deletion, check_nonhiv, check_scramble, check_internal_inversion,
        check_unknown_nucleotides, include_small_orfs, output_csv
    )

if __name__ == "__main__": cli()
