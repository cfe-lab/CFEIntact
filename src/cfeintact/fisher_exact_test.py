from cfeintact.log import logger
from cfeintact.user_error import UserError


def default_fisher_exact(*args, **kwargs):
    raise UserError("Cannot import scipy for fisher_exact test. "
                    "Please check if you have enough RAM for scipy (python -c 'import scipy').")


try:
    import scipy.stats as stats
    fisher_exact = stats.fisher_exact
except ImportError as e:
    logger.warning("Cannot import scipy: %s.", e)
    fisher_exact = default_fisher_exact
except KeyboardInterrupt:
    logger.warning("Cannot import scipy for fisher_exact test. "
                   "Please check if you have enough RAM for scipy (python -c 'import scipy').")
    fisher_exact = default_fisher_exact
