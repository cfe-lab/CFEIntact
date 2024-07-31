from cfeintact.log import logger
from cfeintact.user_error import UserError


def default_fisher_exact(*args, **kwargs):
    raise UserError("Cannot import scipy for fisher_exact test.")


try:
    import scipy.stats as stats
    fisher_exact = stats.fisher_exact
except Exception as e:
    logger.warning("Cannot import scipy: %s.", e)
    fisher_exact = default_fisher_exact
