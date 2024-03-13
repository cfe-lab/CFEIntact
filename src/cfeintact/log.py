import logging

FORMAT = '%(asctime)s - %(name)s - %(levelname)s - %(message)s'

logger = logging.getLogger()
logging.basicConfig(level=logging.INFO, format=FORMAT)


def info(info):
    logger.info(info)
