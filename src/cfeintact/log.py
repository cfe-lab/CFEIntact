import logging

FORMAT = '%(asctime)s - %(name)s - %(levelname)s - %(message)s'

log = logging.getLogger()
logging.basicConfig(level=logging.INFO, format=FORMAT)
