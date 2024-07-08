import logging

FORMAT = '%(levelname)s: [CFEIntact] %(message)s'

logger = logging.getLogger()
logging.basicConfig(level=logging.INFO, format=FORMAT, datefmt='%Y-%m-%d %H:%M:%S')
