import logging

FORMAT = '%(levelname)s: [CFEIntact] %(message)s'

log = logging.getLogger()
logging.basicConfig(level=logging.INFO, format=FORMAT, datefmt='%Y-%m-%d %H:%M:%S')
