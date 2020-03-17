import logging

LOGGER = logging.getLogger("bilana2.log")

# Creating Streamhandler -- This handler will write to console
SH = logging.StreamHandler()
SH.setLevel(logging.DEBUG)
SH.setFormatter(logging.Formatter('%(asctime)s %(funcName)s - %(levelname)s - %(message)s'))
LOGGER.addHandler(SH)
