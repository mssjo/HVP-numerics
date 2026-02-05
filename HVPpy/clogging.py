# clogging: logging but with more color

import logging
import sys

clogger = logging.getLogger("clogging")

# From https://stackoverflow.com/questions/384076/how-can-i-color-python-logging-output
# ... with some customization
class ColorFormatter(logging.Formatter):

    BLACK, RED, GREEN, YELLOW, BLUE, MAGENTA, CYAN, WHITE = range(8)
    RESET = "\033[0m"
    COLOR = "\033[3%dm"
    BRIGHTCOLOR = "\033[9%dm"
    BOLD = "\033[1m"
    REVERSE = "\033[7m"

    def __init__(self, format):
        super().__init__()


        self.formats = {
            logging.DEBUG: self.COLOR%self.BLUE + format + self.RESET,
            logging.INFO: self.COLOR%self.WHITE + self.BOLD + format + self.RESET,
            logging.WARNING: self.COLOR%self.YELLOW + self.BOLD + format + self.RESET,
            logging.ERROR: self.COLOR%self.RED + self.BOLD + format + self.RESET,
            logging.CRITICAL: self.COLOR%self.RED + self.REVERSE + format + self.RESET,
            None: format
        }

    def format(self, record):
        return logging.Formatter(fmt=self.formats.get(record.levelno, self.formats[None])).format(record)


def init_clogging(level=logging.INFO, **kwargs):
    format="[%(levelname)s] %(message)s"
    handler = logging.StreamHandler()
    if sys.stdout.isatty():
        handler.setFormatter(ColorFormatter(format))
    else:
        handler.setFormatter(logging.Formatter(format))

    logging.basicConfig(level=level, handlers=[handler], **kwargs)
