import logging

# Configure the logger
logging.basicConfig(
    level=logging.INFO,  # Set the desired logging level
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[
        logging.StreamHandler(),  # Output logs to the console
    ]
)

# Create a logger instance
logger = logging.getLogger(__name__)