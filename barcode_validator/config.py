import os
import sys
import yaml
import logging


def exception_handler(exc_type, exc_value, exc_traceback):
    """
    Custom exception handler to log uncaught exceptions
    """
    if issubclass(exc_type, KeyboardInterrupt):
        # Don't log keyboard interrupt (Ctrl+C)
        sys.__excepthook__(exc_type, exc_value, exc_traceback)
        return

    # Log the exception
    logging.error("Uncaught exception", exc_info=(exc_type, exc_value, exc_traceback))


class Config:
    _instance = None
    _initialized = False

    def __new__(cls):
        if cls._instance is None:
            cls._instance = super(Config, cls).__new__(cls)
        return cls._instance

    def __init__(self):
        if not self._initialized:
            self.config_data = None
            self.config_path = None
            self._initialized = True

    def load_config(self, config_path):
        """
        Load the configuration from a YAML file.

        Args:
        config_path (str): Path to the configuration YAML file.

        Raises:
        FileNotFoundError: If the config file doesn't exist.
        yaml.YAMLError: If there's an error parsing the YAML file.
        """
        self.config_path = config_path
        if not os.path.exists(config_path):
            raise FileNotFoundError(f"Config file not found: {config_path}")

        with open(config_path, 'r') as config_file:
            try:
                self.config_data = yaml.safe_load(config_file)
            except yaml.YAMLError as e:
                raise yaml.YAMLError(f"Error parsing config file: {e}")

        # Convert relative paths to absolute paths
        self._process_relative_paths(os.path.dirname(os.path.abspath(config_path)))

    def _process_relative_paths(self, config_dir):
        """
        Convert relative paths in the config to absolute paths.

        Args:
        config_dir (str): Directory of the config file.
        """
        for key, value in self.config_data.items():
            if isinstance(value, str) and not os.path.isabs(value) and key.endswith('_file'):
                self.config_data[key] = os.path.join(config_dir, value)

    def setup_logging(self, cmd_log_level=None):
        """
        Set up logging based on the configuration file and command line argument.

        Args:
        cmd_log_level (str, optional): Log level from command line, overrides config file if provided.
        """
        if self.config_data is None:
            raise ValueError("Configuration has not been loaded. Call load_config first.")

        log_level = self.config_data.get('log_level', 'INFO')
        if cmd_log_level:
            log_level = cmd_log_level

        numeric_level = getattr(logging, log_level.upper(), None)
        if not isinstance(numeric_level, int):
            raise ValueError(f'Invalid log level: {log_level}')

        logging.basicConfig(level=numeric_level,
                            format='%(asctime)s - %(levelname)s - %(message)s',
                            filename=self.config_data.get('log_file'),
                            filemode='a')

        sys.excepthook = exception_handler

    def get(self, key, default=None):
        """
        Get a configuration value.

        Args:
        key (str): The configuration key to retrieve.
        default: The default value to return if the key is not found.

        Returns:
        The value associated with the key, or the default value if not found.
        """
        return self.config_data.get(key, default) if self.config_data else default

    def __getitem__(self, key):
        """
        Allow dictionary-style access to configuration values.

        Args:
        key (str): The configuration key to retrieve.

        Returns:
        The value associated with the key.

        Raises:
        KeyError: If the key is not found in the configuration.
        """
        if self.config_data is None:
            raise ValueError("Configuration has not been loaded. Call load_config first.")
        return self.config_data[key]

    def __contains__(self, key):
        """
        Check if a key exists in the configuration.

        Args:
        key (str): The configuration key to check.

        Returns:
        bool: True if the key exists, False otherwise.
        """
        return self.config_data is not None and key in self.config_data
