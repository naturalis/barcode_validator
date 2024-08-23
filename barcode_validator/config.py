import os
import sys
import yaml
import logging


def exception_handler(exc_type, exc_value, exc_traceback):
    """
    Custom exception handler to log uncaught exceptions.
    :param exc_type: Exception type
    :param exc_value: Exception value
    :param exc_traceback: Exception traceback
    """
    if issubclass(exc_type, KeyboardInterrupt):
        # Don't log keyboard interrupt (Ctrl+C)
        sys.__excepthook__(exc_type, exc_value, exc_traceback)
        return

    # Log the exception
    logging.error("Uncaught exception", exc_info=(exc_type, exc_value, exc_traceback))


class Config:
    _instance = None

    def __new__(cls):
        if cls._instance is None:
            cls._instance = super(Config, cls).__new__(cls)
            cls._instance._initialize()
        return cls._instance

    def _initialize(self):
        self.config_data = None
        self.config_path = None
        self.initialized = False

    @classmethod
    def reset(cls):
        """Reset the singleton instance for testing purposes."""
        cls._instance = None

    def load_config(self, config_path):
        """
        Load the configuration from a YAML file.
        :param config_path: Path to the configuration YAML file.
        :raise FileNotFoundError: If the config file doesn't exist.
        :raise yaml.YAMLError: If there's an error parsing the YAML file.
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
        self.initialized = True

    def _process_relative_paths(self, config_dir):
        """
        Convert relative paths in the config to absolute paths.
        :param: config_dir (str): Directory of the config file.
        """
        for key, value in self.config_data.items():
            if isinstance(value, str) and not os.path.isabs(value) and key.endswith('_file'):
                self.config_data[key] = os.path.join(config_dir, value)

    def setup_logging(self, cmd_log_level=None):
        """
        Set up logging based on the configuration file and command line argument.
        :param: cmd_log_level (str, optional): Log level from command line, overrides config file if provided.
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
        :param: key (str): The configuration key to retrieve.
        :return: The value associated with the key, or the default value if not found.
        """
        return self.config_data.get(key, default) if self.config_data else default

    def __getitem__(self, key):
        """
        Allow dictionary-style access to configuration values.
        :param: key (str): The configuration key to retrieve.
        :return: The value associated with the key.
        :raises: KeyError: If the key is not found in the configuration.
        """
        if self.config_data is None:
            raise ValueError("Configuration has not been loaded. Call load_config first.")
        return self.config_data[key]

    def __contains__(self, key):
        """
        Check if a key exists in the configuration.
        :param: key (str): The configuration key to check.
        :return: bool: True if the key exists, False otherwise.
        """
        return self.config_data is not None and key in self.config_data

    def __str__(self):
        """
        Return a string representation of the Config object.
        :return: str: A string representation of the Config object.
        """
        return (f"Config Object:\n"
                f"  Initialized: {self.initialized}\n"
                f"  Config Path: {self.config_path}\n"
                f"  Config Data: {self.config_data}")

    def __repr__(self):
        """
        Return a string representation of the Config object.
        :return: str: A string representation of the Config object.
        """
        return (f"Config(_initialized={self.initialized}, "
                f"config_path='{self.config_path}', "
                f"config_data={self.config_data})")
