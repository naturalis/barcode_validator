import os
import sys
import yaml
import logging


class Config:
    _instance = None

    def __new__(cls):
        if cls._instance is None:
            cls._instance = super(Config, cls).__new__(cls)
        return cls._instance

    def __init__(self):
        if not hasattr(self, 'initialized') or not self.initialized:
            self.config_data = None
            self.config_path = None
            self.initialized = False

    def load_config(self, config_path):
        if self.initialized:
            return  # Config already loaded, do nothing

        self.config_path = config_path
        if not os.path.exists(config_path):
            raise FileNotFoundError(f"Config file not found: {config_path}")

        with open(config_path, 'r') as config_file:
            try:
                self.config_data = yaml.safe_load(config_file)
            except yaml.YAMLError as e:
                raise yaml.YAMLError(f"Error parsing config file: {e}")

        self._process_relative_paths(os.path.dirname(os.path.abspath(config_path)))
        self.initialized = True

    def _process_relative_paths(self, config_dir):
        for key, value in self.config_data.items():
            if isinstance(value, str) and not os.path.isabs(value) and key.endswith('_file'):
                self.config_data[key] = os.path.join(config_dir, value)

    def get(self, key, default=None):
        if not self.initialized:
            raise RuntimeError("Configuration not loaded. Call load_config first.")
        return self.config_data.get(key, default)

    def __getitem__(self, key):
        if not self.initialized:
            raise RuntimeError("Configuration not loaded. Call load_config first.")
        return self.config_data[key]

    def __contains__(self, key):
        if not self.initialized:
            raise RuntimeError("Configuration not loaded. Call load_config first.")
        return key in self.config_data

    def __str__(self):
        status = "Initialized" if self.initialized else "Not Initialized"
        return f"Config Object ({status}):\n  Config Path: {self.config_path}\n  Config Data: {self.config_data}"

    def __repr__(self):
        return f"Config(initialized={self.initialized}, config_path='{self.config_path}')"

    @classmethod
    def reset(cls):
        if cls._instance is not None:
            cls._instance = None

    def setup_logging(self, cmd_log_level=None):
        if not self.initialized:
            raise RuntimeError("Configuration not loaded. Call load_config first.")

        log_level = self.get('log_level', 'INFO')
        if cmd_log_level:
            log_level = cmd_log_level

        numeric_level = getattr(logging, log_level.upper(), None)
        if not isinstance(numeric_level, int):
            raise ValueError(f'Invalid log level: {log_level}')

        logging.basicConfig(level=numeric_level,
                            format='%(asctime)s - %(levelname)s - %(message)s',
                            filename=self.get('log_file'),
                            filemode='a')

        sys.excepthook = exception_handler


def exception_handler(exc_type, exc_value, exc_traceback):
    if issubclass(exc_type, KeyboardInterrupt):
        sys.__excepthook__(exc_type, exc_value, exc_traceback)
        return
    logging.error("Uncaught exception", exc_info=(exc_type, exc_value, exc_traceback))