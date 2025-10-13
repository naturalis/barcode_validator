import argparse
import yaml
import importlib.resources
from pathlib import Path
from typing import Any, Dict, List, Optional, Union
from nbitk.config import Config


class ValidationError(Exception):
    """Exception raised when configuration validation fails."""
    pass


class SchemaConfig(Config):
    """
    Schema-driven configuration class that extends the base Config class.

    This class loads a YAML schema that defines all configuration parameters,
    their types, defaults, choices, and help text. It can automatically
    generate argparse parsers and validate configuration values.

    Examples:
        >>> config = SchemaConfig()
        >>> parser = argparse.ArgumentParser()
        >>> config.populate_argparse(parser)
        >>> config.set('marker', 'COI-5P')
        >>> config.get('marker')
        'COI-5P'
    """

    def __init__(self, schema_package: str = "barcode_validator.config"):
        """
        Initialize the schema-driven configuration.

        :param schema_package: Package containing the schema.yaml file
        """
        # Initialize parent without loading any config
        super().__init__()
        self.schema_package = schema_package
        self.schema: Dict[str, Any] = {}
        self._load_schema()
        self._initialize_with_defaults()

    def _load_schema(self) -> None:
        """Load the schema from the package resources."""
        try:
            # Use importlib.resources to load schema.yaml from the package
            with importlib.resources.open_text(self.schema_package, "schema.yaml") as f:
                self.schema = yaml.safe_load(f)
        except FileNotFoundError:
            raise FileNotFoundError(f"Schema file not found in package {self.schema_package}")
        except yaml.YAMLError as e:
            raise yaml.YAMLError(f"Error parsing schema file: {e}")

        # Validate the schema structure
        self._validate_schema()

    def _validate_schema(self) -> None:
        """Validate that the loaded schema has the correct structure."""
        for key, spec in self.schema.items():
            if not isinstance(spec, dict):
                raise ValidationError(f"Schema entry '{key}' must be a dictionary")

            # Check if this is a container (has nested parameters but no type)
            is_container = self._is_nested_config(spec)

            if is_container:
                # This is a container - validate its nested parameters
                self._validate_container_schema(key, spec)
            else:
                # This is a regular parameter - validate it has required fields
                self._validate_parameter_schema(key, spec)

    def _validate_container_schema(self, container_key: str, spec: Dict[str, Any]) -> None:
        """Validate a container schema entry and its nested parameters."""
        standard_keys = {'type', 'default', 'choices', 'required', 'help'}

        # Container should not have a type
        if 'type' in spec:
            raise ValidationError(f"Container schema entry '{container_key}' should not have a 'type' field")

        # Validate each nested parameter
        for nested_key, nested_spec in spec.items():
            if nested_key in standard_keys:
                continue  # Skip standard keys like 'help'

            if not isinstance(nested_spec, dict):
                raise ValidationError(f"Nested parameter '{container_key}.{nested_key}' must be a dictionary")

            # Validate the nested parameter
            self._validate_parameter_schema(f"{container_key}.{nested_key}", nested_spec)

    def _validate_parameter_schema(self, param_key: str, spec: Dict[str, Any]) -> None:
        """Validate a parameter schema entry."""
        if 'type' not in spec:
            raise ValidationError(f"Parameter schema entry '{param_key}' missing required 'type' field")

        valid_types = ['str', 'int', 'float', 'bool', 'Path']
        if spec['type'] not in valid_types:
            raise ValidationError(f"Parameter schema entry '{param_key}' has invalid type '{spec['type']}'. "
                                  f"Valid types: {valid_types}")

        # Validate choices if present
        if 'choices' in spec and not isinstance(spec['choices'], list):
            raise ValidationError(f"Parameter schema entry '{param_key}' choices must be a list")

    def _initialize_with_defaults(self) -> None:
        """Initialize configuration with default values from schema."""
        self.config_data = {}
        self.initialized = True

        # Set defaults for all schema entries
        for key, spec in self.schema.items():
            if isinstance(spec, dict) and 'default' in spec:
                # Handle nested configurations
                if self._is_nested_config(spec):
                    nested_defaults = {}
                    for nested_key, nested_spec in spec.items():
                        if nested_key not in ['type', 'help'] and isinstance(nested_spec, dict):
                            if 'default' in nested_spec:
                                nested_defaults[nested_key] = nested_spec['default']
                    if nested_defaults:
                        self.config_data[key] = nested_defaults
                else:
                    self.config_data[key] = spec['default']

    def _is_nested_config(self, spec: Dict[str, Any]) -> bool:
        """Check if a schema specification represents a nested configuration."""
        # A spec is nested if it contains keys other than the standard ones
        # and those keys are themselves dictionaries with type/default/etc.
        standard_keys = {'type', 'default', 'choices', 'required', 'help'}
        other_keys = set(spec.keys()) - standard_keys

        for key in other_keys:
            if isinstance(spec[key], dict) and 'type' in spec[key]:
                return True
        return False

    def load_config(self, config_path: str) -> None:
        """
        Override parent method to prevent loading external config files.

        :param config_path: Path to config file (not used)
        :raises ValidationError: Always, since external config loading is disabled
        """
        raise ValidationError("Loading external configuration files is not supported. "
                              "Use command line arguments or set() method instead.")

    def set(self, key: str, value: Any) -> None:
        """
        Set a configuration value with validation.

        :param key: Configuration key (supports dot notation for nested keys)
        :param value: Configuration value
        :raises ValidationError: If key is unknown or value is invalid
        """
        # Handle nested keys (e.g., "taxonomic_validation.min_identity")
        if '.' in key:
            main_key, nested_key = key.split('.', 1)
            if main_key not in self.schema:
                raise ValidationError(f"Unknown configuration key: {main_key}")

            # Get or create nested dict
            if main_key not in self.config_data:
                self.config_data[main_key] = {}

            # Validate nested key exists in schema
            if not self._is_nested_config(self.schema[main_key]):
                raise ValidationError(f"Key '{main_key}' does not support nested configuration")

            if nested_key not in self.schema[main_key]:
                raise ValidationError(f"Unknown nested configuration key: {key}")

            # Validate the nested value
            nested_spec = self.schema[main_key][nested_key]
            validated_value = self._validate_value(nested_key, value, nested_spec)
            self.config_data[main_key][nested_key] = validated_value

        else:
            # Handle top-level keys
            if key not in self.schema:
                raise ValidationError(f"Unknown configuration key: {key}")

            spec = self.schema[key]
            validated_value = self._validate_value(key, value, spec)
            self.config_data[key] = validated_value

        # Handle interdependencies
        self._handle_dependencies(key, value)

    def _validate_value(self, key: str, value: Any, spec: Dict[str, Any]) -> Any:
        """
        Validate a configuration value against its schema specification.

        :param key: Configuration key name
        :param value: Value to validate
        :param spec: Schema specification for the key
        :return: Validated and converted value
        :raises ValidationError: If validation fails
        """
        if value is None and spec.get('required', False):
            raise ValidationError(f"Required configuration key '{key}' cannot be None")

        if value is None:
            return None

        # Type conversion and validation
        expected_type = spec['type']

        try:
            if expected_type == 'str':
                converted_value = str(value)
            elif expected_type == 'int':
                converted_value = int(value)
            elif expected_type == 'float':
                converted_value = float(value)
            elif expected_type == 'bool':
                if isinstance(value, str):
                    converted_value = value.lower() in ('true', 'yes', '1', 'on')
                else:
                    converted_value = bool(value)
            elif expected_type == 'Path':
                converted_value = Path(value)
            else:
                raise ValidationError(f"Unknown type '{expected_type}' for key '{key}'")

        except (ValueError, TypeError) as e:
            raise ValidationError(f"Cannot convert value '{value}' to type '{expected_type}' "
                                  f"for key '{key}': {e}")

        # Validate choices if specified
        if 'choices' in spec and converted_value not in spec['choices']:
            raise ValidationError(f"Invalid value '{converted_value}' for key '{key}'. "
                                  f"Valid choices: {spec['choices']}")

        return converted_value

    def _handle_dependencies(self, key: str, value: Any) -> None:
        """
        Handle configuration dependencies and side effects.

        :param key: Configuration key that was set
        :param value: Value that was set
        """
        # Handle mode dependencies
        if key == 'mode':
            if value == 'structural':
                self.config_data['validate_structure'] = True
                self.config_data['validate_taxonomy'] = False
            elif value == 'taxonomic':
                self.config_data['validate_structure'] = False
                self.config_data['validate_taxonomy'] = True
            elif value == 'both':
                self.config_data['validate_structure'] = True
                self.config_data['validate_taxonomy'] = True

    def get(self, key: str, default: Any = None) -> Any:
        """
        Get a configuration value, supporting dot notation for nested keys.

        :param key: Configuration key (supports dot notation)
        :param default: Default value if key not found
        :return: Configuration value
        """
        if '.' in key:
            main_key, nested_key = key.split('.', 1)
            if main_key in self.config_data and isinstance(self.config_data[main_key], dict):
                return self.config_data[main_key].get(nested_key, default)
            return default
        else:
            return self.config_data.get(key, default)

    def populate_argparse(self, parser: argparse.ArgumentParser) -> None:
        """
        Populate an argparse parser with arguments from the schema.

        :param parser: ArgumentParser to populate
        """
        for key, spec in self.schema.items():
            if self._is_nested_config(spec):
                # Handle nested configurations
                self._add_nested_arguments(parser, key, spec)
            else:
                # Handle regular arguments
                self._add_argument(parser, key, spec)

    def _add_nested_arguments(self, parser: argparse.ArgumentParser,
                              main_key: str, spec: Dict[str, Any]) -> None:
        """Add nested configuration arguments to parser."""
        help_text = spec.get('help', f'Configuration for {main_key}')

        # Create a custom action for nested parameters
        class NestedConfigAction(argparse.Action):
            def __init__(self, option_strings, dest, **kwargs):
                super().__init__(option_strings, dest, **kwargs)
                self.main_key = main_key
                self.spec = spec

            def __call__(self, parser, namespace, values, option_string=None):
                # Parse key=value format
                if '=' not in values:
                    raise argparse.ArgumentTypeError(
                        f"Expected format: nested_key=value, got: {values}")

                nested_key, nested_value = values.split('=', 1)

                # Validate nested key exists in schema
                if nested_key not in self.spec:
                    available_keys = [k for k in self.spec.keys()
                                      if k not in ['type', 'help']]
                    raise argparse.ArgumentTypeError(
                        f"Unknown nested key '{nested_key}' for --{self.main_key}. "
                        f"Available keys: {available_keys}")

                # Initialize nested dict if needed
                if getattr(namespace, self.main_key) is None:
                    setattr(namespace, self.main_key, {})

                nested_dict = getattr(namespace, self.main_key)
                nested_dict[nested_key] = nested_value

        parser.add_argument(
            f'--{main_key.replace("_", "-")}',
            action=NestedConfigAction,
            dest=main_key,  # Use original key as destination
            help=help_text,
            metavar='KEY=VALUE'
        )

    def _add_argument(self, parser: argparse.ArgumentParser,
                      key: str, spec: Dict[str, Any]) -> None:
        """Add a single argument to the parser."""
        arg_name = f'--{key.replace("_", "-")}'
        kwargs = {
            'help': spec.get('help', f'Set {key}'),
            'dest': key
        }

        # Handle different types
        if spec['type'] == 'bool':
            if spec.get('default', False):
                # If default is True, create a --no-flag option
                kwargs['action'] = 'store_false'
                arg_name = f'--no-{key.replace("_", "-")}'
            else:
                kwargs['action'] = 'store_true'
        else:
            kwargs['type'] = self._get_argparse_type(spec['type'])
            if 'default' in spec:
                kwargs['default'] = spec['default']
            if spec.get('required', False):
                kwargs['required'] = True

        # Handle choices
        if 'choices' in spec:
            kwargs['choices'] = spec['choices']

        parser.add_argument(arg_name, **kwargs)

    def _get_argparse_type(self, schema_type: str):
        """Convert schema type to argparse type."""
        if schema_type == 'str':
            return str
        elif schema_type == 'int':
            return int
        elif schema_type == 'float':
            return float
        elif schema_type == 'Path':
            return Path
        else:
            return str

    def update_from_args(self, args: argparse.Namespace) -> None:
        """
        Update configuration from parsed command line arguments.

        :param args: Parsed arguments from argparse
        """
        # Update regular arguments
        for key in self.schema.keys():
            if hasattr(args, key) and getattr(args, key) is not None:
                if self._is_nested_config(self.schema[key]):
                    # Handle nested configurations
                    nested_dict = getattr(args, key, {})
                    for nested_key, nested_value in nested_dict.items():
                        self.set(f'{key}.{nested_key}', nested_value)
                else:
                    self.set(key, getattr(args, key))

    def validate_required_fields(self) -> None:
        """
        Validate that all required fields have been set.

        :raises ValidationError: If required fields are missing
        """
        missing_fields = []

        for key, spec in self.schema.items():
            if spec.get('required', False):
                if key not in self.config_data or self.config_data[key] is None:
                    missing_fields.append(key)

            # Check nested required fields
            if self._is_nested_config(spec):
                for nested_key, nested_spec in spec.items():
                    if (isinstance(nested_spec, dict) and
                            nested_spec.get('required', False)):
                        full_key = f'{key}.{nested_key}'
                        if self.get(full_key) is None:
                            missing_fields.append(full_key)

        if missing_fields:
            raise ValidationError(f"Required fields missing: {', '.join(missing_fields)}")

    def get_schema_help(self, key: str) -> str:
        """
        Get help text for a configuration key.

        :param key: Configuration key
        :return: Help text
        """
        if '.' in key:
            main_key, nested_key = key.split('.', 1)
            if (main_key in self.schema and
                    self._is_nested_config(self.schema[main_key]) and
                    nested_key in self.schema[main_key]):
                return self.schema[main_key][nested_key].get('help', 'No help available')
        elif key in self.schema:
            return self.schema[key].get('help', 'No help available')

        return 'Unknown configuration key'