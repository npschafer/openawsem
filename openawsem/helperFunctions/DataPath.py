import configparser
from pathlib import Path

class DataPath:
    def __init__(self, default_location:Path, default_config_path: Path, custom_config_path: Path):
        self.default_location=default_location
        self._load_config(default_config_path, custom_config_path)    

    def _load_config(self, default_config_path: Path, custom_config_path: Path):
        config = configparser.ConfigParser()

        # Load default paths
        config.read(default_config_path)
        for data_type in config['Data Paths']:
            default_path = self.default_location / config.get('Data Paths', data_type)
            setattr(self, data_type, default_path)

        # Override with custom paths if available
        if custom_config_path.exists():
            print('Using Custom config file')
            config.read(custom_config_path)
            for data_type in config['Data Paths']:
                custom_path = Path(config.get('Data Paths', data_type))
                if custom_path.exists():
                    setattr(self, data_type, custom_path)
