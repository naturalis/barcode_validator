from typing import Optional
from nbitk.config import Config
from barcode_validator.idservices.idservice import IDService
from barcode_validator.constants import RefDB


class IDServiceFactory:

    @staticmethod
    def create_idservice(config: Config, refdb: RefDB) -> Optional[IDService]:
        from barcode_validator.idservices.ncbi import NCBI
        if refdb.value == "ncbi":
            return NCBI(config)