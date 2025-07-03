from typing import Optional
from nbitk.config import Config
from barcode_validator.idservices.idservice import IDService
from barcode_validator.constants import RefDB
from barcode_validator.idservices.ncbi import NCBI
from barcode_validator.idservices.bold import BOLD

class IDServiceFactory:

    @staticmethod
    def create_idservice(config: Config, refdb: RefDB) -> Optional[IDService]:
        if refdb == RefDB.NCBI:
            return NCBI(config)
        elif refdb == RefDB.BOLD:
            return BOLD(config)
        return None