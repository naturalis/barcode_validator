from typing import Optional
from nbitk.config import Config
from barcode_validator.idservices.idservice import IDService
from barcode_validator.constants import RefDB
from barcode_validator.idservices.blast import BLAST as LocalBLAST
from barcode_validator.idservices.bold import BOLD as BOLDIDService
from barcode_validator.idservices.galaxy_blast import GalaxyBLAST
from barcode_validator.idservices.bolddistilled import BLAST as BOLDDistilled

class IDServiceFactory:

    @staticmethod
    def create_idservice(config: Config, refdb: RefDB) -> Optional[IDService]:
        if refdb == RefDB.BLAST:
            return LocalBLAST(config)
        elif refdb == RefDB.BOLD:
            return BOLDIDService(config)
        elif refdb == RefDB.GALAXY:
            return GalaxyBLAST(config)
        elif refdb == RefDB.BOLDDISTILLED:
            return BOLDDistilled(config)
        return None