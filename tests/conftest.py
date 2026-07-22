"""
Shared pytest fixtures and helpers for the barcode_validator test suite.

Several tests exercise taxonomic validation through the Galaxy backend, which submits sequences to the
'Identify reads with blastn and find taxonomy' tool on a remote Galaxy instance. That tool, and the curated
reference databases it searches, are currently deployed on galaxy.naturalis.nl only. Such tests therefore
cannot pass against any other instance, however valid the API key for it may be. The helpers below detect
that situation up front and skip with an explanatory message, rather than letting the run fail deep inside
the validation pipeline with an opaque error.

See the 'Galaxy Integration' section of README.md for the user-facing description of this restriction.
"""
import os

import pytest
import requests

# The Galaxy instance that hosts the tool on which taxonomic validation depends, and the timeout in seconds
# applied to the preflight requests made against it.
SUPPORTED_GALAXY_DOMAIN = 'galaxy.naturalis.nl'
GALAXY_TIMEOUT = 30


def skip_unless_supported_domain(domain: str) -> None:
    """
    Skip the calling test unless the configured Galaxy instance is the one hosting the required tool.

    :param domain: Value of the GALAXY_DOMAIN environment variable
    :return: None
    """
    if domain != SUPPORTED_GALAXY_DOMAIN:
        pytest.skip(
            f"GALAXY_DOMAIN is '{domain}', but the BLASTN tool that taxonomic validation depends on is "
            f"currently deployed only on '{SUPPORTED_GALAXY_DOMAIN}'. Public instances such as "
            f"usegalaxy.org and usegalaxy.eu do not host it, so an API key for those servers will not make "
            f"these tests pass. Use '--taxon-validation method=bold' for a backend that needs no Galaxy "
            f"account; see the 'Galaxy Integration' section of README.md."
        )


def skip_unless_valid_api_key(domain: str, api_key: str) -> None:
    """
    Skip the calling test unless the API key authenticates against the Galaxy instance.

    Galaxy's /api/users/current endpoint returns the authenticated user for a valid key, and an error status
    for one that is absent, malformed, revoked or issued by a different server. This makes it a cheap way to
    establish that the credentials are usable before any sequences are submitted.

    :param domain: Domain of the Galaxy instance, e.g. galaxy.naturalis.nl
    :param api_key: Galaxy API key to verify
    :return: None
    """
    url = f"https://{domain}/api/users/current"
    try:
        response = requests.get(url, headers={'x-api-key': api_key}, timeout=GALAXY_TIMEOUT)
    except requests.RequestException as e:
        pytest.skip(f"Cannot reach Galaxy instance '{domain}': {e}")
    if response.status_code in (400, 401, 403):
        pytest.skip(
            f"GALAXY_API_KEY was rejected by '{domain}' (HTTP {response.status_code}). Generate a fresh key "
            f"under User -> Preferences -> Manage API Key on that instance."
        )
    if not response.ok:
        pytest.skip(f"Galaxy instance '{domain}' returned HTTP {response.status_code} for {url}.")


@pytest.fixture
def galaxy_credentials() -> None:
    """
    Establish that a usable Galaxy instance and API key are available, skipping the test otherwise.

    This fixture is deliberately not autouse: tests that do not touch the Galaxy backend must keep running
    without credentials. Modules that do touch it opt in with a module-level
    `pytestmark = pytest.mark.usefixtures('galaxy_credentials')`.

    :return: None
    """
    domain = os.environ.get('GALAXY_DOMAIN')
    if not domain:
        pytest.skip("GALAXY_DOMAIN not set in environment")
    api_key = os.environ.get('GALAXY_API_KEY')
    if not api_key:
        pytest.skip("GALAXY_API_KEY not set in environment")
    skip_unless_supported_domain(domain)
    skip_unless_valid_api_key(domain, api_key)
