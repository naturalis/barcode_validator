import argparse
import requests
import time
import os
import subprocess
import logging
import sqlite3
import traceback
import base64
from datetime import datetime, timedelta
from config import Config
from barcode_validator.core import BarcodeValidator

GITHUB_TOKEN = os.environ.get('GITHUB_TOKEN')
POLL_INTERVAL = 300  # 5 minutes

headers = {
    'Authorization': f'token {GITHUB_TOKEN}',
    'Accept': 'application/vnd.github.v3+json'
}


def setup_database(config):
    conn = sqlite3.connect(config.get('pr_db_file'))
    c = conn.cursor()
    c.execute('''CREATE TABLE IF NOT EXISTS prs
                 (pr_number INTEGER PRIMARY KEY, status TEXT, last_updated TIMESTAMP)''')
    conn.commit()
    return conn


def get_open_prs(config):
    """
    Get a list of open pull requests for the repository from the remote location on GitHub via its API. The location
    of the pull requests endpoint is constructed from information held by the Config object.
    :return: A list of pull request
    """
    url = f'https://api.github.com/repos/{config.get("repo_owner")}/{config.get("repo_name")}/pulls'
    response = requests.get(url, headers=headers)
    return response.json()


def get_pr_files(config, pr_number):
    """
    Get a list of files for a given pull request. The location of the files-for-this-pull-request endpoint is
    constructed from the repository owner and name (which are held by the Config object) as well as the pull request
    number, which is passed as an argument.
    :param config: The Config object
    :param pr_number: The pull request number
    :return: A list of file objects
    """
    url = f'https://api.github.com/repos/{config.get("repo_owner")}/{config.get("repo_name")}/pulls/{pr_number}/files'
    response = requests.get(url, headers=headers)
    return response.json()


def run_validation(config, pr_number, branch, validator):
    """
    Run the validation process for a given pull request.
    :param config: The Config object
    :param pr_number: The pull request number
    :param branch: The branch name
    :param validator: A BarcodeValidator object
    :return: A list of validation results
    """
    logging.info(f"Starting validation for PR {pr_number}")

    # Go to local clone and fetch the PR
    repo_location = config.get('repo_location')
    logging.info(f"Changing directory to {repo_location}")
    os.chdir(repo_location)

    # Fetch the latest changes
    logging.info(f"Fetching latest changes for PR {pr_number}")
    run_git_command(['git', 'fetch', 'origin'], "Failed to fetch from origin")

    # Create or reset the PR branch
    pr_branch = f"pr-{pr_number}"
    logging.info(f"Creating/resetting branch {pr_branch}")
    run_git_command(['git', 'checkout', '-B', pr_branch, f'origin/{branch}'],
                    f"Failed to create/reset branch {pr_branch}")

    # Get the FASTA files from the PR
    logging.info(f"Getting files for PR {pr_number}")
    files = get_pr_files(config, pr_number)
    fasta_files = [f for f in files if f['filename'].endswith(('.fasta', '.fa', '.fas'))]
    logging.info(f"Found {len(fasta_files)} FASTA files in PR {pr_number}")

    # Run the validation process for each FASTA file
    all_results = []
    for file in fasta_files:
        logging.info(f"Processing file: {file['filename']}")

        # Fetch file content
        file_url = file['raw_url']
        logging.info(f"Fetching file from {file_url}")
        response = requests.get(file_url, headers=headers)
        if response.status_code == 200:
            # Create directory if it doesn't exist
            os.makedirs(os.path.dirname(file['filename']), exist_ok=True)

            # Save file content
            logging.info(f"Saving file content to {file['filename']}")
            with open(file['filename'], 'wb') as f:
                f.write(response.content)

            # Validate file
            logging.info(f"Validating {file['filename']}...")
            results = validator.validate_fasta(file['filename'])
            logging.info(f"Validation complete for {file['filename']}")
            all_results.extend([(file['filename'], result) for result in results])
        else:
            logging.error(f"Failed to fetch file {file['filename']}: {response.status_code}")

    logging.info(f"Validation complete for PR {pr_number}")
    return all_results


def post_results(config, pr_number, results):
    """
    Post a comment to a pull request with the validation results.
    :param config: The Config object
    :param pr_number: The pull request number
    :param results: A list of validation results
    :return: None
    """
    comment = "## Validation Results\n\n"
    for file, result in results:
        comment += f"## File: {file}\n"
        barcode_rank, full_rank, messages = result.calculate_ranks(verbosity=3)
        status_emoji = "✅" if result.passes_all_checks() else "❗"
        obs_taxon_names = "\n".join(f"    - {taxon.name}" for taxon in result.obs_taxon)
        comment += f"""
### {status_emoji} Process ID: {result.process_id}
- Taxonomic check passed: **{result.check_taxonomy()}**
  - Expected species as registered at BOLD: {result.species}
  - Expected {config.get('level')} as registered at BOLD: {result.exp_taxon}
  - Observed BLAST hits at {config.get('level')} level: 
{obs_taxon_names}  
- BIN-compliant sequence length passed: **{result.check_length()}**
  - Net length aligned to marker region: {result.seq_length}  
  - Full sequence length: {result.full_length}
- Sequence quality check passed: **{result.check_seq_quality()}**
  - Ambiguities in marker region: {result.ambiguities}
  - Ambiguities in full sequence: {result.full_ambiguities}
  - Stop codons: {len(result.stop_codons)}
- Rankings:
  - Barcode rank: {barcode_rank}
  - Full rank: {full_rank}
  - Messages: 
    {messages}
"""
    post_comment(config, pr_number, comment)


def post_comment(config, pr_number, comment):
    """
    Post a comment to a pull request.
    :param config: The Config object
    :param pr_number: The pull request number
    :param comment: The comment to post
    :return: None
    """
    url = f'https://api.github.com/repos/{config.get("repo_owner")}/{config.get("repo_name")}/issues/{pr_number}/comments'
    payload = {'body': comment}
    requests.post(url, headers=headers, json=payload)


def run_git_command(command, error_message):
    """Run a git command and log any errors."""
    result = subprocess.run(command, capture_output=True, text=True)
    if result.returncode != 0:
        logging.error(f"{error_message}: {result.stderr}")
        raise RuntimeError(f"Git command failed: {' '.join(command)}")
    return result.stdout


def process_pr(config, validator, conn, pr_number, branch):
    c = conn.cursor()
    c.execute("SELECT status FROM prs WHERE pr_number = ?", (pr_number,))
    row = c.fetchone()

    if row is None or row[0] == 'pending':
        try:

            # Start processing the PR
            c.execute("INSERT OR REPLACE INTO prs (pr_number, status, last_updated) VALUES (?, 'processing', ?)",
                      (pr_number, datetime.now()))
            conn.commit()
            logging.info(f"Changed status of PR {pr_number} to 'processing'")
            post_comment(config, pr_number, "\U0001F916 - Hi! This is an automated message from the barcode validation "
                                            "robot. I'm going to validate the FASTA files in your request. Please wait "
                                            "while I process the files. This takes about two minutes per "
                                            "sequence.")

            # Run the validation process
            results = run_validation(config, pr_number, branch, validator)
            post_results(config, pr_number, results)

            # Update the PR status
            c.execute("UPDATE prs SET status = 'completed', last_updated = ? WHERE pr_number = ?",
                      (datetime.now(), pr_number))
            conn.commit()
            logging.info(f"Changed status of PR {pr_number} to 'completed'")
            post_comment(config, pr_number, "\U0001F916 - Validation complete. If all looks good, "
                                            "notify @rvosa to merge this PR.")

        except Exception as e:
            error_msg = f"Error processing PR {pr_number}: {str(e)}\n"
            error_msg += "Stack trace:\n"
            error_msg += traceback.format_exc()
            logging.error(error_msg)
            c.execute("UPDATE prs SET status = 'error', last_updated = ? WHERE pr_number = ?",
                      (datetime.now(), pr_number))
            conn.commit()
            post_comment(config, pr_number, f"\U0001F916 - {error_msg}")


def main(config_file, verbosity):
    """
    Main function for the barcode validator daemon.
    :param config_file: Path to the configuration file
    :param verbosity: Logging verbosity level
    :return: None
    """

    # Initialize the Config object and setup logging
    config = Config()
    config.load_config(config_file)
    config.setup_logging(verbosity)
    logging.info("*** Barcode Validator Daemon starting ***")

    # Set up the database connection for tracking PRs
    try:
        logging.info(f"Going to initialize PR database at {config.get('pr_db_file')}")
        conn = setup_database(config)
        logging.info("Database initialized")
    except Exception as e:
        logging.error(f"Error setting up database: {str(e)}")
        exit(1)

    # Initialize the BarcodeValidator object
    logging.info("Initializing BarcodeValidator...")
    validator = BarcodeValidator(config)
    validator.initialize()

    # Start the main loop
    logging.info("Starting main loop")
    while True:

        # Get a list of open PRs and iterate over them
        logging.info("Checking for open PRs...")
        prs = get_open_prs(config)
        logging.info(f"Found {len(prs)} open PRs")

        for pr in prs:

            # Process PRs that contain FASTA files
            pr_number = pr['number']
            logging.info(f"Inspecting files from PR {pr['number']}...")
            files = get_pr_files(config, pr_number)
            logging.info(f"Found {len(files)} files in PR {pr['number']}")

            # Iterate over the files in the PR and process if any are FASTA files
            if any(f['filename'].endswith(('.fasta', '.fa', '.fas')) for f in files):
                logging.info(f"Processing PR {pr['number']}")
                process_pr(config, validator, conn, pr_number, pr['head']['ref'])

        # Clean up old completed PRs
        #        c = conn.cursor()
        #        c.execute("DELETE FROM prs WHERE status = 'completed' AND last_updated < ?",
        #                  (datetime.now() - timedelta(days=7),))
        #        conn.commit()

        time.sleep(POLL_INTERVAL)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="GitHub PR Validation Daemon for barcode_validator")
    parser.add_argument("-c", "--config_file", required=True, help="Path to the configuration YAML file")
    parser.add_argument("-v", "--verbosity", choices=['DEBUG', 'INFO', 'WARNING', 'ERROR', 'CRITICAL'],
                        default='INFO', help="Set the logging verbosity (default: INFO)")
    args = parser.parse_args()

    main(args.config_file, args.verbosity)
