import argparse
import requests
import time
import os
import logging
import sqlite3
import traceback
from datetime import datetime
from typing import Optional
from barcode_validator.config import Config
from barcode_validator.core import BarcodeValidator
from barcode_validator.github import GitHubClient


GITHUB_TOKEN = os.environ.get('GITHUB_TOKEN')
POLL_INTERVAL = 300  # 5 minutes


class ValidationDaemon:

    def __init__(self):
        self.bv: Optional[BarcodeValidator] = None
        self.gc: Optional[GitHubClient] = None
        self.conn: Optional[sqlite3.Connection] = None

    def initialize(self, config: Config):
        self.bv = BarcodeValidator()
        self.bv.initialize(config.get('ncbi_taxonomy'), config.get('bold_sheet_file'))
        self.gc = GitHubClient(config.get('repo_owner'), config.get('repo_name'), GITHUB_TOKEN,
                               config.get('repo_location'))
        try:
            logging.info(f"Going to initialize PR database at {config.get('pr_db_file')}")
            self.conn = self.setup_database(config.get('pr_db_file'))
            logging.info("Database initialized")
        except Exception as e:
            logging.error(f"Error setting up database: {str(e)}")
            exit(1)

    @classmethod
    def setup_database(cls, db_file):
        """
        Set up the SQLite database for storing PR status.
        :param db_file: The path to the database file
        :return: A connection object
        """
        conn = sqlite3.connect(db_file)
        c = conn.cursor()
        c.execute('''CREATE TABLE IF NOT EXISTS prs
                     (pr_number INTEGER PRIMARY KEY, status TEXT, last_updated TIMESTAMP)''')
        conn.commit()
        return conn

    def process_pr(self, config, pr_number, branch):
        """
        Process a pull request:
        - Initialize the PR by changing status to 'processing' and notifying the user
        - Validate the PR by fetching its FASTA files and running the validation process
        - Post the validation results to the PR and push the results to the PR branch
        - Finalize the PR by changing status to 'completed' and notifying the user
        :param config: The Config object
        :param pr_number: The pull request number
        :param branch: The branch name
        """
        c = self.conn.cursor()
        c.execute("SELECT status FROM prs WHERE pr_number = ?", (pr_number,))
        row = c.fetchone()

        if row is None or row[0] == 'pending':
            try:
                # Change pending => processing, notify user
                self.initialize_pr(pr_number)

                # Run the validation process
                results = self.validate_pr(config, pr_number, branch)
                self.post_pr_results(config, pr_number, results)

                # Change processing => completed, notify user
                self.finalize_pr(pr_number)

            except Exception as e:
                error_msg = f"Error processing PR {pr_number}: {str(e)}\n"
                error_msg += "Stack trace:\n"
                error_msg += traceback.format_exc()
                logging.error(error_msg)
                c.execute("UPDATE prs SET status = 'error', last_updated = ? WHERE pr_number = ?",
                          (datetime.now(), pr_number))
                self.conn.commit()
                self.gc.post_comment(pr_number, f"\U0001F916 - {error_msg}")

    def initialize_pr(self, pr_number):
        """
        Initialize a pull request by changing its status to 'processing' and notifying the user.
        :param pr_number: The pull request number
        """
        # Start processing the PR
        c = self.conn.cursor()
        c.execute("INSERT OR REPLACE INTO prs (pr_number, status, last_updated) VALUES (?, 'processing', ?)",
                  (pr_number, datetime.now()))
        self.conn.commit()
        logging.info(f"Changed status of PR {pr_number} to 'processing'")
        self.gc.post_comment(pr_number, "\U0001F916 - Hi! This is an automated message from the barcode validation "
                                        "robot. I'm going to validate the FASTA files in your request. Please wait "
                                        "while I process the files. This takes about two minutes per sequence.")

    def finalize_pr(self, pr_number):
        """
        Finalize a pull request by changing its status to 'completed' and notifying the user.
        :param pr_number: The pull request number
        """
        # Update the PR status
        c = self.conn.cursor()
        c.execute("UPDATE prs SET status = 'completed', last_updated = ? WHERE pr_number = ?",
                  (datetime.now(), pr_number))
        self.conn.commit()
        logging.info(f"Changed status of PR {pr_number} to 'completed'")
        self.gc.post_comment(pr_number, "\U0001F916 - Validation complete. If all looks good, notify @rvosa to merge.")

    def validate_pr(self, config, pr_number, branch):
        """
        Run the validation process for a given pull request.
        :param config: The Config object
        :param pr_number: The pull request number
        :param branch: The branch name
        :return: A list of validation results
        """
        logging.info(f"Starting validation for PR {pr_number}")
        fasta_files = self.fetch_pr_fastas(branch, pr_number)

        # Run the validation process for each FASTA file
        all_results = []
        for file in fasta_files:
            logging.info(f"Processing file: {file['filename']}")

            # Fetch file content
            file_url = file['raw_url']
            logging.info(f"Fetching file from {file_url}")
            response = requests.get(file_url, headers=self.gc.headers)
            if response.status_code == 200:
                # Create directory if it doesn't exist
                os.makedirs(os.path.dirname(file['filename']), exist_ok=True)

                # Save file content
                logging.info(f"Saving file content to {file['filename']}")
                with open(file['filename'], 'wb') as f:
                    f.write(response.content)

                # Validate file
                logging.info(f"Validating {file['filename']}...")
                results = self.bv.validate_fasta(file['filename'], config)
                logging.info(f"Validation complete for {file['filename']}")
                all_results.extend([(file['filename'], result) for result in results])
            else:
                logging.error(f"Failed to fetch file {file['filename']}: {response.status_code}")

        logging.info(f"Validation complete for PR {pr_number}")
        return all_results

    def fetch_pr_fastas(self, branch, pr_number):
        """
        Fetch the FASTA files from a pull request.
        :param branch: The branch name
        :param pr_number: The pull request number
        :return: A list of FASTA files
        """

        # Fetch the latest changes
        logging.info(f"Fetching latest changes for PR {pr_number}")
        self.gc.run_git_command(['git', 'fetch', 'origin'], "Failed to fetch from origin")

        # Create or reset the PR branch
        pr_branch = f"pr-{pr_number}"
        logging.info(f"Creating/resetting branch {pr_branch}")
        self.gc.run_git_command(['git', 'checkout', '-B', pr_branch, f'origin/{branch}'],
                                f"Failed to create/reset branch {pr_branch}")

        # Get the FASTA files from the PR
        logging.info(f"Getting files for PR {pr_number}")
        files = self.gc.get_pr_files(pr_number)
        fasta_files = [f for f in files if f['filename'].endswith(('.fasta', '.fa', '.fas'))]
        logging.info(f"Found {len(fasta_files)} FASTA files in PR {pr_number}")
        return fasta_files

    def post_pr_results(self, config, pr_number, results):
        """
        Post a comment to a pull request with the validation results.
        :param config: The Config object
        :param pr_number: The pull request number
        :param results: A list of validation results
        :return: None
        """
        comment = "# Validation Results\n\n"
        current_file_handle = None
        current_file_name = None
        for file, r in results:

            # If this is true, then current_file_name is None or previous file
            if file != current_file_name:

                # Close the previous file handle if it exists, i.e. current_file_name was not None
                if current_file_handle:
                    # Finalize the file and commit it
                    current_file_handle.close()
                    tsv_name = f"{current_file_name}.tsv"
                    self.gc.commit_file(current_file_name, f"Add validated FASTA file {current_file_name}")
                    self.gc.commit_file(tsv_name, f"Add validation results for {current_file_name}")

                # Open the new file and write the header, happens both when current_file_name is None or previous file
                current_file_handle = open(f"{file}.tsv", 'w', buffering=1)
                hlist = r.result_fields(config.get('level'))
                hlist.append('fasta_file')
                current_file_handle.write('\t'.join(hlist) + '\n')
                current_file_name = file

            # Generate the TSV file, with an extra column ('fasta_file') for the analyzed file name
            # and a specificaton of the taxonomic level ('identification_rank')
            r.level = config.get('level')
            rlist = r.get_values()
            rlist.append(file)
            current_file_handle.write('\t'.join(map(str, rlist)) + '\n')

            # Generate the comment to post
            comment = self.generate_markdown(comment, config, file, r)

        # Post the markdown comment and push the TSV files
        self.gc.post_comment(pr_number, comment)
        current_file_handle.close()
        self.gc.commit_file(current_file_name, f"Add validated FASTA file {current_file_name}")
        self.gc.commit_file(f"{current_file_name}.tsv", f"Add validation results for {current_file_name}")
        self.gc.run_git_command(['git', 'push', 'origin', f"pr-{pr_number}"], f"Failed to push branch pr-{pr_number}")

    @classmethod
    def generate_markdown(cls, comment, config, file, r):
        """
        Generate a markdown comment for a validation result.
        :param comment: The current comment string
        :param config: The Config object
        :param file: The file name
        :param r: A DNAAnalysisResult object
        :return: The updated comment string
        """
        barcode_rank, full_rank, messages = r.calculate_ranks(verbosity=3)
        status_emoji = "✅" if r.passes_all_checks() else "❗"
        obs_taxon_names = "\n".join(f"    - {taxon.name}" for taxon in r.obs_taxon)
        messages_list = messages.split('\n') if isinstance(messages, str) else messages
        message_items = "\n".join(f"    - {m}" for m in messages_list)
        comment += f"""<details>
<summary> {status_emoji} Process ID: {r.process_id} </summary>
    
- File: {file}
- {"✅ No errors" if r.error is None else f"⛔{r.error}"}
- {"✅" if r.check_taxonomy() else "❗"} **Taxonomic check**
  - Expected species as registered at BOLD: {r.species}
  - Expected {config.get('level')} as registered at BOLD: {r.exp_taxon}
  - Observed BLAST hits at {config.get('level')} level: 
{obs_taxon_names}  
- {"✅" if r.check_length() else "❗"} **Sequence length**
  - Net length aligned to marker region: {r.seq_length}  
  - Full sequence length: {r.full_length}
- {"✅" if r.check_seq_quality() else "❗"} **Sequence quality**
  - Ambiguities in marker region: {r.ambiguities}
  - Ambiguities in full sequence: {r.full_ambiguities}
  - Stop codons: {len(r.stop_codons)}
- Rankings:
  - Barcode rank: {barcode_rank}
  - Full rank: {full_rank}
  - Messages: 
{message_items}
    
</details>"""
        return comment


def main(config_file, verbosity):
    """
    Main function for the barcode validator daemon.
    :param config_file: Path to the configuration file
    :param verbosity: Logging verbosity level
    :return: None
    """

    # Initialize the Config object, setup logging
    config = Config()
    config.load_config(config_file)
    config.setup_logging(verbosity)
    logging.info("*** Barcode Validator Daemon starting ***")
    daemon = ValidationDaemon()
    daemon.initialize(config)

    # Start the main loop
    logging.info("Starting main loop")
    while True:

        # Get a list of open PRs and iterate over them
        logging.info("Checking for open PRs...")
        prs = daemon.gc.get_open_prs()
        logging.info(f"Found {len(prs)} open PRs")
        for pr in prs:

            # Process PRs that contain FASTA files
            pr_number = pr['number']
            logging.info(f"Inspecting files from PR {pr['number']}...")
            files = daemon.gc.get_pr_files(pr_number)
            logging.info(f"Found {len(files)} files in PR {pr['number']}")

            # Iterate over the files in the PR and process if any are FASTA files
            if any(f['filename'].endswith(('.fasta', '.fa', '.fas')) for f in files):
                logging.info(f"Processing PR {pr['number']}")
                daemon.process_pr(config, pr_number, pr['head']['ref'])

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
