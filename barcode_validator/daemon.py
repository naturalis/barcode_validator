import argparse
import requests
import time
import os
import csv
import yaml
import sqlite3
import traceback
from datetime import datetime
from typing import Optional
from nbitk.config import Config
from nbitk.logger import get_formatted_logger
from barcode_validator.core import BarcodeValidator
from barcode_validator.github import GitHubClient
from barcode_validator.result import DNAAnalysisResult


GITHUB_TOKEN = os.environ.get('GITHUB_TOKEN')
POLL_INTERVAL = 300  # 5 minutes


class ValidationDaemon:

    def __init__(self):
        self.bv: Optional[BarcodeValidator] = None
        self.gc: Optional[GitHubClient] = None
        self.conn: Optional[sqlite3.Connection] = None
        self.logger = None

    def initialize(self, config: Config):
        class_name = self.__class__.__name__
        self.logger = get_formatted_logger(class_name, config)
        self.bv = BarcodeValidator(config)
        self.bv.initialize()
        self.gc = GitHubClient(config)
        try:
            self.logger.info(f"Going to initialize PR database at {config.get('pr_db_file')}")
            self.conn = self.setup_database(config.get('pr_db_file'))
            self.logger.info("Database initialized")
        except Exception as e:
            self.logger.error(f"Error setting up database: {str(e)}")
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
                self.logger.error(error_msg)
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
        self.logger.info(f"Changed status of PR {pr_number} to 'processing'")
        self.gc.post_comment(pr_number, "\U0001F916 - Hi! This is an automated message from the barcode validation "
                                        "robot. I'm going to validate the FASTA files in your request. This takes " 
                                        "less than two minutes per sequence (about two hours per plate). Subscribe "
                                        "to this thread to get updates. You can safely close this tab.")

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
        self.logger.info(f"Changed status of PR {pr_number} to 'completed'")
        self.gc.post_comment(pr_number, "\U0001F916 - Validation complete. If all looks good, notify @rvosa to merge.")

    def validate_pr(self, config, pr_number, branch):
        """
        Run the validation process for a given pull request.
        :param config: The Config object
        :param pr_number: The pull request number
        :param branch: The branch name
        :return: A dict where keys are FASTA file names and values are lists of DNAAnalysisResult objects
        """
        self.logger.info(f"Starting validation for PR {pr_number}")
        fasta_files = self.gc.fetch_pr_files(branch, pr_number, ['.fasta', '.fa', '.fas'])
        csv_files = self.gc.fetch_pr_files(branch, pr_number, ['.csv'])
        yaml_files = self.gc.fetch_pr_files(branch, pr_number, ['.yaml', '.yml'])

        # Run the validation process for each FASTA file
        all_results = {}
        for file in fasta_files:

            # Validate file, store results
            self.logger.info(f"Validating {file}...")
            results = self.bv.validate_fasta(file, config)
            process_id_to_result = {result.process_id: result for result in results}

            # Join the CSV and YAML file(s) to the results
            self.join_csv_to_result(csv_files, file, process_id_to_result)
            self.join_yaml_to_result(yaml_files, file, results)

            all_results[file] = results
            self.logger.info(f"Validation complete for {file}")

        self.logger.info(f"Validation complete for PR {pr_number}")
        return all_results

    def join_yaml_to_result(self, yaml_files, file, results):
        """
        Join the YAML file to the results.
        :param yaml_files: List of YAML files
        :param file: Processed FASTA file
        :param results: List of DNAAnalysisResult objects
        :return:
        """
        # See if there is a matching Y(A)ML file
        base_name = os.path.splitext(file)[0]
        if f'{base_name}.yaml' in yaml_files:
            matching_file = f'{base_name}.yaml'
        elif f'{base_name}.yml' in yaml_files:
            matching_file = f'{base_name}.yml'
        else:
            matching_file = None

        if matching_file is not None:
            self.logger.info(f"Found YAML file for {file}")

            # Read the YAML file and join it with the results
            with open(matching_file, 'r') as yamlfile:
                yaml_data = yaml.safe_load(yamlfile)
                for result in results:
                    for key, value in yaml_data.items():
                        result.ancillary[key] = value

    def join_csv_to_result(self, csv_files, file, process_id_to_result):
        """
        Join the CSV file to the results.
        :param csv_files: List of CSV files
        :param file: Processed FASTA file
        :param process_id_to_result: Dictionary of process_id to DNAAnalysisResult objects
        :return:
        """
        # See if there is a matching CSV file
        base_name = os.path.splitext(file)[0]
        if f'{base_name}.csv' in csv_files:
            self.logger.info(f"Found CSV file for {file}")

            # Read the CSV file and join it with the results
            with open(f'{base_name}.csv', 'r', newline='') as csvfile:
                reader = csv.DictReader(csvfile)
                for row in reader:
                    process_id = row['Process ID']
                    if process_id in process_id_to_result:
                        result = process_id_to_result[process_id]

                        # Add all fields from CSV to the ancillary dictionary
                        for key, value in row.items():
                            if key != 'Process ID':  # Avoid duplicating the process_id
                                result.ancillary[key] = value

    def post_pr_results(self, config, pr_number, resultset):
        """
        Post a comment to a pull request with the validation results.
        :param config: The Config object
        :param pr_number: The pull request number
        :param resultset: A dict where keys are FASTA file names and values are lists of DNAAnalysisResult objects
        :return: None
        """
        for file, results in resultset.items():

            # Open a new TSV file for each file
            tsv_name = f"{file}.tsv"
            tsv_fh = open(tsv_name, 'w', buffering=1)

            # Configure the header for the TSV file by specifying the taxonomic rank at which we matched obs_taxon
            # and by adding a column that specifies the FASTA file name
            hlist = DNAAnalysisResult.result_fields(config.get('level'))
            hlist.append('fasta_file')
            tsv_fh.write('\t'.join(hlist) + '\n')

            # Write the result objects to the TSV file
            for r in results:
                r.level = config.get('level')  # will be serialized to identification_rank
                rlist = r.get_values()  # will include obs_taxon
                rlist.append(file)  # add the FASTA file name under the 'fasta_file' column
                tsv_fh.write('\t'.join(map(str, rlist)) + '\n')

            # Close the TSV file and commit the files
            tsv_fh.close()
            self.logger.info(f"Going to commit {file} and {tsv_name}")
            self.gc.commit_file(file, f"Validated FASTA file {file} for #PR{pr_number}")
            self.gc.commit_file(tsv_name, f"Results from FASTA file {file} for #PR{pr_number}")

            # Post a comment with the validation results for the file
            self.logger.info(f"Posting comment for {file}")
            comment = f"# Validation Results for {file}\n\n"
            for r in results:
                comment = self.generate_markdown(comment, config, file, r)
            self.gc.post_comment(pr_number, comment)

        # Push the TSV files
        self.logger.info(f"Pushing commits for PR {pr_number}")
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
    config.set('log_level', verbosity)
    logger = get_formatted_logger(__name__, config)
    logger.info("*** Barcode Validator Daemon starting ***")
    daemon = ValidationDaemon()
    daemon.initialize(config)

    # Start the main loop
    logger.info("Starting main loop")
    while True:

        try:
            # Get a list of open PRs and iterate over them
            logger.info("Checking for open PRs...")
            prs = daemon.gc.get_open_prs()
            logger.info(f"Found {len(prs)} open PRs")
            for pr in prs:

                # Process PRs that contain FASTA files
                pr_number = pr['number']
                logger.info(f"Inspecting files from PR {pr['number']}...")
                files = daemon.gc.get_pr_files(pr_number)
                logger.info(f"Found {len(files)} files in PR {pr['number']}")

                # Iterate over the files in the PR and process if any are FASTA files
                if any(f['filename'].endswith(('.fasta', '.fa', '.fas')) for f in files):
                    logger.info(f"Processing PR {pr['number']}")
                    daemon.process_pr(config, pr_number, pr['head']['ref'])

        except Exception as e:
            logger.error(f"Error in main loop: {str(e)}")
            logger.error(traceback.format_exc())
            logger.info("Hopefully this was a transient error. Trying again in 5 minutes...")

        time.sleep(POLL_INTERVAL)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="GitHub PR Validation Daemon for barcode_validator")
    parser.add_argument("-c", "--config_file", required=True, help="Path to the configuration YAML file")
    parser.add_argument("-v", "--verbosity", choices=['DEBUG', 'INFO', 'WARNING', 'ERROR', 'CRITICAL'],
                        default='INFO', help="Set the logging verbosity (default: INFO)")
    args = parser.parse_args()

    main(args.config_file, args.verbosity)
