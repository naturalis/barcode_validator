import subprocess
import requests
import os
from nbitk.config import Config
from nbitk.logger import get_formatted_logger


class GitHubClient:
    def __init__(self, config: Config):
        self.logger = get_formatted_logger(__name__, config)
        self.repo_owner = config.get('repo_owner')
        self.repo_name = config.get('repo_name')
        self.github_token = os.environ.get('GITHUB_TOKEN')
        self.clone_path = config.get('repo_location')
        self.base_url = f"https://api.github.com/repos/{self.repo_owner}/{self.repo_name}"
        self.headers = {
            "Authorization": f"token {self.github_token}",
            "Accept": "application/vnd.github.v3+json"
        }

    def get_open_prs(self):
        """
        Get a list of open pull requests.
        :return: A list of open pull requests
        """
        url = f"{self.base_url}/pulls"
        params = {"state": "open"}
        response = requests.get(url, headers=self.headers, params=params)
        response.raise_for_status()
        return response.json()

    def get_pr_files(self, pr_number):
        """
        Get a list of files in a pull request.
        :param pr_number: The pull request number
        :return: A list of files in the pull request
        """
        url = f"{self.base_url}/pulls/{pr_number}/files"
        response = requests.get(url, headers=self.headers)
        response.raise_for_status()
        return response.json()

    def fetch_pr_files(self, branch, pr_number, extensions):
        """
        Fetch the FASTA files from a pull request.
        :param branch: The branch name
        :param pr_number: The pull request number
        :param extensions: A list of file extensions to filter by, e.g. ['.fasta', '.fa', '.fas']
        :return: A list of FASTA files
        """

        # Fetch the latest changes
        self.logger.info(f"Fetching latest changes for PR {pr_number}")
        self.run_git_command(['git', 'fetch', 'origin'], "Failed to fetch from origin")

        # Create or reset the PR branch
        pr_branch = f"pr-{pr_number}"
        self.logger.info(f"Creating/resetting branch {pr_branch}")
        self.run_git_command(['git', 'checkout', '-B', pr_branch, f'origin/{branch}'],
                                f"Failed to create/reset branch {pr_branch}")

        # Request the JSON info about the files from the PR matching the extension
        self.logger.info(f"Getting files for PR {pr_number}")
        files = self.get_pr_files(pr_number)
        matching_files = [f for f in files if f['filename'].endswith(extensions)]
        self.logger.info(f"Found {len(matching_files)} files with {extensions} in PR {pr_number}")

        # Download the files
        result = []
        for file in matching_files:

            # Fetch file content
            file_url = file['raw_url']
            self.logger.info(f"Fetching file from {file_url}")
            response = requests.get(file_url, headers=self.headers)
            if response.status_code == 200:

                # Create directory if it doesn't exist
                os.makedirs(os.path.dirname(file['filename']), exist_ok=True)

                # Save file content
                self.logger.info(f"Saving file content to {file['filename']}")
                with open(file['filename'], 'wb') as f:
                    f.write(response.content)
                    result.append(file['filename'])

        return result

    def post_comment(self, pr_number, comment):
        """
        Post a comment on a pull request.
        :param pr_number: The pull request number
        :param comment: The comment to post
        :return: The posted comment
        """
        url = f"{self.base_url}/issues/{pr_number}/comments"
        data = {"body": comment}
        response = requests.post(url, headers=self.headers, json=data)
        self.logger.info(response)
        response.raise_for_status()
        return response.json()

    def run_git_command(self, command, error_message):
        """
        Run a git command and log any errors.
        :param command: The command to run
        :param error_message: The error message to log
        :return: The output of the command
        """
        self.ensure_correct_directory()
        result = subprocess.run(command, capture_output=True, text=True)
        if result.returncode != 0:
            self.logger.error(f"{error_message}: {result.stderr}")
            raise RuntimeError(f"Git command failed: {' '.join(command)}")
        return result.stdout

    def commit_file(self, filename, message):
        """
        Commit a file to the local git repository.
        :param filename: The name of the file to commit
        :param message: The commit message
        :return: None
        """
        self.run_git_command(['git', 'add', filename], f"Failed to add {filename}")
        self.run_git_command(['git', 'commit', '-m', message], f"Failed to commit {filename}")

    def ensure_correct_directory(self):
        """
        Ensure that the current working directory is the clone path.
        :return: None
        """
        current_dir = os.getcwd()
        if current_dir != self.clone_path:
            self.logger.debug(f"Changing working directory from {current_dir} to {self.clone_path}")
            os.chdir(self.clone_path)