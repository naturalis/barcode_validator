import logging
import subprocess
import requests
import os


class GitHubClient:
    def __init__(self, repo_owner, repo_name, github_token, clone_path):
        self.repo_owner = repo_owner
        self.repo_name = repo_name
        self.github_token = github_token
        self.clone_path = clone_path
        self.base_url = f"https://api.github.com/repos/{repo_owner}/{repo_name}"
        self.headers = {
            "Authorization": f"token {github_token}",
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
            logging.error(f"{error_message}: {result.stderr}")
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
            logging.debug(f"Changing working directory from {current_dir} to {self.clone_path}")
            os.chdir(self.clone_path)