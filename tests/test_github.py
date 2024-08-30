import pytest
from unittest.mock import Mock, patch
import os
import subprocess
from barcode_validator.github import GitHubClient


@pytest.fixture
def github_client():
    return GitHubClient("test_owner", "test_repo", "test_token", "/test/clone/path")


@patch('barcode_validator.github.requests.get')
def test_get_open_prs(mock_get, github_client):
    mock_response = Mock()
    mock_response.json.return_value = [{"number": 1, "title": "Test PR"}]
    mock_get.return_value = mock_response

    prs = github_client.get_open_prs()

    assert prs == [{"number": 1, "title": "Test PR"}]
    mock_get.assert_called_once_with(
        "https://api.github.com/repos/test_owner/test_repo/pulls",
        headers=github_client.headers,
        params={"state": "open"}
    )


@patch('barcode_validator.github.requests.get')
def test_get_pr_files(mock_get, github_client):
    mock_response = Mock()
    mock_response.json.return_value = [{"filename": "test.py", "status": "modified"}]
    mock_get.return_value = mock_response

    files = github_client.get_pr_files(1)

    assert files == [{"filename": "test.py", "status": "modified"}]
    mock_get.assert_called_once_with(
        "https://api.github.com/repos/test_owner/test_repo/pulls/1/files",
        headers=github_client.headers
    )


@patch('barcode_validator.github.requests.post')
def test_post_comment(mock_post, github_client):
    mock_response = Mock()
    mock_response.json.return_value = {"id": 1, "body": "Test comment"}
    mock_post.return_value = mock_response

    comment = github_client.post_comment(1, "Test comment")

    assert comment == {"id": 1, "body": "Test comment"}
    mock_post.assert_called_once_with(
        "https://api.github.com/repos/test_owner/test_repo/issues/1/comments",
        headers=github_client.headers,
        json={"body": "Test comment"}
    )


@patch('barcode_validator.github.subprocess.run')
@patch('barcode_validator.github.os.getcwd')
@patch('barcode_validator.github.os.chdir')
def test_run_git_command(mock_chdir, mock_getcwd, mock_run, github_client):
    mock_getcwd.return_value = "/test/clone/path"
    mock_result = Mock()
    mock_result.returncode = 0
    mock_result.stdout = "Test output"
    mock_run.return_value = mock_result

    output = github_client.run_git_command(["git", "status"], "Failed to get status")

    assert output == "Test output"
    mock_run.assert_called_once_with(["git", "status"], capture_output=True, text=True)
    mock_getcwd.assert_called_once()
    mock_chdir.assert_not_called()


@patch('barcode_validator.github.GitHubClient.run_git_command')
def test_commit_file(mock_run_git_command, github_client):
    github_client.commit_file("test.py", "Test commit")

    mock_run_git_command.assert_any_call(['git', 'add', 'test.py'], "Failed to add test.py")
    mock_run_git_command.assert_any_call(['git', 'commit', '-m', 'Test commit'], "Failed to commit test.py")


@patch('barcode_validator.github.os.getcwd')
@patch('barcode_validator.github.os.chdir')
def test_ensure_correct_directory(mock_chdir, mock_getcwd, github_client):
    mock_getcwd.return_value = "/wrong/path"

    github_client.ensure_correct_directory()

    mock_getcwd.assert_called_once()
    mock_chdir.assert_called_once_with("/test/clone/path")


@patch('barcode_validator.github.requests.get')
def test_get_open_prs_error(mock_get, github_client):
    mock_response = Mock()
    mock_response.raise_for_status.side_effect = Exception("API Error")
    mock_get.return_value = mock_response

    with pytest.raises(Exception):
        github_client.get_open_prs()


if __name__ == "__main__":
    pytest.main()
