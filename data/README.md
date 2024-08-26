# Crowdsourced validation workflow

## 1. Forking

When you upload files for validation, you must do so on a fork of this repository. So, while you are uploading, you
cannot be doing this on `naturalis/barcode_validator`, but on `yourusername/barcode_validator`. If you haven't done
so already, you can create a fork by visiting [this](https://github.com/naturalis/barcode_validator/fork) page.

It is advisable to keep your fork up-to-date with the upstream repository. If you are not actively developing on your
fork, you can sync your fork by clicking the "Sync fork" button on the fork's page on GitHub. To the extent that you
are doing things on your own fork (presumably because you are experimenting with file uploads), it might keep things
simpler for you to discard any needless changes (also an option under the "Sync fork" button).

## 2. Uploading

You can upload files on github by navigating to the `data` directory in your fork and clicking the "Add File > Upload 
files" button. You can upload files one by one, or multiple files at once. The files you upload must be in the `fasta` 
format, i.e. they must have the `.fa`, `.fas` or `.fasta` extension. If you upload files with a different extension, 
they will be ignored.

## 3. Validation

Once you have uploaded your files, you can initiate the validation process by creating a pull request. You can do this
by navigating to the `Pull requests` tab in your fork and clicking the "New pull request" button. You will be presented
with a page that shows the changes you are proposing to make. You can click the "Create pull request" button to create
the pull request. Once you have created the pull request, the validation process will start automatically. Please note:

- Each **sequence** takes about two minutes to validate, so large files will take a while to complete
- If the bot is busy validating, it will not queue new validation jobs to run concurrently
- Closing the pull request will not stop the validation process
- The validation bot checks for new pull requests every 5 minutes
- So, please have patience: this is not a realtime thing
- You can configure your github settings to watch the PR for updates in a variety of ways (e.g. emails etc.)

## 4. Results

The results of the validation process will be posted as comments on the pull request. In addition, the results will be
written to TSV files, which are pushed to a separate branch of the naturalis repository. The branch is named after the
pull request number: `pr-<number>`, which you can discover by looking at the URL of the pull request or by checking the
[branches](https://github.com/naturalis/barcode_validator/branches) page of the naturalis repository.

## 5. Next steps

For now, we are not filtering the results in any way. We keep everything so that we build up a data set that we can use
to compare assembly pipelines and their paramterizations. In the future, we may want to filter the results to only keep
those of a certain quality, merge those into a `clean` branch, and pass those on to CBG for batch submissions.

We may also want to encapsulate the validation files (fasta and TSV) inside RO-crate archives and push those to Zenodo
for long-term storage of fuller provenance records FAIR-ly identifiable by Zenodo-minted DOIs.

