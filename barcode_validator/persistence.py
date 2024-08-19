import logging
# import gspread
import os
import datetime
from barcode_validator.config import Config
# from oauth2client.service_account import ServiceAccountCredentials


def persist_gspread(result):

    # Get the spreadsheet URL and credentials file path from the config
    config = Config()
    spreadsheet_url = config.get('gspread')
    credentials_file = config.get('credentials_file')
    if not spreadsheet_url:
        raise ValueError("Google Spreadsheet URL not found in config")
    if not credentials_file:
        raise ValueError("Credentials file path not found in config")

    # Get the username from the environment variable
    username = os.environ.get('USER', 'Unknown')

    # Use credentials to create a client to interact with the Google Drive API
    scope = ['https://spreadsheets.google.com/feeds', 'https://www.googleapis.com/auth/drive']
    try:
        creds = ServiceAccountCredentials.from_json_keyfile_name(credentials_file, scope)
        client = gspread.authorize(creds)
    except Exception as e:
        logging.error(f"Failed to authenticate: {str(e)}")
        raise

    # Open the spreadsheet
    try:
        sheet = client.open_by_url(spreadsheet_url).sheet1
    except gspread.exceptions.SpreadsheetNotFound:
        logging.error(f"Spreadsheet with URL {spreadsheet_url} not found")
        raise
    except Exception as e:
        logging.error(f"Failed to open spreadsheet: {str(e)}")
        raise

    # Try to find the row with the matching process_id
    try:
        cell = sheet.find(result.process_id)
        row = cell.row
        action = 'updated'
    except gspread.exceptions.CellNotFound:
        # If not found, append a new row
        row = len(sheet.get_all_values()) + 1
        action = 'added'

    # Prepare the data to be inserted/updated
    data = [
        result.process_id,
        result.seq_length,
        result.obs_taxon,
        result.exp_taxon,
        result.species,
        ','.join(map(str, result.stop_codons)),
        result.ambiguities,
        result.passes_all_checks(),
        datetime.now().strftime('%Y-%m-%d %H:%M:%S'),
        username
    ]

    # Update the row
    try:
        sheet.update(f'A{row}:J{row}', [data])
        logging.info(
            f"Result for process_id {result.process_id} has been {action} in the spreadsheet by user {username}")
    except Exception as e:
        logging.error(f"Failed to update spreadsheet: {str(e)}")
        raise


def persist_result(result):
    logging.info(f"Updating results")
    config = Config()
    persistence = config.get('persist')

    # Persist the results to a Google spreadsheet
    if persistence == 'gspread':
        return persist_gspread(result)

    pass


