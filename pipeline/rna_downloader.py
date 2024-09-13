# -*- coding: utf-8 -*-
import os
from selenium import webdriver
from selenium.webdriver.common.by import By
from selenium.webdriver.support.ui import WebDriverWait
from selenium.webdriver.support import expected_conditions as EC
from selenium.webdriver.chrome.service import Service
from webdriver_manager.chrome import ChromeDriverManager
from selenium.webdriver.chrome.options import Options
import requests
import time
import subprocess
import numpy as np
from concurrent.futures import ThreadPoolExecutor
from seq_encoder import *
from mfold import ssDNA_secondary_structure_prediction_mfold
from rna2ssdna import process_rna_to_dna

def wait_and_get_download_link(driver):
    try:
        time.sleep(5)
        # Wait for the link containing "GetResult" text
        download_link = WebDriverWait(driver, 60).until(
            EC.presence_of_element_located((By.XPATH, '//a[contains(@href, "GetResult")]'))
        )
        # Return the href attribute of the link
        return download_link.get_attribute('href')
    except Exception as e:
        driver.save_screenshot("web_screenshot/page_snapshot_exception.png")
        print("Error occurred while waiting for download link:", e)
        return None
    finally:
        driver.quit()

def get_rna_composition_download_link(submit_content):
    """
    Processes an RNA sequence through RNAComposer, and returns the download link.

    :param submit_content: RNA sequence and structure to process
    :return: Download link for the resulting PDB file
    """
    download_link = None
    driver = None
    try_attempts = 5  # Number of attempts to find the download link
    for attempt in range(try_attempts):
        print("="*10 + f"Attempt {attempt + 1} for RNAComposer" + "="*10)
        try:
            chrome_options = Options()
            chrome_options.add_argument("--headless")
            chrome_options.add_argument('--ignore-certificate-errors')

            # Initialize WebDriver instance with options
            chrome_driver_path = r"C:\Users\j9366\.wdm\drivers\chromedriver\win64\127.0.6533.88\chromedriver-win32\chromedriver.exe"
            service = Service(chrome_driver_path)
            driver = webdriver.Chrome(service=service, options=chrome_options)
            
            # Open RNAComposer
            driver.get("https://rnacomposer.cs.put.poznan.pl/")

            # Wait for the text area to be clickable
            wait = WebDriverWait(driver, 10)
            textarea = wait.until(EC.element_to_be_clickable((By.ID, "input")))
            print("="*10 + "Submitting sequence to RNAComposer" + "="*10)

            # Clear the textarea and enter the RNA sequence
            textarea.clear()
            textarea.send_keys(submit_content)

            # Submit the sequence for processing
            compose_button = driver.find_element(By.ID, "send")
            compose_button.click()

            download_link = wait_and_get_download_link(driver)
            if download_link:
                break
            else:
                print(f"Attempt {attempt + 1}: Download link not found, trying again..." if attempt < try_attempts - 1 else "All attempts failed: Download link not found.")
                time.sleep(30)
        except Exception as e:
            print("Error:", e)
        finally:
            if driver:
                driver.quit()
    return download_link

def download_file(download_url, save_path="/path/to/save/file"):
    """
    Downloads a file from the given URL and saves it to the specified path.

    :param download_url: URL of the file to download
    :param save_path: Path where the file should be saved
    :return: Path to the saved file
    """
    try:
        response = requests.get(download_url)
        response.raise_for_status()  # Raise an error for bad responses

        # Make sure dir exists
        os.makedirs(os.path.dirname(save_path), exist_ok=True)

        with open(save_path, 'wb') as file:
            file.write(response.content)
        return save_path
    except requests.RequestException as e:
        print(f"Error downloading the file: {e}")
        return None

def get_rna_composition_and_download(input_sequence, is_dna=False):
    if '\n' in input_sequence:
        sequence, secondary_structure = input_sequence.split('\n')
        submit_sequence = f">rna\n{sequence}\n{secondary_structure}"
    else:
        sequence = input_sequence
        submit_sequence = f">rna\n{sequence}"
    
    print("Submit:")
    print(submit_sequence)
    print("="*20)
    
    download_link = get_rna_composition_download_link(submit_sequence)

    if download_link:
        prefix = "d_" if is_dna else ""
        save_path = f"./Aptamer_pdb/{prefix}{encode_rna(sequence)}.pdb"
        print(f"File will be saved to: {save_path}")
        file_path = download_file(download_link, save_path)
        if file_path:
            print(f"File successfully downloaded and saved to {file_path}")
            return file_path
        else:
            print("Failed to download the file.")
    else:
        print("Failed to get the download link.")
    return None

def get_dna_composition_and_download(dna_sequence):
    # Convert DNA sequence to RNA sequence
    rna_sequence = dna_sequence.replace('T', 'U')
    
    # Predict secondary structure
    secondary_structure = ssDNA_secondary_structure_prediction_mfold(dna_sequence)
    
    if secondary_structure is None:
        print("Unable to predict secondary structure. Cannot continue processing.")
        return None
    
    # Prepare the sequence for RNAComposer
    rna_composer_input = f"{rna_sequence}\n{secondary_structure}"

    rna_file_path = get_rna_composition_and_download(rna_composer_input, is_dna=True)
    
    if rna_file_path:
        # Use process_rna_to_dna function to convert the file
        try:
            dna_file_path = process_rna_to_dna(rna_file_path)
            print(f"Converted ssDNA file saved as: {dna_file_path}")
            return dna_file_path
        except Exception as e:
            print(f"Error occurred while converting RNA to ssDNA: {e}")
            return rna_file_path  # If conversion fails, return the original RNA file path
    else:
        print("Failed to process DNA sequence.")
        return None

def test_rna_composition_and_download():
    # Test case with sequence and structure
    rna_with_structure = "AGCAGCACAGAGGUCAGAUGAUGGAGGUUGGUCGGGUGGGCAAUCAUUCUCCUAUGCGUGCUACCGUGAA\n((((...((.(((....(((((....................)))))...))).))..))))........"
    result = get_rna_composition_and_download(rna_with_structure)

def test_dna_composition_and_download():
    # Test DNA sequence
    dna_sequence = "AGCAGCACAGAGGTCAGATGATGGAGGTTGGTCGGGTGGGCAATCATTCTCCTATGCGTGCTACCGTGAA"
    print("Testing DNA sequence:")
    print(dna_sequence)
    print("="*20)
    
    result = get_dna_composition_and_download(dna_sequence)
    
    if result:
        print(f"DNA processing successful. File saved as: {result}")
    else:
        print("DNA processing failed.")

def test_rna2ssdna_conversion():
    input_file = "./Aptamer_pdb/dna_ROGERRFODIWHRFIFSHIHOWOYPATISIPESIQ.pdb"
    if os.path.exists(input_file):
        try:
            output_file = process_rna_to_dna(input_file)
            if output_file:
                print(f"Conversion successful. Output file: {output_file}")
            else:
                print("Conversion failed. No output file produced.")
        except Exception as e:
            print(f"Conversion failed with error: {e}")
    else:
        print(f"Input file not found: {input_file}")

if __name__ == "__main__":
    # Uncomment the line below to run the RNA test
    # test_rna_composition_and_download()
    
    # Run the DNA test
    # test_dna_composition_and_download()

    # Test RNA to ssDNA conversion
    test_rna2ssdna_conversion()