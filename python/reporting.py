#!/usr/bin/env python
"""Functions related to reporting."""

# Imports
from pathlib import Path
from time import sleep

from munch import Munch
from tabulate import tabulate
from jinja2 import Template
import pandas as pd

from selenium import webdriver
from selenium.webdriver.common.keys import Keys

# Metadata
__author__ = "Gus Dunn"
__email__ = "w.gus.dunn@gmail.com"


# Functions
def df_to_markdown(df):
    """Return markdown representation of a dataframe's table."""
    return tabulate(df, headers=list(df.columns.values), tablefmt='pipe')


def load_template(source):
    """Return jinja2 Template object instantiated with ``source``."""
    try:
        path = Path(source)
        exists = path.exists()

        return Template(source=path.read_text())

    except OSError:
        return Template(source=source)


def get_abs_path(path):
    """Resolve and return full path for `path`."""
    p = Path(path)
    return str(p.resolve())

def generate_metadata_table(**kwargs):
    """Convert key-value pairs into a dataframe and return."""
    return pd.DataFrame(kwargs, index=['']).T


def get_gene_info(gnames):
    """Collect, organize, and return information for each gene for use in generating report."""
    info = Munch()

    for name in gnames:
        info.name = name
        info.region

# Getting genome region images
def init_browser(prefs=None, file_paths=None):
    """Return selenium Chrome browser object instatiated with contents of ``prefs`` and loaded with correct genome browser track options and data."""
    if prefs is None:
        prefs = {}

    if file_paths is None:
        file_paths = []

    try:
        prefs["download.default_directory"]
    except KeyError:
        prefs["download.default_directory"] = get_abs_path('.')

    opts = webdriver.ChromeOptions()
    opts.add_experimental_option("prefs", prefs)

    browser = webdriver.Chrome(chrome_options=opts)
    browser.get('https://www.vectorbase.org/Glossina_fuscipes/Location/View?r=Scaffold1:1-10000;share_ref=conf-s1953407-05bddb91d84e07835baacc651c327c88')

    for path in file_paths:
        add_data_file_to_vb(browser, path)


    return browser

def add_data_file_to_vb(browser, path):
    """Add a single data file to the genome browser session at VB."""
    path = Path(path)

    browser.get('https://www.vectorbase.org/Glossina_fuscipes/UserData/SelectFile')

    f = browser.find_element_by_name('file')
    f.send_keys(str(path.resolve()))
    submit = browser.find_element_by_name('submit_button')
    submit.click()

def enable_all_personal_data_tracks(browser):
    """Activate all of the tracks we added data for."""
    browser.get('https://www.vectorbase.org/Glossina_fuscipes/Location/View?r=Scaffold1:1-10000')

    browser.find_element_by_link_text('Configure this page').click()
    sleep(2)
    browser.find_element_by_xpath('//*[@id="modal_config_viewbottom"]/div[1]/ul/li[5]/a').click() # "Your Data"
    browser.find_element_by_xpath('//*[@id="location_viewbottom_configuration"]/div[1]/div/div/strong').click() # "Enable/disable all tracks"
    browser.find_element_by_xpath('//*[@id="location_viewbottom_configuration"]/div[1]/div/div/ul/li[3]').click() # "Normal"
    browser.find_element_by_xpath('//*[@id="modal_panel"]/div[1]/div[2]').click() # "Save and close"



def get_genome_image(browser, coords, name=None, url_tmpl=None):
    """Save genome region image for coords provided.

    coords should follow this structure:
        {'chr':'Scaffold104','start':593216,'end':597969}
    """
    if name is None:
        name = '{chr}_{start}_{end}.png'.format(**coords)

    if url_tmpl is None:
        url_tmpl = 'https://www.vectorbase.org/Glossina_fuscipes/ImageExport/ImageFormats?component=ViewBottom;data_type=Location;db=core;r={chr}:{start}-{end}'

    url = url_tmpl.format(**coords)

    browser.get(url)

    filename = browser.find_element_by_name('filename')
    filename.clear()
    filename.send_keys(name, Keys.RETURN)

def get_link_to_vb_prot_summary(gname):
    """Return link to protein summary page on vectorbase."""
    return "https://www.vectorbase.org/Glossina_fuscipes/Transcript/ProteinSummary?db=core;g={gname}".format(gname=gname)

def get_vb_prot_domain_table(gname):
    """Return markdown table-ified version of vectorbase's "Download whole table" as csv."""
    link = "https://www.vectorbase.org/Glossina_fuscipes/Transcript/Domains?db=core;g={gname}".format(gname=gname)

def get_orthodb_table(gname):
    """Return markdown table-ified version of orthoDB v9 query for `gname` vs metazoa."""
    link = "http://www.orthodb.org/v9/tab?query={gname}&level=33208&species=33208".format(gname=gname)

    df = pd.read_csv(link, sep='\t')

    if not df.empty:
        return df_to_markdown(df)
    else:
        return None
