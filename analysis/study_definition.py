from cohortextractor import (
    StudyDefinition, 
    patients, 
    filter_codes_by_category
)

# Import codelists.py script
from codelists import *

import pandas as pd

from functions import *

# import the variables for deriving JCVI groups
from elig_definition import (
    elig_variables, 
    study_parameters
)

# set seed so that dummy data can be reproduced
import numpy as np
np.random.seed(study_parameters["seed"])

# regions
regions = pd.read_csv(
    filepath_or_buffer='./analysis/lib/regions.csv',
    dtype=str
)
ratio_regions = { regions['region'][i] : float(regions['ratio'][i]) for i in regions.index }

study=StudyDefinition(

    default_expectations = {
        "date": {"earliest": "2020-01-01", "latest": "2021-12-31"},
        "rate": "uniform",
        "incidence": 0.9,
    },

    ##########################
    ### ELIGIBILITY GROUPS ###
    ##########################
    # population = patients.all(),
    **elig_variables, 

    # ##################
    # ### POPULATION ###
    # ##################

    # extract all for consort diagram
    population = patients.all(),

    has_follow_up = patients.registered_with_one_practice_between(
            start_date = "2020-01-01", # start of 2020 
            end_date = "elig_date - 1 day", 
        ),
    
    ######################
    ### COVID VACCINES ###
    ######################

    # first 2 covid vaccine dates
    **vaccination_date_X(
        name = "covid_vax_disease",
        index_date = "1900-01-01", # set implausibly early to catch errors
        n = 2,
        target_disease_matches="SARS-2 CORONAVIRUS"
    ),

    #############################
    ### DEMOGRAPHIC VARIABLES ###
    #############################

    # STATIC
    # patient sex
    sex = patients.sex(
        return_expectations = {
        "rate": "universal",
        "category": {"ratios": {"M": 0.49, "F": 0.51}},
        "incidence": 0.99,
        }
    ),

    # 1 DEC 2021
    # region - NHS England 9 regions
    region = patients.registered_practice_as_of(
        "elig_date - 1 day",
        returning = "nuts1_region_name",
        return_expectations = {
            "rate": "universal",
            "category": {
                "ratios": ratio_regions,
            },
            "incidence": 0.99
        },
    ),

    # ethnicity (6 categories)
    ethnicity = patients.categorised_as(
        {
        "Unknown": "DEFAULT",
        "White": "eth6='1'",
        "Mixed": "eth6='2'",
        "Asian or Asian British": "eth6='3'",
        "Black or Black British": "eth6='4'",
        "Other": "eth6='5'",
        },
        eth6 = patients.with_these_clinical_events(
            ethnicity_codes_6,
            returning = "category",
            on_or_before = "elig_date - 1 day",
            find_last_match_in_period = True,
            include_date_of_match = False,
            return_expectations = {
                "incidence": 0.75,
                "category": {
                "ratios": { "1": 0.30, "2": 0.20, "3": 0.20, "4": 0.20, "5": 0.05, "6": 0.05, },
                },
            },
        ),
        return_expectations = {
            "rate": "universal",
            "category": {
            "ratios": {
                "White": 0.30,
                "Mixed": 0.20,
                "Asian or Asian British": 0.20,
                "Black or Black British": 0.20,
                "Other": 0.05,
                "Unknown": 0.05,
                },
            },
        },
    ),

    # any death
    death_date=patients.died_from_any_cause(
        returning="date_of_death",
        date_format="YYYY-MM-DD",
        return_expectations = { "incidence": 0.01},
    ),

    dereg_date=patients.date_deregistered_from_all_supported_practices(
        on_or_after="elig_date",
        date_format="YYYY-MM-DD",
        return_expectations = {"incidence": 0.01},
    ),

    ## IMD - quintile
    imd_Q5=patients.categorised_as(
        {
          "Unknown": "DEFAULT",
          "1 (most deprived)": "imd >= 0 AND imd < 32844*1/5",
          "2": "imd >= 32844*1/5 AND imd < 32844*2/5",
          "3": "imd >= 32844*2/5 AND imd < 32844*3/5",
          "4": "imd >= 32844*3/5 AND imd < 32844*4/5",
          "5 (least deprived)": "imd >= 32844*4/5 AND imd <= 32844",
        },
        return_expectations={
          "rate": "universal",
          "category": {"ratios": {"Unknown": 0.02, "1 (most deprived)": 0.18, "2": 0.2, "3": 0.2, "4": 0.2, "5 (least deprived)": 0.2}},
        },
        imd=patients.address_as_of(
          "elig_date - 1 day",
          returning="index_of_multiple_deprivation",
          round_to_nearest=100,
          return_expectations={
          "category": {"ratios": {c: 1/320 for c in range(100, 32100, 100)}}
          }
        ),
    ),


)
