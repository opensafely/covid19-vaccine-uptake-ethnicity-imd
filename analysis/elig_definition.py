from datetime import date

from cohortextractor import (
    patients, 
)

# import codelists.py script
from codelists import *

from functions import *

# import json module
import json

import pandas as pd

### import groups and dates
# atrisk_group
atrisk_group = pd.read_csv(
    filepath_or_buffer='./analysis/lib/atrisk_group.csv',
    dtype=str
)
str_atrisk = atrisk_group['atrisk_group'][0]

# jcvi_groups
jcvi_groups = pd.read_csv(
    filepath_or_buffer='./analysis/lib/jcvi_groups.csv',
    dtype=str
)
dict_jcvi = { jcvi_groups['group'][i] : jcvi_groups['definition'][i] for i in jcvi_groups.index }
ratio_jcvi = { jcvi_groups['group'][i] : 1/len(jcvi_groups.index) for i in jcvi_groups.index }

# elig_dates
elig_dates = pd.read_csv(
    filepath_or_buffer='./analysis/lib/elig_dates.csv',
    dtype=str
)
dict_elig = { elig_dates['date'][i] : elig_dates['description'][i] for i in elig_dates.index }
ratio_elig = { elig_dates['date'][i] : 1/len(elig_dates.index) for i in elig_dates.index }

#study_parameters
with open("./analysis/lib/study_parameters.json") as f:
  study_parameters = json.load(f)

# set seed so that dummy data can be reproduced
import numpy as np
seed=int(study_parameters["seed"])
np.random.seed(seed)

# define variables explicitly
ref_age_1 = study_parameters["ref_age_1"] # reference date for calculating age for phase 1 groups
ref_age_2 = study_parameters["ref_age_2"] # reference date for calculating age for phase 2 groups
ref_cev = study_parameters["ref_cev"] # reference date for calculating clinically extremely vulnerable group
ref_ar = study_parameters["ref_ar"] #reference date for caluclating at risk group
start_date = study_parameters["start_date"] # start of phase 1
pandemic_start = study_parameters["pandemic_start"]

# Notes:
# for inequalities in the study definition, an extra expression is added to align with the comparison definitions in https://github.com/opensafely/covid19-vaccine-coverage-tpp-emis/blob/master/analysis/comparisons.py
# variables that define JCVI group membership MUST NOT be dependent on elig_date (index_date), this is for selecting the population based on registration dates and for deriving descriptive covariates
# JCVI groups are derived using ref_age_1, ref_age_2, ref_cev and ref_ar

elig_variables = dict(
  
  # age on phase 1 reference date
    age_1 = patients.age_as_of(
        ref_age_1,
        return_expectations = {
            "int": {"distribution": "population_ages"},
            "rate": "universal",
        },
    ),

    # age on phase 2 reference date
    age_2 = patients.age_as_of(
        ref_age_2,
        return_expectations = {
            "int": {"distribution": "population_ages"},
            "rate": "universal",
        },
    ),

    # clinically extremely vulnerable group
    cev_group = patients.satisfying(
        "severely_clinically_vulnerable AND NOT less_vulnerable",

        # shielding - first flag all patients with "high risk" codes
        severely_clinically_vulnerable = patients.with_these_clinical_events(
            shield_primis,
            returning = "binary_flag",
            on_or_before = days(ref_cev, -1),
            find_last_match_in_period = True,
        ),

        # find date at which the high risk code was added
        severely_clinically_vulnerable_date = patients.date_of(
            "severely_clinically_vulnerable",
            date_format = "YYYY-MM-DD",
        ),

        # not shielding (medium and low risk) - only flag if later than 'shielded'
        less_vulnerable = patients.with_these_clinical_events(
            nonshield_primis,
            between = ["severely_clinically_vulnerable_date + 1 day", days(ref_cev, -1)],
        ),

        return_expectations={"incidence": 0.01},

    ),

    ## at-risk group variables
    # asthma
    asthma_group=patients.satisfying(
        """
        astadm OR
        (ast AND astrxm1 AND astrxm2 AND astrxm3)
        """,
        # asthma admission codes in past 24 months
        astadm = patients.with_these_clinical_events(
            astadm_primis,
            returning = "binary_flag",
            between = [days(ref_ar, -2*365), days(ref_ar, -1)],
        ),
        # asthma diagnosis code
        ast = patients.with_these_clinical_events(
            ast_primis,
            returning="binary_flag",
            on_or_before = days(ref_ar, -1),
        ),
        # asthma systemic steroid prescription code in month 1
        astrxm1 = patients.with_these_medications(
            astrx_primis,
            returning = "binary_flag",
            between = [days(ref_ar, -31), days(ref_ar, -1)],
        ),
        # asthma systemic steroid prescription code in month 2
        astrxm2 = patients.with_these_medications(
            astrx_primis,
            returning = "binary_flag",
            between = [days(ref_ar, -61), days(ref_ar, -32)],
        ),
        # asthma systemic steroid prescription code in month 3
        astrxm3 = patients.with_these_medications(
            astrx_primis,
            returning = "binary_flag",
            between = [days(ref_ar, -91), days(ref_ar, -62)],
        ),
        return_expectations={"incidence": 0.1},
    ),

    # chronic respiratory disease other than asthma
    resp_group = patients.with_these_clinical_events(
        resp_primis,
        returning = "binary_flag",
        on_or_before = days(ref_ar, -1),
        return_expectations = {"incidence": 0.02},
    ),

    # chronic neurological disease including significant learning disorder
    cns_group = patients.with_these_clinical_events(
        cns_primis,
        returning = "binary_flag",
        on_or_before = days(ref_ar, -1),
        return_expectations = {"incidence": 0.01},
    ),

    # diabetes
    diab_group = patients.satisfying(
        """
        (NOT dmres_date AND diab_date) OR
        (dmres_date < diab_date)
        """,
        diab_date = patients.with_these_clinical_events(
            diab_primis,
            returning = "date",
            find_last_match_in_period = True,
            on_or_before = days(ref_ar, -1),
            date_format = "YYYY-MM-DD",
        ),
        dmres_date = patients.with_these_clinical_events(
            dmres_primis,
            returning = "date",
            find_last_match_in_period = True,
            on_or_before = days(ref_ar, -1),
            date_format = "YYYY-MM-DD",
        ),
        return_expectations = {"incidence": 0.01},
    ),

    # severe mental illness codes
    sevment_group = patients.satisfying(
        """
        (NOT smhres_date AND sev_mental_date) OR
        smhres_date < sev_mental_date
        """,
        # severe mental illness codes
        sev_mental_date = patients.with_these_clinical_events(
            sev_mental_primis,
            returning = "date",
            find_last_match_in_period = True,
            on_or_before = days(ref_ar, -1),
            date_format = "YYYY-MM-DD",
        ),
        # remission codes relating to severe mental illness
        smhres_date = patients.with_these_clinical_events(
            smhres_primis,
            returning = "date",
            find_last_match_in_period = True,
            on_or_before = days(ref_ar, -1),
            date_format = "YYYY-MM-DD",
        ),
        return_expectations = {"incidence": 0.01},
    ),

    # chronic heart disease codes
    chd_group = patients.with_these_clinical_events(
        chd_primis,
        returning = "binary_flag",
        on_or_before = days(ref_ar, -1),
        return_expectations = {"incidence": 0.01},
    ),

    # chronic kidney disease diagnostic codes
    ckd_group = patients.satisfying(
        """
        ckd OR
        (ckd15_date AND ckd35_date AND (ckd35_date >= ckd15_date)) OR
        (ckd35_date AND NOT ckd15_date)
        """,
        # chronic kidney disease codes - all stages
        ckd15_date = patients.with_these_clinical_events(
            ckd15_primis,
            returning = "date",
            find_last_match_in_period = True,
            on_or_before = days(ref_ar, -1),
            date_format = "YYYY-MM-DD",
        ),
        # chronic kidney disease codes-stages 3 - 5
        ckd35_date = patients.with_these_clinical_events(
            ckd35_primis,
            returning = "date",
            find_last_match_in_period = True,
            on_or_before = days(ref_ar, -1),
            date_format = "YYYY-MM-DD",
        ),
        # chronic kidney disease diagnostic codes
        ckd = patients.with_these_clinical_events(
            ckd_primis,
            returning = "binary_flag",
            on_or_before = days(ref_ar, -1),
        ),
        return_expectations = {"incidence": 0.01},
    ),

    # chronic Liver disease codes
    cld_group = patients.with_these_clinical_events(
        cld_primis,
        returning = "binary_flag",
        on_or_before = days(ref_ar, -1),
        return_expectations = {"incidence": 0.01},
    ),

    # immunosuppressed
    immuno_group = patients.satisfying(
        "immrx OR immdx", 
        # immunosuppression diagnosis codes
        immdx = patients.with_these_clinical_events(
            immdx_primis,
            returning = "binary_flag",
            on_or_before = days(ref_ar, -1),
        ),
        # Immunosuppression medication codes
        immrx = patients.with_these_medications(
            immrx_primis,
            returning = "binary_flag",
            between = [days(ref_ar, -6*30), days(ref_ar, -1)],
        ),
        return_expectations = {"incidence": 0.01},
    ),

    # asplenia or dysfunction of the spleen codes
    spln_group = patients.with_these_clinical_events(
        spln_primis,
        returning = "binary_flag",
        on_or_before = days(ref_ar, -1),
        return_expectations = {"incidence": 0.01},
    ),

    # wider learning disability
    learndis_group = patients.with_these_clinical_events(
        learndis_primis,
        returning = "binary_flag",
        on_or_before = days(ref_ar, -1),
        return_expectations = {"incidence": 0.01},
    ),

    # severe obesity
    sevobese_group = patients.satisfying(
        """
        (sev_obesity_date AND NOT bmi_date) OR
        (sev_obesity_date > bmi_date) OR
        (bmi_value_temp >= 40)
        """,
        bmi_stage_date = patients.with_these_clinical_events(
            bmi_stage_primis,
            returning = "date",
            find_last_match_in_period = True,
            on_or_before = days(ref_ar, -1),
            date_format = "YYYY-MM-DD",
        ),
        sev_obesity_date = patients.with_these_clinical_events(
            sev_obesity_primis,
            returning = "date",
            find_last_match_in_period = True,
            ignore_missing_values = True,
            between = ["bmi_stage_date", days(ref_ar, -1)],
            date_format ="YYYY-MM-DD",
        ),
        bmi_date = patients.with_these_clinical_events(
            bmi_primis,
            returning = "date",
            ignore_missing_values = True,
            find_last_match_in_period = True,
            on_or_before = days(ref_ar, -1),
            date_format = "YYYY-MM-DD",
        ),
        bmi_value_temp = patients.with_these_clinical_events(
            bmi_primis,
            returning = "numeric_value",
            ignore_missing_values = True,
            find_last_match_in_period = True,
            on_or_before = days(ref_ar, -1),
            return_expectations = {
                "float": {"distribution": "normal", "mean": 25, "stddev": 5},
            },
        ),
        return_expectations = {"incidence": 0.01},
    ),

    # at risk group
    atrisk_group=patients.satisfying(
        str_atrisk,
        return_expectations = {
        "incidence": 0.1,
        },
    ),

    # patients in long-stay nursing and residential care before start of phase 1
    longres_group = patients.with_these_clinical_events(
        longres_primis,
        returning = "binary_flag",
        on_or_before = days(start_date, -1),
        return_expectations = {"incidence": 0.01},
    ),

    # # check if pregnant on each of the reference dates
    # **pregnancy("preg_cev_group", ref_cev, type = "value"),
    # **pregnancy("preg_ar_group", ref_ar, type = "value"),

    # derive JCVI group
    jcvi_group = patients.categorised_as(
        dict_jcvi,
        return_expectations = {
            "rate": "universal",
            "incidence": 1,
            "category": { 
                "ratios": ratio_jcvi 
                }
        },
    ),

    # derive vaccine eligibility dates
    elig_date = patients.categorised_as(
       dict_elig,
        return_expectations = {
            "category": {"ratios": 
            ratio_elig
            },
            "incidence": 1,
        },
    ),

)