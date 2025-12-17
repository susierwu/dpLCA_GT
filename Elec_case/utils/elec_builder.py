"""
Electricity Foreground builder for multiple locatios, mapped to the right premise/year
Purpose:
    - Build (not yet dynamic) electricity foreground nodes.
    - Support multiple IAMs (IMAGE / REMIND ) as long as you have a mapped PREMISE database.
    - Support multiple SSPs and years.
    - Create activities, production exchanges, technosphere exchanges.
    - Return the technosphere exchange for later temporal distribution assignment.
    - Store all generated nodes under the same foreground DB.
    - Build time distribution (td) for the tech exc (just one high voltage elec input)
    - Map the the right background premise_db ad build final datetime 
"""

import bw2data
from bw_temporalis import TemporalDistribution, easy_timedelta_distribution
from datetime import datetime
from bw_timex import TimexLCA


# ============================================================
# HELPER 1 — Build clean foreground activity CODE
# ============================================================
def build_activity_code(location, pathway, year, ref_name):
    """
    pathway can be e.g.:
        'SSP1-26_IMAGE'
        'SSP2-Base_REMIND'
        'SSP5-HighEff_MESSAGE'

    Output example:
        'elec_highvol_CAN_SSP1-26_IMAGE_2030_1kwh'
    """
    return f"{ref_name}_{location}_{pathway}_{year}_1kwh"


# ============================================================
# HELPER 2 — Build human-readable NAME
# ============================================================

def build_activity_name(location, pathway, year, ref_name):
    """
    Output example:
        'electricity mix, high voltage, CAN, SSP1-26_IMAGE, 2030'
    """
    return f"{ref_name}, {location}, {pathway}, {year}"


# ============================================================
# HELPER 3 — Retrieve IAM/premise high-voltage electricity node
# ============================================================
def resolve_premise_db(ssp_variant: str, year: int, prefix="ei_cutoff_3.11_image"):
    """
    Returns the *full* premise database name [str]:
        ei_cutoff_3.11_image_{SSP}-{Variant}_{YEAR} {timestamp}

    Example:
        resolve_premise_db("SSP2-M", 2030)
        -> "ei_cutoff_3.11_image_SSP2-M_2030 2025-11-22"
    Raises:
        ValueError if zero or multiple matches found.
    """
    search_prefix = f"{prefix}_{ssp_variant}_{year}"

    matches = [
        dbname for dbname in bw2data.databases
        if dbname.startswith(search_prefix)
    ]

    if len(matches) != 1:
        raise ValueError(
            f"Multiple or no DB {len(matches)} matches found for pattern '{search_prefix}': {matches}"
        )
    return matches[0]



def get_premise_electricity_node(premise_db,
                                 location = None,
                                 elec_act = None):    
    
    act = premise_db.get(name=elec_act, location=location)
    #DO NOT do len(act) check, if sucessful, returned value is bw2data.backends.proxies.Activity 
    return act



# ============================================================
# HELPER 4 — Create a fresh foreground electricity activity
# ============================================================
def create_foreground_electricity_node(
    foreground_db,
    location,
    pathway,
    year,
    ref_product_prefix
):
    """
    Creates:
        - A foreground activity
        - A production exchange from the matched premise-db
    Returns
        foreground_activity
    """

    code = build_activity_code(location, pathway, year, ref_product_prefix)
    name = build_activity_name(location, pathway, year, ref_product_prefix)

    ref_product = f"{ref_product_prefix}, {location}, {pathway}, {year}"
    # Create new activity
    act = foreground_db.new_node(
        code,
        name=name,
        location=location,
        unit="kWh",
    )
    act["reference product"] = ref_product
   
    # Production exchange
    act.new_edge(input=act, amount=1, type="production").save()
    
    return act


# ============================================================
# HELPER 5 — Add technosphere exchange and RETURN IT
# ============================================================
def add_technosphere_to_foreground(
    fg_activity,
    premise_node,
    amount=1,
):
    """
    Creates and saves a technosphere exchange.
    because temporal_distribution MUST be attached to the exchange, not the activity
    Returns
        exchange_object (so to assign temporal_distribution)
    """
    exc = fg_activity.new_edge(
        input=premise_node,
        amount=amount,
        type="technosphere",
        name=premise_node["name"],
    )
    exc.save()
    fg_activity.save()

    return exc



# ============================================================
# MAIN BUILDER — One call to build a elec FG node
# ============================================================
def build_dynamic_electricity_all(
    locations,
    pathways,
    years,
    elec_act,
    ref_name,
    fg_db_name="elec_foreground",
    flush_fg_db = False
):
    """
    Creates ONE foreground database, then iterates through
    all (location × pathway × year) combinations.
    For each combination:
        • resolve correct premise DB
        • get ONE IAM HV electricity node
        • create 1 FG activity
        • create 1 exchange
        • return results[(loc, pw, yr)] = (fg_act, exchange)
    """

    print("pathways should match the premise database names, e.g., SSP1-VLLO, SSP2-M, SSP5-H \n also IAM locations for premise database can be found through calling: \n from premise import geomap -> premise.geomap.Geomap('image').iam_regions " )
        
    # ---- Create or reuse foreground DB ----
    if flush_fg_db:
        if fg_db_name in bw2data.databases:
            print(f"Flushing existing foreground DB: '{fg_db_name}'")
            del bw2data.databases[fg_db_name]
        foreground = bw2data.Database(fg_db_name)
        foreground.register()
        print(f"Created fresh foreground DB: '{fg_db_name}'")
    else:
        if fg_db_name in bw2data.databases:
            foreground = bw2data.Database(fg_db_name)
            print(
                f"Using foreground DB: {fg_db_name}, with {len(foreground)} activities currently."
            )
        else:
            # DB does not exist yet → create automatically
            foreground = bw2data.Database(fg_db_name)
            foreground.register()
            print(
                f"No existing DB found; created new foreground DB: '{fg_db_name}'"
            )


    results = {}
    # ---- Main loop ----
    for loc in locations:
        for pw in pathways:
            for yr in years:

                # Map premise DB name
                prem_name = resolve_premise_db(pw, yr)
                prem_db = bw2data.Database(prem_name)

                # HV electricity IAM node
                iam_node = get_premise_electricity_node(premise_db = prem_db, 
                                                        location = loc, elec_act = elec_act)

                # New foreground electricity activity
                fg_act = create_foreground_electricity_node(
                    foreground_db = foreground,
                    location = loc,
                    pathway = pw,
                    year = yr,
                    ref_product_prefix = ref_name
                )

                # Add technosphere link and return exchange
                exc = add_technosphere_to_foreground(
                    fg_act,
                    iam_node,
                    amount=1,
                )

                results[(loc, pw, yr)] = (fg_act, exc)

    return results



def assign_td_from_results_dict(
    results_dict,
    elec_td_year=10,
    resolution="Y",
    kind="uniform",
):
    """
    After building all FG activities + exchanges, assign time distribution to technosphere exchange.
    Parameters
    ----------
    results_dict : dict
        Output of build_dynamic_electricity_all()
        keys: (loc, pathway, year)
        values: (fg_act, exchange)
    elec_td_year : int
        Lifetime / horizon for electricity. Default = 10 years.
    resolution : str
        "Y" for year, "M" for month, etc. Must match timex expectations.
    kind : str
        Distribution type ("uniform", "lognormal", etc.)
    """
    steps = elec_td_year + 1
    # SAME distribution for all the elec exchange for all elec market
    td_elec = easy_timedelta_distribution(
        start=0,
        end=elec_td_year,
        resolution=resolution,
        steps=steps,
        kind=kind,
    )
    # assign to each exchange
    for key, (fg_act, exc) in results_dict.items():
        exc["temporal_distribution"] = td_elec
        exc.save()

    return td_elec  # return to plot via td_elec.graph()




def assign_td_from_foreground_db(
    select_act,
    fg_db_name="elec_foreground",
    elec_td_year=10,
    resolution="Y",
    kind="uniform",
    verbose=False,
):
    """
    Assign temporal distribution to ONE selected foreground activity inside the 
    existing 'elec_foreground' database.

    Parameters
    ----------
    select_act : BW2 ActivityProxy
        A single Brightway activity. Required.
    fg_db_name : str
        Foreground DB name.
    elec_td_year : int
        Lifetime horizon.
    resolution : str
        "Y" or "M".
    kind : str
        'uniform', 'increasing', etc.
    """

    # ------------------------------------------
    # 1. Validate DB existence
    # ------------------------------------------
    if fg_db_name not in bw2data.databases:
        raise ValueError(f"Foreground DB '{fg_db_name}' not found.")

    # ------------------------------------------
    # 2. Validate that select_act is an activity
    # ------------------------------------------
    # Brightway ActivityProxy has attributes: .key, .technosphere(), etc.
    if not hasattr(select_act, "technosphere"):
        raise TypeError(
            f"select_act must be a Brightway activity. Got: {type(select_act)}\n"
            f"Value: {select_act}"
        )

    act = select_act  # rename for clarity

    # ------------------------------------------
    # 3. Build TD once
    # ------------------------------------------
    td = easy_timedelta_distribution(
        start=0,
        end=elec_td_year,
        resolution=resolution,
        steps=elec_td_year + 1,
        kind=kind,
    )

    # ------------------------------------------
    # 4. Extract technosphere exchange
    # ------------------------------------------
    tech_excs = list(act.technosphere())
    if len(tech_excs) != 1:
        raise ValueError(
            f"Activity '{act['name']}' has {len(tech_excs)} "
            f"technosphere exchanges. Expected exactly 1."
        )
    exc = tech_excs[0]

    # ------------------------------------------
    # 5. Assign TD
    # ------------------------------------------
    exc["temporal_distribution"] = td
    exc.save()

    if verbose:
        print(f"TD applied to '{act['name']}' → {exc}")

    return exc

    


def map_foreground_to_one_background_build_datetime(
    pathway: str,
    year: int,
    fg_db_name="elec_foreground",
):
    """
    Map a foreground electricity node (pathway + year) to the ONE correct premise background DB.
    If only one mapped background database is selected, runs into KeyError when tlca.build_timeline()
    Returns:
        {
            prem_db_name: datetime(year,1,1),
            fg_db_name: "dynamic"
        }
    """
    # resolve correct premise DB name
    prem_db_name = resolve_premise_db(pathway, year)  
    # build final mapping dictionary
    # flag databases that should be temporally distributed with "dynamic"
    mapping = {
        prem_db_name: datetime.strptime(str(year), "%Y"),
        fg_db_name: "dynamic"
    }
    return mapping

    



def find_dpGWP100_method(pathway, year, 
                    method_prefix = "Climate Change prospective GWP100",
                    method_suffix = "pGWP100 with dp-AGWPCO2"):
 
    """
    Given a pathway ('SSP5-H', 'SSP2-M', 'SSP1-VLLO') and year (2030/2040/2050),
    return the correct dpGWP100 LCIA method tuple.
    Required method format:
        ('Climate Change prospective GWP100',
         SSPcode,       # e.g. 'SSP585'
         MYcode,        # e.g. 'MY2040'
         'pGWP100 with dp-AGWPCO2')
    """

    SSP_TO_METHOD = {
    "SSP5-H":   "SSP585",
    "SSP2-M":   "SSP245",
    "SSP1-VLLO":"SSP119",
    }
    # -------------------------------
    # 1. Map pathway → LCIA SSP code
    # -------------------------------
    if pathway not in SSP_TO_METHOD:
        raise ValueError(
            f"No LCIA method mapping exists for pathway '{pathway}'. "
            f"Supported: {list(SSP_TO_METHOD.keys())}"
        )

    method_ssp_code = SSP_TO_METHOD[pathway]   # e.g. "SSP585"
    method_year_code = f"MY{year}"             # e.g. "MY2040"

    # -------------------------------
    # 2. Search for a matching method
    # -------------------------------
    for method in bw2data.methods:
        if (
            method[0] == method_prefix
            and method[1] == method_ssp_code
            and method[2] == method_year_code
            and method[3] == method_suffix
        ):
            
            return method

    raise ValueError(
        f"Could not find dpGWP LCIA method for pathway={pathway}, year={year}. "
        f"Searched for SSP='{method_ssp_code}', MY='{method_year_code}'."
    )




def run_dp_timex_lca(
    foreground_act,
    pathway = None,
    year = None,
    method = None,
    database_dates = None,
    temporal_grouping="year", 
    method_prefix = "Climate Change prospective GWP100",
    method_suffix = "pGWP100 - dp-AGWPCO2", 
    fg_db_name = 'elec_foreground'
):
    """
    Find correct LCIA method, initialize TimexLCA,
    build timeline and return tlca instance.
    """
    name = foreground_act.get("name")
    if name is None:
        raise ValueError("not a BW2 activity")
    name_parts = [p.strip() for p in name.split(",")]

    if pathway is None: 
        pathway = name_parts[-2]
    if year is None: 
        year_str = name_parts[-1]
        year = int(year_str)
        
    # 1. map the correct dpGWP100 method
    if method is None:
        method = find_dpGWP100_method(pathway, year, method_prefix, method_suffix)
    # 2. map to premise_db and build datetime
    if database_dates is None: 
        database_dates = map_foreground_to_one_background_build_datetime(pathway, year, fg_db_name)

    print(
        f"for the activity {foreground_act},it's under SSP-{pathway}, year-{year} \n" 
        f" we'll use LCIA {method} "
        f" with background database.datetime = {database_dates}\n"
    )
    # 3. Initialize TimexLCA
    tlca = TimexLCA({foreground_act: 1}, method, database_dates)
    # 4. Build timeline
    tlca.build_timeline(temporal_grouping=temporal_grouping)

    return tlca

