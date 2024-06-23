import os
from lca_algebraic import initProject
import pandas as pd
import numpy as np
import lca_algebraic as agb
from lca_algebraic import *
from lca_algebraic.stats import * 
import bw2data, bw2io
import bw2calc
from premise import *
import premise as prm
from premise.geomap import Geomap
geomap=Geomap(model='image')

def findActwhRPSSP(name=None, loc=None, multiplerp=False, rp=None, SSP_premise=None, 
                   in_name=None, code=None, categories=None, category=None, 
                   db_name=None, single=True, case_sensitive=False, unit=None, limit=1500) -> ActivityExtended :
    
    geomap = Geomap(model='image') 
    if name and '*' in name:
        in_name = name.replace("*", "")
        name = None

    if not case_sensitive:
        if name :
            name = name.lower()
        if in_name:
            in_name = in_name.lower()

    if SSP_premise == "ssp":
        # for electricity, if not market group, the same location as in ecoinvent is used in premise
        # geomap.ecoinvent_to_iam_location('RoW') will convert RoW to World, but we still need "RoW" for "heat, natural gas"
        """
        delete this if as no need to map UPRs to premise? 
        """
        if loc is not None: 
            if name.startswith("market for electricity"): 
                loc = loc
            else:
                if loc == "GLO":
                    pass
                elif loc == "RoW":
                    pass
                # elif loc.startswith('CA-') and 'QC' not in loc:
                #    loc = 'Canada without Quebec'
                else: 
                    loc = geomap.ecoinvent_to_iam_location(loc)
                    print("changed ecoinvent location to premise:", loc)
        else: 
            print ('location is not specified, using a default one.')
    else:
        pass
        #print("not SSP scenario, using ecoinvent default location:", loc)

    def act_filter(act):
        act_name = act['name']
        if not case_sensitive :
            act_name = act_name.lower()
        if name and not name == act_name:
            return False
        if loc and not loc == act['location']:
            return False
        if unit and not unit == act['unit']:
            return False
        if category and not category in act['categories']:
            return False
        if categories and not tuple(categories) == tuple(act['categories']):
            return False
        return True

    search = name if name is not None else in_name
    search = search.lower()
    search = search.replace(',', ' ')
    
    candidates = bw2data.Database(db_name).search(search, limit=limit)
    #print(candidates[0],type(candidates[0]))
    # Exact match
    acts = list(filter(act_filter, candidates))

    # if multiple reference product for one activity
    if multiplerp == True: 
        candidates_rp = []
        for c in candidates: 
            if c.as_dict()["reference product"] == rp: 
                print(c,type(c))
                candidates_rp.append(c)
        acts = list(filter(act_filter, candidates_rp))

    if single and len(acts) == 0:
        any_name = name if name else in_name
        raise Exception("No activity found in '%s' with name '%s' and location '%s'" % (db_name, any_name, loc))
    if single and len(acts) > 1:
        raise Exception("Several activity found in '%s' with name '%s' and location '%s':\n%s" % (
            db_name, name, loc, "\n".join(str(act) for act in acts)))
    if len(acts) == 1:
        return acts[0]
    else:
        return acts


def findSpecificElecwhSSP(voltage = "medium", loc = None, marketgroup = False, SSP_premise = None, db_name = None, limit=1500 ) :
    def act_filter(act):
        act_name = act['name']
        if name and not name == act_name:
            return False
        if loc and not loc == act['location']:
            return False
        return True

    if marketgroup == True:
        if voltage == "medium":
            name = "market group for electricity, medium voltage"
        elif voltage == "low": 
            name = "market group for electricity, low voltage"
    elif marketgroup == False:
        if voltage == "medium":
            name = "market for electricity, medium voltage"
        elif voltage == "low": 
            name = "market for electricity, low voltage"
             
    elecs =  bw2data.Database(db_name).search(name, limit=limit)
    # Exact match
    acts = list(filter(act_filter, elecs))

    if len(acts) == 0:
        any_name = name if name else in_name
        raise Exception("Check electricity market group vs. non-market, No activity found in '%s' with name '%s' and location '%s'" % (db_name, any_name, loc))
    if len(acts) > 1:
        raise Exception("Several activity found in '%s' with name '%s' and location '%s':\n%s" % (
            db_name, name, loc, "\n".join(str(act) for act in acts)))
    if len(acts) == 1:
        return acts[0]
    else:
        return acts 



def findallElecwhSSP(voltage = "medium", marketgroup = True, SSP_premise = None, db_name = None, limit=1500): 
    """
    to search all available electricity mix (either market group or specific-grid) from ecoinvent / premise-ecoinvent
    return two outputs: 
    1. all available location 
    2. all activity name 
    """
    def act_filter(act):
        act_name = act['name']
        if name and not name == act_name:
            return False
        return True

    if marketgroup == True:
        if voltage == "medium":
            name = "market group for electricity, medium voltage"
        elif voltage == "low": 
            name = "market group for electricity, low voltage"
    elif marketgroup == False:
        if voltage == "medium":
            name = "market for electricity, medium voltage"
        elif voltage == "low": 
            name = "market for electricity, low voltage"
    elecs = bw2data.Database(db_name).search(name, limit=limit)
    # Exact match
    acts = list(filter(act_filter, elecs))
    
    avail_loc = []
    for e in acts: 
        print(e.as_dict()["location"])
        avail_loc.append(e.as_dict()["location"]) 
              
    return (avail_loc, acts)
